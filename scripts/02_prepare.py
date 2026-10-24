"""Step 2: clean the raw downloads into small analysis files under data/.

Outputs
- data/tri_facilities_air_2022_2024.csv  facility-level on-site air releases
- data/tracts_study_area.gpkg            study-area tracts with SVI 2022 fields
- data/aqs_pm25_annual_2022_2024.csv     PM2.5 annual monitor means in study areas
- data/airtoxscreen2020_tract_cancer_risk.csv  (if the AirToxScreen files exist)
"""
import os
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

RAW = Path(os.environ.get("RAW_DIR", "data_raw")).resolve()
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
YEARS = [2022, 2023, 2024]
CRS_M = 32615  # UTM zone 15N, meters; covers both regions

# Study regions (county / parish FIPS)
REGIONS = {
    "Houston (Harris Co.)": ["48201"],
    "Louisiana corridor": ["22033", "22121", "22047", "22005", "22093",
                           "22095", "22089", "22051", "22071"],
}
# Facilities are kept if they lie within this distance of a study region,
# so tracts near a region edge still see facilities across the line.
FAC_MARGIN_M = 10_000

SVI_COLS = ["FIPS", "E_TOTPOP", "RPL_THEMES", "RPL_THEME1", "RPL_THEME2",
            "RPL_THEME3", "RPL_THEME4", "EP_POV150", "EP_UNEMP", "EP_NOHSDP",
            "EP_UNINSUR", "EP_AGE65", "EP_AGE17", "EP_DISABL", "EP_LIMENG",
            "EP_MINRTY", "EP_NOVEH", "EP_CROWD", "EP_MOBILE"]


def region_of(geoid):
    for name, fips in REGIONS.items():
        if geoid[:5] in fips:
            return name
    return None


def load_tracts():
    parts = []
    for st in ("48", "22"):
        g = gpd.read_file(f"zip://{RAW / 'tiger' / f'tl_2022_{st}_tract.zip'}")
        parts.append(g[["GEOID", "ALAND", "AWATER", "geometry"]])
    g = pd.concat(parts, ignore_index=True)
    g["region"] = g["GEOID"].map(region_of)
    g = g[g["region"].notna()].copy()
    svi = pd.concat([
        pd.read_csv(RAW / "svi" / f"SVI_2022_{s}.csv", dtype={"FIPS": str}, usecols=SVI_COLS)
        for s in ("Texas", "Louisiana")])
    svi["FIPS"] = svi["FIPS"].str.zfill(11)
    num = [c for c in SVI_COLS if c != "FIPS"]
    svi[num] = svi[num].where(svi[num] > -999)  # -999 = missing in SVI
    g = g.merge(svi, left_on="GEOID", right_on="FIPS", how="left").drop(columns="FIPS")
    return gpd.GeoDataFrame(g, geometry="geometry", crs=parts[0].crs).to_crs(CRS_M)


def load_tri(tracts):
    frames = []
    for y in YEARS:
        for st in ("TX", "LA"):
            d = pd.read_csv(RAW / "tri" / f"tri_{y}_{st}.csv", low_memory=False)
            d.columns = [c.split(". ", 1)[1] if ". " in c else c for c in d.columns]
            frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    for c in ("5.1 - FUGITIVE AIR", "5.2 - STACK AIR"):
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0.0)
    # Dioxin-type chemicals are reported in grams; convert to pounds.
    to_lb = np.where(d["UNIT OF MEASURE"].str.lower().str.startswith("gram"), 1 / 453.592, 1.0)
    d["air_lb"] = (d["5.1 - FUGITIVE AIR"] + d["5.2 - STACK AIR"]) * to_lb
    d["air_carc_lb"] = np.where(d["CARCINOGEN"].astype(str).str.upper() == "YES", d["air_lb"], 0.0)
    # Ammonia is the largest single chemical by mass in this area; keep a version without it
    d["air_nonh3_lb"] = np.where(d["CHEMICAL"].astype(str).str.strip().str.lower() == "ammonia", 0.0, d["air_lb"])
    n_rows = len(d)
    by_fy = (d.groupby(["TRIFD", "YEAR"], as_index=False)
             .agg(air_lb=("air_lb", "sum"), air_carc_lb=("air_carc_lb", "sum"),
                  air_nonh3_lb=("air_nonh3_lb", "sum"),
                  n_chem=("CHEMICAL", "nunique")))
    meta = (d.sort_values("YEAR").groupby("TRIFD", as_index=False)
            .agg(name=("FACILITY NAME", "last"), city=("CITY", "last"), county=("COUNTY", "last"),
                 st=("ST", "last"), lat=("LATITUDE", "last"), lon=("LONGITUDE", "last"),
                 sector=("INDUSTRY SECTOR", "last"), naics=("PRIMARY NAICS", "last")))
    fac = (by_fy.groupby("TRIFD", as_index=False)
           .agg(years_reported=("YEAR", "nunique"),
                air_lb_sum=("air_lb", "sum"), air_carc_lb_sum=("air_carc_lb", "sum"),
                air_nonh3_lb_sum=("air_nonh3_lb", "sum"),
                n_chem_max=("n_chem", "max")))
    # Average over the three years (a year with no report counts as zero)
    fac["air_lb_mean"] = fac["air_lb_sum"] / len(YEARS)
    fac["air_carc_lb_mean"] = fac["air_carc_lb_sum"] / len(YEARS)
    fac["air_nonh3_lb_mean"] = fac["air_nonh3_lb_sum"] / len(YEARS)
    for y in YEARS:
        s = by_fy[by_fy["YEAR"] == y].set_index("TRIFD")["air_lb"]
        fac[f"air_lb_{y}"] = fac["TRIFD"].map(s).fillna(0.0)
    fac = fac.merge(meta, on="TRIFD")
    fac = fac.dropna(subset=["lat", "lon"])
    gf = gpd.GeoDataFrame(fac, geometry=gpd.points_from_xy(fac["lon"], fac["lat"]),
                          crs=4326).to_crs(CRS_M)
    # Keep facilities inside or within the margin of each study region
    keep = []
    for name in REGIONS:
        poly = tracts.loc[tracts["region"] == name].union_all().buffer(FAC_MARGIN_M)
        inside = gf.within(poly)
        sub = gf[inside].copy()
        sub["region"] = name
        sub["in_region"] = sub.within(tracts.loc[tracts["region"] == name].union_all())
        keep.append(sub)
    gf = pd.concat(keep, ignore_index=True)
    print(f"TRI rows read: {n_rows:,}; facilities kept: {len(gf)}")
    chem = (d.assign(TRIFD=d["TRIFD"]).merge(gf[["TRIFD", "region"]].drop_duplicates(), on="TRIFD")
            .groupby(["region", "CHEMICAL"], as_index=False)["air_lb"].sum())
    chem["air_lb_mean_per_yr"] = chem["air_lb"] / len(YEARS)
    chem = chem.sort_values(["region", "air_lb"], ascending=[True, False]).groupby("region").head(10)
    chem.drop(columns="air_lb").to_csv(DATA / "tri_top10_chemicals_air.csv", index=False)
    return gf, n_rows


def load_aqs(tracts):
    rows = []
    fips = sum(REGIONS.values(), [])
    for y in YEARS:
        with zipfile.ZipFile(RAW / "aqs" / f"annual_conc_by_monitor_{y}.zip") as z:
            name = [n for n in z.namelist() if n.endswith(".csv")][0]
            a = pd.read_csv(z.open(name), low_memory=False,
                            dtype={"State Code": str, "County Code": str})
        a = a[(a["Parameter Code"] == 88101) &
              (a["Pollutant Standard"] == "PM25 Annual 2012")]
        a["cfips"] = a["State Code"].str.zfill(2) + a["County Code"].str.zfill(3)
        a = a[a["cfips"].isin(fips)]
        a["site_id"] = a["cfips"] + "-" + a["Site Num"].astype(str).str.zfill(4)
        # one value per site-year: prefer the record with most observations
        a = a.sort_values("Observation Count").groupby(["site_id", "Year"], as_index=False).last()
        rows.append(a[["site_id", "Year", "cfips", "Latitude", "Longitude", "Local Site Name",
                       "Arithmetic Mean", "Observation Count", "Completeness Indicator"]])
    out = pd.concat(rows, ignore_index=True)
    out["region"] = out["cfips"].map(region_of)
    return out


def load_airtox():
    """Population-weighted tract means of AirToxScreen 2020 block-level cancer risk.

    EPA rounds the published "Total Cancer Risk" column (for example 30), so the
    tract total is rebuilt as the sum of the unrounded source-group columns.
    """
    files = sorted((RAW / "airtox").glob("Region6*_CancerRisk_by_block_srcgrp.xlsx"))
    if not files:
        return None
    import openpyxl
    fips = tuple(sum(REGIONS.values(), []))
    recs = []
    for f in files:
        wb = openpyxl.load_workbook(f, read_only=True)
        for ws in wb.worksheets:
            it = ws.iter_rows(values_only=True)
            header = [str(v) for v in next(it)]
            if "Block" not in header:
                continue
            i_blk, i_pop = header.index("Block"), header.index("Population")
            i_pt = header.index("PT-StationaryPoint Cancer Risk (per million)")
            i_src = [i for i, h in enumerate(header)
                     if h.endswith("Cancer Risk (per million)") and not h.startswith("Total")]
            for row in it:
                b = row[i_blk]
                if b is None:
                    continue
                b = str(b).split(".")[0].zfill(15)
                if not b.startswith(fips):
                    continue
                tot = sum(float(row[i] or 0) for i in i_src)
                recs.append((b[:11], row[i_pop], tot, row[i_pt]))
    if not recs:
        return None
    b = pd.DataFrame(recs, columns=["GEOID", "pop", "risk_total", "risk_point"])
    b[["pop", "risk_total", "risk_point"]] = b[["pop", "risk_total", "risk_point"]].apply(pd.to_numeric, errors="coerce")
    b = b[b["pop"] > 0]
    b["wt"], b["wp"] = b["risk_total"] * b["pop"], b["risk_point"] * b["pop"]
    agg = b.groupby("GEOID").agg(ats_pop=("pop", "sum"), wt=("wt", "sum"), wp=("wp", "sum"),
                                 ats_blocks=("pop", "size")).reset_index()
    agg["ats_cancer_risk_total"] = agg["wt"] / agg["ats_pop"]
    agg["ats_cancer_risk_point"] = agg["wp"] / agg["ats_pop"]
    print(f"AirToxScreen populated blocks used: {len(b):,}; tracts: {len(agg)}")
    return agg.drop(columns=["wt", "wp"])


def main():
    tracts = load_tracts()
    print(tracts.groupby("region").size())
    fac, _ = load_tri(tracts)
    fac.drop(columns="geometry").to_csv(DATA / "tri_facilities_air_2022_2024.csv", index=False)
    aqs = load_aqs(tracts)
    aqs.to_csv(DATA / "aqs_pm25_annual_2022_2024.csv", index=False)
    print(f"AQS PM2.5 site-years: {len(aqs)}, sites: {aqs['site_id'].nunique()}")
    ats = load_airtox()
    if ats is not None:
        ats.to_csv(DATA / "airtoxscreen2020_tract_cancer_risk.csv", index=False)
    g = tracts.copy()
    g["geometry"] = g.geometry.simplify(20)
    g.to_file(DATA / "tracts_study_area.gpkg", layer="tracts", driver="GPKG")


if __name__ == "__main__":
    main()
