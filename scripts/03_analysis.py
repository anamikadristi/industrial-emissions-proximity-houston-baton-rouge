"""Step 3: proximity and exposure-index analysis.

Exposure here means residential proximity to TRI facilities and their
reported air releases. It is not a measured dose and says nothing about
health outcomes.

Outputs go to outputs/tables and outputs/gis.
"""
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from esda.getisord import G_Local
from libpysal.weights import KNN
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA, TAB, GIS = ROOT / "data", ROOT / "outputs" / "tables", ROOT / "outputs" / "gis"
TAB.mkdir(parents=True, exist_ok=True)
GIS.mkdir(parents=True, exist_ok=True)
CRS_M = 32615
BUFFERS_KM = [1, 3, 5]
PRIMARY_KM = 3
DECAY_KM = 2.0     # e-folding distance of the exponential distance weight
CUTOFF_KM = 10.0   # facilities farther than this contribute zero
KNN_K = 8
SEED = 20261023

COMPARE_VARS = {
    "RPL_THEMES": "SVI overall percentile (state-relative)",
    "EP_MINRTY": "% racial/ethnic minority",
    "EP_POV150": "% below 150% poverty",
    "EP_NOHSDP": "% no high school diploma",
    "EP_UNINSUR": "% uninsured",
    "EP_LIMENG": "% limited English",
    "EP_NOVEH": "% households with no vehicle",
    "EP_AGE17": "% age 17 and under",
    "EP_AGE65": "% age 65 and over",
}


def load():
    t = gpd.read_file(DATA / "tracts_study_area.gpkg", layer="tracts")
    f = pd.read_csv(DATA / "tri_facilities_air_2022_2024.csv", dtype={"TRIFD": str})
    f = gpd.GeoDataFrame(f, geometry=gpd.points_from_xy(f["lon"], f["lat"]), crs=4326).to_crs(CRS_M)
    ats_path = DATA / "airtoxscreen2020_tract_cancer_risk.csv"
    if ats_path.exists():
        ats = pd.read_csv(ats_path, dtype={"GEOID": str})
        t = t.merge(ats[["GEOID", "ats_cancer_risk_total", "ats_cancer_risk_point"]], on="GEOID", how="left")
        COMPARE_VARS["ats_cancer_risk_total"] = "AirToxScreen 2020 modeled cancer risk, all sources (per million)"
        COMPARE_VARS["ats_cancer_risk_point"] = "AirToxScreen 2020 modeled cancer risk, point sources (per million)"
    return t, f


def exposure_index(pts, fac, col):
    """Sum over facilities of air release x exp(-d / DECAY), d <= CUTOFF."""
    px = np.c_[pts.x.values, pts.y.values]
    fx = np.c_[fac.geometry.x.values, fac.geometry.y.values]
    d = np.sqrt(((px[:, None, :] - fx[None, :, :]) ** 2).sum(-1)) / 1000.0
    w = np.exp(-d / DECAY_KM) * (d <= CUTOFF_KM)
    return w @ fac[col].values, d.min(axis=1)


def rank_biserial(x, y):
    u = stats.mannwhitneyu(x, y, alternative="two-sided")
    r = 2 * u.statistic / (len(x) * len(y)) - 1  # equals Cliff's delta
    return u.statistic, u.pvalue, r


def bh(p):
    p = np.asarray(p, float)
    n = np.sum(~np.isnan(p))
    order = np.argsort(np.where(np.isnan(p), np.inf, p))
    q = np.full_like(p, np.nan)
    prev = 1.0
    for rank, i in reversed(list(enumerate(order[:n], start=1))):
        prev = min(prev, p[i] * n / rank)
        q[i] = prev
    return q


def main():
    t, f = load()
    f_air = f[f["air_lb_mean"] > 0].copy()
    regions = list(t["region"].unique())
    pts = t.geometry.representative_point()

    # ---- exposure index and nearest-facility distance per tract
    t["ei_all"] = 0.0
    t["ei_carc"] = 0.0
    t["near_km"] = np.nan
    for r in regions:
        m = (t["region"] == r).values
        fr = f_air[f_air["region"] == r]
        ei, dmin = exposure_index(pts[m], fr, "air_lb_mean")
        eic, _ = exposure_index(pts[m], fr, "air_carc_lb_mean")
        ein, _ = exposure_index(pts[m], fr, "air_nonh3_lb_mean")
        t.loc[m, "ei_all"], t.loc[m, "ei_carc"], t.loc[m, "near_km"] = ei, eic, dmin
        t.loc[m, "ei_nonh3"] = ein
    t["log_ei_all"] = np.log10(1 + t["ei_all"])
    t["log_ei_carc"] = np.log10(1 + t["ei_carc"])
    t["log_ei_nonh3"] = np.log10(1 + t["ei_nonh3"])
    valid = (t["E_TOTPOP"] > 0) & t["RPL_THEMES"].notna()

    # ---- facility summary
    fs = (f_air.groupby(["region", "in_region"])
          .agg(n_facilities=("TRIFD", "size"), air_lb_mean_total=("air_lb_mean", "sum"),
               air_carc_lb_mean_total=("air_carc_lb_mean", "sum"),
               air_nonh3_lb_mean_total=("air_nonh3_lb_mean", "sum"),
               air_lb_2022=("air_lb_2022", "sum"), air_lb_2023=("air_lb_2023", "sum"),
               air_lb_2024=("air_lb_2024", "sum")).reset_index())
    fs.to_csv(TAB / "facility_summary.csv", index=False)
    top = (f_air[f_air["in_region"]].sort_values("air_lb_mean", ascending=False)
           .groupby("region").head(10)[["region", "name", "city", "sector", "air_lb_mean", "air_carc_lb_mean"]])
    top.to_csv(TAB / "top10_facilities_air_releases.csv", index=False)

    # ---- buffer analysis (areal weighting of tract population)
    t["area"] = t.geometry.area
    rows, buffers = [], []
    for r in regions:
        tr = t[(t["region"] == r)]
        fr = f_air[f_air["region"] == r]
        for km in BUFFERS_KM:
            buf = fr.buffer(km * 1000).union_all()
            buffers.append({"region": r, "km": km, "geometry": buf})
            frac = (tr.geometry.intersection(buf).area / tr["area"]).clip(0, 1)
            t.loc[tr.index, f"frac_{km}km"] = frac.values
            pop = tr["E_TOTPOP"].fillna(0)
            w_in, w_out = pop * frac, pop * (1 - frac)
            ok = tr["RPL_THEMES"].notna()
            row = {"region": r, "buffer_km": km, "n_facilities": len(fr),
                   "pop_total": pop.sum(), "pop_in_buffer": w_in.sum(),
                   "pct_pop_in_buffer": 100 * w_in.sum() / pop.sum()}
            for v in ["RPL_THEMES", "EP_MINRTY", "EP_POV150", "EP_NOHSDP", "EP_LIMENG"]:
                ok = tr[v].notna()
                row[f"{v}_popwt_in"] = np.average(tr.loc[ok, v], weights=w_in[ok]) if w_in[ok].sum() > 0 else np.nan
                row[f"{v}_popwt_out"] = np.average(tr.loc[ok, v], weights=w_out[ok]) if w_out[ok].sum() > 0 else np.nan
            rows.append(row)
    pd.DataFrame(rows).to_csv(TAB / "buffer_population_svi.csv", index=False)
    buffers = gpd.GeoDataFrame(buffers, geometry="geometry", crs=CRS_M)

    # ---- exposed vs unexposed tracts (Wilcoxon rank-sum, rank-biserial r)
    rows = []
    for r in regions:
        for km in BUFFERS_KM:
            sub = t[valid & (t["region"] == r)]
            exp_ = sub["near_km"] <= km
            for v, lab in COMPARE_VARS.items():
                x, y = sub.loc[exp_, v].dropna(), sub.loc[~exp_, v].dropna()
                if len(x) < 5 or len(y) < 5:
                    continue
                U, p, rb = rank_biserial(x, y)
                rows.append({"region": r, "exposed_def_km": km, "variable": v, "label": lab,
                             "n_exposed": len(x), "n_unexposed": len(y),
                             "median_exposed": x.median(), "median_unexposed": y.median(),
                             "median_diff": x.median() - y.median(), "U": U, "p": p,
                             "rank_biserial_r": rb})
    wt = pd.DataFrame(rows)
    wt["q_bh"] = bh(wt["p"])
    wt.to_csv(TAB / "exposed_vs_unexposed_wilcoxon.csv", index=False)

    # ---- correlations of the exposure index with SVI and demographics
    rows = []
    for r in regions:
        sub = t[valid & (t["region"] == r)]
        for ei in ["log_ei_all", "log_ei_carc", "log_ei_nonh3"]:
            for v in COMPARE_VARS:
                s = sub[[ei, v]].dropna()
                rho, p = stats.spearmanr(s[ei], s[v])
                rows.append({"region": r, "index": ei, "variable": v, "n": len(s), "spearman_rho": rho, "p": p})
    ct = pd.DataFrame(rows)
    ct["q_bh"] = bh(ct["p"])
    ct.to_csv(TAB / "exposure_index_spearman.csv", index=False)

    # ---- region-level contrasts (raw estimates; SVI ranks are state-relative)
    rows = []
    a, b = [t[valid & (t["region"] == r)] for r in regions]
    for v in ["log_ei_all", "log_ei_carc", "log_ei_nonh3", "near_km"] + [k for k in COMPARE_VARS if k != "RPL_THEMES"]:
        x, y = a[v].dropna(), b[v].dropna()
        U, p, rb = rank_biserial(x, y)
        rows.append({"variable": v, "region_a": regions[0], "region_b": regions[1],
                     "median_a": x.median(), "median_b": y.median(), "n_a": len(x), "n_b": len(y),
                     "p": p, "rank_biserial_r_a_vs_b": rb})
    rc = pd.DataFrame(rows)
    rc["q_bh"] = bh(rc["p"])
    rc.to_csv(TAB / "region_contrast.csv", index=False)

    # ---- Getis-Ord Gi* hot spots (KNN weights on tract points, per region)
    rows = []
    for c in ["gi_ei_class", "gi_svi_class"]:
        t[c] = "not tested"
    for r in regions:
        sub = t[valid & (t["region"] == r)]
        w = KNN.from_array(np.c_[pts[sub.index].x, pts[sub.index].y], k=KNN_K)
        w.transform = "B"
        for v, c in [("log_ei_all", "gi_ei"), ("RPL_THEMES", "gi_svi")]:
            g = G_Local(sub[v].values, w, star=True, permutations=999, seed=SEED)
            q = bh(g.p_sim)
            cls = np.where(g.p_sim >= 0.05, "not significant", np.where(g.Zs > 0, "hot spot", "cold spot"))
            t.loc[sub.index, f"{c}_z"] = g.Zs
            t.loc[sub.index, f"{c}_p"] = g.p_sim
            t.loc[sub.index, f"{c}_q"] = q
            t.loc[sub.index, f"{c}_class"] = cls
            for k in ["hot spot", "cold spot", "not significant"]:
                rows.append({"region": r, "variable": v, "class": k, "n_tracts": int((cls == k).sum()),
                             "n_fdr_q_lt_0.05": int(((cls == k) & (q < 0.05)).sum())})
        both = (t.loc[sub.index, "gi_ei_class"] == "hot spot") & (t.loc[sub.index, "gi_svi_class"] == "hot spot")
        rows.append({"region": r, "variable": "both", "class": "hot spot for both", "n_tracts": int(both.sum()),
                     "n_fdr_q_lt_0.05": np.nan})
        t.loc[sub.index, "gi_both_hot"] = both.values
    pd.DataFrame(rows).to_csv(TAB / "getis_ord_summary.csv", index=False)

    # ---- AQS PM2.5 monitors: exposure index at the monitor location
    aqs = pd.read_csv(DATA / "aqs_pm25_annual_2022_2024.csv", dtype={"cfips": str})
    site = (aqs.groupby(["site_id", "region", "Local Site Name"], as_index=False)
            .agg(lat=("Latitude", "first"), lon=("Longitude", "first"),
                 pm25_mean_2022_2024=("Arithmetic Mean", "mean"), n_years=("Year", "nunique")))
    sg = gpd.GeoDataFrame(site, geometry=gpd.points_from_xy(site["lon"], site["lat"]), crs=4326).to_crs(CRS_M)
    for r in regions:
        m = (sg["region"] == r).values
        if m.any():
            ei, dmin = exposure_index(sg.geometry[m], f_air[f_air["region"] == r], "air_lb_mean")
            sg.loc[m, "ei_all"], sg.loc[m, "near_km"] = ei, dmin
    sg["log_ei_all"] = np.log10(1 + sg["ei_all"])
    sg.drop(columns="geometry").to_csv(TAB / "aqs_pm25_sites_with_exposure_index.csv", index=False)
    rows = []
    for r, sub in [("both regions", sg)] + [(r, sg[sg["region"] == r]) for r in regions]:
        rho, p = stats.spearmanr(sub["log_ei_all"], sub["pm25_mean_2022_2024"])
        rows.append({"region": r, "n_sites": len(sub), "pm25_median": sub["pm25_mean_2022_2024"].median(),
                     "spearman_rho": rho, "p": p})
    pd.DataFrame(rows).to_csv(TAB / "aqs_pm25_vs_index.csv", index=False)

    # ---- per-tract table and GIS layers
    keep = ["GEOID", "region", "E_TOTPOP", "RPL_THEMES", "EP_MINRTY", "EP_POV150", "near_km",
            "ei_all", "ei_carc", "ei_nonh3", "log_ei_all", "log_ei_carc", "log_ei_nonh3", "frac_1km", "frac_3km", "frac_5km",
            "gi_ei_z", "gi_ei_p", "gi_ei_class", "gi_svi_z", "gi_svi_p", "gi_svi_class", "gi_both_hot"]
    keep += [c for c in ["ats_cancer_risk_total", "ats_cancer_risk_point"] if c in t]
    t[keep].to_csv(TAB / "tract_exposure_results.csv", index=False)
    gpkg = GIS / "industrial_proximity.gpkg"
    if gpkg.exists():
        gpkg.unlink()
    tg = t[keep + ["geometry"]].copy()
    tg["gi_both_hot"] = tg["gi_both_hot"].fillna(False).astype(int)
    tg.to_file(gpkg, layer="tracts", driver="GPKG")
    f_air[["TRIFD", "name", "city", "region", "in_region", "sector", "air_lb_mean", "air_carc_lb_mean",
           "air_lb_2022", "air_lb_2023", "air_lb_2024", "geometry"]].to_file(gpkg, layer="tri_facilities", driver="GPKG")
    buffers.to_file(gpkg, layer="facility_buffers", driver="GPKG")
    sg.to_file(gpkg, layer="aqs_pm25_sites", driver="GPKG")
    kml = f_air[["name", "region", "air_lb_mean", "geometry"]].to_crs(4326).rename(columns={"name": "Name"})
    kml["Description"] = kml.apply(lambda x: f"{x['region']}; mean on-site air releases 2022-2024: {x['air_lb_mean']:,.0f} lb/yr", axis=1)
    kml[["Name", "Description", "geometry"]].to_file(GIS / "tri_facilities.kml", driver="KML")
    print("analysis done")


if __name__ == "__main__":
    main()
