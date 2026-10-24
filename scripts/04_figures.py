"""Step 4: maps and charts from the outputs of step 3."""
from pathlib import Path

import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB, FIG = ROOT / "outputs" / "tables", ROOT / "outputs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
GPKG = ROOT / "outputs" / "gis" / "industrial_proximity.gpkg"

BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORANGE = ["#fde3d6", "#f8bfa4", "#f39a73", "#eb6834", "#c9501f", "#9c3c15", "#6e2a0e"]
SEQ_BLUE = LinearSegmentedColormap.from_list("blue", BLUE)
SEQ_ORANGE = LinearSegmentedColormap.from_list("orange", ORANGE)
HOT, COLD, NS = "#e34948", "#2a78d6", "#e4e3df"
REGION_COL = {"Houston (Harris Co.)": "#2a78d6", "Louisiana corridor": "#eb6834"}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 10,
                     "axes.titleweight": "bold", "figure.dpi": 150})


def base(ax, tr):
    tr.plot(ax=ax, color="none", edgecolor="#bdbcb7", linewidth=0.15)
    ax.set_axis_off()


def fac_sizes(v):
    return 4 + 60 * np.sqrt(v / v.max())


def main():
    t = gpd.read_file(GPKG, layer="tracts")
    f = gpd.read_file(GPKG, layer="tri_facilities")
    b = gpd.read_file(GPKG, layer="facility_buffers")
    s = gpd.read_file(GPKG, layer="aqs_pm25_sites")
    regions = list(t["region"].unique())

    # Fig 1: study areas, facilities and 3 km buffers
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2))
    for ax, r in zip(axes, regions):
        tr = t[t["region"] == r]
        base(ax, tr)
        b[(b["region"] == r) & (b["km"] == 3)].plot(ax=ax, color="#eb6834", alpha=0.18, edgecolor="none")
        b[(b["region"] == r) & (b["km"] == 1)].plot(ax=ax, color="#eb6834", alpha=0.35, edgecolor="none")
        fr = f[(f["region"] == r)]
        fr.plot(ax=ax, markersize=fac_sizes(fr["air_lb_mean"].clip(lower=1)), color=INK,
                alpha=0.75, edgecolor="white", linewidth=0.3)
        ax.set_xlim(tr.total_bounds[[0, 2]] + [-3000, 3000])
        ax.set_ylim(tr.total_bounds[[1, 3]] + [-3000, 3000])
        ax.set_title(f"{r}: {len(fr)} TRI facilities with air releases")
    handles = [Patch(color="#eb6834", alpha=0.35, label="within 1 km"),
               Patch(color="#eb6834", alpha=0.18, label="within 3 km"),
               Line2D([], [], marker="o", ls="", color=INK, label="TRI facility (size = mean air lb/yr)")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False)
    fig.suptitle("TRI facilities reporting on-site air releases, 2022-2024", color=INK)
    fig.tight_layout(rect=[0, 0.05, 1, 0.97])
    fig.savefig(FIG / "fig1_facilities_buffers.png")
    plt.close(fig)

    # Fig 2: exposure index; Fig 3: SVI
    for col, cmap, title, fn, lab in [
        ("log_ei_all", SEQ_BLUE, "Distance-weighted air release index, log10(1 + EI)",
         "fig2_exposure_index.png", "log10(1 + lb/yr, distance weighted)"),
        ("RPL_THEMES", SEQ_ORANGE, "CDC/ATSDR SVI 2022 overall percentile (state-relative)",
         "fig3_svi.png", "SVI percentile")]:
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.4))
        vmax = t[col].quantile(0.99)
        for ax, r in zip(axes, regions):
            tr = t[t["region"] == r]
            tr.plot(ax=ax, column=col, cmap=cmap, vmin=0, vmax=vmax, linewidth=0.1, edgecolor="white",
                    missing_kwds={"color": "#f0efec"})
            fr = f[f["region"] == r]
            fr.plot(ax=ax, markersize=3, color=INK, alpha=0.6)
            ax.set_axis_off()
            ax.set_title(r)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, vmax))
        cb = fig.colorbar(sm, ax=axes, orientation="horizontal", fraction=0.04, pad=0.02, shrink=0.5)
        cb.set_label(lab)
        cb.outline.set_visible(False)
        fig.suptitle(title + " (black dots: TRI facilities)", color=INK)
        fig.savefig(FIG / fn, bbox_inches="tight")
        plt.close(fig)

    # Fig 4: Gi* hot spots
    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    cmap = {"hot spot": HOT, "cold spot": COLD, "not significant": NS, "not tested": "#ffffff"}
    for j, (col, lab) in enumerate([("gi_ei_class", "Exposure index"), ("gi_svi_class", "SVI")]):
        for i, r in enumerate(regions):
            ax = axes[j, i]
            tr = t[t["region"] == r]
            tr.plot(ax=ax, color=tr[col].map(cmap), edgecolor="white", linewidth=0.1)
            ax.set_axis_off()
            ax.set_title(f"{lab}: {r}")
    fig.legend(handles=[Patch(color=HOT, label="Gi* hot spot (p < 0.05)"),
                        Patch(color=COLD, label="Gi* cold spot (p < 0.05)"),
                        Patch(color=NS, label="not significant"),
                        Patch(facecolor="white", edgecolor=GRID, label="no population / no SVI")],
               loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Getis-Ord Gi* hot spots (8 nearest neighbours, 999 permutations)", color=INK)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(FIG / "fig4_getis_ord_hotspots.png")
    plt.close(fig)

    # Fig 5: exposed vs unexposed distributions (3 km)
    res = pd.read_csv(TAB / "tract_exposure_results.csv", dtype={"GEOID": str})
    res = res[(res["E_TOTPOP"] > 0) & res["RPL_THEMES"].notna()]
    vars_ = [("RPL_THEMES", "SVI percentile"), ("EP_MINRTY", "% minority"), ("EP_POV150", "% < 150% poverty")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    for ax, (v, lab) in zip(axes, vars_):
        data, pos, cols = [], [], []
        for k, r in enumerate(regions):
            sub = res[res["region"] == r]
            near = sub["near_km"] <= 3
            data += [sub.loc[near, v].dropna(), sub.loc[~near, v].dropna()]
            pos += [k * 3 + 1, k * 3 + 2]
            cols += [REGION_COL[r], "#f0efec"]
        bp = ax.boxplot(data, positions=pos, widths=0.7, patch_artist=True, showfliers=False,
                        medianprops={"color": INK, "linewidth": 1.5})
        for p, c in zip(bp["boxes"], cols):
            p.set_facecolor(c)
            p.set_edgecolor(INK2)
        ax.set_xticks([1, 2, 4, 5], ["<= 3 km", "> 3 km", "<= 3 km", "> 3 km"])
        ax.set_title(lab)
        ax.grid(axis="y", color=GRID, linewidth=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    fig.legend(handles=[Patch(color=REGION_COL[r], label=f"{r}, tract within 3 km of a facility") for r in regions]
               + [Patch(facecolor="#f0efec", edgecolor=INK2, label="tract beyond 3 km")],
               loc="lower center", ncol=3, frameon=False)
    fig.suptitle("Tracts within 3 km of a TRI air-release facility vs other tracts", color=INK)
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(FIG / "fig5_exposed_vs_unexposed.png")
    plt.close(fig)

    # Fig 6: share of population inside buffers
    bp = pd.read_csv(TAB / "buffer_population_svi.csv")
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    w = 0.38
    for k, r in enumerate(regions):
        d = bp[bp["region"] == r]
        x = np.arange(len(d)) + (k - 0.5) * w
        ax.bar(x, d["pct_pop_in_buffer"], width=w - 0.03, color=REGION_COL[r], label=r)
        for xi, v in zip(x, d["pct_pop_in_buffer"]):
            ax.text(xi, v + 1, f"{v:.0f}%", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(np.arange(len(bp["buffer_km"].unique())), [f"{k} km" for k in sorted(bp["buffer_km"].unique())])
    ax.set_ylabel("% of region population (areal weighting)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(frameon=False)
    ax.set_title("Population living within each distance of a TRI air-release facility")
    fig.tight_layout()
    fig.savefig(FIG / "fig6_population_in_buffers.png")
    plt.close(fig)

    # Fig 7: AirToxScreen point-source risk vs exposure index (if present)
    if "ats_cancer_risk_point" in res:
        fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
        for ax, r in zip(axes, regions):
            sub = res[res["region"] == r].dropna(subset=["ats_cancer_risk_point"])
            ax.scatter(sub["log_ei_all"], sub["ats_cancer_risk_point"], s=9, color=REGION_COL[r],
                       alpha=0.6, edgecolor="none")
            ax.set_title(r)
            ax.set_xlabel("log10(1 + exposure index)")
            ax.grid(color=GRID, linewidth=0.5)
            ax.spines[["top", "right"]].set_visible(False)
        axes[0].set_ylabel("AirToxScreen 2020 point-source\nmodeled cancer risk (per million)")
        fig.suptitle("TRI proximity index vs EPA AirToxScreen modeled point-source risk (screening estimate)", color=INK)
        fig.tight_layout()
        fig.savefig(FIG / "fig7_index_vs_airtoxscreen.png")
        plt.close(fig)
    print("figures done")


if __name__ == "__main__":
    main()
