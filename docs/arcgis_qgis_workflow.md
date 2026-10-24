# Reproducing the GIS steps in ArcGIS Pro or QGIS

The analysis in this repository was run in Python with open-source tools
(GeoPandas, Shapely, PySAL/esda). ArcGIS Pro was **not** run for this
repository. The steps below let you repeat the main GIS operations in
ArcGIS Pro or QGIS and check them against the Python results. Treat any
ArcGIS result as yours only after you have run it yourself.

## Inputs

- `outputs/gis/industrial_proximity.gpkg` (opens directly in ArcGIS Pro 3.x and QGIS)
  - `tracts`: study-area census tracts with SVI fields, exposure index, Gi* results
  - `tri_facilities`: TRI facilities with mean 2022-2024 on-site air releases (lb/yr)
  - `facility_buffers`: dissolved 1, 3 and 5 km buffers per region
  - `aqs_pm25_sites`: AQS PM2.5 monitors with 2022-2024 mean annual concentration
- `outputs/gis/tri_facilities.kml` (Google Earth Pro or ArcGIS "KML To Layer")

All layers are in WGS 84 / UTM zone 15N (EPSG:32615), units in meters,
except the KML (WGS 84 geographic).

## ArcGIS Pro (not run here)

1. Add the GeoPackage layers to a new map (Catalog, right-click the .gpkg, Add To Current Map).
2. Buffers: Analysis > Pairwise Buffer on `tri_facilities`, distances 1, 3 and 5 km,
   Dissolve Type = Dissolve all. Compare with the `facility_buffers` layer.
3. Population within buffers: Tabulate Intersection (zone = buffer, class = `tracts`,
   sum field `E_TOTPOP`), then multiply population by `PERCENTAGE / 100`. This is the
   same areal weighting used in `scripts/03_analysis.py`. Python results are in
   `outputs/tables/buffer_population_svi.csv`.
4. Distance to nearest facility: Near (Analysis) from tract centroids (Feature To Point,
   Inside option) to `tri_facilities`. Compare with field `near_km`.
5. Hot spots: Spatial Statistics > Hot Spot Analysis (Getis-Ord Gi*) on `tracts`,
   field `log_ei_all` (then `RPL_THEMES`), Conceptualization = K nearest neighbors,
   Number of Neighbors = 8, and run once per region (select by `region` first).
   ArcGIS uses analytical z-score p-values (optionally FDR), while the Python run uses
   999 conditional permutations, so expect small differences in which tracts reach
   p < 0.05.
6. Symbolize `gi_ei_class` and `gi_svi_class` to compare with `outputs/figures/fig4_getis_ord_hotspots.png`.

## QGIS (free)

1. Layer > Add Layer > Add Vector Layer, choose the .gpkg and all sublayers.
2. Buffers: Vector > Geoprocessing Tools > Buffer, distance 1000/3000/5000, Dissolve result.
3. Overlay: Processing Toolbox > Overlap analysis (tracts vs buffer) gives the area
   share used for population weighting.
4. Hot spots: install the "Hotspot Analysis" plugin (uses PySAL) and run Gi* with
   k nearest neighbours = 8 on `log_ei_all`.

## Google Earth Pro

Open `outputs/gis/tri_facilities.kml`. Each placemark lists the region and the mean
2022-2024 on-site air releases.
