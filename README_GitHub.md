# Kuamut IFM — contrasting paired plots

Site selection for a second wave of forest biomass plots in the Kuamut IFM
project (Sabah, Malaysia, ~84,000 ha). For every existing plot the code places a
new 30 m-diameter plot 100–300 m away in a **contrasting** biomass condition, so
that each pair samples two ends of the local aboveground carbon density (ACD)
range under near-identical conditions of access, terrain and acquisition date.

The pairs are intended to strengthen the plot-to-remote-sensing calibration: the
original 46 plots were placed for accessibility and to span the biomass range,
which leaves the model weakly constrained where ACD changes sharply over short
distances. A paired design puts two measurements either side of that gradient.

Everything runs offline in R from four local layers. No account, API key or
cloud service is needed.

---

## Contents

- [How it works](#how-it-works)
- [Repository layout](#repository-layout)
- [Requirements](#requirements)
- [Input data](#input-data)
- [Quick start](#quick-start)
- [Selection criteria](#selection-criteria)
- [Scoring and assignment](#scoring-and-assignment)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Results of the delivered run](#results-of-the-delivered-run)
- [Known issues and limitations](#known-issues-and-limitations)

---

## How it works

1. **Direction.** Each existing plot is compared with the project-wide ACD
   median. Above it, the code looks for low biomass nearby; below it, high.
2. **Candidates.** A 10 m grid is laid over the 100–300 m annulus around the
   plot, giving roughly 1,000–2,200 admissible points per plot after the hard
   filters (project area, road distance, ACD floor, local homogeneity,
   separation, nearest-neighbour rule).
3. **Score.** Every candidate is scored on contrast, accessibility, local
   homogeneity and walking distance.
4. **Assignment.** Greedy, most-constrained plot first, then a refinement loop
   that re-picks each plot against all the others until nothing improves.
5. **Verification.** The script re-checks the finished design: pair distances,
   minimum separation, and that every new plot's nearest neighbour really is its
   own reference plot.
6. **Reporting.** An overview map, an eight-panel diagnostics figure, per-pair
   zoom panels and a summary table are written alongside the spatial outputs.

---

## Repository layout

```
.
├── README.md
├── LICENSE
├── .gitignore
├── scripts/
│   └── 02_kuamut_paired_plots.R      the whole pipeline: selection + graphics
├── data/                             inputs — NOT tracked
└── outputs/                          results — NOT tracked
```

Project data (the ACD map, plot coordinates, road network, project boundary) is
not in this repository and should not be committed. Suggested `.gitignore`:

```gitignore
data/
outputs/
*.tif
*.shp
*.shx
*.dbf
*.prj
*.cpg
*.gpkg
.Rhistory
.RData
.Rproj.user/
```

---

## Requirements

R ≥ 4.1 (the script uses the native `|>` pipe) with:

```r
install.packages(c("terra", "sf", "dplyr"))
```

`terra` and `sf` need GDAL, PROJ and GEOS. On Windows the CRAN binaries bundle
them; on Linux install `gdal-bin libgdal-dev libproj-dev libgeos-dev` first.
The figures use base graphics only — no extra plotting packages.

---

## Input data

All layers are handled in **EPSG:32650 (WGS 84 / UTM zone 50N)**; the script
reprojects vectors on read, so the raster must already be in that CRS.

| File | Type | Notes |
|---|---|---|
| `acd_2025_30m.tif` | raster, 30 m | Band 1 = ACD in Mg C ha⁻¹. Bands 2–3 (`ci_low`, `ci_upper`) are the map's own confidence interval and are not used in the selection. |
| `KuamutLocations.shp` | points | Existing plot centres. Needs the plot-number field named in `PLOT_ID_FIELD` (default `Pl`). |
| `Kuamut_ProjectArea.shp` | polygon | Project boundary; candidates must fall inside it. |
| `Permian Roads- all.shp` | lines | Road and skid-trail network used for the accessibility term. |

**No river layer was available** when this was built. If you have one, add a
minimum distance-to-river filter to the candidate mask — see
[Known issues](#known-issues-and-limitations).

---

## Quick start

Set the two folders at the top of the script and source it:

```r
DIR_IN  <- "data"        # folder holding the inputs
DIR_OUT <- "outputs"
```

Everything else runs unattended; expect a few minutes, most of it in
`terra::distance()` on the road raster.

> **Before committing:** `DIR_IN` and `DIR_OUT` are absolute Windows paths in the
> working copy. Replace them with relative paths or with `here::here()` before
> pushing — absolute paths break the script for everyone else and leak local
> directory structure.

Console output ends with the verification block:

```
plots total: 92 | closest pair overall: 98.6 m
pair distances: 100-298 m (target 100-300)
road distance of new plots: median 90 m, max 376 m
nearest-neighbour rule: 0 violations out of 46 new plots
```

---

## Selection criteria

| Criterion | Rule | Why |
|---|---|---|
| Pair distance | 100–300 m from the reference plot | Design requirement: same forest type, same access, different condition. |
| Contrast direction | ACD_ref > project median → target low ACD; otherwise target high | Defines what "contrasting" means for each plot. |
| Contrast magnitude | \|ΔACD\| ≥ 40 Mg C ha⁻¹, relaxed through 30 → 20 → 10 → 0 only if infeasible | Below ~40 the difference is not separable from map error. The fallback level reached is recorded per plot. |
| Separation | ≥ 100 m from every other plot centre, existing and new | Avoids overlapping 15 m footprints and heavy spatial autocorrelation. |
| Nearest-neighbour rule | The new plot must be ≥ 25 m closer to its own reference plot than to any other plot | A pair is only a pair if its members are each other's local context. Without this, a new plot can end up beside a different pair and confound the comparison. |
| Road distance | ≥ 25 m, and ≤ max(200 m, the reference plot's own road distance) | The 15 m footprint must not include the road; access is never worse than its pair's. |
| ACD floor | ≥ 10 Mg C ha⁻¹ | Rules out water, bare rock and road platform masquerading as "low biomass". |
| Local homogeneity | sd(ACD) in the 90 × 90 m window ≤ 35 Mg C ha⁻¹ | The 30 m plot should represent its pixel, not straddle an edge. |
| Project area | Inside the boundary | — |

---

## Scoring and assignment

```
score = 0.50 · min(ΔACD / 100, 1)               contrast, saturating
      + 0.20 · (1 − road_distance / 400)        accessibility
      + 0.15 · (1 − local_sd / 30)              homogeneity
      + 0.15 · (1 − (pair_distance − 100)/200)  walking effort
```

**Contrast saturates at 100 Mg C ha⁻¹ deliberately.** Once a large difference is
achieved, the remaining terms decide. Without saturation the optimiser chases the
single most extreme pixel in the annulus, which is almost always on an edge or at
the least convenient spot on the ground.

**Assignment** (`greedy_assign()`) runs in two stages:

1. *Greedy*, taking the reference plot with the fewest admissible candidates
   first, so the most constrained plots are not squeezed out by later ones.
2. *Refinement*, up to `ROUNDS` passes: each plot, worst contrast first, is
   re-picked against the current choices of all the others and moved only if its
   score improves. This removes the dependence on the arbitrary greedy order.

The nearest-neighbour rule is enforced in both places: against the existing plots
when candidates are generated, and mutually between new plots during assignment —
a candidate is rejected if it would end up closer to an already-placed new plot
than either of them is to its own reference plot.

---

## Parameters

All at the top of `02_kuamut_paired_plots.R`.

| Parameter | Default | Meaning |
|---|---|---|
| `PAIR_MIN`, `PAIR_MAX` | 100, 300 | m, admissible distance to the reference plot |
| `MIN_SEP` | 100 | m, minimum distance to any other plot centre |
| `NN_MARGIN` | 25 | m, how much closer a new plot must be to its own pair than to anything else |
| `ROAD_MIN` | 25 | m, minimum distance to a road |
| `ROAD_CAP_FLOOR` | 200 | m, road-distance allowance for plots that sit on a road |
| `SD_HARD` | 35 | Mg C ha⁻¹, maximum ACD sd in the 90 m window |
| `SD_SCALE` | 30 | Mg C ha⁻¹, scaling of the homogeneity score term |
| `MIN_ACD` | 10 | Mg C ha⁻¹, candidate ACD floor |
| `CONTRAST_SAT` | 100 | Mg C ha⁻¹, contrast saturation point |
| `DELTA_LEVELS` | 40, 30, 20, 10, 0 | contrast ladder, tried in order |
| `GRID` | 10 | m, candidate search grid |
| `ROUNDS` | 8 | maximum refinement passes |
| `PLOT_ID_FIELD` | `"Pl"` | plot-number field in the plot shapefile |
| `EPSG` | 32650 | working CRS |
| `W_CONTRAST`, `W_ACCESS`, `W_HOMOG`, `W_PROX` | 0.50, 0.20, 0.15, 0.15 | score weights, must sum to 1 |

Sensitivity worth knowing: `NN_MARGIN` at 0, 10, 25 and 50 m gives an **identical
design**. The nearest-neighbour rule itself is binding, its tolerance is not, so
25 m is free insurance against GPS error.

---

## Outputs

Written to `DIR_OUT`, all stemmed `Kuamut_paired_plots_2025`:

**Spatial**

| File | Contents |
|---|---|
| `*.gpkg` | layers `paired_plots`, `existing_plots`, `pair_links`, `new_plot_footprints_15m` (UTM 50N) |
| `*.kml` | WGS84 points for Google Earth / Avenza / field navigation |
| `*.csv` | full attribute table with UTM and lon/lat |

Attribute columns: `Pl_new`, `Pl_ref`, `target`, `cls_ref`/`cls_new` (ACD
tercile class), `ACD_ref`, `ACD_new`, `dACD`, `d_pair`, `d_other` (distance to
the nearest other plot), `road_ref`, `road_new`, `acd_sd`, `score`, `x`/`y`,
`lon`/`lat`. Naming convention: the pair of plot `20` is `20B`.

**Figures and report**

| File | Contents |
|---|---|
| `*_map.png` | Overview: ACD raster, project boundary, roads, both plot sets, pair links, labels, scale bar and north arrow |
| `*_diagnostics.png` | Eight panels: (a) where plots sit in the project ACD distribution, (b) per-pair reference → paired shift, (c) tercile shares of area vs plots, (d) plots per ACD decile, (e) contrast achieved per pair, (f) plot-level ACD boxplots, (g) distance-to-road ECDF old vs new, (h) pair distance histogram with local ACD sd overlaid |
| `*_pairs_01.png`, `_02.png`, … | Per-pair zoom panels, 12 per page, each a 500 m window on the ACD map |
| `*_summary.csv` | The numbers behind the figures: tercile shares, decile occupancy, ACD-space coverage, contrast and access statistics |

Two coverage metrics in the summary are worth reading together. **Deciles
occupied** counts how many of the ten equal-area ACD deciles contain at least one
plot. **ACD space covered** is the share of the project area whose ACD sits within
one 5 Mg C ha⁻¹ bin of at least one plot value — a stricter measure of how much of
the biomass range the calibration data actually constrains.

---

## Results of the delivered run

46 existing plots (the source shapefile holds 46, not 49 — `Pl` = 1…46, no
duplicates or nulls), 46 pairs, 92 plots total.

| Metric | Value |
|---|---|
| Pairs placed | 46 / 46 |
| Meeting ΔACD ≥ 40 | 44; plot 14 falls to the 30 level (37.8), plot 1 to the 10 level (18.8) |
| \|ΔACD\| | median 87, mean 80, max 118 Mg C ha⁻¹ |
| Pair distance | 100–298 m, median 190 m |
| Nearest-neighbour rule | 0 violations; margin 184 m at the tightest, 850 m median |
| Road distance, new plots | median 90 m, max 376 m (existing plots: median 173 m, max 498 m) |
| Local ACD sd at new plots | median 15, max 34 Mg C ha⁻¹ |
| Tercile transitions | 32 pairs cross high↔low; 14 involve the middle tercile |

The 14 pairs involving the middle class are reference plots that already sat
mid-distribution, where a full opposite contrast does not exist within 300 m.

Plots 1B and 14B sit in the dense south-east cluster, where the high-contrast
ground is closer to a neighbouring plot than to their own reference plot and is
therefore excluded by the nearest-neighbour rule. Before that rule was added both
reached ΔACD ≥ 40 and the median was 90. If 1B's weak contrast matters more to you
than the rule, that single plot is the place to make a manual exception.

---

## Known issues and limitations

1. **Road distance is enforced on a 30 m raster.** `terra::distance()` runs on
   the ACD grid, so the `ROAD_MIN = 25 m` filter is only accurate to about half a
   pixel. In the delivered run five plots have a true distance below 25 m, the
   closest at 7 m — close enough that the 15 m footprint may clip the verge. Fix
   by post-filtering the chosen points with exact geometry:
   ```r
   new$road_exact <- as.numeric(st_distance(
     st_as_sf(new, coords = c("x","y"), crs = EPSG), st_union(roads)))
   ```
   and re-picking any plot below `ROAD_MIN`, or by computing the distance raster
   at 10 m instead of 30 m.
2. **No river layer.** A low-ACD candidate could be a river bar or an eroding
   bank. Add `river_dist >= 30` to the candidate mask when the layer is available
   — one line alongside the existing road filter.
3. **No terrain screening.** The pipeline takes no DEM, so nothing stops a plot
   landing on a cliff or well above its reference plot. To add it, read a clipped
   DEM with `terra::rast()`, derive slope with `terra::terrain()`, and add
   `elev <= cap`, `slope <= 30` and `abs(elev - elev_ref) <= 80` to the candidate
   mask. A bare-earth DEM is worth the trouble here: surface models still carry
   the canopy over closed forest, so their slope partly describes the treetops
   rather than the ground the crews walk on.
4. **The ACD map is modelled.** `ci_low` / `ci_upper` give the map's own interval
   at each point and are wide in several cases. Check a handful of picks against
   high-resolution imagery before committing a field campaign.
5. **Design objective.** The score maximises local contrast. If the real goal is
   to reduce calibration error over a specific range — say ACD > 150, where plots
   are scarce — replace the contrast term with a weight on under-represented ACD
   bins. That is two lines in the score, and it is worth deciding before fieldwork
   rather than after.

---

## Licence

Add a `LICENSE` file for the code before publishing; the repository has none at
the time of writing. Project data is not distributed here.
