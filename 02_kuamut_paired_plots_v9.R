# =============================================================================
#  Kuamut IFM (Sabah, MY) - contrasting paired plots
#
#  MODE = "local" : full pipeline from the local files (ACD map + roads +
#                   project area + existing plots). No Earth Engine needed.
#                   This reproduces Kuamut_paired_plots_v9.* exactly.
#  MODE = "gee"   : read the candidate table exported by
#                   01_gee_kuamut_paired_plots_v9.js (which adds the AlphaEarth
#                   Satellite Embedding dissimilarity to the score) and run the
#                   final greedy assignment that enforces the 100 m separation
#                   between the new plots.
#
#  Requires: terra, sf, dplyr
# =============================================================================
#install.packages("terra")
library(terra); library(sf); library(dplyr)

SEED <- 123          # nothing here draws random numbers, but the seed is set so
set.seed(SEED)       #   the script stays reproducible if anything is added, and
                     #   every ordering below has an explicit tie-break

MODE <- "local"        # "local" or "gee"

# ------------------------------- paths ---------------------------------------
DIR_IN   <- "C:/Users/JavierRuizRamos/OneDrive - Permian Global Research Limited/Desktop/Kuamut/VegetationPlot campaing_20252026/Kuamut_PlotLocationsCheck_092026/KuamutPlotPairs_Claude/CodebaseData"                     # folder holding the inputs
DIR_OUT  <- "C:/Users/JavierRuizRamos/OneDrive - Permian Global Research Limited/Desktop/Kuamut/VegetationPlot campaing_20252026/Kuamut_PlotLocationsCheck_092026/KuamutPlotPairs_Claude/CodeRResults"
# one ACD map per year, all on the same 30 m grid
F_ACD_YEARS <- c("2021" = file.path(DIR_IN, "acd_2021_30.tif"),
                 "2022" = file.path(DIR_IN, "acd_2022_30.tif"),
                 "2023" = file.path(DIR_IN, "acd_2023_30.tif"),
                 "2024" = file.path(DIR_IN, "acd_2024_30m.tif"),
                 "2025" = file.path(DIR_IN, "acd_2025_30m.tif"))
F_ACD    <- F_ACD_YEARS[["2025"]]      # grid reference and the display map
F_PLOTS  <- file.path(DIR_IN, "KuamutLocations.shp")
F_AOI    <- file.path(DIR_IN, "Kuamut_ProjectArea.shp")
F_ROADS  <- file.path(DIR_IN, "Permian Roads- all.shp")
F_NROAD  <- file.path(DIR_IN, "kuamutnorthroad.shp")   # the accessible north road
F_DEM    <- file.path(DIR_IN, "Kuamut_Elevation.tif")  # elevation, any CRS
NROAD_EPSG <- 32650                                    # the file ships with no .prj
F_GEECSV <- file.path(DIR_IN, "Kuamut_paired_plot_candidates_v9.csv")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

# ----------------------------- parameters ------------------------------------
PAIR_MIN       <- 100    # m   distance to the reference plot
PAIR_MAX       <- 300    # m
PLOT_RADIUS    <- 30     # m   plot radius (60 m diameter). Everything that
                         #     depends on plot size is derived from this.
MIN_SEP        <- 100    # m   min distance to any other plot centre; leaves
                         #     MIN_SEP - 2*PLOT_RADIUS between footprint edges
ROAD_MIN       <- PLOT_RADIUS + 10   # m  keep the footprint clear of the road
ROAD_CAP_FLOOR <- 200    # m   road-distance allowance for roadside plots
SD_HARD        <- 35     # Mg C/ha  max ACD sd in the 90 m neighbourhood
SD_SCALE       <- 30
MIN_ACD        <- 10     # Mg C/ha  floor: avoids water / bare / road pixels
CONTRAST_SAT   <- 100    # Mg C/ha  contrast reward saturates here
DELTA_LEVELS   <- c(40, 30, 20, 10, 0)   # required |dACD|, relaxed in turn
GRID           <- 10     # m   candidate search grid
NN_MARGIN      <- 25     # m   a new plot must be at least this much closer to
                         #     its own reference plot than to any other plot
ROUNDS         <- 12     # refinement passes after the greedy pass
TOL            <- 0.01   # min score gain to accept a move in refinement; stops
                         #   the coordinate descent oscillating on near-ties
# --- v3: spread the WHOLE sample across the ACD range, tails first -----------
DIRECTION      <- "free" # "free"   : contrast either way, the spread objective
                         #            decides which side of the reference plot
                         # "median" : v1/v2 rule, references above the project
                         #            median always look for lower biomass
# --- v4: temporal stability across the five annual maps ---------------------
RS_MAX         <- 0.20   # max 1.4826*MAD/median across years (relative spread)
SLOPE_MAX      <- 6      # Mg C/ha/yr, max |linear trend| 2021-2025: excludes
                         #   ground that was logged or is regrowing
RS_FLOOR       <- 20     # Mg C/ha, floor in the relative spread so genuinely
                         #   low-biomass ground is not penalised for small wobble
DELTA_YEAR_MIN <- 20     # Mg C/ha, the contrast must reach this in EVERY year,
                         #   with a consistent sign
# --- v6: brothers along the north road for the hard-to-reach plots ----------
HARD_PLOTS     <- c(23, 24, 32, 33, 37, 38, 41, 42, 45, 46)
NROAD_MAX      <- 1500   # m, a brother must lie within this of the north road
NROAD_PREF     <- 500    # m, preferred stand-off from the road
NROAD_PREF_SD  <- 400    # m, width of that preference
MATCH_LEVELS   <- c(5, 10, 15, 20, 30)  # Mg C/ha, ACD match ladder
MATCH_YEAR_TOL <- 30     # Mg C/ha, max difference in any single annual map
MATCH_SCALE    <- 30     # Mg C/ha, scaling of the match score
BROTHER_TRIES  <- 60     # replacement candidates tried, best first, until one
                         #   of them admits a valid pair
MAX_PAIR_EXACT  <- 400   # annulus candidates that get the exact footprint when
                         #   pairing a replacement: the exact extraction
                         #   dominates the runtime and pixel values are enough
                         #   to rule out the hopeless ones
MAX_TRIES_TOTAL <- 25    # hard cap on pair searches per replacement, so one
                         #   awkward plot cannot stall the whole script
MAX_PER_BAND   <- 200    # corridor candidates per road band that get the exact
                         #   footprint treatment. That extraction is the slow
                         #   step and only the best ~60 are ever used, so the
                         #   field is ranked on cheap pixel values first.
# --- v7: terrain and spread for the brothers --------------------------------
SLOPE_LEVELS       <- c(15, 20, 25)  # deg, footprint-mean slope ladder
CLIMB_ROAD_MAX     <- 60   # m, |height above the nearest north-road point|
APPROACH_SLOPE_MAX <- 25   # deg, steepest ground on the straight line from that
                           #   road point to the plot
CLIMB_PAIR_MAX     <- 40   # m, |height difference| brother -> its pair
NODATA_BUFFER      <- 90   # m, keep every plot this far from any cell that is
                           #   nodata in ANY year: the ACD maps carry no data
                           #   over open water, so this is the river mask, and
                           #   it also drops the artefact-prone banks beside it
N_BANDS            <- 10   # road split into this many northing bands; one
                           #   brother per band is the target
W_MATCH <- 0.30; W_NROAD <- 0.15; W_SLOPE <- 0.15; W_CLIMB <- 0.10
W_NORTH <- 0.10; W_SPREAD <- 0.10; W_BHOMOG <- 0.05; W_BSTAB <- 0.05
N_BINS         <- 10     # strata = deciles of the project ACD
TAIL_BOOST     <- 2      # extreme strata are (1+TAIL_BOOST)x as valuable as the
                         #   central ones
PLOT_ID_FIELD  <- "Pl"
EPSG           <- 32650

W_CONTRAST <- 0.30; W_COVER <- 0.30; W_ACCESS <- 0.15; W_HOMOG <- 0.10; W_PROX <- 0.15

# with embeddings (MODE = "gee") the GEE script already applies:
# 0.40 contrast / 0.20 embedding / 0.15 access / 0.10 homogeneity / 0.15 proximity

# footprints must not touch
stopifnot(MIN_SEP > 2 * PLOT_RADIUS)

# =============================================================================
#  shared helper: greedy assignment + refinement, most-constrained plot first
# =============================================================================
#  cand: data.frame with Pl_ref, x, y, d_pair, score, and everything to carry
#        through. Rows must already satisfy every per-plot constraint,
#        including the nearest-neighbour rule against the EXISTING plots.
#  The rule enforced here is the one between the NEW plots: no new plot may end
#  up closer to another plot than to its own reference plot.
greedy_assign <- function(cand, counts0, min_sep = MIN_SEP,
                          nn_margin = NN_MARGIN, rounds = ROUNDS) {
  by_ref <- split(cand, cand$Pl_ref)
  # deterministic: fewest options first, ties broken by reference plot number
  by_ref <- by_ref[order(sapply(by_ref, nrow),
                         as.numeric(names(by_ref)))]

  # how much a candidate in each stratum helps the tail-weighted spread
  cover <- function(bins, counts) {
    deficit <- pmax(TARGET - counts, 0) * U
    if (max(deficit) <= 0) return(rep(0, length(bins)))
    deficit[bins] / max(deficit)
  }
  full_score <- function(d, counts) d$base + W_COVER * cover(d$bin, counts)

  pick_best <- function(d, taken, counts) {
    if (nrow(taken) > 0) {
      ok <- rep(TRUE, nrow(d))
      for (k in seq_len(nrow(taken))) {
        dk <- sqrt((d$x - taken$x[k])^2 + (d$y - taken$y[k])^2)
        ok <- ok & dk >= min_sep &
              dk >= pmax(d$d_pair, taken$d_pair[k]) + nn_margin
      }
      d <- d[ok, , drop = FALSE]
    }
    if (nrow(d) == 0) return(NULL)
    d$score <- full_score(d, counts)
    d[which.max(d$score), , drop = FALSE]
  }

  counts <- counts0
  picked <- list()
  for (nm in names(by_ref)) {
    taken <- if (length(picked)) do.call(rbind, picked)[, c("x", "y", "d_pair")]
             else data.frame(x = numeric(0), y = numeric(0), d_pair = numeric(0))
    best <- pick_best(by_ref[[nm]], taken, counts)
    if (is.null(best)) { message("no candidate left for reference plot ", nm); next }
    picked[[nm]] <- best
    counts[best$bin] <- counts[best$bin] + 1
  }

  # the greedy order is arbitrary: re-pick each plot against the others, with
  # its own stratum taken back out of the counts, until nothing moves
  for (it in seq_len(rounds)) {
    improved <- 0
    worst_first <- names(picked)[order(sapply(picked, function(p) p$score),
                                       as.numeric(names(picked)))]
    for (nm in worst_first) {
      cur <- picked[[nm]]
      others <- picked[names(picked) != nm]
      taken <- if (length(others)) do.call(rbind, others)[, c("x", "y", "d_pair")]
               else data.frame(x = numeric(0), y = numeric(0), d_pair = numeric(0))
      counts_wo <- counts; counts_wo[cur$bin] <- counts_wo[cur$bin] - 1
      best <- pick_best(by_ref[[nm]], taken, counts_wo)
      if (is.null(best)) next
      cur_score <- full_score(cur, counts_wo)      # fair comparison
      if (best$score > cur_score + TOL) {
        picked[[nm]] <- best
        counts <- counts_wo; counts[best$bin] <- counts[best$bin] + 1
        improved <- improved + 1
      } else {
        counts <- counts_wo; counts[cur$bin] <- counts[cur$bin] + 1
      }
    }
    message("refinement round ", it, ": ", improved, " plots moved")
    if (improved == 0) break
  }
  message("final stratum counts: ", paste(counts, collapse = " "),
          "  | target: ", paste(round(TARGET, 1), collapse = " "))
  do.call(rbind, picked)
}

write_outputs <- function(new_pts, plots, stem) {
  g <- st_as_sf(new_pts, coords = c("x", "y"), crs = EPSG, remove = FALSE)
  ll <- st_transform(g, 4326)
  g$lon <- round(st_coordinates(ll)[, 1], 6)
  g$lat <- round(st_coordinates(ll)[, 2], 6)
  st_write(g, file.path(DIR_OUT, paste0(stem, ".gpkg")), layer = "paired_plots",
           delete_dsn = TRUE, quiet = TRUE)
  st_write(plots, file.path(DIR_OUT, paste0(stem, ".gpkg")),
           layer = "existing_plots", append = TRUE, quiet = TRUE)
  st_write(st_buffer(g, PLOT_RADIUS), file.path(DIR_OUT, paste0(stem, ".gpkg")),
           layer = paste0("new_plot_footprints_", PLOT_RADIUS, "m"),
           append = TRUE, quiet = TRUE)
  kml <- file.path(DIR_OUT, paste0(stem, ".kml"))
  st_write(ll, kml, delete_dsn = file.exists(kml), quiet = TRUE)
  write.csv(st_drop_geometry(g), file.path(DIR_OUT, paste0(stem, ".csv")),
            row.names = FALSE)
  # The links layer is a convenience for QGIS, not a deliverable, so a failure
  # here must not cost the run: the points, footprints, KML and CSV are already
  # on disk by this line.
  tryCatch({
    # pair links, useful for a quick visual check in QGIS. Each row is joined to
    # the plot it is paired WITH, which is not always an existing plot: a
    # replacement (role "road_brother") links to its own pair, and that pair links
    # back to the replacement. Everything below works on plain numeric vectors,
    # because `g` is an sf object and its geometry column is sticky.
    gx   <- as.numeric(g$x);  gy <- as.numeric(g$y)
    pid  <- as.character(as.integer(g$Pl_ref))
    pnew <- as.character(g$Pl_new)
    role <- if (is.null(g$role)) rep("pair", nrow(g)) else as.character(g$role)

    ex_id <- as.character(as.integer(plots[[PLOT_ID_FIELD]]))
    ex_xy <- st_coordinates(plots)
    ax <- rep(NA_real_, nrow(g)); ay <- rep(NA_real_, nrow(g))

    i <- which(role == "pair")                    # -> its existing plot
    if (length(i)) {
      k <- match(pid[i], ex_id)
      ax[i] <- ex_xy[k, 1]; ay[i] <- ex_xy[k, 2]
    }
    i <- which(role == "brother_pair")            # -> its replacement plot
    if (length(i)) {
      k <- match(pid[i], pnew)
      ax[i] <- gx[k]; ay[i] <- gy[k]
    }
    i <- which(role == "road_brother")            # -> its own pair
    if (length(i)) {
      k <- match(paste0(pid[i], "B"), pnew)
      ax[i] <- gx[k]; ay[i] <- gy[k]
    }

    ok <- !is.na(ax) & !is.na(ay)
    if (any(!ok))
      message(sum(!ok), " plot(s) had no partner to draw a link to: ",
              paste(pnew[!ok], collapse = ", "))
    if (any(ok)) {
      links <- st_sfc(lapply(which(ok), function(j)
        st_linestring(matrix(c(ax[j], gx[j], ay[j], gy[j]), ncol = 2))), crs = EPSG)
      st_write(st_sf(Pl_new = pnew[ok], Pl_ref = g$Pl_ref[ok], role = role[ok],
                     geometry = links),
               file.path(DIR_OUT, paste0(stem, ".gpkg")), layer = "pair_links",
               append = TRUE, quiet = TRUE)
    }
  }, error = function(e)
    message("pair_links layer skipped: ", conditionMessage(e)))

  message("wrote ", nrow(g), " paired plots to ", DIR_OUT)
  invisible(g)
}

# =============================================================================
#  MODE = "local"
# =============================================================================
if (MODE == "local") {

  acd   <- rast(F_ACD)[[1]]; names(acd) <- "acd"
  plots <- st_read(F_PLOTS, quiet = TRUE) |> st_transform(EPSG)
  aoi   <- st_read(F_AOI,   quiet = TRUE) |> st_transform(EPSG)
  roads <- st_read(F_ROADS, quiet = TRUE) |> st_transform(EPSG)

  # --- clip ACD to the project area -----------------------------------------
  acd <- mask(acd, vect(aoi))

  # --- v4: five annual maps, footprint mean per year, median across years ----
  # A 60 m plot spans about four 30 m pixels, so the value a field crew records
  # is the footprint mean, not the centre pixel. And a single year's map carries
  # the error of one Landsat / SAR composite, so selection works from the median
  # across 2021-2025 and rejects ground that is not stable.
  YEARS  <- as.numeric(names(F_ACD_YEARS))
  w_plot <- focalMat(acd, PLOT_RADIUS, "circle")
  fp_one <- function(f) {
    r <- terra::mask(rast(f)[[1]], vect(aoi))
    terra::focal(r, w = w_plot, fun = "sum", na.rm = TRUE) /
      terra::focal(!is.na(r), w = w_plot, fun = "sum", na.rm = TRUE)
  }
  FP <- rast(lapply(F_ACD_YEARS, fp_one)); names(FP) <- paste0("y", YEARS)

  acd_plot <- app(FP, median, na.rm = TRUE); names(acd_plot) <- "acd"
  mad_t    <- 1.4826 * app(abs(FP - acd_plot), median, na.rm = TRUE)
  rs_t     <- mad_t / max(acd_plot, RS_FLOOR); names(rs_t) <- "rs"
  yc       <- YEARS - mean(YEARS)
  slope_t  <- app(FP - acd_plot, function(v) sum(v * yc) / sum(yc^2))
  names(slope_t) <- "slope"

  # NOTE: rs_t / slope_t / stable are pixel-centred and are used only for the
  # map layers and the project-wide statistics below. Every plot and candidate
  # value is recomputed exactly from its own footprint (foot_exact) further down.
  # --- v8: open water and its banks ----------------------------------------
  # The annual ACD maps are nodata over water, so the river network is already
  # in them. Anything within NODATA_BUFFER of a nodata cell in any year is out.
  nodata_any <- app(FP, function(v) any(is.na(v)))
  d_nodata   <- distance(ifel(nodata_any, 1, NA))
  dry        <- d_nodata >= NODATA_BUFFER
  in_aoi <- !is.na(acd_plot)          # the raster is wider than the project
  message(sprintf("excluded as water or water edge: %.1f%% of the project area",
                  100 * (1 - global(mask(dry, in_aoi, maskvalues = c(FALSE, NA)),
                                    "mean", na.rm = TRUE)[[1]])))
  acd_plot <- mask(acd_plot, dry, maskvalues = c(FALSE, NA))

  stable <- (rs_t <= RS_MAX) & (abs(slope_t) <= SLOPE_MAX) & dry
  names(stable) <- "stable"
  message(sprintf("temporally stable: %.0f%% of the project area",
                  100 * global(mask(stable, in_aoi, maskvalues = c(FALSE, NA)),
                               "mean", na.rm = TRUE)[[1]]))

  # --- ACD heterogeneity over the 90 m (3 x 3) neighbourhood ----------------
  # 90 m covers the 60 m footprint plus a 15 m margin, so this now tests
  # whether the plot itself straddles an edge
  acd_sd <- focal(acd_plot, w = 3, fun = sd, na.rm = TRUE, na.policy = "omit")
  names(acd_sd) <- "acd_sd"

  # --- distance to the nearest road -----------------------------------------
  road_r <- rasterize(vect(roads), rast(acd), field = 1, touches = TRUE)
  road_d <- distance(road_r)            # distance to the nearest non-NA cell
  names(road_d) <- "road_d"

  # --- distance to the north road only: the corridor the crews work from ----
  nroad  <- st_read(F_NROAD, quiet = TRUE)
  if (is.na(st_crs(nroad))) nroad <- st_set_crs(nroad, NROAD_EPSG)
  nroad  <- st_transform(nroad, EPSG)
  nroad_r <- rasterize(vect(nroad), rast(acd), field = 1, touches = TRUE)
  nroad_d <- distance(nroad_r); names(nroad_d) <- "nroad_d"

  # --- v7: terrain ----------------------------------------------------------
  DEM   <- project(rast(F_DEM), rast(acd), method = "bilinear"); names(DEM) <- "elev"
  SLOPE <- terrain(DEM, v = "slope", unit = "degrees"); names(SLOPE) <- "slope_deg"
  # elevation of the nearest point on the north road, and the height above it
  nroad_elev <- mask(DEM, nroad_r)
  NR_ELEV <- focal(nroad_elev, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
  for (k in 1:40)   # spread the road elevation outwards to fill the corridor
    NR_ELEV <- focal(NR_ELEV, w = 3, fun = "mean", na.rm = TRUE, na.policy = "only")
  CLIMB_ROAD <- DEM - NR_ELEV; names(CLIMB_ROAD) <- "climb_road"

  # exact footprint mean of a single raster, and the steepest ground on the
  # straight approach from the road
  foot_one <- function(xy, r) {
    b <- vect(st_buffer(st_as_sf(as.data.frame(xy), coords = 1:2, crs = EPSG),
                        PLOT_RADIUS))
    terra::extract(r, b, fun = "mean", exact = TRUE, ID = FALSE)[, 1]
  }
  nroad_u <- st_union(nroad)          # built once, not once per candidate
  approach_max_slope <- function(px, py) {
    nr <- st_nearest_points(st_sfc(st_point(c(px, py)), crs = EPSG), nroad_u)
    rp <- st_coordinates(nr)[1, 1:2]
    n  <- max(2, ceiling(sqrt(sum((c(px, py) - rp)^2)) / res(acd)[1]) + 1)
    t  <- seq(0, 1, length.out = n)
    max(terra::extract(SLOPE, cbind(rp[1] + t * (px - rp[1]),
                                    rp[2] + t * (py - rp[2])))[, 1], na.rm = TRUE)
  }
  roads_all <- st_union(st_union(roads), st_union(nroad))
  road_dist_exact <- function(px, py)
    as.numeric(st_distance(st_sfc(st_point(c(px, py)), crs = EPSG), roads_all))

  # --- project-wide ACD distribution ----------------------------------------
  # per-pixel, because this describes the area; plot values are footprint means
  # and so sit slightly closer to the centre of this distribution
  v_all <- values(acd_plot, mat = FALSE); v_all <- v_all[!is.na(v_all)]
  qs  <- quantile(v_all, c(1/3, 0.5, 2/3))
  t33 <- qs[[1]]; med <- qs[[2]]; t67 <- qs[[3]]
  message(sprintf("ACD terciles: <%.1f | %.1f-%.1f | >%.1f Mg C/ha",
                  t33, t33, t67, t67))
  acd_class <- function(v) ifelse(v < t33, "low", ifelse(v > t67, "high", "mid"))

  # --- v9: the hard plots are retired ---------------------------------------
  # They stay only as the ACD targets their north-road replacements must match:
  # they leave the design, stop constraining it, and are not written out.
  # Defined here because the strata targets below already need it.
  plots$active <- !(plots[[PLOT_ID_FIELD]] %in% HARD_PLOTS)
  ACT <- which(plots$active)
  message(sprintf("retired as hard to reach: %s -> %d existing plots kept, %d replaced",
                  paste(sort(HARD_PLOTS), collapse = ", "), length(ACT),
                  length(HARD_PLOTS)))

  # --- strata for the spread objective --------------------------------------
  # Deciles of the footprint-scale ACD, so each stratum holds 10% of the project
  # area. Tail strata carry more weight: that is where the biomass-to-remote-
  # sensing relationship is least constrained by the existing plots.
  v_fp  <- values(acd_plot, mat = FALSE); v_fp <- v_fp[!is.na(v_fp)]
  EDGES <- quantile(v_fp, seq(0, 1, length.out = N_BINS + 1))
  EDGES[1] <- -Inf; EDGES[length(EDGES)] <- Inf
  acd_bin <- function(v) as.integer(cut(v, EDGES, labels = FALSE,
                                        include.lowest = TRUE))
  U      <<- 1 + TAIL_BOOST * abs(seq_len(N_BINS) - (N_BINS + 1) / 2) /
                              ((N_BINS - 1) / 2)
  TARGET <<- 2 * (length(ACT) + length(HARD_PLOTS)) * U / sum(U)
  message("decile edges: ", paste(round(EDGES[2:N_BINS], 1), collapse = " "))
  message("target counts: ", paste(round(TARGET, 1), collapse = " "))

  # --- reference plot properties --------------------------------------------
  pxy <- st_coordinates(plots)
  plots$road_ref <- as.numeric(st_distance(plots, st_union(roads)))

  # --- v5: EXACT footprint values -------------------------------------------
  # The value of a plot is the area-weighted mean of each annual map over the
  # circle of radius PLOT_RADIUS centred on the PLOT POINT itself. terra's
  # exact=TRUE weights every cell by the fraction of it the circle covers, so
  # this is the true 60 m-diameter footprint and not a pixel-centred stand-in.
  foot_exact <- function(xy) {
    b <- vect(st_buffer(st_as_sf(as.data.frame(xy), coords = 1:2, crs = EPSG),
                        PLOT_RADIUS))
    as.matrix(terra::extract(FP, b, fun = "mean", exact = TRUE, ID = FALSE))
  }
  foot_stats <- function(vy) {                       # vy: n x n_years
    med_ <- apply(vy, 1, median, na.rm = TRUE)
    mad_ <- 1.4826 * apply(abs(vy - med_), 1, median, na.rm = TRUE)
    list(med = med_, rs = mad_ / pmax(med_, RS_FLOOR),
         slope = as.numeric((vy - med_) %*% yc) / sum(yc^2))
  }

  # --- is each EXISTING plot itself stable? ---------------------------------
  REF_Y  <- foot_exact(pxy)                          # n_plots x n_years
  ref_st <- foot_stats(REF_Y)
  plots$ACD_ref   <- ref_st$med                      # exact footprint median
  plots$rs_ref    <- ref_st$rs
  plots$slope_ref <- ref_st$slope
  plots$target    <- ifelse(plots$ACD_ref > med, "low", "high")

  # --- v9: the hard plots are retired ---------------------------------------
  # They stay only as the ACD targets their north-road replacements must match:
  # they leave the design, stop constraining it, and are not written out.
  PXY_A <- pxy[ACT, , drop = FALSE]
  plots$stable_ref <- plots$rs_ref <= RS_MAX & abs(plots$slope_ref) <= SLOPE_MAX
  rep <- data.frame(Pl = plots[[PLOT_ID_FIELD]], ACD = round(plots$ACD_ref, 1),
                    rs = round(plots$rs_ref, 3), slope = round(plots$slope_ref, 2),
                    stable = plots$stable_ref, round(REF_Y, 1))
  names(rep)[6:(5 + length(YEARS))] <- as.character(YEARS)
  message("--- temporal stability of the existing plots ---")
  print(rep[order(-rep$rs), ], row.names = FALSE)
  message(sum(!plots$stable_ref), " of ", nrow(plots),
          " existing plots are NOT temporally stable")

  # --- candidate offsets: a 10 m grid inside the 100-300 m annulus ----------
  off <- seq(-PAIR_MAX, PAIR_MAX, by = GRID)
  gr  <- expand.grid(dx = off, dy = off)
  gr$d <- sqrt(gr$dx^2 + gr$dy^2)
  gr  <- gr[gr$d >= PAIR_MIN & gr$d <= PAIR_MAX, ]

  cand_list <- vector("list", nrow(plots))
  for (i in seq_len(nrow(plots))) {
    if (!plots$active[i]) next
    cx <- pxy[i, 1] + gr$dx; cy <- pxy[i, 2] + gr$dy
    xy <- cbind(cx, cy)
    # cheap cell-level filters first; the exact footprint extraction is the
    # expensive step, so only survivors reach it
    v  <- terra::extract(c(acd_plot, acd_sd, road_d), xy)
    road_cap <- max(ROAD_CAP_FLOOR, plots$road_ref[i])
    keep <- !is.na(v$acd) & !is.na(v$acd_sd) &
            v$acd_sd <= SD_HARD &
            v$road_d >= ROAD_MIN & v$road_d <= road_cap
    if (!any(keep)) { cand_list[[i]] <- NULL; next }
    xy <- xy[keep, , drop = FALSE]; v <- v[keep, ]; dpk <- gr$d[keep]
    # at least MIN_SEP from every existing plot centre, and - the
    # nearest-neighbour rule - closer to its own reference plot than to any
    # other existing plot
    dall <- apply(xy, 1, function(p)
      sqrt((pxy[ACT, 1] - p[1])^2 + (pxy[ACT, 2] - p[2])^2))    # nkept x ncand
    dmin     <- apply(dall, 2, min)
    dmin_oth <- apply(dall[which(ACT != i), , drop = FALSE], 2, min)
    geo <- dmin >= MIN_SEP & dmin_oth >= dpk + NN_MARGIN
    if (!any(geo)) { cand_list[[i]] <- NULL; next }
    xy <- xy[geo, , drop = FALSE]; v <- v[geo, ]; dpk <- dpk[geo]
    dmin_oth <- dmin_oth[geo]

    # exact 60 m footprint, every year, centred on each candidate point
    vy  <- foot_exact(xy)
    fst <- foot_stats(vy)

    # the contrast must hold in EVERY annual map: same sign, and at least
    # DELTA_YEAR_MIN Mg C/ha each year
    dy        <- sweep(vy, 2, REF_Y[i, ], "-")       # n_cand x n_years
    same_sign <- apply(dy > 0, 1, all) | apply(dy < 0, 1, all)
    dy_min    <- apply(abs(dy), 1, min)
    ok <- !is.na(fst$med) & fst$med >= MIN_ACD &
          fst$rs <= RS_MAX & abs(fst$slope) <= SLOPE_MAX &
          same_sign & dy_min >= DELTA_YEAR_MIN
    if (!any(ok)) { cand_list[[i]] <- NULL; next }
    ref <- plots$ACD_ref[i]
    dacd <- if (DIRECTION == "free") abs(fst$med[ok] - ref) else
            if (plots$target[i] == "low") ref - fst$med[ok] else fst$med[ok] - ref
    cand_list[[i]] <- data.frame(
      Pl_ref = plots[[PLOT_ID_FIELD]][i], target = plots$target[i],
      ACD_ref = ref, ACD_new = fst$med[ok], dACD = dacd,
      acd_sd = v$acd_sd[ok], road_new = v$road_d[ok],
      road_ref = plots$road_ref[i], d_pair = dpk[ok], d_other = dmin_oth[ok],
      dACD_min_year = dy_min[ok],
      rs_new = fst$rs[ok], slope_new = fst$slope[ok],
      rs_ref = plots$rs_ref[i], slope_ref = plots$slope_ref[i],
      ref_stable = plots$stable_ref[i],
      x = xy[ok, 1], y = xy[ok, 2])
  }
  cand <- do.call(rbind, cand_list)
  cand$n_cand <- ave(cand$dACD, cand$Pl_ref, FUN = length)

  # --- keep, per reference plot, the strictest contrast level that is met ---
  cand <- do.call(rbind, lapply(split(cand, cand$Pl_ref), function(d) {
    lev <- DELTA_LEVELS[which(sapply(DELTA_LEVELS, function(l) any(d$dACD >= l)))[1]]
    d <- d[d$dACD >= lev, ]; d$dACD_req <- lev; d
  }))

  # --- score ----------------------------------------------------------------
  # everything except the spread term, which depends on what the other plots
  # have already taken and so is computed inside greedy_assign()
  cl <- function(v, lo = 0, hi = 1) pmin(pmax(v, lo), hi)
  cand$bin  <- acd_bin(cand$ACD_new)
  cand$base <- W_CONTRAST * cl(cand$dACD / CONTRAST_SAT) +
               W_ACCESS   * (1 - cl(cand$road_new / 400)) +
               W_HOMOG    * (1 - cl(cand$acd_sd / SD_SCALE)) +
               W_PROX     * (1 - cl((cand$d_pair - PAIR_MIN) / (PAIR_MAX - PAIR_MIN)))

  # --- greedy assignment ----------------------------------------------------
  # the 46 existing plots already occupy strata; the new ones fill the gaps
  n_ref <- tabulate(acd_bin(plots$ACD_ref[ACT]), nbins = N_BINS)
  message("existing plots per stratum: ", paste(n_ref, collapse = " "))
  new <- greedy_assign(cand, n_ref)
  new$target <- ifelse(new$ACD_new > new$ACD_ref, "high", "low")
  new$cls_ref <- acd_class(new$ACD_ref)
  new$cls_new <- acd_class(new$ACD_new)
  new$Pl_new  <- paste0(new$Pl_ref, "B")
  new <- new[order(new$Pl_ref), ]

  # road_d is a 30 m raster, so re-measure the chosen points exactly: with a
  # 30 m radius a half-cell error is enough to put the road inside the plot
  new$road_exact <- as.numeric(st_distance(
    st_as_sf(new, coords = c("x", "y"), crs = EPSG, remove = FALSE),
    st_union(roads)))
  bad_road <- new$road_exact < ROAD_MIN
  if (any(bad_road))
    warning(sum(bad_road), " plot(s) within ", ROAD_MIN, " m of a road: ",
            paste(sprintf("%s (%.0f m)", new$Pl_new[bad_road],
                          new$road_exact[bad_road]), collapse = ", "),
            " - move these by hand or raise ROAD_MIN and re-run")

  print(new[, c("Pl_new", "target", "cls_ref", "cls_new", "ACD_ref", "ACD_new",
                "dACD", "dACD_min_year", "bin", "d_pair", "d_other",
                "road_new", "road_exact", "rs_new", "ref_stable")],
        digits = 4)
  cat(sprintf("\n|dACD| median %.1f, min %.1f Mg C/ha\n",
              median(new$dACD), min(new$dACD)))

  # ==========================================================================
  #  v6: brothers along the north road for the hard-to-reach plots
  # ==========================================================================
  # Each hard plot keeps its own pair. In addition it gets a BROTHER on the
  # north road whose multi-year ACD matches it as closely as possible, and that
  # brother gets its own contrasting pair under exactly the same rules.

  # best contrasting pair for an arbitrary centre point, same rules as above
  pair_for_point <- function(bx, by, ref_years_b, ref_acd_b, occ, occ_dp, counts_now,
                             slope_cap = NULL, elev_b = NULL) {
    cx <- bx + gr$dx; cy <- by + gr$dy
    xy <- cbind(cx, cy)
    v  <- terra::extract(c(acd_plot, acd_sd, road_d), xy)
    keep <- !is.na(v$acd) & !is.na(v$acd_sd) & v$acd_sd <= SD_HARD &
            v$road_d >= ROAD_MIN & v$road_d <= max(ROAD_CAP_FLOOR, 400)
    if (!any(keep)) return(NULL)
    xy <- xy[keep, , drop = FALSE]; v <- v[keep, ]; dpk <- gr$d[keep]

    # cheap pixel-level contrast screen before any exact footprint work: a
    # candidate whose pixel ACD is nowhere near the required contrast will not
    # reach it once measured exactly either. 20 Mg C/ha of slack covers the
    # difference between the pixel value and the footprint mean.
    dpix <- abs(v$acd - ref_acd_b)
    k2 <- dpix >= (min(DELTA_LEVELS[DELTA_LEVELS > 0]) - 20)
    if (!any(k2)) return(NULL)
    if (sum(k2) > MAX_PAIR_EXACT) {
      j2 <- which(k2)
      j2 <- j2[order(-dpix[j2], xy[j2, 1], xy[j2, 2])][seq_len(MAX_PAIR_EXACT)]
      k2 <- rep(FALSE, length(k2)); k2[sort(j2)] <- TRUE
    }
    xy <- xy[k2, , drop = FALSE]; v <- v[k2, ]; dpk <- dpk[k2]

    # separation and the nearest-neighbour rule against everything placed
    ok <- rep(TRUE, nrow(xy))
    for (k in seq_len(nrow(occ))) {
      dk <- sqrt((xy[, 1] - occ[k, 1])^2 + (xy[, 2] - occ[k, 2])^2)
      ok <- ok & dk >= MIN_SEP & dk >= pmax(dpk, occ_dp[k], na.rm = TRUE) + NN_MARGIN
    }
    if (!any(ok)) return(NULL)
    xy <- xy[ok, , drop = FALSE]; v <- v[ok, ]; dpk <- dpk[ok]
    if (!is.null(slope_cap)) {                       # v7: terrain for the pair
      tsl <- foot_one(xy, SLOPE)
      ter <- terra::extract(c(DEM, CLIMB_ROAD), xy)
      tk  <- !is.na(tsl) & tsl <= slope_cap & abs(ter$climb_road) <= CLIMB_ROAD_MAX
      if (!is.null(elev_b)) tk <- tk & abs(ter$elev - elev_b) <= CLIMB_PAIR_MAX
      if (!any(tk)) return(NULL)
      xy <- xy[tk, , drop = FALSE]; v <- v[tk, ]; dpk <- dpk[tk]
      tsl <- tsl[tk]; ter <- ter[tk, ]
    }
    vy  <- foot_exact(xy); fst <- foot_stats(vy)
    dy  <- sweep(vy, 2, ref_years_b, "-")
    good <- !is.na(fst$med) & fst$med >= MIN_ACD &
            fst$rs <= RS_MAX & abs(fst$slope) <= SLOPE_MAX &
            (apply(dy > 0, 1, all) | apply(dy < 0, 1, all)) &
            apply(abs(dy), 1, min) >= DELTA_YEAR_MIN
    if (!any(good)) return(NULL)
    dabs <- abs(fst$med - ref_acd_b)
    for (delta in DELTA_LEVELS) {
      m <- good & dabs >= delta
      if (!any(m)) next
      cl <- function(z, lo = 0, hi = 1) pmin(pmax(z, lo), hi)
      bins <- acd_bin(fst$med[m])
      deficit <- pmax(TARGET - counts_now, 0) * U
      cov <- if (max(deficit) > 0) deficit[bins] / max(deficit) else rep(0, sum(m))
      sc <- W_CONTRAST * cl(dabs[m] / CONTRAST_SAT) + W_COVER * cov +
            W_ACCESS * (1 - cl(v$road_d[m] / 400)) +
            W_HOMOG  * (1 - cl(v$acd_sd[m] / SD_SCALE)) +
            W_PROX   * (1 - cl((dpk[m] - PAIR_MIN) / (PAIR_MAX - PAIR_MIN)))
      j <- which(m)[which.max(sc)]
      if (!is.null(slope_cap)) {
        if (approach_max_slope(xy[j, 1], xy[j, 2]) > APPROACH_SLOPE_MAX) next
        if (road_dist_exact(xy[j, 1], xy[j, 2]) < ROAD_MIN) next
      }
      return(list(x = xy[j, 1], y = xy[j, 2], acd = fst$med[j], yr = vy[j, ],
                  terr_slope = if (is.null(slope_cap)) NA_real_ else tsl[j],
                  elev = if (is.null(slope_cap)) NA_real_ else ter$elev[j],
                  climb_road = if (is.null(slope_cap)) NA_real_ else ter$climb_road[j],
                  rs = fst$rs[j], acd_slope = fst$slope[j], sd3 = v$acd_sd[j],
                  road = v$road_d[j], dpair = dpk[j],
                  dy_min = min(abs(dy[j, ])), bin = acd_bin(fst$med[j]),
                  delta_req = delta, score = max(sc), d = dabs[j]))
    }
    NULL
  }

  # corridor cells: in the project area, within NROAD_MAX of the north road,
  # clear of any road, homogeneous enough to hold a plot
  corr <- !is.na(acd_plot) & dry & stable & nroad_d <= NROAD_MAX & road_d >= ROAD_MIN &
          !is.na(acd_sd) & acd_sd <= SD_HARD &
          !is.na(CLIMB_ROAD) & abs(CLIMB_ROAD) <= CLIMB_ROAD_MAX
  corr_xy <- xyFromCell(corr, which(values(corr, mat = FALSE) %in% TRUE))
  message(sprintf("north-road corridor: %d candidate cells (%.0f ha)",
                  nrow(corr_xy), nrow(corr_xy) * prod(res(acd)) / 1e4))
  corr_v <- terra::extract(c(acd_plot, acd_sd, road_d, nroad_d, DEM, CLIMB_ROAD),
                           corr_xy)
  # northing bands along the road: one brother per band is the target
  NY0 <- st_bbox(nroad)[["ymin"]]; NY1 <- st_bbox(nroad)[["ymax"]]
  road_band <- function(y) pmin(pmax(floor((y - NY0) / (NY1 - NY0) * N_BANDS), 0),
                               N_BANDS - 1) + 1
  band_used <- integer(N_BANDS)

  occ    <- rbind(PXY_A, as.matrix(new[, c("x", "y")]))
  occ_dp <- c(rep(NA_real_, nrow(PXY_A)), new$d_pair)
  counts_now <- n_ref + tabulate(new$bin, nbins = N_BINS)
  extra <- list()
  # replacements continue the plot numbering: 47, 48, ... in hard-plot order
  NEXT_ID <- max(plots[[PLOT_ID_FIELD]]) + 1
  brother_id <- setNames(NEXT_ID + seq_along(sort(HARD_PLOTS)) - 1, sort(HARD_PLOTS))

  for (h in sort(HARD_PLOTS)) {
    i <- which(plots[[PLOT_ID_FIELD]] == h)
    ref_acd <- plots$ACD_ref[i]; ref_yrs <- REF_Y[i, ]
    placed <- NULL
    tries  <- 0
    # The exact footprint work is the expensive step, so do it ONCE per plot on
    # the widest candidate set the ladder could ever reach, then subset.
    pre <- which(!is.na(corr_v$acd) &
                 abs(corr_v$acd - ref_acd) <= max(MATCH_LEVELS) + 10)
    if (!length(pre)) { message("!! no corridor candidate for hard plot ", h); next }
    # rank on cheap pixel-level values and keep the best MAX_PER_BAND in each
    # road band, so every band still gets a fair chance at the exact stage
    cl0 <- function(z, lo = 0, hi = 1) pmin(pmax(z, lo), hi)
    sl_px  <- terra::extract(SLOPE, corr_xy[pre, , drop = FALSE])[, 1]
    pre_sc <- W_MATCH * (1 - cl0(abs(corr_v$acd[pre] - ref_acd) / MATCH_SCALE)) +
              W_NROAD * exp(-((corr_v$nroad_d[pre] - NROAD_PREF) / NROAD_PREF_SD)^2) +
              W_SLOPE * (1 - cl0(sl_px / max(SLOPE_LEVELS))) +
              W_NORTH * cl0((corr_xy[pre, 2] - NY0) / (NY1 - NY0))
    pre_bd <- road_band(corr_xy[pre, 2])
    keep_i <- unlist(lapply(sort(unique(pre_bd)), function(bd) {
      j <- which(pre_bd == bd)
      j[order(-pre_sc[j], corr_xy[pre[j], 1], corr_xy[pre[j], 2])][
        seq_len(min(MAX_PER_BAND, length(j)))]
    }))
    pre <- pre[sort(keep_i)]

    bxy <- corr_xy[pre, , drop = FALSE]
    by_ <- foot_exact(bxy); bst <- foot_stats(by_)
    bsl <- foot_one(bxy, SLOPE)
    base_ok <- !is.na(bst$med) & bst$med >= MIN_ACD & !is.na(bsl) &
               bst$rs <= RS_MAX & abs(bst$slope) <= SLOPE_MAX &
               apply(abs(sweep(by_, 2, ref_yrs, "-")), 1, max) <= MATCH_YEAR_TOL
    message(sprintf("  plot %d: %d candidates kept for exact checks, %d pass",
                    h, length(pre), sum(base_ok)))
    for (slope_cap in SLOPE_LEVELS) {
     for (tol in MATCH_LEVELS) {
      sel <- which(base_ok & bsl <= slope_cap & abs(bst$med - ref_acd) <= tol)
      if (!length(sel)) next
      cl <- function(z, lo = 0, hi = 1) pmin(pmax(z, lo), hi)
      dnv <- corr_v$nroad_d[pre][sel]; sdv <- corr_v$acd_sd[pre][sel]
      byv <- bxy[sel, 2]; bnd <- road_band(byv)
      sc <- W_MATCH  * (1 - cl(abs(bst$med[sel] - ref_acd) / MATCH_SCALE)) +
            W_NROAD  * exp(-((dnv - NROAD_PREF) / NROAD_PREF_SD)^2) +
            W_SLOPE  * (1 - cl(bsl[sel] / slope_cap)) +
            W_CLIMB  * (1 - cl(abs(corr_v$climb_road[pre][sel]) / CLIMB_ROAD_MAX)) +
            W_NORTH  * cl((byv - NY0) / (NY1 - NY0)) +
            W_SPREAD * as.numeric(band_used[bnd] == 0) +
            W_BHOMOG * (1 - cl(sdv / SD_SCALE)) +
            W_BSTAB  * (1 - cl(bst$rs[sel] / RS_MAX))
      # deterministic: score, then easting, then northing
      ord <- order(-sc, bxy[sel, 1], bxy[sel, 2])
      for (t in head(ord, BROTHER_TRIES)) {
        if (tries >= MAX_TRIES_TOTAL) break
        k <- sel[t]; px <- bxy[k, 1]; py <- bxy[k, 2]
        dmin_occ <- min(sqrt((occ[, 1] - px)^2 + (occ[, 2] - py)^2))
        if (dmin_occ < MIN_SEP) next
        if (approach_max_slope(px, py) > APPROACH_SLOPE_MAX) next
        if (road_dist_exact(px, py) < ROAD_MIN) next
        tries <- tries + 1
        pr <- pair_for_point(px, py, by_[k, ], bst$med[k], occ, occ_dp, counts_now,
                             slope_cap = slope_cap,
                             elev_b = corr_v$elev[pre][sel][t])
        if (is.null(pr)) next
        if (dmin_occ < pr$dpair + NN_MARGIN) next      # brother nearest its pair
        placed <- list(px = px, py = py, k = k, pr = pr, tol = tol,
                       dn = dnv[t], dr = corr_v$road_d[pre][sel][t],
                       sd = sdv[t], med = bst$med[k], rs = bst$rs[k],
                       slope = bst$slope[k], yr = by_[k, ], n = length(sel),
                       tsl = bsl[k], elev = corr_v$elev[pre][sel][t],
                       climb = corr_v$climb_road[pre][sel][t], band = bnd[t])
        break
      }
      if (!is.null(placed) || tries >= MAX_TRIES_TOTAL) break
     }
     if (!is.null(placed) || tries >= MAX_TRIES_TOTAL) break
    }
    if (is.null(placed)) { message("!! no road brother for hard plot ", h); next }
    pr <- placed$pr
    occ    <- rbind(occ, c(placed$px, placed$py), c(pr$x, pr$y))
    occ_dp <- c(occ_dp, pr$dpair, pr$dpair)
    counts_now[acd_bin(placed$med)] <- counts_now[acd_bin(placed$med)] + 1
    counts_now[pr$bin] <- counts_now[pr$bin] + 1
    band_used[placed$band] <- band_used[placed$band] + 1
    message(sprintf("plot %2d -> replacement %d  ACD %6.1f vs %6.1f (tol %2.0f)  %5.0f m from the north road | pair %dB dACD %5.1f",
                    h, brother_id[[as.character(h)]], ref_acd, placed$med,
                    placed$tol, placed$dn, brother_id[[as.character(h)]], pr$d))
    extra[[length(extra) + 1]] <- data.frame(
      Pl_ref = brother_id[[as.character(h)]],
      Pl_new = as.character(brother_id[[as.character(h)]]), replaces = h,
      role = "road_brother",
      ACD_ref = ref_acd, ACD_new = placed$med, dACD = abs(placed$med - ref_acd),
      d_pair = pr$dpair, d_other = dmin_occ, road_new = placed$dr,
      nroad_new = placed$dn, acd_sd = placed$sd, bin = acd_bin(placed$med),
      rs_new = placed$rs, slope_new = placed$slope,
      terr_slope_deg = placed$tsl, elev_m = placed$elev,
      climb_road_m = placed$climb, road_band = placed$band,
      dACD_min_year = NA_real_,
      match_max_year = max(abs(placed$yr - ref_yrs)),
      x = placed$px, y = placed$py, score = NA_real_)
    extra[[length(extra) + 1]] <- data.frame(
      Pl_ref = brother_id[[as.character(h)]],
      Pl_new = paste0(brother_id[[as.character(h)]], "B"), replaces = h,
      role = "brother_pair",
      ACD_ref = placed$med, ACD_new = pr$acd, dACD = pr$d,
      d_pair = pr$dpair, d_other = NA_real_, road_new = pr$road,
      nroad_new = terra::extract(nroad_d, cbind(pr$x, pr$y))[, 1],
      acd_sd = pr$sd3, bin = pr$bin, rs_new = pr$rs, slope_new = pr$acd_slope,
      terr_slope_deg = pr$terr_slope, elev_m = pr$elev,
      climb_road_m = pr$climb_road, road_band = NA_integer_,
      dACD_min_year = pr$dy_min, match_max_year = NA_real_,
      x = pr$x, y = pr$y, score = pr$score)
  }
  brothers <- do.call(rbind, extra)
  new$role <- "pair"; new$nroad_new <- NA_real_; new$match_max_year <- NA_real_
  new$replaces <- NA_real_
  new <- dplyr::bind_rows(new, brothers)
  message("final design: ", sum(plots$active), " existing + ", sum(new$role == "pair"),
          " pairs + ", sum(new$role != "pair"), " road brothers and their pairs")

  plots <- plots[plots$active, ]        # retired plots leave the design
  write_outputs(new, plots, "Kuamut_paired_plots_v9")
}

# =============================================================================
#  MODE = "gee" - final assignment on the exported candidate table
# =============================================================================
if (MODE == "gee") {

  plots <- st_read(F_PLOTS, quiet = TRUE) |> st_transform(EPSG)
  cand  <- read.csv(F_GEECSV)
  stopifnot(all(c("Pl_ref", "lon", "lat", "score") %in% names(cand)))

  xy <- st_coordinates(st_transform(
    st_as_sf(cand, coords = c("lon", "lat"), crs = 4326), EPSG))
  cand$x <- xy[, 1]; cand$y <- xy[, 2]
  cand <- cand[is.finite(cand$score), ]

  # the spread objective needs the ACD raster, which this branch does not read,
  # so the GEE score is used as-is and W_COVER has no effect here
  cand$base <- cand$score; cand$bin <- 1L
  new <- greedy_assign(cand, rep(0L, N_BINS))
  new$Pl_new <- paste0(new$Pl_ref, "B")
  new <- new[order(new$Pl_ref), ]

  print(new[, c("Pl_new", "target", "ACD_ref", "ACD_new", "dACD",
                "emb_dist", "d_pair", "d_other", "road_new", "score")], digits = 4)
  cat(sprintf("\nranks used: %s\n",
              paste(names(table(new$rank)), table(new$rank),
                    sep = "x", collapse = " ")))
  write_outputs(new, plots, "Kuamut_paired_plots_embeddings_v9")
}

# =============================================================================
#  final sanity checks - run on whichever mode produced `new`
# =============================================================================
# Which plot is each new row paired WITH? An ordinary pair points at its
# existing plot, a replacement at its own pair, and that pair back at the
# replacement. Used by the checks below and by the figures.
NEW_ROLE <- if (is.null(new$role)) rep("pair", nrow(new)) else as.character(new$role)
NEW_PID  <- as.character(as.integer(new$Pl_ref))
NEW_NAME <- as.character(new$Pl_new)
EX_ID    <- as.character(as.integer(plots[[PLOT_ID_FIELD]]))
EX_XY    <- st_coordinates(plots)

anchor_xy <- matrix(NA_real_, nrow(new), 2)
i <- which(NEW_ROLE == "pair")
if (length(i)) anchor_xy[i, ] <- EX_XY[match(NEW_PID[i], EX_ID), , drop = FALSE]
i <- which(NEW_ROLE == "brother_pair")
if (length(i)) anchor_xy[i, ] <- cbind(new$x, new$y)[match(NEW_PID[i], NEW_NAME), , drop = FALSE]
i <- which(NEW_ROLE == "road_brother")
if (length(i)) anchor_xy[i, ] <- cbind(new$x, new$y)[match(paste0(NEW_PID[i], "B"),
                                                           NEW_NAME), , drop = FALSE]

all_xy <- rbind(EX_XY, as.matrix(new[, c("x", "y")]))
dm <- as.matrix(dist(all_xy)); diag(dm) <- Inf
cat(sprintf("plots total: %d | closest pair overall: %.1f m (%.1f m between %d m footprints)\n",
            nrow(all_xy), min(dm), min(dm) - 2 * PLOT_RADIUS, PLOT_RADIUS))
cat(sprintf("pair distances: %.0f-%.0f m (target %d-%d)\n",
            min(new$d_pair), max(new$d_pair), PAIR_MIN, PAIR_MAX))
cat(sprintf("road distance of new plots: median %.0f m, max %.0f m\n",
            median(new$road_new), max(new$road_new)))

# nearest-neighbour rule: for every new plot the closest plot must be its pair
lab     <- c(paste0("P", EX_ID), NEW_NAME)
own_lab <- ifelse(NEW_ROLE == "pair", paste0("P", NEW_PID),
           ifelse(NEW_ROLE == "brother_pair", NEW_PID, paste0(NEW_PID, "B")))
own_ix  <- match(own_lab, lab)
dm2 <- dm; viol <- 0
for (k in seq_len(nrow(new))) {
  idx <- nrow(EX_XY) + k
  own <- own_ix[k]
  if (is.na(own)) next
  if (which.min(dm2[idx, ]) != own) {
    viol <- viol + 1
    cat(sprintf("  VIOLATION %s: closest plot is %s (%.0f m) not its pair (%.0f m)\n",
                lab[idx], lab[which.min(dm2[idx, ])], min(dm2[idx, ]), dm2[idx, own]))
  }
}
cat(sprintf("nearest-neighbour rule: %d violations out of %d new plots\n",
            viol, nrow(new)))

## Graphics - Reference information

# =============================================================================
#  Kuamut IFM - results map and diagnostic graphics   (MODE = "local")
#
#  Append this to 02_kuamut_paired_plots_v9.R, after the final sanity checks.
#  It reuses objects built in the "local" block:
#     acd, acd_sd, road_d, aoi, roads, plots, new, t33, med, t67, v_all,
#     acd_class, PLOT_ID_FIELD, DIR_OUT, PAIR_MIN, PAIR_MAX, DELTA_LEVELS
#
#  Outputs in DIR_OUT:
#     Kuamut_paired_plots_v9_map.png          overview map
#     Kuamut_paired_plots_v9_diagnostics.png  8-panel figure
#     Kuamut_paired_plots_v9_pairs_NN.png     per-pair zoom panels
#     Kuamut_paired_plots_v9_summary.csv      the numbers behind the figures
# =============================================================================
if (MODE == "local") {
  
  STEM    <- "Kuamut_paired_plots_v9"
  COL_REF <- "#1B4F72"                                  # existing plots
  COL_NEW <- "#C0392B"                                  # new paired plots
  RAMP    <- hcl.colors(100, "Greens", rev = TRUE)      # light = low ACD
  
  open_png <- function(f, w = 11, h = 9, res = 200)
    png(file.path(DIR_OUT, paste0(STEM, f)),
        width = w, height = h, units = "in", res = res)
  
  ref_xy     <- EX_XY
  new_ref_xy <- anchor_xy          # role-aware: see the sanity-check block
  
  # ---------------------------------------------------------------------------
  #  1. overview map
  # ---------------------------------------------------------------------------
  open_png("_map.png", 11, 9)
  plot(acd, col = RAMP, mar = c(3, 3, 3, 5.5),
       plg = list(title = "ACD\nMg C/ha", title.cex = 0.8, cex = 0.8),
       main = "Kuamut IFM - existing plots and new contrasting pairs")
  plot(st_geometry(aoi),   add = TRUE, col = NA, border = "grey15", lwd = 1.4)
  plot(st_geometry(roads), add = TRUE, col = "grey35", lwd = 0.7)
  segments(new_ref_xy[, 1], new_ref_xy[, 2], new$x, new$y,
           col = "black", lwd = 1.1)
  points(ref_xy, pch = 21, bg = COL_REF, col = "white", cex = 1.1, lwd = 0.6)
  points(new$x, new$y, pch = 24, bg = COL_NEW, col = "white", cex = 1.1, lwd = 0.6)
  text(new$x, new$y, new$Pl_new, pos = 4, cex = 0.45, col = "grey10")
  legend("topleft", bty = "n", cex = 0.8, pt.cex = 1.2,
         pch = c(21, 24, NA, NA), pt.bg = c(COL_REF, COL_NEW, NA, NA),
         lty = c(NA, NA, 1, 1), col = c("white", "white", "black", "grey35"),
         legend = c("existing plot", "new paired plot", "pair link", "road"))
  try(terra::sbar(5000, xy = "bottomleft", type = "bar", divs = 2,
                  below = "m", cex = 0.7), silent = TRUE)
  try(terra::north(xy = "bottomright", type = 2), silent = TRUE)
  dev.off()
  
  # ---------------------------------------------------------------------------
  #  2. numbers behind the graphics
  # ---------------------------------------------------------------------------
  acd_ref_v <- as.numeric(na.omit(plots$ACD_ref))
  acd_new_v <- as.numeric(new$ACD_new)
  acd_all_v <- c(acd_ref_v, acd_new_v)
  vs <- if (length(v_all) > 2e5) sample(v_all, 2e5) else v_all   # for plotting
  
  lev  <- c("low", "mid", "high")
  shr  <- function(x) 100 * as.numeric(table(factor(x, levels = lev))) / length(x)
  area_pct <- 100 * c(mean(v_all <  t33),
                      mean(v_all >= t33 & v_all <= t67),
                      mean(v_all >  t67))
  M <- rbind(`project area`  = area_pct,
             `existing only` = shr(acd_class(acd_ref_v)),
             `new only`      = shr(acd_class(acd_new_v)),
             `existing+new`  = shr(acd_class(acd_all_v)))
  colnames(M) <- lev
  
  # decile occupancy: each decile of the project ACD holds 10% of the area
  dec_br <- quantile(v_all, seq(0, 1, 0.1)); dec_br[1] <- -Inf
  dec_br[length(dec_br)] <- Inf
  dec_of <- function(x) cut(x, dec_br, labels = FALSE, include.lowest = TRUE)
  d_ref  <- as.numeric(table(factor(dec_of(acd_ref_v), levels = 1:10)))
  d_new  <- as.numeric(table(factor(dec_of(acd_new_v), levels = 1:10)))
  
  # ACD-space coverage: share of the project area whose ACD sits within ~1 bin
  # (5 Mg C/ha) of at least one plot value
  bw  <- 5
  brk <- seq(floor(min(v_all) / bw) * bw, ceiling(max(v_all) / bw) * bw + bw, bw)
  bin <- function(x) cut(x, brk, labels = FALSE, include.lowest = TRUE)
  expand1 <- function(o) sort(unique(c(o - 1, o, o + 1)))
  cov_ref <- 100 * mean(bin(v_all) %in% expand1(unique(bin(acd_ref_v))))
  cov_all <- 100 * mean(bin(v_all) %in% expand1(unique(bin(acd_all_v))))
  
  # ---------------------------------------------------------------------------
  #  3. eight-panel diagnostics figure
  # ---------------------------------------------------------------------------
  open_png("_diagnostics.png", 14, 8)
  par(mfrow = c(2, 4), mar = c(4.2, 4.2, 3, 1.2), mgp = c(2.4, 0.7, 0),
      cex.main = 1, font.main = 1, las = 1)
  
  ## (a) where the plots sit in the project ACD distribution
  h <- hist(vs, breaks = 40, plot = FALSE)
  plot(h, col = "grey88", border = "white", freq = TRUE,
       main = "(a) ACD distribution", xlab = "ACD (Mg C/ha)", ylab = "pixels")
  abline(v = c(t33, t67), lty = 2, col = "grey30")
  rug(acd_ref_v, col = COL_REF, lwd = 1.5)
  rug(acd_new_v, col = COL_NEW, lwd = 1.5, side = 3)
  legend("topright", bty = "n", cex = 0.75, lty = 1, lwd = 2,
         col = c(COL_REF, COL_NEW), legend = c("existing (below)", "new (above)"))
  
  ## (b) paired shift in ACD
  o <- order(new$ACD_ref); d <- new[o, ]; yy <- seq_len(nrow(d))
  plot(range(c(d$ACD_ref, d$ACD_new)), range(yy), type = "n", yaxt = "n",
       main = "(b) reference -> paired plot", xlab = "ACD (Mg C/ha)", ylab = "")
  abline(v = c(t33, t67), lty = 3, col = "grey60")
  arrows(d$ACD_ref, yy, d$ACD_new, yy, length = 0.05, col = "grey45")
  points(d$ACD_ref, yy, pch = 19, col = COL_REF, cex = 0.7)
  points(d$ACD_new, yy, pch = 17, col = COL_NEW, cex = 0.7)
  axis(2, at = yy, labels = d$Pl_ref, cex.axis = 0.45, tick = FALSE)
  
  ## (c) share of the ACD terciles
  bp <- barplot(M, beside = TRUE, ylim = c(0, max(M) * 1.25),
                col = c("grey75", COL_REF, COL_NEW, "#7D3C98"),
                border = NA, main = "(c) ACD tercile shares",
                xlab = "ACD class", ylab = "% of area / % of plots")
  legend("topleft", bty = "n", cex = 0.7, fill = c("grey75", COL_REF, COL_NEW,
                                                   "#7D3C98"), border = NA, legend = rownames(M))
  text(bp, M + max(M) * 0.03, sprintf("%.0f", M), cex = 0.55, col = "grey20")
  
  ## (d) occupancy of the ACD deciles
  barplot(rbind(d_ref, d_new), beside = FALSE, border = NA,
          col = c(COL_REF, COL_NEW), names.arg = 1:10,
          main = "(d) plots per ACD decile", xlab = "project ACD decile (10% of area each)",
          ylab = "number of plots")
  legend("topright", bty = "n", cex = 0.7, fill = c(COL_REF, COL_NEW),
         border = NA, legend = c("existing", "new"))
  mtext(sprintf("deciles occupied: %d/10 existing -> %d/10 with new plots",
                sum(d_ref > 0), sum((d_ref + d_new) > 0)),
        side = 3, line = -1.1, cex = 0.6, col = "grey25")
  
  ## (e) achieved contrast
  o2 <- order(new$dACD)
  barplot(new$dACD[o2], border = NA, space = 0.15,
          col = ifelse(new$target[o2] == "high", "#2E86C1", "#E67E22"),
          main = "(e) contrast achieved", ylab = "|dACD| (Mg C/ha)", xlab = "pair")
  abline(h = DELTA_LEVELS[DELTA_LEVELS > 0], lty = 3, col = "grey40")
  legend("topleft", bty = "n", cex = 0.7, fill = c("#2E86C1", "#E67E22"),
         border = NA, legend = c("target higher ACD", "target lower ACD"))
  
  ## (f) ACD of existing vs new plots
  boxplot(list(existing = acd_ref_v, new = acd_new_v), col = c(COL_REF, COL_NEW),
          border = "grey25", outline = FALSE, main = "(f) plot-level ACD",
          ylab = "ACD (Mg C/ha)")
  set.seed(1)
  points(jitter(rep(1, length(acd_ref_v)), 8), acd_ref_v, pch = 16, cex = 0.5,
         col = adjustcolor("white", 0.9))
  points(jitter(rep(2, length(acd_new_v)), 8), acd_new_v, pch = 16, cex = 0.5,
         col = adjustcolor("white", 0.9))
  abline(h = c(t33, t67), lty = 3, col = "grey50")
  
  ## (g) access: distance to the nearest road
  plot(ecdf(as.numeric(na.omit(plots$road_ref))), col = COL_REF, lwd = 2,
       do.points = FALSE, verticals = TRUE, main = "(g) distance to road",
       xlab = "distance (m)", ylab = "cumulative share of plots",
       xlim = c(0, max(c(plots$road_ref, new$road_new), na.rm = TRUE)))
  plot(ecdf(new$road_new), col = COL_NEW, lwd = 2, do.points = FALSE,
       verticals = TRUE, add = TRUE)
  legend("bottomright", bty = "n", cex = 0.75, lty = 1, lwd = 2,
         col = c(COL_REF, COL_NEW), legend = c("existing", "new"))
  
  ## (h) pair geometry and local homogeneity
  hist(new$d_pair, breaks = seq(PAIR_MIN, PAIR_MAX, by = 20), col = "grey80",
       border = "white", main = "(h) pair distance / ACD sd",
       xlab = "distance to reference plot (m)", ylab = "pairs")
  box()
  par(new = TRUE)
  plot(new$d_pair, new$acd_sd, axes = FALSE, xlab = "", ylab = "", pch = 21,
       bg = adjustcolor(COL_NEW, 0.6), col = "white", cex = 0.9,
       xlim = c(PAIR_MIN, PAIR_MAX))
  axis(4, col.axis = COL_NEW, col = COL_NEW, cex.axis = 0.8)
  mtext("ACD sd in 90 m window (Mg C/ha)", side = 4, line = 2, cex = 0.6,
        col = COL_NEW, las = 0)
  
  mtext(sprintf("Kuamut IFM paired plots - %d new plots for %d existing plots  |  ACD terciles %.0f / %.0f Mg C/ha",
                nrow(new), nrow(plots), t33, t67),
        outer = TRUE, line = -1.4, cex = 0.8)
  dev.off()
  
  # ---------------------------------------------------------------------------
  #  4. per-pair zoom panels
  # ---------------------------------------------------------------------------
  per_page <- 12
  pages <- split(seq_len(nrow(new)), ceiling(seq_len(nrow(new)) / per_page))
  for (p in seq_along(pages)) {
    png(file.path(DIR_OUT, sprintf("%s_pairs_%02d.png", STEM, p)),
        width = 11, height = 8.5, units = "in", res = 200)
    par(mfrow = c(3, 4), oma = c(0, 0, 2, 0))
    for (k in pages[[p]]) {
      rx <- new_ref_xy[k, ]
      if (anyNA(rx)) next
      e  <- ext(min(rx[1], new$x[k]) - 250, max(rx[1], new$x[k]) + 250,
                min(rx[2], new$y[k]) - 250, max(rx[2], new$y[k]) + 250)
      r  <- crop(acd, e)
      plot(r, col = RAMP, legend = FALSE, axes = FALSE, mar = c(0.5, 0.5, 2, 0.5),
           main = sprintf("%s  %.0f -> %.0f  (d=%.0f m)", new$Pl_new[k],
                          new$ACD_ref[k], new$ACD_new[k], new$d_pair[k]),
           cex.main = 0.85)
      plot(st_geometry(roads), add = TRUE, col = "grey30", lwd = 0.8)
      segments(rx[1], rx[2], new$x[k], new$y[k], lwd = 1.2)
      points(rx[1], rx[2], pch = 21, bg = COL_REF, col = "white", cex = 1.3)
      points(new$x[k], new$y[k], pch = 24, bg = COL_NEW, col = "white", cex = 1.3)
      box(col = "grey70")
    }
    mtext("Pair detail - reference plot (circle) and new plot (triangle) on ACD",
          outer = TRUE, line = 0.3, cex = 0.9)
    dev.off()
  }
  
  # ---------------------------------------------------------------------------
  #  5. summary table + console report
  # ---------------------------------------------------------------------------
  summ <- data.frame(
    metric = c("existing plots", "new plots",
               "project area low / mid / high (%)",
               "existing plots low / mid / high (%)",
               "existing+new low / mid / high (%)",
               "ACD deciles occupied - existing",
               "ACD deciles occupied - existing+new",
               "ACD space covered - existing (%)",
               "ACD space covered - existing+new (%)",
               "median |dACD| (Mg C/ha)", "min |dACD| (Mg C/ha)",
               "mean ACD existing / new (Mg C/ha)",
               "median road distance new (m)", "pair distance range (m)"),
    value = c(nrow(plots), nrow(new),
              paste(sprintf("%.0f", M["project area", ]), collapse = " / "),
              paste(sprintf("%.0f", M["existing only", ]), collapse = " / "),
              paste(sprintf("%.0f", M["existing+new", ]), collapse = " / "),
              sprintf("%d/10", sum(d_ref > 0)),
              sprintf("%d/10", sum((d_ref + d_new) > 0)),
              sprintf("%.0f", cov_ref), sprintf("%.0f", cov_all),
              sprintf("%.1f", median(new$dACD)), sprintf("%.1f", min(new$dACD)),
              sprintf("%.0f / %.0f", mean(acd_ref_v), mean(acd_new_v)),
              sprintf("%.0f", median(new$road_new)),
              sprintf("%.0f-%.0f", min(new$d_pair), max(new$d_pair))),
    stringsAsFactors = FALSE)
  write.csv(summ, file.path(DIR_OUT, paste0(STEM, "_summary.csv")),
            row.names = FALSE)
  cat("\n--- ACD representation summary ---\n")
  print(summ, right = FALSE, row.names = FALSE)
  message("figures written to ", DIR_OUT)
}