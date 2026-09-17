# =============================================================================
#  Kuamut IFM (Sabah, MY) - contrasting paired plots
#
#  MODE = "local" : full pipeline from the local files (ACD map + roads +
#                   project area + existing plots). No Earth Engine needed.
#                   This reproduces Kuamut_paired_plots_2025_v3.* exactly.
#  MODE = "gee"   : read the candidate table exported by
#                   01_gee_kuamut_paired_plots_v3.js (which adds the AlphaEarth
#                   Satellite Embedding dissimilarity to the score) and run the
#                   final greedy assignment that enforces the 100 m separation
#                   between the new plots.
#
#  Requires: terra, sf, dplyr
# =============================================================================
#install.packages("terra")
library(terra); library(sf); library(dplyr)

MODE <- "local"        # "local" or "gee"

# ------------------------------- paths ---------------------------------------
DIR_IN   <- "C:/Users/JavierRuizRamos/OneDrive - Permian Global Research Limited/Desktop/Kuamut/VegetationPlot campaing_20252026/Kuamut_PlotLocationsCheck_092026/KuamutPlotPairs_Claude/CodebaseData"                     # folder holding the inputs
DIR_OUT  <- "C:/Users/JavierRuizRamos/OneDrive - Permian Global Research Limited/Desktop/Kuamut/VegetationPlot campaing_20252026/Kuamut_PlotLocationsCheck_092026/KuamutPlotPairs_Claude/CodeRResults"
F_ACD    <- file.path(DIR_IN, "acd_2025_30m.tif")
F_PLOTS  <- file.path(DIR_IN, "KuamutLocations.shp")
F_AOI    <- file.path(DIR_IN, "Kuamut_ProjectArea.shp")
F_ROADS  <- file.path(DIR_IN, "Permian Roads- all.shp")
F_GEECSV <- file.path(DIR_IN, "Kuamut_paired_plot_candidates.csv")
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
  by_ref <- by_ref[order(sapply(by_ref, nrow))]   # fewest options go first

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
    worst_first <- names(picked)[order(sapply(picked, function(p) p$score))]
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
  # pair links, useful for a quick visual check in QGIS
  ref_xy <- st_coordinates(plots)[match(g$Pl_ref, plots[[PLOT_ID_FIELD]]), ]
  links <- st_sfc(lapply(seq_len(nrow(g)), function(i)
    st_linestring(rbind(ref_xy[i, ], c(g$x[i], g$y[i])))), crs = EPSG)
  st_write(st_sf(Pl_ref = g$Pl_ref, geometry = links),
           file.path(DIR_OUT, paste0(stem, ".gpkg")), layer = "pair_links",
           append = TRUE, quiet = TRUE)
  st_write(ll, file.path(DIR_OUT, paste0(stem, ".kml")),
           delete_dsn = TRUE, quiet = TRUE)
  write.csv(st_drop_geometry(g), file.path(DIR_OUT, paste0(stem, ".csv")),
            row.names = FALSE)
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

  # --- ACD as the plot will actually measure it -------------------------------
  # A 60 m plot spans about four 30 m pixels, so the value a field crew will
  # record is the footprint mean, not the centre pixel. Weights come from a
  # circular kernel of radius PLOT_RADIUS and are renormalised by the valid
  # weight so cells at the project edge are not dragged towards zero.
  w_plot   <- focalMat(acd, PLOT_RADIUS, "circle")
  acd_sum  <- focal(acd, w = w_plot, fun = "sum", na.rm = TRUE)
  wgt_sum  <- focal(!is.na(acd), w = w_plot, fun = "sum", na.rm = TRUE)
  acd_plot <- acd_sum / wgt_sum
  names(acd_plot) <- "acd"          # keeps the downstream column name `acd`

  # --- ACD heterogeneity over the 90 m (3 x 3) neighbourhood ----------------
  # 90 m covers the 60 m footprint plus a 15 m margin, so this now tests
  # whether the plot itself straddles an edge
  acd_sd <- focal(acd, w = 3, fun = sd, na.rm = TRUE, na.policy = "omit")
  names(acd_sd) <- "acd_sd"

  # --- distance to the nearest road -----------------------------------------
  road_r <- rasterize(vect(roads), rast(acd), field = 1, touches = TRUE)
  road_d <- distance(road_r)            # distance to the nearest non-NA cell
  names(road_d) <- "road_d"

  # --- project-wide ACD distribution ----------------------------------------
  # per-pixel, because this describes the area; plot values are footprint means
  # and so sit slightly closer to the centre of this distribution
  v_all <- values(acd, mat = FALSE); v_all <- v_all[!is.na(v_all)]
  qs  <- quantile(v_all, c(1/3, 0.5, 2/3))
  t33 <- qs[[1]]; med <- qs[[2]]; t67 <- qs[[3]]
  message(sprintf("ACD terciles: <%.1f | %.1f-%.1f | >%.1f Mg C/ha",
                  t33, t33, t67, t67))
  acd_class <- function(v) ifelse(v < t33, "low", ifelse(v > t67, "high", "mid"))

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
  TARGET <<- 2 * nrow(plots) * U / sum(U)
  message("decile edges: ", paste(round(EDGES[2:N_BINS], 1), collapse = " "))
  message("target counts: ", paste(round(TARGET, 1), collapse = " "))

  # --- reference plot properties --------------------------------------------
  pxy <- st_coordinates(plots)
  plots$ACD_ref  <- terra::extract(acd_plot, pxy)[, 1]   # footprint mean
  plots$road_ref <- as.numeric(st_distance(plots, st_union(roads)))
  plots$target   <- ifelse(plots$ACD_ref > med, "low", "high")

  # --- candidate offsets: a 10 m grid inside the 100-300 m annulus ----------
  off <- seq(-PAIR_MAX, PAIR_MAX, by = GRID)
  gr  <- expand.grid(dx = off, dy = off)
  gr$d <- sqrt(gr$dx^2 + gr$dy^2)
  gr  <- gr[gr$d >= PAIR_MIN & gr$d <= PAIR_MAX, ]

  cand_list <- vector("list", nrow(plots))
  for (i in seq_len(nrow(plots))) {
    cx <- pxy[i, 1] + gr$dx; cy <- pxy[i, 2] + gr$dy
    xy <- cbind(cx, cy)
    v  <- terra::extract(c(acd_plot, acd_sd, road_d), xy)
    road_cap <- max(ROAD_CAP_FLOOR, plots$road_ref[i])
    keep <- !is.na(v$acd) & !is.na(v$acd_sd) &
            v$acd >= MIN_ACD & v$acd_sd <= SD_HARD &
            v$road_d >= ROAD_MIN & v$road_d <= road_cap
    if (!any(keep)) { cand_list[[i]] <- NULL; next }
    xy <- xy[keep, , drop = FALSE]; v <- v[keep, ]; dpk <- gr$d[keep]
    # at least MIN_SEP from every existing plot centre, and - the
    # nearest-neighbour rule - closer to its own reference plot than to any
    # other existing plot
    dall <- apply(xy, 1, function(p)
      sqrt((pxy[, 1] - p[1])^2 + (pxy[, 2] - p[2])^2))          # nplots x ncand
    dmin     <- apply(dall, 2, min)
    dmin_oth <- apply(dall[-i, , drop = FALSE], 2, min)
    ok <- dmin >= MIN_SEP & dmin_oth >= dpk + NN_MARGIN
    if (!any(ok)) { cand_list[[i]] <- NULL; next }
    ref <- plots$ACD_ref[i]
    dacd <- if (DIRECTION == "free") abs(v$acd[ok] - ref) else
            if (plots$target[i] == "low") ref - v$acd[ok] else v$acd[ok] - ref
    cand_list[[i]] <- data.frame(
      Pl_ref = plots[[PLOT_ID_FIELD]][i], target = plots$target[i],
      ACD_ref = ref, ACD_new = v$acd[ok], dACD = dacd,
      acd_sd = v$acd_sd[ok], road_new = v$road_d[ok],
      road_ref = plots$road_ref[i], d_pair = dpk[ok], d_other = dmin_oth[ok],
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
  n_ref <- tabulate(acd_bin(plots$ACD_ref), nbins = N_BINS)
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
                "dACD", "bin", "d_pair", "d_other", "road_new", "road_exact")],
        digits = 4)
  cat(sprintf("\n|dACD| median %.1f, min %.1f Mg C/ha\n",
              median(new$dACD), min(new$dACD)))
  write_outputs(new, plots, "Kuamut_paired_plots_2025_v3")
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
  write_outputs(new, plots, "Kuamut_paired_plots_2025_embeddings_v3")
}

# =============================================================================
#  final sanity checks - run on whichever mode produced `new`
# =============================================================================
all_xy <- rbind(st_coordinates(plots), as.matrix(new[, c("x", "y")]))
dm <- as.matrix(dist(all_xy)); diag(dm) <- Inf
cat(sprintf("plots total: %d | closest pair overall: %.1f m (%.1f m between %d m footprints)\n",
            nrow(all_xy), min(dm), min(dm) - 2 * PLOT_RADIUS, PLOT_RADIUS))
cat(sprintf("pair distances: %.0f-%.0f m (target %d-%d)\n",
            min(new$d_pair), max(new$d_pair), PAIR_MIN, PAIR_MAX))
cat(sprintf("road distance of new plots: median %.0f m, max %.0f m\n",
            median(new$road_new), max(new$road_new)))

# nearest-neighbour rule: for every new plot the closest plot must be its pair
ref_xy <- st_coordinates(plots)
lab <- c(paste0("P", plots[[PLOT_ID_FIELD]]), paste0(new$Pl_ref, "B"))
dm2 <- as.matrix(dist(all_xy)); diag(dm2) <- Inf
viol <- 0
for (k in seq_len(nrow(new))) {
  idx <- nrow(ref_xy) + k
  own <- match(new$Pl_ref[k], plots[[PLOT_ID_FIELD]])
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
#  Append this to 02_kuamut_paired_plots_v3.R, after the final sanity checks.
#  It reuses objects built in the "local" block:
#     acd, acd_sd, road_d, aoi, roads, plots, new, t33, med, t67, v_all,
#     acd_class, PLOT_ID_FIELD, DIR_OUT, PAIR_MIN, PAIR_MAX, DELTA_LEVELS
#
#  Outputs in DIR_OUT:
#     Kuamut_paired_plots_2025_v3_map.png          overview map
#     Kuamut_paired_plots_2025_v3_diagnostics.png  8-panel figure
#     Kuamut_paired_plots_2025_v3_pairs_NN.png     per-pair zoom panels
#     Kuamut_paired_plots_2025_v3_summary.csv      the numbers behind the figures
# =============================================================================
if (MODE == "local") {
  
  STEM    <- "Kuamut_paired_plots_2025_v3"
  COL_REF <- "#1B4F72"                                  # existing plots
  COL_NEW <- "#C0392B"                                  # new paired plots
  RAMP    <- hcl.colors(100, "Greens", rev = TRUE)      # light = low ACD
  
  open_png <- function(f, w = 11, h = 9, res = 200)
    png(file.path(DIR_OUT, paste0(STEM, f)),
        width = w, height = h, units = "in", res = res)
  
  ref_xy     <- st_coordinates(plots)
  ref_i      <- match(new$Pl_ref, plots[[PLOT_ID_FIELD]])
  new_ref_xy <- ref_xy[ref_i, , drop = FALSE]
  
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