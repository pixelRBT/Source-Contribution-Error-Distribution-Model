#!/usr/bin/env Rscript
# =============================================================================
# SCED: Source Contribution Error Distribution | Prediction tool v2.1 (NoSmin)
# =============================================================================
# PURPOSE
#   Use ONE case manifest CSV (or XLSX with readxl installed), listing the
#   mixture/source CSV and completed Rjags .Rdata for every independent case.
#   Use the frozen SCED Beta model to predict source-specific AE distributions.
#   AE = absolute estimation error in source contribution (0–1 scale).
#   MixSIAR remains the underlying isotope mixing model; SCED predicts its error risk.
#   Discrimination CSV is NOT used as a predictor or input.
#
# ASSUMPTIONS: paired d2H/d18O, 3–7 source groups, one ungrouped MixSIAR fit
#   per case, zero isotope discrimination MEAN AND SD, no concentration dependence.
#   The zero-discrimination assumption is user responsibility, because the
#   discrimination file is intentionally not inspected by this tool.
#
# MANIFEST: Case_ID, Mixture_File, Source_File, Rjags_File
#   Paths are relative to the manifest's folder, or may be absolute.
#   Source files can be reused across multiple case rows.
#
# QUICK START:
#   1) Save Case_List.csv with the four columns above.
#   2) Set CASE_LIST_FILE, TRACER_D2H, TRACER_D18O below.
#   3) Run: source("SCED_Prediction_Tool.R")
#
# COMMAND LINE:
#   Rscript SCED_Prediction_Tool.R Case_List.csv [output_dir]
#
# OUTPUT: single batch HTML overview + one combined source CSV + individual
#   illustrated HTML reports. Invalid cases are logged without losing others.
#
# REQUIRES: base R only for CSV; readxl is optional for an XLSX manifest.
# =============================================================================

# =================== USER SETTINGS (edit here) =============================
CASE_LIST_FILE <- "Case_List.csv"
TRACER_D2H <- "H"       # Actual column name of deuterium tracer, e.g. H
TRACER_D18O <- "O"      # Actual column name of oxygen-18 tracer, e.g. O
AE_THRESHOLD <- 0.20    # Probability AE exceeds this contribution fraction
OUTPUT_DIR <- "SCED_Batch_Report"
AUTO_RUN <- TRUE

# The frozen model was calibrated with zero isotope discrimination and no
# concentration dependence. Nonzero discrimination input is NOT supported.
# There is no discrimination file input: users must check this assumption.
# ==========================================================================

# =============================================================================
# FROZEN SCED MODEL (model C; NoSmin OOF full fit, 2026-11-06)
# 998 complete scenarios / 4,994 source-level records; scenario-equal weights.
# Identical conditional mean equation used by the fixed-precision model in 3.3.
# =============================================================================
.SCED_MODEL <- list(
  name = "Source Contribution Error Distribution (SCED; NoSmin, fixed Beta precision)",
  n_scenarios = 998L,
  n_sources = 4994L,
  precision = 11.7368841924365,
  coef = c(
    "(Intercept)" = -2.1360669420759502,
    "NOS4" = 0.118483766158651,
    "NOS5" = 0.018184622692306401,
    "NOS6" = -0.0015444053704498,
    "NOS7" = -0.0118335362242532,
    "z_w" = 0.182819829903792,
    "z_near" = -0.070411388817612205,
    "z_rpos" = 0.0703416914263503,
    "z_dfm" = -0.0928172547826914,
    "z_w:z_rpos" = 0.0712064212209364
  ),
  center = c(
    z_near = -0.36547124335366599,
    z_rpos = 0.50592921420320902,
    z_dfm = 0.34989353490362501,
    z_w = 0.43945515648526701
  ),
  scale = c(
    z_near = 0.61986561076790103,
    z_rpos = 0.25708102859283699,
    z_dfm = 0.20196545195860199,
    z_w = 0.142508100552103
  )
)

# =============================================================================
# Input/validation helpers
# =============================================================================
.require_file <- function(path, label) {
  if (is.na(path) || !is.character(path) || length(path) != 1L || !nzchar(path))
    stop("Please specify the file for: ", label, call. = FALSE)
  if (!file.exists(path)) stop("File not found (", label, "): ", path, call. = FALSE)
  normalizePath(path, mustWork = TRUE)
}
.read_csv <- function(path, label) {
  x <- tryCatch(utils::read.csv(path, check.names = FALSE,
                                 stringsAsFactors = FALSE,
                                 fileEncoding = "UTF-8-BOM"),
                error = function(e) stop("Cannot read ", label, ": ",
                                         conditionMessage(e), call. = FALSE))
  if (!nrow(x)) stop(label, " contains no rows", call. = FALSE)
  x
}
.key <- function(x) gsub("[^a-z0-9]", "", tolower(x))
.col <- function(data, labels, what, must = TRUE) {
  target <- .key(labels)
  ix <- match(target, .key(names(data)))
  ix <- ix[!is.na(ix)]
  if (length(ix)) return(names(data)[ix[1L]])
  if (must) stop("Missing ", what, ". Expected one of: ",
                 paste(labels, collapse = ", "), ". Available: ",
                 paste(names(data), collapse = ", "), call. = FALSE)
  NULL
}
# A SINGLE public tracer-name interface replaces guessed alias lists.
# The user-specified base name is used consistently in three MixSIAR inputs:
#   mixture/raw-source: tracer
#   source-means/discrimination: Mean<tracer>, SD<tracer>
# Dots/underscores/spacing in header spellings are normalized, but the
# name itself MUST come from the two user-supplied declarations.
.tracer_col <- function(data, tracer, kind = c("raw", "mean", "sd"),
                        what = "tracer field", must = TRUE) {
  kind <- match.arg(kind)
  labels <- switch(kind,
    raw  = c(tracer),
    mean = c(paste0("Mean", tracer), paste0(tracer, "Mean")),
    sd   = c(paste0("SD", tracer), paste0(tracer, "SD"))
  )
  .col(data, labels, what, must = must)
}
.validate_tracers <- function(tracer_d2h, tracer_d18o) {
  values <- c(d2H = tracer_d2h, d18O = tracer_d18o)
  if (!is.character(values) || length(values) != 2L ||
      anyNA(values) || any(!nzchar(trimws(values))) ||
      length(unique(.key(values))) != 2L) {
    stop("Set two distinct nonempty names at TRACER_D2H and TRACER_D18O. ",
         "For standard MixSIAR water-isotope files: 'H' and 'O'.", call. = FALSE)
  }
  values
}
.to_num <- function(v, what) {
  x <- suppressWarnings(as.numeric(as.character(v)))
  if (any(!is.finite(x))) stop("Missing or nonnumeric values in ", what,
                               call. = FALSE)
  x
}
.escape <- function(s) {
  x <- as.character(s)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}
.pct <- function(x, digits = 1L) sprintf(paste0("%.", digits, "f%%"), 100*x)
.pp <- function(x) sprintf("%.1f pp", 100*x)
.svg_fmt <- function(x) formatC(as.numeric(x), digits = 2L, format = "f")

.read_mixture <- function(file, tracers) {
  d <- .read_csv(file, "mixture CSV")
  o <- .tracer_col(d, tracers["d18O"], "raw", "mixture d18O")
  h <- .tracer_col(d, tracers["d2H"], "raw", "mixture d2H")
  message(sprintf("Mixture columns: d18O <- '%s'; d2H <- '%s'", o, h))
  # Complete, paired observations, exactly as in the MixSIAR input.
  x <- .to_num(d[[o]], "mixture d18O")
  y <- .to_num(d[[h]], "mixture d2H")
  if (length(x) < 2L) stop("At least two paired mixture observations are required")
  data.frame(d18O = x, d2H = y)
}
.read_sources <- function(file, tracers) {
  d <- .read_csv(file, "source CSV")
  name_col <- .col(d, c("Source", "Sources", "SourceName", "Source_Name"),
                   "source names", must = FALSE)
  # MixSIAR explicitly interprets the first column as source labels.
  if (is.null(name_col)) name_col <- names(d)[1L]
  names_raw <- trimws(as.character(d[[name_col]]))
  if (anyNA(names_raw) || any(!nzchar(names_raw)))
    stop("Source names cannot be blank")
  # Match MixSIAR's alphabetical source ordering, not original file row order.
  ord_names <- sort(unique(names_raw))
  K <- length(ord_names)
  if (!(K %in% 3:7)) stop("The calibrated model supports 3 to 7 sources; found ", K)

  mo <- .tracer_col(d, tracers["d18O"], "mean", "source mean d18O", must = FALSE)
  mh <- .tracer_col(d, tracers["d2H"], "mean", "source mean d2H", must = FALSE)
  so <- .tracer_col(d, tracers["d18O"], "sd", "source SD d18O", must = FALSE)
  sh <- .tracer_col(d, tracers["d2H"], "sd", "source SD d2H", must = FALSE)
  n_mean_fields <- sum(!vapply(list(mo,mh,so,sh), is.null, logical(1)))
  if (n_mean_fields > 0L && n_mean_fields < 4L)
    stop("Incomplete source means/SD fields for configured tracers ",
         tracers["d2H"], " and ", tracers["d18O"],
         ". Expected Mean",tracers["d2H"], ", SD",tracers["d2H"],
         ", Mean",tracers["d18O"], ", SD",tracers["d18O"], ".")
  nn <- .col(d, c("n", "N", "SampleN"), "source sample size", must = FALSE)

  if (any(grepl("^conc", names(d), ignore.case = TRUE)))
    stop("Concentration dependence is not supported by this SCED implementation.")

  if (!is.null(mo) && !is.null(mh) && !is.null(so) && !is.null(sh)) {
    if (is.null(nn)) stop("Means-format source CSV requires the MixSIAR n column")
    if (anyDuplicated(names_raw)) stop(
      "Means-format source CSV has repeated names. This tool supports one source profile per source, not source-by-factor fits.")
    i <- match(ord_names, names_raw)
    message(sprintf("Recognized source means: d18O <- '%s'; d2H <- '%s'; SD18O <- '%s'; SD2H <- '%s'",
                    mo, mh, so, sh))
    ans <- data.frame(
      source = ord_names,
      d18O = .to_num(d[[mo]][i], "source mean d18O"),
      d2H = .to_num(d[[mh]][i], "source mean d2H"),
      sd18 = .to_num(d[[so]][i], "source SD d18O"),
      sd2H = .to_num(d[[sh]][i], "source SD d2H"),
      n = .to_num(d[[nn]][i], "source n"),
      stringsAsFactors = FALSE
    )
    mode <- "MixSIAR source means + SD + n"
  } else {
    ro <- .tracer_col(d, tracers["d18O"], "raw", "raw source d18O")
    rh <- .tracer_col(d, tracers["d2H"], "raw", "raw source d2H")
    message(sprintf("Recognized raw source isotopes: d18O <- '%s'; d2H <- '%s'", ro, rh))
    vals <- data.frame(source = names_raw,
                       d18O = .to_num(d[[ro]], "raw source d18O"),
                       d2H = .to_num(d[[rh]], "raw source d2H"))
    grp <- split(seq_len(nrow(vals)), factor(vals$source, levels = ord_names))
    ans <- do.call(rbind, lapply(ord_names, function(nm) {
      v <- vals[grp[[nm]], , drop = FALSE]
      if (nrow(v) < 2L) stop("Source '", nm,
                            "' has fewer than 2 raw paired observations")
      data.frame(source = nm, d18O = mean(v$d18O), d2H = mean(v$d2H),
                 sd18 = stats::sd(v$d18O), sd2H = stats::sd(v$d2H),
                 n = nrow(v), stringsAsFactors = FALSE)
    }))
    rownames(ans) <- NULL
    mode <- "MixSIAR raw source measurements"
  }
  if (any(ans$sd18 <= 0 | ans$sd2H <= 0 | ans$n < 2 |
          ans$n != as.integer(ans$n)))
    stop("Each source requires positive SDs and an integer sample size >= 2")
  list(data = ans, mode = mode)
}
.get_rjags <- function(file, K) {
  env <- new.env(parent = baseenv())
  objects <- load(file, envir = env)
  if (!length(objects)) stop("The .Rdata file has no R objects")
  candidates <- objects[vapply(objects, function(nm) {
    x <- get(nm, envir = env)
    is.list(x) && is.list(x$BUGSoutput) &&
      !is.null(x$BUGSoutput$sims.list$p.global)
  }, logical(1))]
  if ("run" %in% candidates) {
    chosen <- "run"
  } else if (length(candidates) == 1L) {
    chosen <- candidates[1L]
  } else {
    stop("Cannot uniquely find the MixSIAR R2jags 'run' object with p.global in the .Rdata file. Found: ",
         paste(candidates, collapse = ", "))
  }
  run <- get(chosen, envir = env)
  sim <- run$BUGSoutput$sims.list
  if (any(grepl("^p\\.fac|^p\\.ind|^p\\.both", names(sim)))) stop(
    "A grouped or individual-effect MixSIAR posterior was detected. This tool predicts the single ungrouped p.global scenario only; group-specific deployment requires separate calibration.")
  p <- sim$p.global
  dp <- dim(p)
  if (is.null(dp)) stop("p.global is not a matrix/array of posterior draws")
  if (length(dp) == 2L) {
    p <- as.matrix(p)
  } else if (length(dp) == 3L && dp[3L] == K) {
    p <- matrix(as.numeric(p), ncol = K)
  } else {
    stop("Unsupported p.global dimensions: ", paste(dp, collapse = " x "))
  }
  if (nrow(p) < 100L || ncol(p) != K || any(!is.finite(p)) ||
      any(p < -1e-7 | p > 1+1e-7) ||
      max(abs(rowSums(p)-1)) > 1e-4) stop(
    "p.global posterior draws are incomplete, have wrong source count, or do not sum to one")
  col_ids <- colnames(p)
  if (!is.null(col_ids) && length(col_ids) == K &&
      all(grepl("^p\\.global\\[[0-9]+\\]$", col_ids))) {
    indices <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", col_ids))
    if (!setequal(indices, seq_len(K))) stop("Unexpected p.global source indices")
    p <- p[, order(indices), drop = FALSE]
  }
  post <- data.frame(
    mean = colMeans(p),
    lo = apply(p, 2L, stats::quantile, probs = 0.025, names = FALSE),
    hi = apply(p, 2L, stats::quantile, probs = 0.975, names = FALSE)
  )
  if (any(post$hi <= post$lo | post$mean < 0 | post$mean > 1))
    stop("Invalid source contribution posterior quantiles")
  su <- run$BUGSoutput$summary
  rhat <- NA_real_; ess <- NA_real_
  if (is.matrix(su) && !is.null(rownames(su))) {
    ix <- grep("^p\\.global\\[[[:space:]]*[0-9]+[[:space:]]*\\]$", rownames(su))
    if (length(ix) == K) {
      cn <- tolower(colnames(su))
      ir <- match("rhat", cn)
      ie <- match(TRUE, cn %in% c("n.eff", "neff", "ess"))
      if (!is.na(ir)) {
        z <- suppressWarnings(as.numeric(su[ix, ir]));
        if (all(is.finite(z))) rhat <- max(z)
      }
      if (!is.na(ie)) {
        z <- suppressWarnings(as.numeric(su[ix, ie]));
        if (all(is.finite(z))) ess <- min(z)
      }
    }
  }
  list(post = post, n_draws = nrow(p), rhat_max = rhat,
       ess_min = ess, object = chosen)
}

# =============================================================================
# Geometry consistent with the updated NoSmin model and Methods 2.3.
# S_ij = sqrt((d18_i-d18_j)^2/(SD18_i^2+SD18_j^2) +
#             (d2H_i-d2H_j)^2/(SD2H_i^2+SD2H_j^2))
# Rpos from arithmetic mean of polygon VERTICES and mean mixture point,
# not the polygon area centroid. DFM is axis-range-scaled distance / sqrt(2).
# =============================================================================
.cross <- function(a,b) a[1L]*b[2L]-a[2L]*b[1L]
.geometry <- function(src, mix) {
  K <- nrow(src)
  m <- c(mean(mix$d18O), mean(mix$d2H))
  M <- as.matrix(src[, c("d18O", "d2H")]); storage.mode(M) <- "double"
  distance <- matrix(0, K, K)
  for (i in seq_len(K-1L)) for (j in (i+1L):K) {
    distance[i,j] <- distance[j,i] <- sqrt(
      (src$d18O[i]-src$d18O[j])^2/(src$sd18[i]^2+src$sd18[j]^2) +
      (src$d2H[i]-src$d2H[j])^2/(src$sd2H[i]^2+src$sd2H[j]^2))
  }
  if (any(!is.finite(distance)) || any(distance[upper.tri(distance)] <= 0))
    stop("Two isotope sources have identical means or invalid pairwise distances")
  nearest <- vapply(seq_len(K), function(i) min(distance[i,-i]), numeric(1))
  hull_order <- grDevices::chull(M[,1], M[,2])
  h <- M[hull_order,,drop=FALSE]
  if (nrow(h) < 3L) stop("All source means are collinear: 2-D geometry undefined")
  j <- c(2:nrow(h),1L)
  area2 <- sum(h[,1]*h[j,2] - h[,2]*h[j,1])
  if (abs(area2) < 1e-10) stop("Isotope mixing-space hull is degenerate")
  center <- colMeans(h)
  vec <- m-center
  if (sqrt(sum(vec^2)) < 1e-13) {
    rpos <- 0
  } else {
    tvals <- numeric(0)
    for (z in seq_len(nrow(h))) {
      aa <- as.numeric(h[z,]); bb <- as.numeric(h[j[z],]); edge <- bb-aa
      den <- .cross(vec,edge)
      if (abs(den) < 1e-12) next
      t <- .cross(aa-center,edge)/den
      u <- .cross(aa-center,vec)/den
      if (is.finite(t) && t>1e-12 && u>=-1e-9 && u<=1+1e-9)
        tvals <- c(tvals,t)
    }
    if (!length(tvals)) stop("Could not compute mixture position in isotope hull")
    rpos <- 1/min(tvals)
  }
  if (!is.finite(rpos) || rpos > 1+1e-7 || rpos < -1e-7)
    stop(sprintf("Mixture mean is outside source-means convex hull (Rpos=%.3f). Prediction withheld; check source completeness and isotope assumptions.", rpos))
  rpos <- max(0,min(1,rpos))
  spans <- apply(M, 2L, function(x) diff(range(x)))
  if (any(spans <= 0)) stop("No isotope separation along one axis")
  dfm <- sqrt(((src$d18O-m[1])/spans[1])^2+
                ((src$d2H-m[2])/spans[2])^2)/sqrt(2)
  if (any(!is.finite(dfm))) stop("DFM calculation failed")
  list(nearest = nearest,
       Rpos = rpos, DFM = dfm, mean_mix = m, polygon = h, src_points=M)
}

# =============================================================================
# Responsive SVG charts. Inline and self-contained in HTML: no dependencies.
# =============================================================================
.palette <- c("#256b8e", "#d98248", "#36866b", "#7659a3",
              "#b04d6c", "#a58028", "#476977")
.polyline <- function(x,y,col,width=2.6) {
  pts <- paste0(.svg_fmt(x),",",.svg_fmt(y),collapse=" ")
  sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="%s" stroke-linejoin="round" stroke-linecap="round"/>',
          pts,col,.svg_fmt(width))
}
.line <- function(x1,y1,x2,y2,color="#dfe6ee",dash="") {
  sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="%s" %s/>',
          .svg_fmt(x1),.svg_fmt(y1),.svg_fmt(x2),.svg_fmt(y2),color,
          if (nzchar(dash)) paste0('stroke-dasharray="',dash,'"') else "")
}
.svgtxt <- function(x,y,txt,size=11,color="#344258",anchor="start",weight="normal") {
  sprintf('<text x="%s" y="%s" font-size="%s" fill="%s" text-anchor="%s" font-weight="%s" font-family="Arial, sans-serif">%s</text>',
          .svg_fmt(x),.svg_fmt(y),.svg_fmt(size),color,anchor,weight,.escape(txt))
}
.plot_space_svg <- function(geo, src) {
  W<-570; H<-335; L<-72; R<-30; T<-32; B<-60
  pts <- rbind(geo$src_points, geo$mean_mix)
  xr<-range(pts[,1]); yr<-range(pts[,2]);
  padx <- max(diff(xr)*.10,.1); pady<-max(diff(yr)*.10,.5)
  xr<-xr+c(-padx,padx); yr<-yr+c(-pady,pady)
  X<-function(x) L+(x-xr[1])/diff(xr)*(W-L-R)
  Y<-function(y) H-B-(y-yr[1])/diff(yr)*(H-B-T)
  out<-c(sprintf('<svg viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Source isotope space and observed mixture mean">',W,H))
  for(t in seq(0,1,length.out=5)) {
    xx<-L+t*(W-L-R);yy<-T+t*(H-T-B)
    out<-c(out,.line(xx,T,xx,H-B),.line(L,yy,W-R,yy))
    out<-c(out,.svgtxt(xx,H-B+22,round(xr[1]+t*diff(xr),1),10,anchor="middle"),
           .svgtxt(L-11,yy+3,round(yr[2]-t*diff(yr),1),10,anchor="end"))
  }
  out<-c(out,.line(L,T,L,H-B,"#738296"),.line(L,H-B,W-R,H-B,"#738296"))
  h<-geo$polygon; hj<-c(seq_len(nrow(h)),1L)
  pts_h<-paste0(.svg_fmt(X(h[hj,1])),",",.svg_fmt(Y(h[hj,2])),collapse=" ")
  out<-c(out,sprintf('<polygon points="%s" fill="#d9edf3" fill-opacity=".43" stroke="#6e9fb0" stroke-width="1.7"/>',pts_h))
  for(i in seq_len(nrow(src))) {
    x<-X(src$d18O[i]); y<-Y(src$d2H[i])
    out<-c(out,sprintf('<circle cx="%s" cy="%s" r="5.2" fill="%s" stroke="#fff" stroke-width="1.3"/>',.svg_fmt(x),.svg_fmt(y),.palette[i]),
           .svgtxt(x+9,y-9,src$source[i],11,.palette[i],weight="bold"))
  }
  mx<-X(geo$mean_mix[1]);my<-Y(geo$mean_mix[2])
  diamond<-paste0(.svg_fmt(c(mx,mx+7,mx,mx-7)),",",.svg_fmt(c(my-7,my,my+7,my)),collapse=" ")
  out<-c(out,sprintf('<polygon points="%s" fill="#1c2938" stroke="#fff" stroke-width="1.8"/>',diamond),
         .svgtxt(mx+12,my+17,"Mixture mean",11,"#1c2938",weight="bold"),
         .svgtxt((L+W-R)/2,H-12,"d18O (per mil)",12,anchor="middle",weight="bold"),
         .svgtxt(15,20,"d2H (per mil)",12,weight="bold"),"</svg>")
  paste(out,collapse="\n")
}
.plot_cdf_svg <- function(a,b,q95,threshold,src) {
  W<-570;H<-335;L<-65;R<-25;T<-32;B<-60
  xmax<-min(1,max(.30,threshold*1.35,max(q95)*1.10))
  X<-function(x)L+x/xmax*(W-L-R)
  Y<-function(y)H-B-y*(H-B-T)
  out<-c(sprintf('<svg viewBox="0 0 %d %d" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="Predicted AE cumulative probability distributions">',W,H))
  for(y in seq(0,1,.2)) {
    yy<-Y(y)
    out<-c(out,.line(L,yy,W-R,yy),.svgtxt(L-10,yy+4,sprintf("%d%%",round(100*y)),10,anchor="end"))
  }
  for(x in seq(0,xmax,length.out=5)) {
    xx<-X(x)
    out<-c(out,.line(xx,T,xx,H-B),.svgtxt(xx,H-B+21,sprintf("%.2f",x),10,anchor="middle"))
  }
  out<-c(out,.line(L,T,L,H-B,"#738296"),.line(L,H-B,W-R,H-B,"#738296"))
  out<-c(out,.line(L,Y(.90),W-R,Y(.90),"#8995a7","5 5"),
         .svgtxt(W-R-4,Y(.90)-7,"90%",10,"#5c687a",anchor="end"))
  if (threshold <= xmax) out<-c(out,.line(X(threshold),T,X(threshold),H-B,"#a3aab5","4 5"),
                                     .svgtxt(X(threshold)+5,T+11,paste0("AE = ",threshold),10,"#5c687a"))
  xx<-seq(0,xmax,length.out=260)
  for(i in seq_along(a)) {
    yy<-stats::pbeta(xx,a[i],b[i])
    out<-c(out,.polyline(X(xx),Y(yy),.palette[i],2.4))
  }
  out<-c(out,.svgtxt((L+W-R)/2,H-11,"Absolute error in source contribution (AE)",12,anchor="middle",weight="bold"),
         .svgtxt(L,19,"Probability AE does not exceed x",12,weight="bold"),"</svg>")
  paste(out,collapse="\n")
}

# =============================================================================
# SCED single-case prediction. Called once per row of the case manifest.
# =============================================================================
predict_sced_case <- function(mixture_file, source_file, rdata_file,
                              output_dir = "SCED_Case_Report",
                              case_id = "Case_01",
                              ae_threshold = 0.20,
                              tracer_d2h = TRACER_D2H,
                              tracer_d18o = TRACER_D18O) {
  tracers <- .validate_tracers(tracer_d2h, tracer_d18o)
  message(sprintf("Configured tracers: d2H = '%s', d18O = '%s'",
                  tracers["d2H"], tracers["d18O"]))
  files <- c(
    mixture = .require_file(mixture_file,"mixture CSV"),
    source = .require_file(source_file,"source CSV"),
    posterior = .require_file(rdata_file,"MixSIAR rjags Rdata")
  )
  if (!grepl("\\.[Rr]data$", files["posterior"],ignore.case=TRUE))
    stop("Posterior must be the saved MixSIAR Rjags .Rdata file")
  if (!is.finite(ae_threshold) || ae_threshold <= 0 || ae_threshold >= 1)
    stop("ae_threshold must be strictly between 0 and 1")

  mix <- .read_mixture(files["mixture"], tracers)
  src_result <- .read_sources(files["source"], tracers)
  src <- src_result$data
  posterior <- .get_rjags(files["posterior"],nrow(src))
  geo <- .geometry(src,mix)
  p <- posterior$post
  K <- nrow(src)

  inputs <- data.frame(
    NOS = rep(K,K),
    W95 = p$hi-p$lo,
    log10_S_nearest = log10(geo$nearest),
    Rpos_obs = rep(geo$Rpos,K),
    DFM_scaled = geo$DFM
  )
  if (any(inputs$W95 <= 0 | inputs$W95 >= 1))
    stop("Posterior W95 must be in (0,1). Check convergence or a degenerate posterior")
  z <- list(
    z_near=(inputs$log10_S_nearest-.SCED_MODEL$center["z_near"])/.SCED_MODEL$scale["z_near"],
    z_rpos=(inputs$Rpos_obs-.SCED_MODEL$center["z_rpos"])/.SCED_MODEL$scale["z_rpos"],
    z_dfm=(inputs$DFM_scaled-.SCED_MODEL$center["z_dfm"])/.SCED_MODEL$scale["z_dfm"],
    z_w=(inputs$W95-.SCED_MODEL$center["z_w"])/.SCED_MODEL$scale["z_w"]
  )
  bc <- .SCED_MODEL$coef
  eta <- rep(unname(bc["(Intercept)"]),K)
  if (K > 3) eta <- eta + unname(bc[paste0("NOS",K)])
  for(nm in c("z_w","z_near","z_rpos","z_dfm"))
    eta <- eta + unname(bc[nm])*z[[nm]]
  eta <- eta + unname(bc["z_w:z_rpos"])*z$z_w*z$z_rpos
  mu <- stats::plogis(eta)
  alpha <- mu*.SCED_MODEL$precision
  beta <- (1-mu)*.SCED_MODEL$precision
  q90 <- stats::qbeta(.90,alpha,beta)
  q95 <- stats::qbeta(.95,alpha,beta)
  pex <- stats::pbeta(ae_threshold,alpha,beta,lower.tail=FALSE)
  if (any(!is.finite(c(mu,q90,q95,pex)))) stop("Beta prediction returned nonfinite values")

  # The compact scientific output: one row per source, no extraneous audit tabs.
  ans <- data.frame(
    Source = src$source,
    MixSIAR_mean = p$mean,
    MixSIAR_q025 = p$lo,
    MixSIAR_q975 = p$hi,
    W95 = inputs$W95,
    Predicted_mean_AE = mu,
    AE_Q90 = q90,
    AE_Q95 = q95,
    P_AE_gt_threshold = pex,
    AE_threshold = rep(ae_threshold,K),
    stringsAsFactors = FALSE
  )
  warnings <- character()
  if (nrow(mix) != 6L) warnings <- c(warnings,paste0("The model was developed with 6 mixture observations; this case has ",nrow(mix),"."))
  if (any(src$n != 6L)) warnings <- c(warnings,"Source sample sizes differ from nominal n = 6 used in the simulations.")
  if (is.finite(posterior$rhat_max) && posterior$rhat_max > 1.05)
    warnings <- c(warnings,sprintf("MCMC convergence warning: maximum Rhat = %.3f (target <= 1.05).",posterior$rhat_max))
  if (is.finite(posterior$ess_min) && posterior$ess_min < 400)
    warnings <- c(warnings,sprintf("MCMC effective sample size warning: minimum ESS = %.0f (target >= 400).",posterior$ess_min))
  zmax <- max(unlist(lapply(z,abs)))
  if (zmax > 3) warnings <- c(warnings,sprintf("At least one predictor is %.2f SD from its training mean; extrapolation is possible.",zmax))
  if (!is.finite(posterior$rhat_max) || !is.finite(posterior$ess_min))
    warnings <- c(warnings,"The Rjags summary lacks verifiable Rhat or effective sample size for all source proportions.")
  warnings <- c(warnings,
    "Zero isotope discrimination is assumed for this case; the discrimination file is not read or checked.",
    "SCED was evaluated on simulated water mixtures; actual source errors are unknown in this case.")

  # Write TWO primary deliverables and TWO ready-to-reuse vector figures.
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  csv_file <- file.path(output_dir,"SCED_Source_Predictions.csv")
  html_file <- file.path(output_dir,"SCED_Probability_Report.html")
  utils::write.csv(ans,csv_file,row.names=FALSE,fileEncoding="UTF-8")

  colors <- .palette[seq_len(K)]
  legend <- paste(vapply(seq_len(K),function(i) sprintf(
    '<span class="legenditem"><span class="swatch" style="background:%s"></span><strong>%s</strong></span>',
    colors[i],.escape(src$source[i])),character(1)),collapse="\n")

  # Each SVG is vector artwork embedded within ONE self-contained HTML file.
  space_svg <- .plot_space_svg(geo,src)
  cdf_svg <- .plot_cdf_svg(alpha,beta,q95,ae_threshold,src)
  writeLines(space_svg,file.path(output_dir,"Mixing_Space.svg"),useBytes=TRUE)
  writeLines(cdf_svg,file.path(output_dir,"SCED_Probability_Curves.svg"),useBytes=TRUE)

  rows <- vapply(seq_len(K),function(i) {
    riskbar <- sprintf('<span class="risktrack"><span style="width:%.1f%%;background:%s"></span></span>',
                       100*pex[i],colors[i])
    sprintf(paste0(
      '<tr><td><span class="dot" style="background:%s"></span><strong>%s</strong></td>',
      '<td>%s</td><td class="mono">%s</td><td class="mono">%s</td>',
      '<td class="mono"><strong>%s</strong></td>',
      '<td class="mono">%s %s</td></tr>'),
      colors[i],.escape(src$source[i]),
      .pct(p$mean[i]),paste0(.pct(p$lo[i]),"–",.pct(p$hi[i])),
      .pp(inputs$W95[i]),.pp(q90[i]),.pct(pex[i]),riskbar)
  },character(1))

  note_html <- paste0('<li>',.escape(warnings),'</li>',collapse="\n")
  parameter_rows <- paste0('<tr><td>',.escape(names(bc)),
                           '</td><td>',sprintf('%.8f',unname(bc)),'</td></tr>',collapse="\n")
  scale_rows <- paste0('<tr><td>',.escape(names(.SCED_MODEL$center)),
                        '</td><td>',sprintf('%.8f',unname(.SCED_MODEL$center)),
                        '</td><td>',sprintf('%.8f',unname(.SCED_MODEL$scale)),'</td></tr>',collapse="\n")
  status <- if (length(warnings) > 1L) "Check applicability notes" else "Within basic geometry checks"
  html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
'<meta name="viewport" content="width=device-width,initial-scale=1">',
'<title>SCED | Source contribution error report</title>',
'<style>
:root{--ink:#16293e;--muted:#59697c;--border:#dde5ec;--teal:#1d6d87;--paper:#ffffff;--bg:#f4f7fa}
*{box-sizing:border-box}body{font-family:Arial,"Helvetica Neue",sans-serif;background:var(--bg);color:var(--ink);margin:0}
main{max-width:1190px;margin:30px auto;background:#fff;padding:32px 40px 36px;box-shadow:0 14px 55px #213c5520;border-radius:12px}
.kicker{text-transform:uppercase;letter-spacing:1.8px;color:var(--teal);font-size:11px;font-weight:800}
h1{font-size:29px;line-height:1.15;margin:9px 0 8px;letter-spacing:-.6px}h2{font-size:17px;margin:0 0 4px}
p{color:#415368;line-height:1.5;margin:5px 0}.muted{color:var(--muted)}
.header{display:flex;gap:20px;justify-content:space-between;align-items:start;border-bottom:2px solid var(--ink);padding-bottom:20px}
.tag{border:1px solid var(--border);border-radius:20px;padding:7px 13px;font-size:11px;color:#465b71;white-space:nowrap}
.meta{display:flex;gap:12px;flex-wrap:wrap;margin:20px 0 24px}
.metric{background:#f4f7fa;border-radius:9px;padding:13px 17px;flex:1;min-width:175px}
.metric label{display:block;text-transform:uppercase;letter-spacing:.8px;font-size:10px;font-weight:700;color:#607085;margin-bottom:7px}
.metric strong{font-size:21px;font-variant-numeric:tabular-nums}.metric small{font-size:11px;color:#647487}
.plots{display:grid;grid-template-columns:1fr 1fr;gap:20px}.panel{border:1px solid var(--border);border-radius:10px;padding:16px 17px 9px;min-width:0}
.panel svg{width:100%;height:auto;display:block}.sub{font-size:12px;color:#607085;margin-bottom:4px}
.legend{display:flex;gap:8px 17px;flex-wrap:wrap;margin-top:15px;font-size:11px;color:#405168}
.legenditem{white-space:nowrap}.swatch{width:10px;height:10px;border-radius:3px;display:inline-block;margin-right:5px}
.results{margin-top:24px}table{border-collapse:collapse;width:100%;font-size:12px;margin-top:12px}
th{text-align:left;background:#f2f6f9;color:#455b6e;text-transform:uppercase;font-size:10px;letter-spacing:.4px;padding:11px 9px;border-bottom:1px solid #d7e3eb}
td{padding:12px 9px;border-bottom:1px solid #e9edf2}.mono{font-variant-numeric:tabular-nums;white-space:nowrap}
.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:8px}.sourcelabel{color:#536578;font-size:11px}
.risktrack{display:inline-block;width:57px;height:6px;background:#e9eef3;border-radius:6px;margin-left:8px;vertical-align:middle;overflow:hidden}
.risktrack span{display:block;height:100%;border-radius:6px}
.explainer{padding:13px 18px;background:#eff6f8;border-left:3px solid #3e8293;border-radius:0 7px 7px 0;margin-top:17px;font-size:12px;line-height:1.6;color:#345268}
.notes{margin-top:21px;border-top:1px solid var(--border);padding-top:15px}.notes h2{font-size:14px}
.notes ul{font-size:11px;line-height:1.55;color:#637286;margin-top:8px;padding-left:20px}.notes li{margin:3px 0}
details{margin-top:16px;border-top:1px dashed var(--border);font-size:11px;padding-top:12px;color:#516174}summary{cursor:pointer;font-weight:700;color:#315369}
.params{display:grid;grid-template-columns:1fr 1fr;gap:20px}.params table{font-size:10px}.params td{padding:5px}
.footer{font-size:10px;color:#8090a0;margin-top:19px;border-top:1px solid #e5ebef;padding-top:10px;line-height:1.6}
@media(max-width:760px){main{margin:0;padding:20px;border-radius:0}.header{display:block}.tag{display:inline-block;margin-top:8px}.plots,.params{grid-template-columns:1fr}.results{overflow-x:auto}}
@media print{body{background:white}main{box-shadow:none;margin:0;padding:10px 12px;max-width:none}.panel{break-inside:avoid}.notes,details{break-inside:avoid}.metric{background:#f4f7fa!important;-webkit-print-color-adjust:exact;print-color-adjust:exact}td,th{padding:8px 7px}.plots{gap:12px}.panel{padding:9px}}
</style></head><body><main>',
'<div class="header"><div><div class="kicker">',.escape(case_id),' · SCED · isotope mixing assessment</div>',
'<h1>SCED: Error risk in MixSIAR source estimates</h1>',
'<p>Source-specific predictions of the magnitude and probability of absolute estimation error.</p></div>',
'<div class="tag">',.escape(status),'</div></div>',
'<section class="meta">',
'<div class="metric"><label>Sources</label><strong>',K,'</strong><small> / one mixing case</small></div>',
'<div class="metric"><label>Mixture observations</label><strong>',nrow(mix),'</strong><small> paired isotopes</small></div>',
'<div class="metric"><label>90% AE level · range</label><strong>',
    sprintf('%.1f–%.1f',100*min(q90),100*max(q90)),
    '</strong><small> percentage points</small></div>',
'<div class="metric"><label>P(AE &gt; ',sprintf('%.2f',ae_threshold),') · range</label><strong>',
    sprintf('%.1f–%.1f',100*min(pex),100*max(pex)),
    '</strong><small> % across sources</small></div></section>',
'<div class="plots">',
'<section class="panel"><h2>1 · Mixing space</h2><div class="sub">Paired source isotope means and observed mixture mean</div>',
 space_svg,'</section>',
'<section class="panel"><h2>2 · Probability of AE not exceeding x</h2><div class="sub">One conditional Beta distribution for each source</div>',
 cdf_svg,'</section></div>',
'<div class="legend">',legend,'</div>',
'<section class="results"><h2>3 · Results you can use</h2>',
'<p class="sub">The posterior interval is from your MixSIAR run. AE probabilities are predicted by the fixed SCED model.</p>',
'<table><thead><tr><th>Source</th><th>MixSIAR mean</th><th>MixSIAR 95% interval</th>',
'<th>W95</th><th>90% AE value</th><th>P(AE &gt; ',sprintf('%.2f',ae_threshold),')</th></tr></thead><tbody>',
 paste(rows,collapse="\n"),'</tbody></table></section>',
'<div class="explainer"><strong>How to read this:</strong> A source with a 90% AE value of 20 pp has a model-predicted 90% probability that its absolute contribution error does not exceed 20 percentage points; it is <em>not</em> a guaranteed error bound. ',
'P(AE &gt; ',sprintf('%.2f',ae_threshold),') is the predicted probability that the true contribution differs from the MixSIAR posterior mean by more than ',.pp(ae_threshold),'.</div>',
'<section class="notes"><h2>Applicability and limitations</h2><ul>',note_html,'</ul></section>',
'<details><summary>Model provenance, coefficients and standardization</summary>',
'<p>Frozen Beta-logit Model C; constant precision &phi; = ',sprintf('%.7f',.SCED_MODEL$precision),
'. Trained on ',.SCED_MODEL$n_scenarios,' simulation scenarios (',.SCED_MODEL$n_sources,
' source records). Per-scenario total weights were equal. NOS = 3 is the reference.</p>',
'<div class="params"><table><tr><th>Model term</th><th>Coefficient</th></tr>',parameter_rows,'</table>',
'<table><tr><th>Predictor</th><th>Training mean</th><th>Training SD</th></tr>',scale_rows,'</table></div>',
'<p><strong>Predictors:</strong> NOS, W95, log10(Snearest), Rpos, DFM; with W95 &times; Rpos interaction. Source labels follow MixSIAR alphabetical ordering.</p></details>',
'<div class="footer">Files: mixture=',.escape(basename(files["mixture"])),
'; source=',.escape(basename(files["source"])),
'; posterior=',.escape(basename(files["posterior"])),
'; d2H tracer=',.escape(tracers["d2H"]),
'; d18O tracer=',.escape(tracers["d18O"]),
'. The analysis assumes zero isotope discrimination (not independently checked); it uses no observed true contributions and never refits the model. ',
'Predictions are conditional on the simulation-calibrated model and do not constitute independent field validation.</div>',
'</main></body></html>')
  writeLines(html,html_file,useBytes=TRUE)

  message("\nSCED prediction report created successfully")
  message("  Illustrated HTML report: ", normalizePath(html_file))
  message("  Compact numeric results: ", normalizePath(csv_file))
  message("  Reusable vector figures: Mixing_Space.svg and SCED_Probability_Curves.svg")
  message("  Number of sources: ",K,"; mixture replicates: ",nrow(mix),
          "; Beta precision phi: ",sprintf("%.6f",.SCED_MODEL$precision))
  print(ans[,c("Source","MixSIAR_mean","W95","AE_Q90","P_AE_gt_threshold")],
        row.names=FALSE,digits=3)
  invisible(list(results=ans,html=normalizePath(html_file),
                 csv=normalizePath(csv_file),warnings=warnings,
                 case_metrics=list(NOS=K, mixture_n=nrow(mix),
                   max_Rhat=posterior$rhat_max, min_ESS=posterior$ess_min,
                   max_abs_z=zmax, Rpos_obs=geo$Rpos,
                    prediction_risk_max=max(pex),
                   Q90_min=min(q90), Q90_max=max(q90))))
}


# =============================================================================
# SCED batch interface: a single user-prepared case manifest.
# =============================================================================
.read_case_manifest <- function(path) {
  path <- .require_file(path, "case manifest CSV/XLSX")
  ext <- tolower(tools::file_ext(path))
  d <- if (ext == "csv") {
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8-BOM")
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE))
      stop("To read an Excel manifest, install.packages('readxl'), or save the sheet as CSV.")
    as.data.frame(readxl::read_excel(path), stringsAsFactors=FALSE, check.names=FALSE)
  } else {
    stop("Case list must be .csv or .xlsx")
  }
  columns <- c("Case_ID", "Mixture_File", "Source_File", "Rjags_File")
  if (anyDuplicated(.key(names(d)))) stop("Duplicate manifest column names")
  matched <- match(.key(columns), .key(names(d)))
  if (anyNA(matched))
    stop("Case manifest needs the four columns: ", paste(columns, collapse=", "),
         ". Found: ", paste(names(d), collapse=", "))
  d <- d[, matched, drop=FALSE]
  names(d) <- columns
  if (nrow(d) < 1L) stop("Case manifest contains no cases")
  for (col in columns) {
    d[[col]] <- trimws(as.character(d[[col]]))
    if (anyNA(d[[col]]) || any(!nzchar(d[[col]])))
      stop("Blank/NA entry in manifest column '", col, "'; all four fields are required")
  }
  if (anyDuplicated(d$Case_ID)) stop("Case_ID must be unique in the manifest")
  # Keep all names as entered for presentation. Make safe unique folder names.
  slug <- gsub("[^A-Za-z0-9_-]", "_", d$Case_ID)
  slug[!nzchar(slug)] <- "case"
  d$Folder_ID <- make.unique(slug, sep="_")
  base <- dirname(path)
  absolutize <- function(x) {
    is_abs <- grepl("^(/|[A-Za-z]:[/\\\\]|~)", x)
    if (is_abs) path.expand(x) else file.path(base, x)
  }
  for (col in c("Mixture_File", "Source_File", "Rjags_File"))
    d[[col]] <- vapply(d[[col]], absolutize, character(1))
  attr(d, "manifest_path") <- path
  d
}

.html_batch_summary <- function(log, predictions, output_dir, manifest_path, ae_threshold) {
  .row <- function(x) {
    if (x$status == "SUCCESS") {
      url <- paste0("Cases/",x$folder,"/SCED_Probability_Report.html")
      case_link <- sprintf('<a href="%s">%s</a>', .escape(url), .escape(x$case_id))
      q90 <- sprintf("%.1f–%.1f pp",100*x$q90_min,100*x$q90_max)
      pex <- sprintf("%.1f%%",100*x$max_risk)
      sprintf('<tr><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td><span class="ok">Completed</span></td><td class="explain">%s</td></tr>',
              case_link,x$NOS,q90,pex,.escape(x$note))
    } else {
      sprintf('<tr><td>%s</td><td>—</td><td>—</td><td>—</td><td><span class="bad">Failed</span></td><td class="explain">%s</td></tr>',
              .escape(x$case_id),.escape(x$note))
    }
  }
  ok <- log[log$status == "SUCCESS",,drop=FALSE]
  failed <- log[log$status != "SUCCESS",,drop=FALSE]
  # Plot one descriptive value per scenario, not a statistical fit.
  chart <- '<p>No successful scenarios to chart.</p>'
  if (nrow(ok)) {
    idx <- order(ok$max_risk,decreasing=TRUE)
    a <- ok[idx,,drop=FALSE]
    mx <- max(.25,ceiling(100*max(a$max_risk)/5)*.05)
    cols <- character(nrow(a))
    for (i in seq_len(nrow(a))) {
      w <- 100 * a$max_risk[i] / mx
      cols[i] <- sprintf(paste0('<div class="barline"><div class="barlabel" title="%s">%s</div>',
        '<div class="track"><div class="bar" style="width:%.1f%%"></div></div>',
        '<div class="number">%.1f%%</div></div>'),
        .escape(a$case_id[i]),.escape(a$case_id[i]),w,100*a$max_risk[i])
    }
    chart <- paste(cols,collapse="\n")
  }
  summary_table <- paste(vapply(seq_len(nrow(log)),function(i)
    .row(as.list(log[i,,drop=FALSE])),character(1)),collapse="\n")
  notes <- paste0('<p><strong>Model scope:</strong> two water isotopes (d²H and d¹⁸O), 3–7 sources, ungrouped MixSIAR estimates, zero discrimination and no concentration dependence. Zero discrimination is assumed, not independently verified. SCED outputs are conditional error predictions, not field-observed estimation errors.</p>')
  html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    '<title>SCED | Batch prediction summary</title>',
    '<style>',
    'body{margin:0;background:#f2f5f9;color:#172a3e;font-family:Arial,Helvetica,sans-serif}',
    'main{max-width:1120px;margin:28px auto;background:white;padding:30px 36px;border-radius:13px;box-shadow:0 10px 46px #20394d17}',
    'h1{font-size:28px;margin:7px 0 7px}h2{font-size:17px;margin:25px 0 10px}',
    'p{line-height:1.55;color:#57687b}small{color:#738499}.kicker{font-size:11px;letter-spacing:1.7px;font-weight:800;color:#20768a;text-transform:uppercase}',
    '.stats{display:grid;grid-template-columns:repeat(3,1fr);gap:14px;margin:22px 0}',
    '.stat{padding:17px 19px;background:#f3f7f9;border-radius:10px}.stat strong{font-size:28px;display:block}.stat span{font-size:11px;color:#607389}',
    '.panel{border:1px solid #e0e6ed;padding:20px;border-radius:10px}.barline{display:grid;grid-template-columns:200px 1fr 70px;gap:15px;align-items:center;padding:7px 0}',
    '.barlabel{font-size:12px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.track{height:11px;border-radius:12px;background:#ecf1f5;overflow:hidden}',
    '.bar{height:100%;background:#27768d;border-radius:12px}.number{text-align:right;font-size:12px;font-weight:700;font-variant-numeric:tabular-nums}',
    'table{border-collapse:collapse;width:100%;font-size:12px}th{text-align:left;color:#5b6c7d;background:#f1f6f9;text-transform:uppercase;font-size:10px;letter-spacing:.4px}',
    'td,th{padding:11px 9px;border-bottom:1px solid #e9edf2}a{color:#1f7590;text-decoration:none;font-weight:700}a:hover{text-decoration:underline}',
    '.explain{max-width:310px;overflow-wrap:anywhere;color:#69798a;font-size:11px}.ok{color:#24785f;font-weight:700}.bad{color:#a94e48;font-weight:700}',
    '.note{border-top:1px solid #e5eaf0;margin-top:20px;padding-top:11px;font-size:11px}',
    '@media(max-width:750px){main{margin:0;padding:17px}.stats{grid-template-columns:1fr 1fr}.barline{grid-template-columns:100px 1fr 55px;gap:7px}.panel{overflow-x:auto}}',
    '@media print{main{margin:0;box-shadow:none}.panel{break-inside:avoid}}',
    '</style></head><body><main>',
    '<div class="kicker">SCED | Batch source-water error probability assessment</div>',
    '<h1>SCED prediction | Case overview</h1>',
    '<p>Each case uses its own MixSIAR posterior and isotope source–mixture geometry. Predictions use the same pre-trained SCED Beta model, without refitting.</p>',
    '<div class="stats">',
    sprintf('<div class="stat"><strong>%d</strong><span>Cases supplied</span></div>',nrow(log)),
    sprintf('<div class="stat"><strong>%d</strong><span>Cases completed</span></div>',nrow(ok)),
    sprintf('<div class="stat"><strong>%d</strong><span>Sources predicted</span></div>',nrow(predictions)),
    '</div>',
    sprintf('<h2>Largest source-wise P(AE &gt; %.2f) in each case</h2>',ae_threshold),
    '<div class="panel">',chart,
    '<p style="font-size:11px">Each bar is the maximum across sources within that case; differences do not establish actual estimation errors.</p></div>',
    '<h2>Case results and report links</h2>',
    '<table><thead><tr><th>Case</th><th>Sources</th><th>Q90(AE) range</th>',
    sprintf('<th>Maximum P(AE &gt; %.2f)</th>',ae_threshold),
    '<th>Status</th><th>Notes</th></tr></thead><tbody>',summary_table,'</tbody></table>',
    '<div class="note">',notes,
    '<p>Inputs listed in ',.escape(basename(manifest_path)),
    '. Select a completed case above to open its complete illustrated report. ',
    'Full source-wise predictions are in <strong>SCED_All_Source_Predictions.csv</strong>.</p></div>',
    '</main></body></html>')
  f <- file.path(output_dir, "SCED_Batch_Report.html")
  writeLines(html,f,useBytes=TRUE)
  normalizePath(f)
}

predict_sced_batch <- function(case_list_file=CASE_LIST_FILE,
                                    output_dir=OUTPUT_DIR,
                                    tracer_d2h=TRACER_D2H,
                                    tracer_d18o=TRACER_D18O,
                                    ae_threshold=AE_THRESHOLD) {
  .validate_tracers(tracer_d2h,tracer_d18o)
  if (!is.numeric(ae_threshold) || length(ae_threshold) != 1L ||
      !is.finite(ae_threshold) || ae_threshold <= 0 || ae_threshold >= 1)
    stop("AE_THRESHOLD must be a single value between 0 and 1")
  cases <- .read_case_manifest(case_list_file)
  dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(output_dir,"Cases"), recursive=TRUE, showWarnings=FALSE)
  results <- list()
  logs <- vector("list", nrow(cases))
  for (i in seq_len(nrow(cases))) {
    case_id <- cases$Case_ID[i]
    folder <- cases$Folder_ID[i]
    message(sprintf("\n[%d/%d] Case: %s",i,nrow(cases),case_id))
    case_dir <- file.path(output_dir,"Cases",folder)
    result <- tryCatch(
      predict_sced_case(
        mixture_file=cases$Mixture_File[i],
        source_file=cases$Source_File[i],
        rdata_file=cases$Rjags_File[i],
        output_dir=case_dir,case_id=case_id,
        ae_threshold=ae_threshold,
        tracer_d2h=tracer_d2h,
        tracer_d18o=tracer_d18o
      ),error=function(e)e
    )
    if (inherits(result,"error")) {
      msg <- conditionMessage(result)
      message("  FAILED: ",msg)
      logs[[i]] <- data.frame(case_id=case_id,folder=folder,status="FAILED",
        NOS=NA_integer_,q90_min=NA_real_,q90_max=NA_real_,
        max_risk=NA_real_,note=msg,stringsAsFactors=FALSE)
      next
    }
    x <- result$results
    x <- data.frame(Case_ID=rep(case_id,nrow(x)),x,check.names=FALSE,stringsAsFactors=FALSE)
    results[[length(results)+1L]] <- x
    metric <- result$case_metrics
    note <- if (length(result$warnings) > 2L) "Review applicability notes" else "See case report"
    logs[[i]] <- data.frame(case_id=case_id,folder=folder,status="SUCCESS",
      NOS=metric$NOS,q90_min=metric$Q90_min,q90_max=metric$Q90_max,
      max_risk=metric$prediction_risk_max,note=note,stringsAsFactors=FALSE)
    # One consolidated machine-readable table is enough; keep plots and HTML.
    if (file.exists(result$csv)) unlink(result$csv)
  }
  log <- do.call(rbind,logs)
  combined <- if (length(results)) do.call(rbind,results) else {
    data.frame(Case_ID=character(),Source=character(),
      MixSIAR_mean=numeric(),MixSIAR_q025=numeric(),MixSIAR_q975=numeric(),
      W95=numeric(),Predicted_mean_AE=numeric(),AE_Q90=numeric(),AE_Q95=numeric(),
      P_AE_gt_threshold=numeric(),AE_threshold=numeric(),stringsAsFactors=FALSE)
  }
  combined_file <- file.path(output_dir,"SCED_All_Source_Predictions.csv")
  utils::write.csv(combined,combined_file,row.names=FALSE,fileEncoding="UTF-8")
  report <- .html_batch_summary(log,combined,output_dir,
    attr(cases,"manifest_path"),ae_threshold)
  message(sprintf("\nFinished: %d/%d cases successful; %d source predictions.",
                  sum(log$status=="SUCCESS"),nrow(log),nrow(combined)))
  message("Batch overview: ",report)
  message("Unified source results: ",normalizePath(combined_file))
  invisible(list(predictions=combined,case_log=log,batch_report=report,
                 source_csv=normalizePath(combined_file)))
}

if (isTRUE(AUTO_RUN)) {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) > 2L) stop("Usage: Rscript tool.R Case_List.csv [output_dir]")
  if (length(args) >= 1L) CASE_LIST_FILE <- args[1L]
  if (length(args) == 2L) OUTPUT_DIR <- args[2L]
  predict_sced_batch(
    case_list_file=CASE_LIST_FILE,
    output_dir=OUTPUT_DIR,
    tracer_d2h=TRACER_D2H,
    tracer_d18o=TRACER_D18O,
    ae_threshold=AE_THRESHOLD
  )
}
