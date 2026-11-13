# ============================================================================
# SCED | Reproducible isotope mixing-scenario generation and selection
# Integrated public-release workflow (original scientific methods preserved)
#
# PURPOSE
#   Step 1. Generate dual-isotope mixing scenarios with known source proportions.
#   Step 2. Audit geometry and joint-factor coverage of valid scenarios.
#   Step 3. Select 1,000 scenarios for subsequent MixSIAR fitting.
#
# IMPORTANT SCOPE
#   This script reproduces scenario GENERATION and SELECTION only.
#   It does not fit MixSIAR, examine MCMC convergence, generate the 998 retained
#   scenarios' posterior estimates, or fit the final SCED prediction model.
#   All three original computations run sequentially IN THIS FILE. No source().
#
# REQUIRED DIRECTORY STRUCTURE
#   03_Simulation/
#     SCED_generate_and_select_scenarios.R
#     inputs/
#       contexts.csv
#       pooled_fits.csv
#       residuals.csv
#       source_sd_pairs.csv
#       single_tree_sd_pairs.csv
#
# RUN
#   Rscript SCED_generate_and_select_scenarios.R
#
#   Optional generator settings (default 1000 per study x NOS, seed 20261023):
#   Rscript SCED_generate_and_select_scenarios.R 1000 20261023
#
# RUN FROM ANY WORKING DIRECTORY
#   Output directories are created beside THIS script, not in the current shell
#   directory. Running with Rscript is preferred for reproducible path handling.
#
# REPRODUCIBILITY
#   Original generation V3.2 (seed 20261023), geometry audit V1.0, and
#   selection V1.0 (seed 20261024) are combined without changing their numeric
#   methods, distance definitions, random sampling rules or NOS quotas.
#   Ensure the published inputs and software versions match those of the study.
# ============================================================================

# ---- User-configurable settings (single configuration location) ------------
# Scenarios: 3 empirical studies x 5 NOS levels (3:7) x 1000 = 15,000.
N_PER_STUDY_NOS <- 1000L
GLOBAL_SEED <- 20261023L
N_MIXTURE_REPEATS <- 6L

# Source and xylem-water isotope scenario generation (V3.2).
GEN_OUTPUT_FOLDER <- "outputs_v3_2"
MAX_SOURCE_REDRAWS <- 25L
MAX_PROPOSALS_PER_SCENARIO <- 200000L
CALIB_TOL_NORM <- 0.0005
CALIB_ORDERS <- c(16L, 28L, 40L, 56L, 80L)
CALIB_VERIFY_ORDERS <- c(28L, 40L, 56L, 80L, 112L)
CALIB_MAXIT <- 350L
CALIB_NEWTON_MAXIT <- 100L

# Prescreen geometry audit: descriptive measures, no candidate exclusion.
COVERAGE_OUTPUT_FOLDER <- "coverage_diagnostics"
EXCEL_FALLBACK <- "isotope_candidate_scenarios_V3_2.xlsx"
DANGER_RADIAL_THRESHOLD <- 0.50
POS_BIN_EDGES <- c(1/3, 2/3)
COMMON_Q <- c(0.05, 0.95)

# Fixed, outcome-blind scenario selection; three NOS-balanced subset sizes.
SELECTION_OUTPUT_FOLDER <- "selection_final_V1"
SELECTION_SEED <- 20261024L
PRIMARY_N <- 1000L
COMPARE_N <- c(850L, 1000L, 1150L)
RANDOM_BACKBONE_PER_NOS <- 140L
MIN_PER_NONEMPTY_GRID_CELL <- 3L

# Cross-stage consistency; do not modify independently.
EXPECTED_REPLICATES <- N_MIXTURE_REPEATS
EXPECTED_MIX_REPS <- N_MIXTURE_REPEATS
COVERAGE_FILE <- file.path(COVERAGE_OUTPUT_FOLDER,
                           "scenario_factor_coverage_V3_2.xlsx")
RAW_FOLDER <- GEN_OUTPUT_FOLDER
RAW_EXCEL <- EXCEL_FALLBACK

# Check all packages BEFORE a computationally intensive simulation begins.
for (pkg in c("writexl", "readxl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing R package: ", pkg,
         ". Install once with install.packages('", pkg, "').")
  }
}

stage_message <- function(step, title) {
  message("\n", strrep("=", 74), "\n",
          "SCED workflow | Step ", step, "/3 | ", title, "\n",
          strrep("=", 74))
}

# ---- STEP 1 ---------------------------------------------------------------
# SCIENTIFIC OBJECTIVE
# Generate physically consistent isotope mixtures from known source proportions.
# Source means/residuals come from the empirical-parameter CSVs. For each
# scenario, a Dirichlet draw supplies the known true contributions; the
# theoretical mixture is the contribution-weighted source isotope mean.
# Six xylem-water observations are drawn from a bivariate normal conditional
# on remaining inside the source-isotope convex hull. The conditional mean is
# numerically calibrated to the theoretical mixture, not forced at the sample
# mean level. The original calibration and accept-reject samplers follow.
# OUTPUT: outputs_v3_2/ (candidate scenarios, source values, mixture replicates).
# --------------------------------------------------------------------------
stage_message(1, "Generate candidate isotope mixtures (V3.2)")
OUTPUT_FOLDER <- GEN_OUTPUT_FOLDER
get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly=FALSE), value=TRUE)
  if(length(arg)) return(dirname(normalizePath(sub("^--file=", "", arg[1L]))))
  normalizePath(getwd()) # source() users should setwd() into this script folder
}
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>=1L) N_PER_STUDY_NOS <- as.integer(args[1L])
if(length(args)>=2L) GLOBAL_SEED <- as.integer(args[2L])
if(length(args)>2L) stop("Usage: Rscript SCED_generate_and_select_scenarios.R [n_per_study_nos] [seed]")
if(!is.finite(N_PER_STUDY_NOS) || N_PER_STUDY_NOS<1L ||
   !is.finite(GLOBAL_SEED) || GLOBAL_SEED<0L) stop("Invalid n or global seed")
if(!requireNamespace("writexl",quietly=TRUE))
  stop("Please install Excel writer ONCE: install.packages('writexl')")
if(getRversion() < "3.6.0") stop("R 3.6.0 or later required")
RNGkind(kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
set.seed(GLOBAL_SEED)

ROOT <- get_script_dir()
INPUT <- file.path(ROOT,"inputs")
OUTPUT <- file.path(ROOT,OUTPUT_FOLDER)
if(!dir.exists(INPUT)) stop("Missing inputs/ adjacent to SCED_generate_and_select_scenarios.R")
dir.create(OUTPUT,recursive=TRUE,showWarnings=FALSE)
read_input <- function(name) {
  p<-file.path(INPUT,name)
  if(!file.exists(p)) stop("Missing calibration input: ",p)
  read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8")
}
files<-c("contexts.csv","pooled_fits.csv","residuals.csv",
         "source_sd_pairs.csv","single_tree_sd_pairs.csv")
contexts<-read_input(files[1L]);fits<-read_input(files[2L]);residuals<-read_input(files[3L]);
source_sds<-read_input(files[4L]);trees<-read_input(files[5L])
STUDIES<-c("Rodriguez-Una 2025","WATSON 2025","Wang et al. 2023")
NOS_VALUES<-3:7
N<-length(STUDIES)*length(NOS_VALUES)*N_PER_STUDY_NOS
if(N*N_MIXTURE_REPEATS>1048575L || N*max(NOS_VALUES)>1048575L)
  stop("Excel row limit exceeded; reduce N_PER_STUDY_NOS")
required_ctx<-c("study","context","n_sources","d18O_center","d2H_center",
                "d18O_sd_sources","d18O_min","d18O_max")
if(!all(required_ctx %in% names(contexts))) stop("contexts.csv columns invalid")
if(!all(STUDIES %in% contexts$study) || !all(STUDIES %in% fits$study))
  stop("One of the three expected studies is missing")
if(any(!is.finite(contexts$d18O_sd_sources)) ||
   any(contexts$d18O_max<=contexts$d18O_min)) stop("Invalid background source range")
if(any(!is.finite(fits$slope_within))) stop("Invalid study isotope slopes")
if(any(!is.finite(source_sds$sd18)) || any(!is.finite(source_sds$sd2H)) ||
   any(source_sds$sd18<=0 | source_sds$sd2H<=0)) stop("Invalid source SD values")
if(any(!is.finite(trees$sd18)) || any(!is.finite(trees$sd2H)) ||
   any(trees$sd18<=0 | trees$sd2H<=0) ||
   any(!is.finite(trees$corr)) || any(abs(trees$corr)>=1))
  stop("Invalid single-tree variation/correlation")
ctx_ix<-lapply(STUDIES,function(s) which(contexts$study==s));names(ctx_ix)<-STUDIES
slope<-setNames(fits$slope_within,fits$study)
resid_pool<-lapply(STUDIES,function(s) {
  v<-residuals$deviation_from_pooled_within_line[residuals$study==s]
  v<-v[is.finite(v)]
  if(length(v)<10L) stop("Too few empirical isotope residuals for ",s)
  v-mean(v)
}); names(resid_pool)<-STUDIES
SD_REFERENCE<-c("Rodriguez-Una 2025"="Rodriguez-Una 2025",
                "WATSON 2025"="Rodriguez-Una 2025",
                "Wang et al. 2023"="Wang et al. 2023")
sd_ix<-lapply(STUDIES,function(s)
  which(source_sds$study==unname(SD_REFERENCE[[s]])))
names(sd_ix)<-STUDIES
if(any(lengths(sd_ix)<1L)) stop("Missing paired within-source SDs")

# --------------------------- SOURCE MEANS -------------------------------
draw_oxygen<-function(n,center,sd,minval,maxval) {
  s<-max(sd,(maxval-minval)/20,0.005)
  lo<-pnorm(minval,mean=center,sd=s);hi<-pnorm(maxval,mean=center,sd=s)
  if(!is.finite(lo) || !is.finite(hi) || hi-lo<1e-12)
    return(runif(n,minval,maxval))
  qnorm(runif(n,lo,hi),mean=center,sd=s)
}
make_hull<-function(x,y) {
  dx<-diff(range(x));dy<-diff(range(y))
  if(dx<=1e-12 || dy<=1e-12) return(NULL)
  ux<-(x-min(x))/dx;uy<-(y-min(y))/dy
  h<-chull(ux,uy)
  if(length(h)<3L) return(NULL)
  j<-c(h[-1L],h[1L]);a2<-sum(ux[h]*uy[j]-ux[j]*uy[h])
  if(!is.finite(a2) || abs(a2)<1e-10) return(NULL)
  list(px=ux[h],py=uy[h],orientation=sign(a2),
       xmin=min(x),ymin=min(y),dx=dx,dy=dy,
       normalized_area=abs(a2)/2)
}
inside_hull<-function(xx,yy,h,tol=1e-10) {
  if(is.null(h)) return(rep(FALSE,length(xx)))
  if(length(xx)!=length(yy)) stop("Coordinate lengths differ")
  u<-(xx-h$xmin)/h$dx;v<-(yy-h$ymin)/h$dy
  k<-length(h$px);j<-c(2:k,1L)
  inside<-rep(TRUE,length(xx))
  for(i in seq_len(k)) {
    cross<-(h$px[j[i]]-h$px[i])*(v-h$py[i])-
           (h$py[j[i]]-h$py[i])*(u-h$px[i])
    inside<-inside & (h$orientation*cross>=-tol)
  }
  inside
}

# ---------- DETERMINISTIC POLYGON QUADRATURE AND MEAN MATCHING ----------
# Golub-Welsch Gauss-Legendre rule on [0,1], no external R package.
# Cache the deterministic rules across scenarios (no random-number consumption).
.gauss_cache <- new.env(parent=emptyenv())
gauss_rule<-function(n) {
  key<-as.character(n)
  if(exists(key,envir=.gauss_cache,inherits=FALSE))
    return(get(key,envir=.gauss_cache,inherits=FALSE))
  ii<-seq_len(n-1L)
  off<-ii/sqrt(4*ii^2-1)
  J<-matrix(0,n,n)
  J[cbind(ii,ii+1L)]<-off
  J[cbind(ii+1L,ii)]<-off
  eg<-eigen(J,symmetric=TRUE)
  ord<-order(eg$values)
  result<-list(nodes=(eg$values[ord]+1)/2,weights=eg$vectors[1L,ord]^2)
  assign(key,result,envir=.gauss_cache)
  result
}
# V3.2: triangle fan anchored at the TRUE mixture point (always inside hull).
# Unlike V3.1's fan anchored at a single hull vertex, this distributes
# quadrature nodes near the requested mean, including configurations in which
# the true contribution is close to one source/boundary. It integrates the
# WHOLE convex polygon; it does not clip or shrink the feasible region.
polygon_quadrature<-function(h,order,target) {
  if(is.null(h) || length(target)!=2L || any(!is.finite(target)))
    stop("Invalid polygon or quadrature target")
  rule<-gauss_rule(order)
  tt<-rep(rule$nodes,each=order)
  ss<-rep(rule$nodes,times=order)
  wt<-rep(rule$weights,each=order)*rep(rule$weights,times=order)
  k<-length(h$px)
  if(k<3L) return(NULL)
  pts<-vector("list",k);ww<-vector("list",k)
  for(j in seq_len(k)) {
    nxt<-if(j==k) 1L else j+1L
    ax<-target[1L];ay<-target[2L]
    bx<-h$px[j]-ax;by<-h$py[j]-ay
    cx<-h$px[nxt]-ax;cy<-h$py[nxt]-ay
    determinant<-abs(bx*cy-by*cx)
    if(!is.finite(determinant) || determinant<=0) next
    u<-ax+tt*((1-ss)*bx+ss*cx)
    v<-ay+tt*((1-ss)*by+ss*cy)
    pts[[j]]<-cbind(u,v)
    ww[[j]]<-wt*tt*determinant
  }
  ok<-which(lengths(ww)>0L)
  if(!length(ok)) return(NULL)
  z<-list(points=do.call(rbind,pts[ok]),
          weights=unlist(ww[ok],use.names=FALSE))
  if(any(!is.finite(z$weights)) || any(z$weights<=0) ||
     abs(sum(z$weights)-h$normalized_area) >
       1e-8*max(1,h$normalized_area))
    stop("Quadrature polygon area check failed")
  z
}

# The geometry test is deterministic and never changes the global seed.
smoke_hull <- make_hull(c(-12,-9,-10), c(-85,-72,-76))
if(is.null(smoke_hull)) stop("PRE-FLIGHT FAILED: invalid smoke polygon")
smoke_target<-c(mean(smoke_hull$px),mean(smoke_hull$py))
smoke_quad <- polygon_quadrature(smoke_hull,16L,smoke_target)
if(is.null(smoke_quad) ||
   !all(inside_hull(smoke_hull$xmin+smoke_hull$dx*smoke_quad$points[,1L],
                    smoke_hull$ymin+smoke_hull$dy*smoke_quad$points[,2L],
                    smoke_hull)))
  stop("PRE-FLIGHT FAILED: polygon quadrature points outside polygon")
# A uniform distribution on a triangle has its mean at vertex centroid.
smoke_uniform_mean<-colSums(smoke_quad$points*smoke_quad$weights)/
                    sum(smoke_quad$weights)
if(max(abs(smoke_uniform_mean-smoke_target))>1e-10)
  stop("PRE-FLIGHT FAILED: polygon quadrature first moment incorrect")
message("V3.2 pre-flight target-fan quadrature check: PASS")
# Natural parameter beta = Sigma_u^{-1} (theta_u - true_u).
# f(u|H) proportional to exp(beta' (u-true_u) - .5 (u-true_u)' Q (u-true_u)).
# Solve weighted conditional E[u] = true_u via convex exponential-family dual.
weighted_conditional<-function(beta,quad,target,Q) {
  delta<-sweep(quad$points,2L,target,"-")
  quadterm<-rowSums((delta%*%Q)*delta)
  logw<-log(quad$weights)+as.vector(delta%*%beta)-0.5*quadterm
  maxw<-max(logw)
  w<-exp(logw-maxw);w<-w/sum(w)
  list(mean=target+colSums(delta*w),
       residual=colSums(delta*w),
       logZ=maxw+log(sum(exp(logw-maxw))),
       conditional_sd=sqrt(pmax(0,colSums(sweep(delta,2L,colSums(delta*w),"-")^2*w))))
}
# Numerically minimize the convex dual log-partition function:
#       F(beta)=log Integral_H exp(beta'(u-target)-quadratic(u)) du
# with gradient E_beta[u-target] and Hessian Cov_beta(u-target).
# First use exact-gradient BFGS, then damped Newton if BFGS stalls.
# Independent, higher-order quadrature checks the final expectation.
fit_calibration<-function(h,mu18,mu2H,sd18,sd2H,rho) {
  target<-c((mu18-h$xmin)/h$dx,(mu2H-h$ymin)/h$dy)
  if(!inside_hull(mu18,mu2H,h,tol=1e-8))
    return(list(success=FALSE,reason="target_outside_hull"))
  covU<-matrix(c((sd18/h$dx)^2,rho*sd18*sd2H/(h$dx*h$dy),
                 rho*sd18*sd2H/(h$dx*h$dy),(sd2H/h$dy)^2),2L,2L)
  if(!all(is.finite(covU)) || det(covU)<=1e-14)
    return(list(success=FALSE,reason="invalid_scaled_covariance"))
  Q<-solve(covU)
  beta_start<-c(0,0)
  best_error<-Inf
  for(attempt in seq_along(CALIB_ORDERS)) {
    qfit<-polygon_quadrature(h,CALIB_ORDERS[attempt],target)
    qverify<-polygon_quadrature(h,CALIB_VERIFY_ORDERS[attempt],target)
    if(is.null(qfit) || is.null(qverify))
      return(list(success=FALSE,reason="polygon_quadrature_failed"))
    delta<-sweep(qfit$points,2L,target,"-")
    base<-log(qfit$weights)-0.5*rowSums((delta%*%Q)*delta)
    cache<-new.env(parent=emptyenv())
    eval_beta<-function(beta) {
      beta<-as.numeric(beta)
      if(exists("par",envir=cache,inherits=FALSE) &&
         identical(cache$par,beta)) return(cache$value)
      lw<-base+as.vector(delta%*%beta)
      mx<-max(lw)
      if(!is.finite(mx)) stop("Nonfinite calibration objective")
      ww<-exp(lw-mx);ww<-ww/sum(ww)
      moment<-colSums(delta*ww)
      centered<-sweep(delta,2L,moment,"-")
      hessian<-crossprod(centered,centered*ww)
      value<-list(f=mx+log(sum(exp(lw-mx))),grad=moment,hess=hessian)
      cache$par<-beta;cache$value<-value
      value
    }
    objective<-function(b) eval_beta(b)$f
    gradient<-function(b) eval_beta(b)$grad
    starts<-list(beta_start)
    if(sum(abs(beta_start))>0) starts[[2L]]<-c(0,0)
    best_opt<-NULL
    for(initial in starts) {
      opt<-tryCatch(optim(initial,objective,gradient,method="BFGS",
                          control=list(maxit=CALIB_MAXIT,reltol=1e-12)),
                    error=function(e) NULL)
      if(is.null(opt) || any(!is.finite(opt$par))) next
      g<-gradient(opt$par)
      if(any(!is.finite(g))) next
      if(is.null(best_opt) || max(abs(g))<best_opt$fit_error)
        best_opt<-list(beta=opt$par,fit_error=max(abs(g)),
                       convergence=opt$convergence,
                       n_optim=opt$counts[["function"]],method="BFGS")
      if(max(abs(g))<=CALIB_TOL_NORM/5) break
    }
    if(is.null(best_opt)) {
      best_opt<-list(beta=beta_start,fit_error=Inf,
                     convergence=99L,n_optim=0L,method="BFGS_unavailable")
    }
    # Damped Newton: inverse Hessian is regularized if the support becomes
    # highly concentrated near a vertex. A descent line search is mandatory.
    if(best_opt$fit_error>CALIB_TOL_NORM/5) {
      beta<-best_opt$beta;iter_used<-0L
      for(iter in seq_len(CALIB_NEWTON_MAXIT)) {
        cur<-eval_beta(beta)
        if(max(abs(cur$grad))<=CALIB_TOL_NORM/10) break
        damping<-1e-10
        step<-NULL
        for(j in seq_len(12L)) {
          step<-tryCatch(as.vector(solve(cur$hess+diag(damping,2L),-cur$grad)),
                         error=function(e) NULL)
          if(!is.null(step) && all(is.finite(step)) &&
             sum(step*cur$grad)<0) break
          step<-NULL;damping<-damping*10
        }
        if(is.null(step)) break
        accepted<-FALSE
        multiplier<-1
        for(j in seq_len(30L)) {
          trial<-beta+multiplier*step
          trial_value<-tryCatch(objective(trial),error=function(e) Inf)
          if(is.finite(trial_value) &&
             trial_value<=cur$f+1e-4*multiplier*sum(step*cur$grad)) {
            beta<-trial;accepted<-TRUE;break
          }
          multiplier<-multiplier/2
        }
        if(!accepted) break
        iter_used<-iter
      }
      newerr<-max(abs(gradient(beta)))
      if(is.finite(newerr) && newerr<best_opt$fit_error) {
        best_opt<-list(beta=beta,fit_error=newerr,
                       convergence=if(newerr<=CALIB_TOL_NORM/5) 0L else 1L,
                       n_optim=best_opt$n_optim+iter_used,
                       method="BFGS+Newton")
      }
    }
    beta<-best_opt$beta
    verification<-weighted_conditional(beta,qverify,target,Q)
    err<-verification$residual
    if(all(is.finite(err))) best_error<-min(best_error,max(abs(err)))
    if(all(is.finite(err)) && max(abs(err))<=CALIB_TOL_NORM) {
      theta_u<-target+as.vector(covU%*%beta)
      return(list(success=TRUE,reason="OK",beta=beta,
                  theta_u=theta_u,covU=covU,Q=Q,
                  expected_u=verification$mean,
                  sd_u=verification$conditional_sd,
                  max_error_norm=max(abs(err)),
                  quad_fit=CALIB_ORDERS[attempt],
                  quad_verify=CALIB_VERIFY_ORDERS[attempt],
                  convergence=best_opt$convergence,
                  n_optim=best_opt$n_optim,
                  solver=best_opt$method,
                  verification_quad=qverify))
    }
    # Keep a finite, not catastrophically divergent warm start for refinement.
    # A poor initial quadrature can drive beta towards infinity if the true
    # target lies outside its finite-node convex hull; restart in that case.
    if(all(is.finite(err)) && max(abs(err))<0.05 &&
       all(is.finite(beta)) && max(abs(beta))<1e6)
      beta_start<-beta
    else beta_start<-c(0,0)
  }
  list(success=FALSE,
       reason=paste0("mean_calibration_nonconvergence;best_error_norm=",
                     format(best_error,digits=5)))
}

# Pre-flight: exercise the COMPLETE mean-calibration solver on a broad,
# deterministic three-source polygon before generating any random scenarios.
# This catches numerical/typographical errors that pure quadrature checks miss.
.pre_x<-c(-12,-7,-9)
.pre_y<-c(-88,-48,-76)
.pre_p<-c(0.3,0.4,0.3)
.pre_h<-make_hull(.pre_x,.pre_y)
.pre_cal<-fit_calibration(.pre_h,sum(.pre_x*.pre_p),sum(.pre_y*.pre_p),
                          sd18=0.5,sd2H=2.2,rho=0.7)
if(!isTRUE(.pre_cal$success) ||
   !is.finite(.pre_cal$max_error_norm) ||
   .pre_cal$max_error_norm>CALIB_TOL_NORM)
  stop("PRE-FLIGHT FAILED: mean calibration did not meet its tolerance: ",
       .pre_cal$reason,". No scenarios generated.")
message("V3.2 pre-flight moment calibration check: PASS")
rm(.pre_x,.pre_y,.pre_p,.pre_h,.pre_cal)

# ----------------- REJECTION SAMPLING ON FULL POLYGON -------------------
# Uniform triangle sampling, exactly uniform over the entire hull.
uniform_hull<-function(h,n) {
  k<-length(h$px);B<-2:(k-1L);C<-B+1L
  area<-abs((h$px[B]-h$px[1L])*(h$py[C]-h$py[1L])-
            (h$py[B]-h$py[1L])*(h$px[C]-h$px[1L]))/2
  ix<-sample.int(length(area),n,replace=TRUE,prob=area)
  bb<-B[ix];cc<-C[ix]
  t<-sqrt(runif(n));s<-runif(n)
  cbind(h$px[1L]+t*((1-s)*(h$px[bb]-h$px[1L])+s*(h$px[cc]-h$px[1L])),
        h$py[1L]+t*((1-s)*(h$py[bb]-h$py[1L])+s*(h$py[cc]-h$py[1L])))
}
# Maximum of Gaussian kernel over polygon: Mahalanobis projection of theta
# to polygon edge if theta is outside; mode is theta itself if within hull.
mode_log_kernel<-function(h,theta,Q) {
  theta_raw<-c(h$xmin+h$dx*theta[1L],h$ymin+h$dy*theta[2L])
  if(inside_hull(theta_raw[1L],theta_raw[2L],h)) return(0)
  k<-length(h$px);ii<-c(2:k,1L)
  score<-rep(Inf,k)
  for(j in seq_len(k)) {
    a<-c(h$px[j],h$py[j]);b<-c(h$px[ii[j]],h$py[ii[j]])
    d<-b-a
    den<-as.numeric(crossprod(d,Q%*%d))
    if(den<=0) next
    tt<-as.numeric(crossprod(d,Q%*%(theta-a)))/den
    tt<-max(0,min(1,tt))
    v<-a+tt*d-theta
    score[j]<-as.numeric(crossprod(v,Q%*%v))
  }
  -0.5*min(score)
}
log_kernel<-function(u,theta,Q) {
  d<-sweep(u,2L,theta,"-")
  -0.5*rowSums((d%*%Q)*d)
}
calibrated_sampler<-function(n,h,cal,max_proposals) {
  values<-matrix(NA_real_,n,2L)
  if(!isTRUE(cal$success))
    return(list(values=values,status="calibration_failed",accepted=0L,
                proposed=0L,acceptance=NA_real_,proposal_mode="none"))
  theta<-cal$theta_u;Q<-cal$Q;covU<-cal$covU
  q<-cal$verification_quad
  logf<-log_kernel(q$points,theta,Q)
  m<-max(log(q$weights)+logf)
  logF<-m+log(sum(exp(log(q$weights)+logf-m)))
  logPnormal<-logF-log(2*pi)-0.5*log(det(covU))
  logmax<-mode_log_kernel(h,theta,Q)
  logPuniform<-logF-log(h$normalized_area)-logmax
  if(!is.finite(logF) || !is.finite(logmax) || !is.finite(logPuniform))
    stop("Nonfinite Gaussian envelope normalizer")
  strategy<-if(is.finite(logPnormal) && logPnormal>=logPuniform)
    "normal_rejection" else "polygon_rejection"
  filled<-0L;proposed<-0L;eligible_total<-0L
  while(filled<n && proposed<max_proposals) {
    size<-min(256L,max_proposals-proposed)
    if(strategy=="normal_rejection") {
      z1<-rnorm(size);z2<-rnorm(size)
      sx<-sqrt(covU[1L,1L]);sy<-sqrt(covU[2L,2L])
      rho<-covU[1L,2L]/(sx*sy)
      uu<-cbind(theta[1L]+sx*z1,
                theta[2L]+sy*(rho*z1+sqrt(max(0,1-rho^2))*z2))
      xx<-h$xmin+h$dx*uu[,1L];yy<-h$ymin+h$dy*uu[,2L]
      yes<-which(inside_hull(xx,yy,h))
    } else {
      uu<-uniform_hull(h,size)
      lw<-log_kernel(uu,theta,Q)-logmax
      # The upper bound is calculated analytically. The tiny tolerance
      # only protects floating-point roundoff; it does not change the law.
      if(any(lw>1e-7)) stop("CRITICAL: Gaussian envelope is not a valid upper bound")
      yes<-which(log(runif(size))<pmin(0,lw))
    }
    proposed<-proposed+size;eligible_total<-eligible_total+length(yes)
    if(length(yes)) {
      take<-min(n-filled,length(yes))
      index<-seq.int(filled+1L,filled+take)
      values[index,]<-uu[yes[seq_len(take)],,drop=FALSE]
      filled<-filled+take
    }
  }
  # All accepted values are mapped to original per-mille units.
  values[,1L]<-h$xmin+h$dx*values[,1L]
  values[,2L]<-h$ymin+h$dy*values[,2L]
  list(values=values,accepted=filled,proposed=proposed,
       acceptance=if(proposed>0L) eligible_total/proposed else NA_real_,
       proposal_mode=strategy,
       status=if(filled==n) "complete" else "proposal_limit")
}

separation_stats<-function(x,y,sx,sy) {
  n<-length(x);z<-numeric(n*(n-1L)/2L);k<-0L
  for(i in seq_len(n-1L)) for(j in seq.int(i+1L,n)) {
    k<-k+1L
    z[k]<-sqrt(((x[i]-x[j])/sqrt(sx[i]^2+sx[j]^2))^2+
               ((y[i]-y[j])/sqrt(sy[i]^2+sy[j]^2))^2)
  }
  c(min=min(z),median=median(z))
}

# -------------------------- GENERATE SCENARIOS --------------------------
extra_fields<-c("calibrated_expected_d18O","calibrated_expected_d2H",
                "calibration_bias_d18O","calibration_bias_d2H",
                "calibration_max_error_norm","theta_pretrunc_d18O",
                "theta_pretrunc_d2H","quad_order_fit","quad_order_verify",
                "calibration_convergence","calibration_iterations","source_draw_attempts")
main_fields<-c("scenario_id","NOS","context_source_n","source_d18O_span",
               "source_d2H_span","min_separation_scaled","median_separation_scaled",
               "mix_true_d18O","mix_true_d2H","mix_obs_mean_d18O","mix_obs_mean_d2H",
               "mix_obs_sd_d18O","mix_obs_sd_d2H","mean_shift_d18O","mean_shift_d2H",
               "tree_SD_d18O","tree_SD_d2H","tree_corr","hull_area_normalized",
               "mix_true_inside_hull","all_obs_inside_hull","n_obs_generated",
               "n_proposed","acceptance_rate",paste0("p_true_",1:7))
SM<-matrix(NA_real_,N,length(main_fields),dimnames=list(NULL,main_fields))
EX<-matrix(NA_real_,N,length(extra_fields),dimnames=list(NULL,extra_fields))
SRC<-matrix(NA_real_,N*7L,9L)
colnames(SRC)<-c("scenario_id","source_index","d18O_mean","d2H_mean",
                 "d18O_sd","d2H_sd","source_sample_n","p_true","source_line_residual")
MIX<-matrix(NA_real_,N*N_MIXTURE_REPEATS,4L)
colnames(MIX)<-c("scenario_id","replicate_id","d18O","d2H")
study_col<-character(N);context_col<-character(N)
source_sd_col<-character(N);tree_col<-character(N)
status_col<-character(N);cal_status_col<-character(N);sampler_col<-character(N)
solver_col<-character(N)
id<-0L
message("V3.2 start: ",N," cases, GLOBAL_SEED=",GLOBAL_SEED,
        "; calibration via deterministic polygon quadrature.")
for(study in STUDIES) for(n in NOS_VALUES) for(rep_idx in seq_len(N_PER_STUDY_NOS)) {
  id<-id+1L
  choices<-ctx_ix[[study]]
  bg<-contexts[choices[sample.int(length(choices),1L)],,drop=FALSE]
  sd_choices<-sd_ix[[study]]
  ss<-source_sds[sd_choices[sample.int(length(sd_choices),n,replace=TRUE)],,drop=FALSE]
  x<-rep(NA_real_,n);y<-rep(NA_real_,n);e<-rep(NA_real_,n)
  hull<-NULL;source_attempts<-0L
  for(attempt in seq_len(MAX_SOURCE_REDRAWS)) {
    source_attempts<-attempt
    x<-draw_oxygen(n,bg$d18O_center,bg$d18O_sd_sources,bg$d18O_min,bg$d18O_max)
    e<-sample(resid_pool[[study]],n,replace=TRUE)
    y<-bg$d2H_center+unname(slope[[study]])*(x-bg$d18O_center)+e
    hull<-make_hull(x,y)
    if(!is.null(hull)) break
  }
  z<-rgamma(n,shape=1,rate=1);p<-z/sum(z)
  mu18<-sum(p*x);mu2H<-sum(p*y)
  t<-trees[sample.int(nrow(trees),1L),,drop=FALSE]
  if(is.null(hull)) {
    cal<-list(success=FALSE,reason="source_hull_degenerate_after_retries")
  } else {
    cal<-tryCatch(fit_calibration(hull,mu18,mu2H,t$sd18,t$sd2H,t$corr),
                  error=function(err) list(success=FALSE,reason=paste0("calibration_exception: ",conditionMessage(err))))
  }
  # Unexpected coding errors are not legitimate scenario failures: stop immediately.
  if(!isTRUE(cal$success) &&
     grepl("^calibration_exception:", cal$reason))
    stop("INTERNAL CALIBRATION ERROR at scenario ", id,
         " [", study, ", NOS=", n, "]: ", cal$reason,
         ". Fix the script before rerunning; no Excel will be written.")
  if(isTRUE(cal$success)) {
    draw<-tryCatch(calibrated_sampler(N_MIXTURE_REPEATS,hull,cal,MAX_PROPOSALS_PER_SCENARIO),
                   error=function(err) list(values=matrix(NA_real_,N_MIXTURE_REPEATS,2L),
                      accepted=0L,proposed=0L,acceptance=NA_real_,
                      proposal_mode="none",status="sampler_exception",
                      error_message=conditionMessage(err)))
  } else {
    draw<-list(values=matrix(NA_real_,N_MIXTURE_REPEATS,2L),
               accepted=0L,proposed=0L,acceptance=NA_real_,
               proposal_mode="none",status="calibration_failed")
  }
  if(identical(draw$status,"sampler_exception"))
    stop("INTERNAL SAMPLER ERROR at scenario ",id,": ",draw$error_message,
         ". Fix the script before rerunning; no Excel will be written.")
  complete<-identical(draw$status,"complete")
  obs18<-draw$values[,1L];obs2H<-draw$values[,2L]
  m18<-if(complete) mean(obs18) else NA_real_
  m2H<-if(complete) mean(obs2H) else NA_real_
  sample_sd18<-if(complete) sd(obs18) else NA_real_
  sample_sd2H<-if(complete) sd(obs2H) else NA_real_
  sep<-separation_stats(x,y,ss$sd18,ss$sd2H)
  true_inside<-if(is.null(hull)) NA else inside_hull(mu18,mu2H,hull)
  obs_inside<-if(complete) all(inside_hull(obs18,obs2H,hull)) else NA
  SM[id,]<-c(id,n,bg$n_sources,diff(range(x)),diff(range(y)),
             sep[["min"]],sep[["median"]],mu18,mu2H,m18,m2H,
             sample_sd18,sample_sd2H,m18-mu18,m2H-mu2H,
             t$sd18,t$sd2H,t$corr,
             if(is.null(hull)) NA_real_ else hull$normalized_area,
             as.numeric(true_inside),as.numeric(obs_inside),draw$accepted,
             draw$proposed,draw$acceptance,c(p,rep(NA_real_,7L-n)))
  if(isTRUE(cal$success)) {
    EX[id,]<-c(hull$xmin+hull$dx*cal$expected_u[1L],
                hull$ymin+hull$dy*cal$expected_u[2L],
                hull$dx*(cal$expected_u[1L]-(mu18-hull$xmin)/hull$dx),
                hull$dy*(cal$expected_u[2L]-(mu2H-hull$ymin)/hull$dy),
                cal$max_error_norm,
                hull$xmin+hull$dx*cal$theta_u[1L],
                hull$ymin+hull$dy*cal$theta_u[2L],
                cal$quad_fit,cal$quad_verify,cal$convergence,
                cal$n_optim,source_attempts)
  } else {
    EX[id,"source_draw_attempts"]<-source_attempts
  }
  study_col[id]<-study;context_col[id]<-bg$context
  source_sd_col[id]<-unname(SD_REFERENCE[[study]])
  tree_col[id]<-paste(t$species,t$tree,sep="|")
  status_col[id]<-draw$status
  cal_status_col[id]<-if(!is.null(draw$error_message))
    paste(cal$reason,draw$error_message,sep=" | ") else cal$reason
  sampler_col[id]<-draw$proposal_mode
  solver_col[id]<-if(isTRUE(cal$success)) cal$solver else "none"
  SRC[(id-1L)*7L+seq_len(n),]<-cbind(rep(id,n),seq_len(n),x,y,
                                     ss$sd18,ss$sd2H,ss$n,p,e)
  MIX[(id-1L)*N_MIXTURE_REPEATS+seq_len(N_MIXTURE_REPEATS),]<-
    cbind(rep(id,N_MIXTURE_REPEATS),seq_len(N_MIXTURE_REPEATS),obs18,obs2H)
  if(id%%max(1L,min(10L,N_PER_STUDY_NOS))==0L) message("V3.2 progress ",id,"/",N,
                           "; incomplete: ",sum(status_col[seq_len(id)]!="complete"))
}
scenarios<-data.frame(study=study_col,context=context_col,
                      source_SD_reference=source_sd_col,tree_reference=tree_col,
                      generation_status=status_col,
                      calibration_status=cal_status_col,
                      proposal_mode=sampler_col,calibration_solver=solver_col,SM,EX,
                      stringsAsFactors=FALSE,check.names=FALSE)
scenarios$scenario_id<-as.integer(scenarios$scenario_id)
scenarios$NOS<-as.integer(scenarios$NOS)
scenarios$context_source_n<-as.integer(scenarios$context_source_n)
sources<-as.data.frame(SRC[!is.na(SRC[,1L]),,drop=FALSE],check.names=FALSE)
sources$scenario_id<-as.integer(sources$scenario_id)
sources$source_index<-as.integer(sources$source_index)
sources$source_label<-paste0("S",sources$source_index)
sources$study<-study_col[sources$scenario_id]
sources$source_SD_reference<-source_sd_col[sources$scenario_id]
mixtures<-as.data.frame(MIX,check.names=FALSE)
mixtures$scenario_id<-as.integer(mixtures$scenario_id)
mixtures$replicate_id<-as.integer(mixtures$replicate_id)
mixtures$generation_status<-status_col[mixtures$scenario_id]

# ----------------------------- HARD QA ----------------------------------
if(nrow(scenarios)!=N || nrow(mixtures)!=N*N_MIXTURE_REPEATS)
  stop("CRITICAL: generated table row counts do not match config")
if(any(abs(rowsum(sources$p_true,sources$scenario_id)[,1L]-1)>1e-10))
  stop("CRITICAL: true contributions do not sum to 1")
if(any(as.integer(table(sources$scenario_id)) != scenarios$NOS))
  stop("CRITICAL: source count inconsistent with NOS")
if(any(scenarios$mix_true_inside_hull!=1,na.rm=TRUE))
  stop("CRITICAL: theoretical mixture outside convex hull")
if(any(scenarios$all_obs_inside_hull!=1,na.rm=TRUE))
  stop("CRITICAL: accepted observation outside convex hull")
if(any(scenarios$calibration_max_error_norm > CALIB_TOL_NORM,na.rm=TRUE))
  stop("CRITICAL: calibration verification error exceeds tolerance")
verified_mu18<-rowsum(sources$d18O_mean*sources$p_true,sources$scenario_id)[,1L]
verified_mu2H<-rowsum(sources$d2H_mean*sources$p_true,sources$scenario_id)[,1L]
max_mass_error<-max(abs(verified_mu18-scenarios$mix_true_d18O),
                    abs(verified_mu2H-scenarios$mix_true_d2H))
if(!is.finite(max_mass_error) || max_mass_error>1e-9)
  stop("CRITICAL: source contribution mass balance check failed")
valid<-scenarios$generation_status=="complete"
if(!any(valid))
  stop("FATAL QA: zero complete scenarios. Check calibration_status; no output written.")
if(any(!is.finite(scenarios$calibration_max_error_norm[valid])))
  stop("CRITICAL: completed scenario lacks valid calibration diagnostics")
qcs<-do.call(rbind,lapply(STUDIES,function(study) {
  do.call(rbind,lapply(NOS_VALUES,function(n) {
    z<-scenarios[scenarios$study==study & scenarios$NOS==n,,drop=FALSE]
    v<-z$generation_status=="complete"
    data.frame(study=study,NOS=n,n_total=nrow(z),n_complete=sum(v),
               n_incomplete=sum(!v),
               median_mean_error_norm=if(any(v)) median(z$calibration_max_error_norm[v]) else NA_real_,
               max_mean_error_norm=if(any(v)) max(z$calibration_max_error_norm[v]) else NA_real_,
               median_obs_sd_d18O=if(any(v)) median(z$mix_obs_sd_d18O[v]) else NA_real_,
               median_obs_sd_d2H=if(any(v)) median(z$mix_obs_sd_d2H[v]) else NA_real_,
               n_retried_hulls=sum(z$source_draw_attempts>1L,na.rm=TRUE),
               n_normal_rejection=sum(z$proposal_mode=="normal_rejection"),
               n_polygon_rejection=sum(z$proposal_mode=="polygon_rejection"))
  }))
}))
row.names(qcs)<-NULL
n_bad<-sum(!valid)
manifest<-data.frame(item=c("SCRIPT_VERSION","GLOBAL_SEED","R_VERSION","RNG_KIND",
                            "N_PER_STUDY_NOS","N_MIXTURE_REPEATS","N_TOTAL",
                            "N_COMPLETE","N_INCOMPLETE","MAX_SOURCE_REDRAWS",
                            "MAX_PROPOSALS_PER_SCENARIO","CALIB_TOL_NORM",
                            "CALIB_QUAD_ORDERS","CALIB_VERIFY_ORDERS",
                            "CALIB_QUADRATURE_GEOMETRY","CALIB_SOLVER",
                            "MAX_MASS_BALANCE_ERROR","FINAL_SELECTION",
                            "OVERALL_STATUS","SCRIPT_MD5",paste0("INPUT_MD5_",files)),
                     value=c("V3.2",as.character(GLOBAL_SEED),R.version.string,
                             paste(RNGkind(),collapse=";"),
                             as.character(N_PER_STUDY_NOS),as.character(N_MIXTURE_REPEATS),
                             as.character(N),as.character(N-n_bad),as.character(n_bad),
                             as.character(MAX_SOURCE_REDRAWS),
                             as.character(MAX_PROPOSALS_PER_SCENARIO),
                             as.character(CALIB_TOL_NORM),
                             paste(CALIB_ORDERS,collapse=","),
                             paste(CALIB_VERIFY_ORDERS,collapse=","),
                             "true-mixture-anchored full-polygon triangle fan",
                             "BFGS; damped Newton fallback",
                             format(max_mass_error,digits=15),"NO",
                             if(n_bad==0L) "PENDING_SCIENTIFIC_REVIEW" else "INCOMPLETE_NEEDS_REVIEW",
                             unname(tools::md5sum(file.path(ROOT,"SCED_generate_and_select_scenarios.R"))),
                             unname(tools::md5sum(file.path(INPUT,files)))),
                     stringsAsFactors=FALSE)
notes<-data.frame(topic=c("Scope","Source mean model","Source SD","Plant data",
                          "Mean calibration","Independent observations","Hull constraint",
                          "Calibration tolerance","Failure policy","No MixSIAR",
                          "Source data caution"),
                  explanation=c(
                    "3 independent studies x NOS 3..7, each scenario independently generated.",
                    "Study-specific local context and empirical isotope covariance/residuals; no fixed 10-source library.",
                    "Paired observed within-source SD; WATSON uses Rodriguez-Una SD as documented proxy.",
                    "Single-tree Goldsmith branch variation as reference covariance, with a recorded field-label swap assumption.",
                    "Numerically solve the normal center theta so the full-polygon truncated expectation equals the true mixed isotope value; target-anchored quadrature with adaptive orders and Newton fallback.",
                    "Six IID draws conditional on calibrated normal distribution and full convex hull; their SAMPLE mean is NOT forced to truth.",
                    "All accepted observations must lie within the source-MEAN convex hull; this is a study-specific idealization.",
                    "Max coordinate-wise calibrated expectation error divided by source coordinate range; precision independently checked at higher quadrature order.",
                    "Calibration or draw failures retained in outputs with flags; no silent cases removed or alternate statistical models substituted.",
                    "Candidate-generation stage only; final selection is handled by Step 3. No MixSIAR runs.",
                    "All input study values are empirical constraints, not an exhaustive representation of global water isotopes."),
                  stringsAsFactors=FALSE)
write.csv(scenarios,file.path(OUTPUT,"01_candidate_scenarios.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(sources,file.path(OUTPUT,"02_candidate_sources.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(mixtures,file.path(OUTPUT,"03_candidate_mixtures.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
write.csv(qcs,file.path(OUTPUT,"04_QC_by_study_NOS.csv"),row.names=FALSE,na="",fileEncoding="UTF-8")
csv_paths<-file.path(OUTPUT,c("01_candidate_scenarios.csv","02_candidate_sources.csv",
                              "03_candidate_mixtures.csv","04_QC_by_study_NOS.csv"))
hashes<-data.frame(file=basename(csv_paths),md5=unname(tools::md5sum(csv_paths)))
write.csv(hashes,file.path(OUTPUT,"05_output_hashes.csv"),row.names=FALSE)
workbook<-file.path(OUTPUT,"isotope_candidate_scenarios_V3_2.xlsx")
writexl::write_xlsx(list(README=notes,MANIFEST=manifest,
                        QC_by_group=qcs,Candidates=scenarios,
                        Sources=sources,Mixtures=mixtures,Output_hashes=hashes),
                   path=workbook)
message("--------------------------------------------------------")
message("V3.2 R finished. Workbook: ",normalizePath(workbook))
message("Complete ",N-n_bad," / ",N,"; unsuccessful are EXPLICITLY flagged.")
message("Mean calibration validated by higher-order independent quadrature.")
message("Inspect the generated diagnostics before using the selected scenarios in a scientific analysis.")


# ---- STEP 2 ---------------------------------------------------------------
# SCIENTIFIC OBJECTIVE
# Independently characterize the *observed sample-mean* isotope mixing space
# for every valid scenario. Important quantities:
#   min_sep  : smallest pairwise standardized source separation (screening);
#   Rpos_obs : radial position of mean xylem isotope mixture in the source
#              convex hull (0 central reference, 1 boundary);
#   DFM      : source-to-mixture distance (diagnostic at this stage).
# This stage is deterministic, uses NO RNG, and selects NO scenarios.
# OUTPUT: coverage_diagnostics/scenario_factor_coverage_V3_2.xlsx
#         The Scenario_metrics worksheet feeds Step 3.
# --------------------------------------------------------------------------
stage_message(2, "Audit mixing-space geometry and candidate coverage")
INPUT_FOLDER <- GEN_OUTPUT_FOLDER
OUTPUT_FOLDER <- COVERAGE_OUTPUT_FOLDER
if(!requireNamespace("writexl",quietly=TRUE))
  stop("Install Excel writer once: install.packages('writexl')")
ROOT <- get_script_dir()
# Resolve inputs beside this consolidated script; works independently of shell cwd.
input_candidates <- unique(c(file.path(ROOT,INPUT_FOLDER),ROOT))
exists_csv <- function(p) all(file.exists(file.path(p,c(
  "01_candidate_scenarios.csv","02_candidate_sources.csv","03_candidate_mixtures.csv"))))
exists_excel <- function(p) file.exists(file.path(p,EXCEL_FALLBACK))
idx_csv <- which(vapply(input_candidates,exists_csv,logical(1)))
idx_excel <- which(vapply(input_candidates,exists_excel,logical(1)))
if(length(idx_csv)) {
  INPUT <- input_candidates[idx_csv[1L]]
  input_type <- "V3.2_CSV"
  candidates <- read.csv(file.path(INPUT,"01_candidate_scenarios.csv"),
                         stringsAsFactors=FALSE,check.names=FALSE)
  sources <- read.csv(file.path(INPUT,"02_candidate_sources.csv"),
                      stringsAsFactors=FALSE,check.names=FALSE)
  mixtures <- read.csv(file.path(INPUT,"03_candidate_mixtures.csv"),
                       stringsAsFactors=FALSE,check.names=FALSE)
  used_paths <- file.path(INPUT,c("01_candidate_scenarios.csv",
                      "02_candidate_sources.csv","03_candidate_mixtures.csv"))
} else if(length(idx_excel)) {
  if(!requireNamespace("readxl",quietly=TRUE))
    stop("CSV files not found. Either retain the original V3.2 CSVs, or install.packages('readxl')")
  INPUT <- input_candidates[idx_excel[1L]]
  input_type <- "V3.2_Excel"
  p <- file.path(INPUT,EXCEL_FALLBACK)
  candidates <- as.data.frame(readxl::read_excel(p,sheet="Candidates"))
  sources <- as.data.frame(readxl::read_excel(p,sheet="Sources"))
  mixtures <- as.data.frame(readxl::read_excel(p,sheet="Mixtures"))
  used_paths <- p
} else stop("No V3.2 data found. Put this script beside outputs_v3_2/ (containing the 3 CSVs), or put the V3.2 Excel beside it.")

require_fields <- function(x,need,label) {
  miss <- setdiff(need,names(x))
  if(length(miss)) stop(label," missing columns: ",paste(miss,collapse=", "))
}
require_fields(candidates,c("scenario_id","study","NOS","context",
  "generation_status","calibration_status","n_obs_generated",
  "mix_obs_mean_d18O","mix_obs_mean_d2H","mix_true_d18O",
  "mix_true_d2H","min_separation_scaled","median_separation_scaled",
  "calibration_max_error_norm"),"Candidates")
require_fields(sources,c("scenario_id","source_index","d18O_mean","d2H_mean",
  "d18O_sd","d2H_sd","p_true"),"Sources")
require_fields(mixtures,c("scenario_id","replicate_id","d18O","d2H"),"Mixtures")
for(nm in c("scenario_id","NOS","n_obs_generated"))
  candidates[[nm]] <- as.integer(candidates[[nm]])
for(nm in c("scenario_id","source_index")) sources[[nm]]<-as.integer(sources[[nm]])
mixtures$scenario_id<-as.integer(mixtures$scenario_id)
mixtures$replicate_id<-as.integer(mixtures$replicate_id)
if(anyNA(candidates$scenario_id) || anyDuplicated(candidates$scenario_id))
  stop("Missing/duplicate scenario IDs in Candidates")
if(anyNA(sources$scenario_id) || anyDuplicated(sources[,c("scenario_id","source_index")]))
  stop("Invalid source IDs or duplicated scenario/source pairs")
if(!all(sources$scenario_id %in% candidates$scenario_id) ||
   !all(mixtures$scenario_id %in% candidates$scenario_id))
  stop("Unmatched scenario IDs across tables")

complete <- candidates$generation_status=="complete" &
            candidates$calibration_status=="OK" &
            candidates$n_obs_generated==EXPECTED_REPLICATES &
            is.finite(candidates$mix_obs_mean_d18O) &
            is.finite(candidates$mix_obs_mean_d2H)
complete[is.na(complete)] <- FALSE
valid_ids<-candidates$scenario_id[complete]
if(!length(valid_ids)) stop("No complete candidate scenarios")
flagged <- candidates[!complete,c("scenario_id","study","NOS","context",
             "generation_status","calibration_status","n_obs_generated"),drop=FALSE]
if(!nrow(flagged)) flagged <- flagged[FALSE,,drop=FALSE]

# Isotope axes are separately rescaled by scenario source ranges for
# 2D positional / DFC / DFM geometry, avoiding domination by delta2H units.
# Radial position itself is affine-invariant; the rescaling improves arithmetic.
cross2 <- function(a,b) a[1L]*b[2L]-a[2L]*b[1L]
hull_from_sources <- function(x,y) {
  rangex<-range(x);rangey<-range(y)
  dx<-diff(rangex);dy<-diff(rangey)
  if(!is.finite(dx) || !is.finite(dy) || dx<1e-12 || dy<1e-12)
    stop("Source isotope span too small for geometry")
  xy<-cbind((x-rangex[1L])/dx,(y-rangey[1L])/dy)
  ix<-chull(xy[,1L],xy[,2L])
  if(length(ix)<3L) stop("Degenerate source convex hull")
  poly<-xy[ix,,drop=FALSE]
  j<-c(seq.int(2L,nrow(poly)),1L)
  twice_area<-sum(poly[,1L]*poly[j,2L]-poly[j,1L]*poly[,2L])
  if(twice_area<0) {ix<-rev(ix);poly<-xy[ix,,drop=FALSE];twice_area<- -twice_area}
  if(!is.finite(twice_area) || twice_area<1e-12) stop("Near-zero hull area")
  list(poly=poly,xy=xy,indices=ix,center=colMeans(poly),
       dx=dx,dy=dy,xmin=rangex[1L],ymin=rangey[1L],area=twice_area/2)
}
# For CCW polygons, intersection of ray C + t*(P-C) with the hull boundary.
# True/observed mixture P corresponds to t=1; radial score = 1/t_boundary.
radial_position <- function(pt,h,tol=1e-9) {
  P<-c((pt[1L]-h$xmin)/h$dx,(pt[2L]-h$ymin)/h$dy)
  C<-h$center;direction<-P-C
  if(sum(direction^2)<1e-24) return(0)
  poly<-h$poly;k<-nrow(poly)
  hits<-numeric(0)
  for(i in seq_len(k)) {
    j<-if(i==k) 1L else i+1L
    a<-poly[i,];edge<-poly[j,]-a
    left_at_center<-cross2(edge,C-a)
    derivative<-cross2(edge,direction)
    if(derivative < -1e-13)
      hits<-c(hits,left_at_center/(-derivative))
  }
  if(!length(hits)) stop("Could not intersect center-to-mixture ray with hull")
  tmax<-min(hits)
  if(!is.finite(tmax) || tmax<=0) stop("Invalid radial intersection")
  r<-1/tmax
  if(r>1+tol) stop("Mixture lies outside source hull; radial=",r)
  max(0,min(1,r))
}
# Strict study-compatible min/median standardized Euclidean pairwise distance.
sep_dist <- function(x,y,sx,sy) {
  n<-length(x);out<-numeric(n*(n-1L)/2L);k<-0L
  for(i in seq_len(n-1L)) for(j in seq.int(i+1L,n)) {
    k<-k+1L
    out[k]<-sqrt((x[i]-x[j])^2/(sx[i]^2+sx[j]^2)+
                 (y[i]-y[j])^2/(sy[i]^2+sy[j]^2))
  }
  c(min=min(out),median=median(out))
}

SRC<-split(sources,sources$scenario_id)
MIX<-split(mixtures,mixtures$scenario_id)
rows<-vector("list",length(valid_ids))
dfm_rows<-vector("list",length(valid_ids))
message("Coverage diagnostic: ",length(valid_ids)," valid scenarios; no new random draws")
for(z in seq_along(valid_ids)) {
  id<-valid_ids[z]
  a<-candidates[match(id,candidates$scenario_id),,drop=FALSE]
  src<-SRC[[as.character(id)]]
  mix<-MIX[[as.character(id)]]
  if(is.null(src) || nrow(src)!=a$NOS ||
     any(!is.finite(as.matrix(src[,c("d18O_mean","d2H_mean","d18O_sd","d2H_sd","p_true")]))))
    stop("Bad sources for scenario ",id)
  if(is.null(mix) || sum(is.finite(mix$d18O) & is.finite(mix$d2H))!=EXPECTED_REPLICATES)
    stop("Missing mixture replicates for complete scenario ",id)
  x<-src$d18O_mean;y<-src$d2H_mean
  h<-hull_from_sources(x,y)
  pobs<-c(a$mix_obs_mean_d18O,a$mix_obs_mean_d2H)
  ptrue<-c(a$mix_true_d18O,a$mix_true_d2H)
  r_obs<-radial_position(pobs,h)
  r_true<-radial_position(ptrue,h)
  centroid_raw<-c(h$xmin+h$center[1L]*h$dx,
                  h$ymin+h$center[2L]*h$dy)
  center_dif_obs<-pobs-centroid_raw
  center_dif_true<-ptrue-centroid_raw
  box_diagonal<-sqrt(h$dx^2+h$dy^2)
  dfc_original_proxy<-sqrt(sum(center_dif_obs^2))/box_diagonal
  dfc_scaled<-sqrt(sum((center_dif_obs/c(h$dx,h$dy))^2))/sqrt(2)
  dfc_true_scaled<-sqrt(sum((center_dif_true/c(h$dx,h$dy))^2))/sqrt(2)
  dm18<-x-pobs[1L];dm2H<-y-pobs[2L]
  dfm_raw<-sqrt(dm18^2+dm2H^2)
  dfm_scaled<-sqrt((dm18/h$dx)^2+(dm2H/h$dy)^2)/sqrt(2)
  sep<-sep_dist(x,y,src$d18O_sd,src$d2H_sd)
  if(max(abs(sep-c(a$min_separation_scaled,a$median_separation_scaled)))>1e-7)
    stop("Independent source separation check failed at scenario ",id)
  if(abs(sum(src$p_true)-1)>1e-8 ||
     abs(sum(src$p_true*x)-ptrue[1L])>1e-8 ||
     abs(sum(src$p_true*y)-ptrue[2L])>1e-8)
    stop("Source weighted truth check failed at scenario ",id)
  mean18<-mean(mix$d18O[is.finite(mix$d18O)])
  mean2H<-mean(mix$d2H[is.finite(mix$d2H)])
  if(abs(mean18-pobs[1L])>1e-8 || abs(mean2H-pobs[2L])>1e-8)
    stop("Mixture mean / Candidates mismatch for scenario ",id)
  rows[[z]]<-data.frame(
    scenario_id=id,study=a$study,context=a$context,NOS=a$NOS,
    min_sep=sep["min"],median_sep=sep["median"],
    Rpos_obs=r_obs,Rpos_true=r_true,
    Rpos_observed_minus_true=r_obs-r_true,
    Fry50_center_proxy_obs=(r_obs<=DANGER_RADIAL_THRESHOLD),
    Fry50_center_proxy_true=(r_true<=DANGER_RADIAL_THRESHOLD),
    DFC_original_proxy=dfc_original_proxy,
    DFC_scaled=dfc_scaled,DFC_scaled_true=dfc_true_scaled,
    DFM_min_scaled=min(dfm_scaled),DFM_median_scaled=median(dfm_scaled),
    DFM_max_scaled=max(dfm_scaled),
    closest_source_index=src$source_index[which.min(dfm_scaled)],
    hull_area_normalized=h$area,hull_n_vertices=nrow(h$poly),
    n_nonvertex_sources=a$NOS-nrow(h$poly),
    source_span_d18O=h$dx,source_span_d2H=h$dy,
    sample_mean_shift_d18O=pobs[1L]-ptrue[1L],
    sample_mean_shift_d2H=pobs[2L]-ptrue[2L],
    calibration_max_error_norm=a$calibration_max_error_norm,
    stringsAsFactors=FALSE)
  dfm_rows[[z]]<-data.frame(
    scenario_id=id,study=a$study,NOS=a$NOS,
    source_index=src$source_index,
    p_true=src$p_true,
    source_d18O=x,source_d2H=y,
    source_sd_d18O=src$d18O_sd,source_sd_d2H=src$d2H_sd,
    DFM_raw_permil=dfm_raw,DFM_scaled=dfm_scaled,
    hull_vertex=(seq_len(nrow(src)) %in% h$indices),
    stringsAsFactors=FALSE)
  if(z%%250L==0L) message("  processed ",z," / ",length(valid_ids))
}
MET<-do.call(rbind,rows)
DFM<-do.call(rbind,dfm_rows)
row.names(MET)<-NULL;row.names(DFM)<-NULL
if(nrow(MET)!=length(valid_ids) || nrow(DFM)!=sum(MET$NOS))
  stop("Output row counts mismatch")

# Within each study × NOS, separator terciles are used ONLY to diagnose
# feasible coverage; final association models retain continuous min_sep.
MET$sep_tercile<-NA_character_
MET$pos_bin<-ifelse(MET$Rpos_obs<=POS_BIN_EDGES[1L],"center_0_1_3",
             ifelse(MET$Rpos_obs<=POS_BIN_EDGES[2L],"middle_1_3_2_3",
                                                       "outer_2_3_1"))
group_keys<-unique(MET[,c("study","NOS")])
group_keys<-group_keys[order(group_keys$study,group_keys$NOS),,drop=FALSE]
sep_thresholds<-vector("list",nrow(group_keys))
for(k in seq_len(nrow(group_keys))) {
  idx<-which(MET$study==group_keys$study[k] & MET$NOS==group_keys$NOS[k])
  bounds<-as.numeric(quantile(MET$min_sep[idx],c(1/3,2/3),na.rm=TRUE,names=FALSE))
  if(bounds[1]>=bounds[2]) {
    warning("Separation quantile ties in group ",group_keys$study[k]," NOS=",group_keys$NOS[k])
  }
  MET$sep_tercile[idx]<-ifelse(MET$min_sep[idx]<=bounds[1],"low",
                   ifelse(MET$min_sep[idx]<=bounds[2],"mid","high"))
  sep_thresholds[[k]]<-data.frame(study=group_keys$study[k],NOS=group_keys$NOS[k],
     sep_q33=bounds[1],sep_q67=bounds[2],
     n_group=length(idx),stringsAsFactors=FALSE)
}
sep_thresholds<-do.call(rbind,sep_thresholds)

takeq <- function(x,q) {
  if(!length(x) || all(!is.finite(x))) return(rep(NA_real_,length(q)))
  as.numeric(quantile(x[is.finite(x)],q,names=FALSE))
}
dist_rows<-vector("list",nrow(group_keys))
coverage_rows<-list();ii<-0L
levels_sep<-c("low","mid","high")
levels_pos<-c("center_0_1_3","middle_1_3_2_3","outer_2_3_1")
for(k in seq_len(nrow(group_keys))) {
  study<-group_keys$study[k];nos<-group_keys$NOS[k]
  v<-MET[MET$study==study & MET$NOS==nos,,drop=FALSE]
  sq<-takeq(v$min_sep,c(.05,.25,.5,.75,.95))
  rq<-takeq(v$Rpos_obs,c(.05,.25,.5,.75,.95))
  dist_rows[[k]]<-data.frame(study=study,NOS=nos,n=nrow(v),
    sep_q05=sq[1],sep_q25=sq[2],sep_q50=sq[3],sep_q75=sq[4],sep_q95=sq[5],
    Rpos_q05=rq[1],Rpos_q25=rq[2],Rpos_q50=rq[3],Rpos_q75=rq[4],Rpos_q95=rq[5],
    Fry50_count=sum(v$Fry50_center_proxy_obs),
    Fry50_fraction=mean(v$Fry50_center_proxy_obs),
    median_DFC_scaled=median(v$DFC_scaled),
    median_DFM_min_scaled=median(v$DFM_min_scaled),
    median_hull_area=median(v$hull_area_normalized),
    median_nonvertex_sources=median(v$n_nonvertex_sources),
    stringsAsFactors=FALSE)
  for(s in levels_sep) for(p in levels_pos) {
    ii<-ii+1L
    coverage_rows[[ii]]<-data.frame(study=study,NOS=nos,
       sep_tercile=s,pos_bin=p,
       n_available=sum(v$sep_tercile==s & v$pos_bin==p),
       group_n=nrow(v),stringsAsFactors=FALSE)
  }
}
DISTRIBUTION<-do.call(rbind,dist_rows)
COVERAGE<-do.call(rbind,coverage_rows)
COVERAGE$frac_of_group<-COVERAGE$n_available/COVERAGE$group_n

# Summarize overlap of CENTRAL 90% ranges, without forcing unsupported
# joint strata. Across NOS within each study and across studies within NOS.
common_rows<-list();rr<-0L
for(axis in c("across_NOS_within_study","across_studies_within_NOS")) {
  anchors<-if(axis=="across_NOS_within_study") unique(MET$study) else sort(unique(MET$NOS))
  for(anchor in anchors) {
    if(axis=="across_NOS_within_study") {
      sub<-MET[MET$study==anchor,,drop=FALSE]; groups<-split(sub,sub$NOS)
    } else {
      sub<-MET[MET$NOS==anchor,,drop=FALSE]; groups<-split(sub,sub$study)
    }
    sl<-lapply(groups,function(d) takeq(d$min_sep,COMMON_Q))
    pl<-lapply(groups,function(d) takeq(d$Rpos_obs,COMMON_Q))
    slo<-max(vapply(sl,`[`,numeric(1),1L));shi<-min(vapply(sl,`[`,numeric(1),2L))
    plo<-max(vapply(pl,`[`,numeric(1),1L));phi<-min(vapply(pl,`[`,numeric(1),2L))
    has_sep<-slo<=shi;has_pos<-plo<=phi
    for(g in names(groups)) {
      v<-groups[[g]]
      rr<-rr+1L
      common_rows[[rr]]<-data.frame(scope=axis,anchor=as.character(anchor),group=g,
        n_group=nrow(v),
        sep_overlap_low=if(has_sep) slo else NA_real_,
        sep_overlap_high=if(has_sep) shi else NA_real_,
        Rpos_overlap_low=if(has_pos) plo else NA_real_,
        Rpos_overlap_high=if(has_pos) phi else NA_real_,
        n_in_sep_overlap=if(has_sep) sum(v$min_sep>=slo & v$min_sep<=shi) else 0L,
        n_in_Rpos_overlap=if(has_pos) sum(v$Rpos_obs>=plo & v$Rpos_obs<=phi) else 0L,
        n_in_joint_overlap=if(has_sep && has_pos)
           sum(v$min_sep>=slo & v$min_sep<=shi &
               v$Rpos_obs>=plo & v$Rpos_obs<=phi) else 0L,
        stringsAsFactors=FALSE)
    }
  }
}
COMMON<-do.call(rbind,common_rows)

safe_cor <- function(a,b) {
  ok<-is.finite(a) & is.finite(b)
  if(sum(ok)<5L || length(unique(a[ok]))<3L || length(unique(b[ok]))<3L)
    return(NA_real_)
  suppressWarnings(cor(a[ok],b[ok],method="spearman"))
}
cor_rows<-list();kk<-0L
append_cors<-function(v,scope,key) {
  data.frame(scope=scope,group=key,n=nrow(v),
    rho_minsep_Rpos=safe_cor(v$min_sep,v$Rpos_obs),
    rho_minsep_DFC=safe_cor(v$min_sep,v$DFC_scaled),
    rho_Rpos_DFC=safe_cor(v$Rpos_obs,v$DFC_scaled),
    rho_NOS_minsep=safe_cor(v$NOS,v$min_sep),
    rho_NOS_Rpos=safe_cor(v$NOS,v$Rpos_obs),
    rho_minsep_DFM_min=safe_cor(v$min_sep,v$DFM_min_scaled),
    stringsAsFactors=FALSE)
}
cor_rows[[length(cor_rows)+1L]]<-append_cors(MET,"all","all")
for(s in unique(MET$study)) cor_rows[[length(cor_rows)+1L]]<-
  append_cors(MET[MET$study==s,,drop=FALSE],"study",s)
for(n in sort(unique(MET$NOS))) cor_rows[[length(cor_rows)+1L]]<-
  append_cors(MET[MET$NOS==n,,drop=FALSE],"NOS",as.character(n))
for(k in seq_len(nrow(group_keys))) {
  s<-group_keys$study[k];n<-group_keys$NOS[k]
  cor_rows[[length(cor_rows)+1L]]<-
    append_cors(MET[MET$study==s & MET$NOS==n,,drop=FALSE],
                "study_x_NOS",paste0(s," | NOS=",n))
}
CORRELATIONS<-do.call(rbind,cor_rows)

# A compact decision summary; no arbitrary claim that a given final n is enough.
SUMMARY<-data.frame(
  item=c("total_generated_candidates","valid_complete_candidates",
         "invalid_or_incomplete_candidates","n_studies","n_NOS_levels",
         "n_study_NOS_groups","n_3x3_cells_total","n_3x3_empty_cells",
         "n_3x3_cells_below_3","n_3x3_cells_below_5",
         "n_Fry50_center_proxy","fraction_Fry50_center_proxy",
         "note"),
  value=as.character(c(nrow(candidates),nrow(MET),nrow(flagged),
           length(unique(MET$study)),length(unique(MET$NOS)),nrow(group_keys),
           nrow(COVERAGE),sum(COVERAGE$n_available==0),
           sum(COVERAGE$n_available<3),sum(COVERAGE$n_available<5),
           sum(MET$Fry50_center_proxy_obs),mean(MET$Fry50_center_proxy_obs),
           "Coverage audit only; no final selection or MixSIAR runs")),
  stringsAsFactors=FALSE)
DEFINITIONS<-data.frame(
  field=c("min_sep","median_sep","Rpos_obs","Rpos_true",
          "DFC_original_proxy","DFC_scaled","DFM_raw_permil","DFM_scaled",
          "Fry50_center_proxy_obs","sep_tercile","pos_bin",
          "COMMON ranges","centroid definition","hull_vertex","source SD note"),
  definition=c(
    "Minimum pairwise sqrt((d18O_i-d18O_j)^2/(sd18_i^2+sd18_j^2)+(d2H_i-d2H_j)^2/(sd2H_i^2+sd2H_j^2)); same as V3.2",
    "Median pairwise standardized distance; supplementary",
    "Observed six-sample mixture MEAN: distance from hull vertex-mean center to mixture / center-to-boundary ray distance; 0 center, 1 boundary",
    "Same radial measure using true mixing isotope value; sensitivity only",
    "Observed mix-to-hull-vertex-mean Euclidean distance in raw delta units / raw bounding rectangle diagonal; legacy-compatible proxy, d2H dominates",
    "Observed mix-to-center Euclidean distance AFTER normalizing d18O and d2H by their respective source spans / sqrt(2)",
    "Raw Euclidean source-to-observed-mixture distance in the isotope delta coordinate system, per mille; dominated by d2H scale",
    "Source-to-observed-mixture Euclidean distance after each tracer is normalized by source-span, divided by sqrt(2)",
    "Rpos_obs <= 0.5; SIMPLIFIED center-50%-scaled hull proxy; NOT a direct implementation of Fry 2013 minmax danger-zone boundaries",
    "Within each study x NOS: min_sep split at the observed 1/3 and 2/3 quantiles; used only for coverage diagnostics, not modeling",
    "Absolute Rpos_obs: center 0..1/3, middle (1/3)..(2/3), outer (2/3)..1; NOT separately quantiled",
    "For each study: intersection of study-specific central-90%-ranges across all NOS; or for each NOS across studies; diagnostic only, not proof of joint support",
    "Arithmetic mean of source convex-hull vertices in normalized isotope coordinates; inside hull; not area-weighted polygon centroid",
    "TRUE if the source is a vertex of the convex hull; internal source may make NOS larger than the number of hull vertices",
    "Separation uses only per-tracer marginal source SD (no cross-tracer covariance); WATSON source SD calibration may rely on external study"),
  stringsAsFactors=FALSE)

OUT <- file.path(if(basename(INPUT)==INPUT_FOLDER) dirname(INPUT) else INPUT,
                 OUTPUT_FOLDER)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
filename<-file.path(OUT,"scenario_factor_coverage_V3_2.xlsx")
manifest<-data.frame(item=c("DIAGNOSTIC_SCRIPT","INPUT_TYPE","INPUT_PATH",
  "DATE_LOCAL","R_VERSION","RNG_USED","SOURCE_GENERATOR","RPOS_BASE",
  "Fry50_DEFINITION","COMMON_Q","INPUT_MD5_FILES","INPUT_MD5_VALUES"),
  value=c("V1.0",input_type,INPUT,as.character(Sys.Date()),R.version.string,
  "NONE: deterministic diagnostic, no set.seed() required","V3.2",
  "Observed six-sample mixture mean; true mixture is retained as sensitivity",
  "Rpos <= 0.50; geometric proxy, NOT Fry's exact source minmax zones",
  paste(COMMON_Q,collapse=","),paste(basename(used_paths),collapse="; "),
  paste(unname(tools::md5sum(used_paths)),collapse="; ")),
  stringsAsFactors=FALSE)
all_sheets <- list(
  Summary=SUMMARY,Manifest=manifest,Definitions=DEFINITIONS,
  Group_distribution=DISTRIBUTION,Separation_thresholds=sep_thresholds,
  Coverage_3x3=COVERAGE,Common_support=COMMON,Correlations=CORRELATIONS,
  Scenario_metrics=MET,Source_DFM=DFM,Excluded_cases=flagged)
writexl::write_xlsx(all_sheets,path=filename)
message("-------------------------------------------------------------")
message("Diagnostic Excel: ",normalizePath(filename,winslash="/"))
message("Valid scenarios: ",nrow(MET),"; flagged: ",nrow(flagged))
message("Empty 3x3 cells: ",sum(COVERAGE$n_available==0),
        " / ",nrow(COVERAGE),"; cells with fewer than 5: ",
        sum(COVERAGE$n_available<5))
message("Geometry-audit stage complete: no random draws, no selection, no MixSIAR.")


# ---- STEP 3 ---------------------------------------------------------------
# SCIENTIFIC OBJECTIVE
# Select a NOS-balanced subset WITHOUT looking at posterior accuracy or AE.
# Within each NOS (=3,...,7):
#   1) Select 140 candidates by simple random sampling (seed 20261024 + NOS).
#   2) Protect nonempty 3 x 3 geometry cells based on source separation
#      terciles and radial position thirds (up to 3 examples per cell).
#   3) Include separation/position extrema.
#   4) Fill remaining slots with greedy maximin (Kennard-Stone style) in
#      the (ranked separation, radial position) plane, producing a fixed
#      priority list. Retain 200 per NOS = 1,000 for the primary sample.
# The 850 and 1,150 alternative subsets share the same priority lists and
# are diagnostics, NOT additional MixSIAR model runs.
# OUTPUT: selection_final_V1/ (selected IDs, sources and mixture replicates).
# --------------------------------------------------------------------------
stage_message(3, "Select 1,000 NOS-balanced MixSIAR input scenarios")
OUTPUT_FOLDER <- SELECTION_OUTPUT_FOLDER
if (!requireNamespace('readxl',quietly=TRUE))
  stop("Install missing package: install.packages('readxl')")
if (!requireNamespace('writexl',quietly=TRUE))
  stop("Install missing package: install.packages('writexl')")

ROOT <- get_script_dir()
# The consolidated script lives in the project root (03_Simulation/).
ROOT_CANDIDATES <- unique(c(ROOT,dirname(ROOT),
                            file.path(ROOT,'isotope_final_selection_R')))
root_ok <- function(p) file.exists(file.path(p,COVERAGE_FILE)) &&
  (all(file.exists(file.path(p,RAW_FOLDER,c('01_candidate_scenarios.csv',
                                           '02_candidate_sources.csv',
                                           '03_candidate_mixtures.csv')))) ||
     file.exists(file.path(p,RAW_FOLDER,RAW_EXCEL)) ||
     file.exists(file.path(p,RAW_EXCEL)))
ok_roots <- ROOT_CANDIDATES[vapply(ROOT_CANDIDATES,root_ok,logical(1))]
if (!length(ok_roots)) stop(
 "Cannot find coverage_diagnostics/scenario_factor_coverage_V3_2.xlsx and V3.2 raw files together.\n",
 "Run from the parent folder containing coverage_diagnostics/ and outputs_v3_2/.")
BASE <- ok_roots[1L]
coverage_path <- file.path(BASE,COVERAGE_FILE)
read_s <- function(file,sheet) as.data.frame(readxl::read_excel(file,sheet=sheet),
                                               stringsAsFactors=FALSE)
M <- read_s(coverage_path,'Scenario_metrics')
needed_metrics <- c('scenario_id','study','NOS','min_sep','Rpos_obs',
                    'DFC_scaled','DFM_min_scaled',
                    'sample_mean_shift_d18O','sample_mean_shift_d2H')
check_cols <- function(d,cols,who) {
  missing <- setdiff(cols,names(d))
  if(length(missing)) stop(who,' missing: ',paste(missing,collapse=', '))
}
check_cols(M,needed_metrics,'Scenario_metrics')
M$scenario_id <- as.integer(M$scenario_id)
M$NOS <- as.integer(M$NOS)
for (nm in c('min_sep','Rpos_obs','DFC_scaled','DFM_min_scaled'))
  M[[nm]] <- as.numeric(M[[nm]])
if(anyNA(M$scenario_id) || anyDuplicated(M$scenario_id))
  stop('Duplicate/invalid scenario_id in Scenario_metrics')
if(any(!is.finite(M$min_sep)) || any(M$min_sep<0) ||
   any(!is.finite(M$Rpos_obs)) || any(M$Rpos_obs < -1e-8 | M$Rpos_obs>1+1e-8))
  stop('Invalid min_sep / Rpos_obs values')
if(anyNA(M$study) || anyNA(M$NOS)) stop('Missing study or NOS')

csv_names <- c('01_candidate_scenarios.csv','02_candidate_sources.csv',
               '03_candidate_mixtures.csv')
raw_csv_paths <- file.path(BASE,RAW_FOLDER,csv_names)
if(all(file.exists(raw_csv_paths))) {
  raw_paths <- raw_csv_paths
  raw_type <- 'three_original_V3.2_CSVs'
  C <- read.csv(raw_csv_paths[1],stringsAsFactors=FALSE,check.names=FALSE)
  S <- read.csv(raw_csv_paths[2],stringsAsFactors=FALSE,check.names=FALSE)
  X <- read.csv(raw_csv_paths[3],stringsAsFactors=FALSE,check.names=FALSE)
} else {
  raw_p <- file.path(BASE,RAW_FOLDER,RAW_EXCEL)
  if(!file.exists(raw_p)) raw_p <- file.path(BASE,RAW_EXCEL)
  if(!file.exists(raw_p)) stop('V3.2 raw Candidates/Sources/Mixtures not found')
  raw_paths <- raw_p
  raw_type <- 'original_V3.2_Excel'
  C <- read_s(raw_p,'Candidates')
  S <- read_s(raw_p,'Sources')
  X <- read_s(raw_p,'Mixtures')
}
check_cols(C,c('scenario_id','NOS','study','generation_status',
               'calibration_status','n_obs_generated',
               'min_separation_scaled','mix_obs_mean_d18O','mix_obs_mean_d2H',
               'mix_true_d18O','mix_true_d2H'),
               'Candidates')
check_cols(S,c('scenario_id','source_index','p_true'),'Sources')
check_cols(X,c('scenario_id','replicate_id','d18O','d2H'),'Mixtures')
C$scenario_id<-as.integer(C$scenario_id); C$NOS<-as.integer(C$NOS)
C$n_obs_generated<-as.integer(C$n_obs_generated)
S$scenario_id<-as.integer(S$scenario_id)
X$scenario_id<-as.integer(X$scenario_id)
if(anyDuplicated(C$scenario_id)) stop('Repeated IDs in raw Candidates')
valid <- C$generation_status=='complete' & C$calibration_status=='OK' &
         C$n_obs_generated == EXPECTED_MIX_REPS
valid[is.na(valid)] <- FALSE
valid_C <- C[valid,,drop=FALSE]
if(!setequal(M$scenario_id,valid_C$scenario_id))
  stop('CRITICAL: Coverage metrics IDs do not match complete scenarios in raw V3.2.\n',
       'Do not mix different V3.2 runs or different Excel versions.')
mc <- match(M$scenario_id,valid_C$scenario_id)
if(any(M$NOS != valid_C$NOS[mc]) ||
   any(as.character(M$study)!=as.character(valid_C$study[mc])))
  stop('CRITICAL: scenario_id joins to different NOS/study across files')
# Guard against different V3.2 reruns that coincidentally share all scenario IDs.
if(any(abs(M$min_sep-as.numeric(valid_C$min_separation_scaled[mc]))>1e-8,na.rm=TRUE) ||
   any(abs(M$sample_mean_shift_d18O-
           (as.numeric(valid_C$mix_obs_mean_d18O[mc])-
            as.numeric(valid_C$mix_true_d18O[mc])))>1e-8,na.rm=TRUE) ||
   any(abs(M$sample_mean_shift_d2H-
           (as.numeric(valid_C$mix_obs_mean_d2H[mc])-
            as.numeric(valid_C$mix_true_d2H[mc])))>1e-8,na.rm=TRUE))
  stop('CRITICAL: coverage metrics and raw candidates are from different V3.2 runs')
ns <- sort(unique(M$NOS))
if(!identical(ns,3:7)) stop('Expected five NOS levels 3:7, got ',paste(ns,collapse=','))
if(anyDuplicated(S[,c('scenario_id','source_index')])) stop('Duplicate source rows')
if(anyDuplicated(X[,c('scenario_id','replicate_id')])) stop('Duplicate replicate rows')

# Source separation terciles are computed WITHIN NOS across all 3 studies.
# Position thirds remain absolute [0,1/3],(1/3,2/3],(2/3,1].
# The grid protects joint combinations; no study, DFC, DFM or true p used.
M$sep_bin <- NA_integer_
M$pos_bin <- ifelse(M$Rpos_obs<=1/3,1L,
                    ifelse(M$Rpos_obs<=2/3,2L,3L))
for(nos in ns) {
  idx<-which(M$NOS==nos)
  q<-as.numeric(quantile(M$min_sep[idx],c(1/3,2/3),names=FALSE))
  if(q[1]>=q[2]) stop('Source separation quantile ties in NOS=',nos)
  M$sep_bin[idx]<-ifelse(M$min_sep[idx]<=q[1],1L,
                         ifelse(M$min_sep[idx]<=q[2],2L,3L))
}
M$joint_cell <- paste0('S',M$sep_bin,'_P',M$pos_bin)

# Within-NOS space filling uses separation PERCENTILE rank and Rpos directly.
# This prevents unusually large source distances from dominating the metric.
# The rank is only used to SELECT samples, never used as the scientific predictor.
priority_parts <- list()
for(nos in ns) {
  G <- M[M$NOS==nos,,drop=FALSE]
  G <- G[order(G$scenario_id),,drop=FALSE]
  n <- nrow(G)
  if(n <= RANDOM_BACKBONE_PER_NOS) stop('Too few candidates for NOS=',nos)
  set.seed(SELECTION_SEED+nos)  # reproducible independent NOS lists
  Z <- cbind(separation_rank=(rank(G$min_sep,ties.method='average')-1)/(n-1),
             mixture_position=G$Rpos_obs)
  D2 <- as.matrix(dist(Z,method='euclidean'))^2
  diag(D2) <- 0
  sel <- integer(0); reasons<-character(0)
  add_one <- function(i, reason) {
    if(!i %in% sel) {sel <<- c(sel,i); reasons <<- c(reasons,reason)}
  }
  choose_cover <- function(candidates) {
    candidates<-setdiff(as.integer(candidates),sel)
    if(!length(candidates)) stop('No unselected scenarios left')
    if(!length(sel)) return(candidates[which.min(G$scenario_id[candidates])])
    # maximize distance from the nearest selected scenario
    ds <- D2[candidates,sel,drop=FALSE]
    nearest <- apply(ds,1L,min)
    # deterministic ID tie-break
    candidates[order(-nearest,G$scenario_id[candidates])[1L]]
  }
  # Stage 1: genuine NOS-stratified simple random baseline.
  for(i in sample.int(n,RANDOM_BACKBONE_PER_NOS,replace=FALSE))
    add_one(i,'random_baseline')
  # Stage 2: preserve at least 3 cases per *existing* 3x3 joint cell.
  # Process the scarcest cells first; never synthesize unsupported cases.
  cell_counts<-sort(table(G$joint_cell),decreasing=FALSE)
  for(cell in names(cell_counts)) {
    cell_members<-which(G$joint_cell==cell)
    goal<-min(MIN_PER_NONEMPTY_GRID_CELL,length(cell_members))
    while(sum(sel %in% cell_members)<goal)
      add_one(choose_cover(cell_members),'joint_cell_protection')
  }
  # Stage 3: retain within-NOS extrema of both main continuous variables.
  for(i in unique(c(which.min(G$min_sep),which.max(G$min_sep),
                    which.min(G$Rpos_obs),which.max(G$Rpos_obs))))
    add_one(i,'range_endpoint')
  # Stage 4: sequential greedy maximin / Kennard-Stone-style augmentation.
  while(length(sel)<n) add_one(choose_cover(seq_len(n)),'joint_space_fill')
  priority_parts[[as.character(nos)]] <- data.frame(
    scenario_id=G$scenario_id[sel],NOS=nos,
    within_NOS_priority=seq_along(sel),selection_reason=reasons,
    stringsAsFactors=FALSE)
}
PRIORITY <- do.call(rbind,priority_parts)
row.names(PRIORITY)<-NULL
if(anyDuplicated(PRIORITY$scenario_id) || !setequal(PRIORITY$scenario_id,M$scenario_id))
  stop('Invalid priority index: duplicates or missing IDs')

nsizes <- sort(unique(c(as.integer(COMPARE_N),as.integer(PRIMARY_N))))
quota_for_n <- function(k) {
  q <- rep(k %/% length(ns),length(ns))
  if(k %% length(ns)) q[seq_len(k%%length(ns))] <- q[seq_len(k%%length(ns))]+1L
  q
}
if(any(nsizes<length(ns)*RANDOM_BACKBONE_PER_NOS) ||
   any(nsizes>nrow(M))) stop('Selection N outside feasible range')
for(k in nsizes) {
  q<-quota_for_n(k)
  if(any(q>as.numeric(table(factor(M$NOS,levels=ns)))))
    stop('N=',k,' exceeds available candidates in at least one NOS')
}

# Diagnostics on discrete joint factor coverage and 2D nearest-case coverage.
make_grid <- function(selected_ids,N) {
  d<-merge(M[,c('scenario_id','NOS','joint_cell')],
           data.frame(scenario_id=selected_ids,selected=TRUE),
           by='scenario_id',all.x=TRUE,sort=FALSE)
  d$selected[is.na(d$selected)]<-FALSE
  tab<-aggregate(cbind(available=rep(1L,nrow(d)),
                       selected=as.integer(d$selected)),
                 by=d[c('NOS','joint_cell')],FUN=sum)
  tab$N<-N
  tab$missed_cell <- tab$selected==0L
  tab$below_protection <- tab$selected<pmin(MIN_PER_NONEMPTY_GRID_CELL,tab$available)
  tab[order(tab$NOS,tab$joint_cell),,drop=FALSE]
}

nn_quantiles <- function(selected_ids) {
  ans<-numeric(0)
  for(nos in ns) {
    G<-M[M$NOS==nos,,drop=FALSE]
    G<-G[order(G$scenario_id),,drop=FALSE]
    ix<-which(G$scenario_id %in% selected_ids)
    zz<-cbind((rank(G$min_sep,ties.method='average')-1)/(nrow(G)-1),
              G$Rpos_obs)
    dd<-as.matrix(dist(zz))
    nn<-apply(dd[,ix,drop=FALSE],1,min)
    ans<-c(ans,nn)
  }
  c(q95=as.numeric(quantile(ans,.95,names=FALSE)),
    max=max(ans),median=median(ans))
}
ks_distance <- function(a,b) {
  v<-sort(unique(c(a,b)))
  max(abs(ecdf(a)(v)-ecdf(b)(v)))
}
summary_list<-list();grid_list<-list();shift_list<-list();study_list<-list()
all_ids<-list()
for(k in nsizes) {
  quota<-quota_for_n(k)
  selected_ids<-unlist(lapply(seq_along(ns),function(j) {
    p<-PRIORITY[PRIORITY$NOS==ns[j],,drop=FALSE]
    p$scenario_id[seq_len(quota[j])]
  }),use.names=FALSE)
  if(length(selected_ids)!=k || anyDuplicated(selected_ids)) stop('N quota mismatch')
  all_ids[[as.character(k)]]<-selected_ids
  grid<-make_grid(selected_ids,k)
  grid_list[[as.character(k)]]<-grid
  nn<-nn_quantiles(selected_ids)
  summary_list[[as.character(k)]]<-data.frame(
    N=k,fraction_retained=k/nrow(M),
    groups_NOS=length(ns),joint_cells_existing=nrow(grid),
    joint_cells_selected=sum(grid$selected>0),
    cells_below_min3=sum(grid$below_protection),
    joint_coverage_pct=100*mean(grid$selected>0),
    q95_nearest_distance=unname(nn['q95']),
    max_nearest_distance=unname(nn['max']),
    median_nearest_distance=unname(nn['median']),
    random_baseline_n=sum(PRIORITY$scenario_id %in% selected_ids &
                        PRIORITY$selection_reason=='random_baseline'),
    mean_source_sep_selected=mean(M$min_sep[M$scenario_id %in% selected_ids]),
    mean_Rpos_selected=mean(M$Rpos_obs[M$scenario_id %in% selected_ids]),
    stringsAsFactors=FALSE)
  for(nos in c('ALL',as.character(ns))) {
    pool<-if(nos=='ALL') M else M[M$NOS==as.integer(nos),,drop=FALSE]
    samp<-pool[pool$scenario_id %in% selected_ids,,drop=FALSE]
    for(v in c('min_sep','Rpos_obs')) {
      qp<-as.numeric(quantile(pool[[v]],c(.05,.5,.95),names=FALSE))
      qs<-as.numeric(quantile(samp[[v]],c(.05,.5,.95),names=FALSE))
      shift_list[[length(shift_list)+1L]]<-data.frame(
        N=k,NOS=nos,predictor=v,n_available=nrow(pool),n_selected=nrow(samp),
        q05_available=qp[1],q05_selected=qs[1],
        median_available=qp[2],median_selected=qs[2],
        q95_available=qp[3],q95_selected=qs[3],
        KS_distance=ks_distance(pool[[v]],samp[[v]]),
        min_retained=min(samp[[v]])<=min(pool[[v]])+1e-12,
        max_retained=max(samp[[v]])>=max(pool[[v]])-1e-12,
        stringsAsFactors=FALSE)
    }
  }
  for(nos in ns) for(study in sort(unique(M$study))) {
    d<-M[M$NOS==nos & M$study==study,,drop=FALSE]
    study_list[[length(study_list)+1L]]<-data.frame(
      N=k,NOS=nos,study=study,n_available=nrow(d),
      n_selected=sum(d$scenario_id %in% selected_ids),stringsAsFactors=FALSE)
  }
}
SUMMARY<-do.call(rbind,summary_list)
GRID<-do.call(rbind,grid_list)
SHIFT<-do.call(rbind,shift_list)
STUDY_AUDIT<-do.call(rbind,study_list)
row.names(SUMMARY)<-NULL;row.names(GRID)<-NULL
row.names(SHIFT)<-NULL;row.names(STUDY_AUDIT)<-NULL

sel_ids<-all_ids[[as.character(PRIMARY_N)]]
if(is.null(sel_ids)) stop('PRIMARY_N not found')
FINAL<-merge(PRIORITY[PRIORITY$scenario_id %in% sel_ids,,drop=FALSE],
             M,by=c('scenario_id','NOS'),sort=FALSE)
FINAL<-FINAL[order(FINAL$NOS,FINAL$within_NOS_priority),,drop=FALSE]
if(nrow(FINAL)!=PRIMARY_N) stop('Final count mismatch')
if(any(FINAL$selection_reason=='random_baseline') &&
   any(table(FINAL$NOS[FINAL$selection_reason=='random_baseline']) !=
       RANDOM_BACKBONE_PER_NOS)) stop('Random-baseline counts inconsistent')

selected_candidates<-C[match(sel_ids,C$scenario_id),,drop=FALSE]
selected_sources<-S[S$scenario_id %in% sel_ids,,drop=FALSE]
selected_mixtures<-X[X$scenario_id %in% sel_ids,,drop=FALSE]
if(anyNA(selected_candidates$scenario_id) ||
   nrow(selected_sources)!=sum(selected_candidates$NOS) ||
   nrow(selected_mixtures)!=PRIMARY_N*EXPECTED_MIX_REPS ||
   any(table(selected_sources$scenario_id)!=selected_candidates$NOS[match(
        as.integer(names(table(selected_sources$scenario_id))),
        selected_candidates$scenario_id)]) ||
   any(table(selected_mixtures$scenario_id)!=EXPECTED_MIX_REPS))
  stop('CRITICAL: selected raw data are incomplete or inconsistent')

OUT<-file.path(BASE,OUTPUT_FOLDER)
dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
write.csv(FINAL,file.path(OUT,'01_selected_IDs_and_factors.csv'),row.names=FALSE,
          fileEncoding='UTF-8')
write.csv(selected_candidates,file.path(OUT,'02_selected_Candidates.csv'),
          row.names=FALSE,fileEncoding='UTF-8')
write.csv(selected_sources,file.path(OUT,'03_selected_Sources.csv'),
          row.names=FALSE,fileEncoding='UTF-8')
write.csv(selected_mixtures,file.path(OUT,'04_selected_Mixtures.csv'),
          row.names=FALSE,fileEncoding='UTF-8')
write.csv(PRIORITY,file.path(OUT,'05_all_candidate_priority.csv'),
          row.names=FALSE,fileEncoding='UTF-8')

# Provenance: selection is now embedded in this single consolidated script.
script_candidates<-c(file.path(BASE,'SCED_generate_and_select_scenarios.R'))
script_file<-script_candidates[file.exists(script_candidates)]
manifest <- data.frame(
  property=c('VERSION','SELECTION_SEED','GENERATION_SEED','PRIMARY_N',
             'COMPARE_N','RANDOM_BACKBONE_PER_NOS','MIN_PER_NONEMPTY_GRID_CELL',
             'SELECTION_VARS','GRID_DEFINITION','STUDY_IN_SELECTION',
             'OUTCOME_VARIABLES_IN_SELECTION','SPATIAL_SCALE','PRIORITY_NESTED',
             'RAW_TYPE','METRICS_SOURCE_MD5','RAW_SOURCE_MD5',
             'SELECT_SCRIPT_MD5','R_VERSION','OUTPUT_PATH'),
  value=c('V1.0',as.character(SELECTION_SEED),as.character(GLOBAL_SEED),as.character(PRIMARY_N),
          paste(nsizes,collapse=','),as.character(RANDOM_BACKBONE_PER_NOS),
          as.character(MIN_PER_NONEMPTY_GRID_CELL),
          'NOS; min_sep; Rpos_obs',
          'NOS-wise min_sep tercile x absolute Rpos thirds; minimum 3 if available',
          'NO (audit only)','NO (MixSIAR has NOT been run)',
          'Within-NOS min_sep empirical rank x Rpos absolute',
          'YES: 850 subset of 1000 subset of 1150',raw_type,
          unname(tools::md5sum(coverage_path)),
          paste(unname(tools::md5sum(raw_paths)),collapse=';'),
          if(length(script_file)) unname(tools::md5sum(script_file[1])) else 'not_found',
          R.version.string,OUT),
  stringsAsFactors=FALSE)
INFO <- data.frame(
  item=c('Interpretation','Unit_of_analysis','Scientific_variables',
         'Random_baseline','Coverage_augmentation','Study_handling',
         'Outcome_blinding','Meaning_of_DFC_DFM_Fry',
         'What_is_not_claimed','Nested_alternatives'),
  definition=c(
    'Candidate subset supports joint-factor analysis; not a random population-prevalence sample',
    'Scenario; one MixSIAR fit per scenario with six mixture observations',
    'NOS, original continuous min_sep, original continuous Rpos_obs',
    'First 140 candidates within each NOS selected by simple random sampling',
    'Protect joint 3x3 cells with at least min(3, available); save four factor extrema; fill by greedy maximin',
    'Three empirical studies only generate parameter diversity; not stratification, quota or objective',
    'No MixSIAR accuracy/uncertainty results, true p, DFC, DFM, or danger zone used in selection',
    'Retained for later source-level or diagnostic analyses, not independent strata',
    'No causal effect inferred from selection; no claim of population-averaged MixSIAR performance',
    'Within each NOS, priority list is fixed; 850 is a subset of 1000 which is a subset of 1150'),
  stringsAsFactors=FALSE)
summary_xlsx <- file.path(OUT,'scenario_final_selection_audit_V1.xlsx')
workbook <- list(
  Summary=SUMMARY,Manifest=manifest,Definitions=INFO,
  Selected_primary=FINAL,Grid_coverage=GRID,
  Quantile_and_KS=SHIFT,Study_audit_only=STUDY_AUDIT,
  Candidate_priority=PRIORITY,
  Selected_850=data.frame(scenario_id=all_ids[['850']]),
  Selected_1000=data.frame(scenario_id=all_ids[['1000']]),
  Selected_1150=data.frame(scenario_id=all_ids[['1150']]))
# Alternatives list names are configurable; omit the standard sheets when absent.
for(n in c(850L,1000L,1150L)) if(!as.character(n) %in% names(all_ids))
  workbook[[paste0('Selected_',n)]] <- NULL
writexl::write_xlsx(workbook,path=summary_xlsx)
cat('\nScenario selection V1 complete.\n')
cat('Valid candidates:',nrow(M),' | PRIMARY_N:',PRIMARY_N,
    ' | Full raw mixture rows:',nrow(selected_mixtures),'\n')
print(SUMMARY[,c('N','joint_cells_existing','joint_cells_selected',
                 'cells_below_min3','q95_nearest_distance','max_nearest_distance')],
      row.names=FALSE)
cat('Saved:',summary_xlsx,'\n')
cat('Study used ONLY as an audit label. No MixSIAR has been run.\n')


# ---- FINAL CHECK ----------------------------------------------------------
# The selection script above already checks all source/mixture record counts,
# scenario identities, and raw-data/coverage consistency. The selected CSVs
# are direct inputs for the SUBSEQUENT, independent MixSIAR fitting workflow.
# This combined script intentionally does not run the posterior inference.
message("\nSCED scenario-generation workflow finished successfully.")
message("Candidate data: ", file.path(get_script_dir(), GEN_OUTPUT_FOLDER))
message("Geometry diagnostics: ", file.path(get_script_dir(), COVERAGE_OUTPUT_FOLDER))
message("Selected MixSIAR scenarios: ", file.path(get_script_dir(), SELECTION_OUTPUT_FOLDER))
message("Next external step: run MixSIAR on the 1,000 selected cases.")
