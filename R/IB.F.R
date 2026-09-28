#' @title Bootstrap confidence bands for fatigue lifetime
#'
#' @description It performs bootstrap confidence bands for fatigue lifetime. The lifetime
#' matrix is calculated by bootstrap resampling by means the above mentioned methodologies
#' (see craks.growth). The confidence bands are estimated by the quantile based method.
#'
#' @param z cracks.growth object.
#' @param nB Number of bootstrap resampling.
#' @param alpha	Confidence level.
#' @param method	Character string showing the distribution estimates method: "SEP-lme_bkde",
#' "SEP-lme_kde" or "PB-nlme.
#' @param seed	Seed for reproducibility.
#' @param ncores Controls how many R processes are launched in parallel.
#'
#' @details IB.F is performed from the output of cracks.growth function.
#'
#' @return
#' \describe{ Return the following values:
#'   \item{\code{Mat.F.B}}{Matrix that contents the fatigue lifetimes corresponding to
#'   each bootstrap resampling.}
#'   \item{\code{I.Bootstrap}}{Data frame that contents the bootstrap confidence bands
#'   for lifetime distribution, at a confidence level of 95 percent (by default). It is
#'   composed by two columns corresponding to the bands limits: low, up.}
#' }
#'
#' @author
#' Antonio Meneses antoniomenesesfreire@hotmail.com, Salvador Naya salva@udc.es,
#' Javier Tarrio-Saavedra jtarrio@udc.es, Ignacio Lopez-Ullibarri ilu@udc.es
#'
#' @references
#' Meeker, W., Escobar, L. (1998) Statistical Methods for Reliability Data. John Wiley & Sons,
#' Inc. New York.
#'
#' Pinheiro JC., Bates DM. (2000) Mixed-effects models in S ans S-plus. Statistics and Computing.
#' Springer-Verlang. New York.
#'
#' Paris, P.C. and Erdogan, F. (1963) A critical analysis of crack propagation laws. J. Basic Eng.,
#' 85, 528.
#'
#' @import KernSmooth
#' @import kerdiest
#' @import nlme
#' @import parallel
#' @importFrom grDevices boxplot.stats colors
#' @importFrom graphics abline lines points
#' @importFrom stats coef coefficients lm na.omit predict quantile resid rnorm selfStart spline var ave fitted
#'
#' @export
#' @examples
#' 
#' ## Using the Alea.A dataset
#' data(Alea.A)
#' x <- Alea.A
#' ## Critical crack length
#' aF <- 1.6
#' ## Censoring time
#' T_c <- 0.12
#' ## cracks.growth function applied to Alea.A data
#' cg <- cracks.growth (x, aF, T_c, method = c("SEP-lme_bkde", "SEP-lme_kde",
#'                     "PB-nlme"), nBKDE = 1000, nKDE = 1000, nMC = 5000)
#' ## z is a cracks.growth object
#' z <- cg
#' ## Number of bootstrap resamplings
#' nB <- 100
#' ## Application of IB.F function to cg object
#' ic.b <- IB.F(z, nB, alpha = 0.05, method = c("SEP-lme_bkde", "SEP-lme_kde",
#'                                              "PB-nlme"))
#' ## ic.b values obtainde by the "SEP-lme_bkde" model
#' names(ic.b)
#' # [1] "Mat.F.B"     "I.Bootstrap"
#' ## Chart with the empirical and estimated distribution functions,
#' ## with bootstrap confidence bands at 95
#' # Observations from which the distribution function is estimated
#' F1.F <- z$F.est[,2]
#' plot( ic.b$I.Bootstrap$low,F1.F, col=2, type="l", lty=2, lwd=2,
#'       xlim=c(0.05,0.18),
#'       main="Plot: distributions of failure times\n  confidence intervals",
#'       xlab="million cycles",  ylab="probability",  cex.lab=1.7,
#'       cex.main=2, las=1)
#' lines(ic.b$I.Bootstrap$up, F1.F, col=2, lty=2, lwd=2)
#' points(z$F.est, pch=20)
#' points(z$F.emp, col=4, pch=20, cex=1.5)
#' legend("topleft", c("Empirical", "Estimated","Bootstrap (95 percent)"),
#'        col=c("blue","black","red"),  lty=c(1,1,1), pch=c(20,20,20),
#'        cex=1.5, bty="n")
#' ## Graph with confidence bands
#' matplot(ic.b$Mat.F.B, F1.F,  main="Bootstrap resampling lines",
#'         type="l", lwd=2, xlim=c(0.05,0.18), xlab="million cycles",
#'         ylab="probability", cex.lab=1.7,  cex.main=2, las=1)
IB.F <- function(z, nB, alpha = 0.05,
                 method = c("SEP-lme_bkde", "SEP-lme_kde", "PB-nlme"),
                 seed = NULL,
                 ncores = min(2L, max(1L, parallel::detectCores() - 1L))) {
  
  method <- match.arg(method)
  
  
  if (!is.null(seed)) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }
  
  
  x  <- z$data
  a0 <- stats::ave(x[[2]], x[[3]], FUN = function(v) v[1])
  x[[4]] <- a0
  colnames(x) <- c("cycles", "cracks", "sample", "a0")
  samples   <- unique(x$sample)
  nS        <- length(samples)
  split.idx <- split(seq_len(nrow(x)), x$sample)
  aF        <- z$a.F
  param.gc  <- as.matrix(z$param)
  
  T_func <- function(a_0, aF, C, m)
    (1 - (aF / a_0)^(-(m - 2) / 2)) / ((a_0^(m/2 - 1)) * (m/2 - 1) * C * pi^(m/2))
  
  run_parallel <- function(FUN, n) {
    if (.Platform$OS.type == "unix") {
      parallel::mclapply(seq_len(n), FUN, mc.cores = ncores, mc.set.seed = TRUE)
    } else {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl))
      parallel::clusterEvalQ(cl, {
        library(mgcv); library(pspline); library(KernSmooth)
        library(kerdiest); library(MASS); library(nlme)
      })
      parallel::clusterExport(cl, varlist = ls(environment()), envir = environment())
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
      parallel::parLapply(cl, seq_len(n), FUN)
    }
  }
  
  
  if (method %in% c("SEP-lme_bkde", "SEP-lme_kde")) {
    
    sigma.rr <- z$sigma
    grid.n   <- if (method == "SEP-lme_bkde") z$nBKDE else z$nKDE
    bw       <- z$bw
    
    one.rep <- function(b) {
      bb      <- sample.int(nS, size = nS, replace = TRUE)
      param.b <- param.gc[bb, , drop = FALSE]
      
      grieta.BB <- unlist(Map(function(idx, i) {
        C <- param.b[i, 1]; m <- param.b[i, 2]
        cyc <- x$cycles[idx]; a_0 <- x$a0[idx]
        (a_0^(1 - m/2) + (1 - m/2) * C * pi^(m/2) * cyc)^(2/(2 - m)) +
          rnorm(length(cyc), 0, sigma.rr)
      }, split.idx, seq_len(nS)), use.names = FALSE)
      
      
      grieta.BB <- unlist(lapply(split.idx, function(idx) {
        v <- grieta.BB[idx]
        v + (x$cracks[idx][1] - v[1])
      }), use.names = FALSE)
      
      
      coef.B <- lapply(split.idx, function(idx) {
        y_i <- grieta.BB[idx] - x$a0[idx]
        sp  <- pspline::smooth.Pspline(x$cycles[idx], y_i, norder = 2, method = 3)
        dpr <- as.vector(predict(sp, x$cycles[idx], nderiv = 1))
        try(coefficients(lm(log(dpr[-1]) ~ log(sqrt(pi * grieta.BB[idx][-1])))),
            silent = TRUE)
      })
      if (any(vapply(coef.B, function(cc) inherits(cc, "try-error"), logical(1))))
        return(rep(NA_real_, 101))
      
      C.B <- vapply(coef.B, function(cc) exp(cc[1]), numeric(1))
      m.B <- vapply(coef.B, function(cc) cc[2], numeric(1))
      
      TT.M <- vapply(seq_len(nS), function(i) {
        a000 <- x$a0[split.idx[[i]]][1]
        T_func(a000, aF, C.B[i], m.B[i])
      }, numeric(1))
      if (any(!is.finite(TT.M))) return(rep(NA_real_, 101))
      
      LL2 <- seq(min(x$cycles), max(TT.M), length.out = 1001)
      
      TT.FALLO.B <- vapply(seq_len(nS), function(i) {
        a0v <- x$a0[split.idx[[i]]][1]
        g <- (a0v^(1 - m.B[i]/2) + (1 - m.B[i]/2) * C.B[i] * pi^(m.B[i]/2) * LL2)^(2/(2 - m.B[i]))
        LL2[which.min(abs(g - aF))]
      }, numeric(1))
      TT.FALLO.BB <- TT.FALLO.B[!is.na(TT.FALLO.B)]
      if (length(TT.FALLO.BB) < 3) return(rep(NA_real_, 101))
      
      if (method == "SEP-lme_bkde") {
        dens <- KernSmooth::bkde(TT.FALLO.BB, gridsize = grid.n,
                                 range.x = range(TT.FALLO.BB), bandwidth = bw)
        bin  <- diff(range(TT.FALLO.BB)) / grid.n
        F.x  <- cumsum(abs(dens$y)) * bin
        keep <- !duplicated(F.x)
        xg   <- dens$x[keep]
        F.x  <- unique(F.x) / max(F.x)
        sp <- spline(F.x, xg, n = 101, method = "fmm", xmin = 0, xmax = 100,
                     xout = seq(0, 1, length.out = 101), ties = mean)
        round(sp$y, 6)
      } else { # SEP-lme_kde
        range.F <- range(TT.FALLO.BB)
        FFF <- kerdiest::kde(type_kernel = "e", vec_data = TT.FALLO.BB,
                             y = seq(range.F[1], range.F[2], length.out = grid.n), bw = bw)
        ii <- order(FFF$grid)
        sp <- spline(FFF$Estimated_values[ii], FFF$grid[ii], n = 101, method = "fmm",
                     xmin = 0, xmax = 100, xout = seq(0, 1, length.out = 101), ties = mean)
        sort(round(sp$y, 6))
      }
    }
    
    reps  <- run_parallel(one.rep, nB)
    FF.BB <- do.call(cbind, reps)
    
    ICB <- t(FF.BB)
    ic  <- apply(ICB, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    I.Bootstrap <- data.frame(low = ic[1, ], up = ic[2, ])
    return(list(Mat.F.B = FF.BB, I.Bootstrap = I.Bootstrap))
  }
  
  
  sigma.r <- z$sigma
  nMC     <- z$nMC
  
  Paris.F <- function(cycles, C, m, a0, F.cte = 1, S = 1)
    (a0^(1 - m/2) + (1 - m/2) * C * F.cte^m * S^m * pi^(m/2) * cycles)^(2/(2 - m))
  
  DerivInit <- function(data, cycles, y) {
    xx <- data$cycles; y <- data$cracks
    derivada <- diff(y) / diff(xx)
    grieta.media <- (y[-1] + y[-length(y)]) / 2
    data.frame(cycles = xx[-1], derivada, grieta.media)
  }
  
  ParisInit <- function(mCall, data, LHS) {
    dd  <- DerivInit(data, mCall[["cycles"]], LHS)
    fit <- lm(log(derivada) ~ log(sqrt(pi * grieta.media)), data = dd)
    res <- c(exp(coefficients(fit)[1]), coefficients(fit)[2])
    names(res) <- mCall[c("C", "m")]
    res
  }
  
  SSParis <- selfStart(Paris.F, initial = ParisInit, parameters = c("C", "m"))
  
  one.rep.PB <- function(b, max_try = 200) {
    for (attempt in seq_len(max_try)) {
      bb <- sample.int(nS, size = nS, replace = TRUE)
      param.b <- param.gc[bb, , drop = FALSE]
      
      grieta.B <- unlist(Map(function(idx, i) {
        Cc <- param.b[i, 1]; mM <- param.b[i, 2]
        cyc <- x$cycles[idx]; a_0 <- x$a0[idx]
        (a_0^(1 - mM/2) + (1 - mM/2) * Cc * pi^(mM/2) * cyc)^(2/(2 - mM)) +
          rnorm(length(cyc), 0, sigma.r)
      }, split.idx, seq_len(nS)), use.names = FALSE)
      
      monotone_ok <- all(vapply(split.idx, function(idx) min(diff(grieta.B[idx])) >= 0, logical(1)))
      if (!monotone_ok) next
      
      grieta.B <- unlist(lapply(split.idx, function(idx) {
        v <- grieta.B[idx]
        v + (x$cracks[idx][1] - v[1])
      }), use.names = FALSE)
      
      x.gD.BB <- data.frame(cycles = x$cycles, cracks = grieta.B,
                            sample = x$sample, a0 = x$a0)
      
      model <- try(nlsList(cracks ~ SSParis(cycles, C, m, a0, F.cte = 1, S = 1) | sample,
                           na.action = na.omit, data = x.gD.BB), silent = TRUE)
      if (inherits(model, "try-error")) next
      if (min(coef(model)[, 1], na.rm = TRUE) <= 0 || min(coef(model)[, 2], na.rm = TRUE) <= 0) next
      
      model.nlme.B <- try(nlme(model, random = C + m ~ 1, na.action = na.omit), silent = TRUE)
      if (inherits(model.nlme.B, "try-error")) next
      
      C.B <- mean(coef(model.nlme.B)[, 1])
      m.B <- mean(coef(model.nlme.B)[, 2])
      COV.B <- var(coef(model.nlme.B))
      
      params <- MASS::mvrnorm(n = nMC, mu = c(C.B, m.B), Sigma = COV.B)
      a_0 <- min(x$a0)
      t_fallo.B <- T_func(a_0, aF, params[, 1], params[, 2])
      
      BOXP.B <- boxplot.stats(t_fallo.B, coef = 3)$stats
      t.fallo.2.B <- t_fallo.B[t_fallo.B > 0 & t_fallo.B <= max(BOXP.B)]
      if (length(t.fallo.2.B) < 3) next
      
      return(as.vector(quantile(t.fallo.2.B, probs = seq(0, 1, 0.01))))
    }
    rep(NA_real_, 101)   
  }
  
  reps    <- run_parallel(one.rep.PB, nB)
  Mat.F.B <- do.call(cbind, reps)
  
  ICB.nlme <- t(Mat.F.B)
  ic.nlme  <- apply(ICB.nlme, 2, quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  I.Bootstrap <- data.frame(low = ic.nlme[1, ], up = ic.nlme[2, ])
  return(list(Mat.F.B = Mat.F.B, I.Bootstrap = I.Bootstrap))
}
