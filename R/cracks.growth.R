#' @title Fatigue crack growth in reliability analysis
#'
#' @description It provides the lifetime distribution of metallic materials that fail due to the crack
#' growth produced by mechanical fatigue efforts. The crack growth trends are fitted by
#' linear or nonlinear mixed effects regression models in order to make predictions
#' about the material lifetime. The lifetime is defined as the time passed before the material
#' does not meet the specification requirements and it is conditioned by a critical crack
#' length that induces the material failure. Three different methods can be applied to estimate the
#' fatigue lifetime distribution: "SEP-lme_bkde" and "SEP-lme_kde" are nonparametric
#' while "PB-nlme" corresponds to the parametric approach proposed by Pinheiro and
#' Bates (2000).
#'
#' @param x Matrix or data frame composed of three columns: times or number of cycles, crack
#' lengths and specimen number.
#' @param aF Critical crack length for which the material failure is produced.
#' @param T_c Censoring time or frequency.
#' @param method A string of characters: "SEP-lme_bkde" (default methodology) indicates that
#' a mixed effects linear regression model is applied to crack growth data and the lifetime
#' density is estimated by bkde method. "SEP-lme_kde" indicates that a mixed effects linear
#' regression model is applied to crack growth data and the lifetime distribution is estimated
#' by kde method. "PB-nlme" indicates that a mixed effects nonlinear regression model is applied
#' to crack growth data and the lifetime the parameters are estimated maximum likelihood methodologies,
#' and the lifetime distribution by Monte Carlo.
#' @param nMC Number of Monte Carlo estimates, by default 5000.
#' @param nBKDE Number of bkde estimates, by default 5000.
#' @param nKDE Number of kde estimates, by default 5000.
#' @param seed Seed for reproducibility.
#'
#' @details This function provides a simultaneous fitting of crack growth data corresponding to
#' different specimens when these are subjected to mechanical fatigue efforts. For this purpose,
#' mixed effects linear models (lme) with B-spline smoothing are applied. Since the failure is
#' defined at a specific critical crack length, predictions of material lifetime are obtained
#' assuming the linearized Paris-Erdogan law and the material lifetime distribution is estimated.
#' There are available three different techniques to estimate the lifetime distribution: the binned
#' kernel density estimate (bkde), the kernel estimator for the distribution function (kde) computed
#' by Quintela del Rio and Estevez-Perez (2012), and in addition the parametric method proposed by
#' Pinheiro and Bates (2000) based on mixed effects nonlinear regression (nlme), maximum likelihood
#' and Monte Carlo simulation.
#'
#' @return
#' \describe{
#'   \item{\code{data}}{Data frame with the data corresponding to number of cycles, crack length,
#' and sample.}
#'   \item{\code{a.F}}{Critical crack length.}
#'   \item{\code{Tc}}{Censoring time.}
#'   \item{\code{param}}{Data frame with the estimates of Paris law parameters: C and m}
#'   \item{\code{crack.est}}{Data frame with time, crack growth estimates, and corresponding sample
#' or specimen.}
#'   \item{\code{sigma}}{Residual standard deviation.}
#'   \item{\code{residuals}}{Residuals resulting from the crack length fitting.}
#'   \item{\code{crack.pred}}{Data frame with time, crack growth predictions out of the experimental
#' time interval, and corresponding sample or specimen.}
#'   \item{\code{F.emp}}{Data frame with the empirical lifetime distribution and the corresponding
#' time: time, Fe.}
#'   \item{\code{bw}}{Bandwidth used in bkde and kde methods.}
#'   \item{\code{F.est}}{Data frame with the estimated lifetime distribution and the corresponding
#' time: time, F.}
#'   \item{\code{nBKDE}}{Number of bkde estimates.}
#'   \item{\code{nKDE}}{Number of kde estimates.}
#'   \item{\code{nMC}}{Number of Monte Carlo estimates.}
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
#' @import KernSmooth
#' @import kerdiest
#' @import nlme
#' @importFrom grDevices boxplot.stats colors
#' @importFrom graphics abline lines points
#' @importFrom stats coef coefficients lm na.omit predict quantile resid rnorm selfStart spline var ave fitted
#'
#'
#' @export
#' @examples
#' 
#' ## Using the Alea.A dataset
#' data(Alea.A)
#' x <-  Alea.A
#' ## Critical crack length
#' aF <- 1.6
#' ## Censoring time
#' T_c <- 0.12
#' ## cracks.growth function applied to Alea.A data
#' cg <- cracks.growth (x, aF, T_c, method = c("SEP-lme_bkde", "SEP-lme_kde",
#'                     "PB-nlme"), nBKDE = 1000, nKDE = 1000, nMC = 5000)
#' ## cracks.growth values using the "SEP-lme_bkde" by default method.
#' names(cg)
#' # [1]	 "data"       "a.F"        "Tc"         "param"      "crack.est"
#' # [6] 	"sigma"      "residuals"  "crack.pred" "F.emp"      "bw"
#' #[11]	 "F.est"      "nBKDE"
cracks.growth <- function(x, aF, T_c,
                          method = c("SEP-lme_bkde", "SEP-lme_kde", "PB-nlme"),
                          nMC = 5000, nBKDE = 1000, nKDE = 1000, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  method <- match.arg(method)          
  
  a0 <- stats::ave(x[[2]], x[[3]], FUN = function(v) v[1])
  x[[4]] <- a0
  colnames(x) <- c("cycles", "cracks", "sample", "a0")
  samples <- unique(x$sample)
  nS <- length(samples)
  split.idx <- split(seq_len(nrow(x)), x$sample)
  
  Paris.F <- function(cycles, C, m, a0, F.cte = 1, S = 1)
    (a0^(1 - m/2) + (1 - m/2) * C * F.cte^m * S^m * pi^(m/2) * cycles)^(2/(2 - m))
  
  DerivInit <- function(data, cycles, y) {
    xx <- data$cycles
    y  <- data$cracks
    deriv.ii <- sfsmisc::D1D2(xx[-1], y[-1])
    data.frame(cycles = xx[-1], derivada = deriv.ii$D1, grieta.media = y[-1])
  }
  
  ParisInit <- function(mCall, data, LHS) {
    dd  <- DerivInit(data = data, cycles = mCall[["cycles"]], y = LHS)
    fit <- lm(log(derivada) ~ log(sqrt(pi * grieta.media)), data = dd)
    res <- c(exp(coefficients(fit)[1]), coefficients(fit)[2])
    names(res) <- mCall[c("C", "m")]
    res
  }
  
  SSParis <- selfStart(Paris.F, initial = ParisInit, parameters = c("C", "m"))
  
  T_func <- function(a_0, aF, C, m)
    (1 - (aF / a_0)^(-(m - 2) / 2)) / ((a_0^(m/2 - 1)) * (m/2 - 1) * C * pi^(m/2))
  
  
  if (method == "PB-nlme") {
    
    model      <- nlsList(cracks ~ SSParis(cycles, C, m, a0, F.cte = 1, S = 1) | sample,
                          na.action = na.omit, data = x)
    model.nlme <- nlme(model, random = C + m ~ 1, na.action = na.omit)
    
    COEF    <- coef(model.nlme)
    CC      <- mean(COEF[, 1]); mm <- mean(COEF[, 2])
    COV     <- var(COEF)
    sigma.r <- model.nlme$sigma
    resid1  <- resid(model.nlme)
    
    ord <- order(as.numeric(row.names(COEF)))
    ORDER1.COEF <- data.frame(C = COEF[, 1][ord], m = COEF[, 2][ord])
    
    grieta.P.nlme <- Map(function(idx, i) {
      C <- ORDER1.COEF[i, 1]; m <- ORDER1.COEF[i, 2]
      t <- x$cycles[idx]; a00 <- x$a0[idx][1]
      (a00^(1 - m/2) + (1 - m/2) * C * pi^(m/2) * t)^(2/(2 - m))
    }, split.idx, seq_along(split.idx))
    
    Grieta.P.nlme <- data.frame(time = x$cycles,
                                growt.est = unlist(grieta.P.nlme, use.names = FALSE),
                                sample = x$sample)
    
    T.M <- vapply(seq_len(nS), function(i) {
      a00 <- x$a0[split.idx[[i]]][1]
      T_func(a00, aF, ORDER1.COEF[i, 1], ORDER1.COEF[i, 2])
    }, numeric(1))
    
    L2 <- seq(min(x$cycles), max(T.M), length.out = 1001)
    
    pred.P.nlme <- lapply(seq_len(nS), function(i) {
      C <- ORDER1.COEF[i, 1]; m <- ORDER1.COEF[i, 2]
      a00 <- x$a0[split.idx[[i]]][1]
      g <- (a00^(1 - m/2) + (1 - m/2) * C * pi^(m/2) * L2)^(2/(2 - m))
      data.frame(tt = L2, g.P.nlme = g)
    })
    
    Grieta.P.nlme.pred <- do.call(rbind, Map(function(p, i)
      data.frame(time = p$tt, growt.pred = p$g.P.nlme, sample = i),
      pred.P.nlme, seq_len(nS)))
    
    T.FALLO <- vapply(pred.P.nlme, function(p) p$tt[which.min(abs(p$g.P.nlme - aF))], numeric(1))
    
    params  <- MASS::mvrnorm(n = nMC, mu = c(CC, mm), Sigma = COV)
    valid   <- params[, 1] > 0 & params[, 2] > 0
    a_0     <- min(x$cracks)
    t_fallo <- numeric(nMC)
    t_fallo[valid] <- T_func(a_0, aF, params[valid, 1], params[valid, 2])
    
    BOXP      <- boxplot.stats(t_fallo, coef = 3)$stats
    t.fallo.2 <- t_fallo[t_fallo > 0 & t_fallo <= max(BOXP)]
    F1        <- as.vector(quantile(t.fallo.2, probs = seq(0, 1, 0.01)))
    
    Fe.1      <- round(((seq_len(nS) - 0.5) / nS), 3)
    T.FALLO.1 <- sort(T.FALLO)
    
    return(list(
      data = x[, -4], a.F = aF, Tc = T_c, param = ORDER1.COEF,
      crack.est = Grieta.P.nlme, sigma = sigma.r, residuals = resid1,
      crack.pred = Grieta.P.nlme.pred,
      F.emp = data.frame(time = T.FALLO.1, Fe = Fe.1),
      F.est = data.frame(time = F1, F = seq(0, 1, 0.01)),
      nMC = nMC
    ))
  }
  
  x.gamm <- data.frame(cycles = x$cycles, cracks = x$cracks,
                       sample = as.factor(x$sample), a0 = x$a0)
  
  ajuste.gam <- mgcv::gam((cracks - a0) ~ s(cycles, by = sample, k = 6, bs = "cr") + sample,
                          data = x.gamm, method = "REML")
  
  resid2    <- resid(ajuste.gam)
  sigma.gam <- sqrt(ajuste.gam$sig2)
  pred.lme  <- stats::fitted(ajuste.gam)
  be.pred   <- data.frame(cycles = x$cycles, pred = as.vector(pred.lme),
                          cracks = x$cracks, sample = x$sample)
  
  split.idx2 <- split(seq_len(nrow(be.pred)), be.pred$sample)
  
  coef.lme <- lapply(split.idx2, function(idx) {
    sp   <- pspline::smooth.Pspline(be.pred$cycles[idx], be.pred$pred[idx],
                                    norder = 2, method = 3)
    dpr  <- as.vector(predict(sp, be.pred$cycles[idx], nderiv = 1))
    lin  <- lm(log(dpr[-1]) ~ log(sqrt(pi * be.pred$cracks[idx][-1])))
    coefficients(lin)
  })
  
  COEF.lme <- data.frame(
    C = vapply(coef.lme, function(cc) exp(cc[1]), numeric(1)),
    m = vapply(coef.lme, function(cc) cc[2], numeric(1))
  )
  
  grieta.P.lme <- Map(function(idx, i) {
    Ce <- COEF.lme[i, 1]; me <- COEF.lme[i, 2]
    t   <- be.pred$cycles[idx]; aoo <- x$a0[split.idx2[[i]]][1]
    (aoo^(1 - me/2) + (1 - me/2) * Ce * pi^(me/2) * t)^(2/(2 - me))
  }, split.idx2, seq_along(split.idx2))
  
  Grieta.P.lme <- data.frame(time = x$cycles,
                             growt.est = unlist(grieta.P.lme, use.names = FALSE),
                             sample = x$sample)
  
  TT.M <- vapply(seq_len(nS), function(i) {
    a000 <- x$a0[split.idx2[[i]]][1]
    T_func(a000, aF, COEF.lme[i, 1], COEF.lme[i, 2])
  }, numeric(1))
  
  LL2 <- seq(min(x$cycles), max(TT.M), length.out = 1001)
  
  pred.P.lme <- lapply(seq_len(nS), function(i) {
    C   <- COEF.lme[i, 1]; m <- COEF.lme[i, 2]
    a0v <- x$a0[split.idx2[[i]]][1]
    g   <- (a0v^(1 - m/2) + (1 - m/2) * C * pi^(m/2) * LL2)^(2/(2 - m))
    data.frame(ttt = LL2, g.P.lme = g)
  })
  
  Grieta.P.lme.pred <- do.call(rbind, Map(function(p, i)
    data.frame(time = p$ttt, growt.pred = p$g.P.lme, sample = i),
    pred.P.lme, seq_len(nS)))
  
  TT.FALLO <- vapply(pred.P.lme, function(p) p$ttt[which.min(abs(p$g.P.lme - aF))], numeric(1))
  
  Fe.1x           <- round(((seq_len(nS) - 0.5) / nS), 3)
  TT.FALLO.sorted <- sort(TT.FALLO)
  
  if (method == "SEP-lme_bkde") {
    
    hop      <- KernSmooth::dpik(TT.FALLO)
    dens.lme <- KernSmooth::bkde(TT.FALLO, gridsize = nBKDE,
                                 range.x = range(TT.FALLO), bandwidth = hop)
    bin      <- diff(range(TT.FALLO)) / nBKDE
    F.xwj    <- cumsum(abs(dens.lme$y)) * bin
    keep     <- !duplicated(F.xwj)
    xwj      <- dens.lme$x[keep]
    F.xwj    <- unique(F.xwj) / max(F.xwj)
    
    sp <- spline(F.xwj, xwj, n = 101, method = "fmm", xmin = 0, xmax = 100,
                 xout = seq(0, 1, length.out = 101), ties = mean)
    
    return(list(
      data = x[, -4], a.F = aF, Tc = T_c, param = COEF.lme,
      crack.est = Grieta.P.lme, sigma = sigma.gam,
      residuals = as.numeric(resid2), crack.pred = Grieta.P.lme.pred,
      F.emp = data.frame(time = TT.FALLO.sorted, Fe = Fe.1x),
      bw = hop, F.est = data.frame(time = sp$y, F = sp$x),
      nBKDE = nBKDE
    ))
    
  } else { 
    
    range.F <- range(TT.FALLO)
    FFF <- kerdiest::kde(type_kernel = "e", vec_data = TT.FALLO,
                         y = seq(range.F[1], range.F[2], length.out = nKDE))
    
    sp <- spline(FFF$Estimated_values, FFF$grid, n = 101, method = "fmm",
                 xmin = 0, xmax = 100, xout = seq(0, 1, length.out = 101), ties = mean)
    
    return(list(
      data = x[, -4], a.F = aF, Tc = T_c, param = COEF.lme,
      crack.est = Grieta.P.lme, sigma = sigma.gam,
      residuals = as.numeric(resid2), crack.pred = Grieta.P.lme.pred,
      F.emp = data.frame(time = TT.FALLO.sorted, Fe = Fe.1x),
      bw = FFF$bw, F.est = data.frame(time = sp$y, F = sp$x),
      nKDE = nKDE
    ))
  }
}
