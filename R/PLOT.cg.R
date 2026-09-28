#' @title Fatigue Crack Growth in Reliability plots
#'
#' @description It provides graphical outputs composed of the trends corresponding to the
#' crack length growth due to mechanical fatigue, the crack length estimates by the models,
#' crack length predictions, and lifetime distribution estimates.
#'
#' @param x cracks.growth object.
#'
#' @details Specifically, the following graphs are provided: exploratory dataset graph,
#' plot with the crack length estimates and predictions, residuals graph, empirical
#' and estimated lifetime distribution plot obtained by SEP-lme_bkde, SEP-lme_kde or
#' PB-nlme methods.
#'
#' @return
#' \describe{ Return the following values:
#'   \item{\code{plot.data}}{Exploratory chart.}
#'   \item{\code{plot.pred}}{Plot for fatigue lifetimes estimates and predictions.}
#'   \item{\code{plot.F}}{Plot for the empirical distribution and lifetimes
#'   distribution estimates of fatigue lifetimes.}
#'   \item{\code{plot.resid}}{Residuals chart.}
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
#' @importFrom grDevices boxplot.stats colors
#' @importFrom graphics abline lines points
#' @importFrom stats coef coefficients lm na.omit predict quantile resid rnorm selfStart spline var
#'
#' @export
#'
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
#' cg <- cracks.growth(x, aF, T_c,
#'                    method = c("SEP-lme_bkde", "SEP-lme_kde", "PB-nlme"),
#'                    nBKDE = 1000, nKDE = 1000, nMC = 5000)
#' ## PLOT.cg applied to cg object.
#' PLOT <- PLOT.cg(cg)
#' names(PLOT)
#' ## [1]  "plot.data"  "plot.pred"  "plot.F"     "plot.resid"
#' ## Exploratory chart for the Alea.A dataset
#' PLOT$plot.data(main = "Plot:  crack growth", xlab = "million cycles",
#'                ylab = "cracks(inches)",  cex.lab=1.8,
#'                cex.main = 2)
#' text(0.02, aF + 0.05, "Failure", cex = 1.8)
#' text(0.095, 0.95, "Censoring time->", cex = 1.5)
#' ## Plot for fatigue lifetimes estimates and predictions.
#' PLOT$plot.pred(xlab = "million cycles", ylab = "cracks(inches)",
#'                main = "Plot: crack growth, estimation and prediction\n failure times (red)",
#'                cex.lab = 1.8, cex.main = 1.5)
#' text(0.02,aF+0.05, "Failure", cex = 1.8)
#' text(0.085,0.95, "Censoring time->", cex = 1.5)
#' ## Plot for the empirical distribution and lifetimes distribution estimates
#' ## of  fatigue lifetimes
#' PLOT$plot.F(main = "Plot: distributions of failure times",
#'             xlab = "million cycles", ylab = "probability",
#'             cex.lab = 1.7, cex.main=2)
#' text(0.14, 0.1, "<-Censoring time", cex = 1.5)
#' legend("topleft", c("Empirical", "Estimated"), col = c("blue","black"),
#'        pch=c(20,20), cex=1.5, bty="n")
#' ## Residuals chart.
#' PLOT$plot.resid(main = "Plot: residual", xlab = "fitted", ylab = "residuals",
#'                 cex = 1.5, col = "blue", cex.lab = 1.7, cex.main = 2)
PLOT.cg <- function (x)
{
    if (x$data[1, 2] <= x$data[length(x$data[, 2][x$data[, 3] ==
        unique(x$data[, 3])[1]]), 2]) {
        m = 1
    }
    else {
      message("cracks are not growing")
    }
    plot.data = function(y = x$data, z = x$a.F, u = x$Tc, ...) {
        COL = c(1:7, colors()[82:150])
        plot.data = plot(subset(y[, c(1, 2)], y[, 3] == unique(y[,
            3])[1]), type = "b", xlim = c(min(y[, 1]), max(y[,
            1])), ylim = c(min(y[, 2]), max(y[, 2])), col = 1,
            las = 1, pch = 20, cex = 1.5, ...)
        for (i in 2:length(unique(y[, 3]))) points(subset(y[,
            c(1, 2)], y[, 3] == unique(y[, 3])[i]), type = "b",
            col = COL[i], pch = 20, cex = 1.5)
        abline(h = z, col = "gray30", lwd = 3, lty = 2)
        abline(v = u, col = "gray30", lwd = 3, lty = 2)
    }
    plot.pred = function(y = x$data, z = x$a.F, u = x$Tc, v = x$crack.pred,
        w = x$F.emp, ...) {
        COL = c(1:7, colors()[82:150])
        plot.pred = plot(subset(y[, c(1, 2)], y[, 3] == unique(y[,
            3])[1]), type = "p", xlim = c(min(y[, 1]), max(w[,
            1])), ylim = c(min(y[, 2]), max(y[, 2])), col = 1,
            las = 1, pch = 20, cex = 1.5, ...)
        abline(h = z, col = "gray30", lwd = 3, lty = 2)
        abline(v = u, col = "gray30", lwd = 3, lty = 2)
        for (i in 1:length(unique(v[, 3]))) lines(subset(v[,
            c(1, 2)], v[, 3] == unique(v[, 3])[i]), col = COL[i],
            lwd = 2)
        for (i in 1:length(unique(y[, 3]))) points(subset(y[,
            c(1, 2)], y[, 3] == unique(y[, 3])[i]), type = "p",
            pch = 20, cex = 1.5)
        points(w[, 1], rep(z, length(w[, 1])), cex = 2, col = 2,
            lwd = 2)
    }
    plot.F = function(y = x$F.est, u = x$Tc, w = x$F.emp, ...) {
        plot.F = plot(y, xlim = c(min(y[, 1]), max(w[, 1])),
            col = 1, las = 1, pch = 20, cex = 1.5, ...)
        abline(v = u, col = "gray30", lwd = 3, lty = 2)
        points(w, col = 4, pch = 20, cex = 2)
    }
    plot.resid = function(y = x$crack.est, u = x$residuals, ...) {
        plot.resid = plot(y[, 2], u, pch = 20, ...)
        abline(h = 0, col = "gray30", lwd = 2)
    }
    list(plot.data = plot.data, plot.pred = plot.pred, plot.F = plot.F,
        plot.resid = plot.resid)
}
