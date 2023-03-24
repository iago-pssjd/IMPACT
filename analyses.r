

if(Sys.info()["sysname"] == "Linux"){
	data_path <- "~/Dropbox/juanpablo spanish pain study/"
	dev_lib_path <- "~/R/x86_64-pc-linux-gnu-library/dev"
} else if(Sys.info()["sysname"] == "Windows"){
	data_path <- paste0(grep("^[A-Z]:$", sub(":(.*)", ":",shell("wmic logicaldisk get name", intern = TRUE)), value = TRUE), "/NCext/PSSJD_other/IMPACT")
	data_path <- data_path[file.exists(data_path)]
	data_path <- paste0(data_path, "/")
	dev_lib_path <- .libPaths()
}

options(max.print=99999)

# Libraries ---------------------------------------------------------------

library(lattice)
library(ggplot2)
library(matrixStats)
library(haven)
# library(nlme)
library(lme4)
library(data.table)
library(gimme)



# Auxiliar functions ------------------------------------------------------

prepanel.ci <- function(x, y, se, subscripts, ...) {
        if (is.null(se)) 
            return(list())
        x <- as.numeric(x)
        hw <- 1.96 * as.numeric(se[subscripts])
        list(xlim = range(transf(x - hw), transf(x + hw), finite = TRUE))
}

panel.ci <- function(x, y, se, subscripts, pch = 16, horizontal = TRUE, 
        col = dot.symbol$col, lty = dot.line$lty, lwd = dot.line$lwd, 
        col.line = dot.line$col, levels.fos = unique(y), groups = NULL, 
        ...) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        dot.line <- trellis.par.get("dot.line")
        dot.symbol <- trellis.par.get("dot.symbol")
        sup.symbol <- trellis.par.get("superpose.symbol")
        panel.abline(h = levels.fos, col = col.line, lty = lty, 
            lwd = lwd)
        panel.abline(v = 0, col = col.line, lty = lty, lwd = lwd)
        if (!is.null(se)) {
            se <- as.numeric(se[subscripts])
            panel.segments(transf(x - 1.96 * se), y, transf(x + 
                1.96 * se), y, col = "black")
        }
        panel.xyplot(transf(x), y, pch = pch, ...)
    }

# https://stackoverflow.com/questions/51259346/how-to-get-names-of-dot-dot-dot-arguments-in-r
f <- function(nx, main = TRUE, transf = I, ...) {
  dots <- match.call(expand.dots = FALSE)$...

          ss <- lme4:::asDf0(x, nx)
        mtit <- if (main) 
            nx
        dotplot(.nn ~ values | ind, ss, se = ss$se, prepanel = prepanel.ci, 
            panel = panel.ci, xlab = NULL, main = mtit, ...)
    }


# Data loading --------------------------------------------------------------

impactdt <- read_sav(paste0(data_path, "Dataset_IMPACT-EMA_v2.sav"))
impactdt <- as_factor(impactdt)
setDT(impactdt)



# Global variables --------------------------------------------------------

EMA <- names(impactdt[, .SD, .SDcols = PainIntensity:Inaction])
pOUTCOMES <- c("InterferLeasure", "InterferSocial", "InterferWork")
sOUTCOMES <- c("Sadness", "PainIntensity", "PainControl", "SleepDisturb", "Stress")
outcomes <- c(pOUTCOMES, sOUTCOMES)
processes <- setdiff(EMA, outcomes)
superpose.symbol <- trellis.par.get("superpose.symbol")
superpose.symbol$fill[5] <- 'lightgreen'
trellis.par.set("superpose.symbol", superpose.symbol)


# Data arranging ----------------------------------------------

impactdtres <- impactdt[!ParticipantID %in% impactdt[, lapply(.SD, \(.x) uniqueN(.x) - any(is.na(.x))), by = ParticipantID, .SDcols = EMA][, min_var := rowMins(as.matrix(.SD)), .SDcols = !c('ParticipantID')][min_var <= 1]$ParticipantID]

beta10 <- beta90 <- beta <- chisq <- pchisq <- matrix(nrow = length(processes), ncol = length(outcomes))
rownames(beta10) <- rownames(beta90) <- rownames(beta) <- rownames(chisq) <- rownames(pchisq) <- processes
colnames(beta10) <- colnames(beta90) <- colnames(beta) <- colnames(chisq) <- colnames(pchisq) <- outcomes

# Multilevel analyses -----------------------------------------------------

for(ixout in seq_along(outcomes)){
  for(ixproc in seq_along(processes)){
    iout <- outcomes[ixout]
    iproc <- processes[ixproc]
    # rimodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = ~ 1 | ParticipantID, method = "ML", na.action = na.omit)
    # rsmodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = as.formula(paste("~", iproc, "| ParticipantID")), method = "ML", na.action = na.omit)
    # chisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["L.Ratio"]][2]
    # pchisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["p-value"]][2]
    # beta10[ixproc, ixout] <- quantile(coef(rsmodel)[[iproc]], probs = c(0.1))
    # beta90[ixproc, ixout] <- quantile(coef(rsmodel)[[iproc]], probs = c(0.9))
    # plot(coef(rsmodel), panel = function(...) {panel.dotplot(..., col = 'darkgreen'); panel.abline(v = fixef(rsmodel), lty = 3, col = 'purple')})

    rimodel <- lmer(paste0(iout, " ~ ", iproc, " + (1|ParticipantID)"), data = impactdtres, REML = FALSE)
    rsmodel <- lmer(paste0(iout, " ~ ", iproc, " + (", iproc,"|ParticipantID)"), data = impactdtres, REML = FALSE)
    # ranef(rsmodel)[["ParticipantID"]][,"(Intercept)"] # merMod methods
    # ranef(rsmodel)[["ParticipantID"]][,iproc] # merMod methods
    chisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Chisq"]][2]
    pchisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Pr(>Chisq)"]][2]
    beta[ixproc, ixout] <- fixef(rsmodel)[[iproc]]
    beta10[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.1))
    beta90[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.9))
    # plot(compareFits(coef(rsmodel), coef(rimodel)), mark = fixef(rsmodel))
    dotplot(ranef(rsmodel, condVar = TRUE), col = 'darkgreen')
  }
}

coef(rsmodel)
