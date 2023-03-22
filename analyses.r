

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

library(ggplot2)
library(matrixStats)
library(haven)
library(nlme)
# library(lme4)
library(data.table)
library(gimme)



# Auxiliar functions ------------------------------------------------------



getRESE.ranef.lme <- function(x){
  ss <- stack(x)
  ss$ind <- factor(as.character(ss$ind), levels = colnames(x))
  ss$.nn <- rep.int(reorder(factor(rownames(x)), x[[1]], 
        FUN = mean, sort = sort), ncol(x))
  pv <- attr(x, "postVar")
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
    rimodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = ~ 1 | ParticipantID, method = "ML", na.action = na.omit)
    rsmodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = as.formula(paste("~", iproc, "| ParticipantID")), method = "ML", na.action = na.omit)
    # rimodel <- lmer(paste0(iout, " ~ ", iproc, " + (1|ParticipantID)"), data = impactdtres, REML = FALSE)
    rsmodel2 <- lmer(paste0(iout, " ~ ", iproc, " + (", iproc,"|ParticipantID)"), data = impactdtres, REML = FALSE)
    # ranef(rsmodel)[["ParticipantID"]][,"(Intercept)"] # merMod methods
    # ranef(rsmodel)[["ParticipantID"]][,iproc] # merMod methods
    chisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["L.Ratio"]][2]
    pchisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["p-value"]][2]
    beta[ixproc, ixout] <- fixef(rsmodel)[[iproc]]
    beta10[ixproc, ixout] <- quantile(coef(rsmodel)[[iproc]], probs = c(0.1))
    beta90[ixproc, ixout] <- quantile(coef(rsmodel)[[iproc]], probs = c(0.9))
    # plot(compareFits(coef(rsmodel), coef(rimodel)), mark = fixef(rsmodel))
    plot(coef(rsmodel), panel = function(...) {panel.dotplot(..., col = 'darkgreen'); panel.abline(v = fixef(rsmodel), lty = 3, col = 'purple')})
  }
}

coef(rsmodel)
