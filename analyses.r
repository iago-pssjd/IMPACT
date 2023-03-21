

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

library(haven)
# library(nlme)
library(lme4)
library(data.table)
library(gimme)



# Global variables --------------------------------------------------------

EMA <- names(impactdt[, .SD, .SDcols = PainIntensity:Inaction])
pOUTCOMES <- c("InterferLeasure", "InterferSocial", "InterferWork")
sOUTCOMES <- c("Sadness", "PainIntensity", "PainControl", "SleepDisturb", "Stress")
outcomes <- c(pOUTCOMES, sOUTCOMES)
processes <- setdiff(EMA, outcomes)

# Data loading --------------------------------------------------------------

impactdt <- read_sav(paste0(data_path, "Dataset_IMPACT-EMA_v2.sav"))





# Data arranging ----------------------------------------------

impactdt <- as_factor(impactdt)
setDT(impactdt)
impactdtres <- impactdt[!ParticipantID %in% impactdt[, lapply(.SD, \(.x) uniqueN(.x) - any(is.na(.x))), by = ParticipantID, .SDcols = EMA][, min_var := rowMins(as.matrix(.SD)), .SDcols = !c('ParticipantID')][min_var <= 1]$ParticipantID]

chisq <- pchisq <- matrix(nrow = length(processes), ncol = length(outcomes))
rownames(chisq) <- rownames(pchisq) <- processes
colnames(chisq) <- colnames(pchisq) <- outcomes

# Multilevel analyses -----------------------------------------------------

for(ixout in seq_along(outcomes)){
  for(ixproc in seq_along(processes)){
    iout <- outcomes[ixout]
    iproc <- processes[ixproc]
    # rimodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = ~ 1 | ParticipantID, method = "ML", na.action = na.omit)
    # rsmodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = as.formula(paste("~", iproc, "| ParticipantID")), method = "ML", na.action = na.omit)
    rimodel <- lmer(paste0(iout, " ~ ", iproc, " + (1|ParticipantID)"), data = impactdtres, REML = FALSE)
    rsmodel <- lmer(paste0(iout, " ~ ", iproc, " + (", iproc,"|ParticipantID)"), data = impactdtres, REML = FALSE)
    # ranef(rsmodel)[["ParticipantID"]][,"(Intercept)"]
    # ranef(rsmodel)[["ParticipantID"]][,iproc]
    chisq[ixproc, ixout] <- anova(rsmodel, rimodel, test = "Chisq")[["Chisq"]][2]
    pchisq[ixproc, ixout] <- anova(rsmodel, rimodel, test = "Chisq")[["Pr(>Chisq)"]][2]
  }
}
