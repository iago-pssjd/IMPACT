
#R! Preliminaries

#R!! OS dependencies

if(Sys.info()["sysname"] == "Linux"){
	data_path <- paste0("/media/", 
			    system("whoami", intern = TRUE), 
			    "/", 
			    system(paste0("ls /media/", system("whoami", intern = TRUE)), intern = TRUE), 
			    "/Dropbox/juanpablo spanish pain study/")
	data_path <- data_path[file.exists(data_path)]
	dev_lib_path <- "~/R/x86_64-pc-linux-gnu-library/dev"
} else if(Sys.info()["sysname"] == "Windows"){
	data_path <- paste0(grep("^[A-Z]:$", sub(":(.*)", ":",shell("wmic logicaldisk get name", intern = TRUE)), value = TRUE), "/NCext/PSSJD_other/IMPACT")
	data_path <- data_path[file.exists(data_path)]
	data_path <- paste0(data_path, "/")
	dev_lib_path <- .libPaths()
}

options(max.print=99999)

#R!! Libraries


library(openxlsx) # createWorkbook
library(metafor)
library(forecast)
library(psych)
library(lattice) # dotplot
library(matrixStats) # rowMins
library(haven) # read_sav
library(data.table)
# library(gimme)


#R!! Auxiliar functions

dotplot.lmerMod <- function(x, data, main = TRUE, transf = I, ...){
	xf <- fixef(x)
	x <- ranef(x, condVar = TRUE)
	for(grf in names(x)){
		for(erf in names(x[[grf]]))
			x[[grf]][, erf] <- x[[grf]][, erf] + xf[erf]
	}
	dotplot(x, data = data, main = main, transf = I, v = xf, ...)
}


dotplot.ranef.mer <- function(x, data, main = TRUE, transf = I, ...){
	prepanel.ci <- function(x, y, se, subscripts, ...){
		if (is.null(se)) return(list())
		x <- as.numeric(x)
		hw <- 1.96 * as.numeric(se[subscripts])
		list(xlim = range(transf(x - hw), transf(x + hw), finite = TRUE))
	}
	panel.ci <- function(x, y, se, subscripts, pch = 16, horizontal = TRUE, v = 0,
			     col = dot.symbol$col, lty.h = dot.line$lty, lty.v = dot.line$lty, lwd.h = dot.line$lwd, lwd.v = dot.line$lwd, 
			     col.line.h = dot.line$col, col.line.v = dot.line$col, levels.fos = unique(y), groups = NULL, 
			     ...){
		x <- as.numeric(x)
		y <- as.numeric(y)
		dot.line <- trellis.par.get("dot.line")
		dot.symbol <- trellis.par.get("dot.symbol")
		sup.symbol <- trellis.par.get("superpose.symbol")
		panel.abline(h = levels.fos, col = col.line.h, lty = lty.h, lwd = lwd.h)
		if(length(v) == 1){
			panel.abline(v = v, col = col.line.v, lty = lty.v, lwd = lwd.v)
		}else{
			panel.abline(v = v[panel.number()], col = col.line.v, lty = lty.v, lwd = lwd.v) 
		}
		if (!is.null(se)) {
			se <- as.numeric(se[subscripts])
			panel.segments(transf(x - 1.96 * se), y, transf(x + 
									1.96 * se), y, col = "black")
		}
		panel.xyplot(transf(x), y, pch = pch, col = col, ...)
	}
	f <- function(nx, ...){
		ss <- lme4:::asDf0(x, nx)
		mtit <- if (main) 
			nx
		dotplot(.nn ~ values | ind, ss, se = ss$se, prepanel = prepanel.ci, 
			panel = panel.ci, xlab = NULL, main = mtit, ...)
	}
	setNames(lapply(names(x), f, ...), names(x))
}


#R!! Data loading

impactdt <- read_sav(paste0(data_path, "Dataset_IMPACT-EMA_v2.sav"))
impactdt <- as_factor(impactdt)
setDT(impactdt)
# wb <- createWorkbook()
# wbcwb <- createWorkbook()
wbT1 <- wbmeta <- createWorkbook()



#R!! Global variables

EMA <- names(impactdt[, .SD, .SDcols = PainIntensity:Inaction])
pOUTCOMES <- c("InterferLeasure", "InterferSocial", "InterferWork") # primary outcomes
sOUTCOMES <- c("Sadness", "PainIntensity", "PainControl", "SleepDisturb", "Stress") # secondary outcomes
outcomes <- c(pOUTCOMES, sOUTCOMES)
processes <- setdiff(EMA, outcomes)



#R!! Data arranging

# remove 30 participants exhibiting no variability on at least one EMA item
impactdtres <- impactdt[!ParticipantID %in% impactdt[, lapply(.SD, \(.x) uniqueN(.x) - any(is.na(.x))), by = ParticipantID, .SDcols = EMA][, min_var := rowMins(as.matrix(.SD)), .SDcols = !c('ParticipantID')][min_var <= 1]$ParticipantID]

# impactdtres[, uniqueN(ParticipantID), by = Arm]


#R! Multilevel analyses

#R!! Global analyses

beta10 <- beta90 <- beta <- chisq <- pchisq <- matrix(nrow = length(processes), ncol = length(outcomes))
rownames(beta10) <- rownames(beta90) <- rownames(beta) <- rownames(chisq) <- rownames(pchisq) <- processes
colnames(beta10) <- colnames(beta90) <- colnames(beta) <- colnames(chisq) <- colnames(pchisq) <- outcomes



for(ixout in seq_along(outcomes)){
	for(ixproc in seq_along(processes)){
		iout <- outcomes[ixout]
		iproc <- processes[ixproc]
		rimodel <- lmer(paste0(iout, " ~ ", iproc, " + (1|ParticipantID)"), data = impactdtres, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
		rsmodel <- lmer(paste0(iout, " ~ ", iproc, " + (", iproc,"|ParticipantID)"), data = impactdtres, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
		chisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Chisq"]][2]
		pchisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Pr(>Chisq)"]][2]
		beta[ixproc, ixout] <- fixef(rsmodel)[[iproc]]
		beta10[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.1))
		beta90[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.9))
		if(Sys.info()["sysname"] == "Windows"){
			png(paste0(data_path, "graphics/", iout, "-", iproc, ".png"), bg = "transparent")
			print(dotplot(rsmodel, col = 'darkgreen', col.line.v = 'purple', lty.v = 3, scales = list(y = list(draw = FALSE)), par.settings = list(strip.border = list(alpha = 0)), strip = strip.custom(factor.levels = c("Intercept (Z score)", "Slope (beta)"), bg = "transparent"))[["ParticipantID"]])
			dev.off()
		}
	}
}

# saving matrices
sheet_list <- list("chi-square" = chisq, "chi-square p-value" = pchisq, "beta" = beta, "beta10" = beta10, "beta90" = beta90)
for(sheetname in names(sheet_list)){
	addWorksheet(wb, sheetName = sheetname)
	writeData(wb, sheet = sheetname, sheet_list[[sheetname]], rowNames = TRUE)
}


#R!! Multilevel analyses per arm


for(arm in levels(impactdtres$Arm)){
	beta10 <- beta90 <- beta <- chisq <- pchisq <- matrix(nrow = length(processes), ncol = length(outcomes))
	rownames(beta10) <- rownames(beta90) <- rownames(beta) <- rownames(chisq) <- rownames(pchisq) <- processes
	colnames(beta10) <- colnames(beta90) <- colnames(beta) <- colnames(chisq) <- colnames(pchisq) <- outcomes
	sarm <- sub("\\+", "", arm, "/")
	dir.create(paste0(data_path, "graphics/", sarm))

	for(ixout in seq_along(outcomes)){
		for(ixproc in seq_along(processes)){
			iout <- outcomes[ixout]
			iproc <- processes[ixproc]
			rimodel <- lmer(paste0(iout, " ~ ", iproc, " + (1|ParticipantID)"), data = impactdtres, REML = FALSE, control = lmerControl(optimizer = "bobyqa"), subset = Arm == arm)      # print(summary(rimodel))
			rsmodel <- lmer(paste0(iout, " ~ ", iproc, " + (", iproc,"|ParticipantID)"), data = impactdtres, REML = FALSE, control = lmerControl(optimizer = "bobyqa"), subset = Arm == arm)
			chisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Chisq"]][2]
			pchisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Pr(>Chisq)"]][2]
			beta[ixproc, ixout] <- fixef(rsmodel)[[iproc]]
			beta10[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.1))
			beta90[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.9))
			if(Sys.info()["sysname"] == "Windows"){
				png(paste0(data_path, "graphics/", sarm, "/", iout, "-", iproc, ".png"), bg = "transparent")
				print(dotplot(rsmodel, col = 'darkgreen', col.line.v = 'purple', lty.v = 3, scales = list(y = list(draw = FALSE)), par.settings = list(strip.border = list(alpha = 0)), strip = strip.custom(factor.levels = c("Intercept (Z score)", "Slope (beta)"), bg = "transparent"))[["ParticipantID"]])
				dev.off()
			}
		}
	}

	# saving matrices
	sheet_list <- list("chi-square" = chisq, "chi-square p-value" = pchisq, "beta" = beta, "beta10" = beta10, "beta90" = beta90)
	names(sheet_list) <- paste(names(sheet_list), sarm, sep = "-")
	for(sheetname in names(sheet_list)){
		addWorksheet(wb, sheetName = sheetname)
		writeData(wb, sheet = sheetname, sheet_list[[sheetname]], rowNames = TRUE)
	}
}



#R! Correlation within-between

#R!! Global correlations

cwbEMA <- statsBy(impactdtres[, .SD, .SDcols = c("ParticipantID", EMA)], group = "ParticipantID", cors = FALSE)

# saving matrices
sheet_list <- list("Within correlations" = cwbEMA$rwg, "Within correlation p-values" = cwbEMA$pwg, "Between correlations" = cwbEMA$rbg, "Between correlation p-values" = cwbEMA$pbg)
for(sheetname in names(sheet_list)){
	addWorksheet(wbcwb, sheetName = sheetname)
	writeData(wbcwb, sheet = sheetname, sheet_list[[sheetname]], rowNames = TRUE)
}


#R!! Correlation within-between per arm


for(arm in levels(impactdtres$Arm)){
	cwbEMA <- statsBy(impactdtres[Arm == arm, .SD, .SDcols = c("ParticipantID", EMA)], group = "ParticipantID")
	sarm <- sub("\\+", "", arm, "/")

	# saving matrices
	sheet_list <- list("Within correlations" = cwbEMA$rwg, "Within p-values" = cwbEMA$pwg, "Between correlations" = cwbEMA$rbg, "Between p-values" = cwbEMA$pbg)

	names(sheet_list) <- paste(names(sheet_list), sarm, sep = "-")

	for(sheetname in names(sheet_list)){
		addWorksheet(wbcwb, sheetName = sheetname)
		writeData(wbcwb, sheet = sheetname, sheet_list[[sheetname]], rowNames = TRUE)
	}

}

#R! ARIMAX models
impactIDs <- unique(impactdtres[, .(ParticipantID, Arm)])

#R!! Global models


#R!!! Only outcome models (Table 1)

T1 <- matrix(nrow = 7*3, ncol = length(outcomes))
rownames(T1) <- paste(rep(c("p", "d", "q", "s", "P", "D", "Q"), each = 3), ":", c("None", "One", "Over One"))
colnames(T1) <- outcomes
T1list <- list()

for(ixout in seq_along(outcomes)){
  # ixout <- 1
  iout <- outcomes[ixout]
    opDT <- impactdtres[, c(.(ParticipantID = ParticipantID, Time = Time, Time_Series = Time_Series), .SD), .SDcols = c(iout)]
    iDout <- matrix(nrow = nrow(impactIDs), ncol = 7)
    rownames(iDout) <- impactIDs[, ParticipantID]
    colnames(iDout) <- c("p", "d", "q", "s", "P", "D", "Q")
    for(ixid in seq_len(nrow(impactIDs))){
      # ixid <- 1
      iD <- impactIDs[ixid, ParticipantID]
      iDdt <- opDT[ParticipantID == iD]
      dt2ts <- ts(iDdt[, iout, with = FALSE])
      iDfit <- auto.arima(dt2ts)
      iDout[ixid, c("p", "d", "q", "s", "P", "D", "Q")] <- iDfit$arma[c(1, 6, 2, 5, 3, 7, 4)]
      rm(dt2ts, iDfit)
    }
    iDout2 <- apply(iDout, 2, cut, breaks = c(-Inf, 0.5, 1.5, Inf), labels = c("None", "One", "Over One"))
    for(rn in rownames(T1)){
      T1[rn, iout] <- sum(iDout2[, sub("\\s:.*$", "", rn)] == sub("^.*:\\s", "", rn))/nrow(iDout2)
    }
    T1list[[iout]] <- as.data.table(iDout)[, ARorder := interaction(p, d, q, s, P, D, Q, drop = TRUE)]
    rm(iDout, iDout2)
    gc()
}

T11 <- rbindlist(T1list, idcol = "outcome")[, .N, by = ARorder][, p:= round(N * 100/sum(N), 2)][]
addWorksheet(wbT1, sheetName = "iARIMA order frequencies")
writeData(wbT1, sheet = "iARIMA order frequencies", round(T1*100, 2), rowNames = TRUE, colNames = TRUE)
addWorksheet(wbT1, sheetName = "iARIMA pattern frequencies")
writeData(wbT1, sheet = "iARIMA pattern frequencies", T11, rowNames = FALSE, colNames = TRUE)
saveWorkbook(wbT1, paste0(data_path,"iArimaT1.xlsx"), overwrite = TRUE)

#R!!! Outcome-process models

beta_bands <- c("(-Inf, -0.3)", "[-0.3, -0.2)", "[-0.2, -0.1)", "[-0.1, 0.1]", "(0.1, 0.2]", "(0.2, 0.3]", "(0.3, +Inf)")
pooledOut <- list()
mainOut <- matrix(nrow = length(processes) * length(outcomes), ncol = 4 + length(beta_bands))
rownames(mainOut) <- paste(rep(outcomes, each = length(processes)), processes, sep = "-")
colnames(mainOut) <- c("beta", "SE", "I2", "Q", beta_bands)


for(ixout in seq_along(outcomes)){
  # ixout <- 1
  iout <- outcomes[ixout]
  for(ixproc in seq_along(processes)){
    # ixproc <- 1
    iproc <- processes[ixproc]
    wbAR <- createWorkbook()
    opDT <- impactdtres[, c(.(ParticipantID = ParticipantID, Time = Time, Time_Series = Time_Series), .SD), .SDcols = c(iout, iproc)]
    iDout <- matrix(nrow = nrow(impactIDs), ncol = 11)
    rownames(iDout) <- impactIDs[, ParticipantID]
    colnames(iDout) <- c("beta", "SE", "p", "d", "q", "s", "P", "D", "Q", "beta_reg", "SE_reg")
    for(ixid in seq_len(nrow(impactIDs))){
      # ixid <- 1
      iD <- impactIDs[ixid, ParticipantID]
      iDdt <- opDT[ParticipantID == iD]
      dt2ts <- ts(iDdt[, iout, with = FALSE])
      dt2xreg <- as.matrix(iDdt[, iproc, with = FALSE])
      iDfit <- auto.arima(dt2ts, xreg = dt2xreg)
      iDregression <- lm(reformulate(iproc, response = iout), data = iDdt)
      iDout[ixid, "beta"] <- coef(iDfit)[[iproc]]
      iDout[ixid, "SE"] <- sqrt(iDfit$var.coef[iproc, iproc])
      iDout[ixid, c("p", "d", "q", "s", "P", "D", "Q")] <- iDfit$arma[c(1, 6, 2, 5, 3, 7, 4)]
      iDout[ixid, "beta_reg"] <- coef(iDregression)[[iproc]]
      iDout[ixid, "SE_reg"] <- sqrt(vcov(iDregression)[iproc, iproc])
      rm(dt2ts, dt2xreg, iDfit)
    }
    
    # save matrix of individual models outputs
    addWorksheet(wbAR, sheetName = "Individual models output")
    writeData(wbAR, sheet = "Individual models output", iDout, rowNames = TRUE, colNames = TRUE)
    saveWorkbook(wbAR, paste0(data_path, "ARIMAX/",iout,"-",iproc,".xlsx"), overwrite = TRUE)
    
    # meta-analysis for ARIMAX models
    res.nomod <- rma(yi = iDout[,"beta"], sei = iDout[,"SE"], measure = "GEN", method = "REML")
    
    # meta-analysis for regression models
    res.reg <- rma(yi = iDout[,"beta_reg"], sei = iDout[,"SE_reg"], measure = "GEN", method = "REML")
    
    # meta-analysis for ARIMAX models adjusted per study arm
    # res.mod <- rma(yi = iDout[,"beta"], sei = iDout[,"SE"], mods = impactIDs[, Arm], measure = "GEN", method = "REML")
    
    mainOut[paste(iout, iproc, sep = "-"), "beta"] <- res.nomod$beta[1]
    mainOut[paste(iout, iproc, sep = "-"), c("SE", "I2", "Q")] <- unlist(res.nomod[c("se", "I2", "QE")])

    
    iDout <- as.data.table(iDout)
    pooledOut[[paste(iout, iproc, sep = "-")]] <- iDout
    mainOut[paste(iout, iproc, sep = "-"), beta_bands] <- iDout[, betaq := factor(fcase(between(beta, -.1, .1, incbounds = TRUE), 0, between(beta, .1, .2, incbounds = TRUE), 1, between(beta, .2, .3, incbounds = TRUE), 2, beta > .3, 3, between(beta, -.2, -.1, incbounds = TRUE), -1, between(beta, -.3, -.2, incbounds = TRUE), -2, beta < -.3, -3), levels = seq(-3, 3, 1), labels = beta_bands)][, .N, keyby = betaq][, round(100 * N/sum(N), 2)]
    
    rm(wbAR, iDout)
    gc()
  }
}

T2 <- rbindlist(pooledOut, idcol = "OPinteraction")[, c("outcome", "process") := tstrsplit(OPinteraction, split = "_", fixed = TRUE)][, .(reg_avg = mean(beta_reg), iARIMAX_avg = mean(beta), cor_reg_iARIMAX = cor(beta, beta_reg), reg_avg_SE = mean(SE_reg), iARIMAX_avg_SE = mean(SE)), by = outcome]
# forest(res.nomod, slab = paste0(impactIDs[, ParticipantID], " (", impactIDs[, Arm], ")"))


#R! Save


saveWorkbook(wb, paste0(data_path, "multilevel-models.xlsx"), overwrite = TRUE)
saveWorkbook(wbcwb, paste0(data_path, "correlations-within-between.xlsx"), overwrite = TRUE)

