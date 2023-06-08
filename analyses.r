
#R! Preliminaries

#R!! OS dependencies

if(Sys.info()["sysname"] == "Linux"){
	data_path <- paste0("/media/", 
			    system("whoami", intern = TRUE), 
			    "/", 
			    system(paste0("ls /media/", system("whoami", intern = TRUE)), intern = TRUE), 
			    "/Dropbox/juanpablo spanish pain study/")
	data_NCpath <- paste0("/media/", 
			    system("whoami", intern = TRUE), 
			    "/", 
			    system(paste0("ls /media/", system("whoami", intern = TRUE)), intern = TRUE), 
			    "/ExtNC/PSSJD_other/IMPACT/")
	data_path <- data_path[file.exists(data_path)]
	data_NCpath <- data_NCpath[file.exists(data_NCpath)]
	dev_lib_path <- "~/R/x86_64-pc-linux-gnu-library/dev"
} else if(Sys.info()["sysname"] == "Windows"){
	data_path <- paste0(grep("^[A-Z]:$", sub(":(.*)", ":",shell("wmic logicaldisk get name", intern = TRUE)), value = TRUE), "/NCext/PSSJD_other/IMPACT")
	data_path <- data_path[file.exists(data_path)]
	data_NCpath <- data_path <- paste0(data_path, "/")
	dev_lib_path <- .libPaths()
}

options(max.print=99999)

#R!! Libraries


library(openxlsx) # createWorkbook
library(cluster)
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

respondersdt <- read_sav(paste0(data_NCpath, "Responders (wide form).sav"))
respondersdt <- as_factor(respondersdt)
setDT(respondersdt)

impactdt <- merge(impactdt, respondersdt, by.x = "ParticipantID", by.y = "ID", all.x = TRUE, all.y = FALSE)

# wb <- createWorkbook()
# wbcwb <- createWorkbook()
wbmeta <- createWorkbook()



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
	sarm <- sub("\\+", "", arm)
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
	sarm <- sub("\\+", "", arm)

	# saving matrices
	sheet_list <- list("Within correlations" = cwbEMA$rwg, "Within p-values" = cwbEMA$pwg, "Between correlations" = cwbEMA$rbg, "Between p-values" = cwbEMA$pbg)

	names(sheet_list) <- paste(names(sheet_list), sarm, sep = "-")

	for(sheetname in names(sheet_list)){
		addWorksheet(wbcwb, sheetName = sheetname)
		writeData(wbcwb, sheet = sheetname, sheet_list[[sheetname]], rowNames = TRUE)
	}

}

#R! iARIMAX models

impactIDs <- unique(impactdtres[, .(ParticipantID, Arm, Wave, N_sessions, Completers, Years_diagnosis, Age, Gender, Depression_diagnosis, (.SD)), .SDcols = BPI_Pre:Responder_FollowUp])[, `:=` (age = droplevels(cut(Age, breaks = c(0,50,55,60,+Inf), right = TRUE)), 
																							    Gender = droplevels(Gender),
																							    years_diagnosis = droplevels(cut(Years_diagnosis, breaks = c(0,5,10,20,+Inf), right = TRUE)))]
ARIMA_parameters <- c("p", "d", "q", "s", "P", "D", "Q")
beta_bands <- c("(-Inf, -0.3)", "[-0.3, -0.2)", "[-0.2, -0.1)", "[-0.1, 0.1]", "(0.1, 0.2]", "(0.2, 0.3]", "(0.3, +Inf)")
# moderators <- c("Gender", "age", "Depression_diagnosis", "BPI_Pre", "NRS_Pre", "DASS21_Ans_Pre", "DASS21_Dep_Pre", "DASS21_St_Pre", "PCS_Pre", "CPAQ_Pre", "BADSSF_Pre", "PIPS_Pre")
# dput(attr(res.mod$b, "dimnames")[[1]])
# moderators <- c("intrcpt", "GenderFemale", "age(40,45]", "age(45,50]", "age(50,55]", "age(55,60]", "age(60,65]", "age(65,70]", "Depression_diagnosisDepression", "BPI_Pre", "NRS_Pre", "DASS21_Ans_Pre", "DASS21_Dep_Pre", "DASS21_St_Pre", "PCS_Pre", "CPAQ_Pre", "BADSSF_Pre", "PIPS_Pre")
# beta_moderators <- paste("betamod", moderators, sep = "_")
# SE_moderators <- paste("SEmod", moderators, sep = "_")
# p_moderators <- paste("pmod", moderators, sep = "_")
moderators <- c("Gender", "age", "Depression_diagnosis", "years_diagnosis")


T1 <- matrix(data = NA, nrow = length(ARIMA_parameters)*3, ncol = length(outcomes))
rownames(T1) <- paste(rep(ARIMA_parameters, each = 3), ":", c("None", "One", "Over One"))
colnames(T1) <- outcomes
T1list <- list()

mainOut <- matrix(data = NA, nrow = length(processes) * length(outcomes), ncol = 6 + length(beta_bands) + 5 * 2)
rownames(mainOut) <- paste(rep(outcomes, each = length(processes)), processes, sep = "-")
colnames(mainOut) <- c("beta", "SE", "p", "I2", "Q", "pQ", beta_bands, "I2mod", "QEmod", "pQEmod", "QMmod", "pQMmod", "I2arm", "QEarm", "pQEarm", "QMarm", "pQMarm")
pooledOut <- list()
modOut <- list()
modArmOut <- list()

sarm <- sub("\\+", "", levels(impactdtres$Arm))

moderatorsArm <- list(moderators, c(moderators, "Completers"), c(moderators, "Completers"))
names(moderatorsArm) <- sarm
T1arm <- list(T1, T1, T1)
names(T1arm) <- sarm
T11arm <- vector(mode = "list", length = length(sarm))
names(T11arm) <- sarm

mainOutArm <- list(mainOut[, -grep("arm", colnames(mainOut))], mainOut[, -grep("arm", colnames(mainOut))], mainOut[, -grep("arm", colnames(mainOut))])
names(mainOutArm) <- sarm
T2arm <- vector(mode = "list", length = length(sarm))
names(T2arm) <- sarm
modOutArm <- vector(mode = "list", length = length(sarm))
names(modOutArm) <- sarm




#R!! Only outcome models (Table 1)

for(ixout in seq_along(outcomes)){
	# ixout <- 1
	iout <- outcomes[ixout]
	iDout <- matrix(data = NA, nrow = nrow(impactIDs), ncol = 7)
	rownames(iDout) <- impactIDs[, ParticipantID]
	colnames(iDout) <- ARIMA_parameters
	for(ixid in seq_len(nrow(impactIDs))){
		# ixid <- 1
		iD <- impactIDs[ixid, ParticipantID]
		iDdt <- impactdtres[ParticipantID == iD]
		dt2ts <- ts(iDdt[, iout, with = FALSE])
		iDfit <- try(auto.arima(dt2ts), silent = TRUE)
		if(!inherits(iDfit, what = "try-error", which = FALSE)){
			iDout[ixid, ARIMA_parameters] <- iDfit$arma[c(1, 6, 2, 5, 3, 7, 4)]
		}
		rm(dt2ts, iDfit)
	}
	iDout2 <- apply(iDout, 2, cut, breaks = c(-Inf, 0.5, 1.5, Inf), labels = c("None", "One", "Over One"))
	for(rn in rownames(T1)){
		T1[rn, iout] <- sum(iDout2[, sub("\\s:.*$", "", rn)] == sub("^.*:\\s", "", rn))/nrow(iDout2)
	}
	T1list[[iout]] <- as.data.table(iDout, keep.rownames = "ParticipantID")[, `:=` (ARorder = interaction(p, d, q, s, P, D, Q, drop = TRUE), 
											Arm = impactIDs[, Arm])]

	for(arm in levels(impactdtres$Arm)){
		sarm <- sub("\\+", "", arm)
		iDoutArm <- iDout[impactIDs[Arm == arm, ParticipantID],]
		iDout2 <- apply(iDoutArm, 2, cut, breaks = c(-Inf, 0.5, 1.5, Inf), labels = c("None", "One", "Over One"))
		for(rn in rownames(T1arm[[sarm]])){
			T1arm[[sarm]][rn, iout] <- sum(iDout2[, sub("\\s:.*$", "", rn)] == sub("^.*:\\s", "", rn))/nrow(iDout2)
		}
	}

	rm(iDout, iDout2)
	gc()
}

T11 <- rbindlist(T1list, idcol = "outcome")[, .N, by = ARorder][, p := round(N * 100/sum(N), 2)][]

for(arm in levels(impactdtres$Arm)){
	sarm <- sub("\\+", "", arm)
	T11arm[[sarm]] <- rbindlist(T1list, idcol = "outcome")[Arm == arm][, .N, by = ARorder][, p := round(N * 100/sum(N), 2)][]
}

#R!! Outcome-process models



for(ixout in seq_along(outcomes)){
	# ixout <- 5
	iout <- outcomes[ixout]
	for(ixproc in seq_along(processes)){
		# ixproc <- 3
		iproc <- processes[ixproc]
		iopString <- paste(iout, iproc, sep = "-")
		iDout <- matrix(data = NA, nrow = nrow(impactIDs), ncol = 4 + length(ARIMA_parameters))
		rownames(iDout) <- impactIDs[, ParticipantID]
		colnames(iDout) <- c("beta", "SE", ARIMA_parameters, "beta_reg", "SE_reg")
		for(ixid in seq_len(nrow(impactIDs))){
			# ixid <- 56
			iD <- impactIDs[ixid, ParticipantID]
			iDdt <- impactdtres[ParticipantID == iD]
			dt2ts <- ts(iDdt[, iout, with = FALSE])
			dt2xreg <- as.matrix(iDdt[, iproc, with = FALSE])
			iDfit <- try(auto.arima(dt2ts, xreg = dt2xreg), silent = TRUE)
			if(!inherits(iDfit, what = "try-error", which = FALSE)){
				iDout[ixid, "beta"] <- coef(iDfit)[[iproc]]
				iDout[ixid, "SE"] <- sqrt(iDfit$var.coef[iproc, iproc])
				iDout[ixid, ARIMA_parameters] <- iDfit$arma[c(1, 6, 2, 5, 3, 7, 4)]
			}
			iDregression <- try(lm(reformulate(iproc, response = iout), data = iDdt), silent = TRUE)
			if(!inherits(iDregression, what = "try-error", which = FALSE)){
				iDout[ixid, "beta_reg"] <- coef(iDregression)[[iproc]]
				iDout[ixid, "SE_reg"] <- sqrt(vcov(iDregression)[iproc, iproc])
			}
			rm(dt2ts, dt2xreg, iDfit, iDregression)
		}

		iDout <- merge(as.data.table(iDout, keep.rownames = "ParticipantID"), impactIDs, by = "ParticipantID")
		pooledOut[[iopString]] <- iDout


		# meta-analysis for ARIMAX models
		res.nomod <- try(rma(yi = beta, sei = SE, data = iDout, measure = "GEN", method = "REML"), silent = TRUE)
		res.mod <- try(rma(yi = beta, sei = SE, mods = reformulate(termlabels = moderators), data = iDout, measure = "GEN", method = "REML"), silent = TRUE)
		
		if(!inherits(res.mod, what = "try-error", which = FALSE)){
			modOut[[iopString]] <- as.data.table(do.call(cbind, res.mod[c("beta", "se", "pval")]), keep.rownames = "IV")
			mainOut[iopString, c("I2mod", "QEmod", "pQEmod", "QMmod", "pQMmod")] <- unlist(res.mod[c("I2", "QE", "QEp", "QM", "QMp")])
		}
		# meta-analysis for regression models
		# res.reg <- rma(yi = beta_reg, sei = SE_reg, data = iDout, measure = "GEN", method = "REML")

		# meta-analysis for ARIMAX models adjusted per study arm
		res.modArm <- try(rma(yi = beta, sei = SE, data = iDout, mods = ~ Arm, measure = "GEN", method = "REML"), silent = TRUE)

		if(!inherits(res.modArm, what = "try-error", which = FALSE)){
			modArmOut[[iopString]] <- as.data.table(do.call(cbind, res.modArm[c("beta", "se", "pval")]), keep.rownames = "IV")
			mainOut[iopString, c("I2arm", "QEarm", "pQEarm", "QMarm", "pQMarm")] <- unlist(res.modArm[c("I2", "QE", "QEp", "QM", "QMp")])
		}
		
		if(!inherits(res.nomod, what = "try-error", which = FALSE)){
			png(paste0(data_NCpath, "graphics/iARIMAX/forestplots/", iout, "-", iproc, ".png"), bg = "transparent", width = 4800, height = 6400, units = "px", res = 320, type = "cairo")
			print(forest(res.nomod, annotate = FALSE, slab = NA, order = "obs", lwd = 0.5))
			dev.off()
			mainOut[iopString, "beta"] <- res.nomod$beta[1]
			mainOut[iopString, c("SE", "p", "I2", "Q", "pQ")] <- unlist(res.nomod[c("se", "pval", "I2", "QE", "QEp")])
		}



		mainOut[iopString, beta_bands] <- iDout[, betaq := factor(fcase(between(beta, -.1, .1, incbounds = TRUE), 0, 
										between(beta, .1, .2, incbounds = TRUE), 1, 
										between(beta, .2, .3, incbounds = TRUE), 2, 
										beta > .3, 3, 
										between(beta, -.2, -.1, incbounds = TRUE), -1, 
										between(beta, -.3, -.2, incbounds = TRUE), -2, 
										beta < -.3, -3), 
									  levels = seq(-3, 3, 1), labels = beta_bands)
							][levels(betaq), on = "betaq", .N, by = .EACHI
							][, round(100 * N/sum(N), 2)]

		rm(res.nomod, res.mod, res.modArm)


		for(arm in levels(impactdtres$Arm)){
			sarm <- sub("\\+", "", arm)
			
			# meta-analysis for ARIMAX models
			res.nomod <- try(rma(yi = beta, sei = SE, data = iDout, measure = "GEN", method = "REML", subset = Arm == arm), silent = TRUE)
			res.mod <- try(rma(yi = beta, sei = SE, mods = reformulate(termlabels = moderatorsArm[[sarm]]), data = iDout, measure = "GEN", method = "REML", subset = Arm == arm), silent = TRUE)
			if(!inherits(res.mod, what = "try-error", which = FALSE)){
				modOutArm[[sarm]][[iopString]] <- as.data.table(do.call(cbind, res.mod[c("beta", "se", "pval")]), keep.rownames = "IV")
				mainOutArm[[sarm]][iopString, c("I2mod", "QEmod", "pQEmod", "QMmod", "pQMmod")] <- unlist(res.mod[c("I2", "QE", "QEp", "QM", "QMp")])
			}

			if(!inherits(res.nomod, what = "try-error", which = FALSE)){
				png(paste0(data_NCpath, "graphics/iARIMAX/forestplots/", sarm, "/", iout, "-", iproc, ".png"), bg = "transparent", width = 4800, height = 6400, units = "px", res = 320, type = "cairo")
				print(forest(res.nomod, annotate = FALSE, slab = NA, order = "obs", lwd = 0.5))
				dev.off()
				mainOutArm[[sarm]][iopString, "beta"] <- res.nomod$beta[1]
				mainOutArm[[sarm]][iopString, c("SE", "p", "I2", "Q", "pQ")] <- unlist(res.nomod[c("se", "pval", "I2", "QE", "QEp")])
			}


			mainOutArm[[sarm]][iopString, beta_bands] <- iDout[Arm == arm
									   ][levels(betaq), on = "betaq", .N, by = .EACHI
									   ][, round(100 * N/sum(N), 2)]

			rm(res.nomod, res.mod)
		}

		rm(iDout)
		gc()
	}
}



iDout <- rbindlist(pooledOut, idcol = "OPinteraction")[, c("outcome", "process") := tstrsplit(OPinteraction, split = "-", fixed = TRUE)]

modOut <- rbindlist(modOut, idcol = "OPinteraction")[, c("outcome", "process") := tstrsplit(OPinteraction, split = "-", fixed = TRUE)]
modArmOut <- rbindlist(modArmOut, idcol = "OPinteraction")[, c("outcome", "process") := tstrsplit(OPinteraction, split = "-", fixed = TRUE)]

setnames(modOut, old = c("IV", "V1", "se", "pval"), new = c("covariate", "beta", "SE", "p"))
setnames(modArmOut, old = c("IV", "V1", "se", "pval"), new = c("covariate", "beta", "SE", "p"))

T2 <- iDout[, .(reg_avg = mean(abs(beta_reg), na.rm = TRUE),
		iARIMAX_avg = mean(abs(beta), na.rm = TRUE),
		cor_reg_iARIMAX = cor(abs(beta), abs(beta_reg), use = "complete.obs"),
		reg_avg_SE = mean(abs(SE_reg), na.rm = TRUE),
		iARIMAX_avg_SE = mean(abs(SE), na.rm = TRUE)),
	    by = outcome]

T2 <- transpose(T2, make.names = "outcome", keep.names = "Strength of relationship measure")

# forest(res.nomod, slab = paste0(impactIDs[, ParticipantID], " (", impactIDs[, Arm], ")"))
# forest(res.nomod, annotate = FALSE, slab = NA, order = "obs")

for(arm in levels(impactdtres$Arm)){
	sarm <- sub("\\+", "", arm)

	modOutArm[[sarm]] <- rbindlist(modOutArm[[sarm]], idcol = "OPinteraction")[, c("outcome", "process") := tstrsplit(OPinteraction, split = "-", fixed = TRUE)]
	setnames(modOutArm[[sarm]], old = c("IV", "V1", "se", "pval"), new = c("covariate", "beta", "SE", "p"))
	
	T2arm[[sarm]] <- iDout[Arm == arm, .(reg_avg = mean(abs(beta_reg), na.rm = TRUE), 
					     iARIMAX_avg = mean(abs(beta), na.rm = TRUE), 
					     cor_reg_iARIMAX = cor(abs(beta), abs(beta_reg), use = "complete.obs"),
					     reg_avg_SE = mean(abs(SE_reg), na.rm = TRUE),
					     iARIMAX_avg_SE = mean(abs(SE), na.rm = TRUE)), 
			       by = outcome]

	T2arm[[sarm]] <- transpose(T2arm[[sarm]], make.names = "outcome", keep.names = "Strength of relationship measure")
}

# saving matrices
sheet_list <- list("iARIMA order frequencies" = list(round(T1*100, 2), TRUE), 
		   "iARIMA pattern frequencies" = list(T11, FALSE), 
		   "Individual models output" = list(iDout, FALSE), # save matrix of individual models outputs
		   "iARIMAX-regression comparison" = list(T2, FALSE), # save T2 
		   "iARIMAX meta-analysis" = list(mainOut, TRUE), # save T3
		   "iARIMAX meta-an moderators" = list(modOut, FALSE),
		   "iARIMAX meta-an ARM mod" = list(modArmOut, FALSE))
for(sheetname in names(sheet_list)){
	addWorksheet(wbmeta, sheetName = sheetname)
	writeData(wbmeta, sheet = sheetname, sheet_list[[sheetname]][[1]], rowNames = sheet_list[[sheetname]][[2]], colNames = TRUE)
}


for(arm in levels(impactdtres$Arm)){
	sarm <- sub("\\+", "", arm)

	# saving matrices
	sheet_list <- list("iARIMA T1" = list(round(T1arm[[sarm]]*100, 2), TRUE), 
			   "iARIMA patterns" = list(T11arm[[sarm]], FALSE), 
			   "iARIMAX T2" = list(T2arm[[sarm]], FALSE), # save T2 
			   "iARIMAX meta" = list(mainOutArm[[sarm]], TRUE), # save T3
			   "iARIMAX moderators" = list(modOutArm[[sarm]], FALSE))

	names(sheet_list) <- paste(names(sheet_list), sarm, sep = "-")

	for(sheetname in names(sheet_list)){
		addWorksheet(wbmeta, sheetName = sheetname)
		writeData(wbmeta, sheet = sheetname, sheet_list[[sheetname]][[1]], rowNames = sheet_list[[sheetname]][[2]], colNames = TRUE)
	}

}



#R! Save

saveWorkbook(wb, paste0(data_path, "multilevel-models.xlsx"), overwrite = TRUE)
saveWorkbook(wbcwb, paste0(data_path, "correlations-within-between.xlsx"), overwrite = TRUE)
saveWorkbook(wbmeta, paste0(data_NCpath, "iARIMAX.xlsx"), overwrite = TRUE)
save(T1, T1arm, T1list, iDout, modOut, modOutArm, modArmOut, mainOut, mainOutArm, file = paste0(data_NCpath, "iARIMAX.rdata"))

#R! Extra analysis

#R!! Covariables descriptives

wbcov <- loadWorkbook(file = paste0(data_NCpath, "iARIMAX.xlsx"))

addWorksheet(wbcov, sheetName = "Global covariates N")
writeData(wbcov, sheet = "Global covariates N", as.data.frame(unclass(summary(impactIDs[, moderators, with = FALSE]))), rowNames = FALSE, colNames = TRUE)


for(arm in levels(impactdtres$Arm)){
	sarm <- sub("\\+", "", arm)
	addWorksheet(wbcov, sheetName = paste(sarm,"covariates N"))
	writeData(wbcov, sheet = paste(sarm, "covariates N"), as.data.frame(unclass(summary(impactIDs[Arm == arm, moderatorsArm[[sarm]], with = FALSE]))), rowNames = FALSE, colNames = TRUE)

}


saveWorkbook(wbcov, paste0(data_NCpath, "iARIMAX.xlsx"), overwrite = TRUE)


#R!! Id strength plots


load(paste0(data_NCpath, "iARIMAX.rdata"))

for(idx in seq_len(nrow(impactIDs))){
	idsp <- impactIDs[idx, "ParticipantID"]
	for(ixout in seq_along(outcomes)){
		iout <- outcomes[ixout]
		idout <-iDout[ParticipantID == idsp & outcome == iout, .(ParticipantID, outcome, beta, SE, process)]
		# res.id <- try(rma(yi = beta, sei = SE, data = idout, measure = "GEN", method = "REML"), silent = TRUE)
		# if(!inherits(res.id, what = "try-error", which = FALSE)){
		png(paste0(data_NCpath, "graphics/iARIMAX/idplots/", iout, "-id", idsp, ".png"), bg = "transparent", width = 4800, height = 4800, units = "px", res = 320, type = "cairo")
		# print(forest(x = res.id, slab = idout$process, main = paste(iout, idsp)))
		ci.lb <- with(idout, beta - SE * qnorm(0.05/2, lower.tail = FALSE))
		ci.ub <- with(idout, beta + SE * qnorm(0.05/2, lower.tail = FALSE))
		print(idforest <- forest(x = idout$beta, sei = idout$SE, slab = rep("", nrow(idout)), col = ifelse(between(0, ci.lb, ci.ub), "red", "green"), annotate = FALSE, main = paste(iout, idsp)))
		print(text(idforest$textpos[1], seq(nrow(idout), 1, -1), idout$process, pos = 4))
		dev.off()
		# }
		# rm(res.id)
	}
}



#R!! Cluster analysis


clusterdt <- dcast(iDout[, .(ParticipantID, outcome, beta, SE, process, Arm, Responder_Post, Responder_FollowUp)], ParticipantID + outcome + Arm +Responder_Post + Responder_FollowUp ~ process, value.var = c("beta", "SE"))



sink(paste0(data_NCpath, "kmedoidsclustering.txt"), split = TRUE)

for(ixout in seq_along(outcomes)){
	iout <- outcomes[ixout]
	print(iout)
	tmpdt <- clusterdt[outcome == iout, .SD, .SDcols = patterns("^(beta|SE)_")]
	bnkmd <- c(NA, -Inf)
	for(nkmd in 2:9){
		print(paste(nkmd, "clusters"))
		cat("\n")
		pam.res <- pam(tmpdt, nkmd)
		print("Summary of silhouette widths per cluster")
		print(summary(pam.res$silinfo$clus.avg.widths))
		# print("(Weighted) Average silhouette width of total data set")
		# print(pam.res$silinfo$avg.width)
		print("Summary of k-medoids clustering")
		print(summary(pam.res))
		if(bnkmd[2] <= pam.res$silinfo$avg.width){
			bnkmd[1] <- nkmd
			bnkmd[2] <- pam.res$silinfo$avg.width
		}
		cat("\n")
	}
	pam.res <- pam(tmpdt, bnkmd[1])
	# pam.res <- pam(tmpdt, 3)
	png(paste0(data_NCpath, "graphics/iARIMAX/clusters/", iout, ".png"), bg = "transparent", width = 6400, height = 4800, units = "px", res = 320, type = "cairo")
	par(mfcol = c(1,2))
	plot(pam.res, which.plots = 1, main = paste("PAM clustering plot for", iout), col.p = clusterdt$Responder_FollowUp)
	plot(pam.res, which.plots = 2, main = paste("PAM clustering silhouette plot for", iout))
	dev.off()
	cat("\n\n")


}
sink()
