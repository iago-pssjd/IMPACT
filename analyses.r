

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

# Libraries ---------------------------------------------------------------

library(openxlsx) # createWorkbook
# library(optimx)
library(lattice) # dotplot
# library(ggplot2)
library(matrixStats) # rowMins
library(haven) # read_sav
# library(nlme)
library(lme4) # lmer
library(data.table)
# library(gimme)


# Auxiliar functions ------------------------------------------------------

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


# Data loading --------------------------------------------------------------

impactdt <- read_sav(paste0(data_path, "Dataset_IMPACT-EMA_v2.sav"))
impactdt <- as_factor(impactdt)
setDT(impactdt)
wb <- createWorkbook()


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
		# rimodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = ~ 1 | ParticipantID, method = "ML", na.action = na.omit)
		# rsmodel <- lme(reformulate(iproc, response = iout), data = impactdtres, random = as.formula(paste("~", iproc, "| ParticipantID")), method = "ML", na.action = na.omit)
		# chisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["L.Ratio"]][2]
		# pchisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["p-value"]][2]
		# beta10[ixproc, ixout] <- quantile(coef(rsmodel)[[iproc]], probs = c(0.1))
		# beta90[ixproc, ixout] <- quantile(coef(rsmodel)[[iproc]], probs = c(0.9))
		# plot(coef(rsmodel), panel = function(...) {panel.dotplot(..., col = 'darkgreen'); panel.abline(v = fixef(rsmodel), lty = 3, col = 'purple')})

		rimodel <- lmer(paste0(iout, " ~ ", iproc, " + (1|ParticipantID)"), data = impactdtres, REML = FALSE, lmerControl(optimizer = "Nelder_Mead"))
		rsmodel <- lmer(paste0(iout, " ~ ", iproc, " + (", iproc,"|ParticipantID)"), data = impactdtres, REML = FALSE, lmerControl(optimizer = "Nelder_Mead"))
		# ranef(rsmodel)[["ParticipantID"]][,"(Intercept)"] # merMod methods
		# ranef(rsmodel)[["ParticipantID"]][,iproc] # merMod methods
		chisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Chisq"]][2]
		pchisq[ixproc, ixout] <- anova(rsmodel, rimodel)[["Pr(>Chisq)"]][2]
		beta[ixproc, ixout] <- fixef(rsmodel)[[iproc]]
		beta10[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.1))
		beta90[ixproc, ixout] <- quantile(coef(rsmodel)[["ParticipantID"]][[iproc]], probs = c(0.9))
		# plot(compareFits(coef(rsmodel), coef(rimodel)), mark = fixef(rsmodel))
		if(Sys.info()["sysname"] == "Windows"){
		  png(paste0(data_path, "graphics/", iout, "-", iproc, ".png"), bg = "transparent")
		  print(dotplot(rsmodel, col = 'darkgreen', col.line.v = 'purple', lty.v = 3, scales = list(x = list(relation = 'free'), y = list(draw = FALSE)), par.settings = list(strip.border = list(alpha = 0)), strip = strip.custom(factor.levels = c("Intercept (Z score)", "Slope (beta)"), bg = "transparent"))[["ParticipantID"]])
		  dev.off()
		}
	}
}
sheet_list <- list("chi-square" = chisq, "chi-square p-value" = pchisq, "beta" = beta, "beta10" = beta10, "beta90" = beta90)
for(sheetname in names(sheet_list)){
  addWorksheet(wb, sheetName = sheetname)
  writeData(wb, sheet = sheetname, sheet_list[[sheetname]], rowNames = TRUE)
}

saveWorkbook(wb, paste0(data_path, "multilevel-models.xlsx"), overwrite = TRUE)

