# preview.r


# OS dependencies ----------------------------------------------------------


if(Sys.info()["sysname"] == "Linux"){
	data_path <- "~/Dropbox/juanpablo spanish pain study/"
	dev_lib_path <- "~/R/x86_64-pc-linux-gnu-library/dev"
} else if(Sys.info()["sysname"] == "Windows"){
	data_path <- paste0(grep("^[A-Z]:$", sub(":(.*)", ":",shell("wmic logicaldisk get name", intern = TRUE)), value = TRUE), "/NCext/PSSJD_other/IMPACT")
	data_path <- data_path[file.exists(data_path)]
	data_path <- paste0(data_path, "/")
	dev_lib_path <- .libPaths()
}



# Libraries ---------------------------------------------------------------

library(ggplot2)
library(haven)
library(matrixStats)
library(data.table)
library(gimme)




# Data loading --------------------------------------------------------------

impactdt <- read_sav(paste0(data_path, "Dataset_IMPACT-EMA_v2.sav"))




# Data arranging and merging ----------------------------------------------


impactdt <- as_factor(impactdt)
setDT(impactdt)





# Data check ----------------------------------------------


impactdt
uniqueN(impactdt, by = "ParticipantID")
impactdt[, .N, by = ParticipantID][, table(N)]
112*70 # observations = participants x timepoints
impactdt[, lapply(.SD, \(.x) sum(is.na(.x)))]
names(some_participants_missings) <- some_participants_missings <- c("N_sessions", "Completers", "BPI_Pre", "BPI_Post", "BPI_FollowUp", "NRS_Pre", "NRS_Post", "NRS_FollowUp", "DASS21_Ans_Pre", "DASS21_Ans_Post", "DASS21_Ans_FollowUp", "DASS21_Dep_Pre", "DASS21_Dep_Post", "DASS21_Dep_FollowUp", "DASS21_St_Pre", "PCS_FollowUp", "CPAQ_Pre", "CPAQ_Post", "CPAQ_FollowUp", "BADSSF_Pre", "BADSSF_Post", "BADSSF_FollowUp", "PIPS_Pre", "PIPS_Post", "PIPS_FollowUp")
names(unregular_missings) <- unregular_missings <- c("PainIntensity", "PainControl", "InterferLeasure", "InterferSocial", "InterferWork", "SleepDisturb", "Activity", "Sadness", "Stress", "ExperientialAvoid", "NoContactMoment", "SelfAsContent", "Fusion", "NoContactValues", "Inaction")
lapply(impactdt[, lapply(.SD, \(.x) sum(is.na(.x))), by = ParticipantID][, .SD, .SDcols = some_participants_missings], table)
impactdt[, .N, keyby = .(ParticipantID, Gender, Arm, Completers)][, .N, keyby = .(Gender, Arm, Completers)]
impactdt[, lapply(.SD, \(.x) sum(is.na(.x))), by = .(Gender, Arm, Completers), .SDcols = some_participants_missings]
2800/70 # N_sessions and Completers missing participants
37+37+38 # n for each of the 3 arms
impactdt[, lapply(.SD, \(.x) sum(is.na(.x))), keyby = .(ParticipantID, Gender, Arm, Completers)][, lapply(.SD, \(.x) sum(.x == 70)), keyby = .(Gender, Arm, Completers), .SDcols = some_participants_missings]

lapply(unregular_missings, \(.x) median(impactdt[, lapply(.SD, \(.x) sum(is.na(.x))), by = .(ParticipantID)][[.x]]))

lapply(unregular_missings, \(.y) ggplot(impactdt[, tp := seq_len(.N), by = ParticipantID][, lapply(.SD, \(.x) sum(is.na(.x))), .SDcols = unregular_missings, keyby = tp], aes_string(x = "tp", y = .y)) + geom_bar(stat="identity"))


impactdt[, lapply(.SD, \(.x) uniqueN(.x) - any(is.na(.x))), by = ParticipantID, .SDcols = unregular_missings][, mc := rowMins(as.matrix(.SD)), .SDcols = !c('ParticipantID')][mc <= 1]


# Data saving ----------------------------------------------





setDF(impactdt)
save(impactdt, file = paste0(data_path, "Dataset_IMPACT-EMA_v2.rdata"))

