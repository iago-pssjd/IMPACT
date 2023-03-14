

if(Sys.info()["sysname"] == "Linux"){
	data_path <- "~/Dropbox/juanpablo spanish pain study/"
	dev_lib_path <- "~/R/x86_64-pc-linux-gnu-library/dev"
} else if(Sys.info()["sysname"] == "Windows"){
	data_path <- paste0(grep("^[A-Z]:$", sub(":(.*)", ":",shell("wmic logicaldisk get name", intern = TRUE)), value = TRUE), "/PSSJD/IMPACT/data")
	data_path <- data_path[file.exists(data_path)]
	data_path <- paste0(data_path, "/")
	dev_lib_path <- .libPaths()
}



# Libraries ---------------------------------------------------------------

library(haven)
library(data.table)
library(gimme)



# Data loading --------------------------------------------------------------

impactdt <- read_sav(paste0(data_path, "Dataset_IMPACT-EMA_v2.sav"))

