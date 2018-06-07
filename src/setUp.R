# This module loads all required packages.
# If a package is installed, it will only be loaded.
# Otherwise, it will be installed and then loaded.
# Packages that are from Bioconductor will be installed using "biocLite".
# Packages from CRAN will be installed with "install.packages".

# Load biocLite from bioconductor.org
# try http:// if https:// URLs are not supported
tryCatch(source("https://bioconductor.org/biocLite.R"),
	error=function(e){
		message("https:// URLs are not supported so http:// URL is used.")
		source("http://bioconductor.org/biocLite.R")
	},
	silent=TRUE
	)

# Install packages
if(!require("limma")){
	biocLite("limma")
	load("limma")
}

if(!require("gage")){
	biocLite("gage")
	load("gage")
}

if(!require("globaltest")){
	biocLite("globaltest")
	load("globaltest")
}

if(!require("GSA")){
	install.packages("GSA")
	load("GSA")
}

if(!require("GSVA")){
	biocLite("GSVA")
	load("GSVA")
}

if(!require("PADOG")){
	biocLite("PADOG")
	load("PADOG")
}

if(!require("safe")){
	biocLite("safe")
	load("safe")
}
