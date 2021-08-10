.onLoad <- function(...){
	packageStartupMessage("\n")
	packageStartupMessage("Welcome to RNMF.")
	packageStartupMessage("\n")
	packageStartupMessage("Version: ",utils::packageDescription('RNMF')$Version)
	packageStartupMessage("\n")
	packageStartupMessage("If this is your first time running RNMF you should see ?RNMF")
	packageStartupMessage("\n")
	}