suppressPackageStartupMessages(library(biomAid))
s <- simTrialData(nvar = 12, nsite = 4, nrep = 1, seed = 1, verbose = FALSE,
                   sim.options = list(variety_prefix = "Line",
                                       site_prefix    = "Env"))
cat("Variety levels:\n"); print(levels(factor(s$data$Variety)))
cat("\nSite levels:\n");    print(levels(factor(s$data$Site)))
cat("\nColumn names:\n");   print(names(s$data))
cat("\nHead:\n");           print(head(s$data))
cat("\nClass of Variety and Site columns:\n")
cat("  Variety:", class(s$data$Variety), "\n")
cat("  Site   :", class(s$data$Site), "\n")
