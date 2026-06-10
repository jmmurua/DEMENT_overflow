# To run via bash, include line 3 and exclude 5

args <- commandArgs(trailingOnly=T)

#args <- c("time slot", "id") #YYMMDDHHSS, integer

job.time <- args[1]
task.ID <- as.numeric(args[2])

source("DEMENT.0.7.6.R")

dir.create("outputs",showWarnings=F)

out<- TraitModel(job.time,task.ID)
filename <- paste("outputs/", out[[1]]$timestamp, ".RData", sep = "")
save.image(file=filename)
