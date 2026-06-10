
library(tidyverse)
library(lubridate)

# Create folders for the range of litter C:N values, one folder per value
Treatments <- c("1_Overflow", "2_Stoichiometry", "3_Enzyme", "4_Enzyme-Uptake", "5_All-NoUptake", "6_All")

LCN <- paste("CN", 10:90, sep = "")


litter.cn <- function(mech, range){
  
 lcn <- paste(mech, range, sep = "/")
 
 for(i in 1:length(lcn)){
   
   dir.create(lcn[i], recursive = TRUE)
   
 }
  
}

litter.cn(mech = Treatments[1], range = LCN)
litter.cn(mech = Treatments[2], range = LCN)
litter.cn(mech = Treatments[3], range = LCN)
litter.cn(mech = Treatments[4], range = LCN)
litter.cn(mech = Treatments[5], range = LCN)
litter.cn(mech = Treatments[6], range = LCN)



# ------------------------------------------------------------------------------


# Creating parameter and output directories
Combinations <- expand.grid(LCN, Treatments)
Dirs <- paste(Combinations[,2], Combinations[,1], sep = "/")
params.dir <- paste(Dirs, "params", sep = "/")
outputs.dir <- paste(Dirs, "outputs", sep = "/")
substrate.dir <- paste(Dirs, "substrates.txt", sep = "/")

for(i in 1:length(params.dir)){
  
  dir.create(params.dir[i])
  dir.create(outputs.dir[i])
  
}


# Putting common files in all of the folders
files.to.copy <- paste("Common_files", list.files("Common_files"), sep = "/")

for(i in 1:length(Dirs)){
  
  file.copy(from = files.to.copy, to = Dirs[i])
  
}


# Putting the corresponding parameter files in the appropriate folder

# 1. Overflow
param.files1 <- paste(Treatments[1], list.files(Treatments[1], ".txt"), sep = "/")

params.dir1 <- params.dir[str_detect(params.dir, "Overflow")]

for(i in 1:length(params.dir1)){
  
  file.copy(from = param.files1, to = params.dir1[i], overwrite = TRUE)
  
}


# 2. Stoichiometry
param.files2 <- paste(Treatments[2], list.files(Treatments[2], ".txt"), sep = "/")

params.dir2 <- params.dir[str_detect(params.dir, "Stoichiometry")]

for(i in 1:length(params.dir2)){
  
  file.copy(from = param.files2, to = params.dir2[i], overwrite = TRUE)
  
}


# 3. Enzyme
param.files3 <- paste(Treatments[3], list.files(Treatments[3], ".txt"), sep = "/")

params.dir3 <- params.dir[str_detect(params.dir, "Enzyme")]

for(i in 1:length(params.dir3)){
  
  file.copy(from = param.files3, to = params.dir3[i], overwrite = TRUE)
  
}


# 4. Enzyme-Uptake
param.files4 <- paste(Treatments[4], list.files(Treatments[4], ".txt"), sep = "/")

params.dir4 <- params.dir[str_detect(params.dir, "Enzyme-Uptake")]

for(i in 1:length(params.dir4)){
  
  file.copy(from = param.files4, to = params.dir4[i], overwrite = TRUE)
  
}


# 5. All-NoUptake
param.files5 <- paste(Treatments[5], list.files(Treatments[5], ".txt"), sep = "/")

params.dir5 <- params.dir[str_detect(params.dir, "All-NoUptake")]

for(i in 1:length(params.dir5)){
  
  file.copy(from = param.files5, to = params.dir5[i], overwrite = TRUE)
  
}


# 6. All
param.files6 <- paste(Treatments[6], list.files(Treatments[6], ".txt"), sep = "/")

params.dir6 <- params.dir[str_detect(params.dir, "6_All")]

for(i in 1:length(params.dir6)){
  
  file.copy(from = param.files6, to = params.dir6[i], overwrite = TRUE)
  
}



# ------------------------------------------------------------------------------
#                      LITTER CHEMISTRY FILES


# Original DEMENT litter chemistry file
substrates <- read.table("substrates.txt")

# Since phosphatases cleave phosphate groups that can be selectively taken up,
# I decided to eliminate the carbon in the organic phosphorus polymers and put
# it in cellulose, hemicellulose, and starch

substrates[c("Cellulose", "Hemicellulose", "Starch"), "C"] <- 
  substrates[c("Cellulose", "Hemicellulose", "Starch"), "C"] +
  sum(substrates[c("OrgP1", "OrgP2") , "C"])/3

substrates[c("Chitin", "Lignin", "Protein1", "Protein2", "Protein3"), "N"] <- 
  substrates[c("Chitin", "Lignin", "Protein1", "Protein2", "Protein3"), "N"] +
  substrates["OrgP2", "N"]/5

substrates[c("OrgP1", "OrgP2") , "C"] <- 0
substrates["OrgP2" , "N"] <- 0


substrates$CN <- substrates$C/substrates$N
substrates$CP <- substrates$C/substrates$P
substrates$NP <- substrates$N/substrates$P

# First calculate the amount of total nitrogen in litter across the range of 
# litter C:N

# 1. I maintain the original quantity of carbon, 344.5584

# 2. The range of litter C:N goes from 10 to 90
#    N = 100 / Litter C:N

# 3. Calculate the original proportion of substrates

substrates$Tot.C <- colSums(substrates)["C"]
substrates$Tot.N <- colSums(substrates[!rownames(substrates) %in% "Lignin",])["N"] # I exclude lignin from total nitrogen for reasons stated below
substrates <- substrates %>% mutate(C.prop = C/Tot.C, N.prop = N/Tot.N)
substrates["Lignin", "N.prop"] <- 0


Litter <- 
  data.frame(Litter.CN = 10:90) %>% mutate(Litter.N = 344.5584/Litter.CN) %>% 
  mutate(Litter.P = Litter.N/7) %>%
  mutate(Litter.NP = Litter.N/Litter.P)



# I fixed the amount of lignin across all litter C:N values. Consequently, 
# first I must subtract lignin N from total litter N and then calculate the 
# amount of polymeric N
# Then use polymeric N to calculate C of nitrogenous compounds

C.comp.l <- vector(mode = "list", length = dim(Litter)[1])
N.comp.l <- vector(mode = "list", length = dim(Litter)[1])
OrgP1.l <- vector(mode = "list", length = dim(Litter)[1])
CHS.C.total.l <- vector(mode = "list", length = dim(Litter)[1])
Litters.l <- vector(mode = "list", length = dim(Litter)[1])

names(Litters.l) <- paste("CN", Litter$Litter.CN, sep = "")


for(i in 1:dim(Litter)[1]){
  
  # Use the original proportion of litter N to calculate each quantity of polymer N
  # Then use each polymer N to calculate polymer C  
  N.comp.l[[i]] <- (substrates$N.prop*(Litter$Litter.N[i] - substrates["Lignin", "N"]))
  N.comp.l[[i]][7] <- substrates["Lignin", "N"]  # I reintroduce the lignin N
  C.comp.l[[i]] <- N.comp.l[[i]]*substrates$CN
  
  # Calculate orgP1 to satisfy N:P = 7
  OrgP1.l[[i]] <- (Litter$Litter.P[i]/2)*(substrates$P >0)        # OLD Litter$Litter.P[i] - (C.comp.l[[i]]*(1/substrates$CP))[12]
  # OLD C.comp.l[[i]][11] <- OrgP1.l[[i]]*substrates$CP[11]
  
  # Finally, calculate cellulose, hemicellulose, and starch to complete the carbon
  # Use proportion between these three polymers only
  CHS.C.total.l[[i]] <-  344.5584 - sum(C.comp.l[[i]], na.rm = TRUE)
  C.comp.l[[i]][3:5] <- CHS.C.total.l[[i]]*(substrates[3:5, "C"]/sum(substrates[3:5, "C"]))  
  C.comp.l[[i]][c(1:2, 11:12)] <- 0
  
  # Putting everything in matrix format
  Litters.l[[i]] <- cbind(C.comp.l[[i]], N.comp.l[[i]], OrgP1.l[[i]])
  Litters.l[[i]][1:2,] <- 0
  rownames(Litters.l[[i]]) <- rownames(substrates)
  colnames(Litters.l[[i]]) <- c("C", "N", "P")
  
}


# Checking resulting ratios
sapply(Litters.l, function(x) colSums(x)[1]/colSums(x)[2]) # Litter C:N
sapply(Litters.l, function(x) colSums(x)[2]/colSums(x)[3]) # Litter N:P
sapply(Litters.l, function(x) colSums(x)[1])  # Total Litter C


# Putting the litter chemistry tables in the right place

for(i in 1:length(LCN)){
  
  Litter.table <- Litters.l[[LCN[i]]]
  Litter.ID <- LCN[i]
  
  for(i in 1:length(substrate.dir[str_detect(substrate.dir, Litter.ID)])){
  write.table(Litter.table,
              substrate.dir[str_detect(substrate.dir, Litter.ID)][i],
              quote = FALSE,
              sep = "\t")
  }
  
}



# ------------------------------------------------------------------------------

# Generating list of directories for R files
Batch.files <- list.files(recursive = TRUE)[str_detect(list.files(recursive = TRUE), "Batch")]
Batch.files2 <- Batch.files[!str_detect(Batch.files, "Common_files")]
Batch.files3 <- Batch.files2[!str_detect(Batch.files2, "Batch_files")]

# Replace this line with the correct directory in the High Performance Computing environment
Batch.files4 <- paste("/pub/jmuruaro/DEMENT/EGU25/HPC3_5", Batch.files3, sep = "/")

# Exporting txt file with names of R files to be executed
write.table(Batch.files4, file = "Batch_files.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Creating list of directories for the HPC3
Dirs2 <- paste("/pub/jmuruaro/DEMENT/EGU25/HPC3_5", Dirs, sep = "/")

# Exporting txt file with names of directories list of directories
write.table(Dirs2, file = "Directories.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Making a function to print the time stamp automatically
number.of.sims <- length(Batch.files4)

make.job.time <- function(){
  now <- Sys.time()
  Times <- ymd_hms(now) + minutes(seq(from = 0, length.out = number.of.sims))
  
  as_tibble(Times) %>% mutate(Year = substr(value, 3, 4),
                              month = substr(value, 6, 7),
                              day = substr(value, 9, 10),
                              hour = substr(value, 12, 13),
                              minute = substr(value, 15, 16),
                              second = "00") %>%
    unite("Timestamp", Year:second, sep = "") %>%
    select(Timestamp) %>% unlist() %>% unname()
  
}


# Generating list of timestamps, the nubers must be in format YYMMDDHHMMSS. Set the seconds to 00, so that the task ID can be added at the end of the simulation
Sim.IDs <- make.job.time()

# Exporting txt file with timestamps for each simulation to be executed
write.table(Sim.IDs, file = "Timestamps.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)








