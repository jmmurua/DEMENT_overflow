
#     Processing the data from simulations with high redundancy and redirected overflow to DOC
#     ------------------------------------



library(tidyverse)


# Loading the data
all.Data <- paste(list.files(recursive = TRUE, pattern = ".RData"),
                  sep = "/")

# Including only simulations with high redundancy and redirected overflow
all.Data <- all.Data[str_detect(all.Data, "03.R")]



# Loading substrate files

subs.dirs <- # Directories for substrate files
  paste(list.files(recursive = TRUE, pattern = "substrates"), sep = "/")[
    str_detect(paste(list.files(recursive = TRUE, pattern = "substrates"), sep = "/"), "Overflow")
  ]


all.subs <- vector(length = length(subs.dirs), mode = "list")

for(i in 1:length(subs.dirs)){
  
  all.subs[[i]] <- read.table(subs.dirs[i])
  
}

# Naming the substrate files according to C:N
names(all.subs) <- str_extract(subs.dirs, "CN([0-9]{2})")

all.subs.rat <- lapply(all.subs, function(x) mutate(x, NC = replace_na(N/C, 0), 
                                                    PC = replace_na(P/C, 0)))



# Function to load multiple RData files without overwriting identically named objects
loadEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}


# Creating the names of the environments where the data will be stored
envs <- 
  str_c(
    str_extract(str_split(all.Data, "/", simplify = TRUE)[,1], "[A-Za-z]+([-]*)([A-za-z]*)"), "_",
    str_extract(str_split(all.Data, "/", simplify = TRUE)[,2], "CN[0-9]{2}"),
    rep(".env", length(all.Data)))
Envs <- str_c(envs, unlist(sapply(table(envs), FUN = function(x){1:x[1]}, simplify = TRUE)))



# Saving .RData files in separate environments
for(i in 1:length(all.Data)){
  
  assign(Envs[i], loadEnvironment(all.Data[i]))
}


# Naming the simulations
Sims <- str_remove(Envs, ".env[0-9]")
XSims <- rep(NA, length(Sims))

for(i in 1:length(XSims)){
  XSims[i] <- str_glue("X", Sims[i])
}


# Choose pulse
p <- 1

# Loading biomass series for each simulation
Mic_Sum <- vector(length = length(Envs), mode = "list")

for (i in 1:length(Envs)){
  Mic_Sum[[i]] <- get(Envs[i])$out[[p]]$Mic_Sum
  
}

names(Mic_Sum) <- Sims

Biomass <- lapply(Mic_Sum, function(x) mutate(as_tibble(x), Biomass = rowSums(x)))
Biomass <- lapply(Biomass, function(x) rbind(x, matrix(ncol = 4, nrow = 1000 - nrow(x), dimnames = list(NULL, c("C", "N", "P", "Biomass"))))) # Filling simulations with less than 1000 days
Biomass <- lapply(Biomass, function(x) cbind(x, days = 1:nrow(x)))

# Adding simulation ID as a column
for (i in 1:length(Biomass)){
  Biomass[[i]] <- cbind(Biomass[[i]], Simulation = rep(names(Biomass)[i], nrow(Biomass[[i]])))
  
}


# Loading biomass growth rate per element
Mic.growth <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)){
  
  Mic.growth[[i]] <- get(Envs[i])$out[[p]]$Mic.growthSeries
  
}

names(Mic.growth) <- Sims

Growth <- lapply(Mic.growth, function(x) rbind(x, matrix(ncol = 3, nrow = 1000 - nrow(x))))
Growth <- lapply(Growth, function(x) cbind(x, days = 1:nrow(x)))

# Adding simulation ID as a column
for (i in 1:length(Growth)){
  Growth[[i]] <- cbind(Growth[[i]], Simulation = rep(names(Growth)[i], nrow(Growth[[i]])))
  
}

# Loading taxa population dynamics for each simulation
MicrobesSeries <- vector(length = length(Envs), mode = 'list')

for(i in 1:length(Envs)){
  
  MicrobesSeries[[i]]<- get(Envs[i])$out[[p]]$MicrobesSeries
}

names(MicrobesSeries) <- Sims

MicrobesSeries <- lapply(MicrobesSeries, function(x) rbind(x, matrix(ncol = ncol(x), nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
MicrobesSeries <- lapply(MicrobesSeries, function(x){cbind(x, days = 1:nrow(x))})
MicrobesSeries <- lapply(MicrobesSeries, as.data.frame)

# Adding simulation ID as a column
for (i in 1:length(MicrobesSeries)){
  MicrobesSeries[[i]] <- cbind(MicrobesSeries[[i]], Simulation = rep(names(MicrobesSeries)[i], nrow(MicrobesSeries[[i]])))
  
}



# Loading enzyme secretion series for each simulation
EnzSer <- vector(length = length(Envs), mode = "list")

for (i in 1:length(Envs)){
  EnzSer[[i]] <- get(Envs[i])$out[[p]]$EnzymesSeries
  rownames(EnzSer[[i]]) <- NULL
}

names(EnzSer) <- Sims

EnzSer <- lapply(EnzSer, function(x) rbind(x, matrix(ncol = 12, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
EnzSer <- lapply(EnzSer, function(x){cbind(x, days = 1:nrow(x))})
EnzSer <- lapply(EnzSer, as.data.frame)

# Adding simulation ID as a column
for (i in 1:length(EnzSer)){
  EnzSer[[i]] <- cbind(EnzSer[[i]], Simulation = rep(names(EnzSer)[i], nrow(EnzSer[[i]])))
  
}


# Loading enzyme production series for each simulation
EnzProdSer <- vector(length = length(Envs), mode = "list")

for (i in 1:length(Envs)){
  EnzProdSer[[i]] <- get(Envs[i])$out[[p]]$EnzProdSeries
  rownames(EnzProdSer[[i]]) <- NULL
}

names(EnzProdSer) <- Sims

EnzProdSer <- lapply(EnzProdSer, function(x) rbind(x, matrix(ncol = 12, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
EnzProdSer <- lapply(EnzProdSer, function(x){cbind(x, days = 1:nrow(x))})
EnzProdSer <- lapply(EnzProdSer, as.data.frame)

# Adding simulation ID as a column
for (i in 1:length(EnzProdSer)){
  EnzProdSer[[i]] <- cbind(EnzProdSer[[i]], Simulation = rep(names(EnzProdSer)[i], nrow(EnzProdSer[[i]])))
  
}



# Loading NH3 and PO4 series
NH4Ser <- vector(length = length(Envs), mode = "list")
PO4Ser <- vector(length = length(Envs), mode = "list")

for (i in 1:length(Envs)){
  NH4Ser[[i]] <- get(Envs[i])$out[[p]]$NH4Series
  #colnames(NH4Ser[[i]]) <- "NH4"
  
  PO4Ser[[i]] <- get(Envs[i])$out[[p]]$PO4Series
  #colnames(PO4Ser[[i]]) <- "PO4"
  
}

names(NH4Ser) <- Sims
names(PO4Ser) <- Sims

NH4Ser <- lapply(NH4Ser, function(x) rbind(x, matrix(ncol = 1, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
PO4Ser <- lapply(PO4Ser, function(x) rbind(x, matrix(ncol = 1, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days


# Loading mineralization series
N.minSer <- vector(length = length(Envs), mode = "list")
P.minSer <- vector(length = length(Envs), mode = "list")

for (i in 1:length(Envs)){
  N.minSer[[i]] <- get(Envs[i])$out[[p]]$N.mineralSeries
  P.minSer[[i]] <- get(Envs[i])$out[[p]]$P.mineralSeries
  
}

names(N.minSer) <- Sims
names(P.minSer) <- Sims

N.minSer <- lapply(N.minSer, function(x) c(x, rep(NA, 1000 - length(x)))) # Filling simulations with less than 1000 days
P.minSer <- lapply(P.minSer, function(x) c(x, rep(NA, 1000 - length(x)))) # Filling simulations with less than 1000 days



# Loading respiration series for each simulation
RespSeries <- vector(length = length(Envs), mode = "list")
Resp.compSeries <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)){
  RespSeries[[i]] <-  get(Envs[i])$out[[p]]$RespSeries
  Resp.compSeries[[i]] <- get(Envs[i])$out[[p]]$Resp.compSeries
}

names(RespSeries) <- Sims
names(Resp.compSeries) <- Sims

RespSeries <- lapply(RespSeries, function(x) c(x, rep(NA, 1000 - length(x)))) # Filling simulations with less than 1000 days
Resp.compSeries <- lapply(Resp.compSeries, function(x) rbind(x, matrix(ncol = 3, nrow = 1000 - nrow(x))))
Resp.compSeries <- lapply(Resp.compSeries, function(x){cbind(x, days = 1:nrow(x))})

# Adding simulation ID as a column
for (i in 1:length(Resp.compSeries)){
  Resp.compSeries[[i]] <- cbind(Resp.compSeries[[i]], Simulation = rep(names(Resp.compSeries)[i], nrow(Resp.compSeries[[i]])))
  
}



# Loading cummulative substrate series for each simulation
Cum_SubstrateSeries <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)){
  
  Cum_SubstrateSeries[[i]] <- get(Envs[[i]])$out[[p]]$Cum_SubstrateSeries
}

names(Cum_SubstrateSeries) <- Sims

Cum_SubstrateSeries <- lapply(Cum_SubstrateSeries, function(x) rbind(x, matrix(ncol = 3, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days


# Loading substrate series for each simulation
SubstratesSeries <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)) {
  
  SubstratesSeries[[i]] <- get(Envs[[i]])$out[[p]]$SubstratesSeries
  
}

names(SubstratesSeries) <- Sims

SubstratesSeries <- lapply(SubstratesSeries, function(x) rbind(x, matrix(ncol = 12, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
SubstratesSeries <- lapply(SubstratesSeries, function(x){cbind(x, days = 1:nrow(x))})
SubstratesSeries <- lapply(SubstratesSeries, as.data.frame)


# Adding simulation ID as a column
for (i in 1:length(SubstratesSeries)){
  SubstratesSeries[[i]] <- cbind(SubstratesSeries[[i]], Simulation = rep(names(SubstratesSeries)[i], nrow(SubstratesSeries[[i]])))
  
}


# Litter mass time-series for each simulation
Litter <- lapply(Cum_SubstrateSeries, rowSums)


# Litter C:N time-series for each simulation
L.CN <- lapply(Cum_SubstrateSeries, FUN = function(x){x[, "C"]/x[, "N"]})


# Litter C:P time-series for each simulation
L.CP <- lapply(Cum_SubstrateSeries, FUN = function(x){x[, "C"]/x[, "P"]})


# Litter C:P time-series for each simulation
L.NP <- lapply(Cum_SubstrateSeries, FUN = function(x){x[, "N"]/x[, "P"]})



# Loading enzyme identity
# ReqEnz <- get(Envs[1])$out[[1]]$ReqEnz[[1]]
# rownames(ReqEnz) <- rownames(substrates)

ReqEnz <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)){
  
  ReqEnz[[i]] <- get(Envs[[i]])$out[[p]]$ReqEnz[[1]]
  rownames(ReqEnz[[i]]) <- rownames(all.subs[[1]])
  
}

names(ReqEnz) <- Sims




# Series of litter associated substrates only
#Litter.subs <- lapply(SubstratesSeries, function(x) x[, c(rownames(substrates)[-1:-2], "Simulation", "days")])

Litter.subs <- vector(length = length(Envs), mode = "list")

for(i in 1:length(SubstratesSeries)){
  
  Litter.subs[[i]] <- SubstratesSeries[[i]][, c(rownames(all.subs[[str_extract(names(SubstratesSeries)[i], "CN[0-9]+")]])[-1:-2], "Simulation", "days")]
  
}



# Loading monomer series
C.MonomersSeries <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)) {
  
  C.MonomersSeries[[i]] <- get(Envs[[i]])$out[[p]]$C.MonomersSeries
  
}

names(C.MonomersSeries) <- Sims

C.MonomersSeries <- lapply(C.MonomersSeries, function(x) rbind(x, matrix(ncol = 14, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
C.MonomersSeries <- lapply(C.MonomersSeries, function(x){cbind(x, days = 1:nrow(x))})
C.MonomersSeries <- lapply(C.MonomersSeries, as.data.frame)

# Adding simulation ID as a column
for (i in 1:length(C.MonomersSeries)){
  C.MonomersSeries[[i]] <- cbind(C.MonomersSeries[[i]], Simulation = rep(names(C.MonomersSeries)[i], nrow(C.MonomersSeries[[i]])))
  
}


N.MonomersSeries <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)) {
  
  N.MonomersSeries[[i]] <- get(Envs[[i]])$out[[p]]$N.MonomersSeries
  
}

names(N.MonomersSeries) <- Sims

N.MonomersSeries <- lapply(N.MonomersSeries, function(x) rbind(x, matrix(ncol = 14, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
N.MonomersSeries <- lapply(N.MonomersSeries, function(x){cbind(x, days = 1:nrow(x))})
N.MonomersSeries <- lapply(N.MonomersSeries, as.data.frame)

# Adding simulation ID as a column
for (i in 1:length(N.MonomersSeries)){
  N.MonomersSeries[[i]] <- cbind(N.MonomersSeries[[i]], Simulation = rep(names(N.MonomersSeries)[i], nrow(N.MonomersSeries[[i]])))
  
}


# Loading monomer production series
MonProdSer <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)) {
  
  MonProdSer[[i]] <- get(Envs[[i]])$out[[p]]$MonProdSeries
  
}

names(MonProdSer) <- Sims

MonProdSer <- lapply(MonProdSer, function(x) rbind(x, matrix(ncol = 12, nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
MonProdSer <- lapply(MonProdSer, function(x){cbind(x, days = 1:nrow(x))})
MonProdSer <- lapply(MonProdSer, as.data.frame)

# Adding simulation ID as a column
for (i in 1:length(MonProdSer)){
  MonProdSer[[i]] <- cbind(MonProdSer[[i]], Simulation = rep(names(MonProdSer)[i], nrow(MonProdSer[[i]])))
  
}


# Loading carbon, nitrogen, and phosphorus uptake series
C.UptSer <- vector(length = length(Envs), mode = "list")
N.UptSer <- vector(length = length(Envs), mode = "list")
P.UptSer <- vector(length = length(Envs), mode = "list")
N.MonomerUptakeSeries <- vector(length = length(Envs), mode = "list")

for (i in 1:length(Envs)){
  
  C.UptSer[[i]] <- get(Envs[i])$out[[p]]$C.UptakeSeries
  rownames(C.UptSer[[i]]) <- NULL
  
  N.UptSer[[i]] <- get(Envs[i])$out[[p]]$N.UptakeSeries
  rownames(N.UptSer[[i]]) <- NULL
  
  P.UptSer[[i]] <- get(Envs[i])$out[[p]]$P.UptakeSeries
  rownames(P.UptSer[[i]]) <- NULL
  
  N.MonomerUptakeSeries[[i]] <- get(Envs[i])$out[[p]]$N.MonomerUptakeSeries
  
}

names(C.UptSer) <- Sims
names(N.UptSer) <- Sims
names(P.UptSer) <- Sims
names(N.MonomerUptakeSeries) <- Sims

C.UptSer <- lapply(C.UptSer, function(x) rbind(x, matrix(ncol = ncol(x), nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
N.UptSer <- lapply(N.UptSer, function(x) rbind(x, matrix(ncol = ncol(x), nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
P.UptSer <- lapply(P.UptSer, function(x) rbind(x, matrix(ncol = ncol(x), nrow = 1000 - nrow(x)))) # Filling simulations with less than 1000 days
N.MonomerUptakeSeries <- lapply(N.MonomerUptakeSeries, function(x) rbind(x, matrix(ncol = 14, nrow = 1000 - nrow(x))))

C.UptSer <- lapply(C.UptSer, function(x) as.data.frame(cbind(rowSums(x), days = 1:nrow(x))))
N.UptSer <- lapply(N.UptSer, function(x) as.data.frame(cbind(rowSums(x), days = 1:nrow(x))))
P.UptSer <- lapply(P.UptSer, function(x) as.data.frame(cbind(rowSums(x), days = 1:nrow(x))))
N.MonomerUptakeSeries <- lapply(N.MonomerUptakeSeries, function(x) as.data.frame(cbind(x, days = 1:nrow(x))))


# N.UptSer <- lapply(N.UptSer, cbind, days = 1:nrow(N.UptSer[[1]]))
# N.UptSer <- lapply(N.UptSer, as.data.frame)
# 
# P.UptSer <- lapply(P.UptSer, cbind, days = 1:nrow(P.UptSer[[1]]))
# P.UptSer <- lapply(P.UptSer, as.data.frame)

# Adding simulation ID as a column
for (i in 1:length(C.UptSer)){
  C.UptSer[[i]] <- cbind(C.UptSer[[i]], Simulation = rep(names(C.UptSer)[i], nrow(C.UptSer[[i]])))
}

for (i in 1:length(N.UptSer)){
  N.UptSer[[i]] <- cbind(N.UptSer[[i]], Simulation = rep(names(N.UptSer)[i], nrow(N.UptSer[[i]])))
}

for (i in 1:length(P.UptSer)){
  P.UptSer[[i]] <- cbind(P.UptSer[[i]], Simulation = rep(names(P.UptSer)[i], nrow(P.UptSer[[i]])))
}

for (i in 1:length(N.MonomerUptakeSeries)){
  N.MonomerUptakeSeries[[i]] <- cbind(N.MonomerUptakeSeries[[i]], Simulation = rep(names(N.MonomerUptakeSeries)[i], nrow(N.MonomerUptakeSeries[[i]])))
}


# Loading net communitary carbon use efficiency for each simulation
Net.CUE_Series <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)){
  
  Net.CUE_Series[[i]] <- get(Envs[[i]])$out[[p]]$Net.CUE_Series
}

names(Net.CUE_Series) <- Sims

Net.CUE_Series <- lapply(Net.CUE_Series, function(x) c(x, rep(NA, 1000 - length(x)))) # Filling simulations with less than 1000 days

for(i in 1:length(Envs)){   # Setting NAs to 0
  
  Net.CUE_Series[[i]][lapply(Net.CUE_Series, is.na)[[i]]] <- 0
  
}




# Identifying the enzymes that degrade each substrate
Cellulases <- vector(length = length(Envs), mode = "list")
Hemicellulases <- vector(length = length(Envs), mode = "list")
Amilases <- vector(length = length(Envs), mode = "list")
Chitinases <- vector(length = length(Envs), mode = "list")
Ligninases <- vector(length = length(Envs), mode = "list")


for(i in 1:length(ReqEnz)){
  
  Cellulases[[i]] <- colnames(ReqEnz[[i]]["Cellulose", ReqEnz[[i]]["Cellulose",] == 1, drop = FALSE])
  Hemicellulases[[i]] <- colnames(ReqEnz[[i]]["Hemicellulose", ReqEnz[[i]]["Hemicellulose",] == 1, drop = FALSE])
  Amilases[[i]] <- colnames(ReqEnz[[i]]["Starch", ReqEnz[[i]]["Starch",] == 1, drop = FALSE])
  Chitinases[[i]] <- colnames(ReqEnz[[i]]["Chitin", ReqEnz[[i]]["Chitin",] == 1, drop = FALSE])
  Ligninases[[i]] <- colnames(ReqEnz[[i]]["Lignin", ReqEnz[[i]]["Lignin",] == 1, drop = FALSE])
  
}



# Making for loop for the proteases
logic.v <- vector(length = dim(ReqEnz[[1]])[2])

# logic.v <- 
# lapply(ReqEnz, function(x) dim(x)[2]) %>% 
#   lapply(function(x) vector(length = x)) #%>%
#   #lapply(function(x) apply(x, 1, any(x[c("Protein1", "Protein2", "Protein3"),] == 1)))
# 
# lapply(ReqEnz, function(x) apply(x, 2, any(x[c("Protein1", "Protein2", "Protein3"),] == 1)))

for(i in 1:dim(ReqEnz[[1]])[2]){
  logic.v[i] <- any(ReqEnz[[1]][c("Protein1", "Protein2", "Protein3"), i] == 1)
}

Proteases <- colnames(ReqEnz[[1]])[logic.v]


# Making for loop for the phosphatases1
logic.v2 <- vector(length = dim(ReqEnz[[1]])[2])

for(i in 1:dim(ReqEnz[[1]])[2]){
  logic.v2[i] <- any(ReqEnz[[1]][c("OrgP1"), i] == 1)
}

Phosphatases1 <- colnames(ReqEnz[[1]])[logic.v2]


# Making for loop for the phosphatases2
logic.v3 <- vector(length = dim(ReqEnz[[1]])[2])

for(i in 1:dim(ReqEnz[[1]])[2]){
  logic.v3[i] <- any(ReqEnz[[1]][c("OrgP2"), i] == 1)
}

Phosphatases2 <- colnames(ReqEnz[[1]])[logic.v3]

# Making for loop for necro-enzymes
logic.v4 <- vector(length = dim(ReqEnz[[1]])[2])

for(i in 1:dim(ReqEnz[[1]])[2]){
  logic.v4[i] <- any(ReqEnz[[1]][c("DeadMic", "DeadEnz"), i] == 1)
}

Necro.enzymes <- colnames(ReqEnz[[1]])[logic.v4]

# Pooling the enzymes identities into elements
C.enzymes <- unique(c(Cellulases[[1]], Hemicellulases[[1]], Amilases[[1]]))
N.enzymes <- unique(c(Chitinases[[1]], Proteases, Ligninases[[1]]))
P.enzymes <- unique(c(Phosphatases1, Phosphatases2))
Recyc.enzymes <- unique(c(Necro.enzymes))


# Making look up table
Look.Enz <- 
  c(
    setNames(rep("C", length(C.enzymes)), C.enzymes),
    setNames(rep("N", length(N.enzymes)), N.enzymes),
    setNames(rep("P", length(P.enzymes)), P.enzymes),
    setNames(rep("Recyc", length(Recyc.enzymes)), Recyc.enzymes)
  )







# Loading enzymes per taxa
EnzGenes <- vector(length = length(Envs), mode = "list")

for(i in 1:length(Envs)) {
  
  EnzGenes[[i]] <- get(Envs[[i]])$out[[p]]$EnzGenes
  
}

names(EnzGenes) <- Sims

EnzGenes <- lapply(EnzGenes, as.data.frame)


# Adding simulation ID as a column
for (i in 1:length(EnzGenes)){
  EnzGenes[[i]] <- cbind(EnzGenes[[i]], Simulation = rep(names(EnzGenes)[i], nrow(EnzGenes[[i]])))
  
}

EnzGenes <- lapply(EnzGenes, rownames_to_column, "Taxa")
EnzGenes <- lapply(EnzGenes, function(x) pivot_longer(x, cols = starts_with("Enz"),
                                                      names_to = "Enzyme",
                                                      values_to = "EnzGene"))

# Adding column for enzyme type
EnzGenes <- lapply(EnzGenes, function(x) mutate(x, Enz.type = unname(Look.Enz[x$Enzyme])))


# Characterizing taxa by the enzymes they possess
EnzGenes.df <- do.call("rbind", EnzGenes)
EnzGenes.df.w <- 
  EnzGenes.df %>% #filter(EnzGene == 1) %>% 
  group_by(Simulation, Taxa, Enz.type) %>% 
  summarize(Sum.enz = sum(EnzGene)) %>%
  pivot_wider(names_from = Enz.type, values_from = Sum.enz) %>%
  mutate(C = replace_na(C, 0), 
         N = replace_na(N, 0),
         P = replace_na(P, 0),
         Recyc = replace_na(Recyc, 0),
         Tot.Enz = sum(C, N, P, Recyc, na.rm = TRUE),
         Elem.cover = sum(C > 0, N > 0, P > 0),
         Recycler = Recyc > 0) %>% ungroup()












# ------------------------------------------------------------------------------
#                            Putting all together


Biomass.df <- do.call("rbind", Biomass) %>%
  separate(Simulation, into = c("Mechanism", "Litter.CN"), 
           sep = "_", remove = FALSE)
rownames(Biomass.df) <- NULL
names(Biomass.df)[1:3] <- c("C.biomass", "N.biomass", "P.biomass")
write_tsv(Biomass.df, "Biomass_HRo.txt")

if(ncol(MicrobesSeries[[1]]) > 3){
  Taxa.df <- do.call("rbind", MicrobesSeries)
  rownames(Taxa.df) <- NULL
  Taxa.df.l <- Taxa.df %>% pivot_longer(cols = starts_with("Tax"), names_to = "Taxa", values_to = "Abundance") %>%
    separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
  
  TaxaGenes <- inner_join(EnzGenes.df.w, Taxa.df.l)
  write_tsv(Taxa.df, "Taxa_HRo.txt"); write_tsv(Taxa.df.l, "TaxaL_HRo.txt"); write_tsv(TaxaGenes, "TaxaGenes_HRo.txt")
  
}



Growth.df <- do.call("rbind", Growth) %>% as_tibble() %>% 
  mutate(across(c("C", "N", "P", "days"), as.numeric)) %>%
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE) %>%
  pivot_longer(cols = c("C", "N", "P"), names_to = "Element", values_to = "Growth")
write_tsv(Growth.df, "Growth_HRo.txt")


Respiration.df <- as.data.frame(RespSeries) %>%  cbind(days = 1:length(RespSeries[[1]])) %>%
  pivot_longer(cols = 1:length(RespSeries), names_to = "Simulation", 
               values_to = "Respiration") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(Respiration.df, "Respiration_HRo.txt")


Resp.comp.df <- do.call("rbind", Resp.compSeries) %>% as_tibble() %>% 
  mutate(across(c("Maint", "Growth", "Overflow", "days"), as.numeric)) %>%
  separate(Simulation, into = c("Mechanism", "Litter.CN"), 
           sep = "_", remove = FALSE) %>%
  pivot_longer(cols = c("Maint", "Growth", "Overflow"), names_to = "Source", values_to = "Resp")
write_tsv(Resp.comp.df, "RespComp_HRo.txt")


Necromass.df <- do.call("rbind", SubstratesSeries) %>% as_tibble() %>%
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE) %>%
  select(Simulation, Mechanism, Litter.CN, days, DeadMic, DeadEnz)
write_tsv(Necromass.df, "Necromass_HRo.txt")


Litter.df <- as.data.frame(Litter) %>%  cbind(days = 1:length(Litter[[1]])) %>%
  pivot_longer(cols = 1:length(Litter), names_to = "Simulation", 
               values_to = "Litter") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(Litter.df, "Litter_HRo.txt")

Litter.subs.df <- do.call("rbind", Litter.subs) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(Litter.subs.df) <- NULL
Litter.subs.df$Litter.C <- rowSums(Litter.subs.df[, rownames(all.subs[[1]])[c(-1:-2, -11:-12)]])
Litter.subs.df$Litter.N <- rowSums(Litter.subs.df[, rownames(all.subs[[1]])[c(6:10)]])
write_tsv(Litter.subs.df, "LitterSubs_HRo.txt")

# Making matrix of substrate C:N values
subsNC.m <- matrix(all.subs.rat[[1]][-1:-2,"NC"], byrow = TRUE, ncol = 10, nrow = dim(Litter.subs.df)[1])
colnames(subsNC.m) <- rownames(all.subs[[1]][-1:-2,])

# Making matrix of substrate C:P values
subsPC.m <- matrix(all.subs.rat[[1]][-1:-2,"PC"], byrow = TRUE, ncol = 10, nrow = dim(Litter.subs.df)[1])
colnames(subsPC.m) <- rownames(all.subs[[1]][-1:-2,])


# Substrate series in nitrogen units
Litter.subsN.df <- Litter.subs.df[,rownames(all.subs[[1]][-1:-2,])]*subsNC.m
Litter.subsN.df$Litter.N <- rowSums(Litter.subsN.df)
Litter.subsN.df <- cbind(Litter.subsN.df, Litter.subs.df[,c("Simulation", "days")])
write_tsv(Litter.subsN.df, "LitterSubsN_HRo.txt")

# Substrate series in phosphorus units
Litter.subsP.df <- Litter.subs.df[,rownames(all.subs[[1]][-1:-2,])]*subsPC.m
Litter.subsP.df$Litter.P <- rowSums(Litter.subsP.df)
Litter.subsP.df <- cbind(Litter.subsP.df, Litter.subs.df[,c("Simulation", "days")])
write_tsv(Litter.subsP.df, "LitterSubsP_HRo.txt")


Litter.CN <- data.frame(Litter.CN = Litter.subs.df$Litter.C/Litter.subsN.df$Litter.N,
                        Simulation = Litter.subs.df$Simulation,
                        days = Litter.subs.df$days) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN.i"), sep = "_", remove = FALSE)
write_tsv(Litter.CN, "LitterCN_HRo.txt")

Litter.CP <- data.frame(Litter.CP = Litter.subs.df$Litter.C/Litter.subsP.df$Litter.P,
                        Simulation = Litter.subs.df$Simulation,
                        days = Litter.subs.df$days) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN.i"), sep = "_", remove = FALSE)
write_tsv(Litter.CP, "LitterCP_HRo.txt")


Litter.subs.df.l <- 
  Litter.subs.df %>%
  pivot_longer(cols = rownames(all.subs[[1]])[-1:-2], names_to = "Substrate", values_to = "Abundance")

# Adding variable for the element that is regulating degradation of each substrate
subs.track.look <- data.frame(element = c("NP", "N", "C", "C", "C", "N", "N", "N", "N", "N", "P", "P"),
                              row.names = rownames((all.subs[[1]])))
Litter.subs.df.l$Subs.type <- subs.track.look[Litter.subs.df.l$Substrate,]
write_tsv(Litter.subs.df.l, "LitterSubsL_HRo.txt")


# Litter C:N and C:P time series
L.CN.df <- as.data.frame(L.CN) %>%  cbind(days = 1:length(L.CN[[1]])) %>%
  pivot_longer(cols = 1:length(L.CN), names_to = "Simulation", 
               values_to = "L.CN") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(L.CN.df, "Lcn_HRo.txt")

L.CP.df <- as.data.frame(L.CP) %>%  cbind(days = 1:length(L.CP[[1]])) %>%
  pivot_longer(cols = 1:length(L.CP), names_to = "Simulation", 
               values_to = "L.CP") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(L.CP.df, "Lcp_HRo.txt")

L.NP.df <- as.data.frame(L.NP) %>%  cbind(days = 1:length(L.NP[[1]])) %>%
  pivot_longer(cols = 1:length(L.NP), names_to = "Simulation", 
               values_to = "L.NP") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(L.NP.df, "Lnp_HRo.txt")


# Monomer time series
C.MonSer.df <- do.call("rbind", C.MonomersSeries) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(C.MonSer.df) <- NULL
C.MonSer.df$C.MonSer <- rowSums(C.MonSer.df[, c("DeadMic", "DeadEnz", "Cellulose", 
                                                "Hemicellulose", "Starch", "Chitin", 
                                                "Lignin", "Protein1", "Protein2", 
                                                "Protein3", "OrgP1", "OrgP2")])
write_tsv(C.MonSer.df, "CMonSer_HRo.txt")

C.MonSerLitt.df <- C.MonSer.df[,c("days", "Simulation", "Mechanism", "Litter.CN")]

C.MonSerLitt.df$C.MonSerLitt <- 
  rowSums(C.MonSer.df[, c("Cellulose", "Hemicellulose", "Starch", "Chitin", 
                          "Lignin", "Protein1", "Protein2", 
                          "Protein3", "OrgP1", "OrgP2")])
write_tsv(C.MonSerLitt.df, "CMonSerLitt_HRo.txt")


N.MonSer.df <- do.call("rbind", N.MonomersSeries) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(N.MonSer.df) <- NULL
N.MonSer.df$N.MonSer <- rowSums(N.MonSer.df[, c("NH4", "DeadMic", "DeadEnz", "Cellulose", 
                                                "Hemicellulose", "Starch", "Chitin", 
                                                "Lignin", "Protein1", "Protein2", 
                                                "Protein3", "OrgP1", "OrgP2")])
write_tsv(N.MonSer.df, "NMonSer_HRo.txt")

N.MonSerLitt.df <- N.MonSer.df[,c("days", "Simulation", "Mechanism", "Litter.CN")]


N.MonSerLitt.df$N.MonSerLitt <- rowSums(N.MonSer.df[, c("Cellulose", 
                                                        "Hemicellulose", "Starch", "Chitin", 
                                                        "Lignin", "Protein1", "Protein2", 
                                                        "Protein3", "OrgP1", "OrgP2")])
write_tsv(N.MonSerLitt.df, "NMonSerLitt_HRo.txt")


MonSer.df <- left_join(C.MonSer.df, N.MonSer.df, by = c("Simulation", "days", "Mechanism", "Litter.CN")) %>%
  select(Simulation, Mechanism, Litter.CN, days, C.MonSer, N.MonSer)
write_tsv(MonSer.df, "MonSer_HRo.txt")

MonSerLitt.df <- left_join(C.MonSerLitt.df, N.MonSerLitt.df, by = c("Simulation", "days", "Mechanism", "Litter.CN")) %>%
  select(Simulation, Mechanism, Litter.CN, days, C.MonSerLitt, N.MonSerLitt)
write_tsv(MonSerLitt.df, "MonSerLitt_HRo.txt")



NH4Ser.df <- cbind(as.data.frame(NH4Ser), days = 1:length(NH4Ser[[1]])) %>%
  pivot_longer(cols = 1:length(NH4Ser), names_to = "Simulation",
               values_to = "NH4") %>% 
  mutate(Simulation = str_replace(Simulation, pattern = "\\.", replacement = "-")) %>%
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)


PO4Ser.df <- cbind(as.data.frame(PO4Ser), days = 1:length(PO4Ser[[1]])) %>%
  pivot_longer(cols = 1:length(PO4Ser), names_to = "Simulation",
               values_to = "PO4") %>% 
  mutate(Simulation = str_replace(Simulation, pattern = "\\.", replacement = "-")) %>%
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)


NutSer.df <- inner_join(NH4Ser.df, PO4Ser.df)
NutSer.df.l <- NutSer.df %>% pivot_longer(cols = c("NH4", "PO4"), 
                                          names_to = "Nutrient", 
                                          values_to = "Amount") %>%
  inner_join(Biomass.df)
write_tsv(NutSer.df, "NutSer_HRo.txt"); write_tsv(NutSer.df.l, "NutSerL_HRo.txt")


N.minSer.df <- as.data.frame(N.minSer) %>%  cbind(days = 1:length(N.minSer[[1]])) %>%
  pivot_longer(cols = 1:length(N.minSer), names_to = "Simulation", 
               values_to = "N.mineral") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)


P.minSer.df <- as.data.frame(P.minSer) %>%  cbind(days = 1:length(P.minSer[[1]])) %>%
  pivot_longer(cols = 1:length(P.minSer), names_to = "Simulation", 
               values_to = "P.mineral") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

MineralSer.df <- N.minSer.df %>% left_join(P.minSer.df)
write_tsv(MineralSer.df, "MineralSer_HRo.txt")


# Monomer production series
# Substrate elemental fractions
Subs.frac <- 
  all.subs.rat[[1]] %>% mutate(Cf = C/(C+N+P), Nf = N/(C+N+P), Pf = P/(C+N+P)) %>%
  select(Cf, Nf, Pf)

MonProdSer_C.df <- do.call("rbind", MonProdSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE) %>%
  mutate(Chitin = Chitin*Subs.frac["Chitin", "Cf"],
         Lignin = Lignin*Subs.frac["Lignin", "Cf"],
         Protein1 = Protein1*Subs.frac["Protein1", "Cf"],
         Protein2 = Protein2*Subs.frac["Protein2", "Cf"],
         Protein3 = Protein3*Subs.frac["Protein3", "Cf"]) %>%
  select(-OrgP1, -OrgP2, -DeadMic, -DeadEnz) %>%
  mutate(Tot.C = rowSums(across(Cellulose:Protein3)))

rownames(MonProdSer_C.df) <- NULL


MonProdSer_N.df <-
  do.call("rbind", MonProdSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE) %>%
  mutate(Chitin = Chitin*Subs.frac["Chitin", "Nf"],
         Lignin = Lignin*Subs.frac["Lignin", "Nf"],
         Protein1 = Protein1*Subs.frac["Protein1", "Nf"],
         Protein2 = Protein2*Subs.frac["Protein2", "Nf"],
         Protein3 = Protein3*Subs.frac["Protein3", "Nf"]) %>%
  select(-OrgP1, -OrgP2, -DeadMic, -DeadEnz, -Cellulose, -Hemicellulose, -Starch) %>%
  mutate(Tot.N = rowSums(across(Chitin:Protein3)))

rownames(MonProdSer_N.df) <- NULL

MonProdSer_CN.df <-
  (MonProdSer_C.df %>% select(Tot.C, Simulation, Mechanism, Litter.CN, days)) %>%
  left_join((MonProdSer_N.df %>% select(Tot.N, Simulation, Mechanism, Litter.CN, days))) %>% 
  mutate(Tot.CN = Tot.C/Tot.N)

write_tsv(MonProdSer_CN.df, "MonProdSerCN_HRo.txt")



EnzSer.df <- do.call("rbind", EnzSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(EnzSer.df) <- NULL
#EnzSer.df <- EnzSer.df[, c(names(which(colSums(EnzSer.df[,1:40]) > 0)), "days", "Simulation")] # Leaving only relevant enzymes
write_tsv(EnzSer.df, "EnzSer_HRo.txt")


EnzSer.df.l <- 
  EnzSer.df %>%
  pivot_longer(cols = starts_with("Enz"), 
               names_to = "Enzyme", values_to = "Abundance") %>%
  inner_join(Biomass.df)

# Adding variable for the element that is regulating each enzyme
EnzSer.df.l$Enz.type <-  get(Envs[1])$out[[1]]$Enz.track.look[EnzSer.df.l$Enzyme]
write_tsv(EnzSer.df.l, "EnzSerL_HRo.txt")


# Enzyme production series
EnzProdSer.df <- do.call("rbind", EnzProdSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(EnzProdSer.df) <- NULL


EnzProdSer.df.l <- 
  EnzProdSer.df %>%
  pivot_longer(cols = starts_with("Enz"), 
               names_to = "Enzyme", values_to = "Abundance") %>%
  inner_join(Biomass.df)

# Adding variable for the element that is regulating each enzyme
EnzProdSer.df.l$Enz.type <-  get(Envs[1])$out[[1]]$Enz.track.look[EnzProdSer.df.l$Enzyme]
write_tsv(EnzProdSer.df.l, "EnzProdSerL_HRo.txt")


# Uptake series
C.UptSer.df <- do.call("rbind", C.UptSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(C.UptSer.df) <- NULL
names(C.UptSer.df)[1] <- "Taxon_Uptake_C"

N.UptSer.df <- do.call("rbind", N.UptSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(N.UptSer.df) <- NULL
names(N.UptSer.df)[1] <- "Taxon_Uptake_N"

P.UptSer.df <- do.call("rbind", P.UptSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)

rownames(P.UptSer.df) <- NULL
names(P.UptSer.df)[1] <- "Taxon_Uptake_P"

UptSer.df <- inner_join(C.UptSer.df, N.UptSer.df) %>% inner_join(P.UptSer.df) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(UptSer.df, "UptSer_HRo.txt")


N.MonomerUptake.df <- do.call("rbind", N.MonomerUptakeSeries) %>%
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE) %>%
  mutate(Tot.recycl = NH4 + DeadMic + DeadEnz, 
         Tot.litt = Chitin + Lignin + Protein1 + Protein2 + Protein3 + OrgP2)
rownames(N.MonomerUptake.df) <- NULL
write_tsv(N.MonomerUptake.df, "N_MonomerUptake_HRo.txt")


Net.CUE.df <- as.data.frame(Net.CUE_Series) %>%  cbind(days = 1:1000) %>%
  pivot_longer(cols = 1:length(Net.CUE_Series), names_to = "Simulation", 
               values_to = "CUE") %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(Net.CUE.df, "NetCUE_HRo.txt")




# Integrating biomass in time
library(pracma)


Tot.growth <- # Until 50% of mass loss
  Growth.df %>% left_join(Litter.subs.df) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  filter(Element == "C") %>% 
  summarize(Tot.C.growth = sum(Growth, na.rm = TRUE))
write_tsv(Tot.growth, "Tot_Growth_HRo.txt")

Tot.resp <- # Until 50% of mass loss
  Resp.comp.df %>% left_join(Litter.subs.df) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE))
write_tsv(Tot.resp, "TotResp_HRo.txt")

Tot.resp.comp <- 
  Resp.comp.df %>% 
  group_by(Simulation, Mechanism, Litter.CN, Source) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE))
write_tsv(Tot.resp.comp, "TotRespComp_HRo.txt")

Tot.C.Uptake <- # Until 50% of mass loss
  UptSer.df %>% left_join(Litter.subs.df) %>% 
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.Uptake = sum(Taxon_Uptake_C, na.rm = TRUE))
write_tsv(Tot.C.Uptake, "TotCuptake_HRo.txt")

Tot.CUE <- 
  Tot.growth %>% left_join(Tot.resp) %>% left_join(Tot.C.Uptake) %>% 
  mutate(CUE = (Tot.C.Uptake - Tot.resp)/Tot.C.Uptake, 
         CUE.growth = (Tot.C.growth)/(Tot.C.Uptake))
write_tsv(Tot.CUE, "TotCUE_HRo.txt")


Tot.Biomass <-
  split(Biomass.df, Biomass.df$Simulation) %>% 
  lapply(drop_na) %>%
  lapply(function(x) trapz(x$days, x$Biomass)) %>%
  unlist() %>% 
  data.frame(Tot.Biomass = ., Simulation = names(.), row.names = NULL) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(Tot.Biomass, "TotBiomass_HRo.txt")


# Integrating enzyme production in time
Tot.Enz <- 
  EnzSer.df.l %>% drop_na() %>% group_by(Simulation, Enz.type) %>%
  summarize(Tot.Enz = trapz(days, Abundance)) %>%
  pivot_wider(names_from = Enz.type, values_from = Tot.Enz) %>%
  mutate(CN = C/N, CP = C/P) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(Tot.Enz, "TotEnz_HRo.txt")


# Integrating monomer production in time
Tot.Mon <- 
  MonSer.df %>% drop_na() %>%
  pivot_longer(cols = c("C.MonSer", "N.MonSer"), names_to = "Element", values_to = "Abundance") %>% 
  group_by(Simulation, Element) %>%
  summarize(Tot.Mon = trapz(days, Abundance)) %>%
  pivot_wider(names_from = Element, values_from = Tot.Mon) %>%
  mutate(CN = C.MonSer/N.MonSer) %>% 
  separate(Simulation, into = c("Mechanism", "Litter.CN"), sep = "_", remove = FALSE)
write_tsv(Tot.Mon, "TotMon_HRo.txt")


Tot.mass.loss <- 
  Litter.subs.df %>% group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  summarize(Final.rem = min(Remaining, na.rm = TRUE))
write_tsv(Tot.mass.loss, "TotMassLoss_HRo.txt")

Days.to.50 <- 
  Litter.subs.df %>% group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>% summarize(Days.to.50 = max(days, na.rm = TRUE))
write_tsv(Days.to.50, "DaysTo50_HRo.txt")

Tot.excretion <- 
  NutSer.df %>% left_join(Litter.subs.df) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.NH4.ex = sum(NH4), Tot.PO4.ex = sum(PO4))
write_tsv(Tot.excretion, "TotExcretion_HRo.txt")



# # ------------------------------------------------------------------------------
# 
# #            Plotting
# #      
# 
# 
# # Total Biomass
# pdf("Total_Biomass_HRo.pdf", height = 3, width = 6)
# Tot.Biomass %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Tot.Biomass)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(x = "Litter C:N", y = "Total biomass (mg)") +
#   theme_classic()
# dev.off()
# 
# # Total biomass growth
# pdf("Total_Carbon_Growth_HRo.pdf", height = 3, width = 6)
# Tot.CUE %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Tot.C.growth)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(x = "Litter C:N", y = "Total C growth") +
#   theme_classic()
# dev.off()
# 
# # Total Respiration
# pdf("Total_respiration_HRo.pdf", height = 3, width = 6)
# Tot.resp %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Tot.resp)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(x = "Litter C:N", y = "Total Respiration (mg CO2)") +
#   theme_classic() 
# dev.off()
# 
# # Total Respiration per biomass
# pdf("Total_respiration_per_biomass_HRo.pdf", height = 3, width = 6)
# Tot.resp %>% left_join(Tot.Biomass) %>%
#   mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Tot.resp/Tot.Biomass)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(x = "Litter C:N", y = "Total Respiration (mg CO2) / Total Biomass (mg)") +
#   theme_classic()
# dev.off()
# 
# # Overflow respiration
# pdf("Overflow_respiration_HRo.pdf", height = 3, width = 6)
# Tot.resp.comp %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   filter(Source == "Overflow") %>%
#   ggplot(aes(x = Litter.CN, y = Tot.resp)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(y = "Overflow Respiration (mg CO2)", x = "Litter C:N") +
#   theme_classic()
# dev.off()
# 
# # Maintenance respiration
# pdf("Maintenance_respiration_HRo.pdf", height = 3, width = 6)
# Tot.resp.comp %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   filter(Source == "Maint") %>%
#   ggplot(aes(x = Litter.CN, y = Tot.resp)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(y = "Maintenance Respiration (mg CO2)", x = "Litter C:N") +
#   theme_classic()
# dev.off()
# 
# # Maintenance respiration per biomass
# pdf("Maintenance_per_biomass_HRo.pdf", height = 3, width = 6)
# Tot.resp.comp %>% left_join(Tot.Biomass) %>%
#   mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   filter(Source == "Maint") %>%
#   ggplot(aes(x = Litter.CN, y = Tot.resp/Tot.Biomass)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(y = "Maintenance Respiration (mg CO2) / Biomass (mg)") +
#   theme_classic()
# dev.off()
# 
# # Growth respiration
# pdf("Growth_respiration_HRo.pdf", height = 3, width = 6)
# Tot.resp.comp %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   filter(Source == "Growth") %>%
#   ggplot(aes(x = Litter.CN, y = Tot.resp)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(y = "Growth Respiration (mg CO2)", x = "Litter C:N") +
#   theme_classic()
# dev.off()
# 
# # Growth respiration per biomass
# pdf("Growth_respiration_per_biomass_HRo.pdf", height = 3, width = 6)
# Tot.resp.comp %>% left_join(Tot.Biomass) %>%
#   mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   filter(Source == "Growth") %>%
#   ggplot(aes(x = Litter.CN, y = Tot.resp/Tot.Biomass)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(y = "Growth Respiration (mg CO2) / Biomass (mg)") +
#   theme_classic()
# dev.off()
# 
# # Carbon-use efficiency
# pdf("CUE_HRo.pdf", height = 3, width = 6)
# Tot.CUE %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = CUE)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(y = "Carbon-use efficiency", x = "Litter C:N") +
#   theme_classic() +
#   expand_limits(y = c(0, 0.7))
# dev.off()
# 
# # Years to 50% mass loss
# pdf("Years_to_50_mass_loss_HRo.pdf", height = 3, width = 6)
# Days.to.50 %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Days.to.50/365)) +
#   geom_path(aes(color = Mechanism)) +
#   expand_limits(y = 0) +
#   labs(y = "Years to 50% mass loss", x = "Litter C:N") +
#   theme_classic()
# dev.off()
# 
# # Final remaining mass
# pdf("Final_remaining_mass_HRo.pdf", height = 3, width = 6)
# Tot.mass.loss %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Final.rem)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(y = "Remaining mass (%)", x = "Litter C:N") +
#   theme_classic()
# dev.off()
# 
# # Total excreted NH4
# pdf("Excreted_NH4_HRo.pdf", height = 3, width = 6)
# Tot.excretion %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Tot.NH4.ex)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(x = "Litter C:N", y = "Total excreted NH4") +
#   theme_classic()
# dev.off()
# 
# # Total excreted PO4
# pdf("Excreted_PO4_HRo.pdf", height = 3, width = 6)
# Tot.excretion %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = Tot.PO4.ex)) +
#   geom_path(aes(color = Mechanism)) +
#   labs(x = "Litter C:N", y = "Total excreted PO4") +
#   theme_classic()
# dev.off()
# 
# # Aggregated enzyme allocation
# pdf("Enzyme_CN_HRo.pdf", height = 3, width = 6)
# Tot.Enz %>% mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
#   ggplot(aes(x = Litter.CN, y = CN)) +
#   geom_line(aes(color = Mechanism)) +
#   labs(x = "Litter C:N", y = "Enzyme C:N") +
#   ylim(c(0, 2)) +
#   theme_classic()
# dev.off()
# 
