

#               Figures 


library(tidyverse)
library(patchwork)
library(cowplot)
library(grid)


# Loading data


# Time series
Biomass_HR <- read_tsv("Biomass_HR.txt")
Biomass_HRo <- read_tsv("Biomass_HRo.txt")
Biomass_LR <- read_tsv("Biomass_LR.txt")
Biomass_LRo <- read_tsv("Biomass_LRo.txt")
LitterSubs.L_HR <- read_tsv("LitterSubsL_HR.txt")
LitterSubs.L_HRo <- read_tsv("LitterSubsL_HRo.txt")
LitterSubs.L_LR <- read_tsv("LitterSubsL_LR.txt")
LitterSubs.L_LRo <- read_tsv("LitterSubsL_LRo.txt")
UptSer_HR <- read_tsv("UptSer_HR.txt")
UptSer_HRo <- read_tsv("UptSer_HRo.txt")
UptSer_LR <- read_tsv("UptSer_LR.txt")
UptSer_LRo <- read_tsv("UptSer_LRo.txt")
Growth_HR <- read_tsv("Growth_HR.txt")
Growth_HRo <- read_tsv("Growth_HRo.txt")
Growth_LR <- read_tsv("Growth_LR.txt")
Growth_LRo <- read_tsv("Growth_LRo.txt")
RespComp_HR <- read_tsv("RespComp_HR.txt")
RespComp_HRo <- read_tsv("RespComp_HRo.txt")
RespComp_LR <- read_tsv("RespComp_LR.txt")
RespComp_LRo <- read_tsv("RespComp_LRo.txt")
EnzProdSer_HR <- read_tsv("EnzProdSerL_HR.txt")
EnzProdSer_LR <- read_tsv("EnzProdSerL_LR.txt")
EnzProdSer_HRo <- read_tsv("EnzProdSerL_HRo.txt")
EnzProdSer_LRo <- read_tsv("EnzProdSerL_LRo.txt")
Necromass_HR <- read_tsv("Necromass_HR.txt")
Necromass_HRo <- read_tsv("Necromass_HRo.txt")
Necromass_LR <- read_tsv("Necromass_LR.txt")
Necromass_LRo <- read_tsv("Necromass_LRo.txt")
MonUptaSer_HR <- read_tsv("N_MonomerUptake_HR.txt")
MonUptaSer_HRo <- read_tsv("N_MonomerUptake_HRo.txt")
MonUptaSer_LR <- read_tsv("N_MonomerUptake_LR.txt")
MonUptaSer_LRo <- read_tsv("N_MonomerUptake_LRo.txt")


# Look-up table to change the text in the plot legend
LookMech <- c("Overflow" = "Overflow", "Stoichiometry" = "Overflow + Stoichiometry",
              "Enzyme" = "Overflow + Enzyme", "Enzyme-Uptake" = "Overflow + Enzyme + Uptake",
              "All-NoUptake" = "Overflow + Enzyme + Stoichiometry", 
              "All" = "All")
MechLev <- c("Overflow", "Overflow + Stoichiometry", 
             "Overflow + Enzyme", "Overflow + Enzyme + Uptake",
             "Overflow + Enzyme + Stoichiometry", "All")
StoichLook <- c("TRUE" = "Flexible", "FALSE" = "Fixed")
AllocLook <- c("FALSE_FALSE" = "None", "TRUE_FALSE" = "Enzyme", "TRUE_TRUE" = "Enzyme + Uptake")



# Processing time series
LitterSubs_HR <- 
  LitterSubs.L_HR %>% filter(Subs.type == "C") %>% select(-Subs.type) %>%
  pivot_wider(names_from = Substrate, values_from = Abundance)

LitterSubs_HRo <- 
  LitterSubs.L_HRo %>% filter(Subs.type == "C") %>% select(-Subs.type) %>%
  pivot_wider(names_from = Substrate, values_from = Abundance)

LitterSubs_LR <- 
  LitterSubs.L_LR %>% filter(Subs.type == "C") %>% select(-Subs.type) %>%
  pivot_wider(names_from = Substrate, values_from = Abundance)

LitterSubs_LRo <- 
  LitterSubs.L_LRo %>% filter(Subs.type == "C") %>% select(-Subs.type) %>%
  pivot_wider(names_from = Substrate, values_from = Abundance)


Tot.resp_HR <-
  RespComp_HR %>% left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE), Tot.dec.C = max(Litter.C) - min(Litter.C)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Dec.efficiency = (Tot.dec.C - Tot.resp)/Tot.dec.C,
         Community = "High (1 taxon)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.resp_HRo <-
  RespComp_HRo %>% left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE), Tot.dec.C = max(Litter.C) - min(Litter.C)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Dec.efficiency = (Tot.dec.C - Tot.resp)/Tot.dec.C,
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.resp_LR <-
  RespComp_LR %>% left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE), Tot.dec.C = max(Litter.C) - min(Litter.C)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Dec.efficiency = (Tot.dec.C - Tot.resp)/Tot.dec.C,
         Community = "Low (100 taxa)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.resp_LRo <-
  RespComp_LRo %>% left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE), Tot.dec.C = max(Litter.C) - min(Litter.C)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Dec.efficiency = (Tot.dec.C - Tot.resp)/Tot.dec.C,
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


# Total respiration
Tot.respCUE_HR <- 
  RespComp_HR %>% left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.respCUE_HRo <- 
  RespComp_HRo %>% left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.respCUE_LR <- 
  RespComp_LR %>% left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.respCUE_LRo <- 
  RespComp_LRo %>% left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.resp = sum(Resp, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))


# Total respiration components
Tot.resp.comp_HR <- 
  RespComp_HR %>% pivot_wider(names_from = Source, values_from = Resp) %>%
  left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.Maint = sum(Maint, na.rm = TRUE),
            Tot.Growth = sum(Growth, na.rm = TRUE),
            Tot.Overflow = sum(Overflow, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.resp.comp_HRo <- 
  RespComp_HRo %>% pivot_wider(names_from = Source, values_from = Resp) %>%
  left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.Maint = sum(Maint, na.rm = TRUE),
            Tot.Growth = sum(Growth, na.rm = TRUE),
            Tot.Overflow = sum(Overflow, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.resp.comp_LR <- 
  RespComp_LR %>% pivot_wider(names_from = Source, values_from = Resp) %>%
  left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.Maint = sum(Maint, na.rm = TRUE),
            Tot.Growth = sum(Growth, na.rm = TRUE),
            Tot.Overflow = sum(Overflow, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.resp.comp_LRo <- 
  RespComp_LRo %>% pivot_wider(names_from = Source, values_from = Resp) %>%
  left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.Maint = sum(Maint, na.rm = TRUE),
            Tot.Growth = sum(Growth, na.rm = TRUE),
            Tot.Overflow = sum(Overflow, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


# Total uptake
Tot.uptake_HR <- 
  UptSer_HR %>% left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), 
            Tot.N.upt = sum(Taxon_Uptake_N),
            Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.uptake_HRo <- 
  UptSer_HRo %>% left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), 
            Tot.N.upt = sum(Taxon_Uptake_N),
            Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.uptake_LR <- 
  UptSer_LR %>% left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), 
            Tot.N.upt = sum(Taxon_Uptake_N),
            Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.uptake_LRo <- 
  UptSer_LRo %>% left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), 
            Tot.N.upt = sum(Taxon_Uptake_N),
            Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.uptakeCUE_HR <- 
  UptSer_HR %>% left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), Tot.N.upt = sum(Taxon_Uptake_N), Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.uptakeCUE_HRo <- 
  UptSer_HRo %>% left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), Tot.N.upt = sum(Taxon_Uptake_N), Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.uptakeCUE_LR <- 
  UptSer_LR %>% left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), Tot.N.upt = sum(Taxon_Uptake_N), Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.uptakeCUE_LRo <- 
  UptSer_LRo %>% left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.upt = sum(Taxon_Uptake_C), Tot.N.upt = sum(Taxon_Uptake_N), Tot.P.upt = sum(Taxon_Uptake_P)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))


# Total enzyme production
Tot.EnzProd_HR <- 
  EnzProdSer_HR %>% left_join(LitterSubs_HR) %>% 
  group_by(Simulation, Mechanism, Litter.CN, Enzyme) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>% 
  summarize(Tot.enz = sum(Abundance)) %>%
  rename(EnzymeID = Enzyme) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.EnzProd_LR <- 
  EnzProdSer_LR %>% left_join(LitterSubs_LR) %>% 
  group_by(Simulation, Mechanism, Litter.CN, Enzyme) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>% 
  summarize(Tot.enz = sum(Abundance)) %>%
  rename(EnzymeID = Enzyme) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.EnzProd_HRo <- 
  EnzProdSer_HRo %>% left_join(LitterSubs_HRo) %>% 
  group_by(Simulation, Mechanism, Litter.CN, Enzyme) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>% 
  summarize(Tot.enz = sum(Abundance)) %>%
  rename(EnzymeID = Enzyme) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.EnzProd_LRo <- 
  EnzProdSer_LRo %>% left_join(LitterSubs_LRo) %>% 
  group_by(Simulation, Mechanism, Litter.CN, Enzyme) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>% 
  summarize(Tot.enz = sum(Abundance)) %>%
  rename(EnzymeID = Enzyme) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


# Total necromass
Tot.necro_HR <- 
  Necromass_HR %>% left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.necromass = sum(DeadMic),
            Tot.deadEnz = sum(DeadEnz)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "CO[2]",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.necro_LR <- 
  Necromass_LR %>% left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.necromass = sum(DeadMic),
            Tot.deadEnz = sum(DeadEnz)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "CO[2]",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.necro_HRo <- 
  Necromass_HRo %>% left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.necromass = sum(DeadMic),
            Tot.deadEnz = sum(DeadEnz)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.necro_LRo <- 
  Necromass_LRo %>% left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.necromass = sum(DeadMic),
            Tot.deadEnz = sum(DeadEnz)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


# Total uptake per monomer type
Tot.monUptake_HR <- 
  MonUptaSer_HR %>% select(-Cellulose, -Hemicellulose, -Starch) %>%
  left_join(LitterSubs_HR, by = c("Simulation", "days", "Mechanism", "Litter.CN")) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(NH4 = sum(NH4), 
            DeadMic = sum(DeadMic),
            DeadEnz = sum(DeadEnz),
            Chitin = sum(Chitin),
            Lignin = sum(Lignin),
            Protein1 = sum(Protein1),
            Protein2 = sum(Protein2),
            Protein3 = sum(Protein3),
            Tot.recycl = sum(Tot.recycl),
            Tot.litt = sum(Tot.litt)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "CO[2]",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.monUptake_HRo <- 
  MonUptaSer_HRo %>% select(-Cellulose, -Hemicellulose, -Starch) %>%
  left_join(LitterSubs_HRo, by = c("Simulation", "days", "Mechanism", "Litter.CN")) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(NH4 = sum(NH4), 
            DeadMic = sum(DeadMic),
            DeadEnz = sum(DeadEnz),
            Chitin = sum(Chitin),
            Lignin = sum(Lignin),
            Protein1 = sum(Protein1),
            Protein2 = sum(Protein2),
            Protein3 = sum(Protein3),
            Tot.recycl = sum(Tot.recycl),
            Tot.litt = sum(Tot.litt)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.monUptake_LR <- 
  MonUptaSer_LR %>% select(-Cellulose, -Hemicellulose, -Starch) %>%
  left_join(LitterSubs_LR, by = c("Simulation", "days", "Mechanism", "Litter.CN")) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(NH4 = sum(NH4), 
            DeadMic = sum(DeadMic),
            DeadEnz = sum(DeadEnz),
            Chitin = sum(Chitin),
            Lignin = sum(Lignin),
            Protein1 = sum(Protein1),
            Protein2 = sum(Protein2),
            Protein3 = sum(Protein3),
            Tot.recycl = sum(Tot.recycl),
            Tot.litt = sum(Tot.litt)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "CO[2]",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))

Tot.monUptake_LRo <- 
  MonUptaSer_LRo %>% select(-Cellulose, -Hemicellulose, -Starch) %>%
  left_join(LitterSubs_LRo, by = c("Simulation", "days", "Mechanism", "Litter.CN")) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(NH4 = sum(NH4), 
            DeadMic = sum(DeadMic),
            DeadEnz = sum(DeadEnz),
            Chitin = sum(Chitin),
            Lignin = sum(Lignin),
            Protein1 = sum(Protein1),
            Protein2 = sum(Protein2),
            Protein3 = sum(Protein3),
            Tot.recycl = sum(Tot.recycl),
            Tot.litt = sum(Tot.litt)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


# Total biomass
Tot.biomass_HR <- 
  Biomass_HR %>% left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.biomass = sum(Biomass, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.biomass_LR <- 
  Biomass_LR %>% left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.biomass = sum(Biomass, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.biomass_HRo <- 
  Biomass_HRo %>% left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.biomass = sum(Biomass, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


Tot.biomass_LRo <- 
  Biomass_LRo %>% left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.biomass = sum(Biomass, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")),
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev),
         Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


# Total growth
Tot.growth_HR <- 
  Growth_HR %>% pivot_wider(names_from = Element, values_from = Growth) %>%
  left_join(LitterSubs_HR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.growth = sum(C, na.rm = TRUE),
            Tot.N.growth = sum(N, na.rm = TRUE),
            Tot.P.growth = sum(P, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.growth_HRo <- 
  Growth_HRo %>% pivot_wider(names_from = Element, values_from = Growth) %>%
  left_join(LitterSubs_HRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.growth = sum(C, na.rm = TRUE),
            Tot.N.growth = sum(N, na.rm = TRUE),
            Tot.P.growth = sum(P, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.growth_LR <- 
  Growth_LR %>% pivot_wider(names_from = Element, values_from = Growth) %>%
  left_join(LitterSubs_LR) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.growth = sum(C, na.rm = TRUE),
            Tot.N.growth = sum(N, na.rm = TRUE),
            Tot.P.growth = sum(P, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))

Tot.growth_LRo <- 
  Growth_LRo %>% pivot_wider(names_from = Element, values_from = Growth) %>%
  left_join(LitterSubs_LRo) %>%
  group_by(Simulation, Mechanism, Litter.CN) %>%
  mutate(Remaining = (Litter.C/max(Litter.C, na.rm = TRUE))*100) %>%
  filter(Remaining > 50) %>%
  summarize(Tot.C.growth = sum(C, na.rm = TRUE),
            Tot.N.growth = sum(N, na.rm = TRUE),
            Tot.P.growth = sum(P, na.rm = TRUE)) %>%
  mutate(Litter.CN = as.numeric(str_remove(Litter.CN, "CN")))


# Total CUE
Tot.CUE_HR <- 
  Tot.growth_HR %>% left_join(Tot.respCUE_HR) %>% left_join(Tot.uptakeCUE_HR) %>% 
  mutate(CUE = (Tot.C.upt - Tot.resp)/Tot.C.upt, 
         CUE.growth = (Tot.C.growth)/(Tot.C.upt),
         Community = "High (1 taxon)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev)) %>%
  mutate(Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))



Tot.CUE_HRo <- 
  Tot.growth_HRo %>% left_join(Tot.respCUE_HRo) %>% left_join(Tot.uptakeCUE_HRo) %>% 
  mutate(CUE = (Tot.C.upt - Tot.resp)/Tot.C.upt, 
         CUE.growth = (Tot.C.growth)/(Tot.C.upt),
         Community = "High (1 taxon)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev)) %>%
  mutate(Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))



Tot.CUE_LR <- 
  Tot.growth_LR %>% left_join(Tot.respCUE_LR) %>% left_join(Tot.uptakeCUE_LR) %>% 
  mutate(CUE = (Tot.C.upt - Tot.resp)/Tot.C.upt, 
         CUE.growth = (Tot.C.growth)/(Tot.C.upt),
         Community = "Low (100 taxa)",
         Overflow = "CO2",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev)) %>%
  mutate(Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))



Tot.CUE_LRo <- 
  Tot.growth_LRo %>% left_join(Tot.respCUE_LRo) %>% left_join(Tot.uptakeCUE_LRo) %>% 
  mutate(CUE = (Tot.C.upt - Tot.resp)/Tot.C.upt, 
         CUE.growth = (Tot.C.growth)/(Tot.C.upt),
         Community = "Low (100 taxa)",
         Overflow = "DOC",
         Mechanism = factor(LookMech[Mechanism], levels = MechLev)) %>%
  mutate(Stoichiometry = StoichLook[as.character(Mechanism %in% MechLev[c(2, 5:6)])],
         Enzyme = Mechanism %in% MechLev[3:6],
         Uptake = Mechanism %in% MechLev[c(4, 6)]) %>%
  unite(Enzyme, Uptake, col = "EnzUpt", remove = FALSE) %>%
  mutate(EnzUpt = AllocLook[as.character(EnzUpt)]) %>%
  mutate(EnzUpt = factor(EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake")))


# ------------------------------------------------------------------------------

#            Plotting
#      

# Defining colorblind-friendly palette
palette <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#000000", "#CC79A7", "#F0E442")
palette2 <- c("#0072B2", "#E69F00", "#000000", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#F0E442")

grid.data <-
  data.frame(Stoichiometry = c(rep(c("Fixed", "Flexible"), 2)),
             Community = c(rep("High (1 taxon)", 2), rep("Low (100 taxa)", 2)),
             grid = c("(a)", "(b)", "(c)", "(d)"))

grid.data2 <-
  data.frame(Overflow = c(rep(c("CO[2]", "DOC"), 2)),
             Community = c(rep("High~(1~taxon)", 2), rep("Low~(100~taxa)", 2)),
             grid = c("(a)", "(b)", "(c)", "(d)"))

grid.data3 <-
  data.frame(Overflow = c(rep(c("CO[2]", "DOC"), 2)),
             Community = c(rep("High~(1~taxon)", 2), rep("Low~(100~taxa)", 2)),
             grid = c("(a)", "(c)", "(b)", "(d)"))

grid.data4 <-
  data.frame(Community = c(rep("High~(1~taxon)", 3), rep("Low~(100~taxa)", 3)),
             EnzUpt = rep(c("None", "Enzyme", "Enzyme + Uptake"), 2),
             grid = c("(a)", "(c)", "(e)", "(b)", "(d)", "(f)"))

grid.data4$EnzUpt <- factor(grid.data4$EnzUpt, levels = c("None", "Enzyme", "Enzyme + Uptake"))

Replot <- 
  Tot.CUE_HR %>% filter(Stoichiometry == "Fixed", EnzUpt == "None") %>%
  mutate(Community = "Low~(100~taxa)", Overflow = "CO[2]")

# Carbon-use efficiency
pdf("Figure4.pdf", height = 5, width = 6.5)
rbind(
  Tot.CUE_HR, Tot.CUE_HRo, Tot.CUE_LR, Tot.CUE_LRo
) %>% mutate(Overflow = recode(Overflow, "CO2" = "CO[2]", "DOC" = "DOC"),
             Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)")) %>%
  ggplot(aes(x = Litter.CN, y = CUE)) +
  geom_path(data = Replot, color = "gray") +
  geom_path(aes(linetype = Stoichiometry, color = EnzUpt)) +
  geom_text(data = grid.data2, aes(label = grid), x = 12, y = 0.74) +
  labs(y = "Carbon-use efficiency", x = "Litter C:N", subtitle = "Overflow nature", color = "Nutrient allocation") +
  theme_bw() +
  ylim(c(0.3, 0.75)) +
  scale_color_manual(values = palette2) +
  facet_grid(Community~Overflow, labeller = "label_parsed") +
  theme(legend.position = c(0.74, 0.63),
        legend.background = element_blank(),
        legend.key.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.spacing = unit(0.1, "cm"),
        legend.box = "horizontal",
        panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        plot.margin = unit(c(0.1,1,0.2,0.2), "cm"))

grid.text(label = "Community redundancy", 
          x = unit(0.97, "npc"), y = 0.73, rot = 270, gp = gpar(fontsize = 11))

dev.off()


# Total respiration
pdf("Figure5.pdf", height = 3.3, width = 6.5)
rbind(
  (Tot.resp_HRo %>% mutate(Overflow = "DOC")),
  (Tot.resp_HR %>% mutate(Overflow = "CO2")),
  (Tot.resp_LRo %>% mutate(Overflow = "DOC")),
  (Tot.resp_LR %>% mutate(Overflow = "CO2"))
  
) %>% filter(Stoichiometry == "Fixed") %>%
  mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
  mutate(Overflow = recode(Overflow, "CO2" = "CO[2]", "DOC" = "DOC"),
         Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)")) %>%
  ggplot(aes(x = Litter.CN, y = Tot.resp/10000)) +
  geom_path(aes(linetype = Overflow, color = EnzUpt)) +
  geom_text(data = (grid.data3 %>% filter(Overflow == "CO[2]")), aes(label = grid), x = 12, y = 115) +
  labs(x = "Litter C:N", y = expression(Total~Respiration~(mg~CO[2])), 
       subtitle = "Community redundancy", color = "Nutrient allocation") +
  facet_grid(~Community, labeller = "label_parsed") +
  scale_color_manual(values = palette2) +
  scale_linetype_discrete(labels = c(expression(CO[2]), "DOC")) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2),
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 9),
        legend.spacing = unit(0.2, "cm"),
        legend.text = element_text(size = 8),
        plot.margin = unit(c(0.1,1,0.2,0.2), "cm"),
        panel.grid = element_blank())

dev.off()


# Respiration components
rbind(
  (Tot.resp.comp_HRo %>% mutate(Overflow = "DOC")),
  (Tot.resp.comp_HR %>% mutate(Overflow = "CO2")),
  (Tot.resp.comp_LR %>% mutate(Overflow = "CO2")),
  (Tot.resp.comp_LRo %>% mutate(Overflow = "DOC"))
) %>% pivot_longer(cols = Tot.Maint:Tot.Overflow, names_to = "Source", values_to = "Resp") %>%
  filter(Stoichiometry == "Fixed") %>%
  mutate(Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)"),
         EnzUpta = recode(EnzUpt, "None" = "None", "Enzyme" = "Enzyme", "Enzyme + Uptake" = "Enzyme~+~Uptake")) %>%
  
  ggplot(aes(x = Litter.CN, y = Resp/10000)) +
  geom_path(aes(color = Source, linetype = Overflow)) +
  geom_text(data = grid.data4, aes(label = grid), x = 12, y = 78) +
  coord_cartesian(clip = "off") +
  facet_grid(EnzUpt~Community, labeller = "label_parsed") +
  scale_color_manual(values = palette2, labels = c("Growth", "Maintenance", "Overflow")) +
  scale_linetype_discrete(labels = c(expression(CO[2]), "DOC")) +
  labs(y = expression(Respiration~(mg~CO[2])), 
       subtitle = "Community redundancy", x = "Litter C:N") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 3.5, 0.2, 0.2), "cm"),
        legend.position = c(1.2, 0.3))

grid.text(label = "Nutrient allocation", 
          x = unit(0.83, "npc"), y = 0.8, rot = 270, gp = gpar(fontsize = 11))


P3.3 <- recordPlot()

pdf("Figure6.pdf", height = 6, width = 7)
P3.3
dev.off()


# Overflow respiration 
pdf("FigureS1.pdf", height = 5, width = 6.5)
rbind(
  (Tot.resp.comp_HRo %>% mutate(Overflow = "DOC")),
  (Tot.resp.comp_HR %>% mutate(Overflow = "CO2")),
  (Tot.resp.comp_LR %>% mutate(Overflow = "CO2")),
  (Tot.resp.comp_LRo %>% mutate(Overflow = "DOC"))
) %>% pivot_longer(cols = Tot.Maint:Tot.Overflow, names_to = "Source", values_to = "Resp") %>%
  filter(Source == "Tot.Overflow") %>%
  mutate(EnzUpta = recode(EnzUpt, "None" = "None", "Enzyme" = "Enzyme", "Enzyme + Uptake" = "Enzyme~+~Uptake"),
         Source = recode(Source, "Tot.Growth" = "Growth", "Tot.Maint" = "Maintenance", "Tot.Overflow" = "Overflow"),
         Overflow = recode(Overflow, "DOC" = "DOC", "CO2" = "CO[2]"),
         Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)")) %>%
  
  ggplot(aes(x = Litter.CN, y = Resp/10000)) +
  geom_path(aes(color = EnzUpt, linetype = Stoichiometry)) +
  geom_text(data = grid.data2, aes(label = grid), x = 12, y = 37) +
  scale_color_manual(values = palette2) +
  facet_grid(Community~Overflow, labeller = "label_parsed") +
  labs(x = "Litter C:N", y = expression(Overflow~respiration~(mg~CO[2])), color = "Nutrient allocation", subtitle = "Overflow nature") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        legend.position = c(0.63, 0.75),
        legend.spacing = unit(0.1, "cm"),
        plot.margin = unit(c(0.1, 1, 0.2, 0.2), "cm"))

grid.text(label = "Community redundancy", 
          x = unit(0.97, "npc"), y = 0.72, rot = 270, gp = gpar(fontsize = 11))

dev.off()


P1 <- 
  EnzProdSer_HR %>% filter(Litter.CN == "CN90", Mechanism == "Enzyme") %>% select(-Enz.type) %>%
  pivot_wider(names_from = Enzyme, values_from = Abundance) %>%
  ggplot(aes(x = days, y = Enz003/Enz002)) +
  geom_line() +
  geom_area(data = (RespComp_HR %>% pivot_wider(names_from = Source, values_from = Resp) %>%
                      filter(Litter.CN == "CN90", Mechanism == "Enzyme") %>%
                      mutate(PO = ifelse(Overflow > 0, 100, 0))),
            aes(x = days, y = PO),  alpha = 0.3) +
  labs(x = "Time (days)", y = "Enzyme production C:N", tag = "(a)") +
  theme_bw() +
  theme(panel.grid = element_blank())

P2 <- 
  EnzProdSer_HR %>% filter(Litter.CN == "CN90", Mechanism == "Enzyme") %>% select(-Enz.type) %>%
  pivot_wider(names_from = Enzyme, values_from = Abundance) %>%
  ggplot(aes(x = days, y = Enz003/Enz002)) +
  geom_line() +
  geom_area(data = (RespComp_HR %>% pivot_wider(names_from = Source, values_from = Resp) %>%
                      filter(Litter.CN == "CN90", Mechanism == "Enzyme") %>%
                      mutate(PO = ifelse(Overflow > 0, 2, 0))),
            aes(x = days, y = PO),  alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Time (days)", y = "Enzyme production C:N", tag = "(b)") +
  ylim(c(0, 2)) + xlim(c(0, 500)) +
  theme_bw() +
  theme(panel.grid = element_blank())

pdf("FigureS2.pdf", height = 4, width = 10)
plot_grid(P1, P2, align = "h")
dev.off()


# Maintenance respiration
grid.data6 <- 
  data.frame(Stoichiometry = c(rep("Fixed", 2), rep("Flexible", 2)),
             Overflow = rep(c("CO[2]", "DOC")),
             grid = c("(a)", "(b)", "(c)", "(d)"))

pdf("FigureS3.pdf", height = 5, width = 7.5)
rbind(
  (Tot.resp.comp_HRo %>% mutate(Overflow = "DOC")),
  (Tot.resp.comp_HR %>% mutate(Overflow = "CO2")),
  (Tot.resp.comp_LR %>% mutate(Overflow = "CO2")),
  (Tot.resp.comp_LRo %>% mutate(Overflow = "DOC"))
) %>% pivot_longer(cols = Tot.Maint:Tot.Overflow, names_to = "Source", values_to = "Resp") %>%
  filter(Source == "Tot.Maint") %>%
  mutate(EnzUpta = recode(EnzUpt, "None" = "None", "Enzyme" = "Enzyme", "Enzyme + Uptake" = "Enzyme~+~Uptake"),
         Source = recode(Source, "Tot.Growth" = "Growth", "Tot.Maint" = "Maintenance", "Tot.Overflow" = "Overflow"),
         Overflow = recode(Overflow, "DOC" = "DOC", "CO2" = "CO[2]")) %>%
  
  ggplot(aes(x = Litter.CN, y = Resp/10000)) +
  geom_path(aes(color = EnzUpt, linetype = Community)) +
  geom_text(data = grid.data6, aes(label = grid), x = 12, y = 85) +
  scale_color_manual(values = palette2) +
  facet_grid(Stoichiometry~Overflow, labeller = "label_parsed") +
  labs(x = "Litter C:N", y = expression(Maintenance~respiration~(mg~CO[2])), color = "Nutrient allocation", subtitle = "Overflow nature") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1.2, 0.3),
        plot.margin = unit(c(0.1, 4, 0.2, 0.2), "cm"))

grid.text(label = "Stoichiometry", 
          x = unit(0.81, "npc"), y = 0.79, rot = 270, gp = gpar(fontsize = 11))

dev.off()


# Total carbon uptake
pdf("FigureS4.pdf", height = 3.3, width = 6.5)
rbind(
Tot.uptake_HRo, Tot.uptake_HR, Tot.uptake_LR, Tot.uptake_LRo
) %>% filter(Stoichiometry == "Fixed") %>% 
  mutate(Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)")) %>%
  mutate(Litter.CN = as.numeric(str_extract(Litter.CN, "[0-9]{2}"))) %>%
  ggplot(aes(x = Litter.CN, y = Tot.C.upt/10000)) +
  geom_path(aes(color = EnzUpt, linetype = Overflow)) +
  geom_text(data = (grid.data3 %>% filter(Overflow == "CO[2]")), aes(label = grid), x = 12, y = 365) +
  facet_grid(~Community, labeller = "label_parsed") +
  labs(x = "Litter C:N", y = "Total C uptake (mg)", color = "Nutrient allocation", subtitle = "Community redundancy") +
  ylim(c(170, 370)) +
  scale_color_manual(values = palette2) +
  scale_linetype_discrete(labels = c(expression(CO[2]), "DOC")) +
  theme_bw() +
  theme(legend.position = c(0.1, 0.5),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        legend.spacing = unit(0.2, "cm"),
        legend.key.size = unit(0.4, "cm"),
        panel.grid = element_blank())
dev.off()


# Total enzyme production
pdf("FigureS5.pdf", height = 5, width = 7.5)
rbind(
  Tot.EnzProd_HR, Tot.EnzProd_HRo, Tot.EnzProd_LR, Tot.EnzProd_LRo
) %>% mutate(Overflow = recode(Overflow, "CO2" = "CO[2]"),
             Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)")) %>%
  group_by(Simulation, EnzUpt, Stoichiometry, Overflow, Community, Litter.CN) %>%
  summarize(Tot.enz = sum(Tot.enz)) %>%
  ggplot(aes(x = Litter.CN, y = Tot.enz/10000)) +
  geom_line(aes(color = EnzUpt, linetype = Stoichiometry)) +
  geom_text(data = grid.data2, aes(label = grid), x = 12, y = 3.3) +
  facet_grid(Community~Overflow, labeller = "label_parsed") +
  labs(x = "Litter C:N", y = "Total enzyme production (mg)", subtitle = "Overflow nature", color = "Nutrient allocation") +
  scale_color_manual(values = palette2) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(1.2, 0.25),
        plot.margin = unit(c(0.1,4.1,0.2,0.2), "cm"))

grid.text(label = "Community redundancy", 
          x = unit(0.81, "npc"), y = 0.73, rot = 270, gp = gpar(fontsize = 11))

dev.off()


# Necromass
pdf("FigureS6.pdf", heigh = 5, width = 6)
rbind(
  Tot.necro_HR, Tot.necro_HRo, Tot.necro_LR, Tot.necro_LRo
) %>%
  mutate(Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)")) %>%
  ggplot(aes(x = Litter.CN, y = Tot.necromass/10000000)) +
  geom_text(data = grid.data2, aes(label = grid), x = 12, y = 2.3) +
  geom_line(aes(color = EnzUpt, linetype = Stoichiometry)) +
  labs(x = "Litter C:N", y = "Total microbial necromass (g)", color = "Nutrient allocation", subtitle = "Overflow nature") +
  scale_color_manual(values = palette2) +
  facet_grid(Community~Overflow, labeller = "label_parsed") +
  expand_limits(y = 0) +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        legend.box = "horizontal",
        legend.position = c(0.25, 0.8),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(0.1, 1, 0.2, 0.2), "cm"))

grid.text(label = "Community redundancy", 
          x = unit(0.97, "npc"), y = 0.72, rot = 270, gp = gpar(fontsize = 11))

dev.off()


# Nitrogen uptake origin
pdf("FigureS7.pdf", height = 5, width = 6)
rbind(
  Tot.monUptake_HR, Tot.monUptake_HRo, Tot.monUptake_LR, Tot.monUptake_LRo
) %>% 
  pivot_longer(NH4:Protein3, names_to = "Source", values_to = "Amount") %>%
  mutate(Community = recode(Community, "High (1 taxon)" = "High~(1~taxon)", "Low (100 taxa)" = "Low~(100~taxa)")) %>%
  filter(Source == "DeadMic") %>%
  ggplot(aes(x = Litter.CN, y = Amount/10000)) +
  geom_text(data = grid.data2, aes(label = grid), x = 12, y = 3.4) +
  geom_line(aes(color = EnzUpt, linetype = Stoichiometry)) +
  scale_color_manual(values = palette2) +
  facet_grid(Community~Overflow, labeller = "label_parsed") +
  labs(x = "Litter C:N", y = "N uptake from necromass (mg)", color = "Nutrient allocation", subtitle = "Overflow nature") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.box = "horizontal",
        legend.position = c(0.25, 0.82),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(0.1, 1, 0.2, 0.2), "cm"))

grid.text(label = "Community redundancy", 
          x = unit(0.97, "npc"), y = 0.72, rot = 270, gp = gpar(fontsize = 11))
dev.off()


# Total biomass
pdf("FigureS8.pdf", heigh = 3.5, width = 8)
rbind(
  Tot.biomass_HR, Tot.biomass_HRo, Tot.biomass_LR, Tot.biomass_LRo 
) %>% filter(Stoichiometry == "Fixed") %>%
  mutate(Overflow = recode(Overflow, "DOC" = "DOC", "CO2" = "CO[2]")) %>%
  ggplot(aes(x = Litter.CN, y = Tot.biomass/10000000)) +
  geom_path(aes(color = EnzUpt, linetype = Community)) +
  geom_text(data = grid.data6[1:2,], aes(label = grid), x = 12, y = 9.5) +
  scale_color_manual(values = palette2) +
  facet_grid(~Overflow, labeller = "label_parsed") +
  labs(x = "Litter C:N", y = "Total biomass (g)", color = "Nutrient allocation", subtitle = "Overflow nature") +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()


