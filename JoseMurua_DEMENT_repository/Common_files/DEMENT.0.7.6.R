###############################################################################
sum.grid <- function(x, rn, grid.size) {
  y <- NULL
  for (i in 1:dim(x)[2]) y <- cbind(y, rowSums(matrix(x[, i], ncol = grid.size)))
  dimnames(y) <- list(rn, colnames(x))  
  y
}

###############################################################################
expand <- function(x, grid.size) {
  matrix(rep(t(x), grid.size), 
         ncol = dim(x)[2], 
         byrow = T, 
         dimnames = list(rep(rownames(x), grid.size), colnames(x)))
}

###################################### Plot grid concentrations; grid is a vector with x changing fastest
plot.grid <- function(grid, x, zlim = NULL, title = "") {
  grid.mat <- matrix(grid, byrow = T, ncol = x)
  if(!is.null(zlim)) zlim <- c(0, zlim)
  myImagePlot(grid.mat, zlim = zlim, title = title)
}

# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <- c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1, length = 256),  # Red
                    seq(0,1, length = 256),  # Green
                    seq(1,0, length = 256))  # Blue
  ColorLevels <- seq(min, max, length = length(ColorRamp))
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  # Data Map
  image(1:length(xLabels), 1:length(yLabels), t(x), col = ColorRamp, xlab = "",
        ylab = "", axes = FALSE, zlim = c(min,max))
  if( !is.null(title) ){
    title(main = title)
  }
  #rect(0.5, 0.5, length(xLabels) + 0.5, length(yLabels) + 0.5)
}
############################################################################

MakePlots <- function(out) {
  par(cex.axis = 2, 
      cex.lab = 2.5, 
      cex.main = 3, 
      font.lab = 2, 
      font.axis = 2, 
      lwd = 2, 
      mar = c(5, 7, 4, 2) + 0.1, 
      tcl = 0.4, 
      las = 1, 
      mgp = c(3, 0.75, 0), 
      mfrow = c(2, 5))
  
  RGB.palette <- colorRampPalette(c("red", "yellow", "blue"), space = "rgb")
  
  n_taxa <- dim(out$MicrobesSeries)[2]
  colors.RGB <- RGB.palette(out$params["Enz_per_taxon_max", ]+1)
  plot(0:(length(out$RespSeries) - 1), out$Mic_Sum[, "C"], 
       ylim = c(0, max(out$MicrobesSeries/out$grid.size, na.rm = T)), 
       ylab = "", 
       main = "Microbial C", 
       xlab = "Day", 
       type = "n", 
       axes = F, 
       frame.plot = T)
  for (i in 1:n_taxa){
    lines(0:(length(out$RespSeries) - 1), 
          out$MicrobesSeries[,i]/out$grid.size, 
          col = colors.RGB[(rowSums(out$EnzGenes)+1)[i]], 
          lwd=2)
  }
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  n_subs <- length(dimnames(out$SubstratesSeries)[2][[1]])
  color.sub <- c("brown", "darkblue", "goldenrod", "salmon", "orange", "pink", 
                 "gray", "purple", "blue", "turquoise", "forestgreen", "green")
  plot(0:(length(out$RespSeries) - 1), out$Substrate_Sum,
       ylim = c(0, max(out$SubstratesSeries/out$grid.size,na.rm=T)),
       main = expression(bold(paste("Substrate C (mg ", cm^-3, ")"))),
       ylab = NA, xlab = "Day", type = "n", axes = F, frame.plot = T)
  axis(1, lwd = 0, lwd.ticks = 2); axis(2, lwd = 0, lwd.ticks = 2)
  for (i in 1:n_subs){
    lines(0:(length(out$RespSeries) - 1),
          out$SubstratesSeries[,i]/out$grid.size,
          col = color.sub[i],lwd=2)
  }
  legend("topright", 
         legend = c("Dead microbe", "Inactive enzyme", "Cellulose", 
                    "Hemicellulose", "Starch", "Chitin", "Lignin", 
                    "Protein 1", "Protein 2", "Protein 3", 
                    "Phospholipid", "Nucleic acid"), 
         lty = 1,
         col = color.sub,
         box.lty = 0, 
         cex = 1.1)
  
  plot(0:(length(out$RespSeries) - 1), 
       out$Mic_Sum[, "C"]/out$Mic_Sum[, "P"], 
       ylim = c(0, max(out$MicrobesPSeries, na.rm = T)), 
       ylab = NA, 
       main = "Microbial C:P", 
       xlab = "Day", 
       type = "l", 
       axes = F, 
       frame.plot = T)
  for (i in 1:n_taxa){
    lines(0:(length(out$RespSeries) - 1),
          out$MicrobesPSeries[,i],
          col = colors.RGB[(rowSums(out$EnzGenes)+1)[i]])
  }
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  plot(0:(length(out$RespSeries) - 1),
       out$Mic_Sum[,"C"]/out$Mic_Sum[, "N"], 
       ylim = c(0, max(out$MicrobesNSeries,na.rm=T)), 
       ylab = NA, 
       main = "Microbial C:N", 
       xlab = "Day", 
       type = "l", 
       axes = F, 
       frame.plot = T)
  for (i in 1:n_taxa){
    lines(0:(length(out$RespSeries) - 1), 
          out$MicrobesNSeries[, i], 
          col = colors.RGB[(rowSums(out$EnzGenes)+1)[i]])
  }
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  plot(0:(length(out$RespSeries) - 1), 
       out$Mic_Sum[, "N"]/out$Mic_Sum[, "P"], 
       ylim = c(0, max(out$MicrobesNPSeries, na.rm = T)), 
       ylab = NA, 
       main = "Microbial N:P", 
       xlab = "Day", 
       type = "l", 
       axes = F, 
       frame.plot = T)
  for (i in 1:n_taxa){
    lines(0:(length(out$RespSeries) - 1),
          out$MicrobesNPSeries[, i], 
          col = colors.RGB[(rowSums(out$EnzGenes)+1)[i]])
  }
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  # Color code enzymes according to the main substrate they degrade
  enz.potential <- out$ReqEnz[[1]]*out$Vmax
  colors <- color.sub[apply(enz.potential, 2, which.max)]
  plot(0:(length(out$RespSeries) - 1), 
       out$Enzyme_Sum, 
       ylim = c(0, max(out$EnzymesSeries/out$grid.size, na.rm = T)), 
       ylab = NA, 
       main = "Enzyme", 
       xlab = "Day", 
       type = "n", 
       axes = F, 
       frame.plot = T)
  for (i in 1:length(colors)){
    lines(0:(length(out$RespSeries) - 1), 
          out$EnzymesSeries[,i]/out$grid.size, 
          col = colors[i], 
          lwd = 2)
  }
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  plot(0:(length(out$RespSeries) - 1),
       out$NH4Series,
       ylim = c(0, max(cbind(out$NH4Series/out$grid.size,out$PO4Series/out$grid.size))),
       ylab = NA, 
       main = "NH4 | PO4",
       xlab = "Day",
       type = "n",
       axes = F,
       frame.plot = T)
  lines(0:(length(out$RespSeries) - 1),
        out$NH4Series/out$grid.size,
        col = "orange")
  lines(0:(length(out$RespSeries) - 1),
        out$PO4Series/out$grid.size,
        col = "purple")
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  n_monomers <- length(dimnames(out$C.MonomersSeries)[2][[1]])
  colors <- c("black", "black", color.sub)
  plot(0:(length(out$RespSeries) - 1), 
       out$Monomer_Sum, 
       ylim = c(0, max(out$C.MonomersSeries/out$grid.size, na.rm = T)),
       ylab = NA, 
       main = "Monomer C",
       xlab = "Day",
       type = "n",
       axes = F, 
       frame.plot = T)
  for (i in 1:n_monomers){
    lines(0:(length(out$RespSeries) - 1),
          out$C.MonomersSeries[,i]/out$grid.size,
          col=colors[i])
  }
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  plot(0:(length(out$RespSeries) - 1),
       out$RespSeries/out$grid.size,
       ylim = c(0, max(out$RespSeries/out$grid.size)),
       ylab = NA,
       main = "Respiration",
       xlab = "Day",
       type="l",
       axes = F,
       frame.plot = T)
  axis(1, lwd = 2); axis(2, lwd = 2)
  
  index <- is.finite(log10(apply(out$MicrobesSeries, 2, max)))
  plot(rowSums(out$EnzGenes)[index], 
       log10(apply(out$MicrobesSeries, 2, max)/out$grid.size)[index],
       main = expression(bold(paste(Log[10], " max density"))),
       ylab = NA,
       xlab = "Number of enzyme genes",
       pch = 19,
       cex = 2,
       axes = F,
       frame.plot=T,
       #ylim=c(-2,1),
       col = colors.RGB[1+as.vector(rowSums(out$EnzGenes)[index])])
  axis(1, lwd = 2);  axis(2, lwd = 2)
  #abline(lm(log10(apply(out$MicrobesSeries, 2, max)/out$grid.size)[index]~rowSums(out$EnzGenes)[index]), lwd = 1.75)
  
} # End of MakePlots() ########################################################


###############################################################################
RunPulse <- function(
    params,
    timestamp,
    rng.seed,
    grid.size,
    Microbes,
    Substrates,
    SubInput,
    Enzymes,
    Monomers,
    MonInput,
    MonomersProduced,
    ReprodNew,
    Colonization.reset,
    Ea,
    Vmax0,
    Km0,
    ReqEnz,
    EnzGenes,
    EnzProdInduce,
    EnzProdConstit,
    Enz.track.look,
    Upt.track.look,
    UptakeGenes,
    UptakeGenesForEnz,
    Uptake_ReqEnz,
    EnzAttrib,
    Uptake_Ea,
    Uptake_Vmax0,
    Uptake_Km0,
    CUE.ref,
    OptimalRatios,
    MinRatios,
    RangeRatios,
    fb,
    Temp,
    Psi,
    Tolerance
) {
  # Initialize names and indicies
  n_taxa <- params["n_taxa", ]
  n_substrates <- params["n_substrates", ]
  n_enzymes <- params["n_enzymes", ]
  Mon.names <- c("NH4", "PO4", rownames(Substrates[1:n_substrates, ]))
  n_monomers <- length(Mon.names)
  Mon.names.G <- rownames(Monomers)
  Enz.names <- colnames(Vmax0)
  Sub.names <- rownames(Substrates[1:n_substrates, ])
  Mic.names <- rownames(Microbes[1:n_taxa, ])
  is.NH4 <- which(rownames(Monomers) == "NH4")
  is.PO4 <- which(rownames(Monomers) == "PO4")
  is.Hemi <- which(rownames(Monomers) == "Mon004")
  org <- which(rownames(Monomers) != c("NH4", "PO4"))
  mineral <- which(rownames(Monomers) == c("NH4", "PO4"))
  is.deadMic <- which(rownames(Substrates) == "DeadMic")
  is.deadEnz <- which(rownames(Substrates) == "DeadEnz")
  is.cellulose <- which(rownames(Substrates) == "Cellulose")
  is.lignin <- which(rownames(Substrates) == "Lignin")
  MonomerRatios <- Monomers
  MonomerRatios[, ] <- 0
  MonomerRatios[is.NH4, "N"] <- 1
  MonomerRatios[is.PO4, "P"] <- 1
  Tref <- 293
  
  # Create matrices to be filled in
  tev <- matrix(rep(0, n_substrates*grid.size*n_enzymes), 
                ncol = n_enzymes)
  Max_Uptake1 <- matrix(rep(0, grid.size*n_taxa*n_monomers),
                        nrow = n_taxa)
  Max_Uptake2 <- matrix(Max_Uptake1,
                        ncol = n_taxa,
                        dimnames = list(Mon.names.G, Mic.names))
  TUC.mat <- TUN.mat <- TUP.mat <- matrix(rep(0, grid.size*n_taxa*n_monomers), 
                                          nrow = n_monomers)
  DOC.mat <- DON.mat <- DOP.mat <- matrix(rep(0, grid.size*n_monomers), nrow = n_monomers)
  Microbe.C.cell <- Microbe.N.cell <- Microbe.P.cell <- matrix(rep(0, grid.size*n_taxa), nrow = n_taxa)
  EP.mat <- matrix(rep(0, grid.size*n_taxa*n_enzymes),
                   nrow = n_taxa,
                   dimnames = list(Mic.names, rownames(Enzymes)))
  Death.mat <- matrix(rep(0, grid.size*n_taxa*3), 
                      nrow = n_taxa)
  DeadEnz.mat <- matrix(rep(0, grid.size*n_enzymes*3), 
                        nrow = n_enzymes)
  Resp.comp <- matrix(rep(0, 3), ncol = 3, 
                      dimnames = list(NULL, c("Maint", "Growth", "Overflow")))
  Mic.growth <- matrix(rep(0, 3), ncol = 3, 
                       dimnames = list(NULL, c("C", "N", "P")))
  
  # Create index for matrix manipulation of taxon uptake
  # index1 is the shape of the starting matrix
  TU.index1 <- matrix(1:(grid.size*n_taxa*n_monomers),
                      ncol = n_taxa)
  TU.index <- matrix(TU.index1, 
                     nrow = n_monomers)
  for(i in 1:grid.size) {
    rows <- ((i-1)*n_monomers+1):(i*n_monomers)
    cols <- ((i-1)*n_taxa+1):(i*n_taxa)
    TU.index[, cols] <- TU.index1[rows, ]
  }
  
  # Create index for matrix manipulation of taxon enzyme production
  EP.index1 <- matrix(1:(grid.size*n_taxa*n_enzymes),
                      ncol = n_enzymes)
  EP.index <- matrix(EP.index1,
                     nrow = n_taxa)
  for(i in 1:grid.size) {
    rows <- ((i-1)*n_taxa+1):(i*n_taxa)
    cols <- ((i-1)*n_enzymes+1):(i*n_enzymes)
    EP.index[, cols] <- EP.index1[rows, ]
  }
  
  # Index to manipulate Vmax values
  vm.index1 <- matrix(1:(grid.size*n_substrates*n_enzymes), 
                      ncol = n_substrates)
  vm.index <- matrix(vm.index1,
                     ncol = n_enzymes)
  for(i in 1:grid.size) {
    vm.index[((i-1)*n_substrates+1):(i*n_substrates),] <- t(vm.index1[((i-1)*n_enzymes+1):(i*n_enzymes),])
  }
  
  # Set grid sums to initial values
  Monomers.grid <- sum.grid(Monomers, Mon.names, grid.size)
  Enzymes.grid <- sum.grid(Enzymes, Enz.names, grid.size)
  Substrates.grid <- sum.grid(Substrates, Sub.names, grid.size)
  Microbes.grid <- sum.grid(Microbes, Mic.names, grid.size)
  n_fungi <- sum(fb)/grid.size
  Fungi <- which(fb[1:n_taxa] == 1)
  is.fungi <- fb == 1
  Fungi.count.zero <- Fungi.count <- rep(0, n_taxa)
  Fungi.vec <- Fungi.vec.zero <- matrix(rep(0, n_taxa*grid.size), 
                                        ncol = 1,
                                        dimnames = list(rownames(Microbes),"Count"))
  Fungi.count[Fungi] <- Microbes.grid[Fungi, "C"]/(0.5*params["max_size_f", ])
  z <- rep(1:n_taxa, grid.size)
  kill.reset <- rep(0, grid.size*n_taxa)
  
  # Initialize the time series of data to hold
  RespSeries <- 0
  Resp.compSeries <- 0
  C.UptakeSeries <- 0
  N.UptakeSeries <- 0
  P.UptakeSeries <- 0
  N.MonomerUptakeSeries <- 0
  Mic.growthSeries <- 0
  Net.CUE_Series <- 0
  N.mineralSeries <- 0
  P.mineralSeries <- 0
  Cum_Leaching_N <- 0
  Cum_Leaching_P <- 0
  EnzymesSeries <- t(Enzymes.grid[, "C"])
  EnzProdSeries <- as.data.frame(matrix(0, ncol = n_enzymes, dimnames = list(NULL, Enz.names)))
  SubstratesSeries <- Substrates.grid[, "C"]
  Cum_SubstrateSeries <- colSums(Substrates.grid)
  MonProdSeries <- as.data.frame(matrix(0, ncol = n_substrates, dimnames = list(NULL, Sub.names)))
  C.MonomersSeries <- Monomers.grid[, "C"]
  N.MonomersSeries <- Monomers.grid[, "N"]
  NH4Series <- Monomers.grid["NH4", "N"]
  PO4Series <- Monomers.grid["PO4", "P"]
  MicrobesSeries <- Microbes.grid[, "C"]
  MicrobesNSeries <- Microbes.grid[, "C"]/Microbes.grid[, "N"]
  MicrobesPSeries <- Microbes.grid[, "C"]/Microbes.grid[, "P"]
  MicrobesNPSeries <- Microbes.grid[, "N"]/Microbes.grid[, "P"]
  
  # Sum microbe, substrate, monomer, enzyme pools
  Mic_Sum <- Cum_Microbe <- colSums(Microbes.grid)
  Substrate_Sum <- sum(Substrates.grid[, "C"])
  Cum_Substrate <- colSums(Substrates.grid)
  Monomer_Sum <- sum(Monomers.grid[, "C"])
  Cum_Monomer <- colSums(Monomers.grid)
  Enzyme_Sum <- Cum_Enzyme_C <- sum(Enzymes.grid)
  Cum_Enzyme_N <- sum(Enzymes.grid*EnzAttrib[, "N_cost"])
  Cum_Enzyme_P <- sum(Enzymes.grid*EnzAttrib[, "P_cost"])
  
  # Relative enzyme carbon cost for each enzyme gene
  Induce_Enzyme_C <- t(EnzProdInduce)*EnzAttrib[, "C_cost"]
  Constit_Enzyme_C <- t(EnzProdConstit)*EnzAttrib[, "C_cost"]
  # Relative uptake enzyme carbon cost for each enzyme gene
  Uptake_Cost <- t(UptakeGenes)
  rownames(Uptake_Cost) <- sprintf("%s%03d", "Upt", 1:params["n_uptake",]) # Original was 1:params["n_uptake",]
  
  # Expand variables to the size of the grid
  Induce_Enzyme_C <- expand(t(Induce_Enzyme_C), grid.size)
  Constit_Enzyme_C <- expand(t(Constit_Enzyme_C), grid.size)
  Uptake_Cost <- expand(t(Uptake_Cost), grid.size)
  tUptakeGenes <- expand(t(UptakeGenes), grid.size)
  UptakeGenes <- expand(UptakeGenes, grid.size)
  
  # Create variable to hold grid at maximum enzyme biomass
  grid.max <- NULL
  
  for(i_t in 1:params["end_time", ]) { # Loop through time
    
    
    # Mortality rate increases linearly with -Psi times slope beta
    # which can differ for bacteria vs fungi; tolerance of 1 eliminates
    # dependence on Psi. Fungi and bacteria can also differ in baseline
    # mortality according to Death_Ratio
    r.death <- (1-fb)*params["Death_Rate", ]*(1 - (params["beta.bac", ]*Psi[i_t])*(1 - Tolerance)) +
      fb*params["Death_Rate", ]*((1 - (params["beta.fungi", ]*Psi[i_t])*(1 - Tolerance))*params["Death_Ratio", ])
    
    # Calculate Vmax values with Arrhenius equation
    # Gas constant = 0.008314 kJ/(mol K)
    # Arrhenius equation for Vmax multiplied by exponential decay for Psi sensitivity
    Vmax <- Vmax0*exp((-Ea/0.008314)*(1/(Temp[i_t] + 273) - 1/Tref))*exp(params["Psi.slope.Vmax", ]*Psi[i_t])
    Uptake_Vmax <- Uptake_Vmax0*exp((-Uptake_Ea/0.008314)*(1/(Temp[i_t] + 273) - 1/Tref))*exp(params["Psi.slope.uptake", ]*Psi[i_t])
    # Expand to the size of the grid after calculating
    etVmax <- expand(t(Vmax), grid.size)
    Uptake_Vmax <- expand(Uptake_Vmax, grid.size)
    
    # Recalculate Km values with Arrhenius equation
    Km <- Km0*exp((-params["Km_Ea",]/0.008314)*(1/(Temp[i_t] + 273) - 1/Tref))
    Uptake_Km <- Uptake_Km0*exp((-params["Km_Ea",]/0.008314)*(1/(Temp[i_t] + 273) - 1/Tref))
    
    # Carbon use efficiency dependent on temperature and number of enzyme genes
    CUE <- CUE.ref + (Temp[i_t] - (Tref-273))*params["CUE_temp", ]
    
    # Reset the reproduction matrix
    Reprod <- ReprodNew
    # Reset the colonization matrix
    Colonization <- Colonization.reset
    # Fungal translocation: calculate average biomass within fungal taxa
    Mean.fungi <- Microbes.grid/as.vector(Fungi.count)
    Mean.fungi[!is.finite(Mean.fungi)] <- 0
    # Expand the fungal average across the grid
    eMF <- expand(Mean.fungi, grid.size)
    # Reset the fungal count to zero
    Fungi.count <- Fungi.count.zero
    
    # Select the daughter cells that are fungi versus bacteria
    daughters.b <- which(Reprod[, "C"] > 0 & fb == 0)
    daughters.f <- which(Reprod[, "C"] > 0 & fb == 1)
    num.b <- length(daughters.b)
    num.f <- length(daughters.f)
    shift_x <- shift_y <- rep(0, grid.size*n_taxa)
    shift_x[daughters.b] <- sample(c(-params["dist", ]:params["dist", ]), 
                                   num.b,replace=T) # vector of dispersal movements in x direction
    shift_y[daughters.b] <- sample(c(-params["dist", ]:params["dist", ]), 
                                   num.b,replace=T) # vector of dispersal movements in y direction
    shift_x[daughters.f] <- 1 # Fungi always move positively in x direction
    # vector of dispersal movements in y direction; constrained to one box away determined by probability "direct"
    shift_y[daughters.f] <- sample(c(-1:1),
                                   num.f, 
                                   replace = T,
                                   prob = c(0.5*(1 - params["direct", ]),
                                            params["direct", ],
                                            0.5*(1 - params["direct",])))
    new_x <- (rep(1:params["x", ],
                  each = n_taxa,
                  times = params["y", ]) + shift_x + params["x", ])%%params["x", ] # calculate x coordinates of dispersal destinations
    new_x[new_x == 0] <- params["x", ] # Substitute coordinates when there is no shift
    new_y <- (rep(1:params["y", ],each = params["x", ]*n_taxa) + shift_y + params["y", ])%%params["y", ] # calculate y coordinates of dispersal destinations
    new_y[new_y == 0] <- params["y", ]
    Reprod[daughters.f, ] <- eMF[daughters.f, ] # set all fungi equal to their grid averages for translocation before colonization
    index_col <- n_taxa*((new_y-1)*params["x", ] + (new_x-1)) + z # convert x,y coordinates to a vector of destination locations
    # Transfer cells to new locations and sum when two or more of the same taxa go to same location
    for(i in c(daughters.b,daughters.f)) {
      Colonization[index_col[i], ] <- Colonization[index_col[i], ] + Reprod[i, ]
    }
    
    # Set Monomers to the grid average
    Monomers <- expand(Monomers.grid/grid.size, grid.size)
    
    # Translocate nutrients within fungal taxa
    i <- Microbes[fb == 1, "C"] > 0
    Microbes[fb == 1, ][i, ] <- eMF[fb == 1, ][i, ]
    
    # Colonization of dispersing microbes
    Microbes <- Microbes + Colonization
    
    # Multiply Vmax values for each substrate by the quantity of each enzyme
    tev[, ] <- (as.vector(Enzymes)*etVmax)[vm.index]
    
    # Equation for Michaelis-Menten enzyme catalysis
    rss <- rowSums(Substrates)
    Decay <- tev*rss/(Km + rss)
    
    # Loop that pulls out each batch of required enzymes and sums across redundant enzymes
    DecaySums <- NULL
    for(i in 1:dim(ReqEnz)[3]) {
      DecaySums <- cbind(rowSums(as.matrix(ReqEnz[,,i])*Decay), DecaySums)
    }
    
    # Assess the rate-limiting enzyme and set decay to that rate
    # Compare to substrate available and take the min, allowing for a tolerance of 1e-9
    # Link cellulose degradation to lignin concentration (LCI)
    DecayRates <- pmin(pmin(DecaySums[, 1],
                            DecaySums[, 2],
                            na.rm = T),
                       rss - 1e-9*rss,
                       na.rm = T)
    ss7 <- rowSums(Substrates[is.lignin, ])
    DecayRates[is.cellulose] <- DecayRates[is.cellulose]*(1 + params["LCI_slope", ]*ss7/(Substrates[is.cellulose,1] + ss7))
    # Replacing NAs with 0 to avoid premature end of the simulation
    DecayRates[is.na(DecayRates)] <- 0
    
    
    # Monomer stoichiometry: organic monomers follow substrate
    SubstrateRatios <- Substrates/rss
    SubstrateRatios[rss == 0, ] <- 0
    MonomerRatios[org, ] <- SubstrateRatios
    
    # Update monomer pools
    # Decayed substrate stays in the grid box for the whole iteration
    # Preferential access by resident microbes
    # Any leftovers diffuse (get averaged) across grid at next iteration
    Monomers[org, ] <- Monomers[org, ] + (DecayRates+MonInput[org])*MonomerRatios[org, ]
    Monomers[mineral, ] <- Monomers[mineral, ] + MonInput[mineral]*MonomerRatios[mineral, ]
    # Keep track of mass balance for inputs
    Cum_Monomer <- Cum_Monomer + colSums(as.vector(MonInput)*MonomerRatios)
    
    # Monomer production across grid
    MonProd.grid <- colSums(matrix(DecayRates, ncol = n_substrates, byrow = TRUE, dimnames = list(NULL, Sub.names)))
    
    
    # Recalculate Monomer stoichiometry after changes due to dead microbial biomass
    rsm <- rowSums(Monomers)
    MonomerRatios[org, ] <- Monomers[org, ]/rsm[org]
    MonomerRatios[org, ][rsm[org] == 0, ] <- 0
    
    
    # Section for modulating the uptake allocation based on uptake of the previous
    # time step. The higher the uptake for a given element, the lower the 
    # uptake allocation
    
    Upt.repres <- 1
    
    if(params["Upt_repress",] == 1){
      
      if(i_t > 1){
        # Calculating per capita nutrient uptake per taxa per grid box
        C.up.pc <- (Taxon_Uptake_C/(Microbes[,"C"] + 10^-10))/(params["C_min",])
        N.up.pc <- (Taxon_Uptake_N/(Microbes[,"C"] + 10^-10))/(params["N_min",])
        P.up.pc <- (Taxon_Uptake_P/(Microbes[,"C"] + 10^-10))/(params["P_min",])
        NP.up.pc <- apply(cbind(N.up.pc, P.up.pc), 1, min) # Selecting the lowest between N and P
        CN.up.pc <- apply(cbind(C.up.pc, N.up.pc), 1, min)
        CP.up.pc <- apply(cbind(C.up.pc, P.up.pc), 1, min) 
        CNP.up.pc <- apply(cbind(C.up.pc, N.up.pc, P.up.pc), 1, min)
        
        
        Upt.repres.l <- list(
          C = round((exp(-12*(C.up.pc - 0.75)) + 0.1) / (1 + exp(-12*(C.up.pc - 0.75) )), digits = 4),
          N = round((exp(-12*(N.up.pc - 0.75)) + 0.1) / (1 + exp(-12*(N.up.pc - 0.75) )), digits = 4),
          P = round((exp(-12*(P.up.pc - 0.75)) + 0.1) / (1 + exp(-12*(P.up.pc - 0.75) )), digits = 4),
          NP = round((exp(-12*(NP.up.pc - 0.75)) + 0.1) / (1 + exp(-12*(NP.up.pc - 0.75) )), digits = 4),  # Old rep(1, length(Taxon_Uptake_N))
          CN = round((exp(-12*(CN.up.pc - 0.75)) + 0.1) / (1 + exp(-12*(CN.up.pc - 0.75) )), digits = 4),
          CP = round((exp(-12*(CP.up.pc - 0.75)) + 0.1) / (1 + exp(-12*(CP.up.pc - 0.75) )), digits = 4),
          CNP = round((exp(-12*(CNP.up.pc - 0.75)) + 0.1) / (1 + exp(-12*(CNP.up.pc - 0.75) )), digits = 4)
          )
        
        
        Upt.repres <- do.call("cbind", Upt.repres.l[t(Upt.track.look)])
      }
    }
    
    
    # Multiply microbial biomass by each taxon's uptake allocation to get biomass of each uptake enzyme by taxon
    rsi <- rowSums(Microbes)
    MicCXGenes <- rsi*UptakeGenes*Upt.repres
    
    # Equation for hypothetical potential uptake (per unit of compatible uptake protein)
    Potential_Uptake <- Uptake_ReqEnz*Uptake_Vmax*rsm/(Uptake_Km + rsm)
    
    # Matrix multiplication to get max possible uptake by taxon and monomer
    # Must extract each grid point separately for operation
    for(i in 1:grid.size) {
      i.a <- ((i - 1)*n_monomers + 1):(i*n_monomers)
      i.b <- ((i - 1)*n_taxa + 1):(i*n_taxa)
      Max_Uptake <- MicCXGenes[i.b, , drop = F]%*%t(Potential_Uptake[i.a, , drop = F]) # Subsetting Potential uptake by rows in chunks of 14 (grid cells)
      Max_Uptake1[, i.a] <- Max_Uptake # wide format; mon*grid cols; rows are taxa
      Max_Uptake2[i.a, ] <- t(Max_Uptake) # long format; mon*grid rows; cols are taxa
    }
    # Sum the total potential uptake of each monomer
    csmu <- colSums(Max_Uptake1)
    
    # Take the min of the monomer available and the max potential uptake
    Min_Uptake <- pmin(csmu, rsm, na.rm = T)
    
    # Scale the uptake to what's available
    Uptake <- Max_Uptake2*Min_Uptake/csmu
    Uptake[csmu == 0, ] <- 0
    # Prevent total uptake from getting too close to zero
    Uptake <- Uptake - 1e-9*Uptake
    
    # Calculate total uptake by monomer
    Monomer_Uptake <- rowSums(Uptake)*MonomerRatios
    
    
    # Sum across cell grids for total monomer uptake
    Monomer_Uptake.grid <- sum.grid(Monomer_Uptake, 
                                    c("NH4", "PO4", Sub.names), grid.size)
    
    # Calculate total uptake by taxon
    TUC.mat[,] <- (Uptake*MonomerRatios[, "C"])[TU.index]
    TUN.mat[,] <- (Uptake*MonomerRatios[, "N"])[TU.index]
    TUP.mat[,] <- (Uptake*MonomerRatios[, "P"])[TU.index]
    # Sum across monomers
    Taxon_Uptake_C <- colSums(TUC.mat)
    Taxon_Uptake_N <- colSums(TUN.mat)
    Taxon_Uptake_P <- colSums(TUP.mat)
    # Sum across cell grids for total taxon uptake
    Taxon_Uptake.grid <- sum.grid(cbind(Taxon_Uptake_C, Taxon_Uptake_N, Taxon_Uptake_P), 
                                  Mic.names, grid.size)
    
    
    # Section to adjust enzyme production based on uptake. The higher the uptake, 
    # the lower the enzyme production. This is calculated independently for each enzyme
    
    
    # Calculate sigmoid curves for repression of enzyme production
    Enz.repres <- 1
    
    if(params["Enz_repress",] == 1){
      
      
      # Uptake ratios
      Upt.rat <- cbind(Taxon_Uptake_C, Taxon_Uptake_N, Taxon_Uptake_P)/
                 rowSums(cbind(Taxon_Uptake_C, Taxon_Uptake_N, Taxon_Uptake_P))
      
      # Calculating uptake excess
      Upt.ex <- Upt.rat - MinRatios
      
      Upt.lim <- c("C", "N", "P")[as.integer(apply(Upt.ex, 1, which.min))]
      
      
      # Calculating per capita nutrient uptake per taxa per grid box
      C.up.pc <- (Taxon_Uptake_C/(Microbes[,"C"] + 10^-10))/(params["C_min",])
      N.up.pc <- (Taxon_Uptake_N/(Microbes[,"C"] + 10^-10))/(params["N_min",])
      P.up.pc <- (Taxon_Uptake_P/(Microbes[,"C"] + 10^-10))/(params["P_min",])
      NP.up.pc <- apply(cbind(N.up.pc, P.up.pc), 1, min) # Selecting the lowest between N, and P
      CN.up.pc <- apply(cbind(C.up.pc, N.up.pc), 1, min)
      CP.up.pc <- apply(cbind(C.up.pc, P.up.pc), 1, min) 
      CNP.up.pc <- apply(cbind(C.up.pc, N.up.pc, P.up.pc), 1, min)
      
      Enz.repres.l <- list(
        C = round((exp(-12*(C.up.pc - 0.75)) + 0.01) / (1 + exp(-12*(C.up.pc - 0.75) )), digits = 4),
        N = round((exp(-12*(N.up.pc - 0.75)) + 0.01) / (1 + exp(-12*(N.up.pc - 0.75) )), digits = 4),
        P = round((exp(-12*(P.up.pc - 0.75)) + 0.01) / (1 + exp(-12*(P.up.pc - 0.75) )), digits = 4),
        CN = round((exp(-12*(CN.up.pc - 0.75)) + 0.01) / (1 + exp(-12*(CN.up.pc - 0.75) )), digits = 4),
        CP = round((exp(-12*(CP.up.pc - 0.75)) + 0.01) / (1 + exp(-12*(CP.up.pc - 0.75) )), digits = 4),
        NP = round((exp(-12*(NP.up.pc - 0.75)) + 0.01) / (1 + exp(-12*(NP.up.pc - 0.75) )), digits = 4),  # Old rep(1, length(Taxon_Uptake_N))
        CNP = round((exp(-12*(CNP.up.pc - 0.75)) + 0.01) / (1 + exp(-12*(CNP.up.pc - 0.75) )), digits = 4)
      )
      Enz.repres <- do.call("cbind", Enz.repres.l[Enz.track.look[colnames(Induce_Enzyme_C)]])
      
      
    }
    
    
    
    
    # Enzyme production: rows are taxa, cols are enzyme genes
    # Two alternatives: proportional to biomass or proportional to uptake C
    Taxon_Enzyme_Production <- Taxon_Uptake_C*Induce_Enzyme_C*Enz.repres + rsi*Constit_Enzyme_C*Enz.repres
    
    # Constrain enzyme production if it will cost more N than currently available in microbial biomass
    Enzyme_Cost_N <- colSums(t(Taxon_Enzyme_Production)*EnzAttrib[, "N_cost"])
    # Only select entries with microbes present
    i <- which(Enzyme_Cost_N > 0)
    Taxon_Enzyme_Production[i, ] <- Taxon_Enzyme_Production[i, ]*pmin(Microbes[i, "N"], 
                                                                      Enzyme_Cost_N[i])/Enzyme_Cost_N[i]
    
    # Total enzyme production across all taxa
    EP.mat[, ] <- Taxon_Enzyme_Production[EP.index]
    Enzyme_Production <- colSums(EP.mat)
    
    
    # Total enzyme carbon cost for each taxon
    ttep <- t(Taxon_Enzyme_Production)
    Enzyme_Maint <- colSums(ttep*EnzAttrib[, "Maint_cost"])
    Enzyme_Cost <- colSums(ttep) + Enzyme_Maint
    Enzyme_Cost_N <- colSums(ttep*EnzAttrib[, "N_cost"])
    Enzyme_Cost_P <- colSums(ttep*EnzAttrib[, "P_cost"])
    
    # Uptake enzyme production: cols are taxa, rows are enzymes
    # Total uptake enzyme carbon cost for each taxon
    Uptake_Maint <- colSums(t(rsi*Uptake_Cost*Upt.repres))*params["Uptake_Maint_cost", ]
    
    # Calculate loss rates and random death
    Enzyme_Loss <- params["Enzyme_Loss_Rate", ]*Enzymes
    
    # Update enzyme pools
    Enzymes <- Enzymes + Enzyme_Production - Enzyme_Loss
    
    # Update microbe pools
    Microbes[, "C"] <- Microbes[, "C"] + Taxon_Uptake_C*CUE - Enzyme_Cost - Enzyme_Maint - Uptake_Maint
    Microbes[, "N"] <- Microbes[, "N"] + Taxon_Uptake_N - Enzyme_Cost_N
    Microbes[, "P"] <- Microbes[, "P"] + Taxon_Uptake_P - Enzyme_Cost_P
    
    Mic.growth[, "C"] <- sum(Taxon_Uptake_C*CUE - Enzyme_Cost - Enzyme_Maint - Uptake_Maint)
    Mic.growth[, "N"] <- sum(Taxon_Uptake_N - Enzyme_Cost_N)
    Mic.growth[, "P"] <- sum(Taxon_Uptake_P - Enzyme_Cost_P)
    
    # Kill microbes that are starving and transfer biomass to substrate
    Death <- Colonization.reset
    if (is.na(sum(Microbes))) {break}
    starve_index <- (Microbes[, "C"] < params["C_min", ] | Microbes[, "N"] < params["N_min", ] | Microbes[,"P"] < params["P_min", ]) & Microbes[,"C"] > 0
    Death[starve_index, ] <- Microbes[starve_index, ]
    Microbes[starve_index, ] <- 0
    
    # Determine which locations have microbes
    is.m <- which(Microbes[, "C"] > 0)
    # Zero the death index
    kill <- kill.reset
    # Sample for death events (=1)
    kill[is.m] <- runif(length(is.m)) < r.death[is.m]
    # Kill microbes
    killed <- kill*Microbes
    Death <- Death + killed
    Microbes <- Microbes - killed
    
    # Adjust stoichiometry of microbial biomass
    # Index locations of microbial cells
    mic.index <- as.numeric(which(Microbes[, "C"] > 0))
    # drop=F retains matrix structure if only one row is selected
    Mic.subset <- Microbes[mic.index,,drop = F]
    # Calculate microbial ratios
    MicrobeRatios <- (Mic.subset/rowSums(Mic.subset))
    # Select the corresponding minimum quotas
    MinRat <- MinRatios[mic.index,,drop = F]
    # Index only microbes that have below-minimum quotas
    rat.index <- mic.index[MicrobeRatios[, "C"] < MinRat[, "C"] | MicrobeRatios[, "N"] < MinRat[, "N"] | MicrobeRatios[, "P"] < MinRat[, "P"]]
    Mic.subset <- Microbes[rat.index,,drop = F]
    StartMicrobes <- Mic.subset
    MicrobeRatios <- Excess <- (Mic.subset/rowSums(Mic.subset))
    MinRat <- MinRatios[rat.index,,drop = F]
    
    # Calculate difference between min and actual ratios    
    Deficit <- Deficit.0 <- MicrobeRatios - MinRat
    # Determine the limiting nutrient that will be conserved
    Limiting <- max.col(-Deficit/MinRat,
                        ties.method = "first")
    # Set all deficient ratios to their minima
    MicrobeRatios[Deficit < 0] <- MinRat[Deficit < 0]
    # Reduce the mass fractions for non-deficient elements in proportion to the distance from the minimum
    # Calculate how far above the minimum each non-deficient element is
    Excess[Deficit > 0] <- Deficit[Deficit > 0]
    # Set deficient element fractions to zero
    Excess[Deficit < 0] <- 0
    # Partition the total deficit to the excess element(s) in proportion to their distances from their minima
    Deficit.0[Deficit > 0] <- 0
    MicrobeRatios[Deficit > 0] <- MicrobeRatios[Deficit > 0] + (rowSums(Deficit.0)*Excess/rowSums(Excess))[Deficit> 0]
    
    # Construct hypothetical nutrient quotas for each possible minimum nutrient
    MC <- Mic.subset[, "C"]; MN <- Mic.subset[, "N"]; MP <- Mic.subset[, "P"]
    MRC <- MicrobeRatios[, "C"]; MRN <- MicrobeRatios[, "N"]; MRP <- MicrobeRatios[, "P"]
    new.C <- cbind(MC, MN*MRC/MRN, MP*MRC/MRP)
    new.N <- cbind(MC*MRN/MRC, MN, MP*MRN/MRP)
    new.P <- cbind(MC*MRP/MRC, MN*MRP/MRN, MP)
    # Insert the appropriate set of nutrient quotas scaled to the minimum nutrient
    select <- cbind(1:length(rat.index),Limiting)
    Microbes[rat.index,] <- cbind(new.C[select],
                                  new.N[select],
                                  new.P[select])
    # Sum up the element losses from biomass across whole grid and calculate average loss
    MicLoss <- StartMicrobes - Microbes[rat.index,]
    
    # Update substrate pools
    # Substrate stoichiometry
    rss <- rowSums(Substrates)
    SubstrateRatios <- Substrates/rss
    SubstrateRatios[rss==0,] <- 0
    # Add inputs and remove decay
    Substrates <- Substrates + SubInput - DecayRates*SubstrateRatios
    # Account for inputs in mass balance
    Cum_Substrate <- Cum_Substrate + colSums(SubInput)
    # Add dead microbial biomass
    Death.mat[,] <- Death
    Substrates[is.deadMic,] <- Substrates[is.deadMic,] + colSums(Death.mat)
    # Add dead enzymes
    DeadEnz.mat[,] <- c(Enzyme_Loss,
                        Enzyme_Loss*EnzAttrib[,"N_cost"],
                        Enzyme_Loss*EnzAttrib[,"P_cost"])
    Substrates[is.deadEnz,] <- Substrates[is.deadEnz,] + colSums(DeadEnz.mat)
    
    # Update monomer pools
    Monomers[is.NH4,"N"] <- Monomers[is.NH4,"N"] + sum(MicLoss[,"N"])/grid.size
    Monomers[is.PO4,"P"] <- Monomers[is.PO4,"P"] + sum(MicLoss[,"P"])/grid.size
    Monomers1 <- Monomers # Snapshot of monomers before uptake
    Monomers <- Monomers - Monomer_Uptake
    
    # C Overflow. Goes to respiration or dissolved organic matter. For the latter
    # case, to avoid introducing a new substrate, we can put C overflow in xylose 
    # (hemicellulose monomer)
    if(params["DOC.overflow",] == 1){
    
       Monomers[is.Hemi, "C"] <- Monomers[is.Hemi, "C"] + sum(MicLoss[,"C"])/grid.size
    
    }else{Resp.comp[, "Overflow"] <- sum(MicLoss[,"C"])}
    
    
    # Rates of mineralization. Respiration is separated by components
    Resp.comp[, "Maint"] <- sum(Enzyme_Maint) + sum(Uptake_Maint)
    Resp.comp[, "Growth"] <- sum(Taxon_Uptake_C*(1-CUE))
    N.mineralization <- sum(MicLoss[,"N"])
    P.mineralization <- sum(MicLoss[,"P"])
    
    # Sum microbes prior to reproduction
    Microbes.grid <- sum.grid(Microbes,Mic.names,grid.size)
    # Reset the vector of fungal locations
    Fungi.vec <- Fungi.vec.zero
    # Add one or two fungi to the count vector based on size
    Fungi.vec[fb==1][Microbes[fb==1,"C"]>0] <- 1
    Fungi.vec[fb==1][Microbes[fb==1,"C"]>params["max_size_f",]] <- 2
    # Sum up the count vector
    Fungi.count <- sum.grid(Fungi.vec,Mic.names,grid.size)
    
    # Cell division
    MicrobesBeforeDivision <- Microbes
    Microbes[fb==0,][Microbes[fb==0,"C"]>params["max_size_b",],] <- Microbes[fb==0,][Microbes[fb==0,"C"]>params["max_size_b",],]/2
    Microbes[fb==1,][Microbes[fb==1,"C"]>params["max_size_f",],] <- Microbes[fb==1,][Microbes[fb==1,"C"]>params["max_size_f",],]/2
    
    # Add daughter cells to matrix of new reproduction
    ReprodNew <- MicrobesBeforeDivision-Microbes
    
    # Respiration components per taxon
    Taxon_Resp.grid <- sum.grid(cbind(Base_resp = Taxon_Uptake_C*(1-(CUE)), Enzyme_Maint, Uptake_Maint, Resp = NA), 
                                Mic.names, grid.size)
    Taxon_Resp.grid[, "Resp"] <-  rowSums(Taxon_Resp.grid[, 1:3, drop = FALSE])
    
    # Exclude taxa which uptake is less than respiration
    exclude <- Taxon_Uptake.grid[, "Taxon_Uptake_C"] > Taxon_Resp.grid[,"Resp"] 
    
    # Calculating net growth and CUE
    Net.growth <- sum(Taxon_Uptake.grid[exclude, "Taxon_Uptake_C"]) -
      sum(Taxon_Resp.grid[exclude, "Resp"])
    Net.CUE <- Net.growth/sum(Taxon_Uptake.grid[exclude, "Taxon_Uptake_C"])
    
    # Sum pools across grid
    Respiration <- sum(Taxon_Uptake_C*(1-CUE)) + sum(Enzyme_Maint) + sum(Uptake_Maint) + Resp.comp[, "Overflow"]
    Monomers.grid <- sum.grid(Monomers, Mon.names, grid.size)
    Monomers.grid1 <- sum.grid(Monomers1, Mon.names, grid.size)
    Enzymes.grid <- sum.grid(Enzymes, Enz.names, grid.size)
    Substrates.grid <- sum.grid(Substrates, Sub.names, grid.size)
    
    Enzyme_Production.grid <- colSums(matrix(Enzyme_Production, ncol = n_enzymes, byrow = TRUE, dimnames = list(NULL, Enz.names)))
    
    # Update monomer leaching
    Leaching_N <- Monomers.grid["NH4","N"]*params["Leaching",]*exp(params["Psi.slope.leach",]*Psi[i_t])
    Leaching_P <- Monomers.grid["PO4","P"]*params["Leaching",]*exp(params["Psi.slope.leach",]*Psi[i_t])
    # Keep track of mass balance
    Cum_Leaching_N <- Cum_Leaching_N + Leaching_N
    Cum_Leaching_P <- Cum_Leaching_P + Leaching_P
    Monomers.grid["NH4","N"] <- Monomers.grid["NH4","N"] - Leaching_N
    Monomers.grid["PO4","P"] <- Monomers.grid["PO4","P"] - Leaching_P
    
    # Record total pools for each time step
    RespSeries <- c(RespSeries,Respiration)
    Resp.compSeries <- rbind(Resp.compSeries, Resp.comp)
    C.UptakeSeries <- rbind(C.UptakeSeries, Taxon_Uptake.grid[, "Taxon_Uptake_C"])
    N.UptakeSeries <- rbind(N.UptakeSeries, Taxon_Uptake.grid[, "Taxon_Uptake_N"])
    P.UptakeSeries <- rbind(P.UptakeSeries, Taxon_Uptake.grid[, "Taxon_Uptake_P"])
    N.MonomerUptakeSeries <- rbind(N.MonomerUptakeSeries, Monomer_Uptake.grid[, "N"])
    Mic.growthSeries <- rbind(Mic.growthSeries, Mic.growth)
    Net.CUE_Series <- c(Net.CUE_Series, Net.CUE)
    MonProdSeries <- rbind(MonProdSeries, MonProd.grid)
    C.MonomersSeries <- rbind(C.MonomersSeries, Monomers.grid1[,"C"])
    N.MonomersSeries <- rbind(N.MonomersSeries, Monomers.grid1[,"N"])
    N.mineralSeries <- c(N.mineralSeries, N.mineralization)
    P.mineralSeries <- c(P.mineralSeries, P.mineralization)
    NH4Series <- rbind(NH4Series, Monomers.grid1["NH4","N"])
    PO4Series <- rbind(PO4Series, Monomers.grid1["PO4","P"])
    EnzymesSeries <- rbind(EnzymesSeries,t(Enzymes.grid))
    EnzProdSeries <- rbind(EnzProdSeries, Enzyme_Production.grid)
    SubstratesSeries <- rbind(SubstratesSeries,Substrates.grid[,"C"])
    Cum_SubstrateSeries <- rbind(Cum_SubstrateSeries, colSums(Substrates.grid))
    MicrobesSeries <- rbind(MicrobesSeries,Microbes.grid[,"C"])
    MicrobesNSeries <- rbind(MicrobesNSeries,Microbes.grid[,"C"]/Microbes.grid[,"N"])
    MicrobesPSeries <- rbind(MicrobesPSeries,Microbes.grid[,"C"]/Microbes.grid[,"P"])
    MicrobesNPSeries <- rbind(MicrobesNPSeries,Microbes.grid[,"N"]/Microbes.grid[,"P"])
    
    # Sum microbial pools
    Mic_Sum <- rbind(Mic_Sum, colSums(Microbes.grid))
    
    # Output grid if enzyme mass has peaked
    if (i_t>3 & params["output.mid",] == 1) {
      if (Enzyme_Sum[i_t] < Enzyme_Sum[i_t - 1] & Enzyme_Sum[i_t-1] > Enzyme_Sum[i_t-2]) {
        grid.max <- list("Substrates" = Substrates,
                         "Monomers" = Monomers,
                         "Enzymes" = Enzymes,
                         "Microbes" = Microbes)
        grid.time <- i_t
      }
    }
    
    # Output grid images
    if (i_t == 1 || i_t%%params["print.grid",]==0){
      apply(matrix(1:n_taxa),1,function(i)plot.grid(Microbes[n_taxa*(0:(grid.size-1))+i,"C"],params["x",],zlim=72*fb[1:n_taxa][i]+3,paste(i,":",fb[1:n_taxa][i],":",i_t,sep="")))
    }
    
    # Sum substrate, monomer, enzyme pools
    Substrate_Sum <- c(Substrate_Sum,sum(Substrates.grid[,"C"]))
    Monomer_Sum <- c(Monomer_Sum,sum(Monomers.grid[,"C"]))
    Enzyme_Sum <- c(Enzyme_Sum,sum(Enzymes.grid))
    
    print(i_t)
  }
  
  # Calculate recoveries
  Recovered_C <- (sum(RespSeries) + sum(Substrates.grid[,"C"]) + sum(Monomers.grid[,"C"]) + sum(Enzymes.grid) + sum(Microbes.grid[,"C"]))-(Cum_Microbe["C"] + Cum_Substrate["C"] + Cum_Monomer["C"] + Cum_Enzyme_C)
  Recovered_N <- (Cum_Leaching_N + sum(Substrates.grid[,"N"]) + sum(Monomers.grid[,"N"]) + sum(Enzymes.grid*EnzAttrib[,"N_cost"]) + sum(Microbes.grid[,"N"]))-(Cum_Microbe["N"] + Cum_Substrate["N"] + Cum_Monomer["N"] + Cum_Enzyme_N)
  Recovered_P <- (Cum_Leaching_P + sum(Substrates.grid[,"P"]) + sum(Monomers.grid[,"P"]) + sum(Enzymes.grid*EnzAttrib[,"P_cost"]) + sum(Microbes.grid[,"P"]))-(Cum_Microbe["P"] + Cum_Substrate["P"] + Cum_Monomer["P"] + Cum_Enzyme_P)
  
  index <- is.finite(log10(apply(MicrobesSeries,2,max)))
  #correl <- cor.test(rowSums(EnzGenes)[index],log(apply(MicrobesSeries,2,max)/grid.size)[index])  
  
  # Output grid at end if necessary
  if (is.null(grid.max)) {
    grid.max <- list("Substrates" = Substrates,
                     "Monomers" = Monomers,
                     "Enzymes" = Enzymes,
                     "Microbes" = Microbes)
    grid.time <- params["end_time", ]
  }
  
  # Correlations of taxon abundance from the grid
  corr.vec <- matrix(NA,ncol=2,nrow=n_taxa*n_taxa)
  count <- 1
  taxon.mat <- matrix(grid.max$Microbes[,"C"],ncol=n_taxa,byrow=T)
  for (i in 1:n_taxa) {
    for (j in 1:n_taxa) {
      taxon.a <- taxon.mat[,i]
      taxon.b <- taxon.mat[,j]
      overlap <- which(taxon.a>0 & taxon.b>0)
      taxon.1 <- taxon.a[overlap]
      taxon.2 <- taxon.b[overlap]
      if (length(taxon.1)>2) {
        corr.list <- cor.test(taxon.1,taxon.2)
      } else {
        corr.list <- list("estimate"=NA,"p.value"=NA)
      }
      corr.vec[count,1] <- corr.list$estimate
      corr.vec[count,2] <- corr.list$p.value*((params["n_taxa", ]*params["n_taxa", ]-params["n_taxa", ])/2) # Bonferroni corrected p-values
      count <- count + 1
    }
  }
  
  # Output pulse results
  out.pulse <- list(
    "params"=params,
    "seed"=rng.seed,
    "Ea"=Ea,
    "Uptake_Ea"=Uptake_Ea,
    "Vmax"=Vmax0,
    "Uptake_Vmax"=Uptake_Vmax0,
    "Km"=Km0[1:n_substrates,],
    "Uptake_Km"=Uptake_Km0[1:n_monomers,],
    "ReqEnz"=list(ReqEnz[,,1][1:n_substrates,],ReqEnz[,,2][1:n_substrates,]),
    "Uptake_ReqEnz"=Uptake_ReqEnz[1:n_monomers,],
    "OptimalRatios"=OptimalRatios,
    "RangeRatios"=RangeRatios,
    "EnzGenes"=EnzGenes,
    "EnzProdInduce"=EnzProdInduce,
    "EnzProdConstit"=EnzProdConstit,
    "Enz.track.look"=Enz.track.look,
    "UptakeGenes"=UptakeGenes[1:n_taxa,],
    "UptakeGenesForEnz"=UptakeGenesForEnz,
    "EnzAttrib"=EnzAttrib,
    "MonInput"=MonInput[1:n_monomers],
    "SubInput"=SubInput[1:n_substrates,],
    "CUE"=CUE.ref,
    "RespSeries"=RespSeries,
    "Resp.compSeries"=Resp.compSeries,
    "C.UptakeSeries"=C.UptakeSeries,
    "N.UptakeSeries"=N.UptakeSeries,
    "P.UptakeSeries"=P.UptakeSeries,
    "N.MonomerUptakeSeries"=N.MonomerUptakeSeries,
    "Mic.growthSeries"=Mic.growthSeries,
    "Net.CUE_Series"=Net.CUE_Series,
    "EnzymesSeries"=EnzymesSeries,
    "EnzProdSeries"=EnzProdSeries,
    "SubstratesSeries"=SubstratesSeries,
    "Cum_SubstrateSeries"=Cum_SubstrateSeries,
    "MonProdSeries"=MonProdSeries,
    "C.MonomersSeries"=C.MonomersSeries,            # Only for carbon
    "N.MonomersSeries"=N.MonomersSeries,
    "N.mineralSeries"=N.mineralSeries,
    "P.mineralSeries"=P.mineralSeries,
    "NH4Series"=NH4Series,
    "PO4Series"=PO4Series,
    "MicrobesSeries"=MicrobesSeries,
    "MicrobesNSeries"=MicrobesNSeries,
    "MicrobesPSeries"=MicrobesPSeries,
    "MicrobesNPSeries"=MicrobesNPSeries,
    "grid"=grid.max,
    "grid.time"=grid.time,
    "grid.size"=grid.size,
    "corr.vec"=corr.vec,
    # "cor.test"=correl,
    "end_time"=params["end_time",],
    "Mic_Sum"=Mic_Sum,
    "Substrate_Sum"=Substrate_Sum,
    "Monomer_Sum"=Monomer_Sum,
    "Enzyme_Sum"=Enzyme_Sum,
    "timestamp"=timestamp,
    "fb"=fb[1:n_taxa],
    "Temp"=Temp,
    "Psi"=Psi,
    "Tolerance"=Tolerance,
    "RecoveredC_Resp_Sub_Mon_Enz_Mic"=
      c(Recovered_C,
        sum(RespSeries),
        sum(Substrates.grid[,"C"]),
        sum(Monomers.grid[,"C"]),
        sum(Enzymes.grid),
        sum(Microbes.grid[,"C"])),
    "RecoveredN_Leach_Sub_Mon_Enz_Mic"=
      c(Recovered_N,
        Cum_Leaching_N,
        sum(Substrates.grid[,"N"]),
        sum(Monomers.grid[,"N"]),
        sum(Enzymes.grid*EnzAttrib[,"N_cost"]),
        sum(Microbes.grid[,"N"])),
    "RecoveredP_Leach_Sub_Mon_Enz_Mic"=
      c(Recovered_P,
        Cum_Leaching_P,
        sum(Substrates.grid[,"P"]),
        sum(Monomers.grid[,"P"]),
        sum(Enzymes.grid*EnzAttrib[,"P_cost"]),
        sum(Microbes.grid[,"P"])),
    "UptakeperTaxon" = TUC.mat
  )
  out.pulse
} # End of RunPulse() #########################################################

###############################################################################
TraitModel <- function(job.time,task.ID){
  
  # Read in parameters
  params.file <- paste("params/params", task.ID,".txt", sep = "") 
  params1 <- read.table(params.file, 
                        header = F, 
                        row.names = 1, 
                        sep = "\t", 
                        stringsAsFactors = F)
  params <- params1
  grid.size <- params["x", ]*params["y", ]
  
  seed.val <- as.POSIXct(job.time,format="%y%m%d%H%M%S") + task.ID # convert to date integer, add task id, and use as seed
  set.seed(seed.val) 
  timestamp <- format(seed.val, "%y%m%d%H%M%S") # convert date seed to 12 digit number
  rng.seed <- as.numeric(timestamp)
  if(params["set.seed",]==1) {
    set.seed(as.POSIXct(as.character(params["seed",]),format="%y%m%d%H%M%S"))
    rng.seed <- params["seed",]
  }
  
  n_monomers <- params["n_substrates",]+2
  n_genes <- params["n_enzymes",]
  n_upgenes <- params["n_uptake",]
  
  # Read in activation energies for substrates
  Ea.frame <- read.table("Ea.txt", 
                         header = TRUE, 
                         row.names = 1, 
                         sep = "\t",
                         stringsAsFactors = F)
  
  # Choose taxa with "fungal" strategy, specifically ability to translocate C and nutrients
  fb <- sample(c(1,0), 
               params["n_taxa",], 
               replace = T, 
               prob = c(params["fb",],(1-params["fb",])))
  
  # Set stress tolerance allocation for each taxon
  Tolerance <- c(runif(params["n_taxa",],
                       params["Tol_min",],
                       params["Tol_max",]))
  if(params["Tol_min",]==params["Tol_max",]) {dummy <- runif(params["n_taxa",])}
  
  # Pre-exponential constants for enzymes
  # Rows are substrates; cols are enzymes
  Vmax0 <- c(runif(params["n_substrates",]*params["n_enzymes",],
                   params["Vmax0_min",],params["Vmax0_max",]))
  Vmax0 <- matrix(Vmax0,
                  nrow = params["n_substrates",],
                  ncol = params["n_enzymes",],
                  byrow = T,
                  dimnames = list(sprintf("%s%03d","Sub",1:params["n_substrates",]),
                                  sprintf("%s%03d","Enz",1:params["n_enzymes",])))
  
  # Pre-exponential constants for uptake
  # Rows are monomers; cols are uptake enzymes
  Uptake_Vmax0 <- c(runif(params["n_uptake",]*n_monomers,
                          params["Uptake_Vmax0_min",],
                          params["Uptake_Vmax0_max",]))
  Uptake_Vmax0 <- matrix(Uptake_Vmax0,
                         nrow = n_monomers,
                         ncol = params["n_uptake",],
                         byrow = T,
                         dimnames = list(sprintf("%s%03d","Mon",1:n_monomers),
                                         sprintf("%s%03d","Upt",1:params["n_uptake",])))
  
  # Enzyme specificity matrix of activation energies
  # Rows are substrates; cols are enzymes
  Ea <- sapply(seq(len = dim(Ea.frame)[1]),
               function(i)runif(params["n_enzymes",],
                                Ea.frame$Ea_min[i],
                                Ea.frame$Ea_max[i]))
  Ea <- matrix(Ea,
               nrow = params["n_substrates", ],
               ncol = params["n_enzymes", ],
               byrow = T,
               dimnames = list(sprintf("%s%03d","Sub",1:params["n_substrates",]),
                               sprintf("%s%03d","Enz",1:params["n_enzymes",])))
  
  # Uptake specificity matrix of activation energies
  # Rows are monomers; cols are uptake enzymes
  Uptake_Ea <- c(runif(params["n_uptake",]*n_monomers,
                       params["Uptake_Ea_min",],
                       params["Uptake_Ea_max",]))
  Uptake_Ea <- matrix(Uptake_Ea,
                      nrow = n_monomers,
                      ncol = params["n_uptake",],
                      byrow = T,
                      dimnames = list(sprintf("%s%03d","Mon",1:n_monomers),
                                      sprintf("%s%03d","Upt",1:params["n_uptake",])))
  
  # Enzymes required for substrate degradation
  # Rows are substrates; cols are enzymes
  # Same number within row implies redundancy
  # Ensures each substrate is degraded by at least 1 enzyme and every enzyme degrades at least 1 substrate
  
  ReqEnz1 <- diag(12)
  
  # Choose some substrates that require multiple enzymes
  probability_vector <- rep(0, params["n_enzymes", ])
  if(params["Avg_extra_req_enz",] > 0) probability_vector[1:params["Avg_extra_req_enz", ]] <- 1
  ReqEnz2 <- sample(probability_vector,
                    params["n_substrates", ]*params["n_enzymes", ],
                    replace = T)
  ReqEnz2[ReqEnz1 == 1] <- 0
  ReqEnz2 <- matrix(ReqEnz2,
                    nrow = params["n_substrates", ],
                    ncol = params["n_enzymes", ],
                    byrow = T)
  for(j in 1:dim(ReqEnz2)[1]) {
    # Put in NAs if the substrate does not require a second enzyme
    if(rowSums(ReqEnz2)[j] == 0) ReqEnz2[j, ] <- NA
  }
  
  # Generate the correct structure for the required enzyme matrices
  ReqEnz <- c(expand(ReqEnz1, grid.size), expand(ReqEnz2, grid.size))
  # rows, cols, stacks; default assemble by cols
  ReqEnz <- array(ReqEnz,
                  dim = c(grid.size*params["n_substrates",],
                          params["n_enzymes",],2),
                  dimnames = list(rep(sprintf("%s%03d","Sub",1:params["n_substrates", ]), grid.size),
                                  sprintf("%s%03d","Enz",1:params["n_enzymes",]),c("set1","set2")))
  
  # Substrate concentrations
  substrates.frame <- read.table("substrates.txt", 
                                 header = TRUE, 
                                 sep = "\t",
                                 stringsAsFactors = F,
                                 row.names = 1)
  Substrates <- data.matrix(substrates.frame)
  
  
  # Making matrix to define which monomers track which elements for downregulating induced enzyme production
  ReqEnz3 <- ReqEnz[,,1][1:params["n_substrates",],]
  rownames(ReqEnz3) <- rownames(substrates.frame)
  colSums(ReqEnz3) == 1 # Checking that each enzyme degrades only one substrate
  Subs.track <- data.frame(C = c(0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0), # Defining which substrates tracks which elements
                           N = c(1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0),
                           P = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1), 
                           row.names = rownames((substrates.frame)))
  Enz.track <- t(Subs.track[rep(rownames(ReqEnz3), params["n_enzymes",])[c(ReqEnz3 == 1)],]) # Table of element that regulates induced production of each enzyme
  colnames(Enz.track) <- colnames(ReqEnz3)
  Subs.track.look <- data.frame(element = c("NP", "N", "C", "C", "C", "CN", "CN", "CN", "CN", "CN", "P", "P"),
                                row.names = rownames((substrates.frame)))
  Enz.track.look <- as.vector(t(Subs.track.look[rep(rownames(ReqEnz3), params["n_enzymes",])[c(ReqEnz3 == 1)],]))
  names(Enz.track.look) <- colnames(ReqEnz3)
  
  
  
  
  # Enzymes used for uptake
  # Rows are monomers; cols are uptake enzymes
  # Same number within row implies redundancy
  # Make sure each monomer is taken up by at least one transporter and every transporter takes up at least one monomer
  probability_vector <- rep(0, params["n_uptake", ])
  probability_vector[1:params["Uptake_per_monomer", ]] <- 1
  Uptake_ReqEnz1 <- NULL
  for(i in 1:n_monomers) {
    Uptake_ReqEnz1 <- c(Uptake_ReqEnz1, 
                        sample(probability_vector,
                               params["n_uptake",],
                               replace = F))
  }
  
  
  Uptake_ReqEnz <- diag(14)
  dimnames(Uptake_ReqEnz) <-  list(sprintf("%s%03d","Mon",1:n_monomers),
                                   sprintf("%s%03d","Upt",1:params["n_uptake",]))
    
  
  
  # Ensuring there are no transporters that take up more than one monomer
  for (i in 1:params["n_uptake", ]) {
    if(sum(Uptake_ReqEnz[, i]) > 1) Uptake_ReqEnz[, i] <- sample(probability_vector,
                                                                 n_monomers,
                                                                 replace = F)
  }
  
  # Ensuring there is at least one transporter for each monomer
  for(i in 1:14){
    if(sum(Uptake_ReqEnz[i,]) < 1){
      
      probability_vector <- rep(0, n_monomers)
      probability_vector[i] <- 1
      
      # Identify monomer that is already "covered" by the most amount of transporters
      Mon.cov <- names(which(rowSums(Uptake_ReqEnz) == max(rowSums(Uptake_ReqEnz))))
      # Which transporters take up this monomer?
      Upt.cov <- names(which(colSums(Uptake_ReqEnz[Mon.cov,,drop=F]) > 0))
      # Choose one of these transporters to re-define it as a inorganic nutrient transporter
      Uptake_ReqEnz[,sample(Upt.cov, 1)] <- probability_vector
    }
  }
  
  
  # Substrate input rates
  SubInputC <- read.table("inputs.txt", 
                          header = TRUE, 
                          sep = "\t",
                          stringsAsFactors = F,
                          row.names = 1)[, 1]
  SubInputN <- SubInputC*Substrates[, "N"]/Substrates[, "C"]
  SubInputN[c("OrgP1", "OrgP2")] <- 0
  SubInputP <- SubInputC*Substrates[, "P"]/Substrates[, "C"]
  SubInputP[c("OrgP1", "OrgP2")] <- 0
  SubInput <- c(SubInputC, SubInputN, SubInputP)
  SubInput <- matrix(SubInput,
                     nrow = params["n_substrates", ],
                     ncol = 3,
                     byrow = F,
                     dimnames = list(sprintf("%s%03d",
                                             "SubIn",
                                             1:params["n_substrates",]),
                                     c("C","N","P")))
  SubInput["SubIn001",] <- SubInput["SubIn002",] <- 0
  Substrates["DeadMic",] <- Substrates["DeadEnz",] <- 0
  
  # Monomers produced by each substrate (binary)
  # Rows are substrates; cols are monomers
  # All rows should sum to 1
  MonomersProduced <- cbind(rep(0, params["n_substrates", ]),
                            rep(0, params["n_substrates", ]),
                            diag(params["n_substrates", ]))
  MonomersProduced <- as.matrix(MonomersProduced)
  dimnames(MonomersProduced) <- list(sprintf("%s%03d","Sub",1:params["n_substrates",]),
                                     sprintf("%s%03d","Mon",-1:(n_monomers-2)))   # Original: sprintf("%s%03d","Mon",-1:(n_monomers-2)))
  dimnames(MonomersProduced)[[2]][1:4]<- c("NH4","PO4","DeadMic","DeadEnz")
  
  Uptake_ReqEnz2 <- Uptake_ReqEnz[1:14,]
  Uptake_ReqEnz3 <- Uptake_ReqEnz[1:14,]
  rownames(Uptake_ReqEnz2) <- c(colnames(MonomersProduced)[1:2], rownames(substrates.frame)[1:12])
  rownames(Uptake_ReqEnz3) <- c(colnames(MonomersProduced)[1:2], rownames(substrates.frame)[1:2], rownames(Uptake_ReqEnz)[5:14])
  
  Subs.track.look2 <- 
    rbind(
      data.frame(element = c("N", "P"), row.names = c("NH4", "PO4")),
      Subs.track.look[1:5,, drop = FALSE], 
      data.frame(element = c(rep("CN", 5), rep("P", 2)), row.names = rownames(Subs.track.look)[6:12]))
  
  Subs.track.look2["DeadEnz",] <- "CN"
  
  Subs.track2 <- rbind(c(0, 1, 0), c(0, 0, 1), Subs.track[1:5,], 
                       matrix(c(rep(c(1, 1, 0), 5), rep(c(0, 0, 1), 2)), byrow = TRUE, 
                              ncol = 3, dimnames = list(rownames(Subs.track)[6:12], c("C", "N", "P"))))
  Subs.track2["DeadEnz",] <- c(1, 1, 0)
  rownames(Subs.track2) <- c("NH4", "PO4", rownames(Subs.track))
  
  
  Upt.track <- t(Subs.track2[rep(rownames(Uptake_ReqEnz2), n_upgenes)[c(Uptake_ReqEnz2 == 1)],]) # Table of element that regulates induced production of each transporter
  
  # Making matrix to define which monomers track which elements for downregulating uptake
  Upt.track.look <- as.vector(t(Subs.track.look2[rep(rownames(Uptake_ReqEnz2), n_upgenes)[c(Uptake_ReqEnz2 == 1)],]))
  names(Upt.track.look) <- colnames(Uptake_ReqEnz2)
  
  
  # Optimal stoichiometry of bacterial taxa
  OptimalRatios <- rep(c(params["Cfrac_b",],
                         params["Nfrac_b",],
                         params["Pfrac_b",]),
                       params["n_taxa",])
  OptimalRatios <- matrix(OptimalRatios,
                          nrow = params["n_taxa",],
                          ncol = 3,
                          byrow = T,
                          dimnames = list(sprintf("%s%03d",
                                                  "Tax",
                                                  1:params["n_taxa",]),
                                          c("C","N","P")))
  # Substitute fungal stoichiometry
  OptimalRatios[fb==1,] <- matrix(rep(c(params["Cfrac_f",],
                                        params["Nfrac_f",],
                                        params["Pfrac_f",]),
                                      sum(fb)),
                                  ncol = 3,
                                  byrow = T)
  
  RangeRatios <- rep(c(params["Crange",],
                       params["Nrange",],
                       params["Prange",]),
                     params["n_taxa",])
  RangeRatios <- matrix(RangeRatios,
                        nrow = params["n_taxa",],
                        ncol = 3,
                        byrow = T,
                        dimnames = list(sprintf("%s%03d",
                                                "Tax",
                                                1:params["n_taxa",]),
                                        c("C","N","P")))
  RangeWeight <- RangeRatios/OptimalRatios
  
  # Calculate minimum cell quotas
  MinRatios <- OptimalRatios - RangeRatios
  
  # Enzyme genes possessed by taxa (binary)
  # Rows are taxa; cols are genes
 
if(params["n_taxa",] == 100){
  # First, define categories of number of enzymes possessed, from 0 to 12
  # Must have a fair representation of all categories
  Categories <- 
    list(
      combn(sprintf("%s%03d","Enz",1:n_genes), 0)[, sample(1:choose(12, 0), 1)],
      matrix(combn(sprintf("%s%03d","Enz",1:n_genes), 1)[, sample(1:choose(12, 1), 12)], ncol = 12),
      combn(sprintf("%s%03d","Enz",1:n_genes), 2)[, sample(1:choose(12, 2), 8)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 3)[, sample(1:choose(12, 3), 8)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 4)[, sample(1:choose(12, 4), 8)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 5)[, sample(1:choose(12, 5), 9)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 6)[, sample(1:choose(12, 6), 9)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 7)[, sample(1:choose(12, 7), 8)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 8)[, sample(1:choose(12, 8), 8)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 9)[, sample(1:choose(12, 9), 8)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 10)[, sample(1:choose(12, 10), 8)],
      combn(sprintf("%s%03d","Enz",1:n_genes), 11)[, sample(1:choose(12, 11), 12)],
      matrix(combn(sprintf("%s%03d","Enz",1:n_genes), 12)[, sample(1:choose(12, 12), 1)], ncol = 1)
    )
  
  # Second, define genes in each taxon
  EnzGenes <- matrix(ncol = n_genes, nrow = params["n_taxa",], dimnames = list(sprintf("%s%03d", "Tax", 1:params["n_taxa",]), 
                                                                       sprintf("%s%03d","Enz",1:n_genes)))
  
  EnzGenes[1,] <- rep(0, 12)
  EnzGenes[2:13,] <- t(apply(Categories[[2]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[14:21,] <- t(apply(Categories[[3]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[22:29,] <- t(apply(Categories[[4]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[30:37,] <- t(apply(Categories[[5]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[38:46,] <- t(apply(Categories[[6]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[47:55,] <- t(apply(Categories[[7]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[56:63,] <- t(apply(Categories[[8]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[64:71,] <- t(apply(Categories[[9]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[72:79,] <- t(apply(Categories[[10]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[80:87,] <- t(apply(Categories[[11]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[88:99,] <- t(apply(Categories[[12]], 2, function(x) as.numeric(sprintf("%s%03d","Enz",1:n_genes) %in% x)))
  EnzGenes[100,] <- rep(1, 12)
  
}
  
if(params["n_taxa",] == 1 & params["Org.recycl",] == 1){
  
  EnzGenes <- matrix(ncol = n_genes, nrow = params["n_taxa",], dimnames = list(sprintf("%s%03d", "Tax", 1:params["n_taxa",]), 
                                                                               sprintf("%s%03d","Enz",1:n_genes)))
  EnzGenes[1,] <- rep(1, 12)
  
}
  
if(params["n_taxa",] == 1 & params["Org.recycl",] == 0){
  
  EnzGenes <- matrix(ncol = n_genes, nrow = params["n_taxa",], dimnames = list(sprintf("%s%03d", "Tax", 1:params["n_taxa",]), 
                                                                               sprintf("%s%03d","Enz",1:n_genes)))
  EnzGenes[1,] <- c(0, 0, rep(1, 10))
  
}
  
  
  Normalize <- 1
  if (params["NormalizeProd",] == 1) {
    Normalize <- rowSums(EnzGenes)/params["Enz_per_taxon_max",]
  }
  
  
  EnzProdInduce <- apply(EnzGenes, 1, function(x)runif(1,
                                                       params["Enz_Prod_min",],
                                                       params["Enz_Prod_max",]))
  EnzProdInduce <- EnzProdInduce*EnzGenes/Normalize
  EnzProdInduce[is.na(EnzProdInduce)] <- 0
  EnzProdConstit <- apply(EnzGenes, 1, function(x)runif(1,
                                                        params["Constit_Prod_min",],
                                                        params["Constit_Prod_max",]))
  EnzProdConstit <- EnzProdConstit*EnzGenes/Normalize
  EnzProdConstit[is.na(EnzProdConstit)] <- 0
  
  # Uptake genes possessed by taxa (binary)
  # Rows are taxa; cols are genes
  # All taxa must have at least one uptake gene
  # For each taxon, check if it has the enz gene associated with the monomer(s) to be targeted by the uptake gene
  # If so, assign that uptake gene to the taxon with probability p
  # Enz-Taxon%*%Enz-Sub%*%Sub-Monomer to get Taxon-Monomer matrix
  # Matrix multiply this by the required uptake enzyme matrix, then assign uptake genes
  RE2 <- ReqEnz[,,2][1:params["n_substrates",],]
  RE2[is.na(RE2)] <- 0
  enz_sub <- ReqEnz[,,1][1:params["n_substrates",],] + RE2
  # # Matrix multiplication to relate taxa to the monomers they can generate with their enzymes
  # UptakeGenes <- EnzGenes%*%t(enz_sub)%*%Uptake_ReqEnz[3:14,]    # Original was MonomersProduced instead of Uptake_ReqEnz[-1:-2,]
  # UptakeGenes[,1:2] <- 1
  # UptakeGenes[UptakeGenes>0] <- 1
  
  # All taxa can take up all monomers
  UptakeGenes <- matrix(1, ncol = 14, nrow = params["n_taxa",], dimnames = list(rownames(EnzGenes), colnames(Uptake_ReqEnz)))
  
  # Make sure every taxon is likely to have an uptake enzyme for at least 1 organic monomer
  # Not guaranteed unless uptake_prob = 1
  probability_vector <- rep(0,n_upgenes-2)
  probability_vector[1] <- 1
  for (i in 1:params["n_taxa",]) {
    if(sum(UptakeGenes[i,3:n_upgenes])==0) UptakeGenes[i,3:n_upgenes] <- sample(probability_vector,
                                                                                n_upgenes - 2,
                                                                                replace = F)
  }
  
  # Give taxa random number of additional uptake genes between the number they have and n_upgenes
  for (i in 1:params["n_taxa", ]) {
    n.zero <- length(UptakeGenes[i, ][UptakeGenes[i, ] == 0])
    probability_vector <- rep(0,n.zero)
    probability_vector[1:sample(0:n.zero,1)] <- 1
    UptakeGenes[i,][UptakeGenes[i,]==0] <- sample(probability_vector, 
                                                  n.zero,
                                                  replace = F)
  }
  UptakeGenesForEnz <- UptakeGenes
  # If true then the uptake potential is normalized to the number of uptake genes
  Normalize <- 1
  if (params["NormalizeUptake",]==1) {
    Normalize <- rowSums(UptakeGenes)/params["n_uptake",]
  }
  UptakeGenes <- UptakeGenes/Normalize
  # Choose total amount of uptake allocation at random
  UptakeProd <- c(runif(params["n_taxa",],
                        params["Uptake_C_cost_min",],
                        params["Uptake_C_cost_max",]))
  UptakeGenes <- UptakeGenes*UptakeProd
  
  # Calculate CUE as a function of the number of enzyme genes, uptake genes, and reference CUE
  CUE.ref <- params["CUE_enz",]*rowSums(EnzGenes)/(params["Enz_per_taxon_max",]) + 
    params["CUE_uptake",]*rowSums(UptakeGenesForEnz)/(params["n_uptake",]) + params["CUE_ref",] + params["Tol_CUE",]*Tolerance
  
  # Enzyme attributes
  # Rows are enzymes; cols are attributes
  EnzAttrib_C <- rep(params["Enz_C_cost",],
                     params["n_enzymes",])
  EnzAttrib <- rbind(EnzAttrib_C,
                     params["Enz_N_cost",],
                     params["Enz_P_cost",],
                     params["Enz_Maint_cost",])
  EnzAttrib <- matrix(EnzAttrib,
                      nrow = params["n_enzymes",],
                      ncol = 4,
                      byrow = T,
                      dimnames = list(sprintf("%s%03d",
                                              "Enz",
                                              1:params["n_enzymes",]),
                                      c("C_cost","N_cost","P_cost","Maint_cost")))
  
  # Monomer pool sizes for all elements
  Monomers <- rbind(c(0,params["Init_NH4",],0),
                    c(0,0,params["Init_PO4",]),
                    Substrates*(params["Monomer_Substrate_Ratio",]))
  Monomers[is.na(Monomers)] <- 0
  Monomers <- as.matrix(Monomers)
  dimnames(Monomers) <- list(c("NH4","PO4","DeadMic","DeadEnz",
                               sprintf("%s%03d","Mon",3:(n_monomers-2))),
                             c("C","N","P"))
  
  # Monomer input rates
  MonInput <- c(params["Input_NH4",],
                params["Input_PO4",],
                read.table("inputs.txt", 
                           header = TRUE, 
                           sep = "\t",
                           stringsAsFactors = F,
                           row.names = 1)[, 2])
  MonInput <- as.matrix(MonInput)
  rownames(MonInput) <- rownames(Monomers)
  
  # Enzyme pool sizes for all elements
  Enzymes <- matrix(c(runif(params["n_enzymes",],
                            params["Enz_min",],
                            params["Enz_max",])),
                    nrow = params["n_enzymes",],
                    ncol = 1,
                    dimnames = list(sprintf("%s%03d","Enz",1:params["n_enzymes",]),"C"))
  
  # Microbial pool sizes for all elements
  BacC <- 0.5*params["max_size_b",]
  BacN <- BacC*params["Nfrac_b",]/params["Cfrac_b",]
  BacP <- BacC*params["Pfrac_b",]/params["Cfrac_b",]
  FunC <- 0.5*params["max_size_f",]
  FunN <- FunC*params["Nfrac_f",]/params["Cfrac_f",]
  FunP <- FunC*params["Pfrac_f",]/params["Cfrac_f",]
  Microbes <- rep(c(BacC,BacN,BacP),params["n_taxa",]*grid.size)
  Microbes <- matrix(Microbes,
                     nrow = params["n_taxa", ]*grid.size,
                     ncol = 3,
                     byrow = T,
                     dimnames = list(rep(sprintf("%s%03d","Tax",1:params["n_taxa",]),
                                         grid.size),
                                     c("C","N","P")))
  fb <- rep(fb, grid.size)
  Microbes[fb == 1, "C"] <- FunC; Microbes[fb == 1, "N"] <- FunN; Microbes[fb == 1, "P"] <- FunP
  # Randomly assign microbes to each grid box
  p.b <- params["taxa_per_box",]
  p.f <- p.b*params["max_size_b",]/params["max_size_f",]
  choose_taxa <- sample(c(1,0),
                        params["n_taxa",]*grid.size,
                        replace = T, 
                        c(p.b,(1 - p.b)))
  choose_taxa[fb == 1] <- sample(c(1, 0), 
                                 sum(fb), 
                                 replace = T,
                                 c(p.f, (1 - p.f)))
  Microbes[choose_taxa == 0, ] <- 0
  # Initialize the reproduction list
  ReprodNew <- Microbes
  ReprodNew[ReprodNew > 0] <- 0
  Colonization.reset <- ReprodNew
  
  # Account for efficiency-specificity tradeoff by dividing Vmax_0 by the number of substrates (or monomers) targeted
  # and multiplied by a specificity factor
  # Leave Km unchanged, effectively increasing it by the factor that Vmax_0 is reduced
  total_substrates <- colSums(ReqEnz[,,1][1:params["n_substrates",],]) + colSums(RE2)
  if (params["Specif_factor",]==0) {total_substrates[total_substrates>1] <- 1} else
  {total_substrates[total_substrates>1] <- total_substrates[total_substrates>1]*params["Specif_factor",]}
  Vmax0 <- t(t(Vmax0)/total_substrates)
  Vmax0[!is.finite(Vmax0)] <- 0
  
  total_monomers <- colSums(Uptake_ReqEnz)
  if (params["Specif_factor",]==0) {total_monomers[total_monomers>1] <- 1} else
  {total_monomers[total_monomers>1] <- total_monomers[total_monomers>1]*params["Specif_factor",]}
  Uptake_Vmax0 <- t(t(Uptake_Vmax0)/total_monomers)
  Uptake_Vmax0[!is.finite(Uptake_Vmax0)] <- 0
  
  # Implement Vmax-Km tradeoff as a direct correlation with slope Vmax_Km, 
  # and error term normally distributed with magnitude Km_error. Minimum Km constrained to Km_min
  Km <- abs(rnorm(Vmax0, 
                  mean = Vmax0*params["Vmax_Km",], 
                  sd = mean(Vmax0)*params["Km_error",]) + params["Vmax_Km_int", ])
  Km[Km < params["Km_min", ]] <- params["Km_min", ]
  Km <- matrix(Km,nrow = params["n_substrates", ],
               ncol = params["n_enzymes", ],
               byrow = F,
               dimnames = list(sprintf("%s%03d","Sub",1:params["n_substrates",]),
                               sprintf("%s%03d","Enz",1:params["n_enzymes",])))
  Uptake_Km <- abs(rnorm(Uptake_Vmax0, 
                         mean = Uptake_Vmax0*params["Uptake_Vmax_Km",], 
                         sd = mean(Uptake_Vmax0)*params["Km_error",]) + params["Uptake_Vmax_Km_int",])
  Uptake_Km[Uptake_Km < params["Uptake_Km_min",]] <- params["Uptake_Km_min",]
  Uptake_Km <- matrix(Uptake_Km,
                      nrow = n_monomers,
                      ncol = params["n_uptake",],
                      byrow = F,
                      dimnames = list(sprintf("%s%03d","Mon",1:n_monomers),
                                      sprintf("%s%03d","Upt",1:params["n_uptake",])))
  
  # Expand matricies to cover the whole grid
  Substrates <- expand(Substrates, grid.size)
  SubInput <- expand(SubInput, grid.size)
  Monomers <- expand(Monomers, grid.size)
  MonInput <- expand(MonInput, grid.size)
  Enzymes <- expand(Enzymes, grid.size)
  Km0 <- expand(Km, grid.size)
  Uptake_Km0 <- expand(Uptake_Km, grid.size)
  Uptake_ReqEnz <- expand(Uptake_ReqEnz, grid.size)
  MinRatios <- expand(MinRatios, grid.size)
  
  climate <- read.table("climate2.txt", 
                        header = TRUE, 
                        sep = "\t", 
                        stringsAsFactors = F)
  
  # Set up output list containing each pulse
  out <- vector("list", params["pulses",])
  for (i_p in 1:params["pulses",]) {
    zeros <- ifelse(i_p < 10, "0", "")
    filename <- paste("outputs/", timestamp, "_", zeros, i_p, "G.png", sep = "")
    rows <- params["end_time",]%/%params["print.grid",] + 1
    png(file=filename,width=300*params["n_taxa",],height=300*rows)
    par(mfrow=c(rows,params["n_taxa",]),mar=c(3,3,3,3))
    Temp <- climate$Temp[(1+params["end_time",]*(i_p-1)):(params["end_time",]*i_p)]
    Psi <- climate$Psi[(1+params["end_time",]*(i_p-1)):(params["end_time",]*i_p)]
    out[[i_p]] <- RunPulse(
      params,
      timestamp,
      rng.seed,
      grid.size,
      Microbes,
      Substrates,
      SubInput,
      Enzymes,
      Monomers,
      MonInput,
      MonomersProduced,
      ReprodNew,
      Colonization.reset,
      Ea,
      Vmax0,
      Km0,
      ReqEnz,
      EnzGenes,
      EnzProdInduce,
      EnzProdConstit,
      Enz.track.look,
      Upt.track.look,
      UptakeGenes,
      UptakeGenesForEnz,
      Uptake_ReqEnz,
      EnzAttrib,
      Uptake_Ea,
      Uptake_Vmax0,
      Uptake_Km0,
      CUE.ref,
      OptimalRatios,
      MinRatios,
      RangeRatios,
      fb,
      Temp,
      Psi,
      Tolerance
    )
    dev.off()
    filename <- paste("outputs/", timestamp, "_", zeros, i_p, ".png", sep = "")
    png(file=filename,width=2000,height=800)
    MakePlots(out[[i_p]])
    dev.off()
    # Recalculate microbial frequencies and repopulate the grid with microbes based on frequencies
    cum.abundance <- apply(out[[i_p]]$MicrobesSeries, 2, sum)
    # Account for different cell sizes of bacteria and fungi
    cum.abundance[out[[i_p]]$fb==1] <- cum.abundance[out[[i_p]]$fb==1]*params["max_size_b",]/params["max_size_f",]
    frequencies <- cum.abundance/sum(cum.abundance)
    
    # Microbial pool sizes for all elements
    BacC <- 0.5*params["max_size_b",]
    BacN <- BacC*params["Nfrac_b",]/params["Cfrac_b",]
    BacP <- BacC*params["Pfrac_b",]/params["Cfrac_b",]
    FunC <- 0.5*params["max_size_f",]
    FunN <- FunC*params["Nfrac_f",]/params["Cfrac_f",]
    FunP <- FunC*params["Pfrac_f",]/params["Cfrac_f",]
    Microbes <- rep(c(BacC,BacN,BacP),params["n_taxa",]*grid.size)
    Microbes <- matrix(Microbes,
                       nrow = params["n_taxa",]*grid.size,
                       ncol = 3,
                       byrow = T,
                       dimnames = list(rep(sprintf("%s%03d","Tax",1:params["n_taxa",]),
                                           grid.size),
                                       c("C","N","P")))
    Microbes[fb == 1, "C"] <- FunC; Microbes[fb == 1, "N"] <- FunN; Microbes[fb == 1, "P"] <- FunP
    # Randomly assign microbes to each grid box based on prior densities
    probs <- matrix(c(frequencies, (1 - frequencies)), ncol = 2)
    choose_taxa <- matrix(rep(0, grid.size*params["n_taxa", ]), nrow = params["n_taxa", ])
    for (i in 1:params["n_taxa", ]) {
      choose_taxa[i,] <- sample(c(1, 0), grid.size, replace = T, probs[i, ])
    }
    Microbes[as.vector(choose_taxa) == 0, ] <- 0
    # Reset reproduction
    ReprodNew[ReprodNew > 0] <- 0  
  }
  out
  
} # End of TraitModel() #######################################################

