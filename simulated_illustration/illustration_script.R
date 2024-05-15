library(here)

# FUNCTIONS

# Package Installation
installer <- function(vector) {
  for(i in vector) {
    if(!require(i, character.only = T)) install.packages(i)
    library(i, character.only = T)
  }
}

# Condensed Formula Object
form <- function(covariates, outcome = NULL, quotes = T) {
  equation <- as.formula(paste0(outcome, " ~ NULL"))
  for(i in covariates) {
    if(substr(i, 1, 1) %in% c(LETTERS, letters) & substr(i, nchar(i), nchar(i)) %in%  c(LETTERS, letters) |
       quotes == T) {
      equation <- update(equation, as.formula(paste0(" ~ . + `", i,"`"))) 
    }
    else{
      equation <- update(equation, as.formula(paste0(" ~ . + ", i))) 
    }
  }  
  equation
}

# BASE PACKAGES
installer(c("tidyverse", "here", "monoreg"))



# DATA PREPARATION
load(here("simulated_illustration", "data", "sim_illustration.RData"))

wide.data <- sim_illustration$`Wide-Format Data`
long.data <- sim_illustration$`Long-Format Data`

# Grid of Doses and Volumes
dose.range <- c(0, 70) # Range of Dose Distribution
dvals <- seq(dose.range[1], dose.range[2], by = 1)  # 71 (dose indices)
gvals <- seq(0, 1, by = 0.02)   # 51 (volume indices)
num.d <- length(dvals) # Number of Dose Indices
num.g <- length(gvals) # Number of Volume Indices

# Visualization/Test Set Purposes
grid.num.d <- 24
grid.num.g <- 24
dvals.grid <- seq(dose.range[1], dose.range[2], length.out = grid.num.d + 2)[-c(1, grid.num.d + 2)]
gvals.grid <- seq(0, 1, length.out = grid.num.g + 2)[-c(1, grid.num.g + 2)]

# Stochastic Intervention for Bladder DVHs
# Volume at 42 Gy <= 50%
q_bladder <- 0.5
d.star_bladder <- 42; d.star_bladder <- dvals[findInterval(d.star_bladder, dvals)]
d.star.data_bladder <- long.data[which(long.data$Dose == d.star_bladder),]


# Stochastic Intervention for Skin DVHs
# Volume at 20 Gy <= 20%
q_skin <- 0.2
d.star_skin <- 20; d.star_skin <- dvals[findInterval(d.star_skin, dvals)]
d.star.data_skin <- long.data[which(long.data$Dose == d.star_skin),]


# FIGURE 1 (SUPPLEMENTARY FIGURE 9 analogue): VISUALIZATION OF DVHs
source(here("simulated_illustration", "figures", "Figure 1", "figure1_code.R"))

# MARGINAL MODELING (based on observed data [1 sample])
source(here("simulated_illustration", "illustration_modeling_bladder.R"))
source(here("simulated_illustration", "illustration_modeling_skin.R"))

# BOOTSTRAP (1000 resamples)
source(here("simulated_illustration", "illustration_bootstrap_bladder.R"))   # Bootstrap Script (Bladder)
source(here("simulated_illustration", "illustration_bootstrap_skin.R"))      # Bootstrap Script (Skin)

# BOOSTRAP VISUALIZATIONS
source(here("simulated_illustration", "illustration_bootstrap_visualization_bladder.R"))   # Bladder (Figure 2, Table 1)
source(here("simulated_illustration", "illustration_bootstrap_visualization_bladder.R"))   # Bladder (Figure 3, Table 1)





