# Plot vertical wiggle plots with shaded positive and negative values for 
# equally spaced virtual source gathers.


# Load functions
source("../functions/plotting_functions.R")

vsg <- as.matrix(read.csv("virtual_gather.csv", header = F) )

m <- dim(vsg)[1] 

time_vec <- seq(-floor( (m-1)/2), ceiling( (m-1)/2 ) )
vsg <- vsg[,-1]


line_fill <- c("red", "white")
k <- 751
ylims <- c(-1500, 1500)
plot_geophone_wiggles( time_vec, ts_norm(vsg, k) , line_fill )
