#!/usr/bin/RScript

library(RSEIS)
library(fields)
library(oce)
library(Hmisc)

# Load stations pairs textfile
stp <- read.delim("station_pairs.txt", header = F, sep = "\t")

# Load the EGFs, assign the time vector and remove it from the data
egf <- read.csv("EmpiricalGreensFunctions.csv", header = F)
t_int <- egf[,1] 
egf <- egf[,-1]

# Load corresponding spectral data, assign the frequency vector and remove it 
# from the data
psd <- read.csv("EGFspectrum.csv", header = F)
frq <- psd[,1] 
psd <- psd[,-1]

# Sort the data based on their positions 
ind_sort <- sort( stp[,3], index.return = T)$ix

stp <- stp[ind_sort,]
egf <- egf[,ind_sort]
psd <- psd[,ind_sort]

# Get the locations. We may need to average duplicates 
locs <- stp[,3] 
lid <- unique(locs) 

n <- length(lid) 

# Check for duplicates
stack_arr <- 0
rm_ind <- 0
for(i in 1:n) 
{
  L <- locs == lid[i]  
  if( sum(L) > 1 )
  {
    stack_arr <- c(stack_arr, lid[i])
  }
}

if(length(stack_arr) > 1)
{
  for( i in 2:length(stack_arr) )
  {
    ind <- which(locs == stack_arr[i])
    egf_tmp <- egf[,ind ]
    # align and stack to the first egf 
    for(j in 2:( dim(egf_tmp)[2]) )
    {
      z <- xcor2(egf_tmp[,1], egf_tmp[,j], DT = 1, LAG = 100, PLOT = T) 
      # shift the second timeseries in the cross correlation
      
      
      if(z$mlag2 < 0)
      {
        tmp <- egf_tmp[,j] 
        tmp <- c(tmp[-1:z$mlag2], rep(0, abs(z$mlag2) ) )
        egf_tmp[,j] <- tmp
      }else if(z$mlag2 > 0) 
      {
        tmp <- egf_tmp[,j] 
        tmp <- c(rep(0, abs(z$mlag2) ), tmp[1:(length(tmp) - z$mlag2)]  )
        egf_tmp[,j] <- tmp
      }else 
      {
        next
      }
    }
    
    egf[,ind[1] ] <- apply(egf_tmp, 1, mean)
    rm_ind <- c( rm_ind, ind[2:length(ind)])
  }
}

if(length(rm_ind) > 1)
{
  rm_ind <- rm_ind[-1]
  egf <- egf[,-rm_ind]
  psd <- psd[,-rm_ind] 
  stp <- stp[-rm_ind,]
}
X <- stp[,4]

# normalize the psd by the first reciever for each line
norm_factor <- max(psd[,11])/max(psd[12])


zmax <- ceiling(max(psd[1:770,]))
zmin <- floor(min(psd[1:770,]))

psd[,1:11] = psd[,1:11]/(norm_factor*max(psd[,11]))
psd[,12:32] = psd[,12:32]/max(psd[,12])

###############################################################################
##                            Plot the spectrum                              ##
###############################################################################
dev.new()
psd <- psd[1:770,] 
frq <- frq[1:770]

log_frq <- c(0, log10(frq[2:length(frq)]) )
ylims <- c(max(log_frq), min(log_frq) )
image(X,frq, t(as.matrix(psd)) ,col = tim.colors(100) , xaxt = "n", xlab = "",
      ylab = "Frequency (Hz)", ylim = c(max(frq), min(frq) ) )
lines(rep(0, length(frq) ), frq)
lines(rep(110, length(frq)), frq)
#lines(c(-140,X,220), rep(450, length(X)+2 ) )
axis(3)
mtext("Distance From Source (m)", side = 3, padj = -4)

col.labels <- c(zmin,zmax)
color.legend(-20,1550,120,1600, c(""), rect.col=tim.colors(100)) 
mtext("Gain (dB)", side = 1, padj = 4.5)
mtext(zmin, side = 1, padj = 2.6, adj = 0.275)
mtext(zmax, side = 1, padj = 2.6, adj = 0.725)
# col = colorRampPalette(c("blue","orange", "red4"))( (100) )
axis(3, at = seq(-110, 210, by = 10), lwd.ticks = 0.5, labels = F)


###############################################################################
##                              Plot the EGF's                               ##
###############################################################################

pm_polygon_fill <- function(x, y, line_fill) 
{
  ind <- which(y >= 0)
  pos_x <- c(x[ind], rev(x[ind]) )
  pos_y <- c(y[ind], rep(0, length(ind) ) )
  
  neg_x <- c(x[-ind], rev(x[-ind]) )
  neg_y <- c(y[-ind], rep(0, length(y[-ind]) ) )
  
  polygon(pos_y, pos_x, col = line_fill[1])
  polygon(neg_y, neg_x, col = line_fill[2])
  
}


# define a few plotting parameters
ylims <- c(max(t_int), min(t_int) )
line_fill <- c("red4", "yellow")

dev.new()
par(mfrow = c(1, length(X) ), mar = c(0.25, 0.10, 0.25, 0.0), mai = c(0, 0, 0, 0)  )

for(i in 1:length(X) )
{
  plot(egf[,i], t_int, ylim = ylims, type = "l",bg = "transparent", 
       xaxt = "n", yaxt = "n", bty = "n")
  pm_polygon_fill(t_int, egf[,i], line_fill)
}
