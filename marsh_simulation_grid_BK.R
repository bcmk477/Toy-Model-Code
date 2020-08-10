setwd("~/Desktop/PhD/Summer 2020 Research/Marsh Modeling")
getwd

library(simecol)

#Practicing simecol state-change matrices ####
wdist<- matrix(c(0.5,0.5,0.5,0.5,0.5,
                 0.5,1.0,1.0,1.0,0.5,
                 0.5,1.0,1.0,1.0,0.5,
                 0.5,1.0,1.0,1.0,0.5,
                 0.5,0.5,0.5,0.5,0.5), nrow = 5)
n <- 20
m <- 20

x <- matrix(rep(0, m*n), nrow = n)

x[10, 10] <- 1
x[1,5] <- 1
x[n, 15] <- 1
x[5,2] <- 1
x[15, m] <- 1
x[n,1] <-1

opar<- par(mfrow= c(2,2))

image(x)
image(matrix(neighbours(x, wdist = wdist, bounds = 0), nrow = n))
image(matrix(neighbours(x, wdist = wdist, bounds = 1), nrow = n))
image(matrix(neighbours(x, wdist = wdist, bounds = c(0, 1, 0, 1)), nrow = n))
par(opar)

#Altered simecol sample code with our parameters####

mycolors <- function(n) {
  col <- c("white", "darkblue")
  if (n>2) col <- c(col, heat.colors(n - 2))
  col
}

prob_bord <- 0.5  # probability of a bare cell becoming a border when neighboring cell is a border
prob_patch <- 0.7  # probability of a bare cell becoming a patch when neighboring cell is a patch

patch_cover <- 0.075  

bare = 1
patch = 2
border = 3

n <- 100
m <- 100
x <- rep(0, m*n)

x[round(runif(100,1,m*n))] <- patch
dim(x)<- c(n, m)

x[38:42, 38:42] <- 2

image(x, col = mycolors(2))


#red = senescence
#orange = adults
#yellow = juveniles
#green = T0

for (i in 1:100){
  ad <- ifelse(x >= patch & x < border, x, 0)
  nb <- neighbours(ad, wdist = wdist)
  genprob <- nb * runif(nb) * patch_cover
  xgen  <- ifelse(x == 0 & genprob >= 1, 1, 0)
  ## growth if neighboring cells are bare
  #xbare <- ifelse(x >= 1 & x < patch & runif(x) <= prob_bord, x+1, 0)
  ## growth if neighboring cells are patches
  xpatch <- ifelse(x >= patch & x < border & runif(x) <= prob_patch, x+1, 0)
  ## growth if neighboring cells are borders
  xbord <- ifelse(x >= border & runif(x) <= prob_bord, x+2, 0)
  ## make resulting grid of salt marsh
  x     <- xgen + xpatch + xbord
  ## plot resulting grid
  image(x, col = mycolors(max(x)), add = T)
  #if (max(x) == 0) stop("extinction", call. = FALSE)
}

#did populate graph previously, but no longer does (not sure why)