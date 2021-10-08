# Supplementary material for:
# Smoothing spatio-temporal data withcomplex missing data patterns
# by Eleonora Arnone, Laura M. Sangalli, Andrea Vicini
# 
# Section 6

# Load the data for the application
# Data available at http://www.laketemp.net/home_ARCLake/targets_phase3.php

load("LakeVictoria.RData")
# locations: the spatial locations where data may be observed
# timelocations: the time instants where data may be observed
# DataMat: the #locations-by-#timelocations matrix containing the LWT data.
#          NAs correspond to unobserved data
# nodes: the starting nodes the create the mesh
# segments: the segments for defining the boundary of the domain
# holes: the holes in the domain
# TimePoints: the time instants used as B-splines nodes


mesh  = create.mesh.2D(nodes = nodes, segments = segments, holes = holes)
meshr = refine.mesh.2D(mesh,minimum_angle = 30,maximum_area = 0.03)
FEMbasis = create.FEM.basis(meshr)

lambdaT = 1e-3
lambdaS = 1e-3

solution = smooth.FEM.time(locations = locations,
                            time_mesh = TimePoints,
                            time_locations = timelocations,
                            observations = DataMat,
                            FEMbasis = FEMbasis,
                            lambdaS = lambdaS,
                            lambdaT = lambdaT)

# Conversion time-to-date
dates=as.Date("1970-01-01 00:00:00")+time+9296.5
dates=format(dates, "%d %b %Y")

# image function to plot the results
library(viridis)
image.fixedTime=function(time, SolutionObj,lambdaSindex=NULL,lambdaTindex=NULL, Nx = 100, Ny = 100, xlim = NA, ylim = NA, zlim = NULL,out=NULL, nlev = 10, ...)
{
  if(is.na(xlim[1]))
  {
    xmin = min(SolutionObj$fit.FEM.time$FEMbasis$mesh$nodes[,1])
    xmax = max(SolutionObj$fit.FEM.time$FEMbasis$mesh$nodes[,1])
  }
  else
  {
    xmin = xlim[1]
    xmax = xlim[2]
  }
  
  if(is.na(ylim[1]))
  {
    ymin = min(SolutionObj$fit.FEM.time$FEMbasis$mesh$nodes[,2])
    ymax = max(SolutionObj$fit.FEM.time$FEMbasis$mesh$nodes[,2])
  }
  else
  {
    ymin = ylim[1]
    ymax = ylim[2]
  }
  
  X    = matrix(seq(xmin, xmax, len=Nx),ncol=1)
  Y    = matrix(seq(ymin, ymax, len=Ny),ncol=1)    
  
  Xmat = X %*% matrix(1,nrow=1,ncol=Ny)
  Ymat = matrix(1,nrow=Nx,ncol=1) %*% t(Y)
  Xvec = NULL
  for (numc in 1:Ny)
  {
    Xvec=c(Xvec,Xmat[,numc])
  }
  Yvec = NULL
  for (numc in 1:Ny)
  {
    Yvec=c(Yvec,Ymat[,numc])
  }
  
  if(is.null(lambdaSindex))
    lambdaSindex=SolutionObj$bestlambda[1]
  if(is.null(lambdaTindex))
    lambdaTindex=SolutionObj$bestlambda[2]
  
  eval_points = cbind(rep(time,length(Xvec)),Xvec, Yvec)
  eval_sol=rep(NA,nrow(eval_points))
  if(!is.null(out))
    eval_points = eval_points[-out,]
  eval_sol_in = fdaPDE::eval.FEM.time(SolutionObj$fit.FEM.time,space.time.locations = eval_points,lambdaS = lambdaSindex,lambdaT = lambdaTindex)
  if(!is.null(out))
    eval_sol[(1:length(eval_sol))[-out]]=eval_sol_in
  else{
    eval_sol=eval_sol_in
  }
  outside=which(is.na(eval_sol))
  evalmat = matrix(eval_sol, nrow=Nx, ncol=Ny, byrow=F)
  
  if(is.null(zlim))
  {
    zlim[1] = min(eval_sol, na.rm = TRUE)
    zlim[2] = max(eval_sol, na.rm = TRUE)
  }
  
  image2D(z=evalmat,x=as.vector(X),y=as.vector(Y),zlim=zlim,...)
  contour2D(z=evalmat,x=as.vector(X),y=as.vector(Y),zlim=zlim, add=T,colkey=F,col="black",nlevels=nlev)
  # outside
}

# plot of the solution in a fixed time instant
temp.ind = 4 # september 1995
image.fixedTime(time[temp.ind],solution,lambdaSindex=NULL,lambdaTindex=NULL, 
                Nx = 65, Ny = 66, zlim = c(296.2766,301.3383), xlab = "", ylab = "", asp=1, 
                main = paste0(dates[temp.ind]),col=viridis(256),axes=F)
      