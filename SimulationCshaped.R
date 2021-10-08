# Supplementary material for:
# Smoothing spatio-temporal data withcomplex missing data patterns
# by Eleonora Arnone, Laura M. Sangalli, Andrea Vicini
# 
# Section 5.4

rm(list=ls())
library(fdaPDE)
# version needed: 1.1-2 or later
# most updated library version always available at https://github.com/fdaPDE/fdaPDE
# periodically released on CRAN: https://cran.r-project.org/web/packages/fdaPDE/index.html
library(mgcv); library(sp); library(spacetime); library(gstat); library(INLA);

source("censoring.R")
data("horseshoe2D")

##################   Preparation   ##################   
schema="a" 
# schema="b"; nRDD=8
# schema="c"; nRDD=12
# schema="d"; nRDD=40
tf=1; p=0.5; numIter=30; beginIter=1

# function to generate data
f<-function(x,y,t)
{
  K <- (y/0.1*as.double((abs(y)<=0.1 & x>-0.5))+as.double((abs(y)>0.1 | x<=-0.5)))^2
  res=numeric(length =length(x))
  for(i in 1:length(x))
  {
    if(x[i]>=0 && y[i]>0) 
      res[i]=cos(t[i]/tf*pi)*(0.25*pi+x[i])+(y[i]-0.5)^2
    if(x[i]>=0 && y[i]<=0) 
      res[i]=cos(2*t[i]/tf*pi)*(-0.25*pi-x[i])+(-y[i]-0.5)^2
    if(x[i]<0 && y[i]>0) 
      res[i]=cos(t[i]/tf*pi)*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
    if(x[i]<0 && y[i]<=0) 
      res[i]=cos(2*t[i]/tf*pi)*(-atan(y[i]/x[i])*0.5)+(sqrt(x[i]^2+y[i]^2)-0.5)^2*K[i]
  }
  res
}

# create the vector of complete time instants and the matrix of complete observations
NumTimeInstants = 41
timelocations = seq(0, tf, length.out = NumTimeInstants) 

meshh = create.mesh.2D(nodes = horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments)
meshh = refine.mesh.2D(meshh,minimum_angle = 30,maximum_area = 0.04)
locations = meshh$nodes[which(meshh$nodesmarkers==0),]
locations[,2]= -locations[,2]
space_time_locations = cbind(rep(timelocations,each=nrow(locations)),rep(locations[,1],NumTimeInstants),rep(locations[,2],NumTimeInstants))
ExactData = f(space_time_locations[,2],space_time_locations[,3],space_time_locations[,1])

# create the locations for the evaluation of the errors
mesh = create.mesh.2D(nodes = horseshoe2D$boundary_nodes, segments = horseshoe2D$boundary_segments) 
mesh = refine.mesh.2D(mesh=mesh,maximum_area = 0.05)
mesh = refine.mesh.2D(mesh=mesh,minimum_angle = 30)
meshr = refine.mesh.2D(mesh,30,0.015)
evaluatePoints_space = meshr$nodes[which(meshr$nodesmarkers==0),]
evaluatePoints_time = seq(0,tf,length.out = 41)
evaluatePoints = cbind(rep(evaluatePoints_time,each=nrow(evaluatePoints_space)),rep(evaluatePoints_space[,1],length(evaluatePoints_time)),rep(evaluatePoints_space[,2],length(evaluatePoints_time)))
exact_solution = f(evaluatePoints[,2],evaluatePoints[,3],evaluatePoints[,1])
RMSE<-function(s1,s2) sqrt(mean((s1-s2)^2))

# create structure for the censoring of the data
mesh_ref = refine.mesh.2D(mesh,minimum_angle = 30,maximum_area = 0.00001)
mesh_ref_time = seq(0,tf,length.out = 100)

# Objects for ST-PDE ####
# Time and space discretization 
TimePoints = seq(0,tf,length.out =21) 
FEMbasis = create.FEM.basis(mesh)

# Range of lambdaS and lambdaT for the different censoring schemes
if(schema=="d" && p==0.5) {
  lambdaS_Sep=10^seq(-2, 0 ,length.out=5)
  lambdaT_Sep=10^seq(-5,-3,length.out=5)
}else if(schema=="b" && p==0.5) { 
  lambdaS_Sep=10^seq(-1, -0.5 ,length.out=5)
  lambdaT_Sep=10^seq(-5.5,-4.5,length.out=5)
}else if(schema=="c" && p==0.5) { 
  lambdaS_Sep=10^seq(-1.1, -0.6 ,length.out=5)
  lambdaT_Sep=10^seq(-5.5,-4.5,length.out=5)
}else {
  lambdaS_Sep=10^seq(-1, -0.5 ,length.out=5)
  lambdaT_Sep=10^seq(-5,-3,length.out=10)[5:10]}

# Objects for TPS ####
TimeData = rep(timelocations,each=nrow(locations))
xData = rep(locations[,1],length(timelocations))
yData = rep(locations[,2],length(timelocations))
knots = data.frame(xData=locations[,1],yData=locations[,2])
NewPoints = data.frame(xData=evaluatePoints[,2],yData=evaluatePoints[,3],TimeData=evaluatePoints[,1])

# Objects for SOAP ####
bound=mesh$nodes[mesh$nodesmarkers == 1,]
bound = rbind(bound, bound[1,])
fsb=list(list(bound[,1],bound[,2]))
names(fsb[[1]]) <- c("xData","yData")

# Objects for KRIG ####
grd <- SpatialPoints(evaluatePoints_space)
tgrd <- as.POSIXct("2010-08-05")+3600*seq(0,tf,length.out = 41)
pred.grid_C <- STF(sp=grd, time=tgrd)
timeinstants <- as.POSIXct("2010-08-05")+3600*timelocations

# Objects for INLA ####
#Boundary
c.bdy = data.frame(x = rev(c(horseshoe2D$boundary_nodes[,1],horseshoe2D$boundary_nodes[1,1])), 
                   y = rev(c(horseshoe2D$boundary_nodes[,2],horseshoe2D$boundary_nodes[1,2])))
# segm.bnd <- inla.mesh.segment(c.bdy)
pol = Polygon(c.bdy)
pol = Polygons(list(pol), ID = "none")
segm.bnd = SpatialPolygons(list(pol))

#Define mesh
max.edge = 0.375
bound.outer = 1
mesh <- inla.mesh.2d(loc.domain = cbind(c.bdy$x, c.bdy$y), 
                     loc = locations, offset = c(max.edge, bound.outer), cutoff = 0.04,
                     boundary = segm.bnd,
                     max.edge = c(1,5)*max.edge, min.angle = 30)
#Save num locations and time instants
n_locations <- nrow(locations)
n_data <- length(ExactData)
n_time <- as.integer(n_data/n_locations)
#Data, to be update in the for loop
c = data.frame(x = space_time_locations[,2],
               y = space_time_locations[,3],
               time = 1+space_time_locations[,1]*40,
               data = ExactData)
coordinates(c) <- ~x+y
#INLA objects independent on the data
A.c =
  inla.spde.make.A(mesh = mesh,
                   loc = coordinates(c),
                   group = c$time,
                   n.group = n_time
  )

#Barrier models
tl = length(mesh$graph$tv[,1])
# - the number of triangles in the mesh
posTri = matrix(0, tl, 2)
for (t in 1:tl){
  temp = mesh$loc[mesh$graph$tv[t, ], ]
  posTri[t,] = colMeans(temp)[c(1,2)] 
}
posTri = SpatialPoints(posTri)
# - the positions of the triangle centres
normal = over(segm.bnd, SpatialPoints(posTri), returnList=T)
# - checking which mesh triangles are inside the normal area
barrier.tri = setdiff(1:tl, unlist(normal))
# - the triangles inside the barrier area
prior.range = c(1, .5)
prior.sigma = c(3, 0.01)
poly.barrier = inla.barrier.polygon(mesh, barrier.triangles = barrier.tri)
barrier.model = inla.barrier.pcmatern(mesh, barrier.triangles = barrier.tri, prior.range = prior.range, prior.sigma = prior.sigma)
#Create data structure
s.index <- inla.spde.make.index(name = "s",
                                n.spde = barrier.model$f$n,
                                n.group = n_time)
hyper.iid = list(prec = list(prior = 'pc.prec', param = prior.sigma)) 

RMSEstpde=RMSEtps=RMSEsoap=RMSEkriging=RMSEinla=rep(0,numIter)
for(iter in beginIter:numIter)
{
  set.seed(254+iter-1)
  Data = ExactData+ rnorm(length(ExactData), mean = 0, sd =  0.05*diff(range(ExactData)))
  DataMat = matrix(data=Data,nrow=nrow(locations),ncol=length(timelocations))
  set.seed(254+iter-1)
  DataMat=createNA(DataMat,schema=schema,p=p,
                   mesh_ref=mesh_ref,mesh_ref_time=mesh_ref_time,locations = locations,
                   timelocations=timelocations,RDD_groups = nRDD)
  
  # ST-PDE
  solutionSTPDE = smooth.FEM.time(locations = locations,
                                  time_mesh = TimePoints,
                                  time_locations = timelocations,
                                  observations = DataMat,
                                  FEMbasis = FEMbasis,
                                  lambdaS = lambdaS_Sep,
                                  lambdaT = lambdaT_Sep,
                                  lambda.selection.lossfunction	= "GCV",
                                  DOF.evaluation = "stochastic")
  
  sol_eval = eval.FEM.time(solutionSTPDE$fit.FEM.time,space.time.locations = evaluatePoints,lambdaS=solutionSTPDE$bestlambda[1],lambdaT=solutionSTPDE$bestlambda[2])
  RMSEstpde[iter] = RMSE(sol_eval,exact_solution)
  
  # TPS
  solutionTPS = gam(as.vector(DataMat)~
                 te(xData,yData,TimeData,k=c(NA,NA),d=c(2,1),bs=c("tp","cr")), 
                 method="GCV.Cp",knots=knots)
  PredictionTPS = as.vector(predict(solutionTPS,NewPoints))
  RMSEtps[iter] = RMSE(PredictionTPS,exact_solution)
  
  # SOAP
  solutionSOAP <-gam(as.vector(DataMat)~
                       te(xData,yData,TimeData,k=c(NA,NA),d=c(2,1),bs=c("sf","cr"),
                          xt=list(bnd=fsb)),method="GCV.Cp",knots=knots)
  PredictionSOAP<-as.vector(predict(solutionSOAP,NewPoints))
  RMSEsoap[iter]<-RMSE(PredictionSOAP,exact_solution)
  
  # KRIG
  STFDF.data <- stConstruct(t(DataMat),
                            space = list(values = 1:nrow(DataMat)),
                            time = timeinstants,
                            SpatialObj = SpatialPoints(locations),
                            interval = TRUE)
  STFDF.data <- as(STFDF.data, "STSDF")
  vv = variogramST(values ~ 1, STFDF.data, na.omit = T)  # the spatio-temporal sample variogram
  Vgm_model <- vgmST("separable", space = vgm(0.002, "Sph", 1.5, 10^(-3)), time = vgm(0.002, "Sph", 1.5*40, 10^(-3)), sill=0.002)  # the desired spatio-temporal model (separable exponential covariance model)
  v <- fit.StVariogram(vv, Vgm_model, method = "L-BFGS-B") # fit a spatio-temporal variogram of the given type (separable exponential covariance model) to spatio-temporal sample variogram (vv).
  
  PredictionKRIG <- krigeST(values ~ 1, STFDF.data, pred.grid_C, v)[[1]]
  RMSEkriging[iter] <- RMSE(PredictionKRIG, exact_solution)
  
  # INLA
  Data = as.vector(DataMat)
  
  #Create data structure
  c.stack <- inla.stack(data  = list(data = Data),
                        A = list(A.c),
                        effects = list(c(s.index, list(Intercept = 1))),
                        tag = "c.data")
  #Fit model
  form <- data ~ -1 + Intercept + f(s, model = barrier.model, 
                                    group=s.group, 
                                    control.group=list(model="ar1"))
  
  solutionINLA <- inla(form, data = inla.stack.data(c.stack),
             family = "gaussian",
             control.predictor = list(A = inla.stack.A(c.stack), compute = TRUE),
             control.compute = list(cpo = TRUE))
  field = solutionINLA$summary.random$s$mean
  #Get predicted data on grid
  A.grid = inla.spde.make.A(mesh, loc=evaluatePoints_space)
  PredictionINLA <- NULL
  for (ind in unique(c$time)) {
    field_ind = field[(1+(ind-1)*mesh$n):(ind*mesh$n)]
    PredictionINLA = c(PredictionINLA, drop(A.grid %*% field_ind) + solutionINLA$summary.fixed$mean)
  }
  RMSEinla[iter] = RMSE(PredictionINLA,exact_solution)
}
