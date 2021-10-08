# Supplementary material for:
# Smoothing spatio-temporal data withcomplex missing data patterns
# by Eleonora Arnone, Laura M. Sangalli, Andrea Vicini
# 
# Section 5.5

rm(list=ls())
library(fdaPDE)
# version needed: 1.1-2 or later
# most updated library version always available at https://github.com/fdaPDE/fdaPDE
# periodically released on CRAN: https://cran.r-project.org/web/packages/fdaPDE/index.html
library(mgcv); library(sp); library(spacetime); library(gstat); library(INLA); library(sinkr);

source("censoring.R")

##################   Preparation   ##################   
schema="a" 
# schema="b"; nRDD=8
# schema="c"; nRDD=12
# schema="d"; nRDD=40
tf=1; p=0.5; numIter=30; beginIter=1

# function to generate data
coe = function(x,y) 1/2*sin(5*pi*x)*exp(-x^2)+1
f = function(x,y,t) 
{
  sin(2*pi*(coe(y,1)*x*cos(9/5*t/tf-2)-y*sin(9/5*t/tf-2)))*cos(2*pi*(coe(y,1)*x*cos(9/5*t/tf-2+pi/2)+coe(x,1)*y*sin((9/5*t/tf-2)*pi/2)))
}

# create the vector of complete time instants and the matrix of complete observations
NumTimeInstants = 41
timelocations = seq(0, tf, length.out = NumTimeInstants)

k=15
eps=0.05
boundary_nodes=unique(rbind(cbind(seq(0,1,length.out = 17),0),cbind(seq(0,1,length.out = 17),1),cbind(0,seq(0,1,length.out = 17)),cbind(1,seq(0,1,length.out = 17))))
locations=cbind(rep(seq(0+eps,1-eps,length.out = k),k),rep(seq(0+eps,1-eps,length.out = k),each=k))
meshr=create.mesh.2D(nodes=rbind(boundary_nodes,locations))
space_time_locations = cbind(rep(timelocations,each=nrow(locations)),rep(locations[,1],NumTimeInstants),rep(locations[,2],NumTimeInstants))
ExactData = f(space_time_locations[,2],space_time_locations[,3],space_time_locations[,1])

# create the locations for the evaluation of the errors
evaluatePoints=space_time_locations
exact_solution=ExactData
RMSE<-function(s1,s2) sqrt(mean((s1-s2)^2))

# create structure for the censoring of the data
mesh_ref_time=seq(0,1*tf,length.out = 100)
mesh_ref=mesh

# Objects for ST-PDE ####
# Time and space discretization 
TimePoints=seq(0,tf,length.out =11)
FEMbasis = create.FEM.basis(meshr)

# Range of lambdaS and lambdaT 
lambdaS_Sep=10^seq(-5.75, -5.25 ,length.out=5)
lambdaT_Sep=10^seq(-5.5, -4.5,length.out=5)

# Objects for TPS and SOAP ####
j = 9
nmax = 100
knots <- data.frame(xData=rep(seq(0+2*eps,1-2*eps,length.out = j),j),
                    yData=rep(seq(0+2*eps,1-2*eps,length.out = j),each=j))

TimeData=rep(timelocations,each=nrow(locations))
xData=rep(locations[,1],length(timelocations))
yData=rep(locations[,2],length(timelocations))
NewPoints<-data.frame(xData=evaluatePoints[,2],yData=evaluatePoints[,3],TimeData=evaluatePoints[,1])
bound=rbind(c(0,0),c(1,0),c(1,1),c(0,1),c(0,0))
fsb=list(list(bound[,1],bound[,2]))
names(fsb[[1]]) <- c("xData","yData")

# Objects for KRIG ####
grd <- SpatialPoints(locations)
tgrd <- as.POSIXct("2010-08-05")+3600*seq(0,tf,length.out = 41)
pred.grid_C <- STF(sp=grd, time=tgrd)
timeinstants <- as.POSIXct("2010-08-05")+3600*timelocations

# Objects for INLA ####
#Boundary
square.bdy = data.frame(x = c(0,1,1,0,0),
                          y = c(0,0,1,1,0))

#Define mesh with border
mesh <- inla.mesh.2d(loc.domain = cbind(square.bdy$x, square.bdy$y), 
                     loc = cbind(rep(seq(0+eps,1-eps,length.out = k-2),k-2),rep(seq(0+eps,1-eps,length.out = k-2),each=k-2)), 
                     max.edge = 0.115, min.angle = 30)

#Save num locations and time instants
n_locations <- nrow(locations)
n_data <- length(ExactData)
n_time <- as.integer(n_data/n_locations)

#Create SPDE
square.spde <- inla.spde2.matern(mesh = mesh, alpha = 2)

#Data, to be update in the for loop
square = data.frame(x = space_time_locations[,2],
                      y = space_time_locations[,3],
                      time = 1+space_time_locations[,1]*40,
                      data = ExactData)
coordinates(square) <- ~x+y

#INLA objects independent on the data
A.square =
  inla.spde.make.A(mesh = mesh,
                   loc = coordinates(square),
                   group = square$time,
                   n.group = n_time
  )

s.index <- inla.spde.make.index(name = "spatial.field",
                                n.spde = square.spde$n.spde,
                                n.group = n_time)

#Create data structure for prediction
repeated.mesh.nodes = cbind(rep(timelocations,each=mesh$n),rep(mesh$loc[,1],NumTimeInstants),rep(mesh$loc[,2],NumTimeInstants))
repeated.mesh.nodes[,1] = 1 + repeated.mesh.nodes[,1]*40

A.pred <- inla.spde.make.A(mesh = mesh,
                           loc = repeated.mesh.nodes[,2:3],
                           group = repeated.mesh.nodes[,1],
                           n.group = n_time
)
square.stack.pred <- inla.stack(data = list(data = NA),
                                  A = list(A.pred),
                                  effects = list(c(s.index, list (Intercept = 1))),
                                  tag = "square.pred")

RMSEstpde=RMSEtps=RMSEsoap=RMSEkriging=RMSEdineof=RMSEinla=rep(0,numIter)

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
  RMSEstpde[iter]=RMSE(sol_eval,exact_solution)
  
  # TPS
  solutionTPS <-gam(as.vector(DataMat)~
                 te(xData,yData,TimeData,k=c(120,10),d=c(2,1),bs=c("tp","cr")),
               method="GCV.Cp",knots=knots)
  PredictionTPS = as.vector(predict(solutionTPS,NewPoints))
  RMSEtps[iter]<-RMSE(PredictionTPS,exact_solution)
  
  # SOAP
  solutionSOAP = gam(as.vector(DataMat) ~ te(xData,yData,TimeData,bs=c("sf","cr"),k=c(25,10),d=c(2,1),
                                    xt=list(list(bnd=fsb,nmax=nmax),NULL))+
              te(xData,yData,TimeData,bs=c("sw","cr"),k=c(25,10),d=c(2,1),
                 xt=list(list(bnd=fsb,nmax=nmax),NULL)),knots=knots)
  
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
  square.stack <- inla.stack(data  = list(data = Data),
                               A = list(A.square),
                               effects = list(c(s.index, list(Intercept = 1))),
                               tag = "square.data")
  #Join stack
  join.stack <- inla.stack(square.stack, square.stack.pred)
  
  #Fit model
  form <- data ~ -1 + Intercept + f(spatial.field, model = spde, 
                                    group=spatial.field.group, 
                                    control.group=list(model="ar1"))
  solutionINLA <- inla(form, data = inla.stack.data(join.stack, spde = square.spde),
             family = "gaussian",
             control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
             control.compute = list(cpo = TRUE))
  
  #Get predicted data on grid
  index.data <- inla.stack.index(join.stack, "square.data")$data
  
  PredictionINLA <- solutionINLA$summary.fitted.values[index.data, "mean"]
  RMSEinla[iter] = RMSE(PredictionINLA,exact_solution)
  
  # DINEOF
  solutionDINEOF   = sinkr::dineof(Xo=DataMat)
  PredictionDINEOF = as.vector(solutionDINEOF$Xa)
  RMSEdineof[iter] = RMSE(PredictionDINEOF,exact_solution)
}
