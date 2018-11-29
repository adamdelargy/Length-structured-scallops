#### Adam Delargy

### Full script 
# - negative loglikelihood, length dependent sigma 
# - alternative survey simulation 


# compiled 24/11/2018


##### Load data 

setwd("M:/My Documents/Stock assessment/length_assess/scallop length assess/Fished area/My modifications")

#### data 
	
Ns = as.matrix(read.csv("Obs Nsurv.csv", header=F))
dimnames(Ns)   <-  list("length"=1:nrow(Ns)+1, "year" = 2012:2016) 

Ns[,4] <- NA # add in the missing data for 2015 

### visualise data 

# let's look at the survey data
par(mfrow=c(1,1))
plot(Ns[,1]~ rownames(Ns), type="l", xlim=c(min(as.numeric(rownames(Ns))),max(as.numeric(rownames(Ns)))), ylim=c(0,max(Ns,na.rm=T)), xlab="Length (cm)")
for(i in 2:dim(Ns)[2]){
	lines(Ns[,i]~ rownames(Ns), col=i)
}
legend("topright", legend=2012:2016, lty=1,col=1:5)

# we can't have zeros in the model
# so we will start the model from 6cm and we will add the 15cm in to the 14cm class 


Ns <- Ns[-c(1:4),]


par(mfrow=c(1,1))
plot(Ns[,1]~rownames(Ns), type="l", xlim=c(min(as.numeric(rownames(Ns))),max(as.numeric(rownames(Ns)))), ylim=c(0,max(Ns,na.rm=T)), xlab="Length (cm)")
for(i in 2:dim(Ns)[2]){
	lines(Ns[,i]~rownames(Ns), col=i)
}
legend("topright", legend=2012:2016, lty=1,col=1:5)


obs<-list(
	Lb = scan("Obs Lb.csv"),
	Ns = Ns)

Eff <- scan("Obs Effort.csv")


# the catch data here is aggregated
par(mfrow=c(1,1))
plot(obs$Lb,type="l")


library(mgcv); library(Matrix); library(lattice)



####


#############################################################################

# FUNCTIONS 


# Author:
#
#   Dr Douglas C. Speirs
#   Department of Mathematics & Statistics
#   University of Strathclyde
#	

non.negative <- function (x){
	# utlity function to truncate a vector at zero from below
	x[x<0]<-0
	return(x)}
	
pgam <- function (x1,x2,a,b) {
   # Function returning the probability over interval (x1,x2) from a Gamma distribution with +ve parameters "a" and "b"
   # mean = a*b, variance = a*b^2
   pgamma(x2,shape=a,scale=b)-pgamma(x1,shape=a,scale=b)}

pdL <- function (L,L0,Linf,k,beta,dL,dt){
	# Probability of growing from length L0 to L over timestep dt...
	MeanDeltaL<-non.negative((Linf-L0)*(1-exp(-k*dt)))
	alphaL<-MeanDeltaL/beta
    DeltaL<-L-L0
	return(pgam(DeltaL,DeltaL+dL,alphaL,beta))}

logistic <- function (L,alpha,beta){
	# logistic function of length used in gear selectivity and discarding fraction
	return(1/(1+alpha*exp(-beta*L)))}



survey <- function (N,L,pA=NULL,alpha,beta,qmax){
   # Description:
   #
   # Samples a population maxtrix of number-at-length by year given   
   # the survey swept area each year, and logistic survey gear selectivity
   # to yield a matrix of survey number-at-length by year.
   #
   #
   # Usage:
   #
   #     survey(N,L,Pa,qmax,alpha,beta)
   #
   #
   # Arguments:
   #               N: an n x maxiter array where N[,j] is the population length distribution
   #                  in timestep j, n is the number of length classes, and maxiter the number of timesteps
   #              pA: a vector of length maxiter containing the proportion of area swept by 
   #                  the survey each timestep (dimensionless)
   #               L: a vector of length n containing the lengths (lower bound) of each length class (cm)
   #            qmax: maximum survey selectivity (dimensionless, normally of order 1)
   #           alpha: survey selectivity parameter
   #            beta: survey selectivity parameter
   #
   # Value:
   #
   #        An n x maxiter array of surveyed number-at-length for each year
   #  
    if (is.null(pA)) {pA<-rep(1,dim(N)[2])}
	q <- qmax*logistic(L,alpha,beta)
	Survey.N <- q%o%pA*N
	return(Survey.N)}


InitLenDist <- function (P,Z,p,s,n,dt,Rhist) {
   #    Function to calculate an initial length distribution, 
   #    consistent with the growth model and a 'historical' 
   #    (i.e. prior to the main model run) recruitment.  The
   #    historical recruitment values required to acheive the 
   #    initial condition are fitting parameters. 
   #    Returns a vector of number at length
   nhist<-length(Rhist)
   N<-rep(0,n)
   S<-diag(exp(-Z*dt))
   for (i in 1:nhist){N<-P%*%S%*%N+p*Rhist[i]}
   return(N)}

MakeParList <- function(fitting.pars,fixed.pars){
    attach(fixed.pars)
    n1<-length(fitting.pars)-2*maxiter-nhist-4
    n2<-n1+maxiter
    n3<-n2+maxiter
    n4<-n3+nhist
    detach(fixed.pars)
    scalar.pars<-as.list(fitting.pars[1:n1])
    vector.pars<-list(Ft=unname(fitting.pars[(n1+1):n2]),
                      Rt=unname(fitting.pars[(n2+1):n3]),
                      Rhist=unname(fitting.pars[(n3+1):n4]),
				qest = unname(fitting.pars[(n4+1):length(fitting.pars)])
							)
    pars<-c(scalar.pars,vector.pars,fixed.pars)	
    return(as.list(pars))}





#len.sim<-function(params,L,dL,W,dt,maxiter,N0=NULL,Discard=FALSE,alphaD=NULL,betaD=NULL,Survey=FALSE,pA=NULL) {
len.sim<-function(par.list,N0=NULL) {
   # Description:
   #
   # Simulates a length-structured fish population using the population model used by 
   # Sullivan et al. (1990) in their catch-at-length stock assessment model. 
   #
   #
   # Usage:
   #
   #    len.sim (par.list)
   #
   #
   # Arguments:
   #
   #    par.list: a list containing the following elements:
   #             Linf: von Bertalanffy asympotic length (cm)
   #                k: von Bertalanffy growth rate (1/year)
   #             beta: gamma distribution scale parameter for growth increments
   #                M: natural mortality rate (1/year)
   #           alphaS: fishing mortality selectivity parameter
   #            betaS: fishing mortality selectivity parameter
   #           alphaR: gamma distribution shape parameter for recruitment length distribution
   #            betaR: gamma distribution scale parameter for recruitment length distribution
   #             qmax: survey maximum selectivity (dimensionless, normally of order 1, only used if Survey=T)
   #           alphaV: survey selectivity parameter (only used if Survey=T)
   #            betaV: survey selectivity parameter (only used if Survey=T)
   #               Ft: vector of time-dependent component to the fishing mortality (must be of length maxiter) (1/year)
   #               Rt: vector of time-dependent recruitment (must be of length maxiter-1)
   #            Rhist: vector of time-dependent recruitment used by function InitLenDist to generate 
   #                   the initial condition if that is not provided separately
   #
   #               dt: timestep (years)
   #          maxiter: number of timesteps
   #               dL: lengh-class width (cm)            
   #                L: a vector of length class sizes - calculations assume equal-width classes (cm)
   #                W: a vector of individual weight at length (kg)
   #           alphaD: discarding selectivity parameter (only used if Discard=T)
   #            betaD: discarding selectivity parameter (only used if Discard=T)
   #               pA: vector of length maxiter giving the proportion of population 
   #                   area swept out by the survey for each timestep (only used if Survey=T)
   #
   # Details:
   #
   #    The vector 'Ft' must be of length 'maxiter' i.e. need fishing mortalities for all years
   #
   #    The vector 'R't must be of length 'maxiter' i.e. need recruitment for all years after the initial year
   #
   #    The vector 'W' must have the same length as 'L'.
   #    
   #    Requires the user-defined function pgam, and depending on settings the functions logistic, survey, and InitLenDist
   #
   # Value:
   #
   #    A list with the following elements:
   #
   #           N: an n x maxiter array where N[,j] is the population number-at-length distribution
   #              in timestep j. 
   #
   #           C: an n x maxiter array where C[,j] is the catch number-at-length distribution
   #              in timestep j.
   #
   #           F: an n x maxiter array where F[,j] is the fishing mortality (per year) at length
   #              in timestep j.
   #         
   #          Ln: an n x maxiter array where Ln[,j] is the landed catch numbers at length
   #              in timestep j.
   #
   #          Lb: a vector of length maxiter containing the landed biomass (tonnes) by timestep
   #
   #          Ns: an n x maxiter array where Ns[,j] is the survey number-at-length distribution
   #              in timestep j. 
   #
   # References:
   # Sullivan, P.J., Lai, H.-L. & Gallucci, V.F. A catch-at-length analysis that incorporates a 
   #   stochastic model of growth. Canadian Journal of Fisheries and Aquatic Sciences 47, 184-198.
   #
   attach(par.list)
   n<-length(L)         # number of length classes
   empty<-matrix(0,nrow=n,ncol=maxiter)	# create matrix of length class by length class
   P<-outer(L,L,pdL,Linf,k,beta,dL,dt)	# here we work out the probability of every age class reaching Linf in a single timestep
   P[n,]<-1-apply(P[-n,],2,sum)		# Here we replace the largest age class and work out the chance (1 - the sum of chance of all other length classes)
   							# I think the zeros are caused because my Linf is lower than some of the size classes
   p<-pgam(L,L+dL,alphaR,betaR)		# this must be where we fit the recruitment-length parameters 
   s<-logistic(L,alphaS,betaS)		# here we fit logistic curve to estimate selectivity
   F<-s%o%Ft					# here we apply selectivity to our fishing mortality estimates to get fishing mortality at size
   Z<-F+M						# add fishing and natural mortality together 
   N<-empty
   N[,1]<-if(is.null(N0)) InitLenDist(P,Z[,1],p,s,n,dt,Rhist) else N0	# here we calculate an initial abundance 
					# this is described in the papers and uses a survivalship, recruitment history guess, probability of increasing to Linf in dt and recruiment-length parameters 
   for (i in 1:(maxiter-1)){		# now we deal with the other length classes 
      S<-diag(exp(-Z[,i]*dt))		# here this is survivalship estimation - simply taking the numbers and applying mortality 
      N[,i+1]<-P%*%S%*%N[,i]+p*Rt[i]}	# here we apply survivalship, recruitment estimation, probability of increasing to Linf in dt and recruiment-length parameters to obtain numbers in each length class
	
   C<-F/Z*(1-exp(-Z*dt))*N			# here we work out the catch from the numbers - Baranov catch equation
   Ln<- if (Discard==T) logistic(L,alphaD,betaD)*C else C	# here we allow discards to be removed using logistic curve if possible
										# otherwise our landings are just our catch
   if (Survey==T) {
		U1q <- c(qest, rep(qest[length(qest)],n - length(qest)))
		Ns <- U1q%*%pA*N*exp(-0.5*(Z))
		} 
	else { Ns <- empty}	# here we use a logistic equation to estimate survey densities - this is the bit I could improve 
   Lb<-apply(W*Ln,2,sum)*1e-3		# convert landing numbers to biomass 
   detach(par.list)
   return(list(N=N,C=C,F=F,Ln=Ln,Lb=Lb,Ns=Ns))}

#############################################################################


#############################################################################


# Create parameter lists 


# Step 1: set up known parameters and guesses for other parameters


maxiter<-5			# number of timesteps 

Lmin <- 6			# minimum length class
Lmax <- 14			# max length class
dL <- 1       	     # length class width (cm)


Lengths  <- seq(Lmin,Lmax,dL)     # vector of class lengths (cm) 

years <- seq(2012,2016,1)

# set seed to keep the random pattern fixed for now

set.seed(1234)

Fscalar <- c(1.104945e+00,2.774086e-01,3.544160e-01,2.139021e-01,1.317979e-01)
#Fscalar <- rep(0.5,maxiter)		# parameters describing Fa,t combined with the design matrix from GAM
#Rec<-1e5*runif(maxiter, 1, 5)		# guess initial value - we are defining recruitment as the new numbers added to the fishable biomass in a given year 
#RecHist<-1e5*runif(maxiter, 1, 5)		# guess
Rec <- c(5.766029e+07, 6.390263e+06,6.301679e+06,5.866751e+06,9.519651e-10)
RecHist <- c(5.892529e-01, 4.640805e+08,8.187937e-05,5.334389e+07,4.681234e+07)
pA<-c(8.409929e-06, 9.810846e-06, 9.836951e-06, NA, 1.920750e-05) 	# script in my modifications folder 
write.table(Rec,"Rtsim.csv",sep=",",row.names=F,col.names=F)
write.table(RecHist,"Rhsim.csv",sep=",",row.names=F,col.names=F)
write.table(Fscalar,"Ftsim.csv",sep=",",row.names=F,col.names=F)
write.table(pA,"pA.csv",sep=",",row.names=F,col.names=F)
rm(maxiter,pA,dL)


a<-0.0003945466                # coefficient for length-weight relationship (kg cm^-b)
b<-2.500605                	# power for length-weight relationship (dimensionless)




# Step 2: sort the parameters in to whether fixed or fitting 


fitting.pars<-c(
	alphaS = 2.034151e+10,
	betaS =  3.076253e+00,
	sigma = c( 2.145226e-05,2.093416e+00,1.072160e+00, 9.862280e-01),		
    Ft		= scan("Ftsim.csv"),
    Rt		= scan("Rtsim.csv"),
    Rhist	= RecHist<-scan("Rhsim.csv"),
	qest = c(2.083695e-01,2.484236e-01,4.105711e-01,6.383062e-01)
	)

fixed.pars<-list(
	dt		= 1,				        # Model timestep (years)
	maxiter	= 5,				        # number of timesteps        
	nhist   = 5,                        # length of recruitment history    
	dL	= 1,
	L	= seq(Lmin,Lmax,1),     
	W       = a*seq(Lmin,Lmax,1)^b, # vector of class weights (kg)
	M       = 0.6670458,                      # natural mortality (1/timestep)
	Survey  = TRUE,                     # 
	pA      = scan("pA.csv"),           # proportion of sea area swept out by survey each year (dimensionless)
	Discard = FALSE,			        # discarding selectivity parameter
	Linf	= 14.39295,				# vbgf parameters
	k		= 0.3422719,
	beta = 1,
	alphaR = 5.031054e+01,
	betaR = 1.972662e-01
	)


# Step 3: Run un-optimised calculations to generate some better starting values

par.list<-MakeParList(fitting.pars,fixed.pars)

x<-len.sim(par.list)
TSB<-apply(par.list$W*x$N,2,sum)*1e-3
plot(2012:2016,TSB,ylim=c(0,max(TSB)),xlab="year",ylab="TSB (1000's tonnes)",type="l")	


######################################################################################

# objective function
# weighted sum of squares 

fun<-function(ln.fitting.pars,fixed.pars,obs){
	pars<-exp(ln.fitting.pars)
	mod<-len.sim(MakeParList(pars,fixed.pars))
	 # calculate -log likelihood for two sources
	
	# length-dependent sigmas 
	sigmaU <- exp(log(pars[4])+log(pars[5])*Lengths+log(pars[6])*Lengths^2)


 	# calculate -log likelihood for two sources
  
  	nll_c <-  -sum(dnorm(log(mod$Lb),log(obs$Lb),pars[3],TRUE),na.rm=T)

  	nll_s <- 0 
  	for(y in 1:length(years)) {
    		    nll_s <- nll_s-sum(dnorm(log(mod$Ns[,y]),log(obs$Ns[,y]),sigmaU,TRUE),na.rm=T)
 	 }
  	return(nll_c + nll_s  )
}

#########################################################################################


# Begin assessment 


# here we allow parameters to vary by 30% 

# let's set seed for now

#set.seed(1234)
             
#CV<-0.3            
#ln.fitting.pars<-log(fitting.pars*(1+CV*rnorm(length(fitting.pars))))


ln.fitting.pars<-log(fitting.pars)		# here we take the logarithm of the fitting parameters
# as it is easier to find minima on a logarithm scale 


pred0<-len.sim(MakeParList(exp(ln.fitting.pars),fixed.pars))	# make another round of un-optimised predictions 

fit0<-fun(ln.fitting.pars,fixed.pars,obs)				# then use objective function without optim


##########################################
 # begin using optim  
# we will run all the different types of routine optim offers and quickly decide
# on best fit 

# Default - Nelder-Mead      
fit<-optim(ln.fitting.pars,fun,gr=NULL,fixed.pars,obs,control=list(maxit=30000,abstol=1e-15), hessian=T)  
          

# so anything included in ln.fitting.pars will be optimised to find 
# the least sum of squares produced in the 'fun' function above
# the fun function itself calls len.sim, which calls consequent functions to make all the model calculations

#what are the Hessian values?
# Remember they have to be all positive
eigen(fit$hessian)$values


pred<-len.sim(MakeParList(exp(fit$par),fixed.pars))


   
par(mfrow=c(2,3))
plot(obs$Lb, ylim=c(min(obs$Lb,pred0$Lb,pred$Lb,na.rm=T),max(obs$Lb,pred0$Lb,pred$Lb,na.rm=T)))
lines(pred0$Lb,lty=2)  
lines(pred$Lb,col="red",lwd=2) 

for (i in 1:5){
plot(fixed.pars$L,obs$Ns[,i],ylim=c(min(obs$Ns,pred0$Ns,pred$Ns, na.rm=T),max(obs$Ns,pred0$Ns,pred$Ns, na.rm=T)))
lines(fixed.pars$L,pred0$Ns[,i],lty=2)         
lines(fixed.pars$L,pred$Ns[,i],col="red",lwd=2)
abline(v=11,lty=2)
}


# pretty damn good 


par(mfrow=c(2,3))

for (i in 1:5){
plot(fixed.pars$L,obs$Ns[,i],ylim=c(min(obs$Ns,pred0$Ns,pred$Ns, na.rm=T),max(obs$Ns,pred0$Ns,pred$Ns, na.rm=T)), ylab="Numbers caught",xlab="Length (cm)", cex.lab=1.25,cex.axis=1.25, main=paste("Year",2011+i))
lines(fixed.pars$L,pred0$Ns[,i],lty=2)         
lines(fixed.pars$L,pred$Ns[,i],col="red",lwd=2)
abline(v=11,lty=2)
}

#### BFGS fit


fit2<-optim(ln.fitting.pars,fun,gr=NULL,fixed.pars,obs, method="BFGS",control=list(trace=T,maxit=30000),hessian=T)            
pred2<-len.sim(MakeParList(exp(fit2$par),fixed.pars))

par(mfrow=c(2,3))
plot(obs$Lb,ylim=c(min(obs$Lb,pred0$Lb,pred2$Lb,na.rm=T),max(obs$Lb,pred0$Lb,pred2$Lb,na.rm=T)))
lines(pred0$Lb,lty=2)  
lines(pred2$Lb,col="red",lwd=2) 

for (i in 1:5){
plot(fixed.pars$L,obs$Ns[,i],ylim=c(min(obs$Ns,pred0$Ns,pred2$Ns, na.rm=T),max(obs$Ns,pred0$Ns,pred2$Ns, na.rm=T)))
lines(fixed.pars$L,pred0$Ns[,i],lty=2)         
lines(fixed.pars$L,pred2$Ns[,i],col="red",lwd=2)
abline(v=11,lty=2)
}

# very nice




##################################################################

#####################################
# Checking results 

#plot logistic curve assumed for selectivity in catch at length
par(mfrow=c(1,1))
plot(logistic(Lengths, exp(fit$par[1]),exp(fit$par[2]))~Lengths, type="b")

# plot the mean Fbar where most catches are 
Fbar  <- apply(pred$F[2:7,],2, mean )
plot(x=years, y=Fbar, xlab="Year",type="b")


TSB<-apply(par.list$W*pred$N,2,sum)*1e-3
Catch <- apply(par.list$W*pred$C,2,sum)*1e-3
Abundance <- apply(pred$'N',2,sum)*1e-3
Fish.mort <- apply(pred$F,2,mean)

par(mfrow=c(2,2))
plot(2012:2016,TSB,ylim=c(0,max(TSB)),xlab="year",ylab="TSB (Tonnes)",type="l")	
plot(2012:2016,Catch,ylim=c(0,max(Catch)),xlab="year",ylab="Catch (Tonnes)",type="l")	
plot(2012:2016,Abundance,ylim=c(0,max(Abundance)),xlab="year",ylab="Abundance (1000's of scallops)",type="l")	
plot(2012:2016,Fish.mort,ylim=c(0,max(Fish.mort)),xlab="year",ylab="Fishing mortality",type="l")	


par(mfrow=c(1,1))
plot(2012:2016,Catch,ylim=c(0,max(Catch,obs$Lb)),xlab="Year",ylab="Catch (Tonnes)",type="b",col=2, cex.lab=1.5, cex.axis=1.5,lwd=2)
lines(2012:2016, obs$Lb,type='l',lwd=2)
legend("topright", legend=c("Observed landings", "Modelled landings"), lty=1,col=1:2, cex=1.5)

plot(2012:2016,Abundance,ylim=c(0,max(Abundance)),xlab="Year",ylab="Abundance (1000's of scallops)",type="b", cex.axis=1.5,cex.lab=1.5,lwd=2)

par(mfrow=c(1,2))
plot(2012:2016,Fish.mort,ylim=c(0,max(Fish.mort)),xlab="Year",ylab="Fishing mortality",type="b",cex.axis=1.5,cex.lab=1.5,lwd=2)	
plot(2012:2016,Eff,ylim=c(0,max(Eff)),xlab="Year",ylab="Effective Effort (kwdays)",type="b",cex.axis=1.5,cex.lab=1.5,lwd=2)

#### density 

# area of Cardigan Bay is total area between 3 and 12nm minus closed area 

Closedarea <- 255.383 + 423.761 - 131.55
# that's in km2 

# below is the total area of our Cardigan Bay fishing ground from 3 - 12nm
Totalarea <- 3773.985 # in km2

Openarea <- Totalarea - Closedarea

Sarea.m <- Openarea*1000000

Sarea.dm <- Sarea.m/100
# that's in square decameters 

par(mfrow=c(1,1))
plot(2012:2016,(Abundance*1000)/Sarea.dm,ylim=c(0,max((Abundance*1000)/Sarea.dm)),xlab="Year",ylab="Density (Count/100m2)",type="b",col=2, cex.axis=1.5,cex.lab=1.5,lwd=2)

survey.sa <- c(271.3372,316.5363,317.3785,NA,619.7090)
Nsum <- apply(obs$Ns,2,sum)
Ndens <- Nsum/survey.sa
lines(2012:2016,Ndens,type='b',lwd=2)
legend("topright", legend=c("Survey","Model abundance"),lty=1,col=1:2)


#plot the recruitment
recruitment  <- pred$'N'[1,]
plot(x=years, y=recruitment, xlab="Year",type="b")

# survey catchability 

q <- exp(c(fit$par[(length(fit$par)-3):length(fit$par)], rep(fit$par[length(fit$par)], 5)))
plot(q~Lengths)
# hmmmm why isn't this working?

# plot aggregated catches
csum <- obs$Lb
chatsum <- Catch
ymin <- min(cbind(csum,chatsum))
ymax <- max(cbind(csum,chatsum))
par(mfrow=c(2,2))
plot(years,csum,ylim=c(ymin,ymax),type='p')
lines(years,chatsum,lty=1)
cresid <- log(csum/chatsum)
qqnorm(cresid)
qqline(cresid)
acf(cresid)

# plot aggregated survey
tuningsum <- apply(obs$Ns,2,sum)
U1sum <- apply(pred$Ns,2,sum)
ymin <- min(cbind(tuningsum,U1sum),na.rm=T)
ymax <- max(cbind(tuningsum,U1sum),na.rm=T)
par(mfrow=c(2,2))
plot(years,tuningsum,ylim=c(ymin,ymax),type='p')
lines(years,U1sum,type='p',col=2)
residsum <- log(tuningsum/U1sum)
qqnorm(residsum)
qqline(residsum)


par(mfrow=c(1,1))
plot(years,tuningsum,ylim=c(ymin,ymax),type='p', cex.lab=1.5, cex.axis=1.5,
	pch=16,cex=1.5,xlab="Year", ylab="Survey Catch (numbers)")
lines(years,U1sum,type='p',col=2, pch=16,cex=1.5)
legend("topright", legend=c("Observed","Modelled"),pch=16,cex=1.5,col=1:2)



