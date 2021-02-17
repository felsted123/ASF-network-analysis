
# Load EpiModel
library(EpiModel)

# Initialize the network
nw <- network.initialize(n = 75, directed = T) 

# Set up a farmtype attribute on the network;
# This function is the long-form version of the %v% function
nw <- set.vertex.attribute(nw, attrname = "farmtype", rep(1:3, c(5,25,70)))

farm<-get.vertex.attribute(nw, "farmtype")
farm
#nw <- set.vertex.attribute(nw, attrname = "farmtype", rbinom(100,1,0.5)

# Parameterizing the model
formation <- ~edges + nodefactor("farmtype") + nodematch("farmtype")

# Target stats for model 1 (the null model)
target.stats <- c(150, 37.5, 70, 63) #the first number refers to the overall degree for the network, which we assumed to be 2, 
#second number is mean degree for farmtype2, lets assume that is 1.5 (25*1.5), and third number is mean degree for farmtype3, which lets assume to be 1 (70*1)
#and the last number represents the proportion of contacts that are within the same farm type, say 90% (70*.9)

# Parameterizing the dissolution model (will be the same for model 2)
coef.diss <- dissolution_coefs(dissolution = ~offset(edges),
                               duration = 40, d.rate = 0.01) # if we dont want the departure, then we can remove the d.rate part, in this code the departure is 
#similar across all node types, however, if we want the departure 
coef.diss

# Fit the model 1
est1 <- netest(nw, formation, target.stats, coef.diss)

# Examine the model coefficients
summary(est1)

# Model diagnostics. We can input other terms to diagnose using nwstats.formula,
# and here we only change the reference category of the nodefactor term (levels = NULL
# means that there is no reference category, so we can see the stats for both farmtype
# 0 and farmtype 1 groups)
dx1 <- netdx(est1, nsims = 5, nsteps = 100, ncores = 4,
             nwstats.formula = ~edges +
                                nodefactor("farmtype", levels = NULL) +
                                nodematch("farmtype"))
dx1
plot(dx1)

# Target stats for model 2: higher mean degree for farmtype 1 and high assortative mixing

inf<-runif(100, 0.1,0.2)# here we are using a uniform distribution 
contact=c(0.1, 0.2, 0.6)# this reflects the weekly contact rate for each of the three farm types, please use different values to 
#make sure that it is reflected in disease spread, in same manner other parameters could also be specified

# Parameterizing the model
param <- param.net(inf.prob = 0.9, act.rate = contact ,
                   a.rate = 0.01, ds.rate = 0.01, di.rate = 0.01)

#param <- param.net(inf.prob = 0.7, inf.prob.g2 = 0.5, act.rate = 0.5, inf.prob.g3=0.9, balance = "g1",
                 #  rec.rate = 1/25, rec.rate.g2 = 1/50, a.rate = 1/100, a.rate.g2 = NA,
                  # ds.rate = 1/100, ds.rate.g2 = 1/100, di.rate = 1/90, di.rate.g2 = 1/90)

# Initial conditions
init <- init.net(farm[25]) #we should be able to choose the farm we want to be index farm, here the 
#second farm which is a farmtype 1 farm is chosen

#the alternate way is to randomly sample using status vector as shown below
#n=75
#status <- sample(c("s", "i"), size = n, replace = TRUE, prob = c(0.99, 0.01))
#init.net(status.vector = status)

# Control settings (reduced number of simulations for computational efficiency)
control <- control.net(type = "SI", nsteps = 10, nsims = 1,
                       epi.by = "farmtype", delete.nodes = TRUE)

# Simulate the two counterfactual models
sim1 <- netsim(est1, param, init, control)


# Examine the model results
sim1

# Post-simulation diagnostics
par(mfrow = c(1, 2))
plot(sim1, type = "formation", stats = "edges")

# Model prevalence overall, using the add argument to plot one model on top of
# the other
plot(sim1, y = "i.num", qnts = 1, mean.col = "steelblue",
     qnts.col = "steelblue", main = "Total Prevalence")
legend("topleft", c("Model 1"), lwd = 3,
       col = c("steelblue", "firebrick"), bty = "n")

# Model results stratified by farmtype
par(mfrow = c(1, 2))
plot(sim1, y = c("i.num.farmtype3", "i.num.farmtype1", "i.num.farmtype2"),  legend = TRUE, qnts = 1,
     ylim = c(0, 200), main = "M1: Disease Prevalence by farmtype")
