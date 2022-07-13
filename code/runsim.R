# Load packages
library(coda)
library(rjags)

print("finished loading libraries")

# source files
source("./data_generation_funcs.R")

args <- commandArgs(trailingOnly = TRUE)

task_id <- strtoi(args[1])

params_matrix <- read.csv("./parameter_matrix.csv")

print("finished loading parameter matrix")

params <- params_matrix[1, ]

# generate data

# number of sites
nsites <- as.integer(params["nsites"])

# number of individuals in capture history
n_total <- as.integer(params["N"])

# number of sampling occasions in capture-recapture
nsample_cap <- as.integer(params["nsample_cap"])

# number of sampling occasions in presence-absence
nsample_pa <- as.integer(params["nsample_pa"])

# amount of movement for the species
tau <- as.double(params["tau"])

# number of camera traps
ntraps <- as.integer(params["ntraps"])

# setting seed
set.seed(as.integer(params["seed"]))

# homeranges for the n individuals
homeranges <- create_homeranges(n_total)

capture_hist <- generate_capture_hist(ntraps, nsample_cap, tau, homeranges)
capture_hist

presence_absence <- generate_presence_absence(nsites, nsample_pa, tau, homeranges)
presence_absence

output_data(presence_absence, capture_hist)

# Supplementary Code by Blanc et al. (2014)

#---------------------- patch occupancy data ----------------------------#
# Presence-absence data structure: rows = sites; columns = occasions
y <- read.csv("./presence-absence-data.csv", header = FALSE) # load presence-absence data
nsites <- dim(y)[1]
nsurvs <- dim(y)[2]
#--------------------- capture-recapture data -----------------------------#
# Capture-recapture data structure : rows = individuals; columns = capture occasions
mydata <- read.table("./capture-recapture-data.txt", header = FALSE) # load capture-recapture data
extra <- 250 # define large number of extra individual capture histories
n <- nrow(mydata) # number of observed individuals
M <- extra + n
xn <- rowSums(mydata)
x <- c(xn, rep(0, extra))
k <- ncol(mydata)
zerouse <- 0
#-------------------- specify model ------------------------------#
sink("CRandPO_HET.txt")
cat("
model{
    # Patch-occupancy (PO) likelihood
for (i in 1:nsites){
    truocc[i] ~ dbern(psi[i]) # true occupancy for site i
    logit(psi[i]) <- mupsi + betapsi[i] # occupancy probability as a function of μψ and a site random effect βs
    betapsi[i] ~ dnorm(0,sigmabeta)

    logit(ppo[i]) <- logitppo[i] # pOs = site-dependant detection probability from PO model
    logitppo[i] ~ dnorm(mean.p0,sigmaeps)

 # observation model
  for (t in 1:nsurvs) {
    y[i,t] ~ dbern(effp[i,t])
    effp[i,t] <- truocc[i]*ppo[i]
  }
}

    # Capture-recapture (CR) likelihood
 for(i in 1:M){
    x[i] ~ dbin(pee[i],k)
    pee[i] <- pcr[i]*step(N-i)  # pcr = individual detection probability pi from CR model
    logit(pcr[i]) <- mumup + xi*etap[i] # integration of an individual random effect ηi
    etap[i] ~ dnorm(0,sigmeta)
 }

    N ~ dpois(lamn)T(0,M)
    psi0 <- 1/(1+exp(-mupsi))  # μψ back transformed
    log(lamn) <- log(-log(1-psi0)) # definition of λ as a function of ψ
    mumup0 <- 1/(1+exp(-mumup)) # μp back transformed
    zerouse ~ dpois(comb)
    comb <- logfact(M)-logfact(N)+logfact((N-n)*step(N-n))

    # priors for occupancy probability
    mupsi ~ dnorm(0,0.2)
    sigmabeta <- 1/(sdbeta*sdbeta)
    sdbeta~ dunif(0,5)

     # priors for detection probability in PO likelihood
    mean.p0 ~ dnorm(0,0.1)
    sigmaeps <- 1/(sdeps*sdeps)
    sdeps ~ dunif(0,10)

     # priors for detection probability in CR likelihood
    mumup ~ dlogis(0,1)
    sigmeta ~ dgamma(1.5,37.5)
    xi ~ dnorm(0,1)
}
    ", fill = TRUE)
sink()

#-------------------- Jags stuff -----------------------#
# List of data
mydatax <- list(x = x, M = M, n = n, k = k, zerouse = zerouse, y = y, nsites = nsites, nsurvs = nsurvs)

# Parameters monitored
parameters <- c("N", "ppo", "mumup0", "psi0", "mupsi", "mean.p0", "sdeps", "taueps", "betapsi", "etap")


# Initial values
init1 <- list(mupsi = runif(1), truocc = as.numeric(apply(y, 1, sum) > 0), xi = rnorm(1, sd = 2), N = sample(n:M, 1), mumup = rnorm(1), taup = rlnorm(1))
init2 <- list(mupsi = runif(1), truocc = as.numeric(apply(y, 1, sum) > 0), xi = rnorm(1, sd = 2), N = sample(n:M, 1), mumup = rnorm(1), taup = rlnorm(1))
init3 <- list(mupsi = runif(1), truocc = as.numeric(apply(y, 1, sum) > 0), xi = rnorm(1, sd = 2), N = sample(n:M, 1), mumup = rnorm(1), taup = rlnorm(1))
inits <- list(init1, init2, init3)

#-------------------- Call jags from R -----------------------#
jmodel <- jags.model("CRandPO_HET.txt", mydatax, inits, n.chains = 3, n.adapt = 2500)
jsample <- coda.samples(jmodel, parameters, n.iter = 15000, thin = 1)
save(jsample, file = paste0("coda_samples_", task_id))

# Get summary statistics
s <- summary(jsample)

write.table(s$statistics, file = paste0("model_statistics_", task_id, ".txt"))
