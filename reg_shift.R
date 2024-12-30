rm(list = ls())
library(spict)
#source("functions.R")



### Set model name

model <- "quar - SSB - toDec2022 - 0 b"



### Set time limit (Calendar)

timeLimit <- 2023.00



### Read data
# Load and format catch

baseCatch <- read.table("quarCatch.csv", sep = ",", header = T)

obsC <- baseCatch[which(!is.na(baseCatch[, "Catch"])), "Catch"]/1e6
timeC <- baseCatch[which(!is.na(baseCatch[, "Catch"])), "calDate"]


# Load and format indices

baseIndice <- read.table("indiceSSB.csv", sep = ",", header = T)
baseIndice <- baseIndice[which(!is.na(baseIndice[, "SSB"])), ]

posiI1 <- which(baseIndice$Source == 1)
posiI2 <- which(baseIndice$Source == 2)
posiI3 <- which(baseIndice$Source == 3 & baseIndice$Season == 1)
posiI4 <- which(baseIndice$Source == 3 & baseIndice$Season == 2)

obsI1 <- baseIndice$SSB[posiI1]
obsI2 <- baseIndice$SSB[posiI2]
obsI3 <- baseIndice$SSB[posiI3]
obsI4 <- baseIndice$SSB[posiI4]
  
timeI1 <- baseIndice$Date[posiI1]
timeI2 <- baseIndice$Date[posiI2]
timeI3 <- baseIndice$Date[posiI3]
timeI4 <- baseIndice$Date[posiI4]



### Format data base to SPiCT

Anch <- list(obsC = obsC, timeC = timeC)

Anch$obsI <- list()
Anch$obsI[[1]] <- obsI1
Anch$obsI[[2]] <- obsI2
Anch$obsI[[3]] <- obsI3
Anch$obsI[[4]] <- obsI4

Anch$timeI <- list()
Anch$timeI[[1]] <- timeI1
Anch$timeI[[2]] <- timeI2
Anch$timeI[[3]] <- timeI3
Anch$timeI[[4]] <- timeI4

#Anch <- check.inp(Anch)



### Set priors

Anch$priors$logK <- c(log(17.1), 0.20, 1) # K (Csirke et al. 1996)

Anch$priors$logq <- list(c(log(1), 0.20, 1), c(log(1), 0.15, 1), c(log(1), 0.10, 1), c(log(1), 0.15, 1)) # q

Anch$priors$logn <- c(log(0.599), 0.342, 1) # from Thorson 2012

Anch$phases$logsdc <- -1

Anch$ini$logsdc <- log(0.05) # sd catch

Anch$priors$logsdi <- list(c(log(0.10), 0.15, 1), c(log(0.10), 0.15, 1), c(log(0.10), 0.15, 1), c(log(0.10), 0.15, 1))

Anch$priors$logbkfrac <- c(log(1), 0.05, 1) # B/K B at beg pretty close to K

Anch$msytype <- "d"

Anch$dteuler <- 1/16
dteuler <- 1/16
## regime shift
nt1 = length(seq(1950,1971,dteuler))
nt2 = length(seq(1971+dteuler,1993,dteuler))
nt3  = length(seq(1993+dteuler,2024,dteuler))
Anch$MSYregime <- as.factor(c(rep(1,nt1), rep(2,nt2), rep(3,nt3) ))
Anch$ini$logqf = log(0.05)
Anch <- check.inp(Anch)



### Fit the model

res  <- fit.spict(Anch)
res



### Estimte osa residuals

osa <- calc.osa.resid(res)



### Estimate cor among pars

corPars <- cov2cor(res$cov.fixed)



### Check ini

# set.seed(123)
# check.ini_res <- check.ini(res, ntrials=30) # 05Ago23 24from30 80%



### Plots default

pdf(file <- paste("preplot - ", model, ".pdf", sep = ""))
plotspict.data(Anch)
dev.off()

pdf(file <- paste("outputs - ", model, ".pdf", sep = ""))
plot(res)
dev.off()

pdf(file <- paste("residuals - ", model, ".pdf", sep = ""))
plotspict.diagnostic(osa)
dev.off()



### Save outputs

save(res, file = paste("outputs - ", model, ".RData", sep = ""))
save(osa, file = paste("residuals - ", model, ".RData", sep = ""))



### Get data

# Get Catch (obs and est)

time_obsC <- timeC

time_estC <- res$inp$timeCpred
time_estC_2Rem <- which(time_estC >= timeLimit)
time_estC <- time_estC[-time_estC_2Rem]

obsC = obsC
estC = as.numeric(get.par(parname = "logCpred", rep = res, exp = T)[-time_estC_2Rem,2])

# Get Indices (obs and est)

time_obsI <- sort(c(timeI1, timeI2, timeI3, timeI4))
obsI = c(obsI1, obsI2, obsI3, obsI4)
obsI = obsI[order(c(timeI1, timeI2, timeI3, timeI4))]

time_estI <- time_obsI
estI = as.numeric(get.par(parname = "logIpred", rep = res, exp = T)[,2])
estI = estI[order(c(timeI1, timeI2, timeI3, timeI4))]


# Get continuous series (B and F)

timeSeries <- as.numeric(row.names(get.par(parname = "logB", rep = res, exp = T)))
timeSeries_2Rem <- which(timeSeries >= timeLimit)
timeSeries <- timeSeries[-timeSeries_2Rem]

estB = as.numeric(get.par(parname = "logB", rep = res, exp = T)[-timeSeries_2Rem, 2])
estB2 = as.numeric(get.par(parname = "logB", rep = res, exp = F)[-timeSeries_2Rem, 2]) # for annual means

estF = as.numeric(get.par(parname = "logFs", rep = res, exp = T)[-timeSeries_2Rem, 2])
estF2 = as.numeric(get.par(parname = "logFs", rep = res, exp = F)[-timeSeries_2Rem, 2]) # for annual means


# Get BRP according to type

n_out <- as.numeric(get.par(parname = "logn", rep = res, exp = T))[2]

MSY_type <- if(n_out >= 1) {MSY_type = "logMSYs"} else {MSY_type = "logMSYd"}
MSY <- as.numeric(get.par(parname = MSY_type, rep = res, exp = T))[2]
MSY_l <- as.numeric(get.par(parname = MSY_type, rep = res, exp = T))[1]
MSY_u <- as.numeric(get.par(parname = MSY_type, rep = res, exp = T))[3] 

FMSY_type = if(n_out >= 1) {FMSY_type = "logFmsys"} else {FMSY_type = "logFmsyd"}
FMSY = as.numeric(get.par(parname = FMSY_type, rep = res, exp = T))[2]
FMSY_l = as.numeric(get.par(parname = FMSY_type, rep = res, exp = T))[1]
FMSY_u = as.numeric(get.par(parname = FMSY_type, rep = res, exp = T))[3] 

BMSY_type = if(n_out >= 1) {BMSY_type = "logBmsys"} else {BMSY_type = "logBmsyd"}
BMSY = as.numeric(get.par(parname = BMSY_type, rep = res, exp = T))[2]
BMSY_l = as.numeric(get.par(parname = BMSY_type, rep = res, exp = T))[1]
BMSY_u = as.numeric(get.par(parname = BMSY_type, rep = res, exp = T))[3] 


# Get relative (to BRP) values 

BBMSY = as.numeric(get.par(parname = "logBBmsy", rep = res, exp = T)[-timeSeries_2Rem, 2])
BBMSY2 = as.numeric(get.par(parname = "logBBmsy", rep = res, exp = F)[-timeSeries_2Rem, 2])

FFMSY = as.numeric(get.par(parname = "logFFmsy", rep = res, exp = T)[-timeSeries_2Rem, 2])
FFMSY2 = as.numeric(get.par(parname = "logFFmsy", rep = res, exp = F)[-timeSeries_2Rem, 2])

# base_estC <- data.frame(calDate = time_estC, estCatch = round(estC*1e6,0))
# base2Exp_C <- merge(baseCatch, base_estC, by = "calDate", all.x = T)
# base2Exp_C <- base2Exp_C[,-2]
# write.table(base2Exp_C, "base2Exp_C.csv", sep = ",", col.names = T, row.names = F)

# base_estI <- data.frame(calDate = time_estI, estI = estI)
# base2Exp_I <- merge(baseIndice, base_estI, by = "calDate", all.x = T)
# base2Exp_I <- base2Exp_I[,-c(2,3,4)]
# write.table(base2Exp_I, "base2Exp_I.csv", sep = ",", col.names = T, row.names = F)

# png(filename = "compareAbsVal.png", width = 750, height = 500, res = 100)
# yearF1 <- colMeans(matrix(estF, ncol = length(estF)/16))
# yearF2 <- exp(colMeans(matrix(estF2, ncol = length(estF2)/16)))
# plot(yearF1, type = "l", xlab = "", ylab = "annual F", xaxt = "n")
# lines(yearF2, type = "l", col = 2)
# legend("topright", legend = c("mean(exp(logF)", "exp(mean(logF))"), lty = c(1, 1), col = c(1, 2), bty = "n")
# axis(side = 1, at = 1:length(yearF1), labels = seq(1953, 2022, 1))
# dev.off()

### Get breakpoints

require(bcp)

yearB = exp(colMeans(matrix(estB2, ncol = length(estB2)/16)))

yearF = exp(colMeans(matrix(estF2, ncol = length(estF2)/16)))

year_relB = exp(colMeans(matrix(BBMSY2, ncol = length(BBMSY2)/16)))
year_relF = exp(colMeans(matrix(FFMSY2, ncol = length(FFMSY2)/16)))

Years =  as.numeric(unique(trunc(timeSeries)))

set.seed(123)
bcp1 = bcp(yearB, return.mcmc = T)

postD = bcp1$posterior.prob
bkp = which(postD %in% sort(x = postD, decreasing = T)[1:2])

bkp1 <- Years[bkp[1]]
bkp2 <- Years[bkp[2]]



# Save outputs

write.table(data.frame(Date = timeSeries, SSB = estB, relSSB = BBMSY), "SSBout.csv", sep = ",", row.names = F)
write.table(data.frame(Date = Years, SSB = yearB, relSSB = year_relB), "annualSSBout.csv", sep = ",", col.names = T, row.names = F)

write.table(data.frame(Date = timeSeries, Fm = estF, relF = FFMSY), "Fout.csv", sep = ",", row.names = F)
write.table(data.frame(Date = Years, Fm = yearF, relF = year_relF), "annualFout.csv", sep = ",", col.names = T, row.names = F)



# Figures
Figure1 = figure_dataUsed(model = model)
Figure2 = figure_goF(model = model)
Figure3 = figure_contSeriesBRP_P(model = model)
Figure4 = figure_regimes(model = model)
Figure5 = figure_priors(model = model)

# Figure3 = figure_contSeries(model = model)
# Figure4 = figure_sd_contSeries(model = model)
# Figure5 = figure_ma_contSeries(model = model)
# Figure6 = figure_sd_maSeries(model = model)
# Figure7 = figure_ciBRP(model = model)
# Figure8 = figure_dataUsed_P(model = model)
# Figure9 = figure_goF_P(model = model)
# Figure10 = figure_contSeries_P(model = model)
# Figure11 = figure_ciBRP_P(model = model)
