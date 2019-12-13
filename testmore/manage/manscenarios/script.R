## Testing new management functionality
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 11/12/2019

## set seed
set.seed(123)

## load spict
require(spict)

## Test functions
out <- function(..., sep = "", append = TRUE)
  cat(..., "\n", file = "res.out", append = append, sep = sep)
get_nchar <- function(...)
  nchar(paste(as.character(unlist(list(...))), collapse = " "))
header <- function(..., append = TRUE)
  out("\n", ..., "\n", rep("=", get_nchar(...)), "\n", append = append)
test_this <- function(title = "", expression) {
  out(title)
  tryCatch(out(capture.output(eval(expression)), sep = "\n"),
           error = function(e) out("Error:", e$message))
}


header("1: General tests", append = FALSE)
######################################################

## load data
inp <- pol$albacore
inp <- check.inp(inp)

## fit spict
rep <- fit.spict(inp)

## manage
test_this("1.1: run manage all scenarios",{
    repman <- manage(rep)
})

test_this("1.2: print manage summary",{
    sumspict.manage(repman)
})

test_this("1.3: use deprecated summary function",{
    mansummary(repman)
})

test_this("1.4: Make selection of management scenarios",{
    sumspict.manage(man.select(repman,c(1,5,6,8)))
})


test_this("1.5: Change base management scenario (manbase and of rep)" ,{
    names(repman$man)
    repman2 <- manbase.select(repman,"ices")
    all.equal(sumspict.predictions(repman2),
              sumspict.predictions(repman$man[[8]]), tolerance = 0.01)
    all.equal(sumspict.predictions(repman2$manbase),
              sumspict.predictions(repman$man[[8]]), tolerance = 0.01)

})


## testing adding scenarios
test_this("1.6: Adding scenarios",{
    repman1 <- man.scenario(repman)
    repman1 <- man.scenario(repman1, ffac = 0.5)
    repman1 <- man.scenario(repman1, ffac = 0.25)
})

## scenarios should not be overwritten:
out(length(repman1$man))

## names should bemeaningful:
out(names(repman1$man))

## sumspict should still work
out(sumspict.manage(repman1))



## change management interval:
test_this("1.7: Change management interval" ,{
    inp <- pol$albacore
    inp$maninterval <- c(1992,1993)
    res <- fit.spict(inp)
    repman <- manage(res, c(1,5))
    ## should be the same as this:
    repman2 <- manage(res, c(1,5), maninterval = c(1992,1993))

    all.equal(sumspict.predictions(repman),
              sumspict.predictions(repman2), tolerance = 0.01)
})




## report tac only
test_this("1.8: Calculate the TAC for all man scenarios" ,{
    man.tac(repman1)
})

## report tac only
test_this("1.8: Calculate TAC for fitted spict object without $man" ,{
    man.scenario(res, getFit=FALSE)
})


test_this("1.9: Add scenario to fitted spict object" ,{
    res <- man.scenario(res, cfac=0.5, scenarioTitle = "Half catch")
})

inp <- check.inp(pol$albacore)
inp$maninterval <- inp$maninterval + 1
inp2 <- check.man.time(inp)
test_this("1.10: Check man times in input list" ,{
    length(inp$seasonindex2) + 1/inp$dteuler == length(inp2$seasonindex2)
})


res$inp$maninterval <- res$inp$maninterval + 1
res2 <- check.man.time(res)
test_this("1.11: Check man times in rep list" ,{
    length(res$inp$logmcovariatein) + 1/res$inp$dteuler == length(res2$inp$logmcovariatein)
})


header("2: Testing man scenarios with annual data", append = FALSE)
######################################################

out("2.1: no intermediate year")
#################################
inp <- pol$lobster
inp$maninterval <- c(1991, 1992)
inp$maneval <- 1992
fit <- fit.spict(inp)
print.man.timeline(fit)

## Fishing at Fmsy
test_this("2.1.1:", {round(man.scenario(fit, getFit = FALSE),3)})
test_this("2.1.2:", {round(man.scenario(fit, fractiles = list(catch=0.2), getFit = FALSE),3)})
test_this("2.1.3:", {round(man.scenario(fit, fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})

## Fishing at Fmsy with PA buffer
test_this("2.1.4:", {round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.8), getFit = FALSE),3)})
test_this("2.1.5:", {round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.9), getFit = FALSE),3)})
test_this("2.1.6:", {round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.93), getFit = FALSE),3)})
test_this("2.1.7:", {round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.95), getFit = FALSE),3)})
test_this("2.1.8:", {round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.6), getFit = FALSE),3)})
test_this("2.1.9:", {round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.65), getFit = FALSE),3)})
test_this("2.1.10:", {round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.7), getFit = FALSE),3)})
test_this("2.1.11:", {round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.9),
                                         fractiles = list(catch=0.2), getFit = FALSE),3)})
test_this("2.1.12:", {round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.9),
                                         fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})

## MSY hockey-stick rule
test_this("2.1.13:", {round(man.scenario(fit, breakpointB = 0.3, getFit = FALSE),3)})
test_this("2.1.14:", {round(man.scenario(fit, breakpointB = 0.5, getFit = FALSE),3)})
test_this("2.1.15:", {round(man.scenario(fit, breakpointB = 0.5, fractiles = list(catch=0.2), getFit = FALSE),3)})
test_this("2.1.16:", {round(man.scenario(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1),
                                         getFit = FALSE),3)})
test_this("2.1.17:", {round(man.scenario(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1, bbmsy=0.1),
                                         getFit = FALSE),3)})

## MSY hockey-stick rule with PA buffer
test_this("2.1.18:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8),
                                         getFit = FALSE),3)})
test_this("2.1.19:", {round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9),
                                         getFit = FALSE),3)})
test_this("2.1.20:", {round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.93),
                                         getFit = FALSE),3)})
test_this("2.1.21:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.95),
                                         getFit = FALSE),3)})
test_this("2.1.22:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.6),
                                         getFit = FALSE),3)})
test_this("2.1.23:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.65),
                                         getFit = FALSE),3)})
test_this("2.1.24:", {round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.5, prob=0.7),
                                         getFit = FALSE),3)})
test_this("2.1.25:", {round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9),
                                         fractiles = list(catch=0.2), getFit = FALSE),3)})
test_this("2.1.26:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                         fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})


out("2.2: intermediate year with constant F")
#################################
inp$maninterval <- c(1992, 1993)
inp$maneval <- 1993
fit <- fit.spict(inp)

## Fishing at Fmsy
test_this("2.2.1:", {round(man.scenario(fit, getFit = FALSE),3)})
test_this("2.2.2:", {round(man.scenario(fit, fractiles = list(catch=0.2), getFit = FALSE),3)})

## Fishing at Fmsy with biomass safeguard
test_this("2.2.3:", {round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.9), getFit = FALSE),3)})
test_this("2.2.4:", {round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.6), getFit = FALSE),3)})

## MSY hockey-stick rule
test_this("2.2.5:", {round(man.scenario(fit, breakpointB = 0.3, getFit = FALSE),3)})
test_this("2.2.6:", {round(man.scenario(fit, breakpointB = 0.5,
                                        fractiles = list(catch=0.2, ffmsy=0.1, limitB=0.1), getFit = FALSE),3)})

## MSY hockey-stick rule with biomass safeguard
test_this("2.2.7:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8),
                                        getFit = FALSE),3)})
test_this("2.2.8:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                        fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})


out("2.3: intermediate year with constant catch")
#################################
inp$maninterval <- c(1992, 1993)
inp$maneval <- 1993
fit <- fit.spict(inp)
lastC <- tail(inp$obsC,1)

## Fishing at Fmsy
test_this("2.3.1:", {round(man.scenario(fit, catchIntermediateYear = lastC, getFit = FALSE),3)})
test_this("2.3.2:", {round(man.scenario(fit, fractiles = list(catch=0.2), catchintermediateYear = lastC,
                                        getFit = FALSE),3)})

## Fishing at Fmsy with biomass safeguard
test_this("2.3.3:", {round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.9), catchIntermediateYear = lastC,
                                        getFit = FALSE),3)})
test_this("2.3.4:", {round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.6), catchIntermediateYear = lastC,
                                        getFit = FALSE),3)})

## MSY hockey-stick rule
test_this("2.3.5:", {round(man.scenario(fit, breakpointB = 0.3, catchIntermediateYear = lastC,
                                        getFit = FALSE),3)})
test_this("2.3.6:", {round(man.scenario(fit, breakpointB = 0.5,
                                        fractiles = list(catch=0.2, ffmsy=0.1, limitB=0.1),
                                        catchIntermediateYear = lastC, getFit = FALSE),3)})

## MSY hockey-stick rule with biomass safeguard
test_this("2.3.7:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8),
                                        catchIntermediateYear = lastC, getFit = FALSE),3)})
test_this("2.3.8:", {round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9),
                                        fractiles = list(catch=0.2, ffmsy = 0.1),
                                        catchIntermediateYear = lastC, getFit = FALSE),3)})


header("3: Testing management scenarios with seasonal data")
#####################################################################

nt <- 50
inp <- list(nseasons = 4, splineorder = 3)
inp$timeC <- seq(0, nt - 1 / inp$nseasons, by = 1 / inp$nseasons)
inp$timeI <- seq(0.1, nt - 1 / inp$nseasons, by = 0.5)
inp$ini <- list(logK = log(1000), logm=log(800), logq = log(1), logn = log(2),
                logbkfrac = log(0.9), logsdf = log(0.3), logF0 = log(0.8),
                logphi = log(c(0.3, 0.5, 1.8)))
inpsim <- sim.spict(inp)
fit2 <- fit.spict(inpsim)

out(fit2$opt$convergence)


out("3.1: standard advice")
###############################
inpsim$maneval <- 31
inpsim$maninterval <- c(30,31)
inp <- check.inp(inpsim)
fit <- fit.spict(inpsim)

## Fishing at Fmsy
out(round(man.scenario(fit, getFit = FALSE),3)})
out(round(man.scenario(fit, fractiles = list(catch=0.2), getFit = FALSE),3)})
out(round(man.scenario(fit, fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})

## Fishing at Fmsy with PA buffer
out(round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.8), getFit = FALSE),3)})
out(round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.65), getFit = FALSE),3)})
out(round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.9), fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})

## MSY hockey-stick rule
out(round(man.scenario(fit, breakpointB = 0.3, getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1, bbmsy=0.1), getFit = FALSE),3)})

## MSY hockey-stick rule with PA buffer
out(round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.3, prob=0.8), getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9), getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9), fractiles = list(catch=0.2), getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9), fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})


out("3.1: in year advice")
###############################
inpsim$maneval <- 31.5
inpsim$maninterval <- c(30.5,31.5)
inp <- check.inp(inpsim)
fit <- fit.spict(inpsim)

## Fishing at Fmsy
out(round(man.scenario(fit, getFit = FALSE),3)})
out(round(man.scenario(fit, fractiles = list(catch=0.2), getFit = FALSE),3)})
out(round(man.scenario(fit, fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})

## Fishing at Fmsy with PA buffer
out(round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.8), getFit = FALSE),3)})
out(round(man.scenario(fit, safeguardB = list(limitB=0.3, prob=0.9), getFit = FALSE),3)})
out(round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.6), getFit = FALSE),3)})
out(round(man.scenario(fit, safeguardB = list(limitB=0.5, prob=0.65), getFit = FALSE),3)})

## MSY hockey-stick rule
out(round(man.scenario(fit, breakpointB = 0.3, getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.5, getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.5, fractiles = list(catch=0.2, ffmsy=0.1), getFit = FALSE),3)})


## MSY hockey-stick rule with PA buffer
out(round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9), getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.5, safeguardB = list(limitB=0.3, prob=0.9), fractiles = list(catch=0.2), getFit = FALSE),3)})
out(round(man.scenario(fit, breakpointB = 0.3, safeguardB = list(limitB=0.5, prob=0.9), fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})

## Fishing at Fmsy with TAC during assessment year
out(round(man.scenario(fit, catch_pred = 4, getFit = FALSE),3)})
out(round(man.scenario(fit, catch_pred = 4, fractiles = list(catch=0.2), getFit = FALSE),3)})
out(round(man.scenario(fit, catch_pred = 4, fractiles = list(catch=0.2, ffmsy = 0.1), getFit = FALSE),3)})




header("4: Testing management scenarios with mixed data")
#####################################################################

nt <- 50
inp <- list(nseasons=4, splineorder=3)
inp$timeC <- c(seq(0, nt/2-1, by=1),seq(nt/2, nt-1/inp$nseasons, by=1/inp$nseasons))
inp$timeI <- seq(0.1, nt-1/inp$nseasons, by=0.5)
inp$ini <- list(logK=log(1000), logm=log(800), logq=log(1), logn=log(2),
                logbkfrac=log(0.9), logsdf=log(0.3), logF0=log(0.8),
                logphi=log(c(0.3, 0.5, 1.8)))
inpsim <- sim.spict(inp)


inp$timeI <- seq(0.6, 29.6, by = 1)
inp$ini <- list(logK = log(100), logm = log(60), logq = log(1),
                logbkfrac = log(1), logsdf = log(0.3), logF0 = log(0.5),
                logphi = log(c(0.05, 0.1, 1.8)))
inpsim <- sim.spict(inp)
fit3 <- fit.spict(inpsim)

out(fit3$opt$convergence)





header("5: Challenging new functions")
#####################################################################

test_this("4.1", {
    sumspict.manage(fit)
})
