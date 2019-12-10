## Testing some new functions and modifications
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 10/12/2019

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


header("1: annual data", append = FALSE)
## load data
inp <- pol$albacore
inp <- check.inp(inp)

## fit spict
rep <- fit.spict(inp)

## manage
test_this("1.1: run manage",{
    repman <- manage(rep, scenarios = c(1,5))
})

out(sumspict.manage(repman))



## testing adding scenarios
test_this("1.1: run manage",{
    repman1 <- man.scenario(repman)
    repman1 <- man.scenario(repman1, ffac = 0.5)
    repman1 <- man.scenario(repman1, ffac = 0.25)
})

## scenarios should not be overwritten:
length(repman1$man)

## names should bemeaningful:
names(repman1$man)


## BUG: ypred=2 in sumspict.manage does not work


## change management interval:

inp <- pol$albacore
inp$maninterval <- c(1992,1993)
res <- fit.spict(inp)
res <- manage(res)
## should be the same as this:
res <- manage(res, maninterval = c(1992,1993))


## report tac only
tac <- man.scenario(rep, ffac = 1, getFit = FALSE)


## use fixed catch in intermediate year



header("2: seasonal data")

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

test_this("2.1: calculate Bmsy/K ratio", {
  round(calc.bmsyk(fit2), 3)
})


test_this("2.2: calculate order of magnitude", {
  round(calc.om(fit2), 3)
})



header("3: Tests with mixed data")

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

test_this("3.1: calculate Bmsy/K ratio", {
  round(calc.bmsyk(fit3), 3)
})

test_this("3.2: calculate order of magnitude", {
  round(calc.om(fit3), 3)
})



header("4: incorrect input")

out("4.1: wrong input")

test_this("4.1.1: Bmsy/K ratio", {
  calc.bmsyk(inp)
})
