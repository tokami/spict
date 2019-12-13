## Testing new function to retape fitted spict objects
## T.K. Mildenberger <t.k.mildenberger@gmail.com>
## 12/12/2019

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


## load data
inp <- pol$albacore
inp <- check.inp(inp)

## fit spict
rep <- fit.spict(inp)




header("1: Changing prediction horizon", append = FALSE)
#####################################################

test_this("1.1: providing input list without change",{
    rept <- spict:::retape.spict(rep, inp)
})

out(all.equal(sumspict.predictions(rep),
              sumspict.predictions(rept),tolerance = 0.01))


inpt <- check.man.time(inp, maninterval = c(1991,1992), verbose = FALSE)
test_this("1.2: extending prediction horizon with maninterval",{
    rept <- spict:::retape.spict(rep, inpt)
})

out(all(rep$inp$maninterval + 1 == rept$inp$maninterval))

out(round(sumspict.predictions(rept),3))


inpt <- check.man.time(inp, maneval = 1996, verbose = FALSE)
test_this("1.3: extending prediction horizon with maneval",{
    rept <- spict:::retape.spict(rep, inpt)
})

out(all(rep$inp$maneval + 6 == rept$inp$maneval))

out(round(sumspict.predictions(rept),3))


inpt$maneval <- 1989
inpt <- check.man.time(inpt, verbose = FALSE)
test_this("1.3: Shortening prediction horizon",{
    rept <- spict:::retape.spict(rep, inpt)
})

out(all(rep$inp$maneval - 1 == rept$inp$maneval))

out(round(sumspict.predictions(rept),3))


inp$maninterval <- c(1990,1990.5)
inpt <- check.man.time(inp, verbose = FALSE)
test_this("1.4: Half year management interval",{
    rept <- spict:::retape.spict(rep, inpt)
})

out(get.par("logCp",rept,exp=TRUE)[2] < get.par("logCp",rep,exp=TRUE)[2])

out(round(sumspict.predictions(rept),3))




header("2: Changing F vector", append = FALSE)
#####################################################

## Testing ffac
inpt <- make.ffacvec(inp, ffac = 0.4)

## old code
plt <- inpt$parlist
datint <- make.datin(inpt)
objt <- make.obj(datint, plt, inpt, phase=1)
objt$retape()
objt$fn(rep$opt$par)
rep1 <- TMB::sdreport(objt)
rep1$inp <- inpt
rep1$obj <- objt
rep1$opt <- rep$opt
class(rep1) <- "spictcls"

## new function
rep2 <- spict:::retape.spict(rep, inpt)

test_this("2.1: ffac == 0.4",{
all.equal(sumspict.predictions(rep1),
          sumspict.predictions(rep2), tolerance = 0.01)
})

## Testing fcon
inpt <- make.fconvec(inp, fcon = 1.3)

## old code
plt <- inpt$parlist
datint <- make.datin(inpt)
objt <- make.obj(datint, plt, inpt, phase=1)
objt$retape()
objt$fn(rep$opt$par)
rep1 <- TMB::sdreport(objt)
rep1$inp <- inpt
rep1$obj <- objt
rep1$opt <- rep$opt
class(rep1) <- "spictcls"

## new function
rep2 <- spict:::retape.spict(rep, inpt)

test_this("2.2: fcon == 1.3",{
all.equal(sumspict.predictions(rep1),
          sumspict.predictions(rep2), tolerance = 0.01)
})





header("3: Forcing errors", append = FALSE)
#####################################################

test_this("3.1: wrong input to make.ffacvec",{
    inpt <- make.ffacvec(pol$albacore, ffac=1.2)
})

test_this("3.2: wrong input to make.fconvec",{
    inpt <- make.ffacvec(pol$albacore, ffac=1.2)
})

test_this("3.2: incorrect rep",{
    rept <- spict:::retape.spict(inpt, inpt)
})

test_this("3.3: non meaningful inp",{
    rept <- spict:::retape.spict(rep, list())
})

inpt <- inp
inpt$maninterval <- c(1991,1992)
test_this("3.4: manually changed maninterval",{
    rept <- spict:::retape.spict(rep, inpt)
})

## Non-converged model
inp <- pol$lobster
inp$obsC <- rnorm(46, 2000, 3)
fit4 <- fit.spict(inp)

out(fit4$opt$convergence)

inpt <- make.ffacvec(check.inp(inp), ffac=1.2)

test_this("3.5: non-converged model",{
    rept <- spict:::retape.spict(fit4, inpt)
})

out(rept$opt$convergence)
