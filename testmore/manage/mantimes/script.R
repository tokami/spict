## Testing order of defining management times
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
test_this <- function(title="", expression) {
  out(title)
  tryCatch(out(capture.output(eval(expression)), sep = "\n"),
           error = function(e) out("Error:", e$message))
}

## data
inpori <- pol$albacore

## Default assessment
inp <- check.inp(inpori)
rep1 <- fit.spict(inp)



header("1: 4 sequences with the same result", append = FALSE)
###################################################################

## Postponing maninterval # 1st way
inpori$maninterval <- c(1991,1992)
inp2 <- check.inp(inpori)
rep2 <- fit.spict(inp2)

## Postponing maninterval # 2nd way
inp3 <- inp
inp3$maninterval <- c(1991,1992)
inp3$ini <- NULL
inp3$logmcovariatein <- NULL
inp3$ffacvec <- NULL
inp3$fconvec <- NULL
inp3 <- check.inp(inp3)
rep3 <- fit.spict(inp3)

## Postponing maninterval # New function
rep4 <- check.man.time(rep1, maninterval = c(1991,1992))

## All approaches should give same parameter estimates
out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep2), tolerance = 0.01))

out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep3), tolerance = 0.01))

out(all.equal(sumspict.parest(rep1),
          sumspict.parest(rep4), tolerance = 0.01))

b1=get.par("logBBmsy",rep1)[,2]
b2=get.par("logBBmsy",rep2)[,2]
b3=get.par("logBBmsy",rep3)[,2]
b4=get.par("logBBmsy",rep4)[,2]

## Last 3 reps are one year longer than first rep
out(all.equal(length(b1) + 1/rep1$inp$dteuler, length(b2), tolerance = 0.01))
out(all.equal(length(b1) + 1/rep1$inp$dteuler, length(b3), tolerance = 0.01))
out(all.equal(length(b1) + 1/rep1$inp$dteuler, length(b4), tolerance = 0.01))

## Last 3 reps should give same predictions
out(all.equal(sumspict.parest(rep2),
              sumspict.parest(rep3), tolerance = 0.01))

out(all.equal(sumspict.parest(rep3),
          sumspict.parest(rep4), tolerance = 0.01))
