## testing order of defining management times
require(spict)

## data
inpori <- pol$albacore

## Default assessment
inp <- check.inp(inpori)
rep1 <- fit.spict(inp)

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
all.equal(sumspict.parest(rep1),
          sumspict.parest(rep2), tolerance = 0.01)

all.equal(sumspict.parest(rep1),
          sumspict.parest(rep3), tolerance = 0.01)

all.equal(sumspict.parest(rep1),
          sumspict.parest(rep4), tolerance = 0.01)

## q.e.d.


## more tests needed?
b1=get.par("logBBmsy",rep1)[,2]
b2=get.par("logBBmsy",rep2)[,2]
b3=get.par("logBBmsy",rep3)[,2]


plot(b2,t='l',col=2,lwd=3)
lines(b1,col=4,lty=3,lwd=2)
lines(b3,col=3,lty=2,lwd=1)

## works! same parameter estimates -> how to make this work with retape()?

str(inp, list.len=length(inp))
str(inp2, list.len=length(inp2))
str(inp3, list.len=length(inp3))
