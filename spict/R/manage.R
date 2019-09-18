# Stochastic surplus Production model in Continuous-Time (SPiCT)
#    Copyright (C) 2015-2016  Martin W. Pedersen, mawp@dtu.dk, wpsgodd@gmail.com
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



#' @name manage
#' @title Calculate predictions under different management scenarios
#' @details Scenarios that are currently implemented include:
#' \itemize{
#'   \item{"1"}{ Keep the catch of the current year (i.e. the last observed catch).}
#'   \item{"2"}{ Keep the F of the current year.}
#'   \item{"3"}{ Fish at Fmsy i.e. F=Fmsy.}
#'   \item{"4"}{ No fishing, reduce to 1\% of current F.}
#'   \item{"5"}{ Reduce F by X\%. Default X = 25.}
#'   \item{"6"}{ Increase F by X\%. Default X = 25.}
#' }
#' @param repin Result list from fit.spict().
#' @param scenarios Vector of integers specifying which scenarios to run. Default: 'all'.
#' @param manstart Year that management should be initiated.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' repman <- manage(rep)
#' mansummary(repman) # To print projections
manage <- function(repin, scenarios='all', manstart=NULL, dbg=0, catch=NULL, catchList=NULL){
    if (scenarios == 'all'){
        scenarios <- 1:6
    }
    if (is.null(manstart)){
        manstart <- repin$inp$manstart
    } else {
        repin$inp$manstart <- manstart
    }
    maninds <- which(repin$inp$time >= manstart)
    inpin <- repin$inp

    timelastobs <- repin$inp$time[repin$inp$indlastobs]
    if (!repin$inp$timepredc < timelastobs){
        # Add prediction horizons
        inpin$timepredc <- repin$inp$timepredc
        inpin$timepredi <- repin$inp$timepredi
        inpin$manstart <- repin$inp$manstart
        repman <- list() # Output list
        attr(repman, "scenarios") <- scenarios 
        if (1 %in% scenarios){
            # 1. Specify the catch, which will be taken each year in the prediction period
            lastyearidxs <- min( which( cumsum(rev(inpin$dtc))>=1 ) ) ## warning: this will not make sense with subannual/mixed data with missing values
            if(is.null(catch)) catch <- sum(tail(inpin$obsC, lastyearidxs))
            repman[[1]] <- take.c(catch, inpin, repin, dbg=dbg, catchList=catchList)
        }
        if (2 %in% scenarios){
            # Keep current F
            fac2 <- 1.0
            repman[[2]] <- prop.F(fac2, inpin, repin, maninds, dbg=dbg)
        }
        if (3 %in% scenarios){
            # Fish at Fmsy
            Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
            Flast <- get.par('logF', repin, exp=TRUE)[repin$inp$indpred[1], 2]
            fac3 <- Fmsy / Flast
            repman[[3]] <- prop.F(fac3, inpin, repin, maninds, dbg=dbg)
        }
        if (4 %in% scenarios){
            # No fishing, reduce to 0.1% of last F
            fac4 <- 0.001
            repman[[4]] <- prop.F(fac4, inpin, repin, maninds, dbg=dbg)
        }
        if (5 %in% scenarios){
            # Reduce F by X%
            fac5 <- 0.75
            repman[[5]] <- prop.F(fac5, inpin, repin, maninds, dbg=dbg)
        }
        if (6 %in% scenarios){
            # Increase F by X%
            fac6 <- 1.25
            repman[[6]] <- prop.F(fac6, inpin, repin, maninds, dbg=dbg)
        }
        repin$man <- repman
        # Create an baseline F trajectory with constant F and store
        repin$manbase <- prop.F(fac=1, inpin, repin, maninds, dbg=dbg)
    } else {
        stop('Error: Could not do management calculations because prediction horizon is too short. Increase inp$timepredc to be at least one timestep into the future.\n')
    }
    return(repin)
}


#' @name prop.F
#' @title Calculate management for changing F by a given factor.
#' @param fac Factor to multiply current F with.
#' @param inpin Input list.
#' @param repin Results list.
#' @param maninds Indices of time vector for which to apply management.
#' @param corF Make correction to F process such that the drift (-0.5*sdf^2*dt) is cancelled and F remains constant in projection mode
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return List containing results of management calculations.
#' @export
prop.F <- function(fac, inpin, repin, maninds, corF=FALSE, dbg=0){
    inpt <- check.inp(inpin)
    inpt <- make.ffacvec(inpt, fac)
    # Make object
    plt <- repin$obj$env$parList(repin$opt$par)
    datint <- make.datin(inpt, dbg=dbg)
    objt <- make.obj(datint, plt, inpt, phase=1)
    objt$fn(repin$opt$par)
    ## repmant <- sdreport(objt)
    verflag <- as.numeric(gsub('[.]', '', as.character(packageVersion('TMB')))) >= 171
    if (verflag) { 
      repmant <- sdreport(objt,
                          getJointPrecision=repin$inp$getJointPrecision,
                          bias.correct=repin$inp$bias.correct,
                          bias.correct.control=repin$inp$bias.correct.control,
                          getReportCovariance=repin$inp$getReportCovariance)
    } else {
      repmant <-sdreport(objt,
                         getJointPrecision=repin$inp$getJointPrecision,
                         bias.correct=repin$inp$bias.correct,
                         bias.correct.control=repin$inp$bias.correct.control)
    }
    repmant$inp <- inpt
    repmant$obj <- objt
    repmant$opt <- list(convergence=0)
    if (!is.null(repmant)){
        class(repmant) <- "spictcls"
    }
    return(repmant)
}


#' @name take.c
#' @title Calculate management when taking a constant catch (proxy for setting a TAC).
#' @param catch Take this catch 'dtpredc' ahead from manstart time 
#' @param inpin Input list.
#' @param repin Results list.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more output.
#' @param sdfac Take catch with this 'stdevfacC' (default = 1e-3) 
#' @return List containing results of management calculations.
#' @export
take.c <- function(catch, inpin, repin, dbg=0, sdfac=1e-3, catchList=NULL){
    
    inpt <- inpin
    if(is.null(catchList)){
        tmpTime <- repin$inp$timeCpred  
        maninds <- which(tmpTime >= inpin$manstart)
        inpt$timeC <- c( inpt$timeC, tmpTime[maninds] )
        inpt$obsC <- c( inpt$obsC, rep(catch, length(maninds)) )
        inpt$stdevfacC <- c(inpt$stdevfacC, rep(sdfac, length(maninds)) )  
        inpt$dtc <- c(inpt$dtc, rep(inpt$dtpredc, length(maninds)) )
    } else {
        inpt$timeC <- c( inpt$timeC, catchList$timeC )
        inpt$obsC <- c( inpt$obsC, catchList$obsC )
        if(is.null(catchList$stdevfacC))
            inpt$stdevfacC <- c(inpt$stdevfacC, rep(sdfac, length(catchList$timeC)) )  else
            inpt$stdevfacC <- c(inpt$stdevfacC, catchList$stdevfacC)
        
        inpt$dtc <- c(inpt$dtc, catchList$dtc )
    }

    
    inpt <- check.inp(inpt)
    # Make TMB data and object
    plt <- repin$obj$env$parList(repin$opt$par)
    datint <- make.datin(inpt, dbg=dbg)
    objt <- make.obj(datint, plt, inpt, phase=1)
    # Get updated sd report
    objt$fn(repin$opt$par)
    repmant <- sdreport(objt)
    repmant$inp <- inpt
    repmant$obj <- objt
    repmant$opt <- list(convergence=0)
    if (!is.null(repmant)){
        class(repmant) <- "spictcls"
    }
    return(repmant)
}


#' @name mansummary
#' @title Print management summary.
#' @param repin Result list as output from manage().
#' @param ypred Show results for ypred years from manstart.
#' @param include.EBinf Include EBinf/Bmsy in the output.
#' @param include.unc Include uncertainty of management quantities.
#' @param verbose Print more details on observed and predicted time intervals.
#' @return Data frame containing management summary.
#' @export
mansummary <- function(repin, ypred=1, include.EBinf=FALSE, include.unc=TRUE, verbose=TRUE){
    if (!'man' %in% names(repin)){
        stop('Management calculations not found, run manage() to include them.')
    } else {
        repman <- repin$man
        rep <- repin$manbase
        # Calculate percent difference.
        get.pdelta <- function(rep, repman, indstart, indnext, parname='logB'){
            val <- get.par(parname, rep, exp=TRUE)[indstart, 2]
            val1 <- get.par(parname, repman, exp=TRUE)[indnext, 2]
            return(round((val1 - val)/val*100, 1))
        }
        indstart <- which(rep$inp$time == rep$inp$manstart)
        #indstart <- rep$inp$indpred[1]-1 # Current time (last time interval of last year)
        #curtime <- rep$inp$time[indstart+1]
        curtime <- rep$inp$manstart
        indnext <- which(rep$inp$time == curtime+ypred) # Current time + ypred
        if (length(indnext) == 1){
            Cn <- paste0('C')
            Bn <- paste0('B')
            Fn <- paste0('F')
            get.cn <- function(nn){
                nc <- nchar(nn)
                tl <- 7 # Total length
                # Add spaces
                #pad <- ifelse(include.unc, paste0(rep(' ', max(0, tl-nc)), collapse=''), '')
                pad <- ''
                return(c(paste0(nn, '.lo'), paste0(pad, nn), paste0(nn, '.hi')))
            }
            BBn <- paste0('BqBmsy') # Should use / instead of q, but / is not accepted in varnames
            FFn <- paste0('FqFmsy')
            EBinfBn <- paste0('EBinfqBmsy')

            scenarios <- attr(repman,"scenarios")
            nsc <- length( scenarios )
            Cnextyear <- matrix(0, nsc, 3)
            colnames(Cnextyear) <- get.cn(Cn)
            Bnextyear <- matrix(0, nsc, 3)
            colnames(Bnextyear) <- get.cn(Bn)
            Fnextyear <- matrix(0, nsc, 3)
            colnames(Fnextyear) <- get.cn(Fn)
            BBnextyear <- matrix(0, nsc, 3)
            colnames(BBnextyear) <- get.cn(BBn)
            FFnextyear <- matrix(0, nsc, 3)
            colnames(FFnextyear) <- get.cn(FFn)
            perc.dB <- numeric(nsc)
            perc.dF <- numeric(nsc)
            EBinf <- numeric(nsc)
            for(i in 1:nsc){
                rp <- repman[[ scenarios[i] ]]  ##repman[[i]]
                EBinf[i] <- get.EBinf(rp)
                perc.dB[i] <- get.pdelta(rep, rp, indstart, indnext, parname='logB')
                perc.dF[i] <- get.pdelta(rep, rp, indstart, indnext, parname='logF')
                indnextC <- which((rp$inp$timeCpred + rp$inp$dtcp) == curtime+ypred)
                Cnextyear[i, ] <- round(get.par('logCpred', rp, exp=TRUE)[indnextC, 1:3], 1)
                Bnextyear[i, ] <- round(get.par('logB', rp, exp=TRUE)[indnext, 1:3], 1)
                Fnextyear[i, ] <- round(get.par('logF', rp, exp=TRUE)[indnext, 1:3], 3)
                BBnextyear[i, ] <- round(get.par('logBBmsy', rp, exp=TRUE)[indnext, 1:3], 3)
                FFnextyear[i, ] <- round(get.par('logFFmsy', rp, exp=TRUE)[indnext, 1:3], 3)
            }
            indnextCrep <- which((rep$inp$timeCpred+rep$inp$dtcp) == curtime+ypred)
            FBtime <- fd(curtime+ypred)
            Ctime1 <- fd(rep$inp$timeCpred[indnextCrep])
            Ctime2 <- fd(rep$inp$timeCpred[indnextCrep]+rep$inp$dtcp[indnextCrep])
            if (!verbose){
                Cn <- paste0('C', Ctime1)
                Bn <- paste0('B', FBtime)
                Fn <- paste0('F', FBtime)
            }
            # Data frame with predictions
            df <- cbind(Cnextyear[, 2], Bnextyear[, 2], Fnextyear[, 2], BBnextyear[, 2],
                        FFnextyear[, 2], perc.dB, perc.dF)
            colnames(df)[1:5] <- c(Cn, Bn, Fn, BBn, FFn)
            qinds <- grep('q', colnames(df))
            colnames(df)[qinds] <- sub('q', '/', colnames(df)[qinds]) # Replace q with /
            # Data frame with uncertainties of absolute predictions
            inds <- c(1, 3)
            dfabs <- cbind(Cnextyear[, inds,drop=FALSE], Bnextyear[, inds,drop=FALSE], Fnextyear[, inds,drop=FALSE])
            colnames(dfabs) <- c(colnames(Cnextyear)[inds], colnames(Bnextyear)[inds],
                                 colnames(Fnextyear)[inds])
            # Data frame with uncertainties of relateive predictions
            dfrel <- cbind(BBnextyear[, inds,drop=FALSE], FFnextyear[, inds,drop=FALSE])
            colnames(dfrel) <- c(colnames(BBnextyear)[inds], colnames(FFnextyear)[inds])
            qinds <- grep('q', colnames(dfrel))
            colnames(dfrel)[qinds] <- sub('q', '/', colnames(dfrel)[qinds]) # Replace q with /
            # Set row names
            scenarios <- attr(repman, "scenarios")
            rn <- c('1. Keep current catch', '2. Keep current F', '3. Fish at Fmsy',
                    '4. No fishing', '5. Reduce F 25%', '6. Increase F 25%')[scenarios]
            
            rownames(df) <- rn
            rownames(dfrel) <- rn
            rownames(dfabs) <- rn
            #cat('Management summary\n')
            timerangeI <- range(unlist(rep$inp$timeI))
            timerangeC <- range(rep$inp$timeC)
            lastcatchseen <- tail(rep$inp$timeC+rep$inp$dtc, 1)
            # Start printing stuff
            if (verbose){ # Time interval information
                cat(paste0('Observed interval, index:  ',
                           fd(timerangeI[1]),
                           ' - ',
                           fd(timerangeI[2]),
                           '\n'))
                cat(paste0('Observed interval, catch:  ',
                           fd(timerangeC[1]),
                           ' - ',
                           fd(lastcatchseen),
                           '\n\n'))
                cat(paste0('Fishing mortality (F) prediction: ',
                           FBtime, '\n'))
                cat(paste0('Biomass (B) prediction:           ',
                           FBtime, '\n'))
                cat(paste0('Catch (C) prediction interval:    ',
                           Ctime1,
                           ' - ',
                           Ctime2,
                           '\n\n'))
                if (rep$inp$catchunit != ''){
                    cat(paste('Catch/biomass unit:', rep$inp$catchunit, '\n\n'))
                }
                cat('Predictions\n')
            }
            print(df)
            if (include.unc){
                cat('\n95% CIs of absolute predictions\n')
                print(dfabs)
                cat('\n95% CIs of relative predictions\n')
                print(dfrel)
            }
            invisible(df)
        } else {
            cat('Warning: Could not show management results because ypred is larger than the calculated management time frame. Reduce ypred or increase inp$timepredc and run fit.spict() and manage() again.\n')
        }
    }
}


#' @name pred.catch
#' @title Predict the catch of the prediction interval specified in inp
#' @param rep Result list as output from fit.spict().
#' @param fmsyfac Projection are made using F = fmsyfac * Fmsy.
#' @param ffac Projection are made using F = ffac * F_last.
#' @param get.sd Get uncertainty of the predicted catch.
#' @param exp If TRUE report exp of log predicted catch.
#' @param dbg Debug flag, dbg=1 some output, dbg=2 more ourput.
#' @return A vector containing predicted catch (possibly with uncertainty).
#' @export
pred.catch <- function(repin, fmsyfac=1, ffac=NULL, MSEmode=TRUE, get.sd=FALSE, exp=FALSE, dbg=0){
    inpin <- list()
    inpin$dteuler <- repin$inp$dteuler
    inpin$timeC <- repin$inp$timeC
    inpin$obsC <- repin$inp$obsC
    inpin$timeI <- repin$inp$timeI
    inpin$obsI <- repin$inp$obsI
    timelastobs <- repin$inp$time[repin$inp$indlastobs]
    # Always predict at least two years
    inpin$timepredc <- repin$inp$timepredc
    inpin$timepredi <- repin$inp$timepredi
    Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
    Flast <- get.par('logF', repin, exp=TRUE)[repin$inp$indpred[1], 2]
    if(!is.null(ffac) & fmsyfac == 1){
        fac <- ffac + 1e-6        
    }else{
        fac <- (fmsyfac + 1e-6) * Fmsy / Flast
    }    
    inpt <- check.inp(inpin)
    # Set F fac
    inpt <- make.ffacvec(inpt, fac)
    # Make object
    datint <- make.datin(inpt, dbg=dbg)
    datint$MSEmode <- MSEmode    
    plt <- repin$obj$env$parList(repin$opt$par)
    objt <- make.obj(datint, plt, inpt, phase=1)
    objt$fn(repin$opt$par)
    if (get.sd){
        repmant <- sdreport(objt)
        Cp <- get.par('logCp', repmant, exp=exp)
    } else {
        Cp <- c(NA, log(objt$report()$Cp), NA, NA, NA)
        if (exp){
            Cp <- exp(Cp)
        }
        names(Cp) <- names(get.par('logK', repin)) # Just to get the names
    }
    return(Cp)
}



#' @name get.TAC
#' 
#' @title Estimate Total Allowable Catch (TAC)
#' 
#' @param repin Result list as output from fit.spict().
#' @param hcr Harvest control rule. Options: 'msy', 'pa', 'dl', '2/3' (more
#'     information under details)
#' @param fractileC The fractile of the catch distribution to be used
#'     for setting the TAC. Default is median (0.5).
#' @param fractileFFmsy The fractile of the distribution of
#'     F/Fmsy. Default is 0.5 (median).
#' @param fractileBBmsy The fractile of the distribution of
#'     B/Bmsy. Default is 0.5 (median).
#' @param pa Logical; indicating if the precautionary approach should
#'     be applied (reduce F if P(B<Blim) < prob). Default is FALSE.
#' @param prob Probability for the precautionary approach (see
#'     argument 'pa', default is 0.95).
#' @param bfrac Fraction of biomass relativ to biomass reference
#'     levels (dependent on \code{quant}), e.g.  fraction of B/Bmsy
#'     which is defined as threshold (Blim = 0.3 Bmsy, Btrigger = 0.5
#'     Bmsy) or fraction of Bp/Bl
#' @param stab Logical; stability clause. If true F multiplication factor is bound between two values set in lower and upper. Default: FALSE.
#' @param lower Lower bound of the stability clause. Default is 0.8, used if uncertaintyCap = TRUE.
#' @param upper Upper bound of the stability clause. Default is 1.2, used if uncertaintyCap = TRUE.
#' @param tcv threshold for CV (for hcr dl3) when using pbb when msy
#' @param getFit Logical; if TRUE the fitted results list with adjusted fsihing mortality value is returned. Default is FALSE.
#'
#' @details The possible harvest control rules are:
#' \itemize{
#'   \item{"msy"}{Standard MSY approach}
#'   \item{"pa"}{Precautionary approach}
#'   \item{"dl"}{Data-limited approach}
#'   \item{"2/3"}{2 over 3 rule = last 2 observations of abundance index over preceeding 3}
#' }
#' 
#' @return A list with estimated TAC based on harvest control rule
#'     settings or the fitted rep list with adjusted fishing mortality
#'     values if getFit = TRUE and a logical value indicating if the
#'     stability clause was hit or not (if in use).
#' 
#' @export
#' 
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' get.TAC(rep)
get.TAC  <- function(repin,
                     hcr = "msy",
                     fractileC = 0.5,
                     fractileFFmsy = 0.5,
                     fractileBBmsy = 0.5,
                     prob = 0.95,
                     bfrac = 0.3,
                     babs = NA,
                     r23_type = "2/3",
                     r23_pa = FALSE,
                     r23_paRed = 0.2,
                     stab = FALSE,
                     lower = 0.8,
                     upper = 1.2,
                     tcv = 0.5,
                     getFit = FALSE,
                     MSEmode = 0){

    ## hack for DLMtool - The number of stochastic samples of the TAC recommendation. 
    reps = 1
    ## inp
    inpin <- repin$inp    
    ## stop if repin not converged
    if((is.null(repin) || is(repin, "try-error") || repin$opt$convergence != 0 ||
       any(is.infinite(repin$sd))) & hcr != "bmed")
        return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "spict"))
    if(hcr %in% c("msy","pa")){
        quant = "logBpBmsynotS"        
        ## get quantities
        logFpFmsy <- get.par("logFpFmsynotS", repin)
        logBpBmsy <- get.par("logBpBmsynotS",repin)
        Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
        Bmsy <- get.par('logBmsy', repin, exp=TRUE)[2]        
        Flast <- get.par('logFnotS', repin, exp=TRUE)[inpin$indpred[1],2]
##        Flast <- get.par('logF', repin, exp=TRUE)[inpin$indpred[1],2]
        ## second non-convergence stop
        if(any(is.null(c(logFpFmsy[2],logBpBmsy[2],Flast,Fmsy))) ||
           !all(is.finite(c(logFpFmsy[2],logBpBmsy[2],Flast,Fmsy))))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "msy/pa"))
            ## F multiplication factor based on uncertainty in F/Fmsy. Default = median        
            fi <- 1-fractileFFmsy
            fm <- exp(qnorm(fi, logFpFmsy[2], logFpFmsy[4]))
            fm5 <- exp(qnorm(0.5, logFpFmsy[2], logFpFmsy[4]))
            ## F multiplication factor based on uncertainty in B/Bmsy. Default = median
            bi <- 2 * exp(qnorm(fractileBBmsy, logBpBmsy[2], logBpBmsy[4]))
            fmult <- fm5 / fm * min(1, bi)         
            fabs <- (fmult + 1e-6) * Fmsy / Flast
            ## precautionary approach
            if(hcr == "pa"){
                repPA <- repin
                inpPA <- make.ffacvec(repPA$inp, fabs)
                repPA$obj$env$data$ffacvec <- inpPA$ffacvec
                repPA$obj$env$data$MSEmode <- 1
                repPA$obj$retape()
                repPA$obj$fn(repin$opt$par)
                sdr <- try(sdreport(repPA$obj),silent=TRUE)
                ## stop if not converged
                if(is.null(sdr) || is(sdr, "try-error"))
                    return(list(TAC=rep(NA, reps), hitSC=FALSE, conv = FALSE, id = "pa"))
                ## get quantities
                logBpBmsyPA <- get.par("logBpBmsynotS",sdr)
                ll <- qnorm(1-prob,logBpBmsyPA[2],logBpBmsyPA[4])
                bbmsyQ5 <- exp(ll)
                ## stop if not finite
                if(is.null(bbmsyQ5) || !is.finite(bbmsyQ5))
                    return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pa"))
                ## if bref as absolute value
                if(is.numeric(babs)){
                    bfrac <- babs / Bmsy
                }
                ## check if precautionary
                if((bbmsyQ5 - bfrac) < -1e-2){
                    tmp <- try(spict:::get.ffac(repin, bfrac=bfrac, prob=prob,
                                                quant=quant, MSEmode = 1))
                    if(is.null(tmp) || is(tmp, "try-error") || !is.finite(tmp))
                        return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pa"))
                    ## debugging:
                    ## if(tmp > fabs) print(paste0("ffacpa",round(tmp,2)," > ffacmsy",round(fabs,2)))
                    print(tmp)
                    fabs <- tmp
                    fmult <- fabs * Flast / Fmsy
                }
            }
        if(is.null(fmult) || !is.finite(fmult))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "msy/pa"))
        ## Stability clause or uncertainty cap (for dl and 2/3)
        if(stab){
            fmult <- spict:::stabilityClause(fmult, lower, upper)
            if(any(fmult < lower) || any(fmult > upper)) hitSC <- TRUE else hitSC <- FALSE
        }else hitSC <- FALSE
        ## convert back to f multiplication factor    
        fabs <- fmult * Fmsy / Flast
        ## predict catch with fabs
        if(is.null(fabs) || !is.finite(fabs))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "msy/pa"))
        TACi <- spict:::get.TACi(repin, fabs, fractileC)
        TAC <- rep(TACi, reps)
        id <- "msy/pa"
    }else if(hcr %in% c("dl")){
        quant = "logBpBl"
        ## get quantities
        logBpBl <- get.par("logBpBl", repin, exp = FALSE)
        logBBl <- get.par("logBBl", repin, exp = FALSE)
        Flast <- get.par('logFnotS', repin, exp=TRUE)[inpin$indpred[1],2]
        ## second non-convergence stop
        if(any(is.null(c(logBpBl[2],logBBl[,2],Flast))) ||
           !all(is.finite(c(logBpBl[2],logBBl[,2],Flast))))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
        ## TAC = Clast by default
        lastyearidxs <- min( which( cumsum(rev(inpin$dtc))>=1 ) ) ## warning: this will not make sense with subannual/mixed data with missing values
        TACi <- sum(tail(inpin$obsC, lastyearidxs))
        ## check if Bpred >= Blast at least x%
        ll <- qnorm(1-prob,logBpBl[2],logBpBl[4])
        bpblQx <- exp(ll)
        ## stop if not finite
        if(is.null(bpblQx) || !is.finite(bpblQx))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
        ## if smaller
        if(abs(bpblQx - bfrac) > 1e-3){
            tmp <- try(spict:::get.ffac(repin, bfrac=bfrac, prob=prob,
                                        quant=quant, MSEmode = 2))
            if(is.null(tmp) || is(tmp, "try-error") || !is.finite(tmp))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
            fabs <- tmp
            fmult <- fabs / Flast
            if(is.null(fmult) || !is.finite(fmult))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
            ## Stability clause or uncertainty cap (for dl and 2/3)
            if(stab){
                fmult <- spict:::stabilityClause(fmult, lower, upper)
                if(any(fmult < lower) || any(fmult > upper)) hitSC <- TRUE else hitSC <- FALSE
            }else hitSC <- FALSE
            ## convert back to f multiplication factor    
            fabs <- fmult * Flast
            ## predict catch with fabs
            if(is.null(fabs) || !is.finite(fabs))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
            TACi <- spict:::get.TACi(repin, fabs, fractileC)
            TAC <- rep(TACi, reps)
        }else{ ## otherwise keep current Catch
            TAC <- rep(TACi, reps)
            hitSC <- FALSE
        }
        id <- "pbb"
    }else if(hcr %in% c("dl2")){
        quant = "logBpBmsy"        
        ## get quantities
        logFpFmsy <- get.par("logFpFmsy", repin)
        logBpBmsy <- get.par("logBpBmsy",repin)
        Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
        Flast <- get.par('logFnotS', repin, exp=TRUE)[inpin$indpred[1],2]
        ## F multiplication factor based on uncertainty in F/Fmsy. Default = median        
        fi <- 1-fractileFFmsy
        fm <- exp(qnorm(fi, logFpFmsy[2], logFpFmsy[4]))
        fm5 <- exp(qnorm(0.5, logFpFmsy[2], logFpFmsy[4]))
        fmult <- fm5 / fm
        fabs <- (fmult + 1e-6) * Fmsy / Flast
        ## check if Bpred >= Bmsy at least x%
        ll <- qnorm(1-prob,logBpBmsy[2],logBpBmsy[4])
        bpblQx <- exp(ll)
        ## if smaller
        if((bpblQx - bfrac) < -1e-3){
            tmp <- try(spict:::get.ffac(repin, bfrac=bfrac, prob=prob,
                                        quant=quant, MSEmode = 1))
            if(is.null(tmp) || is(tmp, "try-error") || !is.finite(tmp))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
            fabs <- tmp
            fmult <- fabs * Flast / Fmsy            
            if(is.null(fmult) || !is.finite(fmult))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
        }
        ## Stability clause or uncertainty cap (for dl and 2/3)
        if(stab){
            fmult <- spict:::stabilityClause(fmult, lower, upper)
            if(any(fmult < lower) || any(fmult > upper)) hitSC <- TRUE else hitSC <- FALSE
        }else hitSC <- FALSE
        ## convert back to f multiplication factor
        fabs <- fmult * Fmsy / Flast                
        ## predict catch with fabs
        if(is.null(fabs) || !is.finite(fabs))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
        TACi <- spict:::get.TACi(repin, fabs, fractileC)
        TAC <- rep(TACi, reps)
        id <- "pbb"

    }else if(hcr %in% c("dl3")){
        quant = "logBpBl"
        ## get quantities
        logBpBl <- get.par("logBpBl", repin, exp = FALSE)
        logBBl <- get.par("logBBl", repin, exp = FALSE)
        Flast <- get.par('logFnotS', repin, exp=TRUE)[inpin$indpred[1],2]
        ## second non-convergence stop
        if(any(is.null(c(logBpBl[2],logBBl[,2]))) ||
           !all(is.finite(c(logBpBl[2],logBBl[,2]))))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
        ## TAC = Clast by default
        lastyearidxs <- min( which( cumsum(rev(inpin$dtc))>=1 ) ) ## warning: this will not make sense with subannual/mixed data with missing values
        TACi <- sum(tail(inpin$obsC, lastyearidxs))
        ## check CV of rel ref levels
        flfmsy <- get.par("logFlFmsy", repin, exp=TRUE)
        blbmsy <- get.par("logBlBmsy", repin, exp=TRUE)
        hitSC <- FALSE
        if(all(c(flfmsy[5], blbmsy[5]) < tcv)){
            quant = "logBpBmsy"        
            ## get quantities
            logFpFmsy <- get.par("logFpFmsy", repin)
            logBpBmsy <- get.par("logBpBmsy",repin)
            Fmsy <- get.par('logFmsy', repin, exp=TRUE)[2]
            Flast <- get.par('logFnotS', repin, exp=TRUE)[inpin$indpred[1],2]
            ## second non-convergence stop
            if(any(is.null(c(logFpFmsy[2],logBpBmsy[2],Flast,Fmsy))) ||
               !all(is.finite(c(logFpFmsy[2],logBpBmsy[2],Flast,Fmsy))))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "msy"))
            ## F multiplication factor based on uncertainty in F/Fmsy. Default = median        
            fi <- 1-fractileFFmsy
            fm <- exp(qnorm(fi, logFpFmsy[2], logFpFmsy[4]))
            fm5 <- exp(qnorm(0.5, logFpFmsy[2], logFpFmsy[4]))
            ## F multiplication factor based on uncertainty in B/Bmsy. Default = median
            bi <- 2 * exp(qnorm(fractileBBmsy, logBpBmsy[2], logBpBmsy[4]))
            fmult <- fm5 / fm * min(1, bi)         
            fabs <- (fmult + 1e-6) * Fmsy / Flast
            if(is.null(fmult) || !is.finite(fmult))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "msy"))
            ## Stability clause or uncertainty cap (for dl and 2/3)
            if(stab){
                fmult <- spict:::stabilityClause(fmult, lower, upper)
                if(any(fmult < lower) || any(fmult > upper)) hitSC <- TRUE else hitSC <- FALSE
            }else hitSC <- FALSE
            ## convert back to f multiplication factor    
            fabs <- fmult * Fmsy / Flast
            ## predict catch with fabs
            if(is.null(fabs) || !is.finite(fabs))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "msy"))
            TACi <- spict:::get.TACi(repin, fabs, fractileC)
            TAC <- rep(TACi, reps)
            id <- "msy"
        }else{
            ## check if Bpred >= Blast at least x%
            ll <- qnorm(1-prob,logBpBl[2],logBpBl[4])
            bpblQx <- exp(ll)
            ## stop if not finite
            if(is.null(bpblQx) || !is.finite(bpblQx))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
            ## if smaller
            if(abs(bpblQx - bfrac) > 1e-2){
                tmp <- try(spict:::get.ffac(repin, bfrac=bfrac, prob=prob,
                                            quant=quant, MSEmode = 2))
                if(is.null(tmp) || is(tmp, "try-error") || !is.finite(tmp))
                    return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
                fabs <- tmp
                fmult <- fabs / Flast
                if(is.null(fmult) || !is.finite(fmult))
                    return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
                ## Stability clause or uncertainty cap (for dl and 2/3)
                if(stab){
                    fmult <- spict:::stabilityClause(fmult, lower, upper)
                    if(any(fmult < lower) || any(fmult > upper)) hitSC <- TRUE else hitSC <- FALSE
                }else hitSC <- FALSE
                ## convert back to f multiplication factor    
                fabs <- fmult * Flast
                ## predict catch with fabs
                if(is.null(fabs) || !is.finite(fabs))
                    return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))    
                TACi <- spict:::get.TACi(repin, fabs, fractileC)
                TAC <- rep(TACi, reps)
            }else{ ## otherwise keep current Catch
                TAC <- rep(TACi, reps)
                hitSC <- FALSE
            }
            id <- "pbb"
        }
    }else if(hcr %in% c("bmed")){
        quant = "logBpBl"
        ## get quantities
        logBpBl <- get.par("logBpBl", repin, exp = FALSE)
        Flast <- get.par('logFnotS', repin, exp=TRUE)[inpin$indpred[1],2]
        ## second non-convergence stop
        if(any(is.null(c(logBpBl[2],Flast))) ||
           !all(is.finite(c(logBpBl[2],Flast))))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "bmed"))
        ## TAC = Clast by default
        lastyearidxs <- min( which( cumsum(rev(inpin$dtc))>=1 ) ) ## warning: this will not make sense with subannual/mixed data with missing values
        TACi <- sum(tail(inpin$obsC, lastyearidxs))
        ## check if Bpred >= Blast at least x%
        ll <- qnorm(1-0.5,logBpBl[2],logBpBl[4])
        bpblQx <- exp(ll)
        fabs <- Flast
        fmult <- fabs/Flast
        ## stop if not finite
        if(is.null(bpblQx) || !is.finite(bpblQx))
            return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
        ## if smaller
        if(abs(bpblQx - bfrac) > 1e-3){
            tmp <- Flast * 0.75
            fabs <- tmp
            fmult <- fabs / Flast
            if(is.null(fmult) || !is.finite(fmult))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
            ## Stability clause or uncertainty cap (for dl and 2/3)
            if(stab){
                fmult <- spict:::stabilityClause(fmult, lower, upper)
                if(any(fmult < lower) || any(fmult > upper)) hitSC <- TRUE else hitSC <- FALSE
            }else hitSC <- FALSE
            ## convert back to f multiplication factor    
            fabs <- fmult * Flast
            ## predict catch with fabs
            if(is.null(fabs) || !is.finite(fabs))
                return(list(TAC=rep(NA, reps),hitSC=FALSE, conv = FALSE, id = "pbb"))
            TACi <- spict:::get.TACi(repin, fabs, fractileC)
            TAC <- rep(TACi, reps)
        }else{ ## otherwise keep current Catch
            TAC <- rep(TACi, reps)
            hitSC <- FALSE
        }
        id <- "bmed"        
    }else if(hcr %in% c("2/3")){
        ## get quantities
        inds <- inpin$obsI
        if(length(inds) > 1){
            ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
            indtab <- do.call(rbind, inds)
            ind <- apply(indtab, 2, mean)
        }else{
            ind <- unlist(inds)
        }
        ninds <- length(ind)
        r23t <- as.numeric(unlist(strsplit(r23_type, "/")))
        inum <- ind[(ninds-(r23t[1]-1)):ninds]
        iden <- ind[(ninds-(r23t[1]+r23t[2]-1)):(ninds-r23t[1])]
        ## inum <- ind[(ninds-1):ninds]
        ## iden <- ind[(ninds-4):(ninds-2)]
        r23 <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
        ## uncertainty cap
        if(stab){
            r23 <- spict:::stabilityClause(r23, lower, upper)
            if(any(r23 < lower) || any(r23 > upper)) hitSC <- TRUE else hitSC <- FALSE
        }else hitSC <- FALSE
        ## account for seasonal and annual catches
        Cl <- sum(tail(inpin$obsC, tail(1/inpin$dtc,1)))
        TACi <- Cl * r23 * 1 * 1  ## Clast * r * f * b
        ## pa buffer
        if(r23_pa) TACi <- TACi * r23_paRed
        TAC <- rep(TACi, reps)        
        if(getFit){
            fit <- try(take.c(catch = Cl, inpin = inpin, repin = repin),silent=TRUE)
            if(!is(fit,"try-error")) return(fit)            
        }
        return(list(TAC=TAC, hitSC=hitSC, conv = NA, id="23"))        
    }
    ## get fitted object
    if(getFit){
        inpt <- make.ffacvec(repin$inp, fmult)
        inpt$MSEmode <- MSEmode
        fit <- try(fit.spict(inpt),silent=TRUE)
        if(!is(fit,"try-error")) return(fit)
    }
    reslist <- list(TAC=TAC, hitSC=hitSC, conv = TRUE, id = id, fabs = fabs, fmult = fmult)
    return(reslist)
}
