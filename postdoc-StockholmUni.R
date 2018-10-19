##########################################################################
# These functions are 
# Copyrigth (C) 2014-2018 V. Miranda, University of Auckland
# All rights reserved. 
######## Edited by Victor Miranda, University of Auckland, 2018-10-14

library("VGAM")
library("VGAMextra")
library("forecast")

##########################################################################
############# Functions implemented for this task ########################

auto.arima.MSE <- function(x = NULL, which.to.predict = NULL, 
                           cross.Valid = FALSE,
                           max.p = 3, max.q = 3,
                           max.d = 1, max.w = 15) {
  
  if (!is.null(dim(x)))
    stop("'x' must be a vector with a univariate time series.")
  
  if ( (x %>% length) && !is.ts(x))
    stop("Wrong input for argument 'x'. Must be a univariate", 
         "\n", "object of class 'ts'.")
  
  if (!is.character(which.to.predict))
    stop("Wrong input for argument 'which.to.predict'.", "\n",
         "Must be a char-string as 'YYYY-MM-DD'.")
  
  if ((c(max.p, max.q) %>% sum) > max.w)
    stop("Function currently implemented for window sizes, 'max.w', ",
         "\n", "larger than 'max.p + max.q'.")
  
  which.to.bis <- which.to.predict %>% as.Date("%Y-%m-%d") %>% as.numeric
  my.Dates <- c("2015-09-01", "2018-08-31") %>% as.Date("%Y-%m-%d")
 length.dat <- which.to.bis - my.Dates %>% min %>% as.numeric
  
  ######### Not executed
  if (FALSE) {
  read.bitcoinDF  <- 
    read.csv("coindesk-bpi-USD-close_data-2015-09-01_2018-08-31.csv")
  
  my.Dates  <-  as.Date(gsub(" 00:00:00", "",  ## Remove characters 00:00:00
                             read.bitcoinDF[, 1, drop = TRUE]), "%Y-%m-%d")
  my.Labels <- as.Date(seq(as.numeric(min(my.Dates)), 
                           as.numeric(max(my.Dates)), by = 90),
                       origin = "1970-01-01")
  #my.Start  <- as.numeric(min(my.Dates))        ## Origin for bitcoinUS
  bitcoinDF <- data.frame(Date = my.Dates,
                          logbitcoin = ts(log(read.bitcoinDF[, 2]),
                                  start = c(2015, 9), frequency = 365))
  if (is.null(x))
    x <- bitcoinDF[, 2, drop = TRUE]
  } ####### Not executed.
  
                 ### 'min.k' may be an argument
  min.k <- 10    ### Minimum data-length for fitting a model 
  if (which.to.bis < as.numeric(min(my.Dates)) + min.k)
    stop("At least 10 observations are needed to estimate arima models",
         "\n", "Enter a valid date at 'which.to.predict'.")
  
  #crit  <- "MSE"    ## May be an argument later
  series.freq  <- x %>% frequency
  serieslength <- x %>% length
  max.p <- c(max.p, series.freq) %>% min
  max.q <- c(max.q, series.freq) %>% min
  
  bestfit <- search.arima.MSE(x = x,
                              which.to.predict = which.to.bis,
                              cross.Valid = cross.Valid,
                              p = max.p, d = max.d,
                              q = max.q, max.w = max.w,
                              min.k = min.k,
                              extra = list(my.Dates = my.Dates))
  
  c("\n\n", "The best prediction for the date", which.to.predict, 
      "based on the MSE is: ", "\n") %>% cat
  my.mat <- matrix(c(bestfit$pre.forecast$mean %>% rep(2), 
                     bestfit$pre.forecast$lower, 
                     bestfit$pre.forecast$upper), 
                   3, 2, by = TRUE) + c(x)[length.dat]
    
  rownames(my.mat) <- c("Predicted", "Lower CI", "Upper CI")
  colnames(my.mat) <- c("80% CI", "95% CI")
  print(my.mat)
  c("\n", "At this window size, the minimum MSE is:", 
    bestfit$mse.min %>% round(8), "\n") %>% cat()
  
  #c("\n", "The observed time--point is",c(x)[] ,".\n") %>% cat
  c("\n", "The order of the corresponding arima model is p =",
    bestfit$order.f[1], ", d =", bestfit$fin.d, ", and q =",
    bestfit$order.f[1], ". \n") %>% cat()
  
  c("\n", "The 'optimal' window size obtained is:", 
    bestfit$fin.window, ".\n") %>% cat()
  
  return(invisible())
}


##
## Search for the best ARIMA moodel based on the MSE.
## The series NOT differenced.
##

search.arima.MSE <- function(x, which.to.predict, cross.Valid = FALSE,
                             p, d, q, max.w, min.k, extra = list()) {
  
  xArima <- x
  to.pred <- TRUE + which.to.predict - as.numeric(min(extra$my.Dates)) 
  out.put <- vector("list", 1)
  xArima <- xArima[1:to.pred]
  tt <- if (length(xArima)) length(xArima) else stop("Enter 'x'.")
  
  for (ii in max.w:min.k) {
    
    x.sub.set <- xArima[(tt - max.w + 1):tt]
    ## Max of 'max.d' non-seasonal differences allowed.
    fin.d <- 0    ## Change for D (later stage)
    if (!is.na(fin.d)) 
      fin.d <- ndiffs(x.sub.set, test = "kpss", max.d = d)
    
    if (fin.d) {
      if (length(cross.Valid) && cross.Valid)
        cross.data <- diff(x[to.pred:(to.pred + max.w - 1)])
      
      if (d >= 2)
        stop("Function currently implemented for series that require one ",
             "\n", "non-seasonal difference to stabilize.")
      cat("Window -- loop w = ", ii, "\n")
      xdiff <- diff(x.sub.set, lag = 1)
      for (dd in 1:fin.d) 
        bestfit <- search.arma.MSE(xarma = xdiff, p = p, q = q,
                                   cross.Valid = cross.Valid,
                                   extra = list(cross.data = cross.data))
    } else {
      warning("Fuction need to be tested yet when d = 0.", "\n")
      if (length(cross.Valid) && cross.Valid)
        cross.data <- x[to.pred:(to.pred + max.w - 1)]
      
      if (d >= 2)
        stop("Currently implemented for series that require one ",
             "\n", "non-seasonal difference to stabilize")
      cat("Window -- loop w = ", ii, "\n")
      for (dd in 1:fin.d) 
        bestfit <- search.arma.MSE(xarma = x.sub.set, p = p, q = q,
                                   cross.Valid = cross.Valid,
                                   extra = list(cross.data = cross.data))
    }
    
    if (ii == max.w) {
      out.put <- bestfit
    } else {
      if (bestfit$mse.min < out.put$mse.min) {
        out.put <- bestfit
        out.put$fin.window <- ii
      }
    }
    xArima <- c(xArima[-1],
                c(bestfit$pre.forecast$mean) + xArima[to.pred - FALSE])
  }
  
  out.put$fin.d <- fin.d
  return(out.put)
  
}



##
## Search for the best ARMA moodel based on the MSE.
## The series may be differenced.
##
search.arma.MSE <- function(xarma, p, q, cross.Valid = FALSE,
                            extra = list()) {
  
  xArma <- xarma; rm(xarma)
  n <- length(xArma)
  comb.in <- my.comb(x = p, y = q)
  all.mse <- apply(comb.in, 1, function(x) {
    sum(Arima(y = xArma, include.mean = TRUE,
              order = c(x[1], 0, x[2]))$residuals^2)/(n - (p + q + 2))
  })
  
  
  min.test <- which(all.mse == min(all.mse))[1]
  c("minimum MSE =", min(all.mse) %>% round(5), "\n\n") %>% cat
  
  best.pred.in <- Arima(xArma, include.mean = TRUE,
                        order = c(comb.in[min.test, ][1], 0,
                                  comb.in[min.test, ][2]))
  
  pre.forecast <- forecast(best.pred.in, h = 1)
  
  if (cross.Valid) {
    cross.data <- extra$cross.data
    save.pred <- c(forecast(best.pred.in, h = 1)$mean)
    xArma <- c(xArma[-1], c(forecast(best.pred.in, h = 1)$mean))
    for (kk in 1:(n - 1)) {
      fit.cross <- Arima(y = xArma, include.mean = TRUE,
                         order = c(comb.in[min.test, ][1], 0, 
                                   comb.in[min.test, ][2])) 
      save.pred <- c(save.pred, c(forecast(fit.cross, h = 1)$mean))
      xArma <- c(xArma[-1], c(forecast(fit.cross, h = 1)$mean))
    }
    fin.mse <- sum( (cross.data - save.pred)^2)/(length(cross.data))
  } else {
    fin.mse <- NULL
  }
  
  structure(list(order.f = comb.in[min.test, ],
                 pre.forecast = pre.forecast,
                 mse.min = min(all.mse), 
                 fin.mse = fin.mse))
}



my.comb <- function(x, y) {
  to.ret <- NULL
    for (ii in 1:x) {
      for (jj in 1:y)
        to.ret <- rbind(to.ret, c(ii, jj))
    }
  to.ret
}





######### Ignore.

if (FALSE) {
## Example.
read.bitcoinDF  <- 
  read.csv("coindesk-bpi-USD-close_data-2015-09-01_2018-08-31.csv")

my.Dates  <-  as.Date(gsub(" 00:00:00", "",  ## Remove characters 00:00:00
                           read.bitcoinDF[, 1, drop = TRUE]), "%Y-%m-%d")
my.Labels <- as.Date(seq(as.numeric(min(my.Dates)), 
                         as.numeric(max(my.Dates)), by = 90),
                     origin = "1970-01-01")
#my.Start  <- as.numeric(min(my.Dates))        ## Origin for bitcoinUS
bitcoinDF <- data.frame(Date = my.Dates,
                        logbitcoin = ts(log(read.bitcoinDF[, 2]),
                                        start = c(2015, 9), frequency = 365))
## The series of interest.
x <- bitcoinDF[, 2, drop = TRUE]
}