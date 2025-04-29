twsBoyce <- function (fit, obs, nclass = 0, window.w = "default", res = 100, 
    PEplot = TRUE, rm.duplicate = TRUE, method = "spearman") 
{

    boycei <- function(interval, obs, fit) {
      ## interval is a vector of the window min and max values
      ## obs is the suitability of all observations
      ## fit is the suitability of all cells

      ## message(round(interval[1], 3), " ", round(interval[2], 3))

      ## pi is the proportion of all observations in the window
      pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2]))/length(obs)
      ## ei is the proportion of all cells in the window (i.e., expected)
      ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2]))/length(fit)

      ## message("PI: ", pi)
      ## message("EI: ", ei)

      return(round(pi/ei, 10))
    }
    if (inherits(fit, "SpatRaster")) {
        if (is.data.frame(obs) || is.matrix(obs)) {
          ## extract fitness values for observations as a vector
          obs <- terra::extract(fit, as.data.frame(obs), ID = FALSE)
          obs <- as.numeric(obs[, 1])
          ## TWS added: drop missing values
          obs <- obs[!is.na(obs)]
        }
        ## extract background fitness values as a vector
        fit <- terra::values(fit, na.rm = T)
        fit <- as.numeric(fit)
    }
    ## define min and max values
    mini <- min(fit, obs)
    maxi <- max(fit, obs)
    if (length(nclass) == 1) {
        if (nclass == 0) {
            if (window.w == "default") {
              ## default to 1/10 the total range
              window.w <- (max(fit) - min(fit))/10
            }
            ## define intervals as res windows spanning mini to maxi
            vec.mov <- seq(from = mini, to = maxi - window.w, length = res)
            interval <- cbind(vec.mov, vec.mov + window.w)
            ## vec.mov <- seq(from = mini, to = maxi - window.w, 
            ##     by = (maxi - mini - window.w)/res)
            ## vec.mov[res + 1] <- vec.mov[res + 1] + 1
            ## interval <- cbind(vec.mov, vec.mov + window.w)
        }
        else {
            vec.mov <- seq(from = mini, to = maxi - window.w, length = nclass)
            interval <- cbind(vec.mov, vec.mov + window.w)
        }
    }
    else { ## user provided multiple values for nclass - not documented?
      ## Looks like the user can use a vector of values to specify the
      ## window boundaries explicitly.
        vec.mov <- c(mini, sort(nclass[!nclass > maxi | nclass < 
            mini]))
        interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
    f <- apply(interval, 1, boycei, obs, fit)
    to.keep <- which(f != "NaN")
    f <- f[to.keep]
    if (length(f) < 2) {
      ## when would this happen?
      b <- NA
    }
    else {
        r <- 1:length(f)
        if (rm.duplicate == TRUE) {
            r <- c(1:length(f))[f != c(f[-1], TRUE)]
        }
        b <- cor(f[r], vec.mov[to.keep][r], method = method)
    }
    HS <- apply(interval, 1, sum)/2 ## midpoint of window
    HS <- HS[to.keep]
    if (PEplot == TRUE) {
        plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", 
            col = "grey", cex = 0.75)
        points(HS[r], f[r], pch = 19, cex = 0.75)
    }
    return(list(F.ratio = f, cor = round(b, 3), HS = HS, interval = interval))
}


## origBoyce <- function (fit, obs, nclass = 0, window.w = "default", res = 100, 
##     PEplot = TRUE, rm.duplicate = TRUE, method = "spearman") {
##     boycei <- function(interval, obs, fit) {
##       message(round(interval[1], 3), " ", round(interval[2], 3))
##       pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2]))/length(obs)
##       ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2]))/length(fit)
##       message("PI: ", pi)
##       message("EI: ", ei)

##       return(round(pi/ei, 10))
##     }
##     if (inherits(fit, "SpatRaster")) {
##         if (is.data.frame(obs) || is.matrix(obs)) {
##             obs <- terra::extract(fit, as.data.frame(obs), ID = FALSE)
##             obs <- as.numeric(obs[, 1])
##         }
##         fit <- terra::values(fit, na.rm = T)
##         fit <- as.numeric(fit)
##     }
##     mini <- min(fit, obs)
##     maxi <- max(fit, obs)
##     if (length(nclass) == 1) {
##         if (nclass == 0) {
##             if (window.w == "default") {
##                 window.w <- (max(fit) - min(fit))/10
##             }
##             vec.mov <- seq(from = mini, to = maxi - window.w, 
##                 by = (maxi - mini - window.w)/res)
##             vec.mov[res + 1] <- vec.mov[res + 1] + 1
##             interval <- cbind(vec.mov, vec.mov + window.w)
##         }
##         else {
##             vec.mov <- seq(from = mini, to = maxi, by = (maxi - 
##                 mini)/nclass)
##             interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
##         }
##     }
##     else {
##         vec.mov <- c(mini, sort(nclass[!nclass > maxi | nclass < 
##             mini]))
##         interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
##     }
##     f <- apply(interval, 1, boycei, obs, fit)
##     to.keep <- which(f != "NaN")
##     f <- f[to.keep]
##     if (length(f) < 2) {
##         b <- NA
##     }
##     else {
##         r <- 1:length(f)
##         if (rm.duplicate == TRUE) {
##             r <- c(1:length(f))[f != c(f[-1], TRUE)]
##         }
##         b <- cor(f[r], vec.mov[to.keep][r], method = method)
##     }
##     HS <- apply(interval, 1, sum)/2
##     if (length(nclass) == 1 & nclass == 0) {
##         HS[length(HS)] <- HS[length(HS)] - 1
##     }
##     HS <- HS[to.keep]
##     if (PEplot == TRUE) {
##         plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", 
##             col = "grey", cex = 0.75)
##         points(HS[r], f[r], pch = 19, cex = 0.75)
##         arrows(interval[, 1], 0.5, interval[, 2], 0.5)
##     }
##     return(list(F.ratio = f, cor = round(b, 3), HS = HS))
## }


