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


tws.plot.niche.dyn <- 
function (z1, z2, intersection = 0, title = "", name.axis1 = "Axis 1", 
    name.axis2 = "Axis 2", interest = 1, col.abn = "lightgreen", 
    col.unf = "green", col.exp = "red", col.stab = "blue", col.pio = "pink", 
    col.NA = "grey", colZ1 = "green3", colZ2 = "red3", transparency = 70, 
    ...) 
{
    t_col <- function(color, percent = 50, name = NULL) {
        rgb.val <- col2rgb(color)
        t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], alpha = (100 - 
            percent) * 255/100, names = name, maxColorValue = 255)
    }
    col.abn = t_col(col.abn, transparency)
    col.unf = t_col(col.unf, transparency)
    col.exp = t_col(col.exp, transparency)
    col.stab = t_col(col.stab, transparency)
    col.pio = t_col(col.pio, transparency)
    col.NA = t_col(col.NA, transparency)
    cat <- ecospat.niche.dyn.index(z1, z2, intersection)$dyn
    if (is.null(z1$y)) {
        R <- length(z1$x)
        x <- z1$x
        xx <- sort(rep(1:length(x), 2))
        y1 <- z1$z.uncor/max(z1$z.uncor)
        Y1 <- z1$Z/max(z1$Z)
        if (intersection > 0) {
            Y1.quant <- quantile(as.matrix(z1$Z)[which(as.matrix(z1$Z) > 
                0)], probs = seq(0, 1, intersection))[2]/max(as.matrix(z1$Z))
        }
        else {
            Y1.quant <- 0
        }
        Y1.quant <- Y1 - Y1.quant
        Y1.quant[Y1.quant < 0] <- 0
        yy1 <- sort(rep(1:length(y1), 2))[-c(1:2, length(y1) * 
            2)]
        YY1 <- sort(rep(1:length(Y1), 2))[-c(1:2, length(Y1) * 
            2)]
        y2 <- z2$z.uncor/max(z2$z.uncor)
        Y2 <- z2$Z/max(z2$Z)
        if (intersection > 0) {
            Y2.quant <- quantile(as.matrix(z2$Z)[which(as.matrix(z2$Z) > 
                0)], probs = seq(0, 1, intersection))[2]/max(as.matrix(z2$Z))
        }
        else {
            Y2.quant <- 0
        }
        Y2.quant <- Y2 - Y2.quant
        Y2.quant[Y2.quant < 0] <- 0
        yy2 <- sort(rep(1:length(y2), 2))[-c(1:2, length(y2) * 
            2)]
        YY2 <- sort(rep(1:length(Y2), 2))[-c(1:2, length(Y2) * 
            2)]
        plot(x, y1, type = "n", xlab = name.axis1, ylab = "density of occurrences", 
            main = title)
        polygon(x[xx], c(0, y1[yy1], 0, 0), col = col.unf, border = 0)
        polygon(x[xx], c(0, y2[yy2], 0, 0), col = col.exp, border = 0)
        polygon(x[xx], c(0, apply(cbind(y2[yy2], y1[yy1]), 1, 
            min, na.exclude = TRUE), 0, 0), col = col.stab, border = 0)
        lines(x[xx], c(0, Y2.quant[YY2], 0, 0), col = colZ2, 
            lty = "dashed")
        lines(x[xx], c(0, Y1.quant[YY1], 0, 0), col = colZ1, 
            lty = "dashed")
        lines(x[xx], c(0, Y2[YY2], 0, 0), col = colZ2)
        lines(x[xx], c(0, Y1[YY1], 0, 0), col = colZ1)
        segments(x0 = min(x[xx]), y0 = 0, x1 = max(x[xx]), y1 = 0, 
            col = "white")
        seg.cat <- function(inter, cat, col.abn, col.unf, col.stab, 
            col.exp, col.pio, col.NA) {
            if (inter[3] == 0) {
                my.col <- 0
            }
            if (inter[3] == 1) {
                my.col <- col.abn
            }
            if (inter[3] == 2) {
                my.col <- col.unf
            }
            if (inter[3] == 3) {
                my.col <- col.stab
            }
            if (inter[3] == 4) {
                my.col <- col.exp
            }
            if (inter[3] == 5) {
                my.col <- col.pio
            }
            if (inter[3] == 6) {
                my.col <- col.NA
            }
            ##segments(x0 = inter[1], y0 = -0.01, y1 = -0.01, x1 = inter[2], 
            ##    col = my.col, lwd = 4, lty = 2)
        }
        inter <- cbind(z1$x[-length(z1$x)], z1$x[-1], cat[-1][])
        apply(inter, 1, seg.cat, col.unf = col.unf, col.exp = col.exp, 
            col.stab = col.stab, col.pio = col.pio, col.abn = col.abn, 
            col.NA = col.NA)
    }
    if (!is.null(z1$y)) {
        col_category <- c("#FFFFFF", col.abn, col.unf, col.stab, 
            col.exp, col.pio, col.NA)[sort(1 + (unique(terra::values(cat))))]
        if (interest == 1) {
            terra::plot(z1$z.uncor, col = gray(100:0/100), legend = FALSE, 
                xlab = name.axis1, ylab = name.axis2, mar = c(3.1, 
                  3.1, 2.1, 3.1))
        }
        if (interest == 2) {
            terra::plot(z2$z.uncor, col = gray(100:0/100), legend = FALSE, 
                xlab = name.axis1, ylab = name.axis2, mar = c(3.1, 
                  3.1, 2.1, 3.1))
        }
        if (interest == 0) {
            terra::plot(cat, col = col_category, legend = FALSE, 
                box = TRUE, xlab = name.axis1, ylab = name.axis2, 
                mar = c(3.1, 3.1, 2.1, 3.1))
        }
        else {
            terra::plot(cat, col = col_category, add = TRUE, 
                legend = FALSE, box = TRUE)
        }
        title(title)
        terra::contour(z1$Z, add = TRUE, levels = quantile(z1$Z[z1$Z > 
            0], c(0, intersection)), drawlabels = FALSE, lty = c(1, 
            2), col = colZ1)
        terra::contour(z2$Z, add = TRUE, levels = quantile(z2$Z[z2$Z > 
            0], c(0, intersection)), drawlabels = FALSE, lty = c(1, 
            2), col = colZ2)
    }
}
