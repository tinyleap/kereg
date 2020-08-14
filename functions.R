getspec <- function(serverid, total) {
  dist.ids <- rep(serverid, length.out=total)
  machineAddresses <- list()
  for(i in serverid){
    machineAddresses[[i]] <- list(host=paste0('circinus-',i,'.ics.uci.edu'),user='ywang47@ics.uci.edu', ncore=sum(dist.ids==i))
  }
  spec <- lapply(machineAddresses[!sapply(machineAddresses, is.null)], function(machine) {rep(list(list(host=machine$host, user=machine$user)), machine$ncore)})
  spec <- unlist(spec,recursive=FALSE)
  spec
}
## Input: p*n matrices, where usually p=1,2
bound_gauss <- function(x) {
  temp <- colSums(x^2)
  exp(-temp/2)*(temp<=qchisq(0.95, dim(x)[1]))
}
gauss <- function(x) {
  temp <- colSums(x^2)
  exp(-temp/2)
}
# gauss <- function(x)mvtnorm::dmvnorm(t(x))
# bikereg.est <- function(X, Y, time, res, t, s, H, kernel=bound_gauss) {
#   kx <- kernel(solve(H, rbind(time-t, res-s)))*X
#   id <- complete.cases(kx)
#   xtx <- t(kx[id,])%*%X[id,]
#   if (qr(xtx)$rank<dim(xtx)[1]) rep(NA, dim(xtx)[1])
#   else solve(xtx, t(kx[id,])%*%Y[id])
# }
# bikereg.cv <- function(X, Y, time, res, H, kernel=bound_gauss, nfolds=5, cl=NULL) {
#   id <- !is.na(res)
#   coords <- paste(time, res, sep=',')
#   uniq.coords <- unique(coords[id])
#   match.coords <- match(coords, uniq.coords)
#   if (!is.list(H) || length(H)==0) stop('H must be a list!')
#   set.seed(1)
#   folds <- split(sample.int(length(Y)), rep(1:nfolds, length.out=length(Y)))
#   foo <- function(x, i, j) {
#     x <- as.numeric(strsplit(x, ',')[[1]])
#     bikereg.est(X[-folds[[j]],], Y[-folds[[j]]], time[-folds[[j]]], res[-folds[[j]]], x[1], x[2], H[[i]])
#   }
#   if (is.null(cl)) {
#     stop('Please create a cluster!')
#     for (i in 1:length(H)) {
#       cat(i)
#       beta.table[[i]] <- vector('list', nfolds)
#       for (j in 1:nfolds) {
#         cat(j)
#         beta.table[[i]][[j]] <- t(sapply(coords, foo))
#       }
#     }
#   } else {
#     clusterExport(cl, c('folds', 'X', 'Y', 'time', 'res', 'coords', 'foo', 'H', 'match.coords'), envir=environment())
#     clusterExport(cl, c('bikereg.est', 'MSE', 'bound_gauss'))
#     job <- function(x) {
#       x <- as.numeric(strsplit(x, ',')[[1]])
#       i <- x[1]
#       j <- x[2]
#       beta <- t(sapply(uniq.coords, function(x)foo(x, i, j)))
#       MSE(X[folds[[j]],], Y[folds[[j]]], beta[match.coords[folds[[j]]],], res[folds[[j]]])
#     }
#     jobs.id <- expand.grid(1:length(H), 1:nfolds)
#     result <- parLapply(cl, paste(jobs.id$Var1, jobs.id$Var2, sep=','), job)
#   }
#   lapply(1:length(H), function(x)do.call(rbind, result[jobs.id$Var1==x]))
# }
# MSE <- function(X, Y, beta, res) {
#   temp <- Y-rowSums(X*beta)
#   result <- c(mean(temp^2, na.rm=T), mean(is.na(res)), sum(is.na(temp))/sum(!is.na(res)))
#   names(result) <- c('MSE', 'NA_censored', 'NA_among_uncensored')
#   result
# }
kereg.est <- function(t, X, Y, time, H, kernel=bound_gauss) {
  coords <- do.call(rbind, time)
  k <- kernel(solve(H, coords-t))
  id <- !is.na(k) & !is.na(Y) ## when Y is residual, this could be NA where k is not
  kx <- k*X
  xkx <- t(kx[id,])%*%X[id,]
  result <- try(solve(xkx, t(kx[id,])%*%Y[id])[,1], silent = T)
  if (inherits(result, 'try-error')) {
    rep(NA, dim(X)[2])
  } else {
    result
  }
}
kereg.var <- function(t, X, Y, time, H, kernel=bound_gauss) {
  coords <- do.call(rbind, time)
  k <- kernel(solve(H, coords-t))
  id <- !is.na(k) & !is.na(Y)
  kx <- k*X
  xkx <- t(kx[id,])%*%X[id,]
  xk2x <- t(kx[id,])%*%kx[id,]
  inv <- try(solve(xkx), T)
  if (inherits(inv, 'try-error')) {
    matrix(NA, dim(X)[2], dim(X)[2])
  } else {
    inv%*%xk2x%*%inv
  }
}
kereg.cv <- function(X, Y, time, H, kernel=bound_gauss, nfolds=5, cl=NULL, loss=MSE) {
  id <- !Reduce('|', lapply(time, is.na))
  coords <- do.call(function(...)paste(sep=',', ...), time)
  uniq.coords <- unique(coords[id])
  match.coords <- match(coords, uniq.coords)
  if (!is.list(H) || length(H)==0) stop('H must be a list!')
  folds <- split(which(id), rep(1:nfolds, length.out=sum(id))) ## only CV on complete cases
  ## Input x is the (t,s) where to evaluate beta, i is the id of H, j is the id of fold
  ## Output beta(t,s)
  foo <- function(x, i, j) {
    x <- as.numeric(strsplit(x, ',')[[1]])
    train <- unlist(folds[-j])
    kereg.est(x, X[train,,drop=F], Y[train], lapply(time, function(x)x[train]), H[[i]], kernel=kernel)
  }
  ## Input x = '1,2' means use the first H and the second fold
  ## Output the loss 
  job <- function(x) {
    x <- as.numeric(strsplit(x, ',')[[1]])
    i <- x[1]
    j <- x[2]
    beta <- lapply(uniq.coords, function(x)foo(x, i, j))
    beta <- do.call(rbind, beta)
    loss(X[folds[[j]],,drop=F], Y[folds[[j]]], beta[match.coords[folds[[j]]],,drop=F])
  }
  jobs.id <- expand.grid(1:length(H), 1:nfolds)
  if (is.null(cl)) {
    stop('Please finish progressbar first!')
    result <- lapply(paste(jobs.id$Var1, jobs.id$Var2, sep=','), job)
  } else {
    clusterExport(cl, c('folds', 'X', 'Y', 'time', 'foo', 'H', 'match.coords', 'loss', 'kernel'), envir=environment())
    clusterExport(cl, c('kereg.est'))
    result <- parLapply(cl, paste(jobs.id$Var1, jobs.id$Var2, sep=','), job)
  }
  lapply(1:length(H), function(x)do.call(rbind, result[jobs.id$Var1==x]))
}
MSE <- function(X, Y, beta) {
  res <- Y-rowSums(X*beta)
  result <- c(mean(res^2, na.rm=T), length(res), sum(is.na(res)))
  names(result) <- c('MSE', 'obs_num', 'NA_num')
  result
}
kereg.predict <- function(X, Y, time, H_beta, kernel_beta=bound_gauss, newX=NULL, newY=NULL, newtime=NULL, cl=NULL) {
  if (is.null(newX) & is.null(newY) & is.null(newtime)) {
    newdata <- F
  } else if (!is.null(newX) & !is.null(newY) & !is.null(newtime)) {
    newdata <- T
  } else {
    stop('Please provide all components of the new data!')
  }
  if (newdata) {
    id <- !Reduce('|', lapply(newtime, is.na))
    coords <- do.call(function(...)paste(sep=',', ...), newtime)
    uniq.coords <- unique(coords[id])
    match.coords <- match(coords, uniq.coords)
  } else {
    id <- !Reduce('|', lapply(time, is.na))
    coords <- do.call(function(...)paste(sep=',', ...), time)
    uniq.coords <- unique(coords[id])
    match.coords <- match(coords, uniq.coords)
  }
  foo <- function(x) {
    x <- as.numeric(strsplit(x, ',')[[1]])
    kereg.est(x, X, Y, time, H_beta, kernel=kernel_beta)
  }
  foo2 <- function(x) {
    x <- as.numeric(strsplit(x, ',')[[1]])
    kereg.var(x, X, Y, time, H_beta, kernel=kernel_beta)
  }
  if (is.null(cl)) {
    if (dim(X)[2]==1) {
      beta <- sapply(uniq.coords, foo)
      variance <- sapply(uniq.coords, foo2)
    } else {
      beta <- t(sapply(uniq.coords, foo))
      variance <- aperm(sapply(uniq.coords, foo2, simplify = "array"))
    }
  } else {
    clusterExport(cl, c('X', 'Y', 'time', 'H_beta', 'kernel_beta'), envir=environment())
    clusterExport(cl, c('kereg.est', 'kereg.var'))
    if (dim(X)[2]==1) {
      beta <- parSapply(cl, uniq.coords, foo)
      variance <- parSapply(cl, uniq.coords, foo2)
    } else {
      beta <- t(parSapply(cl, uniq.coords, foo))
      variance <- aperm(parSapply(cl, uniq.coords, foo2, simplify = "array"))
    }
  }
  if (newdata) {
    if (is.null(dim(beta))) {
      epsilon <- newY-rowSums(newX*beta[match.coords])
    } else {
      epsilon <- newY-rowSums(newX*beta[match.coords,])
    }
    list(coef=beta, var=variance, resid=epsilon, time=newtime)
  } else {
    if (is.null(dim(beta))) {
      epsilon <- Y-rowSums(X*beta[match.coords])
    } else {
      epsilon <- Y-rowSums(X*beta[match.coords,])
    }
    list(coef=beta, var=variance, resid=epsilon, time=time)
  }
}
kereg.inference <- function(fit, H_sigma, kernel_sigma=bound_gauss) {
  coords <- sapply(rownames(fit$coef), function(x)strsplit(x, ',')[[1]][1])
  uniq.coords <- unique(coords)
  match.coords <- match(coords, uniq.coords)
  foo <- function(x) {
    x <- as.numeric(x)
    kereg.est(x, matrix(1, length(fit$resid)), fit$resid^2, fit$time[1], H_sigma, kernel=kernel_sigma)
  }
  sigma <- sapply(uniq.coords, foo)
  fit$variance <- fit$variance*sigma[match.coords]
  fit[1:2]
}
plot_ly2 <- function(x, y, z=NULL) {
  uniq.x <- sort(unique(x))
  match.x <- match(x, uniq.x)
  uniq.y <- sort(unique(y))
  match.y <- match(y, uniq.y)
  if (is.null(z)) {
    plot_ly(x=uniq.x, y=uniq.y)
  } else {
    z.mat <- matrix(NA, length(uniq.y), length(uniq.x))
    for (i in 1:length(x)) {
      z.mat[match.y[i], match.x[i]] <- z[i]
    }
    plot_ly(x=uniq.x, y=uniq.y, z=z.mat)
  }
}
# library(dplyr)
# library(tidyr)
# library(plotly)
# plotly.data <- function(x, y, z) {
#   rec <- list()
#   for (i in 1:dim(z)[2]) {
#     temp <- data.frame(x, y, z[,i]) %>% arrange(x, y) %>% pivot_wider(names_from = 'x', values_from = 'z...i.')
#     if (i==1) {
#       yy <- temp$y
#     }
#     rec[[i]] <- select(temp, -y) %>% as.matrix()
#   }
#   rec$x <- as.numeric(colnames(rec[[1]]))
#   rec$y <- yy
#   rec
# }


