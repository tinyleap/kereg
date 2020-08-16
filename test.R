X <- cbind(1, foo[c('INC_AGE', 'black', 'female', 'tx_cum')])
X <- as.matrix(X)
Y <- foo$final
time <- foo[c('clm_month', 'res_month')]
temp <- kereg.loso(X, Y, time, foo$USRDS_ID, 
                   H=lapply(1:10, function(x)x*diag(2)), cl=1)
