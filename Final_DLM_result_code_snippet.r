GG <- track_model$GG
WW <- track_model$W
VV <- track_model$V

a0 <- xreg_filtered$m %>% tail(1) %>% as.numeric()
R0 <- xreg_filtered %$% dlmSvd2var(U.C %>% tail(1), D.C %>% tail(1)) %>% extract2(1)


np <- 10
a <- matrix(data = NA, nrow = np, ncol=a0 %>% length())
R <- list()
FF <- cbind(1, dfe_test %>% pull(which_track)) %>%  plyr::alply(1,identity)

Q <- list()
f <- matrix(0, nrow=np, ncol=1)


for (i in 1:np){
  if(i==1){
    a[i,] <- GG %*% a0
  } else {
    a[i,] <- a[i-1,]
  }
  f[i,] <- FF[[i]] %*% a[i,]
  if (i==1){
    R[[1]] <- GG %*% R0 %*% t(GG) + WW
  } else {
    R[[i]] <- GG %*% R[[i-1]] %*% t(GG) + WW
  }
  Q[[i]] <- FF[[i]] %*% R[[i]] %>% tcrossprod(FF[[i]]) +VV
}




bpm_43 <- pull(bpm_all, 43) %>% ts(start=1800, deltat=30)
upperCI3 <- f + qnorm(0.975)*sqrt(unlist(Q)) %>% ts(start = 7200, deltat=30)
lowerCI3 <- f - qnorm(0.975)*sqrt(unlist(Q)) %>% ts(start = 7200, deltat=30)

f_strength <- f  %>% ts(start = 7200, deltat=30)

plot(window(bpm_43, start=6500), type = "o", ylab="BPM", ylim=range(c(upperCI3,lowerCI3, bpm_43)), main="ID 62: 10-steps Forecast: Pooled Dynamic Regression")


lines(f_strength, type="o", pch="x", col="red")
lines(upperCI3, lty=c(5), col="red")
lines(lowerCI3, lty=c(5), col="red")

abline(v=time(f_strength)[1],
       lty="dashed")

guide3 <- c("data", "N-step-ahead forecast pooled dynamic regression ", "95% forecast interval")

legend("bottomleft", col=c("black", "blue", "red" ),
       lty = c(1,1,5), pch=c("o","x",NA),
       legend = guide3, cex=0.9)