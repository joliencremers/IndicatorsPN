#Source Simulation code and required packages for analysis
require(Rcpp)
require(bpnreg)
sourceCpp("Simulation/SamplingData.cpp")

#Create vectors with population values
B0I <- cos(seq(0, 355, by = 5)*(pi/180))
B0II <- sin(seq(0, 355, by = 5)*(pi/180))
B1I <- B0I
B1II <- B0II

sim <- 500
N <- 25
p <- 2
mu <- 0
sd <- 1

#Start Data simulation

#set progress bar
total <- length(B0I)
pb <- winProgressBar(title = "progress bar", min = 0,
                     max = total, width = 300)
count <- 0

#start time
start <- Sys.time()

#start loop


for(i in 1:length(B0I)){
 
  count <- count + 1
  
  Aind1 <- rep(NA, 1, sim)
  Aind2 <- rep(NA, 1, sim)
  Aind3 <- rep(NA, 1, sim)
  Aind4 <- rep(NA, 1, sim)
  Aind5 <- rep(NA, 1, sim)
  Aind6 <- rep(NA, 1, sim)
  Aind7 <- rep(NA, 1, sim)
  Aind8 <- rep(NA, 1, sim)
  
  Aind3r <- rep(NA, 1, sim)
  Aind4r <- rep(NA, 1, sim)
  Aind5r <- rep(NA, 1, sim)
  Aind6r <- rep(NA, 1, sim)
  Aind7r <- rep(NA, 1, sim)
  Aind8r <- rep(NA, 1, sim)
  
  Aind1m <- rep(NA, 1, sim)
  Aind2m <- rep(NA, 1, sim)
  Aind3m <- rep(NA, 1, sim)
  Aind4m <- rep(NA, 1, sim)
  Aind5m <- rep(NA, 1, sim)
  Aind6m <- rep(NA, 1, sim)
  Aind7m <- rep(NA, 1, sim)
  Aind8m <- rep(NA, 1, sim)
  
  Aind1hpd <- matrix(NA, 2, sim)
  Aind2hpd <- matrix(NA, 2, sim)
  Aind3hpd <- matrix(NA, 2, sim)
  Aind4hpd <- matrix(NA, 2, sim)
  Aind5hpd <- matrix(NA, 2, sim)
  Aind6hpd <- matrix(NA, 2, sim)
  Aind7hpd <- matrix(NA, 2, sim)
  Aind8hpd <- matrix(NA, 2, sim)
  
  NEind <- rep(NA, 1, sim)
  measreal <- rep(NA, 1, sim)
  measrealax <- rep(NA, 1, sim)
  SDO <- rep(NA, 1, sim)
  spread <- rep(NA, 1, sim)

  Data <- array(NA, dim = c(N,2,sim))
        
        #Set seed
        set.seed(101)
        
        for(s in 1:sim){
          
          #Sample data
          
          dat <- RData(N, p, c(B0I[i], B1I[i]), c(B0II[i], B1II[i]), mu, sd)
          Data[,,s] <- cbind(dat$Theta, dat$X1[,2])
          
        }
        for(s in 1:sim){
        
          data <- as.data.frame(Data[,,s])
          colnames(data) <- c("theta", "x")
          
          #Analyze Datasets and save results
          
          res <- bpnr(theta ~ x, data, its = 5000, burn = 1000, n.lag = 1, seed = 101)
          
          lin.res <- coef_lin(res)
          circ.res <- coef_circ(res)
          
          axreal = -((B0I[i]*B1I[i] + B0II[i]*B1II[i])/((B1I[i])^2+(B1II[i])^2))
          SDO[s] <- (sqrt(((B0I[i])+ B1I[i]*axreal)^2+((B0II[i])+ B1II[i]*axreal)^2))
          spread[s] <- (sqrt((B0I[i]+ B1I[i]*mean(Data[,2,s]))^2+(B0II[i]+ B1II[i]*mean(Data[,2,s]))^2))
          
          numr   <- sqrt((B0I[i])^2 + (B0II[i])^2)*sqrt((B0I[i] + B1I[i])^2 + (B0II[i] + B1II[i])^2)
          numlxr <- sqrt((B0I[i])^2 + (B0II[i])^2)*sqrt((B0I[i]+ B1I[i]*min(data$x))^2 + (B0II[i]+ B1II[i]*min(data$x))^2)
          numhxr <- sqrt((B0I[i])^2 + (B0II[i])^2)*sqrt((B0I[i]+ B1I[i]*max(data$x))^2 + (B0II[i]+ B1II[i]*max(data$x))^2)
          numlhxr <- sqrt((B0I[i] + B1I[i]*min(data$x))^2 + (B0II[i] + B1II[i]*min(data$x))^2)*sqrt((B0I[i]+ B1I[i]*max(data$x))^2 + (B0II[i]+ B1II[i]*max(data$x))^2)
          numaxr <- sqrt((B0I[i])^2 + (B0II[i])^2)*sqrt((B0I[i]+ B1I[i]*axreal)^2 + (B0II[i]+ B1II[i]*axreal)^2)
          numDr   <- sqrt((B0I[i])^2 + (B0II[i])^2)*sqrt((B1I[i])^2 + (B1II[i])^2)
          
          resmatr      <- c(B0I[i], B0II[i], B0I[i] + B1I[i]            , B0II[i] + B1II[i])
          resmatlowxr  <- c(B0I[i], B0II[i], B0I[i] + B1I[i]*min(data$x), B0II[i] + B1II[i]*min(data$x))
          resmathighxr <- c(B0I[i], B0II[i], B0I[i] + B1I[i]*max(data$x), B0II[i] + B1II[i]*max(data$x))
          resmatlowhighxr <- c(B0I[i] + B1I[i]*min(data$x), B0II[i] + B1II[i]*min(data$x), B0I[i] + B1I[i]*max(data$x), B0II[i] + B1II[i]*max(data$x))
          resmataxr    <- c(B0I[i], B0II[i], B0I[i] + B1I[i]*axreal     , B0II[i] + B1II[i]*axreal)
          resmatDr      <- c(B0I[i], B0II[i], B1I[i], B1II[i])
          
          detr <- resmatr[1]*resmatr[4] - resmatr[2]*resmatr[3]
          detlowxr <- resmatlowxr[1]*resmatlowxr[4] - resmatlowxr[2]*resmatlowxr[3]
          dethighxr <- resmathighxr[1]*resmathighxr[4] - resmathighxr[2]*resmathighxr[3]
          detlowhighxr <- resmatlowhighxr[1]*resmathighxr[4] - resmatlowhighxr[2]*resmathighxr[3] 
          detaxr <- resmataxr[1]*resmataxr[4] - resmataxr[2]*resmataxr[3]
          detDr <- resmatDr[1]*resmatDr[4] - resmatDr[2]*resmatDr[3]
          
          Aind3r[s] <- detr/numr
          Aind4r[s] <- detlowxr/numlxr
          Aind5r[s] <- dethighxr/numhxr
          Aind6r[s] <- detaxr/numaxr
          Aind7r[s] <- detlowhighxr/numlhxr
          Aind8r[s] <- detDr/numDr
          
          num <-      sqrt(res$B1[,1]^2 + res$B2[,1]^2) * sqrt((res$B1[,1] + res$B1[,2])^2 + (res$B2[,1] + res$B2[,2])^2)
          numlowx <-  sqrt(res$B1[,1]^2 + res$B2[,1]^2) * sqrt((res$B1[,1] + res$B1[,2]*min(data$x))^2 + (res$B2[,1]+ res$B2[,2]*min(data$x))^2)
          numhighx <- sqrt(res$B1[,1]^2 + res$B2[,1]^2) * sqrt((res$B1[,1] + res$B1[,2]*max(data$x))^2 + (res$B2[,1]+ res$B2[,2]*max(data$x))^2)
          numlowhighx <- sqrt((res$B1[,1] + res$B1[,2]*min(data$x))^2 + (res$B2[,1]+ res$B2[,2]*min(data$x))^2) * sqrt((res$B1[,1] + res$B1[,2]*max(data$x))^2 + (res$B2[,1]+ res$B2[,2]*max(data$x))^2)
          numax <-    sqrt(res$B1[,1]^2 + res$B2[,1]^2) * sqrt((res$B1[,1] + res$B1[,2]*res$a.x)^2 + (res$B2[,1]+ res$B2[,2]*res$a.x)^2)
          numD <-      sqrt(res$B1[,1]^2 + res$B2[,1]^2) * sqrt(res$B1[,2]^2 + res$B2[,2]^2)
          
          resmat      <- cbind(res$B1[,1], res$B2[,1], (res$B1[,1] + res$B1[,2]), (res$B2[,1] + res$B2[,2]))
          resmatlowx  <- cbind(res$B1[,1], res$B2[,1], (res$B1[,1] + res$B1[,2]*min(data$x)), (res$B2[,1] + res$B2[,2]*min(data$x)))
          resmathighx <- cbind(res$B1[,1], res$B2[,1], (res$B1[,1] + res$B1[,2]*max(data$x)), (res$B2[,1] + res$B2[,2]*max(data$x)))
          resmatlowhighx <- cbind((res$B1[,1] + res$B1[,2]*min(data$x)), (res$B2[,1] + res$B2[,2]*min(data$x)), (res$B1[,1] + res$B1[,2]*max(data$x)), (res$B2[,1] + res$B2[,2]*max(data$x)))
          resmatax    <- cbind(res$B1[,1], res$B2[,1], (res$B1[,1] + res$B1[,2]*res$a.x), (res$B2[,1] + res$B2[,2]*res$a.x))
          resmatD      <- cbind(res$B1[,1], res$B2[,1], res$B1[,2], res$B2[,2])
          
          det <- resmat[,4]*resmat[,1] - resmat[,3]*resmat[,2]
          detlowx <- resmatlowx[,4]*resmatlowx[,1] - resmatlowx[,3]*resmatlowx[,2]
          dethighx <- resmathighx[,4]*resmathighx[,1] - resmathighx[,3]*resmathighx[,2]
          detlowhighx <- resmatlowhighx[,4]*resmatlowhighx[,1] - resmatlowhighx[,3]*resmatlowhighx[,2]
          detax <- resmatax[,4]*resmatax[,1] - resmatax[,3]*resmatax[,2]
          detD <- resmatD[,4]*resmatD[,1] - resmatD[,3]*resmatD[,2]
          
          Aind3m[s] <- mode_est(det/num)
          Aind4m[s] <- mode_est(detlowx/numlowx)
          Aind5m[s] <- mode_est(dethighx/numhighx)
          Aind6m[s] <- mode_est(detax/numax)
          Aind7m[s] <- mode_est(detlowhighx/numlowhighx)
          Aind8m[s] <- mode_est(detD/numD)
          
          Aind3hpd[,s] <- hpd_est(det/num)
          Aind4hpd[,s] <- hpd_est(detlowx/numlowx)
          Aind5hpd[,s] <- hpd_est(dethighx/numhighx)
          Aind6hpd[,s] <- hpd_est(detax/numax)
          Aind7hpd[,s] <- hpd_est(detlowhighx/numlowhighx)
          Aind8hpd[,s] <- hpd_est(detD/numD)
          
          if(lin.res[2,4] < 0 & lin.res[2,5] > 0 & lin.res[4,4] < 0 & lin.res[4,5] > 0){
            
            NEind[s] <- 1
            
          }else{
            
            NEind[s] <- 0
            
          }
          
          if(circ.res["x bc", 4] < 0 &  circ.res["x bc", 5] > 0){
            
            Aind1m[s] <- circ.res["x bc",2]
            Aind1hpd[,s] <- circ.res["x bc", 4:5]
            
            Aind1[s] <- 1
            
          }else{
            
            Aind1m[s] <- circ.res["x bc",2]
            Aind1hpd[,s] <- circ.res["x bc", 4:5]
            
            Aind1[s] <- 0
            
          }
          
          if(circ.res["x SSDO", 4] < 0 &  circ.res["x SSDO", 5] > 0){
            
            Aind2m[s] <- circ.res["x SSDO",2]
            Aind2hpd[,s] <- circ.res["x SSDO", 4:5]
            
            Aind2[s] <- 1
            
          }else{
            
            Aind2m[s] <- circ.res["x SSDO",2]
            Aind2hpd[,s] <- circ.res["x SSDO", 4:5]
            
            Aind2[s] <- 0
            
          }
          
          if(Aind3hpd[1,s] < 0 & Aind3hpd[2,s] > 0){
            
            Aind3[s] <- 1
            
          }else{
            
            Aind3[s] <- 0
            
          }
          
          if(Aind4hpd[1,s] < 0 & Aind4hpd[2,s] > 0){
            
            Aind4[s] <- 1
            
          }else{
            
            Aind4[s] <- 0
            
          }
          
          if(Aind5hpd[1,s] < 0 & Aind5hpd[2,s] > 0){
            
            Aind5[s] <- 1
            
          }else{
            
            Aind5[s] <- 0
            
          }
          
          if(Aind6hpd[1,s] < 0 & Aind6hpd[2,s] > 0){
            
            Aind6[s] <- 1
            
          }else{
            
            Aind6[s] <- 0
            
          }
          
          if(Aind7hpd[1,s] < 0 & Aind7hpd[2,s] > 0){
            
            Aind7[s] <- 1
            
          }else{
            
            Aind7[s] <- 0
            
          }
          
          if(Aind8hpd[1,s] < 0 & Aind8hpd[2,s] > 0){
            
            Aind8[s] <- 1
            
          }else{
            
            Aind8[s] <- 0
            
          }
          
          
        }
        
        results <- rbind(Aind1, Aind2, Aind3, Aind4, Aind5, Aind6, Aind7, Aind8,
                         Aind1m, Aind2m, Aind3m, Aind4m, Aind5m, Aind6m, Aind7m, Aind8m,
                         Aind1hpd, Aind2hpd, Aind3hpd, Aind4hpd, Aind5hpd, Aind6hpd, Aind7hpd, Aind8hpd,
                         NEind, Aind3r, Aind4r, Aind5r, Aind6r, Aind7r, Aind8r,
                         SDO, spread, axreal)
        
        save(results, file = paste("Simulation/ResultsAccuracy180911/", B0I[i], B1I[i], B0II[i], B1II[i], "_", N, "_", mu, "_", sd, ".rda", sep = ""))
        
        setWinProgressBar(pb, count, title=paste( round(count/total*100, 0),
                                                  "% done"))
        
      

}
close(pb)

#end time
end <- Sys.time()


