require(bpnreg)

#Create vectors with population values
B0I <- c(cos(80*(pi/180)), cos(10*(pi/180)), cos(30*(pi/180)), cos(60*(pi/180)), cos(45*(pi/180)))
B0II <- c(sin(80*(pi/180)), sin(10*(pi/180)), sin(30*(pi/180)), sin(60*(pi/180)), sin(45*(pi/180)))
B1I <- c(0, 1, 1, 0.5, 0, 0.5, 2, 0, 2)
B1II <- c(1, 0, 1, 0, 0.5, 0.5, 0, 2, 2)

mf <- seq(1, 5, by = 0.5)

sim <- 500
N <- 150
p <- 2
mu <- 0
sd <- 1

#Set counts
count <- 0


for(j in 1:length(mf)){
  for(i in 1:5){
    for(k in 1:9){

        count <- count + 1
        
        if(count == 1){
          
          load(file = paste("Simulation/Results180911/", B0I[i]*mf[j], B1I[k]*mf[j], B0II[i]*mf[j], B1II[k]*mf[j], "_", N, "_", mu, "_", sd, ".rda", sep = ""))
          load(file = paste("Simulation/Results180911/","DATA", B0I[i]*mf[j], B1I[k]*mf[j], B0II[i]*mf[j], B1II[k]*mf[j], "_", N, "_", mu, "_", sd, ".rda", sep = ""))

          Dataset <- results
          
        }else{
          
          load(file = paste("Simulation/Results180911/", B0I[i]*mf[j], B1I[k]*mf[j], B0II[i]*mf[j], B1II[k]*mf[j], "_", N, "_", mu, "_", sd, ".rda", sep = ""))
          load(file = paste("Simulation/Results180911/DATA", B0I[i]*mf[j], B1I[k]*mf[j], B0II[i]*mf[j], B1II[k]*mf[j], "_", N, "_", mu, "_", sd, ".rda", sep = ""))

          Dataset <- rbind(Dataset, results)
          
      }
    }
  }
}



Data <- as.data.frame(cbind(Dataset, rep(1:405, each = 500)))

colnames(Data) <- c("Aind1", "Aind2", "Aind3", "Aind4",
                    "Aind5", "Aind6", "Aind7", "Aind8",
                    "Aind1m", "Aind2m", "Aind3m", "Aind4m",
                    "Aind5m", "Aind6m", "Aind7m", "Aind8m",
                    "LBAind1", "UBAind1","LBAind2", "UBAind2",
                    "LBAind3", "UBAind3", "LBAind4", "UBAind4",
                    "LBAind5", "UBAind5", "LBAind6", "UBAind6",
                    "LBAind7", "UBAind7", "LBAind8", "UBAind8",
                    "NEind", "Aind3r", "Aind4r",
                    "Aind5r", "Aind6r", "Aind7r", "Aind8r",
                    "SDO", "spread", "axreal", "design")

Data <- subset(Data, !is.nan(axreal) & SDO != 0) #delete the designs where there is no effect (accuracy/location at all)

Data0 <- subset(Data, SDO == 0)
Data1 <- subset(Data, SDO > 0 & SDO < 1) 
Data2 <- subset(Data, SDO >= 1 & SDO < 1.99)
Data3 <- subset(Data, SDO >= 2 & SDO < 2.99)
Data4 <- subset(Data, SDO >= 3 & SDO < 3.99)
Data5 <- subset(Data, SDO >=3.99)

hist(unique(Data1$SDO))
hist(unique(Data2$SDO))
hist(unique(Data3$SDO))
hist(unique(Data4$SDO))
hist(unique(Data5$SDO))

unique(Data0$axreal)
unique(Data1$axreal)
unique(Data2$axreal)
unique(Data3$axreal)
unique(Data4$axreal)
unique(Data5$axreal)

length(unique(Data1$design))
length(unique(Data2$design))
length(unique(Data3$design))
length(unique(Data4$design))
length(unique(Data5$design))

#SSDO
mean(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind2"])/500}))
mean(sapply(unique(Data2$design), function(w){sum(subset(Data2, design == w)[,"Aind2"])/500}))
mean(sapply(unique(Data3$design), function(w){sum(subset(Data3, design == w)[,"Aind2"])/500}))
mean(sapply(unique(Data4$design), function(w){sum(subset(Data4, design == w)[,"Aind2"])/500}))
mean(sapply(unique(Data5$design), function(w){sum(subset(Data5, design == w)[,"Aind2"])/500}))

mode_est(subset(Data1)[,"UBAind2"] - subset(Data1)[,"LBAind2"])
mode_est(subset(Data2)[,"UBAind2"] - subset(Data2)[,"LBAind2"])
mode_est(subset(Data3)[,"UBAind2"] - subset(Data3)[,"LBAind2"])
mode_est(subset(Data4)[,"UBAind2"] - subset(Data4)[,"LBAind2"])
mode_est(subset(Data5)[,"UBAind2"] - subset(Data5)[,"LBAind2"])


#sin(alpha)
mean(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind3"])/500}))
mean(sapply(unique(Data2$design), function(w){sum(subset(Data2, design == w)[,"Aind3"])/500}))
mean(sapply(unique(Data3$design), function(w){sum(subset(Data3, design == w)[,"Aind3"])/500}))
mean(sapply(unique(Data4$design), function(w){sum(subset(Data4, design == w)[,"Aind3"])/500}))
mean(sapply(unique(Data5$design), function(w){sum(subset(Data5, design == w)[,"Aind3"])/500}))

mode_est(subset(Data1)[,"UBAind3"] - subset(Data1)[,"LBAind3"])
mode_est(subset(Data2)[,"UBAind3"] - subset(Data2)[,"LBAind3"])
mode_est(subset(Data3)[,"UBAind3"] - subset(Data3)[,"LBAind3"])
mode_est(subset(Data4)[,"UBAind3"] - subset(Data4)[,"LBAind3"])
mode_est(subset(Data5)[,"UBAind3"] - subset(Data5)[,"LBAind3"])

#sin(lambda)
mean(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind4"])/500}))
mean(sapply(unique(Data2$design), function(w){sum(subset(Data2, design == w)[,"Aind4"])/500}))
mean(sapply(unique(Data3$design), function(w){sum(subset(Data3, design == w)[,"Aind4"])/500}))
mean(sapply(unique(Data4$design), function(w){sum(subset(Data4, design == w)[,"Aind4"])/500}))
mean(sapply(unique(Data5$design), function(w){sum(subset(Data5, design == w)[,"Aind4"])/500}))

mode_est(subset(Data1)[,"UBAind4"] - subset(Data1)[,"LBAind4"])
mode_est(subset(Data2)[,"UBAind4"] - subset(Data2)[,"LBAind4"])
mode_est(subset(Data3)[,"UBAind4"] - subset(Data3)[,"LBAind4"])
mode_est(subset(Data4)[,"UBAind4"] - subset(Data4)[,"LBAind4"])
mode_est(subset(Data5)[,"UBAind4"] - subset(Data5)[,"LBAind4"])

#sin(gamma)
mean(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind5"])/500}))
mean(sapply(unique(Data2$design), function(w){sum(subset(Data2, design == w)[,"Aind5"])/500}))
mean(sapply(unique(Data3$design), function(w){sum(subset(Data3, design == w)[,"Aind5"])/500}))
mean(sapply(unique(Data4$design), function(w){sum(subset(Data4, design == w)[,"Aind5"])/500}))
mean(sapply(unique(Data5$design), function(w){sum(subset(Data5, design == w)[,"Aind5"])/500}))

mode_est(subset(Data1)[,"UBAind5"] - subset(Data1)[,"LBAind5"])
mode_est(subset(Data2)[,"UBAind5"] - subset(Data2)[,"LBAind5"])
mode_est(subset(Data3)[,"UBAind5"] - subset(Data3)[,"LBAind5"])
mode_est(subset(Data4)[,"UBAind5"] - subset(Data4)[,"LBAind5"])
mode_est(subset(Data5)[,"UBAind5"] - subset(Data5)[,"LBAind5"])

#sin(lambda+gamma)
mean(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind7"])/500}))
mean(sapply(unique(Data2$design), function(w){sum(subset(Data2, design == w)[,"Aind7"])/500}))
mean(sapply(unique(Data3$design), function(w){sum(subset(Data3, design == w)[,"Aind7"])/500}))
mean(sapply(unique(Data4$design), function(w){sum(subset(Data4, design == w)[,"Aind7"])/500}))
mean(sapply(unique(Data5$design), function(w){sum(subset(Data5, design == w)[,"Aind7"])/500}))

mode_est(subset(Data1)[,"UBAind7"] - subset(Data1)[,"LBAind7"])
mode_est(subset(Data2)[,"UBAind7"] - subset(Data2)[,"LBAind7"])
mode_est(subset(Data3)[,"UBAind7"] - subset(Data3)[,"LBAind7"])
mode_est(subset(Data4)[,"UBAind7"] - subset(Data4)[,"LBAind7"])
mode_est(subset(Data5)[,"UBAind7"] - subset(Data5)[,"LBAind7"])


#Daniel
mean(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind8"])/500}))
mean(sapply(unique(Data2$design), function(w){sum(subset(Data2, design == w)[,"Aind8"])/500}))
mean(sapply(unique(Data3$design), function(w){sum(subset(Data3, design == w)[,"Aind8"])/500}))
mean(sapply(unique(Data4$design), function(w){sum(subset(Data4, design == w)[,"Aind8"])/500}))
mean(sapply(unique(Data5$design), function(w){sum(subset(Data5, design == w)[,"Aind8"])/500}))

mode_est(subset(Data1)[,"UBAind8"] - subset(Data1)[,"LBAind8"])
mode_est(subset(Data2)[,"UBAind8"] - subset(Data2)[,"LBAind8"])
mode_est(subset(Data3)[,"UBAind8"] - subset(Data3)[,"LBAind8"])
mode_est(subset(Data4)[,"UBAind8"] - subset(Data4)[,"LBAind8"])
mode_est(subset(Data5)[,"UBAind8"] - subset(Data5)[,"LBAind8"])

#Compare SSDO and Aind in small SDO designs

sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"SDO"])/500})

plot(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"SDO"])/500})
     , sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind8"])/500}))

points(sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"SDO"])/500})
     , sapply(unique(Data1$design), function(w){sum(subset(Data1, design == w)[,"Aind2"])/500}), col = "red")
