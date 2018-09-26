require(bpnreg)

B0I <- cos(seq(0, 355, by = 5)*(pi/180))
B0II <- sin(seq(0, 355, by = 5)*(pi/180))
B1I <- B0I
B1II <- B0II

sim <- 500
N <- 150
p <- 2
mu <- 0
sd <- 1

#Set counts
count <- 0


for(i in 1:length(B0I)){
    
  count <- count + 1
        
  if(count == 1){
          
    load(file = paste("Simulation/ResultsAccuracy180911/", B0I[i], B1I[i], B0II[i], B1II[i], "_", N, "_", mu, "_", sd, ".rda", sep = ""))
          

    Dataset <- results
          
          
  }else{
          
    load(file = paste("Simulation/ResultsAccuracy180911/", B0I[i], B1I[i], B0II[i], B1II[i], "_", N, "_", mu, "_", sd, ".rda", sep = ""))
      
    Dataset <- rbind(Dataset, results)
          
  }
}
    



Data <- as.data.frame(cbind(Dataset, rep(1:72, each = 500)))

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

length(unique(Data$design))


mean(sapply(unique(Data$design), function(w){sum(subset(Data, design == w)[,"Aind2"])/500}))
mean(sapply(unique(Data$design), function(w){sum(subset(Data, design == w)[,"Aind3"])/500}))
mean(sapply(unique(Data$design), function(w){sum(subset(Data, design == w)[,"Aind4"])/500}))
mean(sapply(unique(Data$design), function(w){sum(subset(Data, design == w)[,"Aind5"])/500}))
mean(sapply(unique(Data$design), function(w){sum(subset(Data, design == w)[,"Aind7"])/500}))
mean(sapply(unique(Data$design), function(w){sum(subset(Data, design == w)[,"Aind8"])/500}))

mode_est(abs(subset(Data)[,"UBAind2"] - subset(Data)[,"LBAind2"]))
mode_est(abs(subset(Data)[,"UBAind3"] - subset(Data)[,"LBAind3"]))
mode_est(abs(subset(Data)[,"UBAind4"] - subset(Data)[,"LBAind4"]))
mode_est(abs(subset(Data)[,"UBAind5"] - subset(Data)[,"LBAind5"]))
mode_est(abs(subset(Data)[,"UBAind7"] - subset(Data)[,"LBAind7"]))
mode_est(abs(subset(Data)[,"UBAind8"] - subset(Data)[,"LBAind8"]))
