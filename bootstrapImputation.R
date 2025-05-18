library(zoo)
data <- read.csv("spiff_data.csv")  #getting time series data

bootstrapmeans <- function(samp, n) {    #Creating n bootstrap samples and computing the mean for each sample
  means <- c()
  for (boot in 1:n) {
    bootsamp <- sample(samp, size = length(samp), replace = TRUE)
    means <- c(means, mean(bootsamp))
  }
  return(means)
}

data <- ts(data)
plot.ts(data)

data[c(36,1194,2836,3430,4118),-c(1,2)] <- NA
for (i in c(36,1194,2836,3430,4118)) {
  data[c(i-1,i,i+1),] <- na.approx(data[c(i-1,i,i+1),])   #replacing outliers by mean
}

plot.ts(data)

for (ts in 3:9) {    #Running imputation procedure for each time series
name <- colnames(data)[ts]
  
tmpdata <- data[,ts]      #Working with a dummy variable to not affect the original data

miss <- which(is.na(tmpdata[1:(length(tmpdata) - 200)]))    #Indeces for missing value gap

#Imputing the missing values by linear interpolation, removing start & end as they are the same as real data
nainterpolated <- approx(c(miss[1]-1, miss[length(miss)] + 1), c(tmpdata[miss[1]-1], tmpdata[miss[length(miss)] + 1]), xout = (miss[1]-1):(miss[length(miss)] + 1))
nostartendinterptmpdata <- nainterpolated$y[-c(1,length(nainterpolated))]

interpresids <- data.frame(NA)
k <- 0
set.seed(1234)
par(mfrow=c(1,1))
plot.ts(tmpdata, xlim = c(0, miss[1]), ylim = c(min(tmpdata[1:(miss[1] - 1)]), max(tmpdata[1:(miss[1] - 1)])), ylab = "interpolated vs real")
interpnames <- c()

#Doing linear interpolations on different parts of data & computing residuals
for (i in c(seq(from = 52, to = miss[length(miss)], by = 51), seq(from = (miss[length(miss)]+1), to = (length(tmpdata) - 200), by = 51) ) ) {
  if (i %in% (miss[1]):(miss[length(miss)] + 51)  ) {
    next
  }
  k <- k + 1
  interp <- approx(c(i-51,i),c(tmpdata[i-51],tmpdata[i]), xout = (i-51):i)
  lines(interp, col = "red")
  interpresids <- cbind(interpresids, (interp$y[-c(1,length(interp$y))] - tmpdata[(i-51 + 1):(i - 1)]))
  interpnames <- c(interpnames, paste0("interpolation ", k))
}

legend(1,max(tmpdata[1:(miss[1] - 1)]), legend = c(paste0("Interpolated  values for ", name)), col = c("red"), lty = c(1))

#Removing first value since we had to include an NA to initialize the empty data frame
interpresids <- interpresids[,-1]
names(interpresids) <- interpnames

bootmeans <- list()
meanofboots <- c()
upperci <- c()
lowerci <- c()
n <- 1000   #Doing 1000 bootstrap samples
for (i in 1:50) {                              #Calculating bootstrap CI and expected value for expected value of residuals
  bootmeans[[i]] <- sort(bootstrapmeans(as.numeric(interpresids[i,]), n))
  meanofboots <- c(meanofboots, mean(bootmeans[[i]]))
  lowerci <- c(lowerci, bootmeans[[i]][0.025*n])
  upperci <- c(upperci, bootmeans[[i]][0.975*n])
  print(paste0("iter", i))
}

#Getting the CI for expected value of the real value
exptrueupper <- nostartendinterptmpdata - lowerci
exptruelower <- nostartendinterptmpdata - upperci

#The new imputation which we now have a confidence interval for
newimputation <- nostartendinterptmpdata- meanofboots
data[miss,ts] <- newimputation #Adding imputed values to data

plot.ts(tmpdata, xlim = c(miss[1] - 150, miss[length(miss)] + 150),              #Plotting the imputed values and CI
        ylim = c(min(tmpdata[(miss[1] - 150):(miss[length(miss)] + 150)], na.rm = TRUE),
                 max(tmpdata[(miss[1] - 150):(miss[length(miss)] + 150)], na.rm = TRUE)),
        ylab = name)
lines(miss,newimputation, col = "red")
lines(miss,exptrueupper, col = "blue", lty = 2)
lines(miss,exptruelower, col = "blue", lty = 2)
legend(miss[1] - 149, max(tmpdata[(miss[1] - 150):(miss[length(miss)] + 150)], na.rm = TRUE), legend = c("New interpolation", "95% confidence interval"), col = c("red","blue"), lty = 1:2)
}

write.csv(data, "bootimputation.txt")   #Writing data to a txt/csv file

