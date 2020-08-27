##########################################################################################
# Author: Odysseas Vlachopoulos
# The code is free software. It canbe redistributed and/or modified it under the terms of  
# the GNU General Public License v3
##########################################################################################
library(maptools)
library(sp)
library(randomForest)
library(raster)
library (rgdal)

directo <- "directory"
setwd(directo)
shapefile <- paste(directo, "shapefile", sep='/')
inputRaster <-'path'

trainShp <- readOGR(dsn = ".", layer = "shapefile") 
inpRas <- brick(inputRaster)

classNums <- c(1,2,3) 
NAvalue(inpRas) <- -10000

features <- data.frame()
response <- numeric()

for (n in 1:3) {
  trainData<- trainShp[trainShp[['class']]==classNums[n],]
  trainAreas <- sapply(slot(trainData, "polygons"), slot, "area")
  samples <- 10000*(trainAreas/sum(trainAreas))
  for (i in 1:nrow(trainData@data)) {
    sampleTrainPoints <- spsample(trainData[i,], type="random", n=samples[i])
    if (i == 1) trainPoints <- sampleTrainPoints
    else trainPoints <- rbind(trainPoints, sampleTrainPoints)
  }
  response <- c(response, over(trainPoints, trainShp)[['class']])
  features <- rbind(features, extract(inpRas, trainPoints))
  
}

training <- na.omit(cbind(response, features))   
randfor <- randomForest(as.factor(response) ~., data=training, importance=TRUE, na.action=na.omit)
## Metrics
print(randfor$confusion)
varImpPlot(randfor)


classRasterName <- paste(directo, "name.tif", sep='/')
classifiedRaster <- raster(inpRas)
classifiedRaster <- writeStart(classifiedRaster, filename=classRasterName, navalue=0, progress='text', format='GTiff', datatype='INT1U', overwrite=TRUE)

bs <- blockSize(inpRas)
for (i in 1:bs$n) {
  imageBlock <-  getValuesBlock(inpRas, row=bs$row[i], nrows=bs$nrows[i])
  predValues <- predict(randfor, imageBlock, type='response')
  classValues <- as.numeric(levels(predValues))[predValues]
  classifiedRaster <- writeValues(classifiedRaster, classValues, bs$row[i])
}
classifiedRaster <- writeStop(classifiedRaster)


