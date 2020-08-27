##########################################################################################
# Author: Odysseas Vlachopoulos
# The code is free software. It canbe redistributed and/or modified it under the terms of  
# the GNU General Public License v3
##########################################################################################
library(raster)
library(data.table)
library(maptools)
library(sp)
library(randomForest)
library(raster)
library(rgdal)

gc()

direct <- "directory"
setwd(direct)
inputRaster <- paste(direct, "name.tif", sep="/")
segmentationRaster <- paste(direct, "name.tif", sep="/")

nodata <- -10000

percMean <- function(nums,...) {
  percentile=quantile(nums, probs=seq(0,1,0.1), na.rm = TRUE)
  mean(nums[which(nums >= percentile[2] & nums <= percentile[10])])
}

percSD <- function(nums,...) {
  percentile=quantile(nums, probs=seq(0,1,0.1), na.rm = TRUE)
  sd(nums[which(nums >= percentile[2] & nums <= percentile[10])])
}

segmRas <- raster(segmentationRaster)
NAvalue(segmRas) <- nodata

inpRas <- brick(inputRaster)
NAvalue(inpRas) <- nodata

vals <- getValues(inpRas)

rasterWithSegm <- data.table(vals, seg = getValues(segmRas))
setkey(rasterWithSegm, seg)

median <- rasterWithSegm[, lapply(.SD, "median"), by=seg]
mean <- rasterWithSegm[, lapply(.SD, "mean", na.rm = TRUE), by=seg]
SD <- rasterWithSegm[, lapply(.SD, "sd", na.rm = TRUE), by=seg]
percMean <- rasterWithSegm[,lapply(.SD, percMean), by=seg]
percSD <- rasterWithSegm[,lapply(.SD, percSD), by=seg]

bands <- sort(c("RedNorm", "GreenNorm", "BlueNorm", "NIRNorm", "RedEdgeNorm", "BRSRNorm"))
means = SDs = percMeans = percSDs = medians = bands 
for (j in 1:length(bands)){means[j]=paste("Mean", bands[j], sep="_")}
for (j in 1:length(bands)){SDs[j]=paste("SD", bands[j], sep="_")}
for (j in 1:length(bands)){percMeans[j]=paste("PercMean", bands[j], sep="_")}
for (j in 1:length(bands)){percSDs[j]=paste("PercSD", bands[j], sep="_")}
for (j in 1:length(bands)){medians[j]=paste("Median", bands[j], sep="_")}

features <- data.frame(mean, SD, percMean, percSD, median)
drop <- c("seg.1", "seg.2", "seg.3", "seg.4")
features <- features[ , !(names(features) %in% drop)]

names(features) <- c("segnumber", means, SDs, percMeans, percSDs, medians)

gc()

shapefile <- paste(direct, "filename", sep='/')
trainShp <- readOGR(dsn = ".", layer = "filename") 

traindt <- trainShp@data

extractSegments <- extract(segmRas, trainShp, cellnumbers=TRUE)
extractSegments <- extractSegments[!sapply(extractSegments, is.null)]

trainingSegments <- matrix(ncol=3, nrow=0)
for (i in 1:length(extractSegments)) {
  lineResponse <- traindt[i,"class"]
  if (is.matrix(extractSegments[[i]])) {segments <- extractSegments[[i]][which(duplicated(extractSegments[[i]][,2]) == FALSE),]}
  else {segments <- extractSegments[[i]]}
  
  if (is.matrix(segments)) {trainingSegments <- rbind(trainingSegments, cbind(lineResponse, segments))}
  else {trainingSegments <- rbind(trainingSegments, cbind(lineResponse, segments[1], segments[2]))}  
}

trainingSegments_clear <- as.data.frame(na.omit(trainingSegments))
names(trainingSegments_clear) <- c("response", "cell", "segment")

features <- na.omit(features)
training <- as.data.frame(na.omit(features[match(trainingSegments_clear$segment, features$segnumber),]))
responseSegm <- trainingSegments_clear[match(training$segnumber, trainingSegments_clear$segment), c(1,3)]

randfor <- randomForest(as.factor(responseSegm$response) ~. , data=training[,-1], importance=TRUE)
## Metrics
print(randfor$confusion)
varImpPlot(randfor)

bs <- blockSize(segmRas)
classRasterName <- paste(direct, "name.tif", sep='/')
img.out <- brick(segmRas)
img.out <- writeStart(img.out, classRasterName, overwrite=TRUE, datatype='INT1U')

predictions <- predict(randfor, features, type='response')
predictionsSegm <- data.frame(features$segnumber, predictions)

for (i in 1:bs$n) {
  img <- getValues(segmRas, row=bs$row[i], nrows=bs$nrows[i])
  outMatrix <- matrix(nrow=length(img), ncol=0)
  is.na(img) <- img == nodata
  outMatrix <- cbind(outMatrix, predictionsSegm$predictions[match(img, predictionsSegm$features.segnumber)])
  writeValues(img.out, outMatrix, bs$row[i])
}

img.out <- writeStop(img.out)

gc()