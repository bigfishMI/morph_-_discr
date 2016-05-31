############################################################################
#
#       Some sample code for the BIGFISH talk on Morphometrics
# 
#       CN 1/06/2016
#
#       Momocs package depends R > 3.2
############################################################################

# remove previous datasets
rm(list = ls())
dev.off()


#install.packages("Momocs")
#install.packages("jpeg")
library(Momocs)
library(jpeg)

getwd()
setwd(".\\Otolith photos BW")

# Get the list of otolith jpps
files <- list.files()
(files <- files[grep(".jpg", files, ignore.case = T)])

#plot a couple of images
img <- readJPEG(files[1])
plot(1:2, type='n')
rasterImage(img, 1.2, 1.27, 1.8, 1.73) #grayscale

img <- readJPEG(files[3])
plot(1:2, type='n')
rasterImage(img, 1.2, 1.27, 1.8, 1.73) #black and white

# Just plotting an image is not useful
# need the outline => use import_jpg() form the Momocs package
jpgs <- import_jpg(files)
class(jpgs)
length(jpgs)
head(jpgs[[1]])
plot(jpgs[[2]][,1], jpgs[[2]][,2])
plot(jpgs[[3]][,1], jpgs[[3]][,2]) 


# The import is just a list of coordinates
# So we convert it to a Momocs object called a "Coo" (Out, Opn)
oto <- Out(jpgs)
oto
names(oto)
panel(oto)
stack(oto)
#coo_plot(oto)

# center all the otoliths on the origin
oto <- coo_center(oto) 
stack(oto)

# scale them
oto <- coo_scale(oto) 
stack(oto)

# aligning, note this rotates AND FLIPS the otolith 
oto <- coo_align(oto) 
stack(oto)

# rotate by 180 degrees
oto <- coo_rotate(oto, theta = pi) 
stack(oto) 

#flip around x axis
oto <- coo_flipx(oto) 
stack(oto) 


#Perimieter
coo_perim(oto[1])
#No. of points
coo_nb(oto)

#remove 3 and 8
oto[3] <- NULL
coo_nb(oto)
oto[7] <- NULL


stack(oto)




#############################################################################
#
# Bigger, cleaner data set with otolith outlines processed externally
# Otoliths outlines traced and quality controlled in TPSDIG
#
# These are samples from spawning herring from two locations, one in Scotland one in Ireland
# The naming convention is Y-SXX-ZZZ.TPS
# Where Y = year (4 = 2004)
# s = Spawning sample
# XX = location (04 = Donegal, 10 = Cape Wrath)
# ZZZ = fish ID
#
############################################################################


setwd("..\\spawning Otolith TPS")

(files <- list.files())
(files <- files[grep(".tps", files, ignore.case = T)])
length(files)


y <- list()
for(j in files){
  lines <- as.numeric(gsub("POINTS=", "", readLines(j)[3]))
  tps <- read.table(j, skip = 3, nrows = lines, colClasses = "numeric", col.names = c("x", "Y"))
  y[[j]]<- data.matrix(tps)
}
y[1]

##### Converting to Coo and normalising (still need to add factors for hauls, years etc)

oto <- Out(y)
oto
panel(oto)
stack(oto, borders="#1A1A1A22", ldk=FALSE)


oto <- coo_center(oto) # center all the otoliths on the origin
oto <- coo_scale(oto) # scale them
oto <- coo_align(oto) # aligning 
oto <- coo_rotate(oto, theta = pi) # rotate by 180 degrees
stack(oto, borders="#1A1A1A22", ldk=FALSE) 


# Everything above is just a list of matrices of x,y points. Very hard to analyse
# So we convert to a useable format with Elliptic Fourier Analysis (function = efourier)


# Plot showing how successive EFA harmonics improve the shape
coo_plot(oto[1], lwd=2, col='#F2F2F2')
for(i in 1:15){
    ef <- efourier(oto[1], i) # get the ith harmonic of elliptic fourier
    efi <- efourier_i(ef) # get the inverse in order to plet it in the same frame
    coo_draw(efi)
    cat ("Press [enter] to continue")
    line <- readline()
}


#plotting the original shape and the 30th EFA harmonic
coo_plot(oto[1], lwd=1, col='#F2F2F2')
coo_draw(efourier_i(efourier(oto[1], nb.h = 30)), lwd = 2)


# Calculating the Elliptic Fourier Analysis (efourier) for the whole data set
efa<-efourier(oto, nb.h = 30, smooth.it = 0, norm = TRUE, start = FALSE) #start = False "does not define the first point of the outlines as an homologous point" 
class(efa)
efa

# Performing a principal component analysis on the EFA 
efa.pca <- PCA(efa)
efa.pca
plot(efa.pca, morpho=TRUE)


# Adding Factors
oto
names(oto)
(fact <- sapply(strsplit(names(oto), "-"), function(x)  x[2]))
class(fact)
                                    #oto$fac <- data.frame(Type = c(rep("A", length(oto)/2), rep("B", length(oto)/2)))

#set the spawning location (i.e. the stock) as the factor
oto$fac <- data.frame(Stock = fact)
oto

#efourier
efa <- efourier(oto, nb.h = 30, smooth.it = 0, norm = TRUE, start = FALSE)

#PCA on the EFA
efa.pca <- PCA(efa)
efa.pca
plot(efa.pca, morpho=TRUE)
plot(efa.pca, 'Stock')
PCcontrib(efa.pca)


# What is the mean shape of otolith in each stock?

mshapes(efa) # the mean (global) shape
ms <- mshapes(efa, 'Stock')
ms$Coe
class(ms$Coe)
ms <- ms$shp
coo_plot(ms$S04A)
coo_draw(ms$S10A, border='forestgreen')
tps_grid(ms$S04A, ms$S10A) # deformation grid
tps_iso(ms$S04A, ms$S10A) # deformation iso


# A better example

data(bot)
botF <- efourier(bot, nb.h = 15)
x <- mshapes(botF, 'type', nb.pts=80)$shp
fr <- x$beer
to <- x$whisky
tps_arr(fr, to, arr.nb=200, palette=col_sari, amp=3)
tps_arr(fr, to, arr.nb=200, palette=col_sari, amp=3, grid=FALSE)
tps_grid(fr, to, amp=3, grid.size=10)



# Linear Discriminant Analysis (LDA)
oto.f <- efourier(oto, 30)
oto.p <- PCA(oto.f)
LDA(oto.p, 'Stock', retain=0.99) # retains 0.99 of the total variance
LDA(oto.p, 'Stock', retain=5) # retain 5 axis
oto.l <- LDA(oto.p, 'Stock', retain=0.99)
oto.l
plot(oto.l)







###########################  shapeR###########################




#install.packages("shapeR")
library(shapeR)


#fourier
data(shape)
shape = stdCoefs(shape,classes="pop","length_cm")
plotFourier(shape,class.name= "pop",useStdcoef=TRUE)
plotFourierShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1)

#wavelets
shape = stdCoefs(shape,classes="pop","length_cm")
plotWavelet(shape,level=5,class.name= "pop",useStdcoef=TRUE)
# plot mean wavelet shape
plotWaveletShape(shape, "pop",show.angle = TRUE,lwd=2,lty=1)

# Calculate the mean otolith area for each fish population
# The results are in square mm since the calibration ('cal') column
# in the data file is in pixels (1 mm/pixel).
tapply(getMeasurements(shape)$otolith.area, getMasterlist(shape)$pop,mean)






