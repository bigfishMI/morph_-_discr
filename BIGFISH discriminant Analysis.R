############################################################################
#
#       Some sample code for the BIGFISH talk on Discriminant Analysis
# 
#       CN 1/06/2016
#
#       
############################################################################

# remove previous datasets
rm(list = ls())
dev.off()

#install.packages("MASS")
#install.packages("klaR")
#install.packages("mvoutlier")
#install.packages("nortest")
#install.packages("rrcov")
library(MASS)
library(nortest)
library(klaR)
library(rrcov)
library(mvoutlier)


# Simple example of Linear Discriminant Analysis
# Discriminant analysis uses a data set of known origin to classify unknown samples

#Make a random data set
Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
                   Sp = rep(c("s","c","v"), rep(50,3)))
Iris
plot(Iris$Petal.L., Iris$Petal.W., col = as.numeric(Iris$Sp))

#Assign a random half of Iris as the training data set
(train <- sample(1:150, 75))
table(Iris$Sp[train])

plot(Iris$Petal.L.[-train], Iris$Petal.W.[-train], col = "blue")
points(Iris$Petal.L.[train], Iris$Petal.W.[train], col = as.numeric(Iris$Sp[train]))

# Calculate the LDA and predict the missing species
(z <- lda(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train))
plot(z)
(x <- predict(z, Iris[-train, ]))
x$class

# make a scatterplot
plot(x$x[,1],x$x[,2], 
     col = as.numeric(Iris$Sp[-train]),
     main="LDA of Iris", 
     pch = as.numeric(Iris$Sp[-train])) 
z

# More plots of the results
# plot(Iris$Petal.L.[-train], Iris$Petal.W.[-train], col = as.numeric(Iris$Sp[-train]))
# points(Iris$Petal.L.[-train], Iris$Petal.W.[-train], col = as.numeric(x$class), pch = 4)
# 
# plot(Iris$Sepal.L.[-train], Iris$Sepal.W.[-train], col = as.numeric(Iris$Sp[-train]))
# points(Iris$Sepal.L.[-train], Iris$Sepal.W.[-train], col = as.numeric(x$class), pch = 4)




# Same example, better plots showing all permutaions and lines of classification 
#install.packages("klaR")

data(iris)
iris
partimat(Species ~ ., data = iris, method = "lda")

partimat(Species ~ ., data = iris, method = "lda", plot.matrix = TRUE, imageplot = FALSE) 

#Different types of analyses
# 1: Quadratic Discriminant Analyses (sensitive to non-normal data)
partimat(Species ~ ., data = iris, method = "qda")

# 2: regularized (supposed to be more robust against multicollinearity in the data)
partimat(Species ~ ., data = iris, method = "rda")




########################### ##############################################################
#
#      These examples are fine for small datasets but what about a large one with >150 variables?
#       



# Read in a bigger dataset. 
# Body and otolith morphometrics from samples of spawning herring in 8 different 
# locations around Ireland and the UK, plus an outgroup from the Baltic

data     <- read.csv(file.path(".", "DA datasets", "Baseline_data_FINAL.csv"),sep=",", header=TRUE, stringsAsFactors=F)
dim(data)
head(data)

# Create indexes for the 167 variables
colnames(data)
stock      <- match("Stock", colnames(data))
body.morph <- 2:38
oto.morph  <- 44:160
oto.meas   <- 39:43
misc.meas  <- 161:164
morph      <- c(body.morph,oto.meas,oto.morph)
all.meas   <- c(body.morph,oto.meas,oto.morph,misc.meas)

head(data[,oto.meas])


# Some quick checks (sensitive to outliers)
boxplot(data[,body.morph])
boxplot(data[,oto.morph])
boxplot(data[,oto.meas])



# make training portion
training.data <- data[1:dim(data)[1] %% 2 == 0, ]
dim(training.data)

#make the 'unknown' portion
test.data <- data[1:dim(data)[1] %% 2 != 0, ]
dim(test.data)

######################################
#
#         Stepwise variable Variable Selection
#         i.e. how to reduce from 168 variables and remove redundants
#
# stepwise forward-backward variable selection with a random starting variable
# improvement = 0.01   stop criteria of 1% improvement.  Results in 8 variables. roughly 3 mins to run
# improvement = 0.005  stop criteria of 0.5% improvement.  16 variables. Roughly 22 mins to run
# improvement = 0.001  stop criteria of 0.1% improvement.  44 variables. Roughly 2 hrs to run

stepwise <- stepclass(Stock ~ ., 
                      data = training.data[,c(morph, stock)], 
                      method = "qda", 
                      direction = "both",  
                      criterion = "AS", 
                      improvement = 0.01, 
                      start.vars = "TH04")

#load("stepwise variable selection.RData")

stepwise
names(stepwise)
stepwise$process
stepwise$result.pm
stepwise$formula
stepwise$model

step.vars <-  match(stepwise$model[,"name"], colnames(data))

#check the indexing
colnames(data)[step.vars]


source("normal functions.R")

# Normality tests on the chosen variables
par(mfrow = c(2,2))
norm.tests(data, step.vars)


#multivariate normality by stock
par(mfrow = c(1,1))
graph.mvn(data, step.vars)
graph.mvn.Stock(data, step.vars)


#    check for and remove  outliers
res <-chisq.plot((data[,step.vars]))
## remove outliers
dim(data)
data <- data[-res$outliers,]
dim(data)
par(mfrow = c(1,1))
graph.mvn(data, step.vars)
graph.mvn.Stock(data, step.vars)

# make training portion again
training.data <- data[1:dim(data)[1] %% 2 == 0, ]
dim(training.data)

#make the 'unknown' portion again
test.data <- data[1:dim(data)[1] %% 2 != 0, ]
dim(test.data)


###########################################################################
#                                       
#                    Linear Discriminat Analysis  
#
###########################################################################

par(mfrow = c(1,1))


colnames(training.data)[c(stock, step.vars)]
training.data$Stock
training.data$Stock <- as.factor(training.data$Stock)

# Apply the linear DA
training.fit <- lda(Stock ~ ., data=training.data[c(stock, step.vars)]) # fit the linear discriminant analysis to the training data, including the stock
training.fit

# Linear Discriminant Analysis with Jacknifed Prediction 
#training.fit <- lda(Stock ~ ., data=training.data[,c(stock,step.vars)], CV=TRUE)

# predict the stock of each fish using the fit
training.predicted <- predict(training.fit, training.data[,step.vars])  

plot(training.predicted$x[,1],training.predicted$x[,2], col = as.numeric(training.data$Stock),main="training data", pch = as.numeric(training.data$Stock)) # make a scatterplot
legend("topright", col = as.numeric(unique(training.predicted$class)), pch = as.numeric(unique(training.predicted$class)), legend = substr(as.character(unique(training.predicted$class)),1,4))


# Assess the accuracy of the prediction
# percent correct for each stock 
(ct <- table(training.data$Stock, training.predicted$class))
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))

#histogram of posterior probabilities
base.hist <- hist(apply(training.predicted$posterior, MARGIN=1, FUN=max),breaks=6, main = "LDA Posterior Probabilities", xlab = "Posterioir probability")


# Now that we know the fitted LDA is good we can use it to predict unknown classes
test.predicted <- predict(training.fit, test.data[,step.vars])  
# make a scatterplot
plot(test.predicted$x[,1],test.predicted$x[,2], 
     col = as.numeric(test.predicted$class),
     main="Predicted Classification", 
     pch = as.numeric(test.predicted$class)) 

legend("topright", 
       col = as.numeric(unique(test.predicted$class)), 
       pch = as.numeric(unique(test.predicted$class)), 
       legend = substr(as.character(unique(test.predicted$class)),1,4))

hist(apply(test.predicted$posterior, MARGIN=1, FUN=max),breaks=6, main = "LDA Posterior Probabilities", xlab = "Posterior probability")


###########################################################################
#                                       
#                    Quadratic Discriminat Analysis  
#
###########################################################################

# the difference between QDA and LDA is that the former permits
# each group distribution to have its own covariance matrix, whilst the latter assumes
# a common covariance matrix for all group distributions

#QDA is more sensitive to non-normal data however

# fit (train) the QDA
training.fit <- qda(Stock ~ ., data=training.data[,c(stock, step.vars)]) 
training.fit
names(training.fit)
class(training.fit)



# with jack knife prediction (leave one out cross-validation)
# training.fit <- qda(Stock ~ ., data=training.data[,c(stock, step.vars)], CV = TRUE)


# Predict the classes (stocks) in the training dataset using the trained QDA
# i.e. assess its accuracy
training.predicted <- predict(training.fit, training.data[,step.vars])  # predict the stock of each fish using the fit
names(training.predicted)
training.predicted$class



# Assess the accuracy of the prediction
# percent correct for each stock 
(ct <- table(training.data$Stock, training.predicted$class))
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))


for(i in unique(training.predicted$class)) {
  hist(training.predicted$posterior[training.predicted$class == i,i], 
       main = i, xlab = "Posterior Probability")
}



test.predicted <- predict(training.fit, test.data[,step.vars])
class(test.predicted)
names(test.predicted)
test.predicted$class

(ct <- table(test.data$Stock, test.predicted$class))
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))





#Plotting QDAs (using rrcov package this time)
training.fit2 <- QdaClassic(x=training.data[,step.vars], grouping=training.data$Stock)
training.predicted2 <- predict(training.fit2, training.data[,step.vars])

head(training.predicted2@x)
plotdat<-training.predicted2@x
plot(plotdat[,1], plotdat[,2], col = 1:8)
plot(plotdat[,2], plotdat[,3], col = 1:8)

plot(plotdat[,1], plotdat[,2], col = as.numeric(training.predicted2@classification))
plot(plotdat[,2], plotdat[,3], col = as.numeric(training.predicted2@classification))

# 3d rotatable plot

#install.packages("rgl")
library(rgl)
plot3d(x = plotdat[,1], 
       y = plotdat[,2], 
       z = plotdat[,3], 
       xlab = "QDA Axis 1", 
       ylab = "QDA Axis 2", 
       zlab = "QDA Axis 3",
       main="training",  
       col = as.numeric(training.predicted2@classification),
       size = 5) 


# Still too cluttered => group the stocks
levels(training.data$Stock)
levels(training.data$Stock) <- list("VIaS" = c("S03 - ROSAMHIL", "S04 - DONEGAL"),
                                    "VIaN" = c("S10A - CAPE WRATH SPRING", "S10B - CAPE WRATH SUMMER"),
                                    "Other" = c("S01 - CELTIC SEA","S02 - SW IRELAND","S06 - IRISH SEA"),
                                    "Outgroup" = "X01A - BALTIC SEA")
training.data$Stock

training.fit3 <- QdaClassic(x=training.data[,step.vars], grouping=training.data$Stock)
training.predicted2 <- predict(training.fit3, training.data[,step.vars])
plotdat<-training.predicted2@x
plot3d(x = plotdat[,1], 
       y = plotdat[,2], 
       z = plotdat[,3], 
       xlab = "QDA Axis 1", 
       ylab = "QDA Axis 2", 
       zlab = "QDA Axis 3",
       main="training",  
       col = as.numeric(training.predicted2@classification),
       size = 5) 

levels(training.data$Stock)
#black = VIaS
# red = VIaN
# green = Other
# blue = Out


