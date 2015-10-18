
######One-way ANOVA ###### 
######
######
###### load tidyr library and import protein data

library(tidyr)

proteinData <- read.csv("protein-expression.csv")
head(proteinData)
dim(proteinData)
summary(proteinData)

###### Show the distributions

library(beeswarm)
library(RColorBrewer)
boxplot(proteinData,xlab="Cell Type",ylab="Protein Expression",main="Protein Expression")
beeswarm(proteinData, add=TRUE,pch=16,col=brewer.pal(5,"Set1"),method="swarm")

###### Use a histogram

par(mfrow=c(2,3))
cols <- brewer.pal(5,"Set1")
hist(proteinData$A,xlab="A",col=cols[1],main="")
hist(proteinData$B,xlab="B",col=cols[2],main="")
hist(proteinData$C,xlab="C",col=cols[3],main="")
hist(proteinData$D,xlab="D",col=cols[4],main="")
hist(proteinData$E,xlab="E",col=cols[5],main="")

###### Log-transform
proteinData.ln <- log(proteinData)
head(proteinData.ln)

###### check histograms again
par(mfrow=c(2,3))
hist(proteinData.ln$A,xlab="A",col=cols[1],main="")
hist(proteinData.ln$B,xlab="B",col=cols[2],main="")
hist(proteinData.ln$C,xlab="C",col=cols[3],main="")
hist(proteinData.ln$D,xlab="D",col=cols[4],main="")
hist(proteinData.ln$E,xlab="E",col=cols[5],main="")


###### fit the model and check

anovaData <- gather(proteinData.ln)
head(anovaData)
mod <- aov(value~ key, data=anovaData)
mod
par(mfrow=c(2,2))
plot(mod)


summary(aov(mod))

###### barlett test for variances

bt <- bartlett.test(value~key,data=anovaData)
bt

###### post-hoc tests
post.tests<- TukeyHSD(mod)
post.tests

######  Kruskal-Wallis test ######
###### Input headache data by-hand
######

headache <- matrix(c(62,74,86,74,91,37, 
                     69,43,100,94,100,98, 
                     50,-120,100,-288,4,-76), ncol=3)

colnames(headache) <- c("Relaxation / response feedback","Relaxation alone","Untreated")
par(mfrow=c(1,1))
boxplot(headache)

###### tidy the data
library(tidyr)
headache <- data.frame(headache)
headache <- gather(headache)
 
###### do the kruskal-wallis test and check the output
kt <- kruskal.test(value~key,data=headache)
kt
names(kt)
kt$statistic
kt$p.value
######

######  Friedman’s test ######
###### input the example data by-hand and visualise
######
rubinPeters <- data.frame(Person = 1:4, 
                          Before=c(22.2,17,14.1,17), 
                          hrs48 = c(5.4,6.3,8.5,10.7), 
                          mnths6 = c(10.6,6.2,9.3,12.3))
boxplot(rubinPeters[,2:4])
beeswarm(rubinPeters[,2:4],add=TRUE)

###### do the friedman's test
fr <-  friedman.test(as.matrix(rubinPeters[,-1]))
fr
######

###### Median Test ######
######
######

###### RVAideMemoire required for median test
library(RVAideMemoire)

###### input data by-hand

mmse <- data.frame(Group1=c(19,7,17,28,21,6,21,19,27,8,25), 
                   Group2 = c(16,22,30,24,22,23,22,28,29,29,0),
                   Group3 = c(4,9,30,29,25,22,25,26,27,18,10)
)

###### tidy the data
mmse <- gather(mmse)

###### do the median test

mood.medtest(value~key,data=mmse,exact=FALSE)

mood.medtest(value~key,data=mmse)

######  Jonckheere-Terpstra Test ######
######
######

###### clinfun is required for J-T test

library(clinfun)
grade <- data.frame(Grade1 = c(1.99,3.01,4.17,7.13,9.82,9.91,NA,NA), 
                    Grade2 = c(4.40,9.82,10.23,11.99,11.99,13.17,13.20,NA),
                    Grade3 =c(6.94,8.04,9.82,15.75,18.30,25.01,26.40,28.17))

###### tidy the data

grade <- gather(grade)
grade$key <- as.numeric(gsub("Grade","",grade$key))

###### do the J-T test
jonckheere.test(grade$value,grade$key)

###### Linear Regression ######
######
######

###### example data is in the csv file
data <- read.csv("lactoferrin.csv")
data

###### scatter plot 
plot(data$conc, data$growth, pch=16, xlab="Concentration", ylab="Growth rate")




###### fit a linear model
model <- lm(growth ~ conc, data=data)
model

###### overlay the model on the scatter plot
plot(data, pch=16, xlab="Concentration", ylab="Growth rate")
abline(model)

######
fitted <- fitted(model)
fitted

######
residuals <- data$growth - fitted
residuals

######
coefficients(model)
residuals(model)
names(model)
model$coefficients
model$residuals

######
mean(data$conc)
mean(data$growth)
model$coefficients[1] + model$coefficients[2] * mean(data$conc)


###### Calculating standard errors in the regression parameters
SSE <- sum(residuals^2)

######
mean_growth <- mean(data$growth)
mean_growth
SSY <- sum((data$growth - mean_growth)^2)
SSY
SSR <- sum((fitted - mean_growth)^2)
SSR
SSE <- SSY - SSR
SSE

######
error_variance <- SSE / (length(data$growth) - 2)
error_variance
F_ratio <- SSR / error_variance
F_ratio

###### Summarizing the linear model

summary(model)

######
summary <- summary(model)
names(summary)
summary$fstatistic
coefficients(summary)
summary$coefficients
summary$coefficients[1,2]

######
summary$sigma
sqrt(error_variance)

######  Confidence intervals for the model parameters
confint(model, level = 0.95)

###### Measuring the degree of fit
r_squared <- SSR / SSY
r_squared


######
cor.test(data$growth, data$conc)
cor(data$growth, data$conc) ^ 2

###### Checking the model assumptions
par(mfrow=c(2,2))
plot(model)


###### Using the model for prediction
concentrations <- data.frame(conc = c(2.5, 6.5, 12.0))
predicted_growth <- predict(model, concentrations)
predicted_growth

###### Beyond Simple Linear Regression
######
######

data <- read.csv("decay.csv")
x <- data$x
y <- data$y
plot(x, y, pch=16, xlim=c(0, 33), ylim=c(0, 155))
model <- lm(y ~ x)
abline(model)

###### Polynomial regression

plot(model, which=c(1))

quadratic <- lm(y ~ x + I(x^2))
summary(quadratic)
######

linear <- lm(y ~ x)
anova(quadratic, linear)

######

plot(x, y, pch=16, xlim=c(0, 39), ylim=c(0, 155))
xv <- seq(0, 40, 0.1)
yv <- predict(quadratic, data.frame(x = xv))
lines(xv, yv)

###### Linearizing the model by transformation

exponential <- lm(log(y) ~ x)
summary(exponential)

