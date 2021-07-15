## ANOVA, Linear Regression comparison

# needed packages and libraries 
install.packages('AICcmodavg')
library(AICcmodavg)


# Make fake data for ANOVA
set.seed(0)
data = data.frame(Sample=1:50,
                  gene_exp = c( rnorm(25, 47.5), # treated with RNAi
                                rnorm(25, 50)), # RNAi-negative control
                  RNAi = c( rep(1,25), # RNAi treated are first 25
                            rep(0,25)),
                  Drug = rep(c(0,1),25)) # drug treatment is not associated with a drop in expression

data$RNAi = factor(data$RNAi)
data$Drug = factor(data$Drug)

print(data)
str(data)

boxplot(gene_exp ~ Drug + RNAi, data)

# Analyze data with ANOVA model
one.way <- aov(gene_exp ~ Drug, data)
summary(one.way)

two.way <- aov(gene_exp ~ Drug + RNAi, data)
summary(two.way)

two.way_interaction <- aov(gene_exp ~ Drug + RNAi + Drug*RNAi, data)
summary(two.way_interaction)

#Find the best-fit ANOVA model
model.set <- list(one.way, two.way, two.way_interaction)
model.names <- c("one.way", "two.way", "interaction")

aictab(model.set, modnames = model.names)

# Check for homoscedasticity
plot(one.way)
plot(two.way)
plot(two.way_interaction)

par(mfrow=c(2,2))
plot(two.way)
par(mfrow=c(1,1))

# Post hoc test
tukey.plot.test <-TukeyHSD(two.way_interaction)
plot(tukey.plot.test, las = 1)
tukey.plot.test

# Make fake data for Linear Regression
set.seed(0)
data2 = data.frame(Sample=1:50,
                  gene_exp = c( rnorm(25, 47.5), # treated with RNAi
                                rnorm(25, 50)), # RNAi-negative control
                  RNAi = seq(from = 0, to = 300, length.out = 50), # RNAi treatment dosage 
                  Drug = seq(from = 0, to = 500, length.out = 50)) # Drug treatment dosage 
print(data2)
str(data2)

plot(data2$Drug, data2$gene_exp, col = 'blue', xlab = 'Predictors', ylab = 'Gene Expression')
points(data2$RNAi, data2$gene_exp, col = 'red')
legend("topright", c("Drug", "RNAi"), pch=c(19,19), col = c("blue", "red"))


# Analyze data with Linear Regression model
Data2_lm1 <- lm(gene_exp ~ RNAi, data2)
summary(Data2_lm1)

Data2_lm2 <- lm(gene_exp ~ RNAi + Drug, data2)
summary(Data2_lm2)

Data2_lm3 <- lm(gene_exp ~ RNAi + Drug + RNAi*Drug, data2)
summary(Data2_lm3)

#Find the best Linear Regression model
ANOVA_fit1 <- anova(Data2_lm1, Data2_lm2)
ANOVA_fit1

ANOVA_fit2 <- anova(Data2_lm1, Data2_lm3)
ANOVA_fit2

ANOVA_fit3 <- anova(Data2_lm2, Data2_lm3)
ANOVA_fit3

ANOVA_fit4 <- anova(Data2_lm1, Data2_lm2, Data2_lm3)
ANOVA_fit4

model.set <- list(Data2_lm1, Data2_lm2, Data2_lm3)
model.names <- c("lm1", "lm2", "lm3")
aictab(model.set, modnames = model.names)

# Check for homoscedasticity
par(mfrow=c(2,2))
plot(Data2_lm1)
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(Data2_lm2)
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(Data2_lm3)
par(mfrow=c(1,1))