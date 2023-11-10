#### 
# Example of the code to estimate the models proposed  in the paper titled 
# "An analysis of the effect of streaming on civic participation through a causal hidden Markov model" written by 
# - F. Bartolucci (University of Perugia, IT)
# - D. Favaro (University of Padova, IT)
# - F. Pennoni (University of Milano-Bicocca, IT)
# - D. Sciulli (University of Chieti-Pescara, IT)
####
#### Load pseudo data as an example  ####
rm(list=ls())
if(!("LMest"%in%installed.packages())) install.packages("nnet")
if(!("LMest"%in%installed.packages())) install.packages("LMest")
#
load("dtt1.Rdata")
names(dtt1)
summary(dtt1)
# ncdsid: id 
# year: year of the survey
# stream: treatment (not streamed,  low-ability, average-ability, high-ability)
# voted: binary response variable (0 for no, 1 for yes)
# org_poli: binary response variable (0 for no, 1 for yes)
# org_volunt: binary response variable (0 for no, 1 for yes)
# female: gender (binary variable, post)
# read: reading score  (pretreatment covariate)
# math: score in maths (pretreatment covariate)
# region1: region of residence (posttreatment covariate)
# married: marital status (0 for no, 1 for yes, none for missing (posttreatment covariate) missing)

### Code to estimate PS weights
require(nnet)
#  create a cross-sectional dataset
ind <- which(dtt1$year==1991)
data2 <- dtt1[ind,]
names(data2) 
n <- dim(data2)[1]
n
resp2 <- multinom(data2$stream ~    
                    data2$math  +
                    data2$read)

# matrix of the predicted values (multinomial probabilities)
P <- summary(resp2)$fitted.values 
# vector of the individual  weights 
w <- 1/P[cbind(1:n,data2$stream)]    
# rescaling weights 
w <- w/sum(w)*n; 
#
summary(w)

### Code employed to estimate the proposed models ###

require(LMest)
####  Model selection: basic HM model with estimated  weights ####
set.seed(143)
mod1w <- lmestSearch(responsesFormula = 
                      voted+
                      org_polit +
                      org_volunt   ~  NULL,
                      weights  = w,
                      index = c("ncdsid","year"),
                      data = dtt1,
                      k = 1:2,
                      nrep = 10,
                      modBasic = 0, 
                      version = "categorical", 
                      seed = 123)
summary(mod1w)

####  HM model estimation with selected number of states, with weights, treatment and posttreatment covariates #####
mod2 <- lmest(voted +
                org_polit +
                org_volunt   ~ NULL,
                latentFormula = ~ 
                I(0 + (stream == 'low-ability')) +
                I(0 + (stream == 'average-ability')) +
                I(0 + (stream == 'high-ability')) +
                I(0 + (married == 1)) + 
                I(0 + (married == "none")) + 
                I(0 + (region1 == 1)), 
                paramLatent = "difflogit",
                index = c("ncdsid","year"), data = dtt1,
                k = 2:2, weights  = w,
                modBasic = 1, out_se = TRUE)


#### Print results ####
summary(mod2)

# Estimated conditional probabilities
round(mod2$Psi,3)

# Estimated  averaged initial probabilities 
round(apply(mod2$Piv,2,mean),3)

# Estimated averaged transition matrix
round(apply(mod2$PI[, , , 2], c(1, 2), mean),3)

# Results on the initial probabilities
# Estimated regression coefficients for the covariates (Be) 
# and standard errors (seBe)
TabBeYY <-  cbind(mod2$Be,mod2$seBe)
colnames(TabBeYY)<-c("Be", "seBe")
round(TabBeYY,3)

# Results on the transition probabilities
# Estimated regression coefficients for the covariates (Ga) 
TabGa1Y <- round(cbind(mod2$Ga[[2]],mod2$seGa[[2]]),3)
colnames(TabGa1Y) <- c("Ga1", "seGa1")
TabGa1Y

# Decoded states 
deco2 <-lmestDecoding(mod2)
head(deco2$Ug)

