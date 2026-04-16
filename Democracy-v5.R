#################################################################################################
#  This is an empirical application of lineal panel data methods to the effect of democracy
#  on economic growth based on Acemoglu, Naidu, Restrepo and Robinson (2005), "Democracy Does
#  Cause Growth," forthcoming JPE
#
#  14.382 L8 MIT.  V. Chernozhukov and I. Fernandez-Val
#
# Data source: Daron Acemoglu (MIT), N = 147 countries, T = 23 years (1987-2009)
# include data for four lags of dependent variable, balanced panel
#
# Description of the data: the sample selection and variable contruction follow
#
# The variables in the data set include:
#
# country_name    = Country name
# wbcode          = World Bank country code
# year            = Year 
# id              = Generated numeric country code
# dem             = Democracy measure by ANRR
# lgdp            = log of GDP per capita in 2000 USD from World Bank

#################################################################################################


# Updated on 04/22/2018


##################################################################################################




library(foreign);
library(xtable);
library(plm);
library(gmm);
library(readstata13);

data <- read.dta13("/Users/abaraban/classes/metrics/14.383_pset4/democracy-balanced-l4.dta")
data <- pdata.frame(data, index = c("id","year"));

attach(data);


########## Descriptive statistics

options(digits=2);
dstat <- cbind(sapply(data[, c(5:6)], mean), sapply(data[,c(5:6)], sd), apply(data[data$dem==1, c(5:6)],2, mean), apply(data[data$dem==0, c(5:6)], 2, mean));
dstat <- rbind(dstat, c(nrow(data), nrow(data), sum(data$dem==1), sum(data$dem==0)));
dimnames(dstat) <- list(c("Democracy", "Log(GDP)", "Number Obs."), c("Mean", "SD", "Dem = 1", "Dem = 0"));
xtable(dstat);

########## OLS estimation

form.ols <- lgdp ~ dem + lag(lgdp, 1:4) + factor(year) -1 ;

ols.fit   <- plm(form.ols, data, model = "pooling", index = c("id","year"));
coefs.ols <- coef(ols.fit); 
se.ols    <- summary(ols.fit)$coef[ ,2];
HCV.coefs <- vcovHC(ols.fit, method = 'arellano', type = 'HC0', cluster = 'group');
cse.ols   <- sqrt(diag(HCV.coefs)); # Clustered std errors
lr.ols     <- coefs.ols[1]/(1-sum(coefs.ols[2:5]));
jac.lr <- c(1,rep(lr.ols,4))/(1-sum(coefs.ols[2:5]));
cse.lr.ols <- sqrt(t(jac.lr) %*% HCV.coefs[1:5,1:5] %*% jac.lr);


########## Fixed Effects estimation

form.fe <- lgdp ~ dem + lag(lgdp, 1:4) - 1;

fe.fit    <- plm(form.fe, data, model = "within", effect = "twoways", index = c("id","year"));
coefs.fe  <- coef(fe.fit); 
se.fe     <- summary(fe.fit)$coef[ ,2];
HCV.coefs <- vcovHC(fe.fit, cluster = 'group');
cse.fe    <- sqrt(diag(HCV.coefs)); # Clustered std errors
lr.fe     <- coefs.fe[1]/(1-sum(coefs.fe[2:5]));
jac.lr <- c(1,rep(lr.fe,4))/(1-sum(coefs.fe[2:5]));
cse.lr.fe <- sqrt(t(jac.lr) %*% HCV.coefs[1:5,1:5] %*% jac.lr);

########## First Differences estimation


form.fd <- lgdp ~ dem  + lag(lgdp, 1:4) + factor(year) - 1;

fd.fit    <- plm(form.fd, data, model = "fd", index = c("id","year"));
coefs.fd  <- coef(fd.fit); 
se.fd     <- summary(fd.fit)$coef[ ,2];
HCV.coefs <- vcovHC(fd.fit, method = 'arellano', type = 'HC0', cluster = 'group');
cse.fd    <- sqrt(diag(HCV.coefs)); # Clustered std errors
lr.fd     <- coefs.fd[1]/(1-sum(coefs.fd[2:5]));
jac.lr <- c(1,rep(lr.fd,4))/(1-sum(coefs.fd[2:5]));
cse.lr.fd <- sqrt(t(jac.lr) %*% HCV.coefs[1:5,1:5] %*% jac.lr);

########## GMM-FD2 ("Arellano-Bond") estimation 


form.ab <- lgdp ~ dem + lag(lgdp, 1:4) | lag(lgdp, 2:99) + lag(dem, 1:99);
ab.fit <- pgmm(form.ab, data, model = "twosteps", effect = "twoways", robust = TRUE );
coefs.ab  <- coef(ab.fit); 
# se.ab     <- summary(ab.fit)$coef[ ,2];
HCV.coefs <- vcovHC(ab.fit, cluster = 'group');
cse.ab    <- sqrt(diag(HCV.coefs)); # Clustered std errors
Jtest.ab  <- sargan(ab.fit)$statistic;
Jdof.ab  <- sargan(ab.fit)$parameter;
lr.ab     <- coefs.ab[1]/(1-sum(coefs.ab[2:5]));
jac.lr <- c(1,rep(lr.ab,4))/(1-sum(coefs.ab[2:5]));
cse.lr.ab <- sqrt(t(jac.lr) %*% HCV.coefs[1:5,1:5] %*% jac.lr);

########## GMM-FD1 ("Anderson-Hsiao") estimation 


form.ah <- lgdp ~ dem + lag(lgdp, 1:4) | lag(lgdp, 2) + lag(dem, 1);
ah.fit <- pgmm(form.ah, data, model = "twosteps", effect = "twoways", robust = TRUE);
coefs.ah  <- coef(ah.fit); 
# se.ah     <- summary(ah.fit)$coef[ ,2];
HCV.coefs <- vcovHC(ah.fit, cluster = 'group');
cse.ah    <- sqrt(diag(HCV.coefs)); # Clustered std errors
Jtest.ah  <- sargan(ah.fit)$statistic;
Jdof.ah  <- sargan(ah.fit)$parameter;
lr.ah     <- coefs.ah[1]/(1-sum(coefs.ah[2:5]));
jac.lr <- c(1,rep(lr.ah,4))/(1-sum(coefs.ah[2:5]));
cse.lr.ah <- sqrt(t(jac.lr) %*% HCV.coefs[1:5,1:5] %*% jac.lr);


########## Panel bootstrap std errors

set.seed(8);
N <- length(unique(id));
T <- length(unique(year));
R <- 500;

coefs.ols.b <- matrix(0, ncol = length(coefs.ols), nrow = R);
coefs.fe.b  <- matrix(0, ncol = length(coefs.fe), nrow = R);
coefs.fd.b  <- matrix(0, ncol = length(coefs.fd), nrow = R);
coefs.ab.b  <- matrix(0, ncol = length(coefs.ab), nrow = R);
coefs.ah.b  <- matrix(0, ncol = length(coefs.ah), nrow = R);

lr.ols.b <- rep(0,nrow=R);
lr.fe.b <- rep(0,nrow=R);
lr.fd.b <- rep(0,nrow=R);
lr.ah.b <- rep(0,nrow=R);
lr.ab.b <- rep(0,nrow=R);

for (r in 1:R) {;
                ids                <- kronecker(sample.int(N, N, replace = TRUE), rep(1,T));
                data.b             <- data[(ids-1)*T + rep(c(1:T),N), ];
                data.b$id          <- kronecker(c(1:N), rep(1,T));
                data.b$year        <- rep(c(1987:2009),N);
                data.b             <- data.frame(data.b);
                data.b             <- pdata.frame(data.b, index = c("id","year"));         # reset indexes of the panel                        

                ols.fit.b           <- plm(form.ols, data.b, model = "pooling", index = c("id","year"));
                coefs.ols.b[r, ]    <- coef(ols.fit.b);
                lr.ols.b[r]         <- coef(ols.fit.b)[1]/(1-sum(coef(ols.fit.b)[2:5]));

                fe.fit.b           <- plm(form.fe, data.b, model = "within", effect = "twoways", index = c("id","year"));
                coefs.fe.b[r, ]    <- coef(fe.fit.b);
                lr.fe.b[r]         <- coef(fe.fit.b)[1]/(1-sum(coef(fe.fit.b)[2:5]));
                
                fd.fit.b           <- plm(form.fd, data.b, model = "fd", index = c("id","year"));
                coefs.fd.b[r, ]    <- coef(fd.fit.b);
                lr.fd.b[r]         <- coef(fd.fit.b)[1]/(1-sum(coef(fd.fit.b)[2:5]));
                
                ah.fit.b           <- pgmm(form.ah, data.b, model = "twosteps", effect = "twoways", robust = TRUE);
                coefs.ah.b[r, ]    <- coef(ah.fit.b);
                lr.ah.b[r]         <- coef(ah.fit.b)[1]/(1-sum(coef(ah.fit.b)[2:5]));
                
                ab.fit.b             <- pgmm(form.ab, data.b, model = "twosteps", effect = "twoways" );
                coefs.ab.b[r, ]      <- coef(ab.fit.b);
                lr.ab.b[r]           <- coef(ab.fit.b)[1]/(1-sum(coef(ab.fit.b)[2:5]));
                
};

rsd <- function(x) { return((quantile(x,.75)-quantile(x,.25))/(qnorm(.75) - qnorm(.25)))} # robust estimator of std deviation based on IQR

# Robust bootstrap std errors;

bse.ols     <- apply(coefs.ols.b, 2, rsd);
bse.fe      <- apply(coefs.fe.b, 2, rsd);
bse.fd      <- apply(coefs.fd.b, 2, rsd);
bse.ah      <- apply(coefs.ah.b, 2, rsd);
bse.ab      <- apply(coefs.ab.b, 2, rsd);

bse.lr.ols    <- rsd(lr.ols.b);
bse.lr.fe     <- rsd(lr.fe.b);
bse.lr.fd     <- rsd(lr.fd.b);
bse.lr.ah     <- rsd(lr.ah.b);
bse.lr.ab     <- rsd(lr.ab.b);


######## Table of results;

options(digits=3);
table.all <- matrix(NA, nrow = 21, ncol = 5, dimnames = list(c("Democracy", "CSE", "BSE", "L1.log(gdp)",  "CSE1", "BSE1", "L2.log(gdp)",  "CSE2", "BSE2","L3.log(gdp)",  "CSE3", "BSE3", "L4.log(gdp)",  "CSE4", "BSE4", "LR-Democracy","CSE5","BSE5","J-test", "p-val","dof"), c("Pooled", "FD", "FE", "GMM-FD1", "GMM-FD2")));

table.all[c(1,4,7,10,13), 1] <- coefs.ols[1:5];
table.all[c(2,5,8,11,14), 1] <- cse.ols[1:5];
table.all[c(3,6,9,12,15), 1] <- bse.ols[1:5];

table.all[c(1,4,7,10,13), 2] <- coefs.fd[1:5];
table.all[c(2,5,8,11,14), 2] <- cse.fd[1:5];
table.all[c(3,6,9,12,15), 2] <- bse.fd[1:5];

table.all[c(1,4,7,10,13), 3] <- coefs.fe[1:5];
table.all[c(2,5,8,11,14), 3] <- cse.fe[1:5];
table.all[c(3,6,9,12,15), 3] <- bse.fe[1:5];

table.all[c(1,4,7,10,13), 4] <- coefs.ah[1:5];
table.all[c(2,5,8,11,14), 4] <- cse.ah[1:5];
table.all[c(3,6,9,12,15), 4] <- bse.ah[1:5];

table.all[c(1,4,7,10,13), 5] <- coefs.ab[1:5];
table.all[c(2,5,8,11,14), 5] <- cse.ab[1:5];
table.all[c(3,6,9,12,15), 5] <- bse.ab[1:5];

table.all[16, ] <- c(lr.ols, lr.fd, lr.fe, lr.ah, lr.ab);
table.all[17, ] <- c(cse.lr.ols, cse.lr.fd, cse.lr.fe, cse.lr.ah, cse.lr.ab);
table.all[18, ] <- c(bse.lr.ols, bse.lr.fd, bse.lr.fe, bse.lr.ah, bse.lr.ab);

table.all[1, ] <- 100 * table.all[1, ];
table.all[2, ] <- 100 * table.all[2, ];
table.all[3, ] <- 100 * table.all[3, ];
table.all[16, ] <- 100 * table.all[16, ];
table.all[17, ] <- 100 * table.all[17, ];
table.all[18, ] <- 100 * table.all[18, ];


table.all[19 ,4] <- Jtest.ah;
table.all[20 ,4] <- 1 - pchisq(Jtest.ah, Jdof.ah);
table.all[21 ,4] <- as.integer(Jdof.ah);

table.all[19 ,5] <- Jtest.ab;
table.all[20 ,5] <- 1 - pchisq(Jtest.ab, Jdof.ab);
table.all[21 ,5] <- as.integer(Jdof.ab);



xtable(table.all, digits=2);



# Long run effects;


# lr.fe <- 100*coefs.fe[1]/(1 - sum(coefs.fe[2:5]));
# lr.ab <- 100*coefs.ab[1]/(1 - sum(coefs.ab[2:5]));
# lr.ah <- 100*coefs.ah[1]/(1 - sum(coefs.ah[2:5]));
# 
# print(c(lr.fe, lr.ab, lr.ah));