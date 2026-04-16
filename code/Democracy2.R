#################################################################################################
#  This is an empirical application of dynamic lineal panel data methods to the effect of democracy
#  on economic growth based on Acemoglu, Naidu, Restrepo and Robinson (2005), "Democracy Does
#  Cause Growth," forthcoming JPE
#
#  14.382 L9 MIT.  V. Chernozhukov and I. Fernandez-Val
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


# Updated on 04/24/2018


##################################################################################################




### set working directory
filepathIvan<-"/Users/Ivan/Dropbox/Shared/14.382"
filepathVictor<-"/Users/VC/Dropbox/TEACHING/14.382"
filepathMert<-"/Users/mertdemirer1/Dropbox (Personal)/14.382"
filepathVira<-"/Users/virasemenova/Dropbox/14.382"


filepathMatej<-"/Users/mcerman/Documents/2 Class/14.383/14.383_pset4/"
setwd(filepathMatej)

###  read-in TestScores data

library(foreign);
library(xtable);
library(plm);
library(gmm);
library(readstata13);

data <- read.dta13("data/democracy-balanced-l4.dta")
data <- pdata.frame(data, index = c("id","year"));

attach(data);


########## Descriptive statistics

options(digits=2);
dstat <- cbind(sapply(data[, c(5:6)], mean), sapply(data[,c(5:6)], sd), apply(data[data$dem==1, c(5:6)],2, mean), apply(data[data$dem==0, c(5:6)], 2, mean));
dstat <- rbind(dstat, c(nrow(data), nrow(data), sum(data$dem==1), sum(data$dem==0)));
dimnames(dstat) <- list(c("Democracy", "Log(GDP)", "Number Obs."), c("Mean", "SD", "Dem = 1", "Dem = 0"));
xtable(dstat)


########## Anderson-Hsiao estimation 


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


## Split sample bias correction (1 partition)

set.seed(383);
S1         <- 1;
N         <- length(levels(id));
acoeff.ah <- 0*coef(ah.fit);
alr.ah    <- 0;
for (s in 1:S1) {
  sample1   <- sample(N, ceiling(N/2), replace = FALSE);
  ah.fit1   <- pgmm(form.ah, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways", robust = TRUE);
  ah.fit2   <- pgmm(form.ah, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways", robust = TRUE);
  lr.ah1     <- coef(ah.fit1)[1]/(1-sum(coef(ah.fit1)[2:5]));
  lr.ah2     <- coef(ah.fit2)[1]/(1-sum(coef(ah.fit2)[2:5]));
  acoeff.ah <- acoeff.ah + ((coef(ah.fit1) + coef(ah.fit2))/2)/S1;
  alr.ah    <- alr.ah + ((lr.ah1 + lr.ah2)/2)/S1;
}
coefs.ah.jbc  <- 2*coef(ah.fit) - acoeff.ah;
lr.ah.jbc     <- 2*lr.ah - alr.ah;

## Split sample bias correction (5 partitions)

S1         <- 5;
N         <- length(levels(id));
acoeff.ah <- 0*coef(ah.fit);
alr.ah    <- 0;
for (s in 1:S1) {
  sample1   <- sample(N, ceiling(N/2), replace = FALSE);
  ah.fit1   <- pgmm(form.ah, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways", robust = TRUE);
  ah.fit2   <- pgmm(form.ah, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways", robust = TRUE);
  lr.ah1     <- coef(ah.fit1)[1]/(1-sum(coef(ah.fit1)[2:5]));
  lr.ah2     <- coef(ah.fit2)[1]/(1-sum(coef(ah.fit2)[2:5]));
  acoeff.ah <- acoeff.ah + ((coef(ah.fit1) + coef(ah.fit2))/2)/S1;
  alr.ah    <- alr.ah + ((lr.ah1 + lr.ah2)/2)/S1;
}
coefs.ah.jbc5  <- 2*coef(ah.fit) - acoeff.ah;
lr.ah.jbc5     <- 2*lr.ah - alr.ah;


########## Arellano-Bond estimation - All lags


form.ab <- lgdp ~ dem + lag(lgdp, 1:4) | lag(lgdp, 2:99) + lag(dem, 1:99);
ab.fit <- pgmm(form.ab, data, model = "twosteps", effect = "twoways" );
coefs.ab  <- coef(ab.fit); 
# se.ab     <- summary(ab.fit)$coef[ ,2];
HCV.coefs <- vcovHC(ab.fit, cluster = 'group');
cse.ab    <- sqrt(diag(HCV.coefs)); # Clustered std errors
Jtest.ab  <- sargan(ab.fit)$statistic;
Jdof.ab  <- sargan(ab.fit)$parameter;
lr.ab     <- coefs.ab[1]/(1-sum(coefs.ab[2:5]));
jac.lr <- c(1,rep(lr.ab,4))/(1-sum(coefs.ab[2:5]));
cse.lr.ab <- sqrt(t(jac.lr) %*% HCV.coefs[1:5,1:5] %*% jac.lr);

## Split sample bias correction (1 partition)

S2 <- 1;
acoeff.ab <- 0*coef(ab.fit);
alr.ab    <- 0;
for (s in 1:S2) {
  sample1   <- sample(N, ceiling(N/2), replace = FALSE);
  ab.fit1 <- pgmm(form.ab, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways" );
  ab.fit2 <- pgmm(form.ab, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways" );
  lr.ab1     <- coef(ab.fit1)[1]/(1-sum(coef(ab.fit1)[2:5]));
  lr.ab2     <- coef(ab.fit2)[1]/(1-sum(coef(ab.fit2)[2:5]));
  acoeff.ab <- acoeff.ab + ((coef(ab.fit1) + coef(ab.fit2))/2)/S2;
  alr.ab    <- alr.ab + ((lr.ab1 + lr.ab2)/2)/S2;
}
coefs.ab.jbc  <- 2*coef(ab.fit) - acoeff.ab;
lr.ab.jbc     <- 2*lr.ab - alr.ab;

## Split sample bias correction (5 partitions)

S2 <- 5;
acoeff.ab <- 0*coef(ab.fit);
alr.ab    <- 0;
for (s in 1:S2) {
  sample1   <- sample(N, ceiling(N/2), replace = FALSE);
  ab.fit1 <- pgmm(form.ab, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways" );
  ab.fit2 <- pgmm(form.ab, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways" );
  lr.ab1     <- coef(ab.fit1)[1]/(1-sum(coef(ab.fit1)[2:5]));
  lr.ab2     <- coef(ab.fit2)[1]/(1-sum(coef(ab.fit2)[2:5]));
  acoeff.ab <- acoeff.ab + ((coef(ab.fit1) + coef(ab.fit2))/2)/S2;
  alr.ab    <- alr.ab + ((lr.ab1 + lr.ab2)/2)/S2;
}
coefs.ab.jbc5  <- 2*coef(ab.fit) - acoeff.ab;
lr.ab.jbc5     <- 2*lr.ab - alr.ab;


########## Panel bootstrap std errors

R <- 500;

# Function to generate bootstrap data sets;

data.rg <- function(data, mle)
{
  N                  <- length(unique(data$id));
  T                  <- length(unique(data$year));
  ids                <- kronecker(sample.int(N, N, replace = TRUE), rep(1,T));
  data.b             <- data[(ids-1)*T + rep(c(1:T),N), ];
  data.b$id          <- kronecker(c(1:N), rep(1,T));
  data.b$year        <- rep(c(1987:2009),N);
  data.b             <- data.frame(data.b);
  data.b             <- pdata.frame(data.b, index = c("id","year"));         # reset indexes of the panel 
  return(data.b);
}


########## statistics to be computed in each bootstrap draw #####################################
boot.SE<- function(data, form.fe, form.ah, form.ab){

  ah.fit <- pgmm(form.ah, data, model = "twosteps", effect = "twoways", robust = TRUE);
  coefs.ah  <- coef(ah.fit); 
  lr.ah     <- coefs.ah[1]/(1-sum(coefs.ah[2:5]));
  
  N         <- length(unique(data$id));
  S1        <- 1;
  acoeff.ah <- 0*coef(ah.fit);
  alr.ah    <- 0;
  for (s in 1:S1) {
    sample1   <- sample(N, ceiling(N/2), replace = FALSE);
    ah.fit1   <- pgmm(form.ah, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways", robust = TRUE);
    ah.fit2   <- pgmm(form.ah, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways", robust = TRUE);
    lr.ah1     <- coef(ah.fit1)[1]/(1-sum(coef(ah.fit1)[2:5]));
    lr.ah2     <- coef(ah.fit2)[1]/(1-sum(coef(ah.fit2)[2:5]));
    acoeff.ah <- acoeff.ah + ((coef(ah.fit1) + coef(ah.fit2))/2)/S1;
    alr.ah    <- alr.ah + ((lr.ah1 + lr.ah2)/2)/S1;
  }
  coefs.ah.jbc  <- 2*coef(ah.fit) - acoeff.ah;
  lr.ah.jbc     <- 2*lr.ah - alr.ah;

  S1        <- 5;
  acoeff.ah <- 0*coef(ah.fit);
  alr.ah    <- 0;
  for (s in 1:S1) {
    sample1   <- sample(N, ceiling(N/2), replace = FALSE);
    ah.fit1   <- pgmm(form.ah, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways", robust = TRUE);
    ah.fit2   <- pgmm(form.ah, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways", robust = TRUE);
    lr.ah1     <- coef(ah.fit1)[1]/(1-sum(coef(ah.fit1)[2:5]));
    lr.ah2     <- coef(ah.fit2)[1]/(1-sum(coef(ah.fit2)[2:5]));
    acoeff.ah <- acoeff.ah + ((coef(ah.fit1) + coef(ah.fit2))/2)/S1;
    alr.ah    <- alr.ah + ((lr.ah1 + lr.ah2)/2)/S1;
  }
  coefs.ah.jbc5  <- 2*coef(ah.fit) - acoeff.ah;
  lr.ah.jbc5     <- 2*lr.ah - alr.ah;
  
  ab.fit <- pgmm(form.ab, data, model = "twosteps", effect = "twoways" );
  coefs.ab  <- coef(ab.fit); 
  lr.ab     <- coefs.ab[1]/(1-sum(coefs.ab[2:5]));

  S2        <- 1;
  acoeff.ab <- 0*coef(ab.fit);
  alr.ab    <- 0;
  for (s in 1:S2) {
    sample1   <- sample(N, ceiling(N/2), replace = FALSE);
    ab.fit1 <- pgmm(form.ab, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways" );
    ab.fit2 <- pgmm(form.ab, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways" );
    lr.ab1     <- coef(ab.fit1)[1]/(1-sum(coef(ab.fit1)[2:5]));
    lr.ab2     <- coef(ab.fit2)[1]/(1-sum(coef(ab.fit2)[2:5]));
    acoeff.ab <- acoeff.ab + ((coef(ab.fit1) + coef(ab.fit2))/2)/S2;
    alr.ab    <- alr.ab + ((lr.ab1 + lr.ab2)/2)/S2;
  }
  coefs.ab.jbc  <- 2*coef(ab.fit) - acoeff.ab;
  lr.ab.jbc     <- 2*lr.ab - alr.ab;
  
  S2        <- 5;
  acoeff.ab <- 0*coef(ab.fit);
  alr.ab    <- 0;
  for (s in 1:S2) {
    sample1   <- sample(N, ceiling(N/2), replace = FALSE);
    ab.fit1 <- pgmm(form.ab, data[as.double(id) %in% sample1, ], model = "twosteps", effect = "twoways" );
    ab.fit2 <- pgmm(form.ab, data[!(as.double(id) %in% sample1), ], model = "twosteps", effect = "twoways" );
    lr.ab1     <- coef(ab.fit1)[1]/(1-sum(coef(ab.fit1)[2:5]));
    lr.ab2     <- coef(ab.fit2)[1]/(1-sum(coef(ab.fit2)[2:5]));
    acoeff.ab <- acoeff.ab + ((coef(ab.fit1) + coef(ab.fit2))/2)/S2;
    alr.ab    <- alr.ab + ((lr.ab1 + lr.ab2)/2)/S2;
  }
  coefs.ab.jbc5  <- 2*coef(ab.fit) - acoeff.ab;
  lr.ab.jbc5     <- 2*lr.ab - alr.ab;
  
  return(c(coefs.ah[1:5], coefs.ah.jbc[1:5], coefs.ah.jbc5[1:5], coefs.ab[1:5], coefs.ab.jbc[1:5], coefs.ab.jbc5[1:5], 
             lr.ah, lr.ah.jbc, lr.ah.jbc5, lr.ab, lr.ab.jbc, lr.ab.jbc5));
}

library(boot); # library to do bootstrap with paralell computing


result.boot.SE <- boot(data = data, statistic=boot.SE, sim = "parametric", ran.gen = data.rg, mle = 0, form.fe = form.fe, form.ah = form.ah, 
                       form.ab = form.ab, parallel="multicore", ncpus = 16, R=R);


rsd <- function(x) { return((quantile(x,.75,na.rm=TRUE)-quantile(x,.25,na.rm=TRUE))/(qnorm(.75) - qnorm(.25)))} # robust estimator of std deviation based on IQR

# Robust bootstrap std errors;

result      <- structure(vapply(result.boot.SE$t, as.double, numeric(1)), dim=dim(result.boot.SE$t)); # transforms "Error in La.svd(x, nu, nv) : error code 1 from Lapack routine 'dgesdd'\n" to NA
bse.ah      <- apply(result[,1:5], 2, rsd);
bse.ah.jbc  <- apply(result[,6:10], 2, rsd);
bse.ah.jbc5 <- apply(result[,11:15], 2, rsd);
bse.ab      <- apply(result[,16:20], 2, rsd);
bse.ab.jbc  <- apply(result[,21:25], 2, rsd);
bse.ab.jbc5 <- apply(result[,26:30], 2, rsd);

bse.lr.ah     <- rsd(result[,31]);
bse.lr.ah.jbc <- rsd(result[,32]);
bse.lr.ah.jbc5 <- rsd(result[,33]);
bse.lr.ab     <- rsd(result[,34]);
bse.lr.ab.jbc <- rsd(result[,35]);
bse.lr.ab.jbc5 <- rsd(result[,36]);


######## Table of results;

options(digits=2);
table.all <- matrix(NA, nrow = 18, ncol = 6, dimnames = list(c("Democracy", "CSE", "BSE", "L1.log(gdp)",  "CSE1", "BSE1", "L2.log(gdp)",  "CSE2", "BSE2","L3.log(gdp)",  "CSE3", "BSE3", "L4.log(gdp)",  "CSE4", "BSE4", "LR-Democracy","CSE5","BSE5"), c("GMM1", "GMM1-BC", "GMM1-BC5", "GMM2", "GMM2-BC", "GMM2-BC5")));


table.all[c(1,4,7,10,13), 1] <- coefs.ah[1:5];
table.all[c(2,5,8,11,14), 1] <- cse.ah[1:5];
table.all[c(3,6,9,12,15), 1] <- bse.ah[1:5];

table.all[c(1,4,7,10,13), 2] <- coefs.ah.jbc[1:5];
#table.all[c(2,5,8,11,14), 2] <- cse.ah[1:5];
table.all[c(3,6,9,12,15), 2] <- bse.ah.jbc[1:5];

table.all[c(1,4,7,10,13), 3] <- coefs.ah.jbc5[1:5];
#table.all[c(2,5,8,11,14), 3] <- cse.ah[1:5];
table.all[c(3,6,9,12,15), 3] <- bse.ah.jbc5[1:5];

table.all[c(1,4,7,10,13), 4] <- coefs.ab[1:5];
table.all[c(2,5,8,11,14), 4] <- cse.ab[1:5];
table.all[c(3,6,9,12,15), 4] <- bse.ab[1:5];

table.all[c(1,4,7,10,13), 5] <- coefs.ab.jbc[1:5];
#table.all[c(2,5,8,11,14), 5] <- cse.ab[1:5];
table.all[c(3,6,9,12,15), 5] <- bse.ab.jbc[1:5];

table.all[c(1,4,7,10,13), 6] <- coefs.ab.jbc5[1:5];
#table.all[c(2,5,8,11,14), 6] <- cse.ab[1:5];
table.all[c(3,6,9,12,15), 6] <- bse.ab.jbc5[1:5];

table.all[16, 1] <- lr.ah;
table.all[17, 1] <- cse.lr.ah;
table.all[18, 1] <- bse.lr.ah;

table.all[16, 2] <- lr.ah.jbc;
#table.all[17, 2] <- cse.lr.ah;
table.all[18, 2] <- bse.lr.ah.jbc;

table.all[16, 3] <- lr.ah.jbc5;
#table.all[17, 3] <- cse.lr.ah;
table.all[18, 3] <- bse.lr.ah.jbc5;

table.all[16, 4] <- lr.ab;
table.all[17, 4] <- cse.lr.ab;
table.all[18, 4] <- bse.lr.ab;

table.all[16, 5] <- lr.ab.jbc;
#table.all[17, 5] <- cse.lr.ab;
table.all[18, 5] <- bse.lr.ab.jbc;

table.all[16, 6] <- lr.ab.jbc5;
#table.all[17, 6] <- cse.lr.ab;
table.all[18, 6] <- bse.lr.ab.jbc5;

table.all[1, ] <- 100 * table.all[1, ];
table.all[2, ] <- 100 * table.all[2, ];
table.all[3, ] <- 100 * table.all[3, ];
table.all[16, ] <- 100 * table.all[16, ];
table.all[17, ] <- 100 * table.all[17, ];
table.all[18, ] <- 100 * table.all[18, ];


results_table = xtable(table.all, digits=2)

print(results_table, file = "outputs/democracy_results.tex")




