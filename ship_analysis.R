library(fda)
#library(fdapace)
load("~/Downloads/df.RData")
names(df)
df2 <- df[, -c(1,19, 21, 22,23,24,25, 26, 29)]

basis <- create.bspline.basis(rangeval=c(0,1),nbasis=100)
treatment <- unique(df2[, c(16, 18)])
df2 <- df2[, -16]
#Xsp <- smooth.basis(argvals=abscissa, y=Xobs0, fdParobj=basis)
#Xsp0bis <- eval.fd(abscissa, Xsp$fd)


abscissa <- seq(0, 1, length.out = 120)


vessels <- NULL
i = 1
 
chosen <- which(table(df2$VN) > 110)
X <- array(0, dim = c(120, length(unique(df2$VN)),  18 ))
df2 <- df2[df2$VN %in% unique(df2$VN)[chosen],]
for (vessel in unique(df2$VN)){
  tryCatch({
    Xi <- array(0, c(120, 18))
    
    for (p in seq(1, 18)[-17]){
      df.aux <- df2[df2$VN == vessel, c(p, 20)]
      if (p == 1){
        Xsp <- smooth.basis(argvals=df.aux$percent_miles_hat,
                            y=unlist(df.aux[,1]), fdParobj=basis)
        Xi[, 1] = as.vector(eval.fd(abscissa, Xsp$fd))
        Xi[, 1]  = ( Xi[, 1] - mean(Xi[, 1])) / var(Xi[, 1])
        Xi[, 18] = as.vector(eval.fd(abscissa, Xsp$fd, Lfdobj = 1))
        Xi[, 18] = (Xi[, 18] - mean(Xi[, 18]) ) / var(Xi[, 18])
      }
      else{
        Xsp <- smooth.basis(argvals=df.aux$percent_miles_hat,
                            y=unlist(df.aux[,1]), fdParobj=basis)
        Xi[ , p] = as.vector(eval.fd(abscissa, Xsp$fd))
        Xi[ , p] = ( Xi[ , p] - mean(Xi[ , p]) ) / var(Xi[ , p])
      }
    }
    vessels <- c(vessels, vessel)
    print(i)
    X[,i,] = Xi
    i  = i + 1
  }, error=function(e){print(paste0(e, i, p))}
  )
}
X <- X[, 1:(i-1), ]
matplot(t(X[,,1]), type='l')
#colnames(X) = abscissa
#rownames(X) = vesses
ships.fd <- Data2fd(y = X, argvals = abscissa, basisobj = basis)
ships.fd$fdnames$reps <- vessels
save(ships.fd, file="shipsfd.RData")
setwd("~/PycharmProjects/SDA-Functional-Control-Charts/")
load("shipsfd.RData")
ships.fd$fdnames$values <- # TODO
# ships.fd$fdnames = as.list(vessels)
library(fda)

x11()
plot.fd(ships.fd[1:30])
ships.fd$coefs
pca_W.1 <- pca.fd(ships.fd[1:60],nharm=50,centerfns=TRUE)
dim(pca_W.1$scores)


plot(pca_W.1$values[1:35],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:35]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

dim(apply(pca_W.1$scores[, 1:7,], MARGIN = 1, FUN = sum))
X <- array(0, c(60,50))
for (i in 1:60){
  X[i,] <- apply(pca_W.1$scores[i,1:50,], MARGIN = 1, FUN = sum)
}

fit <- lm (y[1:60]~ X)
summary(fit)



# prepare cotrol chart and other stuff. 
T2 <- X^2 # square element - wise
for ( i in ncol(T2)){
  T2[,i] = T2[,i] / pca_W.1$values[i]
}
# set up upper control limits. 4 sdev for 93.75% of data
UCLs <- apply(T2, MARGIN = 2, FUN = sd)

# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first 
# N-1=34 are non-null
"N-1 refers to how many observations at each time point"
plot(pca_W.1$values[1:35],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:35]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))

## Test set --------
# PREPARE TEST SET 
load("y_C02_emissions.RData")
indices_train <- 1:60
indices_test <- (1:dim(y)[1])[-indices_train]
y_na <- c(292, 508, 705)
indices_test <- indices_test[-y_na]
keep_pcs = 50
Z_test <- array(0, dim = c(length(indices_test), keep_pcs) )
problematic.ships <- NULL
k = 0
for (i in indices_test){
  tryCatch({
    #print(i)
    cur.curve <- ships.fd[ i ]
    scores <- predictFPCA(pca_W.1, cur.curve, nharm=50, nvar=18)
    Z_test[k, ] = apply(scores[1:keep_pcs,], MARGIN = 1, FUN = sum)
    k = k+1
    
    
  },
  error = function(e){
    print(e)
    problematic.ships <- c(problematic.ships, i)
    print(i)
  })
  
}

Z_test <- as.data.frame(Z_test)
rownames(Z_test) <- indices_test



# Obtain values and fit --------
Z <- pca_W.1$scores[, 1:5]
y <- aggregate(df[,29], by = list(df$VN), FUN = sum)
y = y[y$Group.1 %in% vessels,]

y$VN == ships.fd[-y_na]$fdnames$reps # check
colnames(y) <- c("VN", "CO2.emissions")
y_na <- which(is.na(y$CO2.emissions))
y <- y[-y_na,]
save(y, file = "y_C02_emissions.RData")
fit <- lm(y ~ Z)
summary(fit)






# Appendix -----------
Data2fd(argvals = abscissa, y = samples_sog)

# pass to object 
basis <- create.bspline.basis(rangeval=c(0,1),nbasis=11)
data.ships.fd <- Data2fd(y = df2.2 ,argvals = df2.2[,-1],
                         basisobj = basis, fdnames = c('percent_miles', 'SOG', 'VN'))
plot.fd(data.ships.fd, main="B-splines")

#### CANADA WEATHER ----
data_W <- CanadianWeather$dailyAv
time <- 1:365
basis.3 <- create.fourier.basis(rangeval=c(0,365),nbasis=109)
data_W.fd.3 <- Data2fd(y = data_W,argvals = time,basisobj = basis.3)
plot.fd(data_W.fd.3)

pca_W.1 <- pca.fd(data_W.fd.3,nharm=25,centerfns=TRUE)
pca_W.1$scores
dim(pca_W.1$scores)
Z <- c(pca_W.1$scores[1, 1:6, 1], pca_W.1$scores[1, 1:6, 2], 
  pca_W.1$scores[1, 1:6, 3])
for (i in 2:10){
  Z <- rbind(Z, c(pca_W.1$scores[i, 1:4, 1], pca_W.1$scores[i, 1:4, 2], 
                  pca_W.1$scores[i, 1:3, 3]) )
  
}

df <- data.frame(scores1=pca)
# scree plot
# pca.fd computes all the 365 eigenvalues, but only the first 
# N-1=34 are non-null
"N-1 refers to how many observations at each time point"
plot(pca_W.1$values[1:35],xlab='j',ylab='Eigenvalues')
plot(cumsum(pca_W.1$values)[1:35]/sum(pca_W.1$values),xlab='j',ylab='CPV',ylim=c(0.8,1))


unique(CanadianWeather$region)
a <- glm(CanadianWeather$region=="Continental" ~ pca_W.1$scores[, 1:5] , family = "binomial")
summary(a)
predict(a, data.frame(pca_W.1$scores[, 1:5]))
a <- multinom(as.factor(CanadianWeather$region) ~ pca_W.1$scores[, 1:5])
summary(a)
?multinom

library(fda)

indices_train <- c(9, 23, 26, 20, 4, 17, 25, 19, 1, 18, 8, 6, 34, 13)
data_W <- CanadianWeather$dailyAv
time <- 1:365
basis.3 <- create.fourier.basis(rangeval=c(0,365),nbasis=109)
data_W.fd.3 <- Data2fd(y = data_W,argvals = time,basisobj = basis.3)

pca_W.1 <- pca.fd(data_W.fd.3[indices_train],nharm=25,centerfns=TRUE)

# obtain dataframe with features
i <- indices_train[1]

Z <- cbind(pca_W.1$scores[i, 1:6, 1], pca_W.1$scores[i, 1:6, 2], pca_W.1$scores[i, 1:4, 3])
Z <- apply(Z, MARGIN = 1, FUN = sum)
for (i in 2:length(indices_train)){
  aux <- cbind(pca_W.1$scores[i, 1:6, 1], pca_W.1$scores[i, 1:6, 2], 
                  pca_W.1$scores[i, 1:6, 3]) )  
  Z <- rbind(Z, apply(aux, MARGIN = 1, FUN = sum))
  
}
Z <- as.data.frame(Z)

### PREDICT SCORES ON NEW FUNCTIONAL SAMPLE
predictFPCA <- function(pca_obj, fdobj, nharm=25, nvar=3){
  # Center (crucial!)
  fdobj$coefs <- fdobj$coefs - pca_obj$meanfd$coefs[,1,]
  
  harmscr  <- array(0, c(nharm, nvar))
  coefarray <- fdobj$coefs
  harmcoefarray <- (pca_obj$harmonics)$coefs
  basisobj <- fdobj$basis
  for (j in 1:nvar) {
    fdobjj  <- fd(as.matrix(coefarray[,j]), basisobj)
    harmfdj <- fd(as.matrix(harmcoefarray[,,j]), basisobj)
    harmscr[,j] <- inprod(fdobjj, harmfdj)
  }
  return (harmscr)
}
fdobj <- data_W.fd.3[1]
predictFPCA(pca_W.1, fdobj, nharm=25, nvar=3)

# PREPARE TEST SET
# PREPARE TEST SET 
test_ind <- (1:dim(data_W)[2])[-indices_train]
Z_test <- array(0, c(length(test_ind), 6))
for (i in 1:length(test_ind)){
  cur.curve <- data_W.fd.3[ test_ind[i] ]
  scores <- predictFPCA(pca_W.1, cur.curve, nharm=25, nvar=3)
  Z_test[i, ] = apply(scores[1:6,], MARGIN = 1, FUN = sum)
}
  
Z_test <- as.data.frame(Z_test)
rownames(Z_test) <- test_ind
  
  






## PREDICT SCORES ON NEW FUNCTIONAL SAMPLE
predictFPCA <- function(pca_obj, fdobj, nharm=25, nvar=3){
  # Center
  fdobj <- center.fd(fdobj)
  
  harmscr  <- array(0, c(1, nharm, nvar))
  coefarray <- fdobj$coefs
  harmcoefarray <- pca_obj$harmonics$coefs
  basisobj <- fdobj$basis
  for (j in 1:nvar) {
    fdobjj  <- fd(as.matrix(coefarray[,,j]), basisobj)
    harmfdj <- fd(as.matrix(harmcoefarray[,,j]), basisobj)
    harmscr[,j] <- inprod(fdobjj, harmfdj)
  }
  return (harmscr)
}
predictFPCA(pca_W.1, (ships.fd[1]), nharm=50, nvar=18)[,1]
(pca_W.1$scores[1,,1])

# PREPARE TEST SET 
test_ind <- (1:dim(data_W)[2])[-indices_train]
Z_test <- array(0, c(length(test_ind), 6))
for (i in 1:length(test_ind)){
  cur.curve <- data_W.fd.3[ test_ind[i] ]
  scores <- predictFPCA(pca_W.1, cur.curve, nharm=25, nvar=3)
  Z_test[i, ] = apply(scores[1:6,], MARGIN = 1, FUN = sum)
}

Z_test <- as.data.frame(Z_test)
rownames(Z_test) <- test_ind











