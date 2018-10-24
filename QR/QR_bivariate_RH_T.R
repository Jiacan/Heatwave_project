# this script perform conditional quantile regression based on temperature

library(ncdf4)
library(lattice)
library(quantreg)
library(pbs)
library(splines)
library(pracma)

n_files = 35
days_per_year = 365
day_of_year = 1:365
upper_year = 36
max_time = days_per_year * upper_year
path = '/Users/jy519/Documents/statistic_ES/data/NYC/'
files = list.files("/Users/jy519/Documents/statistic_ES/data/NYC/")[1:n_files]

# get spline basis for day of year (x) and number of years (t)
x = rep(day_of_year, n_files*upper_year)
t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)
x.int.basis = as.matrix(pbs(x, df=3))  # basis functions for interation terms
x.main.basis = as.matrix(pbs(x, df=16))
t.basis = ns(t, df=4)

#---- get linear basis function for temperature: seasonal + long term
# read Tmax
tmax_1 = lapply(files, function(f) {
  path_f = paste(path, f, sep="")
  read.table(path_f, header=FALSE, skip=2)[,2]
})
tmax_2 = unlist(tmax_1)
# for the long-term trend:
temp1 = matrix(tmax_2,nrow = 365, ncol = 1260)
tmax_anm = colMeans(temp1)
tmax_lt = c(t(matrix(rep(tmax_anm, days_per_year), ncol=days_per_year)))

tmax_lt.basis = ns(tmax, df=4)

# for the seasonal cycle:
temp2 = matrix(tmax_2, nrow=13140, ncol=35)
tmax_dtr = matrix(0,dim(temp2)[1], dim(temp2)[2])
for (k in 1:n_files){
  tmax_dtr[,k] = detrend(temp2[,k],'linear')+mean(temp2[,k])
}
tmax_dtr1 = matrix(tmax_dtr, nrow=459900)

# fit detrend tmax with sinusoidal 
tmax_fit=matrix(0L,nrow=13140, ncol=35)
for (i in 1:n_files){
  ssp <- spectrum(tmax_dtr[,i])
  per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
  tt = x[1:13140]
  reslm <- lm(tmax_dtr[,1] ~ sin(2*pi/per*tt)+cos(2*pi/per*tt))
  tmax_fit[,i] = fitted(reslm)
}
tmax_fit1 = matrix(tmax_fit, nrow=459900)

tmax_sc.basis = as.matrix(pbs(tmax_fit1, df=10))

X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis+tmax_sc.basis)
dim(X)

q=c( 0.05, 0.25, 0.5, 0.75, 0.95)
qq = c('5th','25th','50th','75th','95th')
# using Frisch Newton algorithm for quantile regression
coef = matrix(0, dim(X)[2], length(q))
gc()

one_pixel = lapply(files, function(f) {
    path_f = paste(path, f, sep="")
    read.table(path_f, header=FALSE, skip=2)[,4]
})

for (i in 1:length(q)) {
  coef[,i] = rq.fit.pfn(X, y=unlist(one_pixel), tau=q[i],max.bad.fixup=20)$coefficients
  gc()
}

save(X,coef,file='Conditional_QR_QR-Tsmooth_NYC.Rdata')

# plot the time evolution of relative humidity
pt1 = matrix(0L, nrow=n_files, ncol=days_per_year)
pt2 = matrix(0L, nrow=n_files, ncol=days_per_year)
for (j in 1:n_files){
  pt1[j,]=data.matrix(one_pixel[[j]])[1:365]
  pt2[j,]=data.matrix(one_pixel[[j]])[(13140-364):13140]
}

z = apply(X %*% coef, 2, function(b) {
  c(matrix(b, ncol=upper_year))
})
# Reshape z
z1 = t(z)
dim(z1)<-c(5,13140,35)
zm = apply(z1, c(1,2), mean)
zm = t(zm)

matplot(t(pt1), type='p',pch=1,col='gray', ylim=c(10,110), ylab='relative humidity (%)',
        xlab = 'Day of the year')
par(new=TRUE)
matplot(zm[1:365,],type = 'l',pch=1:5,col = 1:5, ylim=c(10,110), ylab='relative humidity (%)',
        xlab = 'Day of the year') #plot
legend(0,0.03, legend = qq, col=1:5, pch=1) # optional legend)
title('Conditional QR of daily relative humidity at NYC in 1990')

# plot bivariate QR fit of RH and Tmax
matplot(t(pt2), type='p',pch=1,col='gray', ylim=c(10,110), ylab='relative humidity (%)',
        xlab = 'Day of the year')
par(new=TRUE)
matplot(zm[(13140-365+1):13140,],type='l',pch=1:5,col = 1:5, ylim=c(10,110),ylab='relative humidity (%)',
        xlab = 'Day of the year') #plot

legend(0,0.03, legend = qq, col=1:5, pch=1) # optional legend)
title('Conditional QR of daily relative humidity at NYC in 2080')

# compute the fitted line for Tmax
load('quantile_estimate_NYC.Rdata')
z_tmax = apply(X %*% coef, 2, function(b) {
  c(matrix(b, ncol=upper_year))
})
# Reshape z
z1_tmax = t(z_tmax)
dim(z1_tmax)<-c(5,13140,35)
zm_tmax = apply(z1_tmax, c(1,2), mean)
zm_tmax = t(zm_tmax)

# plot joint distribution of SH and Tmax
tmax_em = rowMeans(temp2,dim=1)
SH_1 = matrix(unlist(one_pixel),nrow=13140, ncol=35)
SH_em = rowMeans(SH_1, dim=1)
plot(temp2[1:365,],SH_1[1:365,], xlab='Tmax(K)',ylab='relative humidity (%)')
for (h in 1:5){
  lines(zm_tmax[1:365,h], zm[1:365,h],type='l', col=h+1)
}
legend(260,0.02, legend = qq, col=2:6, pch=1)
title('Bivariate QR of RH v.s. Tmax at NYC in 1990')

plot(temp2[(13140-364):13140,],SH_1[(13140-364):13140,], xlab='Tmax(K)',ylab='Specific Humidity (g/Kg)')
for (h in 1:5){
  lines(zm_tmax[(13140-364):13140,h], zm[(13140-364):13140,h],type='l', col=h+1)
}
legend(270,0.02, legend = qq, col=2:6, pch=1)
title('Bivariate QR of RH v.s. Tmax at NYC in 2080')
