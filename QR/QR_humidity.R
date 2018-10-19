# this script perform conditional quantile regression based on temperature

library(ncdf4)
library(lattice)
library(quantreg)
library(pbs)
library(splines)

n_files = 35
days_per_year = 365
day_of_year = 1:365
upper_year = 36
max_time = days_per_year * upper_year
path = '/Users/jy519/Documents/statistic_ES/data/NYC/'
files = list.files("/Users/jy519/Documents/statistic_ES/data/NYC/")[1:n_files]

# get spline bases for day of year (x) and number of years (t)
x = rep(day_of_year, n_files*upper_year)
t = rep(c(t(matrix(rep(1:upper_year, days_per_year), ncol=days_per_year))), n_files)
x.int.basis = as.matrix(pbs(x, df=3))  # basis functions for interation terms
x.main.basis = as.matrix(pbs(x, df=16))
t.basis = ns(t, df=4)
X = model.matrix(~x.main.basis + t.basis + x.int.basis:t.basis)
dim(X)

q=c( 0.05, 0.25, 0.5, 0.75, 0.95)
qq = c('5th','25th','50th','75th','95th')
# using Frisch Newton algorithm for quantile regression
coef = matrix(0, dim(X)[2], length(q))
rm(t,x, t.basis, x.main.basis, x.int.basis)
gc()

one_pixel = lapply(files, function(f) {
    path_f = paste(path, f, sep="")
    read.table(path_f, header=FALSE, skip=2)[,3]
})

for (i in 1:length(q)) {
  coef[,i] = rq.fit.pfn(X, y=unlist(one_pixel), tau=q[i],max.bad.fixup=20)$coefficients
  gc()
}

# save(X,coef,file='quantile_estimate_NYC.Rdata')

# plot quantile coefficient in the first year

pt1 = matrix(0L, nrow=n_files, ncol=days_per_year)
pt2 = matrix(0L, nrow=n_files, ncol=days_per_year)
for (j in 1:n_files){
  pt1[j,]=data.matrix(one_pixel[[j]])[1:365]
  pt2[j,]=data.matrix(one_pixel[[j]])[(max_time-365+1):max_time]
}

z = apply(X %*% coef, 2, function(b) {
  c(matrix(b, ncol=upper_year))
})
matplot(t(pt1), type='p',pch=1,col='gray', ylim=c(0,0.03), ylab='specific humidity',
        xlab = 'Day of the year')
par(new=TRUE)
matplot(z[1:365,],type = 'l',pch=1:5,col = 1:5, ylim=c(0,0.03), ylab='specific humidity',
        xlab = 'Day of the year') #plot
legend(0,0.03, legend = qq, col=1:5, pch=1) # optional legend)
title('quantile fit of daily specific humidity at NYC in 1990')

# plot quantile coefficient in the last year


matplot(t(pt2), type='p',pch=1,col='gray', ylim=c(0,0.03), ylab='specific humidity',
        xlab = 'Day of the year')
par(new=TRUE)
matplot(z[(13140-365+1):13140,],type='l',pch=1:5,col = 1:5, ylim=c(0,0.03),ylab='specific humidity',
        xlab = 'Day of the year') #plot

legend(0,0.03, legend = qq, col=1:5, pch=1) # optional legend)
title('quantile fit of daily specific humidity at NYC in 2080')
