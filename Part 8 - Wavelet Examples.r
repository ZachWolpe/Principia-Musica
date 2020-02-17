 rm(list = ls(all = TRUE))

# Create two functions:
# Stationary:
m1 = function(s, alpha0, alpha, beta,f,sig = 0)
{
   res = alpha0/2 +s*0
   for(i in 1:length(alpha))
   {
     res = res + alpha[i]*cos(2*pi*f[i]*s) + beta[i]*sin(2*pi*f[i]*s)
   }
   # Add some noise:
   res = res+rnorm(length(s),0,sig)
   return(res)
}

# Non - Stationary:
m2 = function(s, alpha0, alpha, beta,f,lower,upper,sig = 0)
{
   res = alpha0/2+s*0
   for(i in 1:length(alpha))
   {
     res = res + (alpha[i]*cos(2*pi*f[i]*s) + beta[i]*sin(2*pi*f[i]*s))*(s>=lower[i])*(s<upper[i])
   }
   # Add some noise:
   res = res+rnorm(length(s),0,sig)
   return(res)
}

s = seq(0,1,1/1000)

y1 = m1(s,0,c(1,1,1,1),c(0,0,0,0),c(10,25,50,100),sig = 1)
y2 = m2(s,0,c(1,1,1,1),c(0,0,0,0),c(10,25,50,100),c(0.8,0.6,0.3,0),c(1,0.8,0.6,0.3),sig =1)

par(mfrow = c(1,1))
plot(y1~s,type= 'l')
plot(y2~s,type= 'l')

res1 = fft(y1,inverse = F)
plot(Re(res1),xlim = c(0,200),type= 'h')
res2 = fft(y2,inverse = F)
plot(Re(res2),xlim = c(0,200),type= 'h')
# Note that we can actually reverse the process to recover the signal used to fit:
freq_to_signal = function(coefs,s)
{
  N   = length(s)
  i   = complex(real = 0, imaginary = 1)
  y   = rep(0,N)
  k   = 1:length(coefs)-1
  for(j in 1:N)
  {
    # Reconstruct series at each time point based on ALL the coefficients.
    y[j] <- sum(coefs*exp(i*2*pi*k*(j/N)))/N
  }
  return(y)
}

yhat1  = freq_to_signal(res1,s)
yhat1. = freq_to_signal(res1*(abs(res1)>200),s)
plot(y1~s,type ='l')
#lines(yhat1~s,col = 'red')
lines(yhat1.~s,lwd = 3, col = 'blue')

yhat2 = freq_to_signal(res2,s)
yhat2. = freq_to_signal(res2*(abs(res1)>200),s)
plot(y2~s,type ='l',col = 'lightgray')
points(yhat2~s,col = 'red')
lines(yhat2.~s,lwd = 3, col = 'blue')
lines(m2(s,0,c(1,1,1,1),c(0,0,0,0),c(10,25,50,100),c(0.8,0.6,0.3,0),c(1,0.8,0.6,0.3),sig =1)~s)

# Wavelet analysis offers a more compelling approach, in that transient behaviour
# can be represented efficiently in the time-freqeuncy domain
library('WaveletComp')
res3 = analyze.wavelet(data.frame(y1))
res4 = analyze.wavelet(data.frame(y2))
#W1 = dwt(y1, n.levels=3)

yhat3 = reconstruct(res3)
yhat4 = reconstruct(res4)

print(yhat3)
#plot(yhat3$y1)
#plot(yhat4$axis.1~yhat4$axis.2)


par(mfrow = c(1,1))
wt.image(res3, color.key="interval", legend.params=list(lab="wavelet power levels"))
wt.image(res4, color.key="interval", legend.params=list(lab="wavelet power levels"))


# Non - Stationary:
m3 = function(s, alpha0, alpha, beta,f,lower,upper,sig = 0)
{
   res = alpha0/2+s*0
   for(i in 1:length(alpha))
   {
     res = res + (alpha[i]*cos(2*pi*f[i]*sqrt(s)) + beta[i]*sin(2*pi*f[i]*s))
   }
   res = res+rnorm(length(s),0,sig)
   return(res)
}

y3 = m3(s,0,c(1,1),c(0,0),c(10,20),sig = 0.5)
plot(y3~s,type= 'l')
res5 = analyze.wavelet(data.frame(y3))

wt.image(res5, color.key="interval", legend.params=list(lab="wavelet power levels"))

yhat5 = reconstruct(res5)
lines(yhat5)

## load

data <- read.csv("2.csv", header = F)
month.load <- data[,7]

plot(month.load,type ='l')

## wavelet
resmonth <- analyze.wavelet(data.frame(month.load))

wt.image(resmonth, color.key="interval", legend.params=list(lab="wavelet power levels"))

yhat5 = reconstruct(resmonth)

##fourier

res1 = fft(month.load,inverse = F)
x = 1:length(month.load)
plot(Re(res1),xlim = c(0,50),type= 'h')

yhat2 = freq_to_signal(res1,x)
yhat2. = freq_to_signal(res1*(abs(res1)>1000),x)
plot(month.load~x,type ='l',col = 'lightgray')
lines(yhat2~x,col = 'red')
lines(yhat2.~x, col = 'blue')






