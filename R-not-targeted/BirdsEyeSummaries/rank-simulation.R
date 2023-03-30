
# What does a null rank distribution look like?
N <- 1000

xx <- runif(N)
yy <- runif(N)

rank.xx <- rank(xx)
rank.yy <- rank(yy)

hist((rank.xx + rank.yy) / 2)

#END
