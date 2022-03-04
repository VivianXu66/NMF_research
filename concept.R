# Set the true W1, H1, W2 and H2 to get the original target matrices O and E
set.seed(6)
i <- 10
j <- 20
r <- 5
W1_T <- matrix(rlnorm(i*r), nrow = i, ncol = r, byrow = T)
H1_T <- matrix(rlnorm(r*j), nrow = r, ncol = j, byrow = T)
## Here, in order to simulate the common features of O and E, I create a common part for W1 and W2
## Also, in order to create a simple identity matrix as the coupling matrix, the rows of O and E have the same dimension
g <- 10
h <- 25
W2_T <- matrix(rlnorm(g*r), nrow = g, ncol = r, byrow = T)
W2_T[ ,1:2] <- W1_T[ ,1:2]
H2_T <- matrix(rlnorm(r*h), nrow = r, ncol = h, byrow = T)
# Get O and E
#O <- W1_T%*%H1_T + matrix(rnorm(i*j, 0, 0.001), nrow = i, ncol = j, byrow = T)
#E <- W2_T%*%H2_T + matrix(rnorm(g*h, 0, 0.001), nrow = g, ncol = h, byrow = T)
O <- W1_T%*%H1_T
E <- W2_T%*%H2_T
# Initialize W and H randomly
## Generate all entries of W1 and H2 randomly from the lognormal distribution
W1 <- matrix(rlnorm(i*r), nrow = i, ncol = r, byrow = T)
H1 <- matrix(rlnorm(r*j), nrow = r, ncol = j, byrow = T)

## Generate W2 and H2 with the same method as above
W2 <- matrix(rlnorm(g*r), nrow = g, ncol = r, byrow = T)
H2 <- matrix(rlnorm(r*h), nrow = r, ncol = h, byrow = T)
#w1 <- W1
#w2 <- W2
#h1 <- H1
#h2 <- H2
#O <- w1%*%h1
#E <- w2%*%h2

# Set identity matrix as an example of the coupling matrix  
A <- diag(x = 1, nrow = g, ncol = i, names = T)

# Use ALS algorithm to update W1, W2, H1 and H2 for 10 times.
for(n in 0:10){
    W1t <- O%*%t(H1)%*%solve(H1%*%t(H1))
    W1 <- ifelse(W1t<0, 0, W1t)
    H1t <- solve(t(W1)%*%W1)%*%t(W1)%*%O
    H1 <- ifelse(H1t<0, 0, H1t)
    W2t <- E%*%t(H2)%*%solve(H2%*%t(H2))
    W2 <- ifelse(W2t<0, 0, W2t)
    H2 <- ifelse((solve(t(W2)%*%W2)%*%t(W2)%*%E)<0, 0, solve(t(W2)%*%W2)%*%t(W2)%*%E)
}

# Set up the objective function
##OF <- sum(rowSums((O-W1%*%H1)^2))/2 + lambda1*sum(rowSums((E-W2%*%H2)^2))/2 - lambda2*sum(diag(t(W2)%*%A%*%W1)) + miu*(sum(rowSums(W1)^2)+sum(rowSums(W2)^2))
lambda1 <- 1
lambda2 <- 0.01
miu <- sum(rowSums((O-W1%*%H1)^2))/(sum(rowSums(W1)^2)+sum(rowSums(W2)^2))
#W1<-W1_T
#W2<-W2_T
# Use MU algorithm to update W1, H1, W2 and H2 until convergence
v <- matrix(data = NA, nrow = 10000, ncol =1, byrow = T)
f1 <- matrix(data = NA, nrow = 10000, ncol = 1, byrow = T)
f2 <- matrix(data = NA, nrow = 10000, ncol = 1, byrow = T)
f3 <- matrix(data = NA, nrow = 10000, ncol = 1, byrow = T)
f4 <- matrix(data = NA, nrow = 10000, ncol = 1, byrow = T)
for(m in 0:10000){
    W1 <- W1*((O%*%t(H1)+lambda2/2*t(A)%*%W2)/(W1%*%H1%*%t(H1)+2*miu*W1))
    W2 <- W2*((E%*%t(H2)+lambda2/(2*lambda1)*t(A)%*%W1)/(W2%*%H2%*%t(H2)+2*miu*W2))
    H1 <- H1*((t(W1)%*%O)/(t(W1)%*%W1%*%H1))
    H2 <- H2*((t(W2)%*%E)/(t(W2)%*%W2%*%H2))
    v[m] <- sum(rowSums((O-W1%*%H1)^2))/2 + lambda1*sum(rowSums((E-W2%*%H2)^2))/2 - lambda2*sum(diag(t(W2)%*%A%*%W1)) + miu*(sum(rowSums(W1)^2)+sum(rowSums(W2)^2))
    f1[m] <- sqrt(sum(rowSums((W1-W1_T)^2)))
    f2[m] <- sqrt(sum(rowSums((W2-W2_T)^2)))
    f3[m] <- sqrt(sum(rowSums((H1-H1_T)^2)))
    f4[m] <- sqrt(sum(rowSums((H2-H2_T)^2)))
}
par(mfrow=c(2,3))
plot(c(1:10000), v)
plot(c(1:10000), f1)
plot(c(1:10000), f2)
plot(c(1:10000), f3)
plot(c(1:10000), f4)


## Anchorfree Algorithm

vol_test <- vol_preprocess(t(O), pfactor = 10)
vol.anchor_test <- AnchorFree(vol_test, n.comp = 5)
str(vol.anchor_test$C)
W1 <- vol.anchor_test$C * vol_test$col.factors
apply(cor(W1_T, W1), 1, max)
H1 <- infer_intensities(W1, O)

vol_test <- vol_preprocess(t(E), pfactor = 10)
vol.anchor_test <- AnchorFree(vol_test, n.comp = 5)
str(vol.anchor_test$C)
W2 <- vol.anchor_test$C * vol_test$col.factors
apply(cor(W2_T, W2), 1, max)
H2 <- infer_intensities(W2, E)
