library(vrnmf)
library(ggplot2)
##Set the true W1, H1, W2 and H2 to get the original target matrices O and E
set.seed(7)
v1 <- 10
D <- 200
f1 <- 5
v2 <- 10
f2 <- 5
m1 <- matrix(data = NA, nrow = 6, ncol = 20)
m2 <- matrix(data = NA, nrow = 6, ncol = 20)
m3 <- matrix(data = NA, nrow = 6, ncol = 20)
m4 <- matrix(data = NA, nrow = 6, ncol = 20)

for (n in 1:20){

C1_T <- matrix(rlnorm(v1*f1), nrow = v1, ncol = f1, byrow = T)
W1_T <- matrix(rlnorm(f1*D), nrow = f1, ncol = D, byrow = T)
## Here, in order to simulate the common features of O and E, I create a common part for W1 and W2
## Also, in order to create a simple identity matrix as the coupling matrix, the rows of O and E have the same dimension
C2_T <- matrix(rlnorm(v2*f2), nrow = v2, ncol = f2, byrow = T)
W2_T <- matrix(rlnorm(f2*D), nrow = f2, ncol = D, byrow = T)
W2_T[1:2,] <- W1_T[1:2,]

# Get O and E
#O <- C1_T%*%W1_T + matrix(rnorm(v1*D, 0, 0.001), nrow = v1, ncol = D, byrow = T)
#E <- C2_T%*%W2_T + matrix(rnorm(v2*D, 0, 0.001), nrow = v2, ncol = D, byrow = T)
O <- C1_T%*%W1_T
E <- C2_T%*%W2_T

for(i in 1:6){
  vol_test_1 <- vol_preprocess(t(O))
  vol.anchor_test_1 <- AnchorFree(vol_test_1,n.comp=(i+4))
  str(vol.anchor_test_1$C)
  C1 <-  vol.anchor_test_1$C * vol_test_1$col.factors
  W1 <- infer_intensities(C1, O)
  m1[i,n] <- max(round(cor(t(W1), t(W1_T)),2)[,1])
}
for(i in 1:6){
  vol_test_1 <- vol_preprocess(t(O))
  vol.anchor_test_1 <- AnchorFree(vol_test_1,n.comp=(i+4))
  str(vol.anchor_test_1$C)
  C1 <-  vol.anchor_test_1$C * vol_test_1$col.factors
  W1 <- infer_intensities(C1, O)
  m2[i,n] <- max(round(cor(t(W1), t(W1_T)),2)[,2])
}
for(i in 1:6){
  vol_test_2 <- vol_preprocess(t(E))
  vol.anchor_test_2 <- AnchorFree(vol_test_2,n.comp=(i+4))
  str(vol.anchor_test_1$C)
  C2 <-  vol.anchor_test_2$C * vol_test_2$col.factors
  W2 <- infer_intensities(C2, E)
  m3[i,n] <- max(round(cor(t(W2), t(W2_T)),2)[,1])
}
for(i in 1:6){
  vol_test_2 <- vol_preprocess(t(E))
  vol.anchor_test_2 <- AnchorFree(vol_test_2,n.comp=(i+4))
  str(vol.anchor_test_2$C)
  C2 <-  vol.anchor_test_2$C * vol_test_2$col.factors
  W2 <- infer_intensities(C2, E)
  m4[i,n] <- max(round(cor(t(W2), t(W2_T)),2)[,2])
}
}

data1 <- data.frame(factors = rep(c(1:6), 20), values = as.vector(m1), simulations = rep(c(1:20), rep(6, 20)))
data1 %>% ggplot( aes(x=factors, y=values, group=simulations, color=factor(simulations))) +
  geom_point() +
  geom_boxplot(aes(group = factors)) +
  scale_color_viridis_d() +
  theme(
    legend.position="none",
    plot.title = element_text(size=14)
  ) +
  ggtitle("Comparison between W1 and the 1st col of true W1")
  

data2 <- data.frame(factors = rep(c(1:6), 20), values = as.vector(m2), simulations = rep(c(1:20), rep(6, 20)))
data2 %>% ggplot( aes(x=factors, y=values, group=simulations, color=factor(simulations))) +
    geom_line() +
    scale_color_viridis_d() +
  theme(
    legend.position="none",
    plot.title = element_text(size=14)
  ) +
    ggtitle("Comparison between W1 and the 2nd col of true W1")

data3 <- data.frame(factors = rep(c(1:6), 20), values = as.vector(m3), simulations = rep(c(1:20), rep(6, 20)))
data3 %>% ggplot( aes(x=factors, y=values, group=simulations, color=factor(simulations))) +
    geom_line() +
    scale_color_viridis_d() +
  theme(
    legend.position="none",
    plot.title = element_text(size=14)
  ) +
    ggtitle("Comparison between W2 and the 1st col of true W2")
  
data4 <- data.frame(factors = rep(c(1:6), 20), values = as.vector(m4), simulations = rep(c(1:20), rep(6, 20)))
data4 %>% ggplot( aes(x=factors, y=values, group=simulations, color=factor(simulations))) +
    geom_line() +
    scale_color_viridis_d() +
  theme(
    legend.position="none",
    plot.title = element_text(size=14)
  ) +
    ggtitle("Comparison between W2 and the 2nd col of true W2")
  
  

vol_test_1 <- vol_preprocess(t(O))
vol.anchor_test_1 <- AnchorFree(vol_test_1,n.comp=5)
str(vol.anchor_test_1$C)
C1 <-  vol.anchor_test_1$C * vol_test_1$col.factors
W1 <- infer_intensities(C1, O)
r1[n] <- max(round(cor(t(W1), t(W1_T)),2)[,1])


vol_test_1 <- vol_preprocess(t(O))
vol.anchor_test_1 <- AnchorFree(vol_test_1,n.comp=5)
str(vol.anchor_test_1$C)
C1 <-  vol.anchor_test_1$C * vol_test_1$col.factors
W1 <- infer_intensities(C1, O)
r2[n] <- max(round(cor(t(W1), t(W1_T)),2)[,2])


vol_test_2 <- vol_preprocess(t(E))
vol.anchor_test_2 <- AnchorFree(vol_test_2,n.comp=5)
str(vol.anchor_test_2$C)
C2 <-  vol.anchor_test_2$C * vol_test_2$col.factors
W2 <- infer_intensities(C2, E)
r3[n] <- max(round(cor(t(W2), t(W2_T)),2)[,1])


vol_test_2 <- vol_preprocess(t(E))
vol.anchor_test_2 <- AnchorFree(vol_test_2,n.comp=5)
str(vol.anchor_test_2$C)
C2 <-  vol.anchor_test_2$C * vol_test_2$col.factors
W2 <- infer_intensities(C2, E)
r4[n] <- max(round(cor(t(W2), t(W2_T)),2)[,2])


X <- rbind(O, E)
q <- c()
for(i in 1:6){
  vol_test_3 <- vol_preprocess(t(X))
  vol.anchor_test_3 <- AnchorFree(vol_test_3,n.comp=(i+4))
  str(vol.anchor_test_3$C)
  R <-  vol.anchor_test_3$C * vol_test_3$col.factors
  U <- infer_intensities(R, X)
  q3[i] <- max(apply(cor(t(U), t(W1_T)), 1, max))
}
plot(c(1:6), q3, ylim = c(0.5,1))


vol_test_3 <- vol_preprocess(t(X))
vol.anchor_test_3 <- AnchorFree(vol_test_3,n.comp=8)
str(vol.anchor_test_3$C)
R <-  vol.anchor_test_3$C * vol_test_3$col.factors
U <- infer_intensities(R, X)

round(cor(t(U), t(W1_T)),2)
round(cor(t(U), t(W2_T)),2)
round(cor(t(W1_T),t( W2_T)),2)

round(cor(R, rbind(C1_T,C2_T)),2)


