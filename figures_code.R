# Code to reproduce figures in Miller and Allesina (2022)
# Habitat heterogeneity, environmental feedbacks, and species coexistence across timescales
# zachmiller@uchicago.edu

library("tidyverse")
library("reshape2")
source("functions.R")


##### FIG 2 #####

set.seed(16)

n <- 5
m <- 1

# sample diagonal elements from a distribution with larger values (to ensure species are habitat specialists)
diag.dist <- function(n) {
  return(runif(n, m, n * m))
}

P <- sample_P(n, check.feas = TRUE, check.stab = FALSE,
         m = 1.2, # sample with larger m to ensure feasibility over a wider range of w for clearer visualization
         distribution = runif, diag.distribution = diag.dist)$P

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

p_vec <- seq(0, 10, by = 0.005) # a range of p values  
out <- tibble(p = c(), sp = c(), x = c(), w = c(), entropy = c(), feasible = c()) # initialize empty tibble to store values
for(p in p_vec){
  
  # sample over a gradient of w values from very even (low p) to very uneven (high p) distribution of patch types
  
  w <- v^p / sum(v^p)
  x <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% w - m * rep(1, n)) # equilibrium values for x_i (Eq. 3)
  
  out <- out %>% add_row(p = rep(p, n),
                         sp = 1:n,
                         x = as.vector(x),
                         w = as.vector(w),
                         entropy = -sum(w * log(w)), # Shannon entropy of the landscape
                         feasible = ifelse(all(x > 0), rep(TRUE, n), rep(FALSE, n))
                         )  
}

# save 300 x 350
out %>%
  ggplot() + 
  aes(x = entropy, y = w, group = sp, color = as.factor(sp)) + 
  geom_line(size = 1) + 
  annotate(geom = "rect",
            xmin = min(out %>% filter(feasible) %>% pull(entropy)),
            xmax = max(out %>% filter(feasible) %>% pull(entropy)),
            ymin = -Inf, ymax = 1,
            color = NA,
            fill = "blue",
            alpha = 0.25) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  #xlim(c(0.5, log(n))) + 
  scale_x_continuous(trans = "exp") + 
  xlab("Landscape entropy") + ylab("Frequency (w)") + 
  theme_classic() + 
  guides(color = "none")

# save 300 x 350
out %>% filter(feasible) %>%
  ggplot() + 
  aes(x = entropy, y = x, group = sp, color = as.factor(sp)) + 
  geom_line(size = 1) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  scale_x_continuous(trans = "exp") + 
  xlab("Landscape entropy") + ylab("Frequency (x)") + 
  theme_classic() + 
  guides(color = "none")


##### FIG 3 #####

set.seed(6)

n <- 5
m <- 1

# sample diagonal elements from a distribution with larger values (to ensure species are habitat specialists)
diag.dist <- function(n) {
  return(runif(n, m, 3 * m))
}

P <- sample_P(n, check.feas = TRUE, check.stab = FALSE, m = m,
              distribution = runif, diag.distribution = diag.dist)$P

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

# save 300 x 300
tidy_P <- melt(P)
colnames(tidy_P) <- c("row", "col", "val")
tidy_P %>% ggplot() + 
  aes(x = row, y = col, fill = val) + 
  geom_tile() + 
  scale_y_reverse() + 
  scale_fill_viridis_c() + 
  guides(fill = "none") +
  theme_void()

# save 300 x 300
x <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% v - m * rep(1, n))
y <- m * rowSums(solve(P))
xy <- as.vector(x) %*% t(as.vector(y))
tidy_xy <- melt(xy)
colnames(tidy_xy) <- c("row", "col", "val")
tidy_xy %>% ggplot() + 
  aes(x = row, y = col, fill = val) + 
  geom_tile() + 
  scale_y_reverse() + 
  scale_fill_viridis_c() + 
  guides(fill = "none") +
  theme_void()

# save 300 x 300
X <- P * xy
tidy_X <- melt(X)
colnames(tidy_X) <- c("row", "col", "val")
tidy_X %>% ggplot() + 
  aes(x = row, y = col, fill = val) + 
  geom_tile() + 
  scale_y_reverse() + 
  scale_fill_viridis_c() + 
  guides(fill = "none") +
  theme_void()


##### FIG 4 #####

set.seed(10)

n <- 3
m <- 0.1

P <- sample_P(n, check.feas = TRUE, check.stab = TRUE,
              symmetric = TRUE, m = 0.5)$P

d <- 0.00001

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

# sample initial w near v to ensure feasiblity of the x(w(0))
x_of_w <- -1 
while(any(x_of_w < 0)) {
  w_init <- v * runif(n, 0.5, 1.5)
  w_init <- (w_init / sum(w_init)) # normalize
  x_of_w <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% w_init - m * rep(1, n))
}

# random initial x and y
y_init <- 0.5 * w_init
x_init <- runif(n)
x_init <- (x_init / sum(x_init)) * (1 - sum(y_init))

times <- seq(0, 200000, 5)
pars <- list(n = n, P = P, m = m, d = d)
init <- c(x_init, y_init, w_init)

# Integrate the general model (Eq. 5)
out <- ode(y = init, times = times, parms = pars, 
           func = general_feedback_model, method = "ode45")

# Define slow system for w
A <- solve(diag(rowSums(solve(P))) %*% t(P)) - diag(n)
b <- -m * rowSums(solve(t(P)) - solve(P))

# Integrate slow system (Eq. 7)
out_slow <- ode(y = w_init, times = times, 
                parms = list(n = n, P = P, m = m, d = d, A = A, b = b),
                func = slow_system, method = "ode45")

# Collect w values from full model and slow system
w_predicted <- cbind(1:nrow(out), out_slow[, 2:(n + 1)])
colnames(w_predicted) <- c("time", 1:n)
tidy_w_predicted <- pivot_longer(as_tibble(w_predicted), -time)
colnames(tidy_w_predicted) <- c("time", "species", "predicted")

w_actual <- cbind(1:nrow(out), out[, (2 * n + 2):(3 * n + 1)])
colnames(w_actual) <- c("time", 1:n)
tidy_w_actual <- pivot_longer(as_tibble(w_actual), -time)
colnames(tidy_w_actual) <- c("time", "species", "actual")

dat_w <- inner_join(tidy_w_actual, tidy_w_predicted) %>% mutate(type = "Patch states (w)")

# Collect x values from full model and predicted from slow system
w <- w_predicted[, -1]
x_predicted <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% t(w) - m)

x_predicted <- cbind(1:nrow(out), t(predicted))
colnames(x_predicted) <- c("time", 1:n)
tidy_x_predicted <- pivot_longer(as_tibble(x_predicted), -time)
colnames(tidy_x_predicted) <- c("time", "species", "predicted")

x_actual <- cbind(1:nrow(out), out[, 2:(n + 1)])
colnames(x_actual) <- c("time", 1:n)
tidy_x_actual <- pivot_longer(as_tibble(x_actual), -time)
colnames(tidy_x_actual) <- c("time", "species", "actual")

dat_x <- inner_join(tidy_x_actual, tidy_x_predicted) %>% mutate(type = "Species (x)")

dat <- dat_w %>% add_row(dat_x)
dat$type = factor(dat$type, levels=c('Species (x)', 'Patch states (w)'))

# save 400 x 350
scale_round <- function(x) sprintf("%.2f", x)
dat %>%
  ggplot() +
  aes(x = time) + 
  geom_line(aes(y = actual, group = species, color = species), size = 1, alpha = 0.9) +
  geom_line(aes(y = predicted, group = species), lty = "dashed") +
  xlab("Time") + 
  ylab("Frequency") + 
  scale_x_log10() + 
  scale_y_continuous(labels = scale_round) + 
  guides(color = "none") + 
  facet_grid(type~., scales = "free") + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white"))


##### FIG 5 #####

### Panel A

set.seed(18)

n <- 3
m <- 0.5

P <- sample_P(n, check.feas = TRUE, check.stab = TRUE, symmetric = TRUE, m = 0.5)$P

d <- 0.00001

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

w_init <- v * runif(n, 0.5, 1.5)
w_init <- (w_init / sum(w_init))
y_init <- 0.8 * w_init
x_init <- runif(n)
x_init <- (x_init / sum(x_init)) * (1 - sum(y_init))

times <- seq(0, 1000000, 20)
pars <- list(n = n, P = P, m = m, d = d)
init <- c(x_init, y_init, w_init)

out <- ode(y = init, times = times, parms = pars, 
           func = general_feedback_model, method = "ode45")

dat <- cbind(1:nrow(out), out[, 2:(n + 1)])
colnames(dat) <- c("time", 1:n)
dat <- pivot_longer(as_tibble(dat), -time)
colnames(dat) <- c("time", "species", "value")

# save 350 x 250
scale_round <- function(x) sprintf("%.2f", x)
dat %>% #filter(time > 1) %>%
  ggplot() +
  aes(x = time) + 
  geom_line(aes(y = value, group = species, color = species), size = 1, alpha = 0.9) +
  xlab("Time") + 
  ylab("Frequency") + 
  scale_y_continuous(labels = scale_round) + 
  guides(color = "none") + 
  theme_bw()

### Panel B

set.seed(30)

n <- 3
m <- 0.5

P <- sample_P(n, check.feas = TRUE, check.stab = FALSE, symmetric = TRUE, m = 0.5)$P
P <- P + diag(n) # increase diagonal to make coexistence unstable

d <- 0.00001

perron <- Re(eigen(diag(rowSums(solve(P))) %*% t(P))$vectors[, 1]) # Perron eigenvector of D(P^{-1} 1) P^T
v <- m * rowSums(solve(P)) + (1 - m * sum(solve(P))) * (perron / sum(perron)) # choose v according to Eq. S28 to ensure feasibility

w_init <- v * runif(n, 0.98, 1.02)
w_init <- (w_init / sum(w_init))

# start dynamics on slow manifold to avoid distracting transients
x_init <- solve(t(P)) %*% (diag(1 / rowSums(solve(P))) %*% w_init - m * rep(1, n))
y_init <- m %*% rowSums(solve(P))
y_init <- (y_init / sum(y_init)) * (1 - sum(x_init))

# integrate for some time before perturbation
times <- seq(0, 20000, 5)
pars <- list(n = n, P = P, m = m, d = d)
init <- c(x_init, y_init, w_init)

out1 <- ode(y = init, times = times, parms = pars, 
            func = general_feedback_model, method = "ode45")

perturb_idx <- 2 # perturb the species that will go extinct
current <- out1[nrow(out1), -1] # take the current state of the system
y <- current[(n+1):(2*n)]
x <- current[perturb_idx] 
perturb_dist <- y * P[perturb_idx, ] # distribution of patches occupied by the perturbed species
current[perturb_idx] <- x * 0.05 # reduce frequency by 95%
current[(n+1):(2*n)] <- y + x * (1 - 0.05) * (perturb_dist / sum(perturb_dist)) # add newly emptied patches to y 

# integrate for longer following perturbation
out2 <- ode(y = current, times = seq(0, 50000, 5), parms = pars, 
            func = general_feedback_model, method = "ode45")

out <- rbind(out1[, -1], out2[, -1]) # combine time series

# reformat for plotting
dat <- cbind(1:nrow(out), out[, 1:n])
colnames(dat) <- c("time", 1:n)
dat <- pivot_longer(as_tibble(dat), -time)
colnames(dat) <- c("time", "species", "value")

# save 350 x 250
scale_round <- function(x) sprintf("%.2f", x)
dat %>% filter(time > 1) %>%
  ggplot() +
  aes(x = time) + 
  geom_line(aes(y = value, group = species, color = species), size = 1, alpha = 0.9) +
  geom_vline(xintercept = length(times), lty = "dashed") + 
  xlab("Time") + 
  ylab("Frequency") + 
  scale_y_continuous(labels = scale_round) + 
  guides(color = "none") + 
  theme_bw()

##### FIG S1 #####

plot_circle <- function(center, radius){
  
  angles <- seq(0, 2 * pi, 0.0001)
  
  out <- tibble(x = radius * cos(angles) + center[1],
                y = radius * sin(angles) + center[2])
  return(out)
}

outer_circ <- plot_circle(c(0,0), 1)
outer_circ %>% ggplot() + aes(x = x, y = y) + 
  geom_vline(xintercept = 0, alpha = 0.5) + 
  geom_hline(yintercept = 0, alpha = 0.5) + 
  geom_polygon(alpha = 0.2, size = 1, fill = "blue", color = "black") + 
  geom_polygon(data = plot_circle(c(1/2, 0), 1/2), 
               alpha = 0.4, size = 1, fill = "red", color = "black") + 
  geom_point(x = 1, y = 0, pch = 4, size = 3) + 
  xlab("Re") + 
  ylab("Im") +
  coord_equal() + 
  theme_bw()