library(MASS)

######################################
## make distributions
# unimodal
x <- seq(-5, 5, length=100)
y <- dnorm(x)
plot(x,y, type = "l", lwd = 2, axes = FALSE)
  
#bimodal
x <- seq(-5, 5, length = 100)
y1 <- dnorm(x, mean = -2, sd = 1)
y2 <- dnorm(x, mean = 2, sd = 1)
y <- y1 + y2
plot(x, y, type = "l", lwd = 2, axes = FALSE)

######################################
## make networks
# Number of nodes
n <- 10

# Create covariance matrices for the two datasets
cov1 <- diag(n) + 0.8  # Set diagonal elements to high value
cov1[lower.tri(cov1)] <- 0.8  # Set off-diagonal elements to high value

# For the second dataset, set most off-diagonal elements to low values to create fewer and weaker connections
cov2 <- diag(n) + 0.8  # Set diagonal elements to high value
cov2[lower.tri(cov2)] <- 0.3  # Set off-diagonal elements to low value

# Simulate data for the two datasets
data1 <- mvrnorm(n = 100, mu = rep(0, n), Sigma = cov1)
data2 <- mvrnorm(n = 100, mu = rep(0, n), Sigma = cov2)

# Convert the data to binary format (-1 and 1)
data1_binary <- ifelse(data1 > 0, -1, 1)
data2_binary <- ifelse(data2 > 0, -1, 1)

# Display the first few rows of the simulated datasets
head(data1_binary)
head(data2_binary)

# give col names
labels <- c(1:10)
colnames(data1_binary) <- labels
colnames(data2_binary) <- labels

# network colours: can be "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
# plot networks
figure.model.1 <- Ising(data1_binary, vars=labels)
figure.network.1 <- getmatrix(figure.model.1, "omega")
qgraph(figure.network.1, layout = 'circle', theme = 'classic')

figure.model.2 <- Ising(data2_binary, vars=labels)
figure.network.2 <- getmatrix(figure.model.2, "omega")
qgraph(figure.network.2, layout = 'circle', theme = 'classic')


####

n <- 10

# Create covariance matrices for the two datasets
cov1 <- diag(n) + 0.8  # Set diagonal elements to high value
cov1[lower.tri(cov1)] <- 0.8  # Set off-diagonal elements to high value

# Simulate data for the two datasets
data1 <- mvrnorm(n = 100, mu = rep(0, n), Sigma = cov1)

# Convert the data to binary format (-1 and 1)
data1_binary <- ifelse(data1 > 0, -1, 1)

# give col names
labels <- c(1:10)
colnames(data1_binary) <- labels

# network colours: can be "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".
# plot networks
figure.model.1 <- Ising(data1_binary, vars=labels)
figure.network.1 <- getmatrix(figure.model.1, "omega")
qgraph(figure.network.1, layout = 'circle', theme = 'classic')


###

# Load necessary libraries
library(qgraph)
library(IsingFit)

# Number of nodes
n <- 10

# Create covariance matrices for the dataset
cov1 <- diag(n) + 0.8  # Set diagonal elements to high value
cov1[lower.tri(cov1)] <- 0.8  # Set off-diagonal elements to high value

# Simulate data for the dataset
data1 <- mvrnorm(n = 100, mu = rep(0, n), Sigma = cov1)

# Convert the data to binary format (-1 and 1)
data1_binary <- ifelse(data1 > 0, 0, 1)

# Give column names
labels <- c(1:10)
colnames(data1_binary) <- labels

# Fit Ising model to the binary data
model <- IsingFit(data1_binary, plot = FALSE)
adj_matrix <- model$weiadj

# Multiply the edge weights
multiplier <- 2  # You can adjust this multiplier as needed
adj_matrix_multiplied <- adj_matrix * multiplier

# Plot the network with multiplied edge weights
qgraph(adj_matrix_multiplied, layout = 'circle', theme = 'classic')


#########

set.seed(123)

# adjacency matrix
adj_matrix <- matrix(runif(100, min = 0, max = 1), nrow = 10, ncol = 10)

# Make the matrix symmetric to represent an undirected network
adj_matrix <- (adj_matrix + t(adj_matrix)) / 2

# diagonal to zero (no self-loops)
diag(adj_matrix) <- 0

# scaled version of the adjacency matrix
scaled_adj_matrix <- adj_matrix * 0.5

## plot
par(mfrow = c(1, 2))

# network colours: can be "classic","colorblind","gray","Hollywood","Borkulo", "gimme","TeamFortress","Reddit","Leuven"or"Fried".

colours = 'classic'

#stronger
qgraph(adj_matrix, layout = "spring", title = "Original Network", edge.labels = FALSE,
       maximum = max(adj_matrix), minimum = 0, theme=colours)

# weaker
qgraph(scaled_adj_matrix, layout = "spring", title = "Scaled Network", edge.labels = FALSE,
       maximum = max(adj_matrix), minimum = 0, theme=colours)

