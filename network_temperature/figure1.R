library(MASS)

## Figure 1 for Network Temperature paper 

######################################

## Distribution curves
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

## Network structures 

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
colours = 'Hollywood'

#stronger
qgraph(adj_matrix, layout = "spring", title = "Original Network", edge.labels = FALSE,
       maximum = max(adj_matrix), minimum = 0, theme=colours)

# weaker
qgraph(scaled_adj_matrix, layout = "spring", title = "Scaled Network", edge.labels = FALSE,
       maximum = max(adj_matrix), minimum = 0, theme=colours)

