rm(list = ls())
devtools::clean_dll()
devtools::document()
devtools::load_all()


unlink("NAMESPACE")
unlink("man/plot.p_hclust.Rd")

devtools::clean_dll()
devtools::document()
devtools::load_all()

set.seed(123)

build.data <- function() {
  g1 <- cbind(
    rnorm(50, mean = 0, sd = 0.5),
    rnorm(50, mean = 0, sd = 0.5)
  )
  g2 <- cbind(
    rnorm(50, mean = 3.7, sd = 0.5),
    rnorm(50, mean = 2.8, sd = 0.5)
  )
  g3 <- cbind(
    rnorm(100, mean = 20, sd = 4),
    rnorm(100, mean = 0, sd = 4)
  )
  X <- rbind(g1, g2, g3)
  as.data.frame(X)
}
X <- build.data()
Y <- build.data()
true_labels <- as.integer(
  rep(c(1, 2, 3), times = c(50, 50, 100))
)

ggplot(data = X, aes(x = X[, 1], y = X[, 2], color = factor(true_labels))) +
  geom_point() +
  theme_bw()

dismat <- dist(X, method = "euclidean")^2

hcl <- fastcluster::hclust(
  dismat,
  method = "ward.D"
)


p_hcl <- p.dendrogram(X = X,verbose = FALSE)

class(p_hcl)
names(p_hcl)
X[p_hcl$split_indices}

plot(
  p_hcl,
  plot_config = list(
    log_height = TRUE,
    labels_pvalues = TRUE,
    strip_height_factor = 0.0005
  ),
  groups = list()
)


# bien prendre en compte
get.partion
labels <- extract_labels(p_hcl) # get.partition
length(labels)
