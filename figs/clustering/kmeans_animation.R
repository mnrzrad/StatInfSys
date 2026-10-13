# Animated GIF of Lloyd's k-means algorithm, used in the Clustering chapter
# (section "K-means Algorithm").
#
# Regenerate from the project root with:
#   Rscript figs/clustering/kmeans_animation.R
# Requires the packages ggplot2 and gifski.

library(ggplot2)
library(gifski)

# ---- Data: three Gaussian groups in two dimensions -------------------------
set.seed(2026)
X <- rbind(
  cbind(rnorm(50, 0, 0.7), rnorm(50, 0.0, 0.7)),
  cbind(rnorm(50, 4, 0.7), rnorm(50, 0.5, 0.7)),
  cbind(rnorm(50, 2, 0.7), rnorm(50, 3.5, 0.7))
)
colnames(X) <- c("x1", "x2")
K <- 3

# ---- Lloyd's algorithm, recording every step -------------------------------
set.seed(48)                                 # initial centroids: K random observations
C <- X[sample(nrow(X), K), ]

assign_to <- function(P, C) {
  d <- sapply(seq_len(nrow(C)), function(k) rowSums(sweep(P, 2, C[k, ])^2))
  max.col(-d, ties.method = "first")
}
wss <- function(cl, C) sum((X - C[cl, ])^2)

steps <- list(list(type = "init", C = C, C_old = NULL, cl = NULL))
cl <- rep(0, nrow(X))
iter <- 0
repeat {
  new_cl <- assign_to(X, C)
  if (all(new_cl == cl)) break
  cl <- new_cl
  iter <- iter + 1
  steps[[length(steps) + 1]] <- list(type = "assign", iter = iter, C = C, C_old = NULL, cl = cl)
  C_old <- C
  C <- t(sapply(seq_len(K), function(k) colMeans(X[cl == k, , drop = FALSE])))
  steps[[length(steps) + 1]] <- list(type = "update", iter = iter, C = C, C_old = C_old, cl = cl)
}
steps[[length(steps) + 1]] <- list(type = "done", iter = iter, C = C, C_old = NULL, cl = cl)

# ---- One frame per step -----------------------------------------------------
pal  <- c("1" = "#E69F00", "2" = "#56B4E9", "3" = "#009E73")   # colour-blind safe
grid <- expand.grid(x1 = seq(-2.5, 6.5, length.out = 180),
                    x2 = seq(-2.5, 6.0, length.out = 170))
lims <- list(x = c(-2.5, 6.5), y = c(-2.5, 6.0))

make_frame <- function(s) {
  pts  <- data.frame(X, cluster = if (is.null(s$cl)) NA else factor(s$cl, levels = 1:K))
  cent <- data.frame(s$C, cluster = factor(1:K, levels = 1:K))
  colnames(cent)[1:2] <- c("x1", "x2")
  # Nearest-centroid regions for the current centroids
  regions <- cbind(grid, cluster = factor(assign_to(as.matrix(grid), s$C), levels = 1:K))

  title <- switch(s$type,
    init   = "Initialization: K = 3 random observations become the centroids",
    assign = sprintf("Iteration %d, assign step: each point joins its nearest centroid", s$iter),
    update = sprintf("Iteration %d, update step: each centroid moves to the mean of its points", s$iter),
    done   = sprintf("Converged after %d iterations: the assignments no longer change", s$iter))
  subtitle <- if (is.null(s$cl)) "Points not yet assigned" else
    sprintf("Total within-cluster sum of squares = %.1f", wss(s$cl, s$C))

  g <- ggplot() +
    geom_raster(data = regions, aes(x1, x2, fill = cluster), alpha = 0.12) +
    scale_fill_manual(values = pal, guide = "none")

  if (is.null(s$cl)) {
    g <- g + geom_point(data = pts, aes(x1, x2), color = "gray45", size = 2)
  } else {
    g <- g + geom_point(data = pts, aes(x1, x2, color = cluster), size = 2)
  }

  if (!is.null(s$C_old)) {
    moves <- data.frame(s$C_old, s$C)
    colnames(moves) <- c("x_from", "y_from", "x_to", "y_to")
    g <- g +
      geom_point(data = setNames(data.frame(s$C_old), c("x1", "x2")),
                 aes(x1, x2), shape = 4, size = 4, stroke = 1.2, color = "gray40") +
      geom_segment(data = moves, aes(x_from, y_from, xend = x_to, yend = y_to),
                   arrow = arrow(length = unit(0.12, "inches")), linewidth = 0.9)
  }

  g +
    geom_point(data = cent, aes(x1, x2, fill = cluster),
               shape = 23, size = 6, stroke = 1.3, color = "black") +
    scale_color_manual(values = pal, guide = "none") +
    coord_equal(xlim = lims$x, ylim = lims$y, expand = FALSE) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(size = 12.5, face = "bold"),
          panel.grid.minor = element_blank()) +
    labs(title = title, subtitle = subtitle, x = expression(x[1]), y = expression(x[2]))
}

# ---- Write the GIF ----------------------------------------------------------
frames_dir <- file.path(tempdir(), "kmeans_frames")
dir.create(frames_dir, showWarnings = FALSE)
files <- character(0)
for (i in seq_along(steps)) {
  f <- file.path(frames_dir, sprintf("frame_%03d.png", i))
  ggsave(f, make_frame(steps[[i]]), width = 7.5, height = 6.4, dpi = 96, bg = "white")
  # Hold the first and last frames longer
  reps  <- if (i == 1 || i == length(steps)) 3 else 1
  files <- c(files, rep(f, reps))
}

gifski(files, gif_file = "figs/clustering/kmeans_lloyd.gif",
       width = 720, height = 614, delay = 1.2)
cat("Iterations:", iter, " Frames:", length(steps), "\n")
