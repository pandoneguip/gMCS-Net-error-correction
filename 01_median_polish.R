# ==============================================================
# FUNCIONES AUXILIARES
# ==============================================================
setwd("G:/Mi unidad/PFG/gMCS-Net_Angel/code/R")

soft_threshold <- function(z, lambda) {
  sign(z) * pmax(abs(z) - lambda, 0)
}

soft_threshold_nn <- function(z, lambda) {
  pmax(z - lambda, 0)
}

# ==============================================================
# APROXIMACIÓN SPARSE DE RANGO 1
# ==============================================================

sparse_rank1 <- function(A, frac_u = 0.3, frac_v = 0.3,
                         max_iter = 1000, tol = 1e-8, verbose = FALSE) {
  m <- nrow(A)
  n <- ncol(A)
  
  svd1 <- svd(A, nu = 1, nv = 1)
  u <- svd1$u[, 1] * svd1$d[1]
  v <- abs(svd1$v[, 1])
  
  v_norm <- sum(v)
  if (v_norm > 1e-15) {
    v <- v / v_norm
    u <- u * v_norm
  }
  
  obj_prev <- 0.5 * sum((A - u %*% t(v))^2)
  
  for (iter in seq_len(max_iter)) {
    
    Av <- as.vector(A %*% v)
    vtv <- sum(v^2)
    if (vtv < 1e-15) break
    
    lmax_u <- max(abs(Av))
    lambda_u <- frac_u * lmax_u
    u <- soft_threshold(Av, lambda_u) / vtv
    
    utu <- sum(u^2)
    if (utu < 1e-15) break
    
    Atu <- as.vector(crossprod(A, u))
    lmax_v <- max(pmax(Atu, 0))
    lambda_v <- frac_v * lmax_v
    v <- soft_threshold_nn(Atu, lambda_v) / utu
    
    v_norm <- sum(v)
    if (v_norm > 1e-15) {
      u <- u * v_norm
      v <- v / v_norm
    } else {
      break
    }
    
    obj <- 0.5 * sum((A - u %*% t(v))^2)
    
    if (verbose && iter %% 100 == 0) {
      cat(sprintf("  Iter %4d | obj=%.4e | nnz_u=%d/%d | nnz_v=%d/%d\n",
                  iter, obj, sum(u != 0), m, sum(v != 0), n))
    }
    
    if (abs(obj_prev - obj) / (abs(obj_prev) + 1e-15) < tol) break
    obj_prev <- obj
  }
  
  list(u = u, v = v,
       frobenius = sqrt(sum((A - u %*% t(v))^2)),
       iterations = iter,
       nnz_u = sum(u != 0), nnz_v = sum(v != 0),
       u_sparsity = mean(u == 0), v_sparsity = mean(v == 0))
}

# ==============================================================
# CRITERIOS DE SELECCIÓN
# ==============================================================

compute_criteria <- function(A, fit) {
  m <- nrow(A)
  n <- ncol(A)
  mn <- m * n
  
  rss <- fit$frobenius^2
  tss <- sum(A^2)
  df  <- fit$nnz_u + fit$nnz_v
  
  n_eff <- m + n
  bic_corrected <- n_eff * log(rss / n_eff) + log(n_eff) * df
  
  gamma <- 0.5
  comb_u <- if (fit$nnz_u > 0 && fit$nnz_u < m) lchoose(m, fit$nnz_u) else 0
  comb_v <- if (fit$nnz_v > 0 && fit$nnz_v < n) lchoose(n, fit$nnz_v) else 0
  
  ebic <- n_eff * log(rss / n_eff) + log(n_eff) * df +
    2 * gamma * (comb_u + comb_v)
  
  aic_corrected <- n_eff * log(rss / n_eff) + 2 * df
  ratio_pen     <- rss / tss + 0.01 * df / n_eff
  var_explained <- 1 - rss / tss
  
  list(bic_corrected = bic_corrected, ebic = ebic,
       aic_corrected = aic_corrected, ratio_pen = ratio_pen,
       var_explained = var_explained, rss = rss, tss = tss, df = df)
}

# ==============================================================
# BÚSQUEDA EN REJILLA
# ==============================================================

select_lambda <- function(A,
                          frac_u_seq = seq(0.01, 0.95, length.out = 20),
                          frac_v_seq = seq(0.01, 0.95, length.out = 20),
                          criterion  = c("ebic", "bic_corrected",
                                         "aic_corrected", "ratio_pen"),
                          verbose    = TRUE) {
  
  criterion <- match.arg(criterion)
  grid <- expand.grid(frac_u = frac_u_seq, frac_v = frac_v_seq)
  grid$value <- NA; grid$nnz_u <- NA; grid$nnz_v <- NA
  grid$frob  <- NA; grid$var_expl <- NA
  
  for (k in seq_len(nrow(grid))) {
    fit <- tryCatch(
      suppressWarnings(sparse_rank1(A, frac_u = grid$frac_u[k],
                                    frac_v = grid$frac_v[k], max_iter = 500)),
      error = function(e) NULL
    )
    if (is.null(fit) || fit$nnz_u == 0 || fit$nnz_v == 0) next
    
    crit <- compute_criteria(A, fit)
    grid$value[k]    <- crit[[criterion]]
    grid$nnz_u[k]    <- fit$nnz_u
    grid$nnz_v[k]    <- fit$nnz_v
    grid$frob[k]     <- fit$frobenius
    grid$var_expl[k] <- crit$var_explained
  }
  
  grid <- grid[!is.na(grid$value), ]
  best <- grid[which.min(grid$value), ]
  
  if (verbose) {
    cat(sprintf("Criterio: %s\n", criterion))
    cat(sprintf("Grid: %d combinaciones evaluadas\n", nrow(grid)))
    cat(sprintf("Mejor: frac_u=%.3f, frac_v=%.3f\n", best$frac_u, best$frac_v))
    cat(sprintf("       nnz_u=%d, nnz_v=%d, ||err||_F=%.4f\n",
                best$nnz_u, best$nnz_v, best$frob))
    cat(sprintf("       Var explicada=%.2f%%\n", 100 * best$var_expl))
  }
  
  list(best_frac_u = best$frac_u, best_frac_v = best$frac_v,
       best_value  = best$value,  best_nnz_u  = best$nnz_u,
       best_nnz_v  = best$nnz_v,  best_frob   = best$frob,
       criterion   = criterion,   grid        = grid)
}

# ==============================================================
# CARGA DE DATOS
# ==============================================================

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(conflicted)
conflicts_prefer(base::intersect)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

fold_1_data_LELO <- read_delim(
  "G:/Mi unidad/PFG/gMCS-Net_Angel/results/fold_2_data_LELO_results.txt",
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)

Expression <- read_delim(
  "G:/Mi unidad/PFG/gMCS-Net_Angel//data/Expression.txt",
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE
)
Expression <- as.matrix(Expression)

CRISPR <- read_delim(
  "G:\\Mi unidad\\PFG\\gMCS-Net_Angel\\data\\CRISPR.txt",
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE
)
CRISPR <- as.matrix(CRISPR)

genes <- read_delim(
  "G:\\Mi unidad\\PFG\\gMCS-Net_Angel\\data\\genes_index.txt",
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE
)

cells <- read_delim(
  "G:\\Mi unidad\\PFG\\gMCS-Net_Angel\\data\\cell_index.txt",
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE
)

rownames(Expression) <- cells$X1
colnames(Expression) <- genes$X1
rownames(CRISPR)     <- cells$X1
colnames(CRISPR)     <- genes$X1

# ==============================================================
# PREPARACIÓN TRAIN / VAL / TEST
# ==============================================================

split_data <- fold_1_data_LELO %>% group_split(group)

create_matrix <- function(df) {
  df %>%
    select(cell_line, gene, pred) %>%
    pivot_wider(names_from = gene, values_from = pred) %>%
    column_to_rownames("cell_line") %>%
    as.matrix()
}

train_mat <- create_matrix(split_data[[2]])
val_mat   <- create_matrix(split_data[[3]])
test_mat  <- create_matrix(split_data[[1]])

train_mat <- train_mat[, -ncol(train_mat)]
val_mat   <- val_mat[,   -ncol(val_mat)]
test_mat  <- test_mat[,  -ncol(test_mat)]

expr_train <- Expression[rownames(train_mat), ]
pred_train <- CRISPR[rownames(train_mat), ]
viab_train <- train_mat

expr_val <- Expression[rownames(val_mat), ]
pred_val <- CRISPR[rownames(val_mat), ]
viab_val <- val_mat

expr_test <- Expression[rownames(test_mat), ]
pred_test <- CRISPR[rownames(test_mat), ]
viab_test <- test_mat

# ==============================================================
# MATRICES DE ERROR
# ==============================================================

error_train <- pred_train - viab_train
error_val   <- pred_val   - viab_val
error_test  <- pred_test  - viab_test

# ==============================================================
# MEDIAN POLISH: ELIMINACIÓN DE SESGOS SISTEMÁTICOS
# ==============================================================
# El median polish descompone la matriz de error en:
#   error = efecto_global + sesgo_cell_line + sesgo_gen + residuo
#
# Ventaja sobre el centrado simple por medias:
#   - usa medianas → más robusto ante outliers
#   - elimina sesgos tanto por gen como por cell line simultáneamente
#   - el residuo tiene mediana ~0 en filas Y columnas a la vez

mp_train <- medpolish(error_train, maxiter = 100,
                      eps = 0.01, trace.iter = FALSE)

# ----------------------------------------------------------
# Verificación: el residuo debe tener medianas ~0
# ----------------------------------------------------------
cat("=== VERIFICACIÓN MEDIAN POLISH ===\n")
cat(sprintf("Efecto global:                  %.4f\n", mp_train$overall))
cat(sprintf("Mediana global del residuo:     %.2e\n",
            median(mp_train$residuals)))
cat(sprintf("Mediana máx por gen:            %.2e\n",
            max(abs(apply(mp_train$residuals, 2, median)))))
cat(sprintf("Mediana máx por cell line:      %.2e\n\n",
            max(abs(apply(mp_train$residuals, 1, median)))))

# ----------------------------------------------------------
# Análisis de sesgos sistemáticos
# ----------------------------------------------------------
# Cell lines con mayor sesgo: siempre se predicen mal
row_effects <- sort(abs(mp_train$row), decreasing = TRUE)
cat("Top 10 cell lines con mayor sesgo sistemático:\n")
print(head(row_effects, 10))

# Genes con mayor sesgo: el modelo siempre falla en ellos
col_effects <- sort(abs(mp_train$col), decreasing = TRUE)
cat("\nTop 10 genes con mayor sesgo sistemático:\n")
print(head(col_effects, 10))

library(ggplot2)
library(patchwork)

p1 <- ggplot(data.frame(x = mp_train$row), aes(x = x)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = "Sesgo sistemático por línea celular",
       x = "Efecto fila (Median Polish)", y = "Frecuencia") +
  theme_minimal(base_size = 11)

p2 <- ggplot(data.frame(x = mp_train$col), aes(x = x)) +
  geom_histogram(bins = 30, fill = "salmon", color = "white") +
  labs(title = "Sesgo sistemático por gen",
       x = "Efecto columna (Median Polish)", y = "Frecuencia") +
  theme_minimal(base_size = 11)

figura <- p1 + p2

ggsave("figura_median_polish_sesgos.png",
       plot   = figura,
       width  = 10,
       height = 4,
       dpi    = 300,
       bg     = "white")
# ----------------------------------------------------------
# El residuo es la matriz que usaremos para la deflación
# ----------------------------------------------------------
error_train_mp <- mp_train$residuals

cat(sprintf("\nEnergía error original:   %.4e\n", sum(error_train^2)))
cat(sprintf("Energía residuo mp:       %.4e\n",  sum(error_train_mp^2)))
cat(sprintf("Energía eliminada por mp: %.2f%%\n\n",
            100 * (1 - sum(error_train_mp^2) / sum(error_train^2))))
# ==============================================================
# PRIMER AJUSTE EXPLORATORIO (sobre residuo mp)
# ==============================================================

gamma <- 10
A <- t(1 / gamma * asinh(gamma * error_train_mp))

fit <- sparse_rank1(A, frac_u = 0.1, frac_v = 0.1, verbose = TRUE)
str(fit)

# Heatmaps
A_approx      <- fit$u %*% t(fit$v)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
breaks_common <- seq(min(A, A_approx), max(A, A_approx), length.out = 101)

p <- pheatmap(A, cluster_rows = TRUE, cluster_cols = TRUE,
              treeheight_row = 0, treeheight_col = 0,
              show_rownames = FALSE, show_colnames = FALSE,
              color = color_palette, breaks = breaks_common,
              silent = TRUE)

row_order <- p$tree_row$order
col_order <- p$tree_col$order

pheatmap(A, cluster_rows = TRUE, cluster_cols = TRUE,
         treeheight_row = 0, treeheight_col = 0,
         show_rownames = FALSE, show_colnames = FALSE,
         color = color_palette, breaks = breaks_common)

pheatmap(A_approx[row_order, col_order],
         cluster_rows = FALSE, cluster_cols = FALSE,
         treeheight_row = 0, treeheight_col = 0,
         show_rownames = FALSE, show_colnames = FALSE,
         color = color_palette, breaks = breaks_common)

# Nombres para análisis exploratorio
u_named <- fit$u;  names(u_named) <- rownames(A)
v_named <- fit$v;  names(v_named) <- colnames(A)

u_pos <- u_named[u_named > 0]
u_neg <- u_named[u_named < 0]

# Barplot genes
n_show  <- 15
plot_df <- rbind(
  data.frame(gene      = names(head(sort(u_neg), n_show)),
             weight    = as.numeric(head(sort(u_neg), n_show)),
             direction = "negative"),
  data.frame(gene      = names(head(sort(u_pos, decreasing = TRUE), n_show)),
             weight    = as.numeric(head(sort(u_pos, decreasing = TRUE), n_show)),
             direction = "positive")
)
plot_df$gene <- factor(plot_df$gene,
                       levels = plot_df$gene[order(plot_df$weight)])

# Sustituir gene por gene_label en el ggplot
p3 <- ggplot(plot_df, aes(x = gene, y = weight, fill = direction)) +
  geom_col() + coord_flip() +
  labs(title = "Genes que dominan el patrón del componente sparse rank-1",
       x = "Gen", y = "Peso en u") +
  theme_minimal(base_size = 12)

# Símbolos génicos
gene_map_neg <- bitr(names(head(sort(u_neg), 20)),
                     fromType = "ENSEMBL", toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
gene_map_neg$weight <- u_neg[gene_map_neg$ENSEMBL]
gene_map_neg <- gene_map_neg[order(gene_map_neg$weight), ]
print(gene_map_neg)

gene_map_pos <- bitr(names(u_pos),
                     fromType = "ENSEMBL", toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
gene_map_pos$weight <- u_pos[gene_map_pos$ENSEMBL]
gene_map_pos <- gene_map_pos[order(-gene_map_pos$weight), ]
print(gene_map_pos)

# Top 30 cell lines

metadata      <- read.csv("Model.csv", stringsAsFactors = FALSE)
v_nonzero     <- sort(v_named[v_named > 0], decreasing = TRUE)
top_cells_df  <- data.frame(ModelID = names(head(v_nonzero, 30)),
                            weight  = as.numeric(head(v_nonzero, 30)))
annotated_cells <- merge(top_cells_df, metadata, by = "ModelID", all.x = TRUE)
annotated_cells <- annotated_cells[order(-annotated_cells$weight), ]
print(annotated_cells[, c("ModelID", "CellLineName",
                          "OncotreeLineage", "OncotreePrimaryDisease",
                          "weight")])
table(annotated_cells$OncotreeLineage)
# ==============================================================
# PANEL A — Genes con mayor peso (negativos y positivos)
# ==============================================================
gene_map_neg_plot <- head(gene_map_neg, 15)
gene_map_neg_plot$direction <- "negative"

gene_map_pos_plot <- head(gene_map_pos, 15)
gene_map_pos_plot$direction <- "positive"

plot_df <- rbind(gene_map_neg_plot, gene_map_pos_plot)
plot_df$SYMBOL <- factor(plot_df$SYMBOL,
                         levels = plot_df$SYMBOL[order(plot_df$weight)])

p_a <- ggplot(plot_df, aes(x = SYMBOL, y = weight, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("negative" = "#1a5276",
                               "positive" = "#f0a070")) +
  labs(title = "Genes con mayor peso en el componente 1",
       x = NULL, y = "Peso en u") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        axis.text.y     = element_text(size = 9),
        plot.title      = element_text(face = "bold", size = 11))

# ==============================================================
# PANEL B — Distribución de linajes
# ==============================================================
tissue_counts <- as.data.frame(table(annotated_cells$OncotreeLineage))
colnames(tissue_counts) <- c("Lineage", "n")
tissue_counts <- tissue_counts[tissue_counts$n > 0, ]
tissue_counts <- tissue_counts[order(-tissue_counts$n), ]
tissue_counts$Lineage <- factor(tissue_counts$Lineage,
                                levels = tissue_counts$Lineage)

p_b <- ggplot(tissue_counts, aes(x = Lineage, y = n)) +
  geom_col(fill = "#1a5276") +
  labs(title = "Distribución de linajes en líneas con mayor peso",
       x = NULL, y = "Número de líneas celulares") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1, size = 9),
        plot.title   = element_text(face = "bold", size = 11)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# ==============================================================
# COMBINAR Y GUARDAR
# ==============================================================
figura <- p_a / p_b +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave("figura_hipotesis_hemato.png",
       plot   = figura,
       width  = 7,
       height = 11,
       dpi    = 300,
       bg     = "white")



# ==============================================================
# DEFLACIÓN RECURSIVA
# ==============================================================

run_recursive_deflation <- function(error_matrix,
                                    n_components    = 20,
                                    gamma           = 2,
                                    frac_u          = 0.3,
                                    frac_v          = 0.3,
                                    use_grid_search = FALSE,
                                    tol_energy      = 0.005,
                                    verbose         = TRUE) {
  
  R        <- t(error_matrix)
  energy_0 <- sum(R^2)
  cat(sprintf("Energía inicial: %.4e\n\n", energy_0))
  
  components <- list()
  
  for (k in seq_len(n_components)) {
    
    cat(sprintf("============ COMPONENTE %d ============\n", k))
    
    A_k      <- (1 / gamma) * asinh(gamma * R)
    energy_k <- sum(A_k^2)
    cat(sprintf("Energía del residuo (escala asinh): %.4e\n", energy_k))
    
    if (energy_k < tol_energy * (1/gamma * asinh(gamma * sqrt(energy_0)))^2) {
      cat("Energía residual muy pequeña, paramos.\n")
      break
    }
    
    if (use_grid_search) {
      sel <- select_lambda(A_k,
                           frac_u_seq = seq(0.3, 0.8, length.out = 8),
                           frac_v_seq = seq(0.3, 0.8, length.out = 8),
                           criterion  = "ebic",
                           verbose    = FALSE)
      fu <- sel$best_frac_u
      fv <- sel$best_frac_v
      
      if (sel$best_nnz_v > 3 * sel$best_nnz_u ||
          sel$best_nnz_u > 3 * sel$best_nnz_v) {
        cat("  [AVISO] Solución desequilibrada, forzando frac simétrico\n")
        frac_simetrico <- max(fu, fv)
        fu <- frac_simetrico
        fv <- frac_simetrico
      }
    } else {
      fu <- frac_u
      fv <- frac_v
    }
    
    fit_k <- sparse_rank1(A_k, frac_u = fu, frac_v = fv,
                          max_iter = 1000, verbose = FALSE)
    
    cat(sprintf("frac_u=%.3f, frac_v=%.3f | nnz_u=%d, nnz_v=%d | ||err||_F=%.4f\n",
                fu, fv, fit_k$nnz_u, fit_k$nnz_v, fit_k$frobenius))
    
    u_k <- fit_k$u;  names(u_k) <- rownames(A_k)
    v_k <- fit_k$v;  names(v_k) <- colnames(A_k)
    
    R <- R - u_k %*% t(v_k)
    
    energy_after  <- sum(R^2)
    pct_explained <- 100 * (1 - energy_after / sum(t(error_matrix)^2))
    
    cat(sprintf("Energía residual (lineal): %.4e  |  Var. acumulada: %.2f%%\n\n",
                energy_after, pct_explained))
    
    components[[k]] <- list(
      k            = k,   u = u_k,       v = v_k,
      frac_u       = fu,  frac_v = fv,
      nnz_u        = fit_k$nnz_u,        nnz_v = fit_k$nnz_v,
      frobenius    = fit_k$frobenius,
      energy_after = energy_after,        var_cumul = pct_explained
    )
    
    if (energy_after / energy_0 < tol_energy) {
      cat("Criterio de parada por energía alcanzado.\n")
      break
    }
  }
  
  cat(sprintf("\nTotal de componentes extraídos: %d\n", length(components)))
  components
}

# --------------------------------------------------------------
# Ejecutar deflación sobre el RESIDUO del median polish
# --------------------------------------------------------------
components_list <- run_recursive_deflation(
  error_matrix    = error_train_mp,   # <-- residuo del median polish
  n_components    = 20,
  gamma           = 2,
  frac_u          = 0.3,
  frac_v          = 0.3,
  use_grid_search = FALSE,
  tol_energy      = 0.005
)

cat("\nnnz_v componente 1:", components_list[[1]]$nnz_v, "\n")


# ==============================================================
# SOLAPAMIENTO: GENES MEDIAN POLISH vs COMPONENTES SPARSE
# ==============================================================
top_bias_genes <- names(head(col_effects, 20))

for (k in seq_along(components_list)) {
  comp         <- components_list[[k]]
  active_genes <- names(comp$u[comp$u != 0])
  overlap      <- intersect(active_genes, top_bias_genes)
  
  if (length(overlap) > 0) {
    cat(sprintf("Componente %d: %d genes de alto sesgo activos -> %s\n",
                k, length(overlap),
                paste(overlap, collapse = ", ")))
  }
}

# Ver qué genes aparecen más frecuentemente en los componentes
all_overlap_genes <- unlist(lapply(seq_along(components_list), function(k) {
  comp         <- components_list[[k]]
  active_genes <- names(comp$u[comp$u != 0])
  intersect(active_genes, top_bias_genes)
}))

freq_table <- sort(table(all_overlap_genes), decreasing = TRUE)
print(freq_table)

# Convertir a símbolo
gene_freq_df <- bitr(names(freq_table),
                     fromType = "ENSEMBL",
                     toType   = "SYMBOL",
                     OrgDb    = org.Hs.eg.db)
gene_freq_df$frecuencia <- as.numeric(freq_table[gene_freq_df$ENSEMBL])
gene_freq_df <- gene_freq_df[order(-gene_freq_df$frecuencia), ]
print(gene_freq_df)

# ==============================================================
# SCREE PLOT
# ==============================================================


scree_df <- do.call(rbind, lapply(components_list, function(comp) {
  data.frame(componente    = comp$k,
             nnz_u         = comp$nnz_u,
             nnz_v         = comp$nnz_v,
             var_acumulada = comp$var_cumul)
}))
scree_df$var_marginal <- c(scree_df$var_acumulada[1],
                           diff(scree_df$var_acumulada))
print(scree_df)

ggplot(scree_df, aes(x = componente, y = var_marginal)) +
  geom_col(fill = "steelblue", alpha = 0.85) +
  geom_line(aes(y = var_acumulada / max(var_acumulada) * max(var_marginal)),
            color = "red", linewidth = 1) +
  scale_y_continuous(
    name = "Varianza marginal explicada (%)",
    sec.axis = sec_axis(
      ~ . * max(scree_df$var_acumulada) / max(scree_df$var_marginal),
      name = "Varianza acumulada (%)"
    )
  ) +
  labs(title = "Scree plot: varianza explicada por componente sparse",
       x = "Componente") +
  theme_minimal(base_size = 12)

# Comparación con techo SVD
svd_full <- svd(t(error_train_mp))
var_svd  <- cumsum(svd_full$d^2) / sum(svd_full$d^2) * 100

cat("=== COMPARACIÓN SPARSE vs SVD ===\n")
comparison_table <- data.frame(
  n_componentes   = c(5, 10, 20),
  varianza_SVD    = c(var_svd[5], var_svd[10], var_svd[20]),
  varianza_sparse = c(scree_df$var_acumulada[5],
                      scree_df$var_acumulada[10],
                      scree_df$var_acumulada[20])
)
comparison_table$pct_recuperado <- round(
  100 * comparison_table$varianza_sparse / comparison_table$varianza_SVD, 1)
print(comparison_table)

# ==============================================================
# ANÁLISIS DE MEDIO DE CULTIVO
# ==============================================================

analyze_growth_pattern <- function(comp, A, metadata) {
  
  v_full        <- comp$v
  names(v_full) <- colnames(A)
  
  df <- data.frame(ModelID = names(v_full),
                   weight  = as.numeric(v_full))
  
  cols_needed <- c("ModelID", "CellLineName", "GrowthPattern",
                   "OncotreeLineage", "OncotreePrimaryDisease")
  df <- merge(df, metadata[, cols_needed], by = "ModelID", all.x = TRUE)
  
  df$in_component <- df$weight > 0
  active          <- df[df$in_component, ]
  
  cat(sprintf("\n=== COMPONENTE %d: DISTRIBUCIÓN GrowthPattern ===\n", comp$k))
  
  prop_comp  <- prop.table(table(active$GrowthPattern)) * 100
  prop_total <- prop.table(table(df$GrowthPattern))     * 100
  
  comp_table <- data.frame(
    GrowthPattern  = names(prop_comp),
    pct_componente = as.numeric(prop_comp),
    pct_universo   = as.numeric(prop_total[names(prop_comp)]),
    n_componente   = as.numeric(table(active$GrowthPattern)),
    n_total        = as.numeric(table(df$GrowthPattern)[names(prop_comp)])
  )
  comp_table$enriquecimiento <- comp_table$pct_componente / comp_table$pct_universo
  print(comp_table[order(-comp_table$enriquecimiento), ])
  
  p <- ggplot(active[!is.na(active$GrowthPattern), ],
              aes(x = reorder(GrowthPattern, -weight, median),
                  y = weight, fill = GrowthPattern)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.8) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    labs(title = sprintf("Componente %d: peso según tipo de cultivo", comp$k),
         x = "Growth Pattern", y = "Peso en v") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  print(p)
  
  cat("\n=== TEST DE FISHER: Suspensión vs. Resto ===\n")
  df$is_suspension <- grepl("Suspension", df$GrowthPattern, ignore.case = TRUE)
  contingency <- table(en_componente = df$in_component,
                       suspension    = df$is_suspension)
  print(contingency)
  fisher_res <- fisher.test(contingency)
  cat(sprintf("Odds Ratio: %.2f  |  p-value: %.4f\n",
              fisher_res$estimate, fisher_res$p.value))
  
  invisible(df)
}

# Ejecutar para componente 1 y componente 20
df_growth_1  <- analyze_growth_pattern(components_list[[1]],  A, metadata)
df_growth_20 <- analyze_growth_pattern(components_list[[20]], A, metadata)

# ==============================================================
# GO ENRICHMENT POR COMPONENTE
# ==============================================================

run_go_enrichment <- function(u_named, all_genes,
                              top_n_pos    = 50,
                              top_n_neg    = 50,
                              ont          = "BP",
                              pvalueCutoff = 0.05) {
  
  u_pos    <- sort(u_named[u_named > 0], decreasing = TRUE)
  u_neg    <- sort(u_named[u_named < 0], decreasing = FALSE)
  top_pos  <- names(u_pos)[seq_len(min(top_n_pos, length(u_pos)))]
  top_neg  <- names(u_neg)[seq_len(min(top_n_neg, length(u_neg)))]
  genes_oi <- unique(c(top_pos, top_neg))
  
  if (length(genes_oi) < 5) return(NULL)
  
  ego_all <- tryCatch(
    enrichGO(gene = genes_oi, universe = all_genes,
             OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
             ont = ont, pAdjustMethod = "BH",
             pvalueCutoff = pvalueCutoff, qvalueCutoff = 0.2,
             readable = TRUE),
    error = function(e) NULL
  )
  
  ego_pos <- NULL
  if (length(top_pos) >= 5)
    ego_pos <- tryCatch(
      enrichGO(gene = top_pos, universe = all_genes,
               OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
               ont = ont, pAdjustMethod = "BH",
               pvalueCutoff = pvalueCutoff, qvalueCutoff = 0.2,
               readable = TRUE),
      error = function(e) NULL)
  
  ego_neg <- NULL
  if (length(top_neg) >= 5)
    ego_neg <- tryCatch(
      enrichGO(gene = top_neg, universe = all_genes,
               OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
               ont = ont, pAdjustMethod = "BH",
               pvalueCutoff = pvalueCutoff, qvalueCutoff = 0.2,
               readable = TRUE),
      error = function(e) NULL)
  
  list(all = ego_all, pos = ego_pos, neg = ego_neg,
       genes_pos = top_pos, genes_neg = top_neg)
}

all_genes_universe <- names(components_list[[1]]$u)

go_by_component <- lapply(seq_along(components_list), function(k) {
  comp   <- components_list[[k]]
  cat(sprintf("GO enrichment para componente %d...\n", k))
  result <- run_go_enrichment(comp$u, all_genes_universe,
                              top_n_pos = 50, top_n_neg = 50,
                              ont = "BP")
  if (!is.null(result$all) && nrow(result$all) > 0) {
    cat(sprintf("  -> %d términos significativos\n", nrow(result$all)))
    cat(sprintf("  -> Top 3: %s\n",
                paste(head(result$all$Description, 3), collapse = "; ")))
  } else {
    cat("  -> Sin enriquecimiento significativo\n")
  }
  result
})

names(go_by_component) <- paste0("comp_", seq_along(go_by_component))



# Dotplots componentes representativos
for (k in c(1, 2, 3, 4)) {
  ego <- go_by_component[[k]]$all
  if (!is.null(ego) && nrow(ego) > 0) {
    p <- dotplot(ego, showCategory = 10,
                 title = sprintf("GO BP - Componente %d", k)) +
      theme(axis.text.y = element_text(size = 8))
    print(p)
  }
}


# ==============================================================
# TABLA RESUMEN FINAL
# ==============================================================

summary_table <- do.call(rbind, lapply(components_list, function(comp) {
  v       <- comp$v
  top_ids <- names(sort(v[v > 0], decreasing = TRUE)[1:min(20, sum(v > 0))])
  ann     <- metadata[metadata$ModelID %in% top_ids, ]
  
  top_tissue <- if (nrow(ann) > 0)
    names(sort(table(ann$OncotreeLineage), decreasing = TRUE))[1] else NA
  top_growth <- if (nrow(ann) > 0)
    names(sort(table(ann$GrowthPattern),  decreasing = TRUE))[1] else NA
  
  ego    <- go_by_component[[comp$k]]$all
  top_go <- if (!is.null(ego) && nrow(ego) > 0) ego$Description[1] else "—"
  
  data.frame(comp = comp$k, nnz_genes = comp$nnz_u, nnz_cells = comp$nnz_v,
             var_acum_pct = round(comp$var_cumul, 2),
             tejido_dominante = top_tissue, cultivo = top_growth,
             top_GO = top_go)
}))



cat("\n========================================\n")
cat("TABLA RESUMEN DE COMPONENTES\n")
cat("========================================\n")
print(summary_table)

# ==============================================================
# EXPRESIÓN DIFERENCIAL GENES COMPONENTE 1
# ==============================================================

u_1        <- components_list[[1]]$u
top_genes_1 <- names(sort(abs(u_1[u_1 != 0]),
                          decreasing = TRUE)[1:20])

suspension_cells <- metadata$ModelID[grepl("Suspension",
                                           metadata$GrowthPattern)]
adherent_cells   <- metadata$ModelID[grepl("Adherent",
                                           metadata$GrowthPattern)]

expr_susp <- Expression[rownames(Expression) %in% suspension_cells,
                        top_genes_1]
expr_adh  <- Expression[rownames(Expression) %in% adherent_cells,
                        top_genes_1]

pvals     <- sapply(top_genes_1, function(g)
  t.test(expr_susp[, g], expr_adh[, g])$p.value)
pvals_adj <- p.adjust(pvals, method = "BH")

results_df <- data.frame(ensembl  = names(pvals),
                         pval_raw = pvals,
                         pval_adj = pvals_adj)

gene_symbols <- bitr(names(pvals), fromType = "ENSEMBL",
                     toType = "SYMBOL", OrgDb = org.Hs.eg.db)
results_df <- merge(results_df, gene_symbols,
                    by.x = "ensembl", by.y = "ENSEMBL")

results_df$mean_susp <- sapply(results_df$ensembl, function(g)
  if (g %in% colnames(expr_susp)) mean(expr_susp[, g]) else NA)
results_df$mean_adh  <- sapply(results_df$ensembl, function(g)
  if (g %in% colnames(expr_adh))  mean(expr_adh[, g])  else NA)
results_df$log2FC    <- log2(results_df$mean_susp / results_df$mean_adh)

results_df <- results_df[order(results_df$pval_adj), ]
print(results_df[, c("SYMBOL", "pval_raw", "pval_adj", "log2FC")])


# ==============================================================
# ANÁLISIS DE SUBTIPOS HEMATOLÓGICOS EN EL COMPONENTE 1
# ==============================================================
# El componente 1 captura un error sistemático asociado a células
# en suspensión (OR=11.71). Aquí desagregamos dentro de los linajes
# hematológicos (Lymphoid/Myeloid) para identificar qué subtipos
# oncológicos concretos concentran el mayor peso en el componente.

u_1_deflacion   <- components_list[[1]]$u
u_neg_comp1     <- u_1_deflacion[u_1_deflacion < 0]
genes_neg_comp1 <- names(u_neg_comp1)

v_k <- components_list[[1]]$v
df  <- data.frame(ModelID = names(v_k), weight = as.numeric(v_k))
df  <- merge(df, metadata[, c("ModelID", "OncotreeLineage",
                              "OncotreeSubtype", "OncotreePrimaryDisease",
                              "PatientSubtypeFeatures")],
             by = "ModelID", all.x = TRUE)
active <- df[df$weight > 0, ]
hemato <- active[active$OncotreeLineage %in% c("Lymphoid", "Myeloid"), ]

cat(sprintf("=== COMPONENTE 1 (n hemato activos = %d) ===\n", nrow(hemato)))

subtipo_summary <- hemato %>%
  group_by(OncotreeSubtype) %>%
  summarise(n = n(), peso_medio = mean(weight), peso_max = max(weight)) %>%
  filter(n >= 2) %>%
  arrange(desc(peso_medio))
print(subtipo_summary)

cat("Error medio suspensión:", mean(error_train[rownames(error_train) %in% suspension_cells, genes_neg_comp1]), "\n")
cat("Error medio adherentes:", mean(error_train[rownames(error_train) %in% adherent_cells,  genes_neg_comp1]), "\n")
# Error general por tipo de cultivo
cat("=== ERROR MEDIO POR GROWTH PATTERN ===\n")
cells_by_growth <- split(metadata$ModelID, metadata$GrowthPattern)

for (pattern in names(cells_by_growth)) {
  cells <- cells_by_growth[[pattern]]
  cells_in_train <- cells[cells %in% rownames(error_train)]
  if (length(cells_in_train) == 0) next
  
  err <- mean(error_train[cells_in_train, ])
  cat(sprintf("%s (n=%d): %.4f\n", pattern, length(cells_in_train), err))
}
# Error general por linaje
cat("\n=== ERROR MEDIO POR ONCOTREE LINEAGE ===\n")
cells_by_lineage <- split(metadata$ModelID, metadata$OncotreeLineage)

for (lineage in names(cells_by_lineage)) {
  cells <- cells_by_lineage[[lineage]]
  cells_in_train <- cells[cells %in% rownames(error_train)]
  if (length(cells_in_train) < 3) next
  
  err <- mean(error_train[cells_in_train, ])
  cat(sprintf("%-30s (n=%d): %.4f\n", lineage, length(cells_in_train), err))
}





#####
library(ggplot2)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)

# Vectores de pesos
u_named <- components_list[[1]]$u
names(u_named) <- rownames(A)

# Limpiar versiones ENSG si las hay
names(u_named) <- sub("\\..*", "", names(u_named))

# Separar positivos y negativos antes de mapear
u_neg <- sort(u_named[u_named < 0])
u_pos <- sort(u_named[u_named > 0], decreasing = TRUE)

# Mapear negativos
gene_map_neg <- bitr(names(head(u_neg, 20)),
                     fromType = "ENSEMBL", toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
gene_map_neg$weight <- u_neg[gene_map_neg$ENSEMBL]
gene_map_neg <- gene_map_neg[order(gene_map_neg$weight), ]
gene_map_neg$direction <- "negative"

# Mapear positivos
gene_map_pos <- bitr(names(head(u_pos, 20)),
                     fromType = "ENSEMBL", toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
gene_map_pos$weight <- u_pos[gene_map_pos$ENSEMBL]
gene_map_pos <- gene_map_pos[order(-gene_map_pos$weight), ]
gene_map_pos$direction <- "positive"

# Tomar top 15 de cada uno
plot_df <- rbind(head(gene_map_neg, 15), head(gene_map_pos, 15))
plot_df$SYMBOL <- factor(plot_df$SYMBOL, levels = plot_df$SYMBOL[order(plot_df$weight)])

# PANEL A
p_a <- ggplot(plot_df, aes(x = SYMBOL, y = weight, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("negative" = "#1a5276", "positive" = "#f0a070")) +
  labs(title = "Genes con mayor peso en el componente 1",
       x = NULL, y = "Peso en u") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 11, face = "bold"))

# PANEL B
v_named <- components_list[[1]]$v
names(v_named) <- colnames(A)
v_nonzero <- sort(v_named[v_named > 0], decreasing = TRUE)
top_cells_df <- data.frame(ModelID = names(head(v_nonzero, 30)),
                           weight = as.numeric(head(v_nonzero, 30)))
annotated_cells <- merge(top_cells_df, metadata, by = "ModelID", all.x = TRUE)

tissue_counts <- as.data.frame(table(annotated_cells$OncotreeLineage))
colnames(tissue_counts) <- c("Lineage", "n")
tissue_counts <- tissue_counts[tissue_counts$n > 0, ]
tissue_counts <- tissue_counts[order(-tissue_counts$n), ]
tissue_counts$Lineage <- factor(tissue_counts$Lineage, levels = tissue_counts$Lineage)

p_b <- ggplot(tissue_counts, aes(x = Lineage, y = n)) +
  geom_col(fill = "#1a5276") +
  labs(title = "Distribución de tejidos en líneas con mayor peso",
       x = NULL, y = "Número de líneas celulares") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        plot.title = element_text(size = 11, face = "bold")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

# COMBINAR
p_a / p_b + plot_annotation(tag_levels = "a")

