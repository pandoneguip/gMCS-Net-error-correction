# ==============================================================
# PIPELINE K-FOLD ANIDADO DENTRO DEL FOLD 2
# ==============================================================
# Partición de las 658 cell lines de train en 5 grupos
# Esquema: 3 train + 1 val + 1 test rotando
#
# Iter 1: train=1,2,3  val=4  test=5
# Iter 2: train=5,1,2  val=3  test=4
# Iter 3: train=4,5,1  val=2  test=3
# Iter 4: train=3,4,5  val=1  test=2
# Iter 5: train=2,3,4  val=5  test=1
#
# ESQUEMA ESTRICTO POR ITERACION:
#   - TEST no se toca hasta el final de cada iteracion
#   - Todos los hiperparametros se eligen en VAL
#   - TEST se evalua una unica vez al final de cada iteracion
#
# PREREQUISITO: train_mat, pred_train, Expression en memoria
# ==============================================================

library(dplyr)
library(tidyr)
library(tibble)
library(conflicted)

conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::union)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

# ==============================================================
# FUNCIONES
# ==============================================================

pinv <- function(A, k = NULL) {
  m <- nrow(A); n <- ncol(A)
  if (is.null(k)) k <- min(m, n)
  k <- min(k, min(m, n))
  svd_A <- svd(A, nu = k, nv = k)
  U <- svd_A$u[, 1:k, drop = FALSE]
  d <- svd_A$d[1:k]
  V <- svd_A$v[, 1:k, drop = FALSE]
  tol <- max(m, n) * .Machine$double.eps * max(d)
  d_inv <- ifelse(d > tol, 1/d, 0)
  V %*% diag(d_inv, nrow = k) %*% t(U)
}

cor_por_celula <- function(mat_real, mat_pred, label = "") {
  cors <- diag(cor(t(mat_real), t(mat_pred)))
  med <- median(cors, na.rm = TRUE)
  mn  <- mean(cors,   na.rm = TRUE)
  if (label != "")
    cat(sprintf("%-40s media=%.4f  mediana=%.4f\n", label, mn, med))
  invisible(cors)
}

cor_por_gen <- function(mat_real, mat_pred, label = "") {
  cors <- diag(cor(mat_real, mat_pred))
  med <- median(cors, na.rm = TRUE)
  mn  <- mean(cors,   na.rm = TRUE)
  if (label != "")
    cat(sprintf("%-40s media=%.4f  mediana=%.4f\n", label, mn, med))
  invisible(cors)
}

# ==============================================================
# PARTICIÓN EN 5 GRUPOS
# ==============================================================

set.seed(42)

cell_lines_train <- rownames(pred_train)
n <- length(cell_lines_train)

grupos <- sample(rep(1:5, length.out = n))
names(grupos) <- cell_lines_train

cat("Distribución de células por grupo:\n")
print(table(grupos))

esquema <- data.frame(
  iter     = 1:5,
  i_test   = c(5, 4, 3, 2, 1),
  i_val    = c(4, 3, 2, 1, 5),
  i_train1 = c(1, 5, 4, 3, 2),
  i_train2 = c(2, 1, 5, 4, 3),
  i_train3 = c(3, 2, 1, 5, 4)
)

cat("\nEsquema de rotación:\n")
print(esquema[, c("iter", "i_train1", "i_train2", "i_train3",
                   "i_val", "i_test")])

# ==============================================================
# K-FOLD
# ==============================================================

resultados_b   <- list()
resultados_oof <- list()

for (iter in 1:5) {

  cat(sprintf("\n%s\n", strrep("=", 55)))
  cat(sprintf("ITERACION %d — Train: %d,%d,%d | Val: %d | Test: %d\n",
              iter,
              esquema$i_train1[iter], esquema$i_train2[iter],
              esquema$i_train3[iter],
              esquema$i_val[iter], esquema$i_test[iter]))
  cat(sprintf("%s\n", strrep("=", 55)))

  # Índices de cada conjunto
  cells_test    <- names(grupos)[grupos == esquema$i_test[iter]]
  cells_val     <- names(grupos)[grupos == esquema$i_val[iter]]
  cells_train_k <- names(grupos)[grupos %in% c(esquema$i_train1[iter],
                                                esquema$i_train2[iter],
                                                esquema$i_train3[iter])]

  cat(sprintf("Train: %d | Val: %d | Test: %d cell lines\n",
              length(cells_train_k), length(cells_val),
              length(cells_test)))

  # Matrices
  viab_train_k <- train_mat[cells_train_k, ]
  viab_val_k   <- train_mat[cells_val,     ]
  viab_test_k  <- train_mat[cells_test,    ]

  pred_train_k <- pred_train[cells_train_k, ]
  pred_val_k   <- pred_train[cells_val,     ]
  pred_test_k  <- pred_train[cells_test,    ]

  error_train_k <- viab_train_k - pred_train_k

  cat(sprintf("Error medio train: %.4f\n",
              mean(abs(error_train_k))))

  # ----------------------------------------------------------
  # CAPA 1: MEDIAN POLISH — ELEGIR ALPHA EN VAL
  # ----------------------------------------------------------
  mp_k <- medpolish(error_train_k, maxiter = 100,
                     eps = 0.01, trace.iter = FALSE)

  sesgo_genes_k <- mp_k$col + mp_k$overall

  alphas_mp <- seq(0, 2, length.out = 100)
  err_mp_val <- sapply(alphas_mp, function(a) {
    pred_new <- t(t(pred_val_k) + a * sesgo_genes_k)
    sum((viab_val_k - pred_new)^2)
  })
  alpha_mp_k <- alphas_mp[which.min(err_mp_val)]

  # Aplicar solo a VAL
  pred_val_mp_k <- t(t(pred_val_k) + alpha_mp_k * sesgo_genes_k)

  cat(sprintf("Capa 1 MP: alpha=%.4f | RMSE val: %.4f -> %.4f\n",
              alpha_mp_k,
              sqrt(mean((viab_val_k - pred_val_k)^2)),
              sqrt(mean((viab_val_k - pred_val_mp_k)^2))))

  # ----------------------------------------------------------
  # CAPA 2: PSEUDOINVERSA + ASINH — ELEGIR k, GAMMA, ALPHA EN VAL
  # ----------------------------------------------------------

  svd_check <- svd(pred_train_k, nu = 0, nv = 0)
  cat(sprintf("Num condicion pred_train: %.2e\n",
              max(svd_check$d) /
                min(svd_check$d[svd_check$d > 0])))

  ks          <- c(10, 20, 50, 100, 150, 200, 250, 300,
                   min(nrow(pred_train_k), ncol(pred_train_k)))
  ks          <- unique(ks)
  gammas      <- c(1, 2, 5, 10)
  alphas_pinv <- seq(0, 2, length.out = 100)

  mejor_err    <- Inf
  gamma_k      <- 1
  alpha_pinv_k <- 0
  k_opt        <- 50

  for (k in ks) {
    for (g in gammas) {
      err_train_asinh <- (1/g) * asinh(g * error_train_k)
      A_val_k     <- pred_val_mp_k %*% pinv(pred_train_k, k = k)
      err_val_est <- A_val_k %*% err_train_asinh

      if (max(abs(err_val_est)) > 10) next

      for (a in alphas_pinv) {
        pred_new <- pred_val_mp_k + a * err_val_est
        err <- sum((viab_val_k - pred_new)^2)
        if (err < mejor_err) {
          mejor_err    <- err
          gamma_k      <- g
          alpha_pinv_k <- a
          k_opt        <- k
        }
      }
    }
    cat(sprintf("k=%3d completado\n", k))
  }

  cat(sprintf("Capa 2 Pinv: k=%d gamma=%d alpha=%.4f\n",
              k_opt, gamma_k, alpha_pinv_k))

  # Aplicar solo a VAL
  err_train_asinh_opt_k <- (1/gamma_k) * asinh(gamma_k * error_train_k)
  A_val_k       <- pred_val_mp_k %*% pinv(pred_train_k, k = k_opt)
  err_val_est_k <- A_val_k %*% err_train_asinh_opt_k
  pred_val_pinv_k <- pred_val_mp_k + alpha_pinv_k * err_val_est_k

  cat(sprintf("RMSE val tras pinv: %.4f\n",
              sqrt(mean((viab_val_k - pred_val_pinv_k)^2))))

  # ----------------------------------------------------------
  # CAPA 3: EXPRESIÓN GÉNICA — ELEGIR ALPHA EN VAL
  # ----------------------------------------------------------

  cells_expr_train_k <- base::intersect(cells_train_k, rownames(Expression))
  cells_expr_val_k   <- base::intersect(cells_val,     rownames(Expression))
  cells_expr_test_k  <- base::intersect(cells_test,    rownames(Expression))

  expr_means_k <- colMeans(Expression[cells_expr_train_k, ], na.rm = TRUE)

  corr_expr_k <- sapply(colnames(error_train_k), function(g) {
    if (!g %in% colnames(Expression)) return(NA)
    cor(Expression[cells_expr_train_k, g],
        error_train_k[cells_expr_train_k, g],
        use = "complete.obs")
  })

  genes_signal_k <- names(corr_expr_k)[
    abs(corr_expr_k) > 0.10 & !is.na(corr_expr_k)]

  cat(sprintf("Genes con señal expresion: %d\n", length(genes_signal_k)))

  beta_expr_k <- sapply(colnames(error_train_k), function(g) {
    if (!g %in% colnames(Expression)) return(NA)
    expr_g  <- Expression[cells_expr_train_k, g]
    error_g <- error_train_k[cells_expr_train_k, g]
    ok <- !is.na(expr_g) & !is.na(error_g)
    if (sum(ok) < 20) return(NA)
    cov(expr_g[ok], error_g[ok]) / var(expr_g[ok])
  })

  aplicar_expr_k <- function(pred_mat, cells_expr, alpha_e) {
    pred_corr <- pred_mat
    for (g in genes_signal_k) {
      if (!g %in% colnames(Expression)) next
      if (is.na(beta_expr_k[g])) next
      expr_g <- Expression[cells_expr, g]
      expr_g[is.na(expr_g)] <- expr_means_k[g]
      correccion <- alpha_e * beta_expr_k[g] *
        (expr_g - expr_means_k[g])
      correccion <- pmax(pmin(correccion, 0.1), -0.1)
      pred_corr[cells_expr, g] <- pred_corr[cells_expr, g] + correccion
    }
    pred_corr
  }

  # Seleccionar alpha en VAL minimizando error cuadratico
  alphas_expr  <- seq(0, 5, length.out = 200)
  err_expr_val <- sapply(alphas_expr, function(a) {
    pred_new <- aplicar_expr_k(pred_val_pinv_k, cells_expr_val_k, a)
    sum((viab_val_k - pred_new)^2)
  })
  alpha_expr_k <- alphas_expr[which.min(err_expr_val)]

  cat(sprintf("Capa 3 Expr: alpha=%.4f\n", alpha_expr_k))

  # Aplicar solo a VAL para verificar
  pred_val_final_k <- aplicar_expr_k(pred_val_pinv_k,
                                      cells_expr_val_k, alpha_expr_k)

  cat(sprintf("RMSE val final: %.4f\n",
              sqrt(mean((viab_val_k - pred_val_final_k)^2))))

  # ----------------------------------------------------------
  # HIPERPARAMETROS SELECCIONADOS EN VAL
  # ----------------------------------------------------------
  cat(sprintf("Hiperparametros: alpha_mp=%.4f k=%d gamma=%d alpha_pinv=%.4f alpha_expr=%.4f\n",
              alpha_mp_k, k_opt, gamma_k, alpha_pinv_k, alpha_expr_k))

  # ----------------------------------------------------------
  # EVALUACION FINAL EN TEST — UNA SOLA VEZ
  # ----------------------------------------------------------
  cat("--- Evaluacion TEST ---\n")

  # Aplicar pipeline completo a TEST con parametros de VAL
  pred_test_mp_k    <- t(t(pred_test_k) + alpha_mp_k * sesgo_genes_k)

  A_test_k          <- pred_test_mp_k %*% pinv(pred_train_k, k = k_opt)
  err_test_est_k    <- A_test_k %*% err_train_asinh_opt_k
  pred_test_pinv_k  <- pred_test_mp_k + alpha_pinv_k * err_test_est_k

  pred_test_final_k <- aplicar_expr_k(pred_test_pinv_k,
                                       cells_expr_test_k, alpha_expr_k)

  # Metricas finales
  rmse_base  <- sqrt(mean((viab_test_k - pred_test_k)^2))
  rmse_final <- sqrt(mean((viab_test_k - pred_test_final_k)^2))

  cor_cel_base  <- mean(diag(cor(t(viab_test_k), t(pred_test_k))))
  cor_cel_final <- mean(diag(cor(t(viab_test_k), t(pred_test_final_k))))

  cor_gen_base  <- median(diag(cor(viab_test_k, pred_test_k)))
  cor_gen_final <- median(diag(cor(viab_test_k, pred_test_final_k)))

  cat(sprintf("RMSE:    %.4f -> %.4f (%.1f%%)\n",
              rmse_base, rmse_final,
              100*(rmse_base - rmse_final)/rmse_base))
  cat(sprintf("Cor cel: %.4f -> %.4f\n",
              cor_cel_base, cor_cel_final))
  cat(sprintf("Cor gen: %.4f -> %.4f\n",
              cor_gen_base, cor_gen_final))

  resultados_b[[iter]] <- list(
    iter          = iter,
    rmse_base     = rmse_base,
    rmse_final    = rmse_final,
    cor_cel_base  = cor_cel_base,
    cor_cel_final = cor_cel_final,
    cor_gen_base  = cor_gen_base,
    cor_gen_final = cor_gen_final,
    alpha_mp      = alpha_mp_k,
    k_opt         = k_opt,
    gamma         = gamma_k,
    alpha_pinv    = alpha_pinv_k,
    alpha_expr    = alpha_expr_k,
    n_genes_signal= length(genes_signal_k)
  )

  # Guardar predicciones out-of-fold para analisis posterior
  resultados_oof[[iter]] <- list(
    cells     = cells_test,
    pred_base = pred_test_k,
    pred      = pred_test_final_k,
    viab      = viab_test_k
  )
}

# ==============================================================
# RESUMEN FINAL
# ==============================================================

cat(sprintf("\n%s\n", strrep("=", 55)))
cat("RESUMEN K-FOLD COMPLETO\n")
cat(sprintf("%s\n\n", strrep("=", 55)))

resumen_b <- do.call(rbind, lapply(resultados_b, function(r) {
  data.frame(
    iter          = r$iter,
    rmse_base     = round(r$rmse_base,     4),
    rmse_final    = round(r$rmse_final,    4),
    mejora_rmse   = round(100*(r$rmse_base - r$rmse_final) /
                            r$rmse_base, 1),
    cor_cel_base  = round(r$cor_cel_base,  4),
    cor_cel_final = round(r$cor_cel_final, 4),
    cor_gen_base  = round(r$cor_gen_base,  4),
    cor_gen_final = round(r$cor_gen_final, 4),
    alpha_mp      = round(r$alpha_mp,      4),
    k_opt         = r$k_opt,
    gamma         = r$gamma,
    alpha_pinv    = round(r$alpha_pinv,    4),
    alpha_expr    = round(r$alpha_expr,    4),
    genes_signal  = r$n_genes_signal
  )
}))

print(resumen_b)

cat(sprintf("\nMejora RMSE:   media=%.2f%%  sd=%.2f%%  IC95=[%.2f%%, %.2f%%]\n",
            mean(resumen_b$mejora_rmse),
            sd(resumen_b$mejora_rmse),
            mean(resumen_b$mejora_rmse) -
              1.96 * sd(resumen_b$mejora_rmse) / sqrt(5),
            mean(resumen_b$mejora_rmse) +
              1.96 * sd(resumen_b$mejora_rmse) / sqrt(5)))

cat(sprintf("Cor gen final: media=%.4f  sd=%.4f  IC95=[%.4f, %.4f]\n",
            mean(resumen_b$cor_gen_final),
            sd(resumen_b$cor_gen_final),
            mean(resumen_b$cor_gen_final) -
              1.96 * sd(resumen_b$cor_gen_final) / sqrt(5),
            mean(resumen_b$cor_gen_final) +
              1.96 * sd(resumen_b$cor_gen_final) / sqrt(5)))

cat(sprintf("Cor cel final: media=%.4f  sd=%.4f\n",
            mean(resumen_b$cor_cel_final),
            sd(resumen_b$cor_cel_final)))

cat(sprintf("\nK optimo por iteracion:   %s\n",
            paste(resumen_b$k_opt, collapse = ", ")))
cat(sprintf("Alpha expr por iteracion: %s\n",
            paste(round(resumen_b$alpha_expr, 3), collapse = ", ")))

# ==============================================================
# ANALISIS OUT-OF-FOLD: 5 TEST CONCATENADOS
# ==============================================================
# Cada celula de las 658 de train aparece exactamente una vez
# como test a lo largo de las 5 iteraciones. Concatenando se
# recupera una prediccion "out-of-fold" para las 658 celulas.
# ==============================================================

cat(sprintf("\n%s\n", strrep("=", 55)))
cat("ANALISIS OUT-OF-FOLD (5 TEST CONCATENADOS)\n")
cat(sprintf("%s\n\n", strrep("=", 55)))

all_cells     <- unlist(lapply(resultados_oof, function(x) x$cells))
all_pred      <- do.call(rbind, lapply(resultados_oof, function(x) x$pred))
all_pred_base <- do.call(rbind, lapply(resultados_oof, function(x) x$pred_base))
all_viab      <- do.call(rbind, lapply(resultados_oof, function(x) x$viab))

rownames(all_pred)      <- all_cells
rownames(all_pred_base) <- all_cells
rownames(all_viab)      <- all_cells

# Reordenar
ord <- order(rownames(all_pred))
all_pred      <- all_pred[ord, ]
all_pred_base <- all_pred_base[ord, ]
all_viab      <- all_viab[ord, ]

cat(sprintf("Total celulas out-of-fold: %d\n", nrow(all_pred)))

# ----------------------------------------------------------
# Metricas globales OOF
# ----------------------------------------------------------
rmse_base_oof  <- sqrt(mean((all_viab - all_pred_base)^2))
rmse_final_oof <- sqrt(mean((all_viab - all_pred)^2))

cor_gen_base_oof  <- diag(cor(all_viab, all_pred_base))
cor_gen_final_oof <- diag(cor(all_viab, all_pred))

cor_cel_base_oof  <- diag(cor(t(all_viab), t(all_pred_base)))
cor_cel_final_oof <- diag(cor(t(all_viab), t(all_pred)))

cat(sprintf("RMSE OOF:        %.4f -> %.4f (%.1f%%)\n",
            rmse_base_oof, rmse_final_oof,
            100*(rmse_base_oof - rmse_final_oof)/rmse_base_oof))
cat(sprintf("Cor gen mediana: %.4f -> %.4f\n",
            median(cor_gen_base_oof, na.rm = TRUE),
            median(cor_gen_final_oof, na.rm = TRUE)))
cat(sprintf("Cor cel media:   %.4f -> %.4f\n",
            mean(cor_cel_base_oof, na.rm = TRUE),
            mean(cor_cel_final_oof, na.rm = TRUE)))

# ----------------------------------------------------------
# Por linaje tumoral (correlacion por celula)
# ----------------------------------------------------------
oof_df <- data.frame(
  ModelID   = rownames(all_pred),
  cor_base  = cor_cel_base_oof,
  cor_final = cor_cel_final_oof,
  mae_base  = rowMeans(abs(all_viab - all_pred_base)),
  mae_final = rowMeans(abs(all_viab - all_pred))
) %>%
  merge(metadata[, c("ModelID", "OncotreeLineage", "FormulationID")],
        by = "ModelID", all.x = TRUE) %>%
  mutate(medio_base = trimws(gsub("\\s*\\+.*", "", FormulationID)))

cat("\n=== CORRELACION POR CELULA POR LINAJE (OOF) ===\n")
oof_df %>%
  group_by(OncotreeLineage) %>%
  summarise(
    n          = n(),
    cor_base   = round(mean(cor_base,  na.rm = TRUE), 4),
    cor_final  = round(mean(cor_final, na.rm = TRUE), 4),
    mae_base   = round(mean(mae_base),  4),
    mae_final  = round(mean(mae_final), 4)
  ) %>%
  arrange(desc(n)) %>%
  print(n = 25)

cat("\n=== CORRELACION POR CELULA POR MEDIO BASE (OOF, n>=10) ===\n")
oof_df %>%
  filter(!is.na(medio_base) & medio_base != "") %>%
  group_by(medio_base) %>%
  summarise(
    n          = n(),
    cor_base   = round(mean(cor_base,  na.rm = TRUE), 4),
    cor_final  = round(mean(cor_final, na.rm = TRUE), 4),
    mae_base   = round(mean(mae_base),  4),
    mae_final  = round(mean(mae_final), 4)
  ) %>%
  filter(n >= 10) %>%
  arrange(desc(n)) %>%
  print()

# ----------------------------------------------------------
# Boxplots comparativos (correlacion por celula)
# ----------------------------------------------------------
library(ggplot2)
library(tidyr)

oof_long <- oof_df %>%
  pivot_longer(cols = c(cor_base, cor_final),
               names_to  = "etapa",
               values_to = "cor") %>%
  mutate(etapa = ifelse(etapa == "cor_base",
                        "Baseline", "Pipeline final"))

# Por linaje (linajes principales)
linajes_principales <- oof_df %>%
  count(OncotreeLineage, sort = TRUE) %>%
  filter(n >= 15) %>%
  pull(OncotreeLineage)

p_linaje <- ggplot(oof_long %>%
         filter(OncotreeLineage %in% linajes_principales),
       aes(x    = reorder(OncotreeLineage, -cor, median),
           y    = cor,
           fill = etapa)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  scale_fill_manual(values = c("Baseline"      = "#D85A30",
                                "Pipeline final" = "#1D9E75")) +
  labs(title = "Correlacion por celula segun linaje (OOF, 658 celulas)",
       x     = "Linaje tumoral",
       y     = "Correlacion por celula",
       fill  = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        legend.position = "top")

print(p_linaje)

# Por medio de cultivo
medios_principales_oof <- oof_df %>%
  filter(!is.na(medio_base) & medio_base != "") %>%
  count(medio_base, sort = TRUE) %>%
  filter(n >= 10) %>%
  pull(medio_base)

p_medio <- ggplot(oof_long %>%
         filter(medio_base %in% medios_principales_oof),
       aes(x    = reorder(medio_base, -cor, median),
           y    = cor,
           fill = etapa)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  scale_fill_manual(values = c("Baseline"       = "#D85A30",
                                "Pipeline final" = "#1D9E75")) +
  labs(title = "Correlacion por celula segun medio de cultivo (OOF, 658 celulas)",
       x     = "Medio base",
       y     = "Correlacion por celula",
       fill  = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        legend.position = "top")

print(p_medio)

# ----------------------------------------------------------
# Brecha hematologico vs solido (OOF)
# ----------------------------------------------------------
oof_df <- oof_df %>%
  mutate(tipo_tumoral = ifelse(OncotreeLineage %in%
                                 c("Lymphoid", "Myeloid"),
                               "Hematológico", "Tumor sólido"))

cat("\n=== HEMATOLOGICO VS SOLIDO (OOF, correlacion por celula) ===\n")
oof_df %>%
  group_by(tipo_tumoral) %>%
  summarise(
    n         = n(),
    cor_base  = round(mean(cor_base,  na.rm = TRUE), 4),
    cor_final = round(mean(cor_final, na.rm = TRUE), 4),
    mae_base  = round(mean(mae_base),  4),
    mae_final = round(mean(mae_final), 4)
  ) %>%
  print()

wt_cor_base  <- wilcox.test(cor_base  ~ tipo_tumoral, data = oof_df)
wt_cor_final <- wilcox.test(cor_final ~ tipo_tumoral, data = oof_df)

cat(sprintf("Wilcoxon cor_base:  p = %.2e\n", wt_cor_base$p.value))
cat(sprintf("Wilcoxon cor_final: p = %.2e\n", wt_cor_final$p.value))

# ==============================================================
# CORRELACION POR GEN DENTRO DE CADA LINAJE (OOF)
# ==============================================================

# Anadir linaje a las matrices OOF
metadata_oof <- data.frame(ModelID = rownames(all_pred)) %>%
  merge(metadata[, c("ModelID", "OncotreeLineage")], by = "ModelID", all.x = TRUE)

# Funcion para calcular cor por gen dentro de un subgrupo
cor_gen_por_linaje <- function(linaje, min_n = 25) {
  cells <- metadata_oof$ModelID[metadata_oof$OncotreeLineage == linaje]
  if (length(cells) < min_n) return(NULL)
  
  viab_sub <- all_viab[cells, ]
  pred_base_sub  <- all_pred_base[cells, ]
  pred_final_sub <- all_pred[cells, ]
  
  cor_base  <- diag(cor(viab_sub, pred_base_sub))
  cor_final <- diag(cor(viab_sub, pred_final_sub))
  
  data.frame(
    linaje = linaje,
    n      = length(cells),
    cor_gen_base_mediana  = median(cor_base,  na.rm = TRUE),
    cor_gen_final_mediana = median(cor_final, na.rm = TRUE)
  )
}

# Linajes con al menos 15 celulas en OOF
linajes_validos <- metadata_oof %>%
  count(OncotreeLineage, sort = TRUE) %>%
  filter(n >= 15) %>%
  pull(OncotreeLineage)

resultado_cor_gen_linaje <- do.call(rbind,
                                    lapply(linajes_validos, cor_gen_por_linaje))

print(resultado_cor_gen_linaje)

