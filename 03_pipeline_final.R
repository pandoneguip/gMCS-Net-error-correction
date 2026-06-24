# ==============================================================
# PIPELINE FINAL — FOLD 2 COMPLETO
# ==============================================================
# Partición oficial del fold 2:
#   Train = 658 cell lines
#   Val   = 222 cell lines  (elegir hiperparámetros)
#   Test  = 223 cell lines  (evaluación final — una sola vez al final)
#
# Capas:
#   1. Median polish         (alpha en VAL)
#   2. Pseudoinversa + asinh (k, gamma, alpha en VAL)
#   3. Corrección expresión  (alpha en VAL)
#
# ESQUEMA ESTRICTO:
#   - TEST no se toca hasta el final
#   - Todos los hiperparámetros se eligen en VAL
#   - TEST se evalúa una única vez al final
#
# PREREQUISITO: train_mat, val_mat, test_mat,
#               pred_train, pred_val, pred_test,
#               Expression en memoria
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
# BASELINE 
# ==============================================================

cat("========================================\n")
cat("BASELINE (VAL)\n")
cat("========================================\n")

cat(sprintf("RMSE val:              %.4f\n",
            sqrt(mean((val_mat - pred_val)^2))))
cat(sprintf("Error cuadratico val:  %.2f\n",
            sum((val_mat - pred_val)^2)))
cor_por_celula(val_mat, pred_val, "Cor celula VAL baseline")
cor_por_gen(val_mat,    pred_val, "Cor gen VAL baseline")

# ==============================================================
# APRENDER PARÁMETROS EN TRAIN
# ==============================================================

error_train <- train_mat - pred_train

# ==============================================================
# CAPA 1: MEDIAN POLISH — ELEGIR ALPHA EN VAL
# ==============================================================

cat("\n========================================\n")
cat("CAPA 1: MEDIAN POLISH\n")
cat("========================================\n")

mp <- medpolish(error_train, maxiter = 100,
                 eps = 0.01, trace.iter = FALSE)

sesgo_genes <- mp$col + mp$overall
cat(sprintf("Efecto global: %.4f\n", mp$overall))
cat(sprintf("Sesgo maximo por gen: %.4f\n", max(abs(mp$col))))

# Seleccionar alpha en VAL
alphas_mp <- seq(0, 2, length.out = 100)
err_mp_val <- sapply(alphas_mp, function(a) {
  pred_new <- t(t(pred_val) + a * sesgo_genes)
  sum((val_mat - pred_new)^2)
})
alpha_mp_opt <- alphas_mp[which.min(err_mp_val)]

cat(sprintf("Alpha optimo MP: %.4f\n", alpha_mp_opt))
cat(sprintf("Error cuadratico VAL: %.2f -> %.2f\n",
            err_mp_val[1], min(err_mp_val)))

# Aplicar solo a VAL para siguiente capa
pred_val_mp <- t(t(pred_val) + alpha_mp_opt * sesgo_genes)

cat(sprintf("RMSE val tras MP: %.4f\n",
            sqrt(mean((val_mat - pred_val_mp)^2))))
cor_por_gen(val_mat, pred_val_mp, "Cor gen VAL tras MP")

# ==============================================================
# CAPA 2: PSEUDOINVERSA + ASINH — ELEGIR k, GAMMA, ALPHA EN VAL
# ==============================================================

cat("\n========================================\n")
cat("CAPA 2: PSEUDOINVERSA + ASINH\n")
cat("========================================\n")

svd_check <- svd(pred_train, nu = 0, nv = 0)
cat(sprintf("Num condicion pred_train: %.2e\n",
            max(svd_check$d) / min(svd_check$d[svd_check$d > 0])))

ks <- c(10, 20, 50, 100, 150, 200, 250, 300,
        min(nrow(pred_train), ncol(pred_train)))
ks          <- unique(ks)
gammas      <- c(1, 2, 5, 10)
alphas_pinv <- seq(0, 2, length.out = 100)

mejor_err      <- Inf
gamma_opt      <- 1
alpha_pinv_opt <- 0
k_opt          <- 50

cat("Grid search k, gamma, alpha en VAL...\n")

for (k in ks) {
  for (g in gammas) {
    err_train_asinh <- (1/g) * asinh(g * error_train)
    A_val       <- pred_val_mp %*% pinv(pred_train, k = k)
    err_val_est <- A_val %*% err_train_asinh

    if (max(abs(err_val_est)) > 10) next

    for (a in alphas_pinv) {
      pred_new <- pred_val_mp + a * err_val_est
      err <- sum((val_mat - pred_new)^2)
      if (err < mejor_err) {
        mejor_err      <- err
        gamma_opt      <- g
        alpha_pinv_opt <- a
        k_opt          <- k
      }
    }
  }
  cat(sprintf("k=%3d completado\n", k))
}

cat(sprintf("\nK optimo:     %d\n",   k_opt))
cat(sprintf("Gamma optimo: %d\n",   gamma_opt))
cat(sprintf("Alpha optimo: %.4f\n", alpha_pinv_opt))
cat(sprintf("Error VAL optimo: %.2f\n", mejor_err))

# Aplicar solo a VAL para siguiente capa
err_train_asinh_opt <- (1/gamma_opt) * asinh(gamma_opt * error_train)
A_val           <- pred_val_mp %*% pinv(pred_train, k = k_opt)
err_val_est     <- A_val %*% err_train_asinh_opt
pred_val_pinv   <- pred_val_mp + alpha_pinv_opt * err_val_est

cat(sprintf("RMSE val tras pinv: %.4f\n",
            sqrt(mean((val_mat - pred_val_pinv)^2))))
cor_por_gen(val_mat, pred_val_pinv, "Cor gen VAL tras pinv")

# ==============================================================
# CAPA 3: EXPRESIÓN GÉNICA — ELEGIR ALPHA EN VAL
# ==============================================================

cat("\n========================================\n")
cat("CAPA 3: EXPRESION GENICA\n")
cat("========================================\n")

cells_expr_train <- base::intersect(rownames(train_mat), rownames(Expression))
cells_expr_val   <- base::intersect(rownames(val_mat),   rownames(Expression))
cells_expr_test  <- base::intersect(rownames(test_mat),  rownames(Expression))

expr_means <- colMeans(Expression[cells_expr_train, ], na.rm = TRUE)

corr_expr <- sapply(colnames(error_train), function(g) {
  if (!g %in% colnames(Expression)) return(NA)
  cor(Expression[cells_expr_train, g],
      error_train[cells_expr_train, g],
      use = "complete.obs")
})

genes_signal <- names(corr_expr)[
  abs(corr_expr) > 0.10 & !is.na(corr_expr)]

cat(sprintf("Genes con señal expresion (|cor|>0.10): %d\n",
            length(genes_signal)))

beta_expr <- sapply(colnames(error_train), function(g) {
  if (!g %in% colnames(Expression)) return(NA)
  expr_g  <- Expression[cells_expr_train, g]
  error_g <- error_train[cells_expr_train, g]
  ok <- !is.na(expr_g) & !is.na(error_g)
  if (sum(ok) < 20) return(NA)
  cov(expr_g[ok], error_g[ok]) / var(expr_g[ok])
})
#Para cada gen con señal, mira cuánto se desvía la expresión de esa célula 
#respecto a la media de train, multiplica por la pendiente y por alpha, 
#limita el resultado a ±0.1, y lo suma a la predicción.
aplicar_expr <- function(pred_mat, cells_expr, alpha_e) {
  pred_corr <- pred_mat
  for (g in genes_signal) {
    if (!g %in% colnames(Expression)) next
    if (is.na(beta_expr[g])) next
    expr_g <- Expression[cells_expr, g]
    expr_g[is.na(expr_g)] <- expr_means[g]
    correccion <- alpha_e * beta_expr[g] * (expr_g - expr_means[g])
    correccion <- pmax(pmin(correccion, 0.1), -0.1)
    pred_corr[cells_expr, g] <- pred_corr[cells_expr, g] + correccion
  }
  pred_corr
}

# Seleccionar alpha en VAL minimizando error cuadratico
alphas_expr  <- seq(0, 5, length.out = 200)
err_expr_val <- sapply(alphas_expr, function(a) {
  pred_new <- aplicar_expr(pred_val_pinv, cells_expr_val, a)
  sum((val_mat - pred_new)^2)
})
alpha_expr_opt <- alphas_expr[which.min(err_expr_val)]

cat(sprintf("Alpha optimo expresion: %.4f\n", alpha_expr_opt))
cat(sprintf("Error VAL antes expr: %.2f -> %.2f\n",
            sum((val_mat - pred_val_pinv)^2), min(err_expr_val)))

# Aplicar solo a VAL para verificar
pred_val_final <- aplicar_expr(pred_val_pinv, cells_expr_val,
                                alpha_expr_opt)

cat(sprintf("RMSE val final: %.4f\n",
            sqrt(mean((val_mat - pred_val_final)^2))))
cor_por_gen(val_mat, pred_val_final, "Cor gen VAL final")

# ==============================================================
# HIPERPARÁMETROS SELECCIONADOS
# ==============================================================

cat("\n========================================\n")
cat("HIPERPARAMETROS SELECCIONADOS EN VAL\n")
cat("========================================\n")
cat(sprintf("Alpha MP:     %.4f\n", alpha_mp_opt))
cat(sprintf("K optimo:     %d\n",   k_opt))
cat(sprintf("Gamma optimo: %d\n",   gamma_opt))
cat(sprintf("Alpha pinv:   %.4f\n", alpha_pinv_opt))
cat(sprintf("Alpha expr:   %.4f\n", alpha_expr_opt))
cat(sprintf("Genes signal: %d\n",   length(genes_signal)))

# ==============================================================
# EVALUACIÓN FINAL EN TEST — UNA SOLA VEZ
# ==============================================================

cat("\n========================================\n")
cat("EVALUACION FINAL EN TEST\n")
cat("========================================\n")

# Aplicar pipeline completo a TEST con parametros de VAL
pred_test_mp    <- t(t(pred_test) + alpha_mp_opt * sesgo_genes)

A_test          <- pred_test_mp %*% pinv(pred_train, k = k_opt)
err_test_est    <- A_test %*% err_train_asinh_opt
pred_test_pinv  <- pred_test_mp + alpha_pinv_opt * err_test_est

pred_test_final <- aplicar_expr(pred_test_pinv, cells_expr_test,
                                 alpha_expr_opt)

# Métricas finales
cat(sprintf("RMSE baseline TEST:  %.4f\n",
            sqrt(mean((test_mat - pred_test)^2))))
cat(sprintf("RMSE final TEST:     %.4f\n",
            sqrt(mean((test_mat - pred_test_final)^2))))
cat(sprintf("Mejora RMSE TEST:    %.2f%%\n",
            100*(sqrt(mean((test_mat - pred_test)^2)) -
                   sqrt(mean((test_mat - pred_test_final)^2))) /
              sqrt(mean((test_mat - pred_test)^2))))

cat(sprintf("\nError cuadratico TEST baseline: %.2f\n",
            sum((test_mat - pred_test)^2)))
cat(sprintf("Error cuadratico TEST final:    %.2f\n",
            sum((test_mat - pred_test_final)^2)))

cat("\n--- Correlacion por gen TEST ---\n")
cors_gen_base  <- cor_por_gen(test_mat, pred_test,
                               "Baseline TEST")
cors_gen_final <- cor_por_gen(test_mat, pred_test_final,
                               "Final TEST")

cat("\n--- Correlacion por celula TEST ---\n")
cors_cel_base  <- cor_por_celula(test_mat, pred_test,
                                  "Baseline TEST")
cors_cel_final <- cor_por_celula(test_mat, pred_test_final,
                                  "Final TEST")

cat("\n--- Test estadistico (por celula) ---\n")
t_test <- t.test(cors_cel_base, cors_cel_final)
cat(sprintf("p=%.2e  %.4f -> %.4f\n",
            t_test$p.value,
            mean(cors_cel_base), mean(cors_cel_final)))

# ==============================================================
# RESUMEN FINAL
# ==============================================================

cat("\n========================================\n")
cat("RESUMEN\n")
cat("========================================\n")

resumen <- data.frame(
  Conjunto = c("VAL", "VAL", "TEST", "TEST"),
  Estado   = c("Baseline", "Final", "Baseline", "Final"),
  RMSE = c(
    sqrt(mean((val_mat  - pred_val)^2)),
    sqrt(mean((val_mat  - pred_val_final)^2)),
    sqrt(mean((test_mat - pred_test)^2)),
    sqrt(mean((test_mat - pred_test_final)^2))
  ),
  Cor_gen_media = c(
    mean(diag(cor(val_mat,  pred_val))),
    mean(diag(cor(val_mat,  pred_val_final))),
    mean(diag(cor(test_mat, pred_test))),
    mean(diag(cor(test_mat, pred_test_final)))
  ),
  Cor_gen_mediana = c(
    median(diag(cor(val_mat,  pred_val))),
    median(diag(cor(val_mat,  pred_val_final))),
    median(diag(cor(test_mat, pred_test))),
    median(diag(cor(test_mat, pred_test_final)))
  ),
  Cor_cel_media = c(
    mean(diag(cor(t(val_mat),  t(pred_val)))),
    mean(diag(cor(t(val_mat),  t(pred_val_final)))),
    mean(diag(cor(t(test_mat), t(pred_test)))),
    mean(diag(cor(t(test_mat), t(pred_test_final))))
  )
)
resumen$mejora_rmse <- round(
  100 * (resumen$RMSE[c(1,1,3,3)] - resumen$RMSE) /
    resumen$RMSE[c(1,1,3,3)], 2)
print(resumen)

# Verificacion overfitting
cat("\n--- Verificacion overfitting ---\n")
pred_train_mp    <- t(t(pred_train) + alpha_mp_opt * sesgo_genes)
A_train          <- pred_train_mp %*% pinv(pred_train, k = k_opt)
err_train_est    <- A_train %*% err_train_asinh_opt
pred_train_pinv  <- pred_train_mp + alpha_pinv_opt * err_train_est
pred_train_final <- aplicar_expr(pred_train_pinv,
                                  cells_expr_train, alpha_expr_opt)

rmse_train <- sqrt(mean((train_mat - pred_train_final)^2))
rmse_val_f <- sqrt(mean((val_mat   - pred_val_final)^2))
rmse_test_f<- sqrt(mean((test_mat  - pred_test_final)^2))

cat(sprintf("RMSE train: %.4f\n", rmse_train))
cat(sprintf("RMSE val:   %.4f\n", rmse_val_f))
cat(sprintf("RMSE test:  %.4f\n", rmse_test_f))
cat(sprintf("Ratio val/train:  %.3f\n", rmse_val_f / rmse_train))
cat(sprintf("Ratio test/val:   %.3f\n", rmse_test_f / rmse_val_f))

# ==============================================================
# ANÁLISIS POST-CORRECCIÓN — CAPÍTULO 4
# ==============================================================

library(dplyr)

# Calcular error final por cell line
# pred_test_final ya está en memoria del pipeline_final.R

error_test_base  <- pred_test  - test_mat
error_test_final <- pred_test_final - test_mat

# MAE por cell line antes y después
mae_base  <- rowMeans(abs(error_test_base),  na.rm = TRUE)
mae_final <- rowMeans(abs(error_test_final), na.rm = TRUE)

# Data frame con métricas y metadata
post_corr <- data.frame(
  ModelID   = rownames(test_mat),
  mae_base  = mae_base,
  mae_final = mae_final,
  mejora    = mae_base - mae_final
) %>%
  merge(metadata[, c("ModelID", "OncotreeLineage",
                     "OncotreePrimaryDisease",
                     "GrowthPattern", "FormulationID",
                     "PediatricModelType")],
        by = "ModelID", all.x = TRUE) %>%
  mutate(
    medio_base = trimws(gsub("\\s*\\+.*", "", FormulationID)),
    tipo_tumoral = ifelse(OncotreeLineage %in%
                            c("Lymphoid", "Myeloid"),
                          "Hematológico", "Tumor sólido")
  )

# ----------------------------------------------------------
# 1. MAE POR LINAJE ANTES Y DESPUÉS
# ----------------------------------------------------------
cat("=== MAE POR LINAJE ANTES Y DESPUÉS ===\n")
post_corr %>%
  group_by(OncotreeLineage) %>%
  summarise(
    n          = n(),
    mae_base   = round(mean(mae_base),  4),
    mae_final  = round(mean(mae_final), 4),
    mejora_pct = round(100*(mean(mae_base) - mean(mae_final)) /
                         mean(mae_base), 1)
  ) %>%
  arrange(desc(mae_base)) %>%
  print(n = 20)

# ----------------------------------------------------------
# 2. HEMATOLÓGICOS VS SÓLIDOS ANTES Y DESPUÉS
# ----------------------------------------------------------
cat("\n=== HEMATOLÓGICOS VS SÓLIDOS ===\n")
post_corr %>%
  group_by(tipo_tumoral) %>%
  summarise(
    n          = n(),
    mae_base   = round(mean(mae_base),  4),
    mae_final  = round(mean(mae_final), 4),
    mejora_pct = round(100*(mean(mae_base) - mean(mae_final)) /
                         mean(mae_base), 1)
  ) %>%
  print()

# Test Wilcoxon: diferencia hemato vs solido ANTES
wt_base <- wilcox.test(
  mae_base ~ tipo_tumoral, data = post_corr)
cat(sprintf("Wilcoxon ANTES: p = %.2e\n", wt_base$p.value))

# Test Wilcoxon: diferencia hemato vs solido DESPUÉS
wt_final <- wilcox.test(
  mae_final ~ tipo_tumoral, data = post_corr)
cat(sprintf("Wilcoxon DESPUÉS: p = %.2e\n", wt_final$p.value))

# ----------------------------------------------------------
# 3. MAE POR MEDIO DE CULTIVO ANTES Y DESPUÉS
# ----------------------------------------------------------
cat("\n=== MAE POR MEDIO BASE ANTES Y DESPUÉS ===\n")
post_corr %>%
  filter(!is.na(medio_base) & medio_base != "") %>%
  group_by(medio_base) %>%
  summarise(
    n          = n(),
    mae_base   = round(mean(mae_base),  4),
    mae_final  = round(mean(mae_final), 4),
    mejora_pct = round(100*(mean(mae_base) - mean(mae_final)) /
                         mean(mae_base), 1)
  ) %>%
  filter(n >= 5) %>%
  arrange(desc(mae_base)) %>%
  print()

# ----------------------------------------------------------
# 4. BOXPLOT COMPARATIVO POR LINAJE
# ----------------------------------------------------------
library(ggplot2)
library(tidyr)

post_corr_long <- post_corr %>%
  pivot_longer(cols = c(mae_base, mae_final),
               names_to  = "etapa",
               values_to = "mae") %>%
  mutate(etapa = ifelse(etapa == "mae_base",
                        "Baseline", "Pipeline final"))

ggplot(post_corr_long %>%
         filter(OncotreeLineage %in%
                  c("Lymphoid", "Myeloid",
                    "Lung", "Skin", "CNS/Brain",
                    "Breast", "Bowel")),
       aes(x    = reorder(OncotreeLineage, -mae, median),
           y    = mae,
           fill = etapa)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  scale_fill_manual(values = c("Baseline"      = "#D85A30",
                               "Pipeline final" = "#1D9E75")) +
  labs(title = "MAE por linaje tumoral: antes y después del pipeline",
       x     = "Linaje tumoral",
       y     = "MAE",
       fill  = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# ----------------------------------------------------------
# 5. BOXPLOT COMPARATIVO POR MEDIO DE CULTIVO
# ----------------------------------------------------------
medios_principales <- c("RPMI", "IMDM", "DMEM",
                        "MEM", "DMEM:F12", "EMEM", "F12")

post_corr_long %>%
  filter(medio_base %in% medios_principales) %>%
  ggplot(aes(x    = reorder(medio_base, -mae, median),
             y    = mae,
             fill = etapa)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  scale_fill_manual(values = c("Baseline"       = "#D85A30",
                               "Pipeline final" = "#1D9E75")) +
  labs(title = "MAE por medio de cultivo: antes y después del pipeline",
       x     = "Medio base",
       y     = "MAE",
       fill  = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        legend.position = "top")

# ----------------------------------------------------------
# 6. ¿SE REDUCE LA BRECHA HEMATOLÓGICO VS SÓLIDO?
# ----------------------------------------------------------
cat("\n=== BRECHA HEMATOLÓGICO VS SÓLIDO ===\n")
resumen_brecha <- post_corr %>%
  group_by(tipo_tumoral) %>%
  summarise(
    mae_base  = mean(mae_base),
    mae_final = mean(mae_final)
  )

brecha_base  <- diff(resumen_brecha$mae_base)
brecha_final <- diff(resumen_brecha$mae_final)

cat(sprintf("Brecha antes:  %.4f\n", abs(brecha_base)))
cat(sprintf("Brecha después: %.4f\n", abs(brecha_final)))
cat(sprintf("Reducción brecha: %.1f%%\n",
            100*(abs(brecha_base) - abs(brecha_final)) /
              abs(brecha_base)))

# Correlación por célula, baseline y final
cor_cel_base  <- diag(cor(t(test_mat), t(pred_test)))
cor_cel_final <- diag(cor(t(test_mat), t(pred_test_final)))

post_corr$cor_base  <- cor_cel_base
post_corr$cor_final <- cor_cel_final

# Por linaje
post_corr %>%
  group_by(OncotreeLineage) %>%
  summarise(n = n(),
            cor_base  = mean(cor_base),
            cor_final = mean(cor_final)) %>%
  arrange(desc(cor_base)) %>%
  print(n=20)

# Por medio
post_corr %>%
  filter(!is.na(medio_base) & medio_base != "") %>%
  group_by(medio_base) %>%
  summarise(n = n(),
            cor_base  = mean(cor_base),
            cor_final = mean(cor_final)) %>%
  filter(n>=5) %>%
  arrange(desc(cor_base)) %>%
  print()

# Boxplots
post_corr_long2 <- post_corr %>%
  pivot_longer(cols = c(cor_base, cor_final),
               names_to = "etapa", values_to = "cor") %>%
  mutate(etapa = ifelse(etapa=="cor_base","Baseline","Pipeline final"))

ggplot(post_corr_long2 %>%
         filter(OncotreeLineage %in% c("Lymphoid","Myeloid","Lung","Skin","CNS/Brain","Breast","Bowel")),
       aes(x=reorder(OncotreeLineage,-cor,median), y=cor, fill=etapa)) +
  geom_boxplot(alpha=0.8) +
  scale_fill_manual(values=c("Baseline"="#D85A30","Pipeline final"="#1D9E75")) +
  labs(title="Correlación por línea celular: antes y después", x="Linaje", y="Correlación") +
  theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="top")

