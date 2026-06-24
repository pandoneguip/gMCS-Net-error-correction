# ==============================================================
# ANÁLISIS DE ERROR POR MEDIO DE CULTIVO
# ==============================================================
# Se usa MAE (Mean Absolute Error) en vez de RMSE por ser más
# consistente con la filosofía robusta del median polish
# que opera sobre medianas.
#
# A la hora de analizar los medios de cultivos primero queremos
# saber a qué células suelen pertenecer:
#   RPMI  --> Células hematológicas (leucemias, linfomas)
#   DMEM  --> Células adherentes de tumores sólidos
#   MEM   --> Medio para tumores sólidos
#   IMDM  --> Hematológicas (células mieloides y médula ósea)
# ==============================================================


# Cargamos el metadata de DepMap con información de cada cell line
# (tejido, tipo de cultivo, medio de cultivo, etc.)
setwd("G:/Mi unidad/PFG/gMCS-Net_Angel/code/R")
metadata <- read.csv("Model.csv", stringsAsFactors = FALSE)

# Comprobación rápida: ver qué columnas tenemos disponibles
# para el análisis de medios de cultivo
head(metadata[, c("ModelID", "FormulationID", "OnboardedMedia",
                  "GrowthPattern")], 10)

# ----------------------------------------------------------
# Calcular métricas de error por cell line
# sobre el residuo del median polish
# ----------------------------------------------------------

# Creamos un data frame con una fila por cell line
# usando los nombres de fila de la matriz de residuos
error_by_media <- data.frame(
  ModelID = rownames(error_train_mp)
)

# MAE (Mean Absolute Error): promedio del valor absoluto del error
# Es la métrica principal porque es consistente con el median polish
# que trabaja con medianas — ambas son robustas ante outliers
error_by_media$mae   <- apply(error_train_mp, 1,
                              function(x) mean(abs(x)))

# MedAE (Median Absolute Error): mediana del valor absoluto del error
# Es la métrica más robusta posible — completamente inmune a outliers
# La usamos como comprobación de que MAE y MedAE cuentan la misma historia
error_by_media$medae <- apply(error_train_mp, 1,
                              function(x) median(abs(x)))

# Media del error sin valor absoluto: detecta si el modelo
# sobreestima (positivo) o subestima (negativo) sistemáticamente
# en cada cell line
error_by_media$mean_error <- apply(error_train_mp, 1, mean)

# Cruzamos las métricas de error con el metadata
# para añadir el FormulationID (medio de cultivo) y otras
# características de cada cell line
error_by_media <- merge(error_by_media,
                        metadata[, c("ModelID", "FormulationID",
                                     "OnboardedMedia", "GrowthPattern",
                                     "OncotreeLineage", "OncotreePrimaryDisease","SerumFreeMedia",
                                     "PlateCoating", "PediatricModelType")],
                        by = "ModelID", all.x = TRUE)

# ----------------------------------------------------------
# LIMPIEZA Y AGRUPACIÓN POR MEDIO BASE
# ----------------------------------------------------------
# El FormulationID tiene nombres muy específicos con suplementos
# por ejemplo: "RPMI + 10% FBS + 2mM Glutamine + 25mM HEPES"
# Esto fragmenta los grupos en decenas de categorías pequeñas.
# Agrupamos por medio base (RPMI, DMEM, MEM...) para tener
# grupos estadísticamente robustos (n suficiente para tests).

# Primero limpiamos los strings vacíos "" que no son NAs
# pero tampoco tienen información útil — los ponemos a NA
error_by_media$FormulationID_clean <- ifelse(
  error_by_media$FormulationID == "" |
    is.na(error_by_media$FormulationID),
  NA,
  error_by_media$FormulationID
)

# Extraemos el medio base eliminando todo lo que viene después del +
# gsub reemplaza el patrón "\\s*\\+.*" (espacio + y todo lo que sigue)
# por nada "", dejando solo el nombre base del medio
error_by_media$medio_base <- gsub("\\s*\\+.*", "",
                                  error_by_media$FormulationID_clean)

# trimws elimina espacios en blanco al inicio y final del string
error_by_media$medio_base <- trimws(error_by_media$medio_base)

# Vemos cuántas cell lines hay por cada medio base
# para decidir el umbral mínimo del análisis
cat("=== DISTRIBUCIÓN POR MEDIO BASE ===\n")
print(sort(table(error_by_media$medio_base), decreasing = TRUE))

# ----------------------------------------------------------
# Resumen estadístico por medio base
# ----------------------------------------------------------
# Calculamos el MAE medio, desviación estándar y sesgo
# para cada medio base, filtrando los que tienen menos de 20
# cell lines porque con tan pocos datos los resultados no son fiables
medios_base_summary <- error_by_media %>%
  filter(!is.na(medio_base)) %>%           # eliminamos los sin medio
  group_by(medio_base) %>%                 # agrupamos por medio base
  summarise(
    n           = n(),                     # cuántas cell lines hay
    mae_medio   = mean(mae),               # error medio del grupo
    mae_sd      = sd(mae),                 # variabilidad del error
    medae_medio = mean(medae),             # error mediano medio
    sesgo_medio = mean(mean_error)         # dirección del sesgo
  ) %>%
  filter(n >= 20) %>%                      # solo medios con n suficiente
  arrange(desc(mae_medio))                 # ordenar de mayor a menor error

cat("\n=== RESUMEN POR MEDIO BASE (n >= 20) ===\n")
print(medios_base_summary)

# ----------------------------------------------------------
# Preparar datos para el boxplot
# ----------------------------------------------------------
# Filtramos el data frame principal para quedarnos solo con
# las cell lines que pertenecen a medios con n >= 20
error_filtrado_base <- error_by_media %>%
  filter(!is.na(medio_base) &
           medio_base %in% medios_base_summary$medio_base)

cat(sprintf("\nMedios base con n >= 20: %d\n",
            nrow(medios_base_summary)))
cat(sprintf("Cell lines incluidas en el boxplot: %d\n",
            nrow(error_filtrado_base)))

# ----------------------------------------------------------
# Boxplot MAE por medio base
# ----------------------------------------------------------
# reorder() ordena los medios de mayor a menor MAE mediano
# para que la visualización sea más informativa
ggplot(error_filtrado_base,
       aes(x    = reorder(medio_base, -mae, median),
           y    = mae,
           fill = medio_base)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +  # caja con outliers visibles
  geom_jitter(width = 0.15, alpha = 0.3, size = 1) + # puntos individuales
  labs(title = "MAE del modelo por medio base de cultivo (n ≥ 20)",
       x     = "Medio base",
       y     = "MAE (residuo median polish)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# ----------------------------------------------------------
# Test estadístico: Kruskal-Wallis + Dunn post-hoc
# ----------------------------------------------------------
# Usamos Kruskal-Wallis en vez de ANOVA por dos razones:
#   1. Los datos no siguen una distribución normal
#   2. Es consistente con el median polish que usa medianas
#
# El Kruskal-Wallis es el "semáforo" — responde si hay
# diferencias reales en algún sitio (p<0.05) o todo es ruido
cat("\n=== TEST KRUSKAL-WALLIS (MAE ~ medio_base) ===\n")
kw_base <- kruskal.test(mae ~ medio_base,
                        data = error_filtrado_base)
print(kw_base)

# Si el Kruskal-Wallis es significativo, hacemos Dunn post-hoc
# El Dunn compara todos los pares posibles (RPMI vs DMEM, etc.)
# y nos dice exactamente QUÉ pares difieren significativamente
# La corrección BH controla los falsos positivos por múltiples comparaciones
if (kw_base$p.value < 0.05) {
  cat("\nDiferencias significativas. Test de Dunn post-hoc:\n")
  if (requireNamespace("dunn.test", quietly = TRUE)) {
    dunn.test::dunn.test(
      x      = error_filtrado_base$mae,
      g      = error_filtrado_base$medio_base,
      method = "BH"   # corrección de Benjamini-Hochberg
    )
  } else {
    cat("Instala dunn.test: install.packages('dunn.test')\n")
  }
}

# ==============================================================
# ANÁLISIS DE TUMORES EN MEDIOS RPMI e IMDM
# ==============================================================
# RPMI e IMDM son los que tienen MAE significativamente mayor
# según el Dunn. Aquí analizamos qué tipos tumorales contienen
# para entender si el error es por el medio o por el tipo tumoral.

# Filtramos solo las cell lines de RPMI e IMDM
# ordenadas de mayor a menor MAE para ver las más problemáticas primero
rpmi_imdm <- error_filtrado_base %>%
  filter(medio_base %in% c("RPMI", "IMDM")) %>%
  arrange(desc(mae))

# Separamos RPMI e IMDM para analizarlos individualmente

# Distribución de tejidos en RPMI
# Queremos saber si RPMI tiene más hematológicos que otros medios
cat("=== TEJIDOS EN RPMI ===\n")
rpmi_cells <- rpmi_imdm[rpmi_imdm$medio_base == "RPMI", ]
print(sort(table(rpmi_cells$OncotreeLineage), decreasing = TRUE))

# Distribución de tejidos en IMDM
# IMDM es más específico — esperamos ver mieloides principalmente
cat("\n=== TEJIDOS EN IMDM ===\n")
imdm_cells <- rpmi_imdm[rpmi_imdm$medio_base == "IMDM", ]
print(sort(table(imdm_cells$OncotreeLineage), decreasing = TRUE))

# Distribución por tipo de cultivo en RPMI
# Confirmamos que RPMI tiene principalmente células en suspensión
cat("\n=== GROWTH PATTERN EN RPMI ===\n")
print(sort(table(rpmi_cells$GrowthPattern), decreasing = TRUE))

# Distribución por tipo de cultivo en IMDM
cat("\n=== GROWTH PATTERN EN IMDM ===\n")
print(sort(table(imdm_cells$GrowthPattern), decreasing = TRUE))

# Top 20 cell lines con mayor MAE en RPMI
# Cruzamos con metadata para ver el nombre, tejido y subtipo de cada una
# Esto nos permite identificar qué tumores específicos fallan más
cat("\n=== TOP 20 CELL LINES CON MAYOR MAE EN RPMI ===\n")
top_rpmi <- merge(
  rpmi_cells[rpmi_cells$medio_base == "RPMI",
             c("ModelID", "mae", "medio_base")],
  metadata[, c("ModelID", "CellLineName", "OncotreeLineage",
               "OncotreeSubtype", "GrowthPattern")],
  by = "ModelID", all.x = TRUE
)
top_rpmi <- top_rpmi[order(-top_rpmi$mae), ]  # ordenar por MAE descendente
print(head(top_rpmi[, c("CellLineName", "OncotreeLineage",
                        "OncotreeSubtype", "GrowthPattern",
                        "mae")], 20))

# Top 20 cell lines con mayor MAE en IMDM
# Igual que RPMI pero para el medio mieloide
cat("\n=== TOP 20 CELL LINES CON MAYOR MAE EN IMDM ===\n")
top_imdm <- merge(
  rpmi_imdm[rpmi_imdm$medio_base == "IMDM",
            c("ModelID", "mae", "medio_base")],
  metadata[, c("ModelID", "CellLineName", "OncotreeLineage",
               "OncotreeSubtype", "GrowthPattern")],
  by = "ModelID", all.x = TRUE
)
top_imdm <- top_imdm[order(-top_imdm$mae), ]
print(head(top_imdm[, c("CellLineName", "OncotreeLineage",
                        "OncotreeSubtype", "GrowthPattern",
                        "mae")], 20))

# ==============================================================
# PREGUNTA CLAVE: ¿el error en RPMI es por el medio o por
# el tipo tumoral?
# ==============================================================
# Si el problema fuera el medio RPMI, todos los que crecen en él
# fallarían igual. Si el problema es biológico, los hematológicos
# en RPMI fallarían MÁS que los tumores sólidos en RPMI.

# Etiquetamos cada cell line como hematológica o tumor sólido
rpmi_hema_vs_solido <- top_rpmi %>%
  mutate(tipo = ifelse(OncotreeLineage %in%
                         c("Lymphoid", "Myeloid"),
                       "Hematológico", "Tumor sólido"))

cat("=== MAE EN RPMI: HEMATOLÓGICO vs SÓLIDO ===\n")
rpmi_hema_vs_solido %>%
  group_by(tipo) %>%
  summarise(
    n           = n(),
    mae_medio   = mean(mae),
    mae_mediana = median(mae)
  ) %>%
  arrange(desc(mae_medio)) %>%
  print()

# Test de Wilcoxon: compara las distribuciones de MAE entre
# hematológicos y tumores sólidos dentro de RPMI
# Si p<0.05, la diferencia es real y el problema es biológico,
# no del medio de cultivo
wilcox.test(mae ~ tipo, data = rpmi_hema_vs_solido)


###Quiero comprobar esto en algo más especifico como IMDM

imdm_hema_vs_solido <- top_imdm %>%
  mutate(tipo = ifelse(OncotreeLineage %in%
                         c("Lymphoid", "Myeloid"),
                       "Hematológico", "Tumor sólido"))

cat("=== MAE EN IMDM: HEMATOLÓGICO vs SÓLIDO ===\n")
imdm_hema_vs_solido %>%
  group_by(tipo) %>%
  summarise(
    n           = n(),
    mae_medio   = mean(mae),
    mae_mediana = median(mae)
  ) %>%
  arrange(desc(mae_medio)) %>%
  print()

# Test de Wilcoxon: compara las distribuciones de MAE entre
# hematológicos y tumores sólidos dentro de IMDM
# Si p<0.05, la diferencia es real y el problema es biológico,
# no del medio de cultivo
wilcox.test(mae ~ tipo, data = imdm_hema_vs_solido)

####Vamos a mirar DMEM que es el otro medio de cultivo más importante
dmem <- error_filtrado_base %>%
  filter(medio_base %in% c("DMEM")) %>%
  arrange(desc(mae))
cat("=== TEJIDOS EN DMEM ===\n")
dmem_cells <- dmem[dmem$medio_base == "DMEM", ]
print(sort(table(dmem_cells$OncotreeLineage), decreasing = TRUE))

cat("\n=== TOP 20 CELL LINES CON MAYOR MAE EN DMEM ===\n")
top_dmem <- merge(
  dmem_cells[dmem_cells$medio_base == "DMEM",
             c("ModelID", "mae", "medio_base")],
  metadata[, c("ModelID", "CellLineName", "OncotreeLineage",
               "OncotreeSubtype", "GrowthPattern")],
  by = "ModelID", all.x = TRUE
)
top_dmem <- top_dmem[order(-top_dmem$mae), ]  # ordenar por MAE descendente
print(head(top_dmem[, c("CellLineName", "OncotreeLineage",
                        "OncotreeSubtype", "GrowthPattern",
                        "mae")], 20))

### Quiero mirar si hay ciertas enfermedades que me generan más problemas que otras
###Primero comprobamos en el medio de cultivo RPMI y despues en el de IMDM
rpmi_cells %>%
  group_by(OncotreePrimaryDisease) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  ) %>%
  filter(n >= 10) %>%
  arrange(desc(mae_medio))

imdm_cells %>%
  group_by(OncotreePrimaryDisease) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  ) %>%
  filter(n >= 3) %>%
  arrange(desc(mae_medio))

#### Solo poniendo las dos principales enfermedades en RPMI Y IMDM
rpmi_grupo <- rpmi_cells %>%
  mutate(grupo = ifelse(OncotreePrimaryDisease %in% 
                          c("Acute Myeloid Leukemia",
                            "Mature B-Cell Neoplasms"),
                        "Hematológicos clave",
                        "Resto"))

# Resumen
rpmi_grupo %>%
  group_by(grupo) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  )

# Test
wilcox.test(mae ~ grupo, data = rpmi_grupo)

#### Solo poniendo las dos principales enfermedades en RPMI Y IMDM
imdm_grupo <- imdm_cells %>%
  mutate(grupo = ifelse(OncotreePrimaryDisease %in% 
                          c("Acute Myeloid Leukemia",
                            "Mature B-Cell Neoplasms"),
                        "Hematológicos clave",
                        "Resto"))
# Resumen
imdm_grupo %>%
  group_by(grupo) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  )

# Test
wilcox.test(mae ~ grupo, data = imdm_grupo)
##Analizamos si el suero es un factor importante a la hora de analizar errores
error_by_media %>%
  group_by(SerumFreeMedia) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  )
### Mirar pediatrico porque puede ser que tengamos peor a chavales que adultos

error_by_media %>%
  group_by(PediatricModelType) %>%
  summarise(
    n = n(),
    
    mae_medio = mean(mae)
  )

wilcox.test(mae ~ PediatricModelType, data = error_by_media)


hema_data <- error_by_media %>%
  filter(OncotreeLineage %in% c("Lymphoid", "Myeloid"))

hema_data %>%
  group_by(PediatricModelType) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  )
wilcox.test(mae ~ PediatricModelType, data = hema_data)

# ==============================================================
# ANÁLISIS: GrowthPattern (tipo de crecimiento)
# ==============================================================

# Distribución
cat("\n=== DISTRIBUCIÓN GrowthPattern ===\n")
print(table(error_by_media$GrowthPattern, useNA = "ifany"))

# Resumen global
cat("\n=== MAE por GrowthPattern ===\n")
error_by_media %>%
  group_by(GrowthPattern) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  ) %>%
  print()

# Test global (más de 2 grupos)
cat("\n=== TEST KRUSKAL-WALLIS GrowthPattern ===\n")
print(kruskal.test(mae ~ GrowthPattern, data = error_by_media))

# Comparación clave: Suspension vs Adherent
subset_gp <- error_by_media %>%
  filter(GrowthPattern %in% c("Suspension", "Adherent"))

cat("\n=== MAE: Suspension vs Adherent ===\n")
subset_gp %>%
  group_by(GrowthPattern) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  ) %>%
  print()

cat("\n=== TEST Wilcoxon: Suspension vs Adherent ===\n")
print(wilcox.test(mae ~ GrowthPattern, data = subset_gp))


###¿El efecto hematológico es solo porque son suspensión?

#Suspension--> Hematologico vs solido

susp_data <- error_by_media %>%
  filter(GrowthPattern == "Suspension")

susp_data %>%
  group_by(OncotreeLineage %in% c("Lymphoid", "Myeloid")) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  )

susp_data <- susp_data %>%
  mutate(tipo = ifelse(OncotreeLineage %in% c("Lymphoid", "Myeloid"),
                       "Hematológico",
                       "No hematológico"))

wilcox.test(mae ~ tipo, data = susp_data)
###Adherent--> Hematológico vs sólido

adh_data <- error_by_media %>%
  filter(GrowthPattern == "Adherent")

adh_data %>%
  group_by(OncotreeLineage %in% c("Lymphoid", "Myeloid")) %>%
  summarise(
    n = n(),
    mae_medio = mean(mae)
  )
adh_data <- adh_data %>%
  mutate(tipo = ifelse(OncotreeLineage %in% c("Lymphoid", "Myeloid"),
                       "Hematológico",
                       "No hematológico"))

wilcox.test(mae ~ tipo, data = adh_data)

