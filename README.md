# gMCS-Net Error Correction Pipeline

Pipeline de corrección de errores sistemáticos para el modelo gMCS-Net 
de predicción de viabilidad celular tras knockout génico.

## Descripción

Este repositorio contiene el código R desarrollado en el Proyecto de Fin 
de Grado (PFG) del Grado en Ingeniería Biomédica (Tecnun — Universidad 
de Navarra, 2026) para caracterizar y corregir los errores sistemáticos 
del modelo gMCS-Net sobre datos de DepMap (fold 2, esquema LELO).

El trabajo demuestra que los errores de gMCS-Net no son aleatorios sino 
sistemáticos, con patrones biológicamente interpretables asociados 
principalmente a los linajes tumorales hematológicos (linfoides y mieloides). 
Se diseña y valida un pipeline modular de tres capas que reduce el RMSE 
un 57% y mejora la mediana de correlación por gen un 261% sobre el 
conjunto de test.

## Estructura del repositorio

- `01_median_polish.R` — Análisis exploratorio del error mediante Median 
  Polish y descomposición en componentes dispersos. Incluye la 
  identificación del sesgo sistemático por gen y por linaje tumoral, 
  la deflación recursiva de rango 1, el análisis por tipo de cultivo 
  y el enriquecimiento funcional (Gene Ontology).

- `02_analisis_medio_cultivo.R` — Análisis estadístico del sesgo por 
  medio de cultivo. Incluye los tests de Kruskal-Wallis y Dunn post-hoc, 
  la comparación dentro de cada medio entre líneas celulares hematológicas 
  y tumores sólidos mediante el test de Wilcoxon, y el análisis del 
  factor pediátrico como posible variable de confusión.

- `03_pipeline_final.R` — Pipeline completo de corrección de errores en 
  tres capas aplicadas secuencialmente: (1) Median Polish para eliminar 
  sesgos aditivos por gen, (2) pseudoinversa de Moore-Penrose con 
  transformación asinh para corregir patrones de interacción 
  gen × línea celular, y (3) corrección por expresión génica para 
  capturar la heterogeneidad molecular residual. Todos los hiperparámetros 
  se seleccionan en validación y el conjunto de test se evalúa una única 
  vez al final. Incluye el análisis post-corrección por linaje tumoral 
  y medio de cultivo.

- `04_kfold_anidado.R` — Validación cruzada anidada de 5 iteraciones 
  sobre el subconjunto de entrenamiento (658 líneas celulares). Permite 
  estimar la variabilidad de los resultados y confirmar que el pipeline 
  no está sobreajustado a la partición oficial. Incluye el análisis 
  out-of-fold de correlación por gen y por línea celular por linaje tumoral.

## Datos necesarios

Los datos utilizados en este trabajo no se incluyen en el repositorio. 
Deben obtenerse de las siguientes fuentes:

**Datos de DepMap** (disponibles en https://depmap.org/portal/):
- `CRISPR.txt` — matriz de viabilidad celular tras knockout génico 
  (puntuaciones Chronos), dimensiones 1103 líneas celulares × 386 genes 
  metabólicos
- `Expression.txt` — matriz de expresión génica (log2 TPM+1), mismas 
  dimensiones
- `genes_index.txt` — índice de los 386 genes metabólicos (identificadores 
  Ensembl)
- `cell_index.txt` — índice de las 1103 líneas celulares (identificadores 
  DepMap)
- `Model.csv` — metadatos de las líneas celulares (linaje tumoral, medio 
  de cultivo, tipo pediátrico, etc.)

**Resultados de gMCS-Net**:
- `fold_2_data_LELO_results.txt` — predicciones de gMCS-Net para el fold 2 
  del esquema LELO, con columnas: cell_line, gene, essentiality, pred, group. 
  Este archivo fue proporcionado por el grupo de investigación. Para 
  obtenerlo contactar con los directores del trabajo.

## Esquema de partición de datos

El trabajo utiliza el fold 2 del esquema LELO (Leave Every Lineage Out):
- **Entrenamiento:** 658 líneas celulares — aprendizaje de todos los 
  parámetros del pipeline
- **Validación:** 222 líneas celulares — selección de hiperparámetros
- **Test:** 223 líneas celulares — evaluación final, evaluado una única vez

## Requisitos

**R** >= 4.0.0

**Paquetes CRAN:**
```r
install.packages(c("dplyr", "tidyr", "tibble", "ggplot2", 
                   "patchwork", "pheatmap", "readr", "conflicted"))
```

**Paquetes Bioconductor:**
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("org.Hs.eg.db", "clusterProfiler"))
```

## Resultados principales

| Métrica | Baseline | Pipeline final | Mejora |
|---|---|---|---|
| RMSE test | 0.2628 | 0.1131 | −57.0% |
| Error cuadrático test | 5945.43 | 1101.74 | −81.5% |
| Mediana correlación por gen | 0.0530 | 0.1915 | +261% |
| Media correlación por línea celular | 0.9144 | 0.9850 | +7.7pp |

La brecha de error entre líneas celulares hematológicas y tumores sólidos 
se reduce un 41.9% (fold 2) y un 72.5% (fold 3), confirmando que el 
pipeline corrige específicamente el sesgo identificado en el análisis 
exploratorio.

## Referencia

Trabajo de Fin de Grado — Ingeniería Biomédica  
**Autor:** Pablo Andonegui  
**Institución:** Tecnun — Universidad de Navarra, 2026  
