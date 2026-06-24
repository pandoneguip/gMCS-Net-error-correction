# gMCS-Net-error-correction
ANÁLISIS Y CORRECCIÓN ESTRUCTURAL DE ERRORES EN LA PREDICCIÓN DE VIABILIDAD CELULAR TRAS EL KNOCK-OUT DE GENES METABÓLICOS. PFG-Ingeniería biomédica-2026
# gMCS-Net Error Correction Pipeline

Pipeline de corrección de errores para el modelo gMCS-Net de predicción 
de viabilidad celular tras knockout génico.

## Descripción

Este repositorio contiene el código R desarrollado en el Trabajo de Fin 
de Grado (TFG) del Grado en Ingeniería Biomédica (Tecnun — Universidad 
de Navarra) para caracterizar y corregir los errores sistemáticos del 
modelo gMCS-Net sobre datos de DepMap (fold 2, esquema LELO).

## Estructura

- `01_median_polish.R` — Análisis exploratorio y descomposición del error 
  mediante Median Polish y componentes dispersos
- `02_analisis_medio_cultivo.R` — Análisis del sesgo por medio de cultivo 
  y linaje tumoral
- `03_pipeline_final.R` — Pipeline completo de corrección en tres capas: 
  Median Polish, pseudoinversa con transformación asinh, y corrección 
  por expresión génica
- `04_kfold_anidado.R` — Validación cruzada anidada de 5 iteraciones 
  sobre el conjunto de entrenamiento

## Datos

Los datos de DepMap utilizados (CRISPR.txt, Expression.txt, 
fold_2_data_LELO_results.txt) no se incluyen en este repositorio 
por su tamaño. Pueden descargarse desde https://depmap.org/portal/

## Requisitos

R >= 4.0.0

Paquetes: dplyr, tidyr, ggplot2, patchwork, org.Hs.eg.db, 
clusterProfiler, pheatmap

## Referencia

Proyecto de Fin de Grado — Ingeniería Biomédica  
Autor: Pablo Andonegui  
Tecnun — Universidad de Navarra, 2026  

