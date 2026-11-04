# Por qué se muestran PLS, Ridge y SVR en las mezclas nuevas

La figura nominal frente a predicha resume **un modelo por objetivo**, elegido por el menor **MAE global observado en las mezclas nuevas** dentro de los tres candidatos de su familia. Se evaluaron nueve modelos en total: tres para NaCl, tres para KCl y tres para el cloruro total usando los modelos entrenados con PBS. Los modelos de distintos objetivos no compiten en un ranking único.

La selección fue **descriptiva y posterior a la evaluación de estas mezclas**, para resumir la figura. No fue una selección independiente fijada antes de observar el desempeño en este conjunto. No se reentrenaron modelos con las mezclas. En el código de las figuras quedó fijada en `SELECCIONADOS` como NaCl/PLS, PBS/Ridge y KCl/SVR; esa lista no se recalcula al dibujar la figura. Este cálculo comprueba que coincide con los mínimos de MAE global de los resultados individuales originales. Los otros seis modelos siguen incluidos en el heatmap y en las tablas.

## Datos y cálculo

La entrada es `figuras_poster/salidas/tablas/predicciones_mezclas.csv`, dentro del paquete [figuras_poster.zip](../figuras_poster.zip). Contiene 2.700 predicciones individuales, correspondientes a nueve modelos y **300 espectros por modelo**. Cada modelo predice las mismas cinco composiciones M1–M5, con **60 observaciones por composición**. No se interpretan como 300 preparaciones independientes.

Los objetivos nominales se verifican en cada fila:

- Familia NaCl: concentración nominal de NaCl, en mM.
- Familia KCl: concentración nominal de KCl, en mM.
- Familia PBS: cloruro total nominal aportado por NaCl y KCl, igual a la suma de sus concentraciones molares, en mM.

Se calcula el error individual `e_i = predicción_i − nominal_i`, con los valores sin redondear. Por modelo:

- `MAE = sum(abs(e_i)) / 300`.
- `RMSE = sqrt(sum(e_i**2) / 300)`.
- `Sesgo = sum(e_i) / 300`; positivo indica sobrepredicción.
- `R² = 1 − sum(e_i**2) / sum((nominal_i − media_nominal)**2)`.

Como cada mezcla tiene exactamente 60 observaciones, el MAE global también equivale a la media simple de los cinco MAE por mezcla: `MAE_global = (MAE_M1 + … + MAE_M5) / 5`. Esta equivalencia depende de los recuentos iguales. No se usa el error de la predicción media ni se promedian valores ya redondeados de la imagen. Tampoco se selecciona por número de mezclas ganadas.

## Comparación global

Todos los modelos tienen n = 300. Se redondea solo esta tabla para lectura.

| Objetivo | Modelo | Configuración | MAE (mM) | RMSE (mM) |
|---|---|---|---:|---:|
| NaCl | **PLS** | Optimizado | **9,073595** | 9,807609 |
| NaCl | Ridge | Optimizado | 9,141982 | 9,842305 |
| NaCl | OLS | Base | 9,636998 | 10,422726 |
| Cloruro, modelo PBS | **Ridge** | Base | **6,837534** | 8,206509 |
| Cloruro, modelo PBS | PLS | Base | 7,951394 | 10,254925 |
| Cloruro, modelo PBS | SVR (RBF) | Optimizado | 9,912564 | 10,863720 |
| KCl | **SVR (RBF)** | Base | **3,094386** | 3,917850 |
| KCl | GPR-PCA | Base | 16,485525 | 16,713626 |
| KCl | PLS | Base | 18,989948 | 20,284737 |

Ventaja respecto al segundo menor MAE de cada objetivo:

- NaCl: PLS reduce el MAE en **0,068387 mM**, aproximadamente **0,748 %**, respecto a Ridge.
- PBS/cloruro: Ridge reduce el MAE en **1,113860 mM**, aproximadamente **14,008 %**, respecto a PLS.
- KCl: SVR reduce el MAE en **13,391139 mM**, aproximadamente **81,230 %**, respecto a GPR-PCA.

Estos porcentajes de reducción comparan MAE entre modelos: `100 × (MAE_segundo − MAE_elegido) / MAE_segundo`. Son distintos de los puntos porcentuales de los heatmaps. El CSV global también incluye esos puntos porcentuales respecto a la madre: `MAE_pp = 100 × MAE_mM / madre_mM`, con 154 mM para NaCl, 10 mM para KCl y 140,7 mM de cloruro para PBS. Esa conversión no cambia el ranking dentro de cada familia.

PLS/NaCl tiene el menor MAE en cuatro de cinco mezclas, Ridge/PBS en dos de cinco y SVR/KCl en cinco de cinco. En PBS, PLS gana tres composiciones, pero las mejoras de Ridge en las otras dos hacen que Ridge tenga menor MAE global. La elección responde al promedio del error, no al número de victorias.

## Alcance de la justificación

La figura resume el mejor desempeño **observado entre los tres candidatos de cada objetivo en estas cinco mezclas**. La comparación completa de nueve modelos se mantiene en el mapa de calor. La ventaja de PLS sobre Ridge para NaCl es pequeña y no demuestra una diferencia estadística. En KCl, el menor MAE de SVR no implica buen ajuste absoluto: su R² es −2,689795 y su sesgo es +2,995201 mM. No se calcularon pruebas de significación a partir de barridos repetidos.

Si estas mezclas se presentan como un conjunto de prueba para un modelo fijado de antemano, no debe describirse esta elección posterior como una selección independiente previa. Puede describirse como selección para visualización, posterior a la evaluación, conservando la comparación completa.

Texto posible:

> Para resumir la visualización nominal frente a predicha, se presentó por cada objetivo el modelo con menor MAE global observado en los 300 espectros de las cinco mezclas. Esta selección fue descriptiva y posterior a la evaluación; los nueve modelos evaluados se conservaron en el mapa de calor y en las tablas de resultados.

## Archivos y reproducción

- [metricas_globales_nueve_modelos.csv](metricas_globales_nueve_modelos.csv): MAE, RMSE, sesgo, R², MAE en puntos porcentuales, modo, recuentos y ranking dentro de cada objetivo.
- [seleccion_por_objetivo.csv](seleccion_por_objetivo.csv): tres modelos mostrados, segundo modelo, diferencias de MAE y recuento de mezclas ganadas.
- [MAE_por_mezcla_nueve_modelos.csv](MAE_por_mezcla_nueve_modelos.csv): las 45 celdas calculadas con las predicciones individuales.
- [ganador_por_mezcla_y_objetivo.csv](ganador_por_mezcla_y_objetivo.csv): ganador específico de cada composición y objetivo.
- [concentraciones_nominales_mezclas.csv](concentraciones_nominales_mezclas.csv): clave M1–M5 con las concentraciones nominales de ambas sales.
- [procedencia_y_calculo.json](procedencia_y_calculo.json): archivo fuente, SHA-256 de sus bytes originales, fórmulas y comprobación de la lista fijada en el código.
- [calcular_metricas.py](calcular_metricas.py): código reproducible; requiere numpy y pandas del entorno ya proporcionado.

Desde la raíz del repositorio, con `figuras_poster.zip` presente:

```bash
python seleccion_modelos_mezclas/calcular_metricas.py
```

También permite `--paquete ruta/al/figuras_poster.zip --salida carpeta`. El script verifica los nominales, identidades, recuentos y correspondencia entre la lista de la figura y los mínimos calculados. No modifica el paquete, las figuras ni las predicciones.
