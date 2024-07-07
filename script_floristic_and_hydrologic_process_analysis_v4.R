
## pacotes
library(vegan)
library(factoextra) 
library(ggpubr)
library(ggplot2)
library(dplyr)
library(PMCMRplus)
library(PerformanceAnalytics)
library(gsubfn)
library(GGally)
library(ggrepel)
library(fields)
library(spdep)
library(openxlsx)
library(FD)
library(reshape2)
library(tidyr)

##############################
############
diretorio <- "C:/Users/pedro.rajao/Desktop/doutorado UFRJ/9. Cap4_CWM-FD role on rainfall interception canopy-litter/1. final_data"
arquivos_txt <- list.files(path = diretorio, pattern = "\\.txt$", full.names = TRUE)

carregar_arquivo <- function(caminho) {
  nome <- tools::file_path_sans_ext(basename(caminho))
  assign(nome, read.table(caminho, header = TRUE, sep = "\t"), envir = .GlobalEnv)
}

lapply(arquivos_txt, carregar_arquivo)



###########
#### process

result <- merge(df_hidroerosion, df_rainfall, by.x = "ID_evento_chuva", by.y = "ID_chuva_INMET")
clean_result1 <- cbind(result[1:17], result[21:28])
result2 <- merge(clean_result1, df_structure, by.x = "etiqueta", by.y = "etiqueta")
clean_result2 <- cbind(result2[1:25], result2[28:38])
result3 <- merge(clean_result2, df_slope, by.x = "etiqueta", by.y = "ID")
clean_result3 <- cbind(result3[1:31], result3[33:36], result3[49:50])
result4 <- merge(clean_result3, df_litter, by.x = "etiqueta", by.y = "Arquivo_ID")
colunas <- names(result4)
cols_numericas <- c("turbidez", "turbidez_max")
result4$turbidez <- gsub(",", ".", result4$turbidez)
result4$turbidez <- gsub("\\.", ".", result4$turbidez, fixed = TRUE)
result4$turbidez <- as.numeric(result4$turbidez)

result4 <- result4 %>%
  mutate(
    across(c(transprecipitacao_mm, escoamento_mm, coef_trans, coef_inter, coef_esc, coef_esc_trans,
             preciptacao_mm_INMET, turbidez_max, chuva_acumulada_mm, intensidade_max_mmh,
             velocidade_vento, umidade_do_ar_media, temperatura_media, abertura, caco, 
             index_caco, gap_fraction, ARVI, LAI, NDVI, media, sd),
           ~ as.numeric(gsub(",", ".", .)))
  )

# Criar o novo dataframe com médias e CVs
novo_dfA <- result4 %>%
  group_by(etiqueta, ID_evento_chuva) %>%
  summarize(
    turbidez = mean(turbidez, na.rm = TRUE),
    media_coef_trans = mean(coef_trans, na.rm = TRUE) * 100,
    cv_coef_trans = sd(coef_trans, na.rm = TRUE) / mean(coef_trans, na.rm = TRUE) * 100,
    media_coef_inter = mean(coef_inter, na.rm = TRUE) * 100,
    cv_coef_inter = sd(coef_inter, na.rm = TRUE) / mean(coef_inter, na.rm = TRUE) * 100,
    media_coef_esc = mean(coef_esc, na.rm = TRUE) * 100,
    cv_coef_esc = sd(coef_esc, na.rm = TRUE) / mean(coef_esc, na.rm = TRUE) * 100,
    media_coef_esc_trans = mean(coef_esc_trans, na.rm = TRUE) * 100,
    cv_coef_esc_trans = sd(coef_esc_trans, na.rm = TRUE) / mean(coef_esc_trans, na.rm = TRUE) * 100,
    media_dias_sem_chuva1mm = mean(dias_sem_chuva1mm, na.rm = TRUE),
    media_dias_sem_chuva02mm = mean(dias_sem_chuva02mm, na.rm = TRUE),
    media_velocidade_vento = mean(velocidade_vento, na.rm = TRUE),
    media_umidade_do_ar_media = mean(umidade_do_ar_media, na.rm = TRUE),
    media_temperatura_media = mean(temperatura_media, na.rm = TRUE),
    media_abertura = mean(abertura, na.rm = TRUE),
    media_caco = mean(caco, na.rm = TRUE),
    media_index_caco = mean(index_caco, na.rm = TRUE),
    media_gap_fraction = mean(gap_fraction, na.rm = TRUE),
    media_ndvi = mean(NDVI, na.rm = TRUE),
    media_arvi = mean(ARVI, na.rm = TRUE),
    media_lai = mean(LAI, na.rm = TRUE),
    media_riqueza = mean(riqueza, na.rm = TRUE),
    media_densidade = mean(densidade, na.rm = TRUE),
    media_media = mean(media, na.rm = TRUE),
    media_sd = mean(sd, na.rm = TRUE),
    media_chuva_acumulada_mm = mean(chuva_acumulada_mm, na.rm = TRUE),
    media_intensidade_max_mmh = mean(intensidade_max_mmh, na.rm = TRUE),
    litter_cover = mean(litter_cover, na.rm = TRUE),
    number_leaf_liter = mean(leafs_numb, na.rm = TRUE),
    CWM_area_cm_quadrados = mean(CWM_area_cm_quadrados, na.rm = TRUE),
    FDrich_area_cm_quadrados = mean(FDrich_area_cm_quadrados, na.rm = TRUE),
    FDdiv_area_cm_quadrados = mean(FDdiv_area_cm_quadrados, na.rm = TRUE),
    CWM_perimeter = mean(CWM_perimeter, na.rm = TRUE),
    FDrich_perimeter = mean(FDrich_perimeter, na.rm = TRUE),
    FDdiv_perimeter = mean(FDdiv_perimeter, na.rm = TRUE),
    CWM_perimeter_area_ratio = mean(CWM_perimeter_area_ratio, na.rm = TRUE),
    FDrich_perimeter_area_ratio = mean(FDrich_perimeter_area_ratio, na.rm = TRUE),
    FDdiv_perimeter_area_ratio = mean(FDdiv_perimeter_area_ratio, na.rm = TRUE),
    .groups = 'drop'
  )


df_coord <- df_coord %>%
  rename(Name = Name.xcoord.ycoord.zcoord)
df_coord <- df_coord %>%
  separate(Name, into = c("etiqueta", "xcoord", "ycoord", "zcoord"), sep = ",")

novo_df <- merge(novo_dfA, df_coord, by.x = "etiqueta", by.y = "etiqueta")


######
#### 
# exploratory: process

#### Throughfall
##
novo_df_throughfall <- novo_df[complete.cases(novo_df$media_coef_trans),]
mean_values <- tapply(novo_df_throughfall$media_coef_trans, novo_df_throughfall$etiqueta, mean)
ordered_levels <- names(sort(mean_values))
cores_etiquetas <- colorRampPalette(c("lightgreen", "darkgreen"))(length(ordered_levels))
highest_mean <- max(mean_values, na.rm = TRUE)
lowest_mean <- min(mean_values, na.rm = TRUE)

# Create the boxplot
boxplot_trans <- ggplot(novo_df_throughfall, aes(x = factor(etiqueta, levels = ordered_levels), y = media_coef_trans, fill = factor(etiqueta, levels = ordered_levels))) +
  geom_boxplot() +
  labs(x = "runoff plots", y = "coef. Throughfall (%)") +
  scale_fill_manual(values = cores_etiquetas) +
  theme_test() +
  theme(legend.text = element_text(size = 8),
        title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(0, max(novo_df_throughfall$media_coef_trans, na.rm = TRUE), by = 20)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = highest_mean, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = lowest_mean, linetype = "dashed", color = "lightgreen") +
  stat_compare_means(method = "anova", label = "p.format", size = 3, label.x = 2.5, label.y = 162)

print(boxplot_trans)
ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/boxplot_throughfall_completo.png", plot = boxplot_trans, width = 10, height = 6, units = "in", dpi = 300)


## boxplots by Evento de chuva
create_boxplots <- function(data, event_id) {
  # Filtrar dados para o ID_evento_chuva específico
  subset_data <- subset(data, ID_evento_chuva == event_id)
  
  # Calcular média por etiqueta
  mean_values <- tapply(subset_data$coef_trans, subset_data$etiqueta, mean, na.rm = TRUE)
  ordered_levels <- names(sort(mean_values))
  
  # Criar o boxplot
  boxplot_trans <- ggplot(subset_data, aes(x = factor(etiqueta, levels = ordered_levels), y = coef_trans*100, fill = factor(etiqueta, levels = ordered_levels))) +
    geom_boxplot() +
    labs(x = "runoff plots", y = "coef. Throughfall (%)") +
    scale_fill_manual(values = colorRampPalette(c("lightgreen", "darkgreen"))(length(ordered_levels))) +
    theme_test() +
    theme(
      legend.text = element_text(size = 8),
      title = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 45, hjust = 1),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = seq(0, max(subset_data$coef_trans*100, na.rm = TRUE), by=20)) +
    geom_hline(yintercept = max(mean_values*100, na.rm = TRUE), linetype = "dashed", color = "darkgreen") +
    geom_hline(yintercept = min(mean_values*100, na.rm = TRUE), linetype = "dashed", color = "lightgreen") +
    ggtitle(paste("Evento de chuva", event_id))
    
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Trans_evento_", event_id, ".png", sep = "")
  print(boxplot_trans)
  ggsave(filename, plot = boxplot_trans, width = 10, height = 6, units = "in", dpi = 300)
  
  return(boxplot_trans)
}

# Obter IDs únicos de eventos de chuva
unique_event_ids <- unique(result4$ID_evento_chuva)
for (event_id in unique_event_ids) {
  boxplot_result <- create_boxplots(result4, event_id)
  print(boxplot_result)
}


#### Runoff
##
novo_df_runoff <- novo_df_throughfall[complete.cases(novo_df_throughfall$media_coef_esc),]
mean_values_runoff2 <- tapply(novo_df_runoff$media_coef_esc, novo_df_runoff$etiqueta, mean)
#ordered_levels <- names(sort(mean_values_runoff2))
cores_etiquetas <- colorRampPalette(c("lightblue", "darkblue"))(length(ordered_levels))

highest_mean <- max(mean_values_runoff2, na.rm = TRUE)
lowest_mean <- min(mean_values_runoff2, na.rm = TRUE)

# Create the boxplot
boxplot_esc <- ggplot(novo_df_runoff, aes(x = factor(etiqueta, levels = ordered_levels), y = media_coef_esc, fill = factor(etiqueta, levels = ordered_levels))) +
  geom_boxplot() +
  labs(x = "runoff plots", y = "coef. Runoff (%)") +
  scale_fill_manual(values = cores_etiquetas) +
  theme_test() +
  theme(legend.text = element_text(size = 8),
        title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none") +
  scale_y_continuous(breaks = seq(0, max(novo_df1$media_coef_esc, na.rm = TRUE), by = 10)) +
  geom_hline(yintercept = highest_mean, linetype = "dashed", color = "darkblue") +
  geom_hline(yintercept = lowest_mean, linetype = "dashed", color = "lightblue") +
  coord_cartesian(ylim = c(0, 50)) +
  stat_compare_means(method = "anova", label = "p.format", size = 3, label.x = 2.1, label.y = 49.5)  # Replace Your_F_Value with the actual F-value

print(boxplot_esc)
ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/boxplot_runoff_completo.png", plot = boxplot_esc, width = 10, height = 6, units = "in", dpi = 300)

## boxplots by Evento de chuva
create_point_plots <- function(data, event_id) {
  # Filtrar dados para o ID_evento_chuva específico
  subset_data <- subset(data, ID_evento_chuva == event_id)
  
  # Calcular média por etiqueta
  mean_values <- tapply(subset_data$coef_esc, subset_data$etiqueta, mean, na.rm = TRUE)
  ordered_levels <- names(sort(mean_values))
  
  # Criar o gráfico de pontos
  point_plot_esc <- ggplot(subset_data, aes(x = factor(etiqueta, levels = ordered_levels), y = coef_esc * 100, fill = factor(etiqueta, levels = ordered_levels))) +
    geom_point(shape = 21, size = 3, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
    labs(x = "runoff plots", y = "coef. Runoff (%)") +
    scale_fill_manual(values = colorRampPalette(c("lightblue", "darkblue"))(length(ordered_levels))) +
    theme_test() +
    theme(
      legend.text = element_text(size = 8),
      title = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 45, hjust = 1),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
    ) +
    scale_y_continuous(breaks = seq(0, max(subset_data$coef_esc * 100, na.rm = TRUE), by = 2)) +
    geom_hline(yintercept = max(mean_values * 100, na.rm = TRUE), linetype = "dashed", color = "darkblue") +
    geom_hline(yintercept = min(mean_values * 100, na.rm = TRUE), linetype = "dashed", color = "lightblue") +
    ggtitle(paste("Evento de chuva", event_id))
  
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Esc_evento_", event_id, ".png", sep = "")
  print(point_plot_esc)
  ggsave(filename, plot = point_plot_esc, width = 10, height = 6, units = "in", dpi = 300)
  
  return(point_plot_esc)
}

# Obter IDs únicos de eventos de chuva
unique_event_ids <- unique(result4$ID_evento_chuva)

# Criar e exibir gráficos de pontos para cada ID_evento_chuva
for (event_id in unique_event_ids) {
  point_plot_result <- create_point_plots(result4, event_id)
  print(point_plot_result)
}


## Turbity
novo_df_tur <- novo_df_throughfall[complete.cases(novo_df_throughfall$turbidez),]
mean_values_tur <- tapply(novo_df_tur$turbidez, novo_df_tur$etiqueta, mean)
ordered_levels <- names(sort(mean_values_tur))
cores_etiquetas <- colorRampPalette(c("#ffb48a", "#5d2417"))(length(ordered_levels))
highest_mean <- max(mean_values_tur, na.rm = TRUE)
lowest_mean <- min(mean_values_tur, na.rm = TRUE)

# Create the boxplot
boxplot_tur <- ggplot(novo_df_tur, aes(x = factor(etiqueta, levels = ordered_levels), y = turbidez, fill = factor(etiqueta, levels = ordered_levels))) +
  geom_boxplot() +
  labs(x = "runoff plots", y = "turbity (NTU)") +
  scale_fill_manual(values = cores_etiquetas) +
  theme_test() +
  theme(legend.text = element_text(size = 8),
        title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none") +
  geom_hline(yintercept = highest_mean, linetype = "dashed", color = "#5d2417") +
  geom_hline(yintercept = lowest_mean, linetype = "dashed", color = "#ffb48a") +
  stat_compare_means(method = "anova", label = "p.format", size = 3, label.x = 2.2, label.y = 3260)

print(boxplot_tur)
ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/boxplot_turbidez_completo.png", plot = boxplot_tur, width = 10, height = 6, units = "in", dpi = 300)



### regressão por evento de chuva

library(ggplot2)
library(ggpubr)

create_regression_plots <- function(data, event_id) {
  # Filtrar dados para o ID_evento_chuva específico
  subset_data <- subset(data, ID_evento_chuva == event_id)
  
  # Obter valores da regressão
  regression_summary <- summary(lm(coef_esc ~ coef_trans, data = subset_data))
  p_value <- signif(regression_summary$coef[, "Pr(>|t|)"], digits = 3)
  r_squared <- signif(regression_summary$adj.r.squared, digits = 3)
  
  # Criar o gráfico de dispersão com linha de regressão
  regression_plot <- ggplot(subset_data, aes(x = coef_trans, y = coef_esc)) +
    geom_point(shape = 21, fill = "darkblue", color = "black", size = 3) +
    geom_smooth(method = "lm", se = TRUE, col = "black", fill = "lightgray", level = 0.95) +
    labs(x = "coef. Throughfall", y = "coef. Runoff") +
    theme_test() +
    theme(
      legend.text = element_text(size = 8),
      title = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 45, hjust = 1),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
    ) +
    ggtitle(paste("Regressão - Evento de chuva", event_id)) +
    annotate("text", x = min(subset_data$coef_trans, na.rm = TRUE), y = max(subset_data$coef_esc, na.rm = TRUE), 
             label = paste("p =", p_value), vjust = 1, hjust = 0) +
    annotate("text", x = min(subset_data$coef_trans, na.rm = TRUE), y = max(subset_data$coef_esc, na.rm = TRUE) - 0.2, 
             label = paste("R²-adjusted =", r_squared), vjust = 1, hjust = 0)
  
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Regressao_evento_", event_id, ".png", sep = "")
  print(regression_plot)
  ggsave(filename, plot = regression_plot, width = 10, height = 6, units = "in", dpi = 300)
  
  return(regression_plot)
}

# Obter IDs únicos de eventos de chuva
unique_event_ids <- unique(result4$ID_evento_chuva)

# Criar e exibir gráficos de regressão para cada ID_evento_chuva
for (event_id in unique_event_ids) {
  regression_plot_result <- create_regression_plots(result4, event_id)
  print(regression_plot_result)
}


create_regression_plots <- function(data, event_id) {
  # Filtrar dados para o ID_evento_chuva específico
  subset_data <- subset(data, ID_evento_chuva == event_id)
  
  # Obter valores da regressão
  regression_summary <- summary(lm(coef_esc ~ coef_trans, data = subset_data))
  p_value <- signif(regression_summary$coef[, "Pr(>|t|)"], digits = 3)
  r_squared <- signif(regression_summary$adj.r.squared, digits = 3)
  
  # Criar o gráfico de dispersão com linha de regressão
  regression_plot <- ggplot(subset_data, aes(x = coef_trans, y = coef_esc)) +
    geom_point(shape = 21, fill = "darkblue", color = "black", size = 3) +
    geom_smooth(method = "lm", se = TRUE, col = "black", fill = "lightgray", level = 0.95) +
    labs(x = "coef. Throughfall", y = "coef. Runoff") +
    theme_test() +
    theme(
      legend.text = element_text(size = 8),
      title = element_text(size = 8),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 45, hjust = 1),
      axis.ticks.length = unit(0.15, "cm"),
      legend.position = "none"
    ) +
    ggtitle(paste("Regressão - Evento de chuva", event_id)) +
    annotate("text", x = min(subset_data$coef_trans, na.rm = TRUE), y = max(subset_data$coef_esc, na.rm = TRUE), 
             label = paste("p =", p_value), vjust = 1, hjust = 0) +
    annotate("text", x = min(subset_data$coef_trans, na.rm = TRUE), y = max(subset_data$coef_esc, na.rm = TRUE) - 0.2, 
             label = paste("R²-adjusted =", r_squared), vjust = 1, hjust = 0)
  
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Regressao_evento_", event_id, ".png", sep = "")
  print(regression_plot)
  ggsave(filename, plot = regression_plot, width = 10, height = 6, units = "in", dpi = 300)
  
  return(regression_plot)
}

# Obter IDs únicos de eventos de chuva
unique_event_ids <- unique(result4$ID_evento_chuva)

# Criar e exibir gráficos de regressão para cada ID_evento_chuva
for (event_id in unique_event_ids) {
  regression_plot_result <- create_regression_plots(result4, event_id)
  print(regression_plot_result)
}

#### plots
## Correlations
novo_df4 <- novo_df[-1]
cor_matrix2 <- cor(novo_df4, method='pearson')
cor_df <- as.data.frame(as.table(cor_matrix2))
names(cor_df) <- c("Var1", "Var2", "value")

COR <- ggplot(data = cor_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  #geom_point(shape = 21, size = 4, color = "white", stroke = 0.5) +
  scale_fill_gradient2(low = "darkred", high = "darkblue", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", name = "Correlation",
                       na.value = "white") +  # Defina a cor para valores NA como branco
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, color = "black"),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )
ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Correlations_ret.png", plot = COR, width = 10, height = 6, units = "in", dpi = 300)


ggpairs(novo_df4, lower = list(continuous = "smooth"))
ggpairs(novo_df4, ggplot2::aes(colour=as.character(ID_evento_chuva)))
ggpairs(novo_df4, lower = list(continuous = "smooth"))
chart.Correlation(novo_df4, method='pearson')





################
########## relações floristicas
library(tidyverse)
library(dplyr)

# Remove columns 10 to 25
df_data_floristic2 <- df_data_floristic[, -c(10:25)]
df_data_floristic2 <- df_data_floristic2 %>%
  mutate(
    cap = as.numeric(gsub(",", ".", cap)),
    dbh = as.numeric(gsub(",", ".", dbh))
  )

# Calculate area_basal
df_data_floristic2 <- df_data_floristic2 %>%
  mutate(area_basal = (dbh / 2)^2 * pi)

# Calculate absolute density, dominance, and frequency
resultados_fito <- df_data_floristic2 %>%
  group_by(sp) %>%
  summarise(
    densidade_absoluta = sum(n) / n_distinct(etiqueta),
    dominancia_absoluta = sum(area_basal),
    frequencia_absoluta = n_distinct(etiqueta)
  ) %>%
  filter(row_number() > 1)  # Remove the first row if needed

# Calculate relative values
resultados_fito <- resultados_fito %>%
  mutate(
    densidade_absoluta_relativa = 100 * (densidade_absoluta / sum(densidade_absoluta)),
    dominancia_absoluta_relativa = 100 * (dominancia_absoluta / sum(dominancia_absoluta)),
    frequencia_absoluta_relativa = 100 * (frequencia_absoluta / sum(frequencia_absoluta)),
    IVC = (dominancia_absoluta_relativa + densidade_absoluta_relativa) / 2,
    IVI = (dominancia_absoluta_relativa + densidade_absoluta_relativa + frequencia_absoluta_relativa) / 3
  )

# Arrange results by descending IVI
resultados_fito <- resultados_fito %>%
  arrange(desc(IVI))

# View the resulting dataframe
print(resultados_fito)

diretorio <- "C:/Users/pedro.rajao/Desktop/doutorado UFRJ/7. Cap2_leaf traits on storage and drainage, canopy and litter (Texto Qualificação)/1. final_data/resultados_fito_final.txt"
write.table(resultados_fito, diretorio, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


############
######### 
# Diversidade Funcional

#df_completoA <- merge(df_data_floristic2, df_leaf_traits, by = "sp")
df_completo <- left_join(df_data_floristic2, df_leaf_traits, by = "sp")
df_completo <- df_completo %>%
  mutate_all(~as.numeric(gsub(",", ".", .)))

# Assuming your dataframe is named df
df_completo$SLA <- df_completo$la_sla / df_completo$peso_seco_sla
df_completo$LDMC <- df_completo$peso_seco_sla / df_completo$peso_0_sla
df_completo$WHC <- (df_completo$peso_1_whc_escorrido - df_completo$peso_0_whc) / df_completo$peso_0_whc
df_completo$WHCmax <- (df_completo$peso_1_whc_seco.papel - df_completo$peso_0_whc) / df_completo$peso_0_whc

###### LA
cwm_vars <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("la_whc"), weighted.mean, w = n, na.rm = TRUE))

fd <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_la_whc = sd(la_whc, na.rm = TRUE)
  )

resultados_funcional <- left_join(cwm_vars, fd, by = "etiqueta")

FD_CWM <- ggplot(resultados_funcional, aes(x = la_whc, y = SD_la_whc)) +
  geom_point() +
  labs(title = "Relação entre FD_la_whc e CWM_la_whc",
       y = "FD_divergence SD la_whc",
       x = "CWM_la_whc") +
  theme_classic() +
  # Definindo os limites para os eixos x e y
  scale_x_continuous(limits = c(0, 300)) +
  scale_y_continuous(limits = c(0, 150))
print(FD_CWM)

resultados_funcional <- resultados_funcional %>%
  mutate(
    grupo_LA = case_when(
      la_whc <= 120 & SD_la_whc <= 50 ~ "LL",
      la_whc <= 120 & SD_la_whc > 50 ~ "LH",
      la_whc > 120 & SD_la_whc > 50 ~ "HH",
      la_whc > 120 & SD_la_whc <= 50 ~ "HL"
    )
  )

FD_CWM_color <- ggplot(resultados_funcional, aes(x = la_whc, y = SD_la_whc, color = grupo_LA, label = etiqueta)) +
  geom_point() +
  geom_text(aes(color = grupo_LA), vjust = -0.5, hjust = 0.5, size=2.5)  +
  labs(title = "Relação entre FD_la_whc e CWM_la_whc",
       y = "FD_la_whc",
       x = "CWM_la_whc") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 300)) +
  scale_y_continuous(limits = c(0, 150)) +
  scale_color_manual(values = c("LL" = "red", "LH" = "orange", "HH" = "darkgreen", "HL" = "darkblue"))

print(FD_CWM_color)

ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap.3/FD_CWM_LA.png", plot = FD_CWM_color, width = 10, height = 6, units = "in", dpi = 300)


###### CAP
cwm_vars_CAP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("cap"), weighted.mean, w = n, na.rm = TRUE))

fd_CAP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_CAP = sd(cap, na.rm = TRUE)
  )
resultados_funcional_CAP <- left_join(cwm_vars_CAP, fd_CAP, by = "etiqueta")

###### H
cwm_vars_H <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("h"), weighted.mean, w = n, na.rm = TRUE))

fd_H <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_H = sd(h, na.rm = TRUE)
  )
resultados_funcional_H <- left_join(cwm_vars_H, fd_H, by = "etiqueta")


###### SLA
cwm_vars_SLA <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("SLA"), weighted.mean, w = n, na.rm = TRUE))

fd_SLA <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_SLA = sd(SLA, na.rm = TRUE)
  )
resultados_funcional_SLA <- left_join(cwm_vars_SLA, fd_SLA, by = "etiqueta")

FD_CWM_SLA <- ggplot(resultados_funcional_SLA, aes(x = SLA, y = SD_SLA)) +
  geom_point() +
  labs(title = "Relação entre FD_SLA e CWM_SLA",
       y = "FD divergence SD SLA",
       x = "CWM_SLA") +
  theme_classic() +
  # Definindo os limites para os eixos x e y
  scale_x_continuous(limits = c(50, 400)) +
  scale_y_continuous(limits = c(0, 200))
print(FD_CWM_SLA)

resultados_funcional_SLA <- resultados_funcional_SLA %>%
  mutate(
    grupo_SLA = case_when(
      SLA <= 200 & SD_SLA <= 80 ~ "LL",
      SLA <= 200 & SD_SLA > 80 ~ "LH",
      SLA > 200 & SD_SLA > 80 ~ "HH",
      SLA > 200 & SD_SLA <= 80 ~ "HL"
    )
  )

FD_CWM_color_SLA <- ggplot(resultados_funcional_SLA, aes(x = SLA, y = SD_SLA, color = grupo_SLA, label = etiqueta)) +
  geom_point() +
  geom_text(aes(color = grupo_SLA), vjust = -0.5, hjust = 0.5, size=2.5) +
  labs(title = "Relação entre FD_SLA e CWM_SLA",
       y = "FD divergence SD SLA",
       x = "CWM_SLA") +
  theme_classic() +
  scale_x_continuous(limits = c(50, 400)) +
  scale_y_continuous(limits = c(0, 150)) +
  scale_color_manual(values = c("LL" = "red", "LH" = "orange", "HH" = "darkgreen", "HL" = "darkblue"))

print(FD_CWM_color_SLA)

ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap.3/FD_CWM_SLA.png", plot = FD_CWM_color_SLA, width = 10, height = 6, units = "in", dpi = 300)

###### LDMC
cwm_vars_LDMC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("LDMC"), weighted.mean, w = n, na.rm = TRUE))

fd_LDMC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_LDMC = sd(LDMC, na.rm = TRUE)
  )
resultados_funcional_LDMC <- left_join(cwm_vars_LDMC, fd_LDMC, by = "etiqueta")

FD_CWM_LDMC <- ggplot(resultados_funcional_LDMC, aes(x = LDMC, y = SD_LDMC)) +
  geom_point() +
  labs(title = "Relação entre FD_LDMC e CWM_LDMC",
       y = "FD_LDMC",
       x = "CWM_LDMC") +
  theme_classic() 
print(FD_CWM_LDMC)

resultados_funcional_LDMC <- resultados_funcional_LDMC %>%
  mutate(
    grupo_LDMC = case_when(
      LDMC <= 0.47 & SD_LDMC <= 0.13 ~ "LL",
      LDMC <= 0.47 & SD_LDMC > 0.13 ~ "LH",
      LDMC > 0.47 & SD_LDMC > 0.13 ~ "HH",
      LDMC > 0.47 & SD_LDMC <= 0.13 ~ "HL"
    )
  )

FD_CWM_color_LDMC <- ggplot(resultados_funcional_LDMC, aes(x = LDMC, y = SD_LDMC, color = grupo_LDMC, label = etiqueta)) +
  geom_point() +
  geom_text(aes(color = grupo_LDMC), vjust = -0.5, hjust = 0.5, size=2.5) +
  labs(title = "Relação entre FD_LDMC e CWM_LDMC",
       y = "FD divergence SD LDMC",
       x = "CWM_LDMC") +
  theme_classic() +
  scale_color_manual(values = c("LL" = "red", "LH" = "orange", "HH" = "darkgreen", "HL" = "darkblue"))

print(FD_CWM_color_LDMC)

ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap.3/FD_CWM_LDMC.png", plot = FD_CWM_color_LDMC, width = 10, height = 6, units = "in", dpi = 300)

###### WHC
cwm_vars_WHC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("WHC"), weighted.mean, w = n, na.rm = TRUE))

fd_WHC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_WHC = sd(WHC, na.rm = TRUE)
  )
resultados_funcional_WHC <- left_join(cwm_vars_WHC, fd_WHC, by = "etiqueta")

FD_CWM_WHC <- ggplot(resultados_funcional_WHC, aes(x = WHC, y = SD_WHC)) +
  geom_point() +
  labs(title = "Relação entre FD_WHC e CWM_WHC",
       y = "FD_WHC",
       x = "CWM_WHC") +
  theme_classic() 
print(FD_CWM_WHC)


###### WHCmax
cwm_vars_WHCmax <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("WHCmax"), weighted.mean, w = n, na.rm = TRUE))

fd_WHCmax <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_WHCmax = sd(WHCmax, na.rm = TRUE)
  )
resultados_funcional_WHCmax <- left_join(cwm_vars_WHCmax, fd_WHCmax, by = "etiqueta")


FD_CWM_color_WHCmax <- ggplot(resultados_funcional_WHC, aes(x = WHCmax, y = SD_WHCmax, color = grupo_WHCmax, label = etiqueta)) +
  geom_point() +
  geom_text(aes(color = grupo_WHCmax), vjust = -0.5, hjust = 0.5, size=2.5) +
  labs(title = "Relação entre FD_WHCmax e CWM_WHCmax",
       y = "FD_WHCmax",
       x = "CWM_WHCmax") +
  theme_classic() +
  scale_color_manual(values = c("LL" = "red", "LH" = "orange", "HH" = "darkgreen", "HL" = "darkblue"))

print(FD_CWM_color_WHCmax)

ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap.3/FD_CWM_WHC.png", plot = FD_CWM_color_WHCmax, width = 10, height = 6, units = "in", dpi = 300)



novo_df$etiqueta <- as.character(novo_df$etiqueta)
resultados_funcional_LDMC$etiqueta <- as.character(resultados_funcional_LDMC$etiqueta)
merged_df1_2 <- merge(novo_df, resultados_funcional_LDMC, by = "etiqueta", all = TRUE)
merged_df1_2_3 <- merge(merged_df1_2, resultados_funcional_SLA, by = "etiqueta", all = TRUE)
merged_df1_2_3_4 <- merge(merged_df1_2_3, resultados_funcional_WHC, by = "etiqueta", all = TRUE)
merged_df1_2_3_4_5 <- merge(merged_df1_2_3_4, resultados_funcional_CAP, by = "etiqueta", all = TRUE)
merged_df1_2_3_4_5_6 <- merge(merged_df1_2_3_4_5, resultados_funcional_H, by = "etiqueta", all = TRUE)
merged_df1_2_3_4_5_6_7 <- merge(merged_df1_2_3_4_5_6, resultados_funcional_WHCmax, by = "etiqueta", all = TRUE)

final_merged_df <- merge(merged_df1_2_3_4_5_6_7, resultados_funcional, by = "etiqueta", all = TRUE)



View(final_merged_df)

final_merged_dfA <- final_merged_df %>%
  group_by(etiqueta) %>%
  summarize(
    runoff = mean(media_coef_esc, na.rm=TRUE),
    trans = mean(media_coef_trans, na.rm = TRUE),
    SD_trans = sd(media_coef_trans, na.rm = TRUE),
    LA = mean(la_whc, na.rm = TRUE),
    CV_LA = mean(SD_la_whc, na.rm = TRUE),
    SLA = mean(SLA, na.rm = TRUE),
    CV_SLA = mean(SD_SLA, na.rm = TRUE),
    LDMC = mean(LDMC, na.rm = TRUE),
    CV_LDMC = mean(SD_LDMC, na.rm = TRUE)
  )

############
### 
ggplot(final_merged_df, aes(x = la_whc , y = litter_cover)) +
  geom_point() +
  theme_classic() 

ggplot(final_merged_df, aes(x = la_whc , y = media_coef_trans, group = as.factor(ID_evento_chuva), color= as.factor(ID_evento_chuva))) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Relação entre Throughfall e sd LA",
       y = "throughfall",
       x = "LA") +
  theme_classic()
  


## SLA
ggplot(final_merged_df, aes(x = SLA, y = media_coef_trans, color=grupo_SLA)) +
  geom_point() +
  labs(title = "Relação entre Throughfall e SLA",
       y = "throughfall",
       x = "SLA") +
  theme_classic() +
  scale_color_manual(values = c("LL" = "red", "LH" = "orange", "HH" = "darkgreen", "HL" = "darkblue"))

ggplot(final_merged_df, aes(x = SD_SLA, y = media_coef_trans, color= grupo_SLA)) +
  geom_point() +
  labs(title = "Relação entre Throughfall e sd SLA",
       y = "throughfall",
       x = "sd SLA") +
  theme_classic() +  
  scale_color_manual(values = c("LL" = "red", "LH" = "orange", "HH" = "darkgreen", "HL" = "darkblue"))


### LDMC
ggplot(final_merged_df, aes(x = LDMC, y = media_coef_trans)) +
  geom_point() +
  labs(title = "Relação entre Throughfall e LDMC",
       y = "throughfall",
       x = "LDMC") +
  theme_classic() +
  # Definindo os limites para os eixos x e y
  scale_y_continuous(limits = c(0, 200))

ggplot(final_merged_df, aes(x = SD_LDMC, y = media_coef_trans)) +
  geom_point() +
  labs(title = "Relação entre Throughfall e sd LDMC",
       y = "throughfall",
       x = "sd LDMC") +
  theme_bw() +
  # Definindo os limites para os eixos x e y
  scale_y_continuous(limits = c(0, 200))

## WHC
ggplot(final_merged_df, aes(x = WHCmax, y = media_coef_trans, color=grupo_WHCmax)) +
  geom_point() +
  labs(title = "Relação entre Throughfall e WHC",
       y = "throughfall",
       x = "WHCmax") +
  theme_classic() +  
  scale_color_manual(values = c("LL" = "red", "LH" = "orange", "HH" = "darkgreen", "HL" = "darkblue"))

ggplot(final_merged_df, aes(x = SD_WHC, y = media_coef_trans)) +
  geom_point() +
  labs(title = "Relação entre Throughfall e sd WHC",
       y = "throughfall",
       x = "sd WHC") +
  theme_bw() +
  # Definindo os limites para os eixos x e y
  scale_y_continuous(limits = c(0, 200))


summary(lm(litter_cover ~ CWM_perimeter + FDrich_perimeter + FDdiv_perimeter + CWM_area_cm_quadrados + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + la_whc + cap + SD_CAP + SD_H + h + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC, data=final_merged_df))
summary(lm(CWM_area_cm_quadrados ~ la_whc + cap + SD_CAP + SD_H + h + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC, data=final_merged_df))
summary(lm(FDrich_area_cm_quadrados ~ la_whc + cap + SD_CAP + SD_H + h + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC, data=final_merged_df))
summary(lm(FDdiv_area_cm_quadrados ~ la_whc + cap + SD_CAP + SD_H + h + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC, data=final_merged_df))




library(lme4)

# Fit the linear model with random effects for ID_evento_chuva
model <- lmer(media_coef_trans ~ media_riqueza + media_densidade + media_abertura + media_lai + cap + SD_CAP + SD_H + h + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC + (1 | ID_evento_chuva), data = final_merged_df)

# Summary of the model
summary(model)


library(car)

# Fit the model without random effects to calculate VIF
model_vif <- lm(media_coef_trans ~ media_riqueza + media_densidade + media_abertura + media_lai + cap + SD_CAP + SD_H + h + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC, data = final_merged_df)

# Calculate VIF
vif_values <- vif(model_vif)
print(vif_values)
library(lme4)
library(car)

# Remove highly collinear predictors and refit the model
refined_model_vif <- lmer(media_coef_trans ~ media_riqueza + media_densidade + media_abertura + media_lai + la_whc + SD_la_whc + cap + SD_CAP + SD_H + h + WHC + SD_WHC + (1 | ID_evento_chuva), data = final_merged_df)

# Summary of the refined model
summary(refined_model_vif)

# Calculate VIF for the refined model
refined_model_vif_noRE <- lm(media_coef_trans ~ media_riqueza + media_densidade + media_abertura + media_lai + cap + SD_CAP + SD_H + h + WHC + SD_WHC, data = final_merged_df)
vif_values_refined <- vif(refined_model_vif_noRE)
print(vif_values_refined)


summary(lm(media_coef_trans ~ media_riqueza  + media_densidade + media_abertura + media_lai + cap + SD_CAP + SD_H + h + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC + (1 | ID_evento_chuva) , data=final_merged_df))
summary(lm(cv_coef_trans ~ media_riqueza  + media_densidade + media_abertura + media_lai + cap + SD_CAP + SD_H + h + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC+ (1 | ID_evento_chuva), data=final_merged_df))


summary(lm(media_abertura ~ cap + SD_CAP + SD_H + h + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC, data=final_merged_df))
summary(lm(media_lai ~ cap + SD_CAP + SD_H + h + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + WHC + SD_WHC, data=final_merged_df))

summary(lm(media_coef_esc ~ litter_cover + CWM_perimeter + FDrich_perimeter + FDdiv_perimeter + CWM_area_cm_quadrados + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados, data=final_merged_df))


############
######### 
# Criar matrizes
df_data_floristic2[is.na(df_data_floristic2)] <- 0
matriz_presenca <- dcast(df_data_floristic2, sp ~ etiqueta, value.var = "ID_sp", fun.aggregate = length)

matriz_soma_dbh <- dcast(df_data_floristic2, etiqueta ~ sp, value.var = "dbh", fun.aggregate = sum)
matriz_soma_dbh <- matriz_soma_dbh[rowSums(matriz_soma_dbh[, -1]) > 0, ]
matriz_soma_dbh <- matriz_soma_dbh[, -1]
matriz_soma_dbh <- matriz_soma_dbh[, -1]

df_data_floristic2$canopy_cover[is.na(df_data_floristic2$canopy_cover)] <- 0
matriz_soma_canopy_cover <- dcast(df_data_floristic2,  etiqueta ~ sp, value.var = "canopy_cover", fun.aggregate = sum)
matriz_soma_canopy_cover <- matriz_soma_canopy_cover[complete.cases(matriz_soma_canopy_cover), ]
matriz_soma_canopy_cover <- matriz_soma_canopy_cover[, -1]
matriz_soma_canopy_cover <- matriz_soma_canopy_cover[, -1]
str(df_data_floristic2)


#########################
################
############
matriz_dissimilaridade <- vegdist(matriz_soma_canopy_cover, method = "bray")
nmds_result <- metaMDS(matriz_dissimilaridade)
plot(nmds_result)

num_clusters <- 5  # Defina o número de clusters desejado
cluster_result <- kmeans(nmds_result$points, centers = num_clusters)
nmds_result$clusters <- factor(cluster_result$cluster)
dados_plot <- as.data.frame(nmds_result$points)
dados_plot$etiqueta <- df_floristic$etiqueta
dados_plot$clusters <- nmds_result$clusters


NMDS <- ggplot(data = dados_plot, aes(x = MDS1, y = MDS2, color = clusters, label = etiqueta)) +
  geom_point() +
  geom_text_repel(data = subset(dados_plot, clusters == 1), aes(label = etiqueta), box.padding = 0.5) +
  geom_text_repel(data = subset(dados_plot, clusters == 2), aes(label = etiqueta), box.padding = 0.5) +
  geom_text_repel(data = subset(dados_plot, clusters == 3), aes(label = etiqueta), box.padding = 0.5) +
  geom_text_repel(data = subset(dados_plot, clusters == 4), aes(label = etiqueta), box.padding = 0.5) +
  geom_text_repel(data = subset(dados_plot, clusters == 5), aes(label = etiqueta), box.padding = 0.5) +
    scale_color_manual(values = c("#b22222", "#0060dd", "darkgreen", "grey", "orange")) +
  theme_classic() +
  ggtitle("NMDS Plot pela Cobertura de Copa") + 
  theme(legend.text = element_text(size = 8),
        title = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1),
        axis.ticks.length = unit(0.15, "cm"),
        legend.position = "none")
ggsave("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/NMDS cobertura de copa.png", plot = NMDS, width = 10, height = 6, units = "in", dpi = 300)






## Spatial Autocorrelation!

novo_df_Autocorrelacao <- novo_df[complete.cases(novo_df$media_coef_trans), ]
novo_df_Autocorrelacao$ID_evento_chuva <- as.factor(novo_df_Autocorrelacao$ID_evento_chuva)
novo_df_Autocorrelacao$xcoord <- as.numeric(novo_df_Autocorrelacao$xcoord)
novo_df_Autocorrelacao$ycoord <- as.numeric(novo_df_Autocorrelacao$ycoord)
str(novo_df_Autocorrelacao)

# basics plot 
rbPal <- colorRampPalette(c("#b22222", "#0060dd"))
for (evento in unique(novo_df_Autocorrelacao$ID_evento_chuva)) {
  # Subconjunto dos dados para o evento atual
  subset_data <- novo_df_Autocorrelacao[novo_df_Autocorrelacao$ID_evento_chuva == evento, ]
  
  # Defina os intervalos desejados para as classes
  intervalos <- seq(min(subset_data$media_coef_trans), max(subset_data$media_coef_trans), length.out = 5)
  
  # Rótulos formatados
  rotulos <- sprintf("%0.1f - %0.1f", head(intervalos, -1), tail(intervalos, -1))
  
  # Crie um gráfico ggplot para o evento atual
  p <- ggplot(subset_data, aes(x = xcoord, y = ycoord, color = as.factor(cut(media_coef_trans, breaks = intervalos)))) +
    geom_point(size = 4) +
    labs(title = paste("Mapa de Cores Coef Trans chuva", evento),
         x = "latitude", y = "longitude",
         color = "media_coef_trans") +
    scale_color_manual(values = rbPal(5),
                       labels = rotulos) +
    theme_test() +
    theme(legend.position = "bottom", legend.justification = "right") +
    guides(color = guide_legend(title = NULL))  + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 8),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12)) 
  
  # Salve o gráfico em um arquivo
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Mapa_Cores_evento_", evento, ".png", sep = "")
  ggsave(filename, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  print(p)
}

# resultados Moran
for (evento in unique(novo_df_Autocorrelacao$ID_evento_chuva)) {
   subset_data <- novo_df_Autocorrelacao[novo_df_Autocorrelacao$ID_evento_chuva == evento, ]
    sp_data <- SpatialPointsDataFrame(coords = subset_data[, c("xcoord", "ycoord")],
                                    data = subset_data,
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))
  w <- dnearneigh(coordinates(sp_data), 0, 300)  
   w <- nb2listw(w, style = "W")
    y_values <- as.numeric(subset_data$media_coef_trans)
    moran_result <- moran.test(y_values, listw = w, alternative = "greater")
  print(paste("ID_evento_chuva:", evento))
  print(moran_result)
}


## plots correlagram de Coef Trans with Moran Tests
for (evento in unique(novo_df_Autocorrelacao$ID_evento_chuva)) {
    subset_data <- novo_df_Autocorrelacao[novo_df_Autocorrelacao$ID_evento_chuva == evento, ]
    sp_data <- SpatialPointsDataFrame(coords = subset_data[, c("xcoord", "ycoord")],
                                    data = subset_data,
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  w <- dnearneigh(coordinates(sp_data), 0, 300)  # Ajuste o valor 300 conforme necessário
  w <- nb2listw(w, style = "W")
  y_values <- as.numeric(subset_data$media_coef_trans)
  moran_result <- moran.test(y_values, listw = w, alternative = "greater")
  moran_label <- sprintf("Moran p-value = %.2f", moran_result$p.value)
  rotulos <- sprintf("%0.1f - %0.1f", head(intervalos, -1), tail(intervalos, -1))
  
  p <- ggplot(subset_data, aes(x = xcoord, y = ycoord, color = as.factor(cut(media_coef_trans, breaks = 5)))) +
    geom_point(size = 4) +
    labs(title = paste("Mapa de Cores Coef Trans chuva", evento),
         x = "latitude", y = "longitude",
         color = "media_coef_trans") +
    scale_color_manual(values = rbPal(5),
                       labels = rotulos) +
    theme_test() +
    theme(legend.position = "bottom", legend.justification = "right") +
    guides(color = guide_legend(title = NULL)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 8),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12)) +
    annotate("text", x = min(subset_data$xcoord) + 0.00005, y = max(subset_data$ycoord) - 0.00001,
             label = moran_label)
  
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Mapa_Cores_Moran_evento_", evento, ".png", sep = "")
  print(p)
  ggsave(filename, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  
}

## plots correlagram de Coef Runoff with Moran Tests
novo_df_Autocorrelacao2 <- novo_df[complete.cases(novo_df$media_coef_esc), ]
novo_df_Autocorrelacao2$ID_evento_chuva <- as.factor(novo_df_Autocorrelacao2$ID_evento_chuva)
novo_df_Autocorrelacao2$xcoord <- as.numeric(novo_df_Autocorrelacao2$xcoord)
novo_df_Autocorrelacao2$ycoord <- as.numeric(novo_df_Autocorrelacao2$ycoord)
str(novo_df_Autocorrelacao2)

# basics plot 
rbPal <- colorRampPalette(c("#b22222", "#0060dd"))


for (evento in unique(novo_df_Autocorrelacao2$ID_evento_chuva)) {
  subset_data <- novo_df_Autocorrelacao2[novo_df_Autocorrelacao2$ID_evento_chuva == evento, ]
  sp_data <- SpatialPointsDataFrame(coords = subset_data[, c("xcoord", "ycoord")],
                                    data = subset_data,
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  w <- dnearneigh(coordinates(sp_data), 0, 300)  # Ajuste o valor 300 conforme necessário
  w <- nb2listw(w, style = "W")
  y_values <- as.numeric(subset_data$media_coef_esc)
  moran_result <- moran.test(y_values, listw = w, alternative = "greater")
  moran_label <- sprintf("Moran p-value = %.2f", moran_result$p.value)
  rotulos <- sprintf("%0.1f - %0.1f", head(intervalos, -1), tail(intervalos, -1))
  
  p <- ggplot(subset_data, aes(x = xcoord, y = ycoord, color = as.factor(cut(media_coef_esc, breaks = 5)))) +
    geom_point(size = 4) +
    labs(title = paste("Mapa de Cores Coef Esc chuva", evento),
         x = "latitude", y = "longitude",
         color = "media_coef_esc") +
    scale_color_manual(values = rbPal(5),
                       labels = rotulos) +
    theme_test() +
    theme(legend.position = "bottom", legend.justification = "right") +
    guides(color = guide_legend(title = NULL)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 8),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12)) +
    annotate("text", x = min(subset_data$xcoord) + 0.00005, y = max(subset_data$ycoord) - 0.00001,
             label = moran_label)
  
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Mapa_Cores_Runoff_Moran_evento_", evento, ".png", sep = "")
  print(p)
  ggsave(filename, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  
}

## plots correlagram de Coef Runoff with Moran Tests
novo_df_Autocorrelacao3 <- novo_df[complete.cases(novo_df$media_abertura), ]
novo_df_Autocorrelacao3$ID_evento_chuva <- as.factor(novo_df_Autocorrelacao2$ID_evento_chuva)
novo_df_Autocorrelacao3$xcoord <- as.numeric(novo_df_Autocorrelacao3$xcoord)
novo_df_Autocorrelacao3$ycoord <- as.numeric(novo_df_Autocorrelacao3$ycoord)
str(novo_df_Autocorrelacao3)

# basics plot 
rbPal <- colorRampPalette(c("#b22222", "#0060dd"))


for (evento in unique(novo_df_Autocorrelacao3$ID_evento_chuva)) {
  subset_data <- novo_df_Autocorrelacao3[novo_df_Autocorrelacao3$ID_evento_chuva == evento, ]
  sp_data <- SpatialPointsDataFrame(coords = subset_data[, c("xcoord", "ycoord")],
                                    data = subset_data,
                                    proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  w <- dnearneigh(coordinates(sp_data), 0, 300)  # Ajuste o valor 300 conforme necessário
  w <- nb2listw(w, style = "W")
  y_values <- as.numeric(subset_data$media_abertura)
  moran_result <- moran.test(y_values, listw = w, alternative = "greater")
  moran_label <- sprintf("Moran p-value = %.2f", moran_result$p.value)
  rotulos <- sprintf("%0.1f - %0.1f", head(intervalos, -1), tail(intervalos, -1))
  
  p <- ggplot(subset_data, aes(x = xcoord, y = ycoord, color = as.factor(cut(media_abertura, breaks = 5)))) +
    geom_point(size = 4) +
    labs(title = paste("Mapa de Cores Abertura", evento),
         x = "latitude", y = "longitude",
         color = "media_coef_esc") +
    scale_color_manual(values = rbPal(5),
                       labels = rotulos) +
    theme_test() +
    theme(legend.position = "bottom", legend.justification = "right") +
    guides(color = guide_legend(title = NULL)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(angle = 45, hjust = 1, size = 8),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12)) +
    annotate("text", x = min(subset_data$xcoord) + 0.00005, y = max(subset_data$ycoord) - 0.00001,
             label = moran_label)
  
  filename <- paste("C:/Users/pedro.rajao/Desktop/doutorado UFRJ/8. Cap3/Mapa_Cores_Abertura_Moran_evento_", evento, ".png", sep = "")
  print(p)
  ggsave(filename, plot = p, width = 10, height = 6, units = "in", dpi = 300)
  
}

######################
#############
#######
## SEM Models


library(piecewiseSEM)

final_merged_dfB<-na.omit(final_merged_df)
str(final_merged_dfB)

# Create a ORIGINAL model
model.H1 <- psem(
  glm(itter_cover ~ CWM_perimeter + FDrich_perimeter + FDdiv_perimeter + CWM_area_cm_quadrados + FDdiv_area_cm_quadrados + 
        CWM_perimeter_area_ratio + FDrich_perimeter_area_ratio + FDdiv_perimeter_area_ratio, data=final_merged_df),
  data=final_merged_df
  
)



?psem
# Create a ORIGINAL model
model.H1 <- psem(
  
  glm(media_coef_esc ~ litter_cover + CWM_perimeter + FDrich_perimeter + FDdiv_perimeter + CWM_area_cm_quadrados + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + media_media + media_sd + media_coef_trans + cv_coef_trans + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_df),
  glm(media_coef_trans ~ SD_LDMC + SD_WHCmax + media_riqueza + SD_SLA + media_ndvi + la_whc + LDMC + WHCmax + SLA + media_intensidade_max_mmh + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_df),
  glm(cv_coef_trans ~ LDMC + la_whc + SLA + SD_LDMC + SD_WHCmax + media_riqueza + SD_la_whc + SD_SLA + media_ndvi + media_intensidade_max_mmh + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_df),
  glm(media_ndvi ~ media_abertura + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + media_riqueza, data=final_merged_df),
  glm(litter_cover ~ number_leaf_liter + litter_LA + litter_PER + + litter_PERcv + litter_LAcv + litter_PER + litter_PER_AL + litter_PER_ALcv, data=final_merged_df),
  data=final_merged_df
  
)

# Summary of the pSEM
summary(model.H1)
plot(model.H1)


# Create a ORIGINAL model
model.H2 <- psem(
  
  lm(media_coef_esc ~ la_whc + SD_la_whc + media_lai + litter_LAcv + media_coef_trans + number_leaf_liter + media_media + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_coef_trans ~ SD_WHCmax + media_riqueza + SD_SLA + media_lai + la_whc + LDMC + litter_LAcv + SLA + media_intensidade_max_mmh, data=final_merged_dfB),
  lm(cv_coef_trans ~ media_riqueza + SD_la_whc + SD_WHCmax + media_lai + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_lai ~ SD_WHCmax + number_leaf_liter + litter_LAcv + media_abertura + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + media_riqueza, data=final_merged_dfB),
  
  data=final_merged_dfB
  
)
summary(model.H2)
plot(model.H2) #Global goodness-of-fit: Chi-Squared = 39.705 with P-value = 0.089 and on 29 degrees of freedom / Fisher's C = 45.562 with P-value = 0.882 and on 58 degrees of freedom


model.H3 <- psem(
  lm(media_coef_esc ~ litter_LA  + litter_LAcv  + WHCmax + SD_SLA + la_whc + SD_la_whc + SLA + media_ndvi + media_coef_trans + cv_coef_trans, data=final_merged_dfB),
  lm(media_coef_trans ~ LDMC + SLA + SD_SLA + SD_LDMC + media_intensidade_max_mmh, data=final_merged_dfB),
  lm(cv_coef_trans ~ media_riqueza + SD_LDMC + SD_la_whc + SD_WHCmax + media_ndvi + media_chuva_acumulada_mm, data=final_merged_dfB),
  lm(media_ndvi ~ media_riqueza + WHCmax + la_whc + SD_la_whc + SLA + SD_SLA + LDMC, data=final_merged_dfB),
  data=final_merged_dfB
)

summary(model.H3)
plot(model.H3)


model.H4 <- psem(
  
  lm(media_coef_esc ~ media_ndvi + SD_la_whc + litter_cover + litter_LA + media_media + media_sd + media_coef_trans + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_coef_trans ~ litter_cover + litter_PER_ALcv + media_ndvi + la_whc + SD_la_whc + media_sd + SD_LDMC + media_riqueza + SD_SLA + SLA + media_intensidade_max_mmh, data=final_merged_dfB),
  lm(cv_coef_trans ~ media_intensidade_max_mmh + litter_PERcv + media_riqueza + SD_LDMC + litter_LA + SD_WHCmax + SD_la_whc + media_ndvi + media_chuva_acumulada_mm, data=final_merged_dfB),
  lm(media_ndvi ~ media_sd + media_sd + litter_cover + litter_PER_ALcv + litter_PERcv + SD_WHCmax + media_media + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + media_riqueza, data=final_merged_dfB),
  lm(litter_cover ~  media_media + LDMC + la_whc + SD_la_whc + SD_WHCmax + SLA + SD_SLA + SD_LDMC + media_sd + number_leaf_liter + litter_LA + litter_PERcv + litter_PER_ALcv, data=final_merged_dfB),
  data=final_merged_dfB
  
)

summary(model.H4)
plot(model.H4)


model.H5 <- psem(
  lm(turbidez ~ litter_LA + litter_PER + + litter_PERcv + litter_LAcv + media_media + media_sd + litter_cover + media_coef_esc + media_intensidade_max_mmh + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_coef_esc ~ litter_cover + litter_LA + litter_PER + + litter_PERcv + litter_LAcv + media_media + media_sd + media_coef_trans + cv_coef_trans + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_coef_trans ~ SD_LDMC + SD_WHCmax + media_riqueza + SD_SLA + media_ndvi + la_whc + LDMC + WHCmax + SLA + media_intensidade_max_mmh + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(cv_coef_trans ~ LDMC + la_whc + SLA + SD_LDMC + SD_WHCmax + media_riqueza + SD_la_whc + SD_SLA + media_ndvi + media_intensidade_max_mmh + media_chuva_acumulada_mm + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_ndvi ~ media_abertura + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + SD_LDMC + media_riqueza, data=final_merged_dfB),
  lm(litter_cover ~ number_leaf_liter + litter_LA + litter_PER + + litter_PERcv + litter_LAcv + litter_PER + litter_PER_AL + litter_PER_ALcv, data=final_merged_dfB),
  data=final_merged_dfB
  
)
summary(model.H5)
plot(model.H5)



model.H6 <- psem(
  lm(turbidez ~ litter_LA + litter_PERcv + litter_cover + media_coef_esc + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_coef_esc ~ media_ndvi + SD_la_whc + litter_cover + litter_LA + media_media + media_sd + media_coef_trans + media_dias_sem_chuva1mm, data=final_merged_dfB),
  lm(media_coef_trans ~ media_ndvi + media_sd + SD_LDMC + media_intensidade_max_mmh, data=final_merged_dfB),
  lm(cv_coef_trans ~ media_intensidade_max_mmh + litter_PER + media_riqueza + SD_LDMC + litter_PERcv + litter_LA + SD_WHCmax + SD_la_whc + media_ndvi + media_chuva_acumulada_mm, data=final_merged_dfB),
  lm(media_ndvi ~ litter_PER + media_sd + litter_cover + litter_PER_ALcv + SD_WHCmax + media_media + litter_PERcv + la_whc + SD_la_whc + SLA + SD_SLA + LDMC + media_riqueza, data=final_merged_dfB),
  lm(litter_cover ~ media_media + LDMC + SD_SLA + SLA + la_whc + SD_la_whc + SD_WHCmax +  SD_LDMC + media_sd + number_leaf_liter + litter_LA + litter_PERcv + litter_PER + litter_PER_ALcv, data=final_merged_dfB),
  data=final_merged_dfB
  
)
summary(model.H6)
plot(model.H6)
model.H6