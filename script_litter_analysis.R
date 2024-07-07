
## pacotes
library(vegan)
library(factoextra) 
library(ggpubr)
library(ggplot2)
library(dplyr)
library(PMCMRplus)
library(PerformanceAnalytics)
library(gsubfn)
library(openxlsx)
library(FD)
library(tidyr)
library(lme4)
library(MuMIn)
library(purrr)
library(piecewiseSEM)
library(tidyverse)


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
    media_transprecipitacao_mm = mean(transprecipitacao_mm, na.rm = TRUE),
    cv_transprecipitacao_mm = sd(transprecipitacao_mm, na.rm = TRUE) / mean(transprecipitacao_mm, na.rm = TRUE)*100,
    media_escoamento_mm = mean(escoamento_mm, na.rm = TRUE) * 100,
    cv_escoamento_mm = sd(escoamento_mm, na.rm = TRUE) / mean(escoamento_mm, na.rm = TRUE),
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
################
########## relações floristicas

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

diretorio <- "C:/Users/pedro.rajao/Desktop/doutorado UFRJ/9. Cap4_CWM-FD role on rainfall interception canopy-litter/resultados_fito_final.txt"
write.table(resultados_fito, diretorio, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


############
######### 
# Diversidade Funcional

#df_completoA <- merge(df_data_floristic2, df_leaf_traits, by = "sp")
df_completo <- left_join(df_data_floristic2, df_leaf_traits4, by = "sp")
df_completo <- df_completo %>%
  mutate_all(~as.numeric(gsub(",", ".", .)))


# List of columns that are supposed to be numeric
numeric_cols <- c("cap", "dbh", "h", "canopy_cover", "area_basal", "la_whc", "la_sla", 
                  "perimeter", "peso_0_sla", "peso_seco_sla", "peso_0_whc", "peso_1_whc_escorrido", 
                  "peso_24_whc_escorrido", "peso_1_whc_seco.papel", "peso_24_whc_seco.papel", 
                  "LRET", "LREP")

# Apply the transformation to the selected columns
df_completo <- df_completo %>%
  mutate(across(all_of(numeric_cols), ~as.numeric(gsub(",", ".", .))))


#TRAITS
df_completo$SLA <- df_completo$la_sla / df_completo$peso_seco_sla
df_completo$LDMC <- df_completo$peso_seco_sla / df_completo$peso_0_sla
df_completo$WHC_RET <- (df_completo$peso_24_whc_escorrido - df_completo$peso_0_whc) / df_completo$peso_0_whc
df_completo$AHC <- (df_completo$peso_24_whc_seco.papel - df_completo$peso_0_whc) / df_completo$peso_0_whc

###### LA
cwm_vars <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("la_whc"), weighted.mean, w = n, na.rm = TRUE))

fd <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_la_whc = sd(la_whc, na.rm = TRUE)
  )

fdRIC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_LA = if_else(all(is.na(la_whc)), 
                      NA_real_, 
                      max(la_whc, na.rm = TRUE) - min(la_whc, na.rm = TRUE))
  )

###### CAP
cwm_vars_CAP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("cap"), weighted.mean, w = n, na.rm = TRUE))

fd_CAP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_CAP = sd(cap, na.rm = TRUE)
  )
fdRIC_CAP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_CAP = if_else(all(is.na(cap)), 
                       NA_real_, 
                       max(cap, na.rm = TRUE) - min(cap, na.rm = TRUE))
  )

###### H
cwm_vars_H <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("h"), weighted.mean, w = n, na.rm = TRUE))

fd_H <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_H = sd(h, na.rm = TRUE)
  )
fdRIC_H <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_H = if_else(all(is.na(h)), 
                     NA_real_, 
                     max(h, na.rm = TRUE) - min(h, na.rm = TRUE))
  )

###### SLA
cwm_vars_SLA <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("SLA"), weighted.mean, w = n, na.rm = TRUE))

fd_SLA <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_SLA = sd(SLA, na.rm = TRUE)
  )
fdRIC_SLA <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_SLA = if_else(all(is.na(SLA)), 
                       NA_real_, 
                       max(SLA, na.rm = TRUE) - min(SLA, na.rm = TRUE))
  )

###### LDMC
cwm_vars_LDMC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("LDMC"), weighted.mean, w = n, na.rm = TRUE))

fd_LDMC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_LDMC = sd(LDMC, na.rm = TRUE)
  )
fdRIC_LDMC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_LDMC = if_else(all(is.na(LDMC)), 
                        NA_real_, 
                        max(LDMC, na.rm = TRUE) - min(LDMC, na.rm = TRUE))
  )

###### WHC_RET
cwm_vars_WHC_RET <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("WHC_RET"), weighted.mean, w = n, na.rm = TRUE))

fd_WHC_RET <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_WHC_RET = sd(WHC_RET, na.rm = TRUE)
  )
fdRIC_WHC_RET <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_WHC_RET = if_else(all(is.na(WHC_RET)), 
                           NA_real_, 
                           max(WHC_RET, na.rm = TRUE) - min(WHC_RET, na.rm = TRUE))
  )

###### WHC
cwm_vars_AHC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("AHC"), weighted.mean, w = n, na.rm = TRUE))

fd_AHC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_AHC = sd(AHC, na.rm = TRUE)
  )
fdRIC_AHC <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_AHC = if_else(all(is.na(AHC)), 
                       NA_real_, 
                       max(AHC, na.rm = TRUE) - min(AHC, na.rm = TRUE))
  )

###### LRET
cwm_vars_LRET <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("LRET"), weighted.mean, w = n, na.rm = TRUE))

fd_LRET <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_LRET = sd(LRET, na.rm = TRUE)
  )
fdRIC_LRET <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_LRET = if_else(all(is.na(LRET)), 
                        NA_real_, 
                        max(LRET, na.rm = TRUE) - min(LRET, na.rm = TRUE))
  )


###### LREP
cwm_vars_LREP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("LREP"), weighted.mean, w = n, na.rm = TRUE))

fd_LREP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_LREP = sd(LREP, na.rm = TRUE)
  )
fdRIC_LREP <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_LREP = if_else(all(is.na(LREP)), 
                        NA_real_, 
                        max(LREP, na.rm = TRUE) - min(LREP, na.rm = TRUE))
  )


# Lista de data frames para cada variável
la_whc_list <- list(cwm_vars, fd, fdRIC)
cap_list <- list(cwm_vars_CAP, fd_CAP, fdRIC_CAP)
h_list <- list(cwm_vars_H, fd_H, fdRIC_H)
sla_list <- list(cwm_vars_SLA, fd_SLA, fdRIC_SLA)
ldmc_list <- list(cwm_vars_LDMC, fd_LDMC, fdRIC_LDMC)
whc_ret_list <- list(cwm_vars_WHC_RET, fd_WHC_RET, fdRIC_WHC_RET)
ahc_list <- list(cwm_vars_AHC, fd_AHC, fdRIC_AHC)
lret_list <- list(cwm_vars_LRET, fd_LRET, fdRIC_LRET)
lrep_list <- list(cwm_vars_LREP, fd_LREP, fdRIC_LREP)

# Função para juntar listas de data frames por "etiqueta"
join_by_etiqueta <- function(df_list) {
  reduce(df_list, left_join, by = "etiqueta")
}

# Juntando cada variável individualmente
resultados_la_whc <- join_by_etiqueta(la_whc_list)
resultados_cap <- join_by_etiqueta(cap_list)
resultados_h <- join_by_etiqueta(h_list)
resultados_sla <- join_by_etiqueta(sla_list)
resultados_ldmc <- join_by_etiqueta(ldmc_list)
resultados_whc_ret <- join_by_etiqueta(whc_ret_list)
resultados_ahc <- join_by_etiqueta(ahc_list)
resultados_lret <- join_by_etiqueta(lret_list)
resultados_lrep <- join_by_etiqueta(lrep_list)

# Lista de resultados finais de cada variável
resultados_list <- list(
  resultados_la_whc, resultados_cap, resultados_h, resultados_sla, 
  resultados_ldmc, resultados_whc_ret, resultados_ahc, resultados_lret, resultados_lrep
)

# Juntando todos os resultados finais em um único data frame
resultados_funcional <- reduce(resultados_list, left_join, by = "etiqueta")
resultados_funcional$etiqueta <- as.integer(resultados_funcional$etiqueta)

final_merged_dfA <- left_join(resultados_funcional, novo_df, by = "etiqueta")

final_merged_df <- final_merged_dfA %>%
  mutate(
    classe_vol_chuva = case_when(
      media_chuva_acumulada_mm < 15 ~ "leve",
      media_chuva_acumulada_mm >= 15 & media_chuva_acumulada_mm < 30 ~ "moderada",
      media_chuva_acumulada_mm >= 30 ~ "forte"
    ),
    classe_intensidade_chuva = case_when(
      media_intensidade_max_mmh < 3 ~ "leve",
      media_intensidade_max_mmh >= 3 & media_intensidade_max_mmh < 15 ~ "media",
      media_intensidade_max_mmh >= 15 ~ "forte"
    )
  )

final_merged_df <- final_merged_df %>%
  slice(-980)
str(final_merged_df)


summary(lmer(media_coef_esc ~ litter_cover + (1 | ID_evento_chuva), data=final_merged_df))



# Fit the mixed-effects model
model <- lmer(media_coef_esc ~ litter_cover + (litter_cover | ID_evento_chuva), data=final_merged_df)
r_squared_values <- r.squaredGLMM(model)
print(r_squared_values)

model <- lmer(media_coef_esc ~ CWM_area_cm_quadrados + (1 | ID_evento_chuva), data=final_merged_df)
r_squared_values <- r.squaredGLMM(model)
print(r_squared_values)
summary(model)

model <- lmer(litter_cover ~ cap + SD_CAP + SIZE_CAP + SD_H + h + SIZE_H + la_whc + SD_la_whc + SIZE_LA + SLA + SD_SLA + SIZE_SLA + LDMC + SD_LDMC + SIZE_LDMC + AHC + SD_AHC + SIZE_AHC + WHC_RET + SD_WHC_RET + SIZE_WHC_RET + (1 | ID_evento_chuva), data=final_merged_df)
r_squared_values <- r.squaredGLMM(model)
print(r_squared_values)

summary(lm(litter_cover ~ cap + SD_CAP + SIZE_CAP + SD_H + h + SIZE_H + la_whc + SIZE_LA + SLA + SD_SLA + SIZE_SLA + LDMC + SIZE_LDMC + AHC + SIZE_AHC + SD_WHC_RET + SIZE_WHC_RET + (1 | ID_evento_chuva), data=final_merged_df))
model <- lmer(litter_cover ~ cap + SD_CAP + SIZE_CAP + SD_H + h + SIZE_H + la_whc + SIZE_LA + SLA + SD_SLA + SIZE_SLA + LDMC + SIZE_LDMC + AHC + SIZE_AHC + SD_WHC_RET + SIZE_WHC_RET + (1 | ID_evento_chuva), data=final_merged_df)
r_squared_values <- r.squaredGLMM(model)
print(r_squared_values)
summary(model)

plot(model)

ggplot(final_merged_df, aes(x = media_abertura, y = litter_cover)) +
geom_point(size = 3) +  
geom_smooth(method = "glm", family = gaussian(link = "gaussian"), se = T, linetype = "dashed") +
theme_classic() +
theme(axis.text.x = element_text(size = 38),
        axis.text.y = element_text(size = 38),
        legend.text = element_text(size = 15)
  )

summary(lm(media_abertura ~ cap + SD_CAP + SIZE_CAP + SD_H + h + SIZE_H + la_whc + SD_la_whc + SIZE_LA + SLA + SD_SLA + SIZE_SLA + LDMC + SD_LDMC + SIZE_LDMC + AHC + SD_AHC + SIZE_AHC + WHC_RET + SD_WHC_RET + SIZE_WHC_RET + (1 | ID_evento_chuva), data=final_merged_df))


df_filtro_chuva_leve <- final_merged_df %>%
  filter(classe_vol_chuva %in% 'leve')
mixed_model <- lmer(media_coef_esc ~ CWM_area_cm_quadrados + (1 | ID_evento_chuva), data = df_filtro_chuva_leve)
  ggplot(df_filtro_chuva_leve, aes(x = CWM_area_cm_quadrados, y = media_coef_esc, color = as.factor(ID_evento_chuva), shape = as.factor(ID_evento_chuva))) +
  geom_point(size = 5) +  # Ajustando o tamanho dos pontos para melhor visibilidade
  geom_smooth(method = "glm", se = FALSE, aes(group = ID_evento_chuva, color = as.factor(ID_evento_chuva)), linetype = "dashed")  +  # Ajustando os formatos de acordo com o número de cores
  theme_classic() 

df2_filtro_chuva_leve <- df_filtro_chuva_leve %>%
    filter(ID_evento_chuva %in% c(2, 21)) 
ggplot(df2_filtro_chuva_leve, aes(x = CWM_area_cm_quadrados, y = media_coef_esc, color = as.factor(ID_evento_chuva), shape = as.factor(ID_evento_chuva))) +
  geom_point(size = 8) +  # Ajustando o tamanho dos pontos para melhor visibilidade
  geom_smooth(method = "glm", se = FALSE, aes(color = as.factor(ID_evento_chuva)), linetype = "dashed")  + # Ajustando os formatos de acordo com o número de cores
  theme_classic() +
  scale_color_manual(values = c('#00BFFF', '#00FFFF')) +
  theme(axis.text.x = element_text(size = 38),
        axis.text.y = element_text(size = 38),
        legend.text = element_text(size = 15)
  )

  
  df_filtro_chuva_moderada <- final_merged_df %>%
    filter(classe_vol_chuva %in% 'moderada')
  mixed_model <- lmer(media_coef_esc ~ CWM_area_cm_quadrados + (1 | ID_evento_chuva), data = df_filtro_chuva_moderada)
  ggplot(df_filtro_chuva_moderada, aes(x = CWM_area_cm_quadrados, y = media_coef_esc, color = as.factor(ID_evento_chuva), shape = as.factor(ID_evento_chuva))) +
    geom_point(size = 3) +  # Ajustando o tamanho dos pontos para melhor visibilidade
    geom_smooth(method = "glm", se = FALSE, aes(group = ID_evento_chuva, color = as.factor(ID_evento_chuva)), linetype = "dashed")  +  # Ajustando os formatos de acordo com o número de cores
    theme_classic() 
  
  df2_filtro_chuva_moderada <- df_filtro_chuva_moderada %>%
    filter(ID_evento_chuva %in% c(25, 26)) 
  ggplot(df2_filtro_chuva_moderada, aes(x = CWM_area_cm_quadrados, y = media_coef_esc, color = as.factor(ID_evento_chuva), shape = as.factor(ID_evento_chuva))) +
    geom_point(size = 8) +  # Ajustando o tamanho dos pontos para melhor visibilidade
     geom_smooth(method = "glm", se = FALSE, aes(color = as.factor(ID_evento_chuva)), linetype = "dashed")  +  # Ajustando os formatos de acordo com o número de cores
    theme_classic() +
    scale_color_manual(values = c('#4169E1', '#1E90FF')) +
    theme(axis.text.x = element_text(size = 38),
          axis.text.y = element_text(size = 38),
          legend.text = element_text(size = 15)
    ) # Definindo a cor dos pontos como 'lightblue' e a linha suavizada como 'blue'
  
  
  df_filtro_chuva_forte <- final_merged_df %>%
    filter(classe_vol_chuva %in% 'forte')
  mixed_model <- lmer(media_coef_esc ~ CWM_area_cm_quadrados + (1 | ID_evento_chuva), data = df_filtro_chuva_forte)
  ggplot(df_filtro_chuva_forte, aes(x = CWM_area_cm_quadrados, y = media_coef_esc, color = as.factor(ID_evento_chuva), shape = as.factor(ID_evento_chuva))) +
    geom_point(size = 3) +  # Ajustando o tamanho dos pontos para melhor visibilidade
    geom_smooth(method = "glm", se = FALSE, aes(group = ID_evento_chuva, color = as.factor(ID_evento_chuva)), linetype = "dashed")  +  # Ajustando os formatos de acordo com o número de cores
    theme_classic() 
  df2_filtro_chuva_forte <- df_filtro_chuva_forte %>%
    filter(ID_evento_chuva %in% c(22, 23, 27)) 
  ggplot(df2_filtro_chuva_forte, aes(x = CWM_area_cm_quadrados, y = media_coef_esc, color = as.factor(ID_evento_chuva), shape = as.factor(ID_evento_chuva))) +
    geom_point(size = 8) +  # Ajustando o tamanho dos pontos para melhor visibilidade
    geom_smooth(method = "glm", se = FALSE, aes(color = as.factor(ID_evento_chuva)), linetype = "dashed")  +  # Ajustando os formatos de acordo com o número de cores
    theme_classic() +
    scale_color_manual(values = c('#191970', '#483D8B', '#008080')) +
    theme(axis.text.x = element_text(size = 38),
          axis.text.y = element_text(size = 38),
          legend.text = element_text(size = 15)
    ) # Definindo a cor dos pontos e da linha suavizada
  
  
  
  final_merged_dfC <- final_merged_df[!is.na(final_merged_df$media_coef_esc), ]

  
  
  # Create a ORIGINAL model
  model.H1 <- psem(
    lmer(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + ((1 + CWM_area_cm_quadrados + CWM_perimeter) | ID_evento_chuva), data = final_merged_df ),
    lmer(media_coef_esc ~ litter_cover + (1 + litter_cover | ID_evento_chuva), data = final_merged_df )
     )
  summary(model.H1)
 
  
  # Create a df_filtro_chuva_leve model
  model.Hleve <- psem(
    lmer(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + (1 + CWM_area_cm_quadrados + CWM_perimeter| ID_evento_chuva), data = df_filtro_chuva_leve),
    lmer(media_coef_esc ~ litter_cover + (1 + litter_cover | ID_evento_chuva), data = df_filtro_chuva_leve)
  )
  summary(model.Hleve)
  plot(model.Hleve)
  
  # Create a df_filtro_chuva_moderada model
  model.Hmoderada <- psem(
    lmer(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + (1 + CWM_area_cm_quadrados + CWM_perimeter| ID_evento_chuva), data = df_filtro_chuva_moderada),
    lmer(media_coef_esc ~ litter_cover + (1 + litter_cover| ID_evento_chuva), data = df_filtro_chuva_moderada)
  )
  summary(model.Hmoderada)
  
  # Create a df_filtro_chuva_forte model
  model.Hforte <- psem(
    lmer(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + (1 + CWM_area_cm_quadrados + CWM_perimeter| ID_evento_chuva), data = df_filtro_chuva_forte),
    lmer(media_coef_esc ~ litter_cover + (1 + litter_cover| ID_evento_chuva), data = df_filtro_chuva_forte)
  )
  summary(model.Hforte)
  
  
  