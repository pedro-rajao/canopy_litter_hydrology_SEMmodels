
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

columns_with_commas <- c("area_basal", "la_whc", "la_sla", "perimeter", 
                         "peso_0_sla", "peso_seco_sla", "peso_0_whc", 
                         "peso_1_whc_escorrido", "peso_24_whc_escorrido", 
                         "peso_1_whc_seco.papel", "peso_24_whc_seco.papel", 
                         "LRET", "LREP")

# Substituir vírgulas por pontos nas colunas identificadas
df_completo[columns_with_commas] <- lapply(df_completo[columns_with_commas], function(x) {
  as.numeric(gsub(",", ".", x))
})

#TRAITS
df_completo$SLA <- df_completo$la_sla / df_completo$peso_seco_sla
df_completo$LDMC <- df_completo$peso_seco_sla / df_completo$peso_0_sla
df_completo$WHC_RET <- (df_completo$peso_24_whc_escorrido - df_completo$peso_0_whc) / df_completo$peso_0_whc
df_completo$AHC <- (df_completo$peso_24_whc_seco.papel - df_completo$peso_0_whc) / df_completo$peso_0_whc
df_completo$LA_PER <- df_completo$la_whc / df_completo$perimeter


str(comm_matrix)
trait_matrix2 <- trait_matrix[-1, ]


###### LA

cwm_vars <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("la_whc"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("cap"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("h"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("SLA"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("LDMC"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("WHC_RET"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("AHC"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("LRET"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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
  summarise(across(starts_with("LREP"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

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

###### LA_PER

cwm_vars_LA_PER <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(across(starts_with("LA_PER"), ~ weighted.mean(., w = canopy_cover, na.rm = TRUE)))

fd_LA_PER <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SD_LA_PER = sd(LA_PER, na.rm = TRUE)
  )

fdRIC_LA_PER <- df_completo %>%
  group_by(etiqueta) %>%
  summarise(
    SIZE_LA_PER = if_else(all(is.na(LA_PER)), 
                      NA_real_, 
                      max(la_whc, na.rm = TRUE) - min(la_whc, na.rm = TRUE))
  )



# Lista de data frames para cada variável
la_per_list <- list(cwm_vars_LA_PER, fd_LA_PER, fdRIC_LA_PER)
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
resultados_la_per <- join_by_etiqueta(la_per_list)
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
  resultados_la_per, resultados_la_whc, resultados_cap, resultados_h, resultados_sla, 
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
final_merged_df <- final_merged_df %>%
  mutate(ID_evento_chuva = as.factor(ID_evento_chuva))
str(final_merged_df)
final_merged_df <- final_merged_df %>%
  mutate(across(where(is.numeric), scale))

traits_cor_SSS <- final_merged_df[, c("la_whc", "SD_la_whc", "SIZE_LA", "LA_PER", "SD_LA_PER", "SIZE_LA_PER",
                                      "cap", "SD_CAP", "SIZE_CAP", "h", "SD_H", "SIZE_H")]

traits_cor_HT <- final_merged_df[, c("WHC_RET", "SD_WHC_RET", "SIZE_WHC_RET", "AHC", "SD_AHC", "SIZE_AHC", 
                                     "LRET", "SD_LRET", "SIZE_LRET", "LREP", "SD_LREP", "SIZE_LREP")] 

traits_cor_RES <- final_merged_df[, c("LDMC", "SIZE_LDMC", "SD_LDMC", "SLA", "SIZE_SLA", "SD_SLA")]

cor_matrix <- cor(traits_cor_RES)
heatmap(cor_matrix)
chart.Correlation(traits_cor_RES)

######################################################### REVISAR OS MODELOS!

model<-lmer(media_escoamento_mm ~ litter_cover + media_chuva_acumulada_mm + media_intensidade_max_mmh + media_transprecipitacao_mm + cv_transprecipitacao_mm + media_media + media_sd + (1 | ID_evento_chuva), data = final_merged_df)
model2<-lmer(media_escoamento_mm ~ litter_cover + media_transprecipitacao_mm +  (1 | ID_evento_chuva), data = final_merged_df)
r_squared_values <- r.squaredGLMM(model2)
print(r_squared_values)
summary(model2)

model<-lmer(media_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_LA_PER + cap+ SD_CAP+ h + SD_H + WHC_RET + SD_WHC_RET + AHC + SD_AHC + LRET + SD_LRET + LREP+ SD_LREP + LDMC + SD_LDMC + SLA + SD_SLA + (1  | ID_evento_chuva), data=final_merged_df)
model2<-lmer(media_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_H + SD_AHC + SD_LRET + SD_LDMC + (1  | ID_evento_chuva), data=final_merged_df)
r_squared_values <- r.squaredGLMM(model2)
print(r_squared_values)
summary(model2)


model<-(lmer(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_LA_PER + cap+ SD_CAP+ h + SD_H + WHC_RET + SD_WHC_RET + AHC + SD_AHC + LRET + SD_LRET + LREP+ SD_LREP + LDMC + SD_LDMC + SLA + SD_SLA + (1  | ID_evento_chuva), data=final_merged_df))
model2<-(lmer(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + la_whc + LA_PER + cap+ SD_CAP+ h + SD_H + WHC_RET + LRET + SD_LRET + LREP+ LDMC + SD_LDMC + SD_SLA + (1  | ID_evento_chuva), data=final_merged_df))
r_squared_values <- r.squaredGLMM(model2)
print(r_squared_values)
summary(model2)

model<-(lmer(media_abertura ~ SD_LDMC + SD_SLA + la_whc + LDMC  + SLA +  SD_la_whc + LA_PER + SD_LA_PER + cap + SD_CAP + h + SD_H  + (1 | ID_evento_chuva), data=final_merged_df))
model2<-(lmer(media_abertura ~ SD_SLA + la_whc + SLA +  SD_la_whc + LA_PER + SD_LA_PER + cap + SD_H  + (1 | ID_evento_chuva), data=final_merged_df))
r_squared_values <- r.squaredGLMM(model2)
print(r_squared_values)
summary(model2)

model<-lmer(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + FDrich_perimeter + FDdiv_perimeter + (1  | ID_evento_chuva), data = final_merged_df)
r_squared_values <- r.squaredGLMM(model)
print(r_squared_values)
summary(model)

model<-lmer(turbidez ~ CWM_area_cm_quadrados + CWM_perimeter + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + FDrich_perimeter + FDdiv_perimeter + litter_cover + media_escoamento_mm + media_chuva_acumulada_mm + media_intensidade_max_mmh + (1  | ID_evento_chuva), data = final_merged_df)
model2<-lmer(turbidez ~ FDdiv_area_cm_quadrados + media_escoamento_mm + (1  | ID_evento_chuva), data = final_merged_df)
r_squared_values <- r.squaredGLMM(model2)
print(r_squared_values)
summary(model2)



##################################
###########################
##### média por tipo de chuva.

means_df <- final_merged_df %>%
  group_by(etiqueta, classe_vol_chuva) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

df_filtro_chuva_forte <- means_df %>%
  filter(classe_vol_chuva %in% 'forte')
df_filtro_chuva_moderada <- means_df %>%
  filter(classe_vol_chuva %in% 'moderada')
df_filtro_chuva_leve <- means_df %>%
  filter(classe_vol_chuva %in% 'leve')


model.MEDIA_1 <- psem(
  glm(media_escoamento_mm ~ litter_cover + media_transprecipitacao_mm, family = gaussian(), data = final_merged_df),
  glm(media_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_H + SD_AHC + SD_LRET + SD_LDMC, family = gaussian(), data = means_df),
  glm(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + la_whc + LA_PER + cap + SD_CAP + h + SD_H + WHC_RET + LRET + LREP + LDMC + SD_LDMC + SD_SLA, family = gaussian(), data = means_df),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc + cap, family = gaussian(), data = means_df),
  glm(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + FDrich_perimeter + FDdiv_perimeter, family = gaussian(), data = means_df)
)
summary(model.MEDIA_1)

model.MEDIA_2 <- psem(
  glm(media_escoamento_mm ~ SD_SLA + la_whc + SD_LA_PER + SLA + LDMC + WHC_RET + h + SD_LDMC + media_chuva_acumulada_mm + litter_cover + media_transprecipitacao_mm, family = gaussian(), data = means_df),
  glm(media_transprecipitacao_mm ~ SD_la_whc + la_whc + LREP + media_chuva_acumulada_mm + LA_PER + SD_H + SD_AHC + SD_LRET + SD_LDMC, family = gaussian(), data = means_df),
  glm(cv_transprecipitacao_mm ~ SD_la_whc + SD_AHC + media_transprecipitacao_mm + media_chuva_acumulada_mm + la_whc + LA_PER + cap + SD_CAP + h + SD_H + WHC_RET + LRET + SD_LRET + LREP + LDMC + SD_LDMC + SD_SLA, family = gaussian(), data = means_df),
  glm(media_abertura ~ SD_SLA + la_whc + SLA + SD_la_whc + LA_PER + SD_LA_PER + cap + SD_H, family = gaussian(), data = means_df),
  glm(litter_cover ~ SD_SLA + LRET + SD_CAP + cap + SLA + SD_LA_PER + LDMC + LREP + WHC_RET + h + SD_LDMC + SD_LRET + SD_H + LA_PER + SD_la_whc + la_whc + CWM_area_cm_quadrados + CWM_perimeter + FDrich_perimeter + FDdiv_perimeter, family = gaussian(), data = means_df)
)
summary(model.MEDIA_2)
#Chi-Squared = 162.858 with P-value = 0 and on 55 degrees of freedom
#Fisher's C = 123.08 with P-value = 0.186 and on 110 degrees of freedom
plot(model.MEDIA_2)



model.MEDIA_LEVE_1 <- psem(
  glm(media_escoamento_mm ~ litter_cover + media_transprecipitacao_mm, family = gaussian(), data = df_filtro_chuva_leve),
  glm(media_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_H + SD_AHC + SD_LRET + SD_LDMC, family = gaussian(), data = df_filtro_chuva_leve),
  glm(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + la_whc + LA_PER + cap + SD_CAP + h + SD_H + WHC_RET + LRET + LREP + LDMC + SD_LDMC + SD_SLA, family = gaussian(), data = df_filtro_chuva_leve),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc + cap, family = gaussian(), data = df_filtro_chuva_leve),
  glm(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + FDrich_perimeter + FDdiv_perimeter, family = gaussian(), data = df_filtro_chuva_leve)
)
summary(model.MEDIA_LEVE_1) 
#Chi-Squared = 220.345 with P-value = 0 and on 96 degrees of freedom
#Fisher's C = 213.471 with P-value = 0.138 and on 192 degrees of freedom


model.MEDIA_LEVE_2 <- psem(
  glm(media_escoamento_mm ~ litter_cover + media_transprecipitacao_mm, family = gaussian(), data = df_filtro_chuva_leve),
  glm(media_transprecipitacao_mm ~ 1, family = gaussian(), data = df_filtro_chuva_leve),
  glm(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + la_whc + SD_CAP + h + SD_H + WHC_RET + LRET + LREP + LDMC, family = gaussian(), data = df_filtro_chuva_leve),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc, family = gaussian(), data = df_filtro_chuva_leve),
  glm(litter_cover ~ LA_PER + SD_H + SD_LRET + SD_LDMC + h + SLA + CWM_area_cm_quadrados , family = gaussian(), data = df_filtro_chuva_leve)
)
summary(model.MEDIA_LEVE_2)   
#Chi-Squared = 60.535 with P-value = 0.125 and on 49 degrees of freedom
#Fisher's C = 79.918 with P-value = 0.909 and on 98 degrees of freedom
plot(model.MEDIA_LEVE_2)   



model.MEDIA_MODERADA_1 <- psem(
  glm(media_escoamento_mm ~ litter_cover + media_transprecipitacao_mm, family = gaussian(), data = df_filtro_chuva_moderada),
  glm(media_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_H + SD_AHC + SD_LRET + SD_LDMC, family = gaussian(), data = df_filtro_chuva_moderada),
  glm(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + la_whc + LA_PER + cap + SD_CAP + h + SD_H + WHC_RET + LRET + LREP + LDMC + SD_LDMC + SD_SLA, family = gaussian(), data = df_filtro_chuva_moderada),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc + cap, family = gaussian(), data = df_filtro_chuva_moderada),
  glm(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + FDrich_perimeter + FDdiv_perimeter, family = gaussian(), data = df_filtro_chuva_moderada)
)
summary(model.MEDIA_MODERADA_1)

model.MEDIA_MODERADA_2 <- psem(
  glm(media_escoamento_mm ~ litter_cover , family = gaussian(), data = df_filtro_chuva_moderada),
  glm(media_transprecipitacao_mm ~ SD_H + SD_LDMC + SD_la_whc + LREP , family = gaussian(), data = df_filtro_chuva_moderada),
  glm(cv_transprecipitacao_mm ~ SD_LDMC + LREP + la_whc + SD_CAP + SD_H + SD_SLA, family = gaussian(), data = df_filtro_chuva_moderada),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc , family = gaussian(), data = df_filtro_chuva_moderada),
  glm(litter_cover ~ SLA +  h + SD_LDMC + SD_LRET + SD_H + LA_PER + CWM_area_cm_quadrados , family = gaussian(), data = df_filtro_chuva_moderada)
)
summary(model.MEDIA_MODERADA_2)
#Chi-Squared = 44.659 with P-value = 0.65 and on 49 degrees of freedom
#Fisher's C = 74.915 with P-value = 0.96 and on 98 degrees of freedom
plot(model.MEDIA_MODERADA_2)



model.MEDIA_FORTE_1 <- psem(
  glm(media_escoamento_mm ~ litter_cover + media_transprecipitacao_mm, family = gaussian(), data = df_filtro_chuva_forte),
  glm(media_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_H + SD_AHC + SD_LRET + SD_LDMC, family = gaussian(), data = df_filtro_chuva_forte),
  glm(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + la_whc + LA_PER + cap + SD_CAP + h + SD_H + WHC_RET + LRET + LREP + LDMC + SD_LDMC + SD_SLA, family = gaussian(), data = df_filtro_chuva_forte),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc + cap, family = gaussian(), data = df_filtro_chuva_forte),
  glm(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + FDrich_perimeter + FDdiv_perimeter, family = gaussian(), data = df_filtro_chuva_forte)
)
summary(model.MEDIA_FORTE_1)


model.MEDIA_FORTE_2 <- psem(
  glm(media_escoamento_mm ~ LDMC + SD_LDMC + litter_cover, family = gaussian(), data = df_filtro_chuva_forte),
  glm(media_transprecipitacao_mm ~ SD_la_whc + la_whc + LA_PER + SD_H + SD_LDMC, family = gaussian(), data = df_filtro_chuva_forte),
  glm(cv_transprecipitacao_mm ~ SD_la_whc + LREP + LDMC + SD_SLA, family = gaussian(), data = df_filtro_chuva_forte),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc , family = gaussian(), data = df_filtro_chuva_forte),
  glm(litter_cover ~ SLA + h + SD_LDMC + SD_LRET + SD_H + LA_PER + CWM_area_cm_quadrados , family = gaussian(), data = df_filtro_chuva_forte)
)
summary(model.MEDIA_FORTE_2)
##Chi-Squared = 44.921 with P-value = 0.6 and on 48 degrees of freedom
#Fisher's C = 68.204 with P-value = 0.986 and on 96 degrees of freedom
plot(model.MEDIA_FORTE_2)






model.MEDIA_1 <- psem(
  glm(turbidez ~ FDdiv_area_cm_quadrados + media_escoamento_mm, family = gaussian(), data = final_merged_df),
  glm(media_escoamento_mm ~ litter_cover + media_transprecipitacao_mm, family = gaussian(), data = final_merged_df),
  glm(media_transprecipitacao_mm ~ media_chuva_acumulada_mm + media_intensidade_max_mmh + la_whc + SD_la_whc + LA_PER + SD_H + SD_AHC + SD_LRET + SD_LDMC, family = gaussian(), data = means_df),
  glm(cv_transprecipitacao_mm ~ media_chuva_acumulada_mm + la_whc + LA_PER + cap + SD_CAP + h + SD_H + WHC_RET + LRET + LREP + LDMC + SD_LDMC + SD_SLA, family = gaussian(), data = means_df),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc + cap, family = gaussian(), data = means_df),
  glm(litter_cover ~ CWM_area_cm_quadrados + CWM_perimeter + FDrich_area_cm_quadrados + FDdiv_area_cm_quadrados + FDrich_perimeter + FDdiv_perimeter, family = gaussian(), data = means_df)
)
summary(model.MEDIA_1)


model.MEDIA_2 <- psem(
  glm(turbidez ~ SD_CAP + SLA + LDMC + LREP + LRET + LA_PER + media_intensidade_max_mmh + media_chuva_acumulada_mm + media_escoamento_mm, family = gaussian(), data = final_merged_df),
  glm(media_escoamento_mm ~ LDMC + WHC_RET + SD_LDMC + LA_PER + media_chuva_acumulada_mm + litter_cover + media_transprecipitacao_mm, family = gaussian(), data = final_merged_df),
  glm(media_transprecipitacao_mm ~ la_whc + SD_la_whc + LREP + media_chuva_acumulada_mm + LA_PER + SD_H + SD_AHC + SD_LDMC, family = gaussian(), data = means_df),
  glm(cv_transprecipitacao_mm ~ media_transprecipitacao_mm + media_chuva_acumulada_mm + la_whc + cap + SD_CAP + h + SD_H + WHC_RET + LRET + LREP + LDMC + SD_LDMC + SD_SLA, family = gaussian(), data = means_df),
  glm(media_abertura ~ la_whc + SLA + SD_la_whc + cap, family = gaussian(), data = means_df),
  glm(litter_cover ~ cap + SD_CAP + SLA + LDMC + LREP + LRET + WHC_RET + h + SD_LDMC + SD_LRET + SD_H + LA_PER + SD_la_whc + la_whc + CWM_area_cm_quadrados + CWM_perimeter + FDrich_perimeter + FDdiv_perimeter, family = gaussian(), data = means_df)
)
summary(model.MEDIA_2)
#Chi-Squared = 3157.873 with P-value = 0 and on 90 degrees of freedom
#Fisher's C = 189.082 with P-value = 0.237 and on 176 degrees of freedom
plot(model.MEDIA_2)



