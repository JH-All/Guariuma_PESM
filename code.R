# Install/Load packages ------------------------------------
packages = c("adespatial", "tidyverse", "vegan", "readxl", "hrbrthemes",
             "viridis", "ggbeeswarm", "ggthemes", "iNEXT",
             "spaa", "cowplot", "FactoMineR", "factoextra",
             "writexl", "fields", "reshape2", "ade4", "readr", "car",
             "MASS", "broom", "ggrepel", "grid")
lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    } 
    library(x, character.only = TRUE)
  }
)

# Getting data ready ----------------------
data = read_excel("dados.xlsx")
com = data[,16:23]
env = data[,c(6:10)]

data$Comp_maximo_m <- as.numeric(data$Comp_maximo_m)
data$Larg_maxima_m <- as.numeric(data$Larg_maxima_m)
data$Prof_maxima_cm <- as.numeric(data$Prof_maxima_cm)
data$depth_m <- data$Prof_maxima_cm/ 100

V = (2/3) * pi * data$Comp_maximo_m * data$Larg_maxima_m  * data$depth_m
env$Volume <- V
data$Distancia_riacho_proximo = as.numeric(data$Distancia_riacho_proximo)
env$Stream_dist= as.numeric(data$Distancia_riacho_proximo)
env$pH = as.numeric(env$pH)
env$OD = as.numeric(env$OD)
env$Temperatura = as.numeric(env$Temperatura)
env$Condutividade = as.numeric(env$Condutividade)
env$Sol_totais = as.numeric(env$Sol_totais)
env$Periodo = data$Periodo
data$Periodo = as.factor(data$Periodo)

data$Abundance = rowSums(com)
data$S = specnumber(com)
data$H = diversity(com)
data$effective = exp(data$H)
env$abund = data$Abundance
env$camarao = data$Camarão
env$S = data$S 
env$effec = data$effective


# PCA -------------------------------
data$Periodo <- dplyr::recode(data$Periodo,
                       "Seco" = "Dry",
                       "Chuvoso" = "Wet")

env$Periodo <- dplyr::recode(env$Periodo,
                       "Seco" = "Dry",
                       "Chuvoso" = "Wet")

env <- env %>%
  rename(
    Temperature = Temperatura,
    DO          = OD,
    TS          = Sol_totais,
    Cond        = Condutividade
  )

pca.p <- PCA(X = env[,c(1:7)],
             scale.unit = TRUE, graph = FALSE)

fviz_screeplot(pca.p, addlabels = TRUE, ylim = c(0, 70), main = "", 
               xlab = "Dimensões",
               ylab = "Porcentagem de variância explicada") 
var_env <- get_pca_var(pca.p)
summary(pca.p)

fig2 = fviz_pca_biplot(pca.p,
                geom.ind = "point", 
                fill.ind = env$Periodo, 
                col.ind = "black",
                alpha.ind = 0.55,
                pointshape = 21,
                pointsize = 7,
                palette = c("#33a02c", "#1f78b4"), 
                col.var = "black",
                invisible = "quali",
                title = NULL) +
  labs(x = "PC1 (38%)", y = "PC2 (21%)", fill = NULL)+
  theme_classic(base_size = 18)+
  scale_x_continuous(limits = c(-4,5))

fig2

ggsave("Figure_2.jpg", fig2)

# iNEXT -------------------------
datlist <- list()
com_inext <- com

com_inext_dp <- decostand(com_inext[1:15, 1:6 ],
                          method = "pa")
com_inext_wp <- decostand(com_inext[16:34, c(1:5, 7:8)],
                          method = "pa")

datlist$WP <-data.frame(t(com_inext_wp))
datlist$DP <-data.frame(t(com_inext_dp))

datlist

result <- iNEXT(datlist,
                q = 0,
                datatype = "incidence_raw",
                endpoint = 30,
                se = TRUE, 
                nboot = 999)


fig3 = ggiNEXT(result, type = 1) +
  scale_color_manual(
    values = c("#33a02c", "#1f78b4"),
    labels = c("Dry", "Wet")
  ) +
  scale_fill_manual(
    values = c("#33a02c", "#1f78b4"),
    labels = c("Dry", "Wet")
  )+
  guides(shape = "none")+
  theme_classic(base_size = 18)+
  labs(x = "Number of sampled pools", y = "Species richness")

fig3

ggsave("Figure_3.jpg", fig3)

# Abundance and diversity  ----------------
fig4_a = data %>% 
  ggplot(aes(x = Periodo, y = Abundance, fill = Periodo))+
  geom_boxplot(width = 0.45, show.legend = F, alpha = 0.8)+
  geom_jitter(width = 0.12, shape = 21, size = 4.5,
              show.legend = F, alpha = 0.6)+
  scale_fill_manual(values = c("#33a02c", "#1f78b4"))+
  theme_classic(base_size = 18)+
  labs(x = NULL)

fig4_a

wilcox.test(Abundance ~ Periodo, data = data) # W = 151, p = 0.78

fig4_b = data %>% 
  ggplot(aes(x = Periodo, y = S, fill = Periodo))+
  geom_boxplot(width = 0.45, show.legend = F, alpha = 0.8)+
  geom_jitter(width = 0.12, shape = 21, size = 4.5,
              show.legend = F, alpha = 0.6)+
  scale_fill_manual(values = c("#33a02c", "#1f78b4"))+
  theme_classic(base_size = 18)+
  labs(x = NULL, y = "Species richness")

fig4_b

wilcox.test(S ~ Periodo, data = data) # W = 168, p = 0.36


fig4_c = data %>% 
  ggplot(aes(x = Periodo, y = effective, fill = Periodo))+
  geom_boxplot(width = 0.45, show.legend = F, alpha = 0.8)+
  geom_jitter(width = 0.12, shape = 21, size = 4.5,
              show.legend = F, alpha = 0.6)+
  scale_fill_manual(values = c("#33a02c", "#1f78b4"))+
  theme_classic(base_size = 18)+
  labs(x = NULL, y = "Effective number of species")

fig4_c

wilcox.test(effective ~ Periodo, data = data) # W = 158.5, p = 0.58

fig_4 = plot_grid(fig4_a, fig4_b, fig4_c, labels = "AUTO", nrow = 1)
fig_4

ggsave("Figure_4.jpg", fig_4, width = 12, height = 4)

# Predictors -----------------
env_dry <- filter(env, Periodo == "Dry")
env_wet <- filter(env, Periodo == "Wet")
preds <- c("Volume", "Stream_dist", "camarao")

vif_dry <- lm(abund ~ Volume + Stream_dist + camarao, data = env_dry)
vif(vif_dry)

vif_wet <- lm(abund ~ Volume + Stream_dist + camarao, data = env_wet)
vif(vif_wet)

m_abund_dry <- glm.nb(abund ~ Volume + Stream_dist + camarao, data = env_dry)
m_abund_wet <- glm.nb(abund ~ Volume + Stream_dist + camarao, data = env_wet)
m_S_dry <- glm.nb(S ~ Volume + Stream_dist + camarao, data = env_dry)
m_S_wet <- glm.nb(S ~ Volume + Stream_dist + camarao, data = env_wet)
m_eff_dry <- glm(effec ~ Volume + Stream_dist + camarao,
                 data = env_dry,
                 family = Gamma(link = "log"))
m_eff_wet <- glm(effec ~ Volume + Stream_dist + camarao,
                 data = env_wet,
                 family = Gamma(link = "log"))

extract_coefs <- function(model, response, period){
  tidy(model) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      response = response,
      Periodo = period,
      sig = ifelse(p.value < 0.05, "significant", "ns")
    )
}

coefs <- bind_rows(
  extract_coefs(m_abund_dry, "Abundance", "Dry"),
  extract_coefs(m_abund_wet, "Abundance", "Wet"),
  extract_coefs(m_S_dry, "Richness", "Dry"),
  extract_coefs(m_S_wet, "Richness", "Wet"),
  extract_coefs(m_eff_dry, "Effective", "Dry"),
  extract_coefs(m_eff_wet, "Effective", "Wet")
)

coefs2 <- coefs %>%
  mutate(
    response = dplyr::recode(response,
                             "Effective" = "ENS"),
    response = factor(response, levels = c("Abundance", "Richness", "ENS")),
    sig = factor(sig, levels = c("ns", "significant"))
  )

coefs2 <- coefs2 %>%
  mutate(
    term = dplyr::recode(term,
                  "Stream_dist" = "Stream distance",
                  "camarao"     = "Shrimp",
                  "Volume"      = "Volume")
  )

fig_5 = ggplot(coefs2, aes(x = estimate, y = term, color = Periodo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error),
                 height = 0.2, linewidth = 0.9) +
  geom_point(aes(fill = ifelse(sig == "significant", as.character(Periodo), "ns")),
             shape = 21, size = 4, stroke = 1.2) +
  scale_color_manual(values = c("Dry" = "#33a02c", "Wet" = "#1f78b4")) +
  scale_fill_manual(values = c("Dry" = "#33a02c",
                               "Wet" = "#1f78b4",
                               "ns" = "white")) +
  coord_cartesian(xlim = c(-0.20, 0.20)) +
  facet_wrap(~ response, scales = "fixed") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "Model coefficient (± SE)",
       y = NULL)

fig_5

ggsave("Figure_5.jpg", fig_5, width = 12, height = 6)

# PERMANOVA --------------------------
data$Periodo <- as.factor(data$Periodo)
levels(data$Periodo)

distancia <- vegdist(com, method = "bray")

resultado_permanova <- adonis2(distancia ~ data$Periodo, permutations = 999)
print(resultado_permanova) 
dispersao <- betadisper(distancia,data$Periodo)
anova(dispersao)  
# PERMANOVA: F = 0.62, p = 0.701; PERMDISP: F = 0.96, p = 0.334

# RDA -------------------------------
env <- env %>%
  rename(Temp = Temperature)

run_rda_period <- function(period, com, env, vars){
  
  env_p <- env %>% filter(Periodo == period)
  
  com_p <- com[env$Periodo == period, , drop = FALSE]
  
  com_p <- com_p[, colSums(com_p, na.rm = TRUE) > 0, drop = FALSE]
  
  com_hel <- decostand(com_p, method = "hellinger")
  
  form <- as.formula(paste("com_hel ~", paste(vars, collapse = " + ")))
  rda_mod <- rda(form, data = env_p)
  
  list(
    rda = rda_mod,
    anova_global = anova(rda_mod, permutations = 999),
    anova_terms  = anova(rda_mod, by = "terms", permutations = 999),
    vif          = vif.cca(rda_mod),
    eigen        = eigenvals(rda_mod),
    var_exp      = eigenvals(rda_mod) / sum(eigenvals(rda_mod)) * 100,
    species_included = colnames(com_p)
  )
}

env_vars <- c("Stream_dist", "DO", "Temp", "Volume", "TS")

rda_dry <- run_rda_period("Dry", com, env, env_vars)
rda_wet <- run_rda_period("Wet", com, env, env_vars)

rda_dry$anova_global
rda_wet$anova_global

anova(rda_dry$rda, by = "axis", permutations = 999)
anova(rda_wet$rda, by = "axis", permutations = 999)

rda_dry$anova_terms
rda_wet$anova_terms

rda_dry$vif
rda_wet$vif

round(rda_dry$var_exp[1:2], 2)
round(rda_wet$var_exp[1:2], 2)

plot_rda <- function(rda_mod, period = c("Dry", "Wet"), title = NULL,
                     show_species_labels = TRUE){
  
  period <- match.arg(period)
  
  period_cols <- c("Dry" = "#33a02c", "Wet" = "#1f78b4")
  arrow_col <- period_cols[period]
  
  sc <- 2
  scr <- scores(rda_mod, scaling = sc)
  
  species <- as.data.frame(scr$species)
  envs    <- as.data.frame(scr$biplot)
  
  species$sp <- rownames(species)
  envs$var   <- rownames(envs)
  
  if (is.null(title)) title <- paste(period, "period")
  
  p <- ggplot() +
    geom_segment(data = envs,
                 aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                 arrow = arrow(length = unit(3, "mm")),
                 linewidth = 0.9,
                 color = arrow_col) +
    geom_text(data = envs,
              aes(x = RDA1 * 1.10, y = RDA2 * 1.10, label = var),
              color = arrow_col, size = 4) +
    geom_point(data = species,
               aes(x = RDA1, y = RDA2),
               shape = 21, fill = "black", alpha = 0.5, size = 2) +
    theme_classic(base_size = 16) +
    labs(x = "RDA1", y = "RDA2", title = title)
  
  if (show_species_labels) {
    p <- p +
      ggrepel::geom_text_repel(data = species,
                               aes(x = RDA1, y = RDA2, label = sp),
                               size = 3,
                               color = "black",
                               max.overlaps = 30,
                               box.padding = 0.25,
                               point.padding = 0.15,
                               seed = 1)
  }
  
  p
}

fig6_a = plot_rda(rda_dry$rda, period = "Dry", title = "Dry period")+
  labs(title = NULL, x = "RDA1 (52.67%)", y = "RDA2 (7.81%)")+
  theme_classic(base_size = 18)
fig6_b = plot_rda(rda_wet$rda, period = "Wet", title = "Wet period")+
  labs(title = NULL, x = "RDA1 (9.59%)", y = "RDA2 (3.81%)")+
  theme_classic(base_size = 18)+
  scale_x_continuous(limits = c(-0.9, 0.9))

fig6_a
fig6_b

fig_6 = plot_grid(fig6_a, fig6_b, labels = "AUTO", nrow = 1)
ggsave("Figure_6.jpg", fig_6, width = 10, height = 4)
