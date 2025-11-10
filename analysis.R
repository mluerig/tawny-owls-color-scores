# 00a - setup --------------------------------------------------------------

# devtools::install_github("https://github.com/alishinski/lavaanPlot/tree/dir")
# devtools::install_github("Sebastien-Le/YesSiR")
# library(YesSiR) # to export a flextable into MS Excel: exportxlsx() function

## "pacman" does automatic package loading
if (!require("pacman")) install.packages("pacman")
pkgs<-c(
  "data.table","mgcv","scales", "lme4", "lmerTest",
  "nlme","gratia","emmeans","nadiv", "sommer",
  "lavaan","lavaanPlot", "qgraph",
  "ggplot2", "cowplot", "ggridges","viridis",
  "flextable","YesSiR", 
  "plot3Drgl", "rgl",
  "stringr",
  "marginaleffects"
)
pacman::p_load(char=pkgs,install=F)
rm(pkgs)

m_r = function(x, digits=10) list(
  MEAN = round(mean(x, na.rm=T), digits)
)
m_se = function(x, digits=10) list(
  MEAN = round(mean(x, na.rm=T),digits), 
  SE = round((sd(x,na.rm=T) / sqrt(length(x))), digits)
  )
outlier_mad <- function(x, threshold = 3) {
  mad_deviate <- (x - stats::median(x, na.rm = T)) / stats::mad(x, na.rm = T)
  abs(mad_deviate) > threshold
}


# 00b - set WD ------------------------------------------------------------------

setwd("D:/git-repos/data-repos/tawny-owls-color-scores")

dir.create("figures")
dir.create("tables")

# 01 - define variables ---------------------------------------------------

## set theme
theme_set(theme_cowplot())

## set global colors
morph_cols = c(
  "gray"="gray",
  "brown"="chocolate4"
  )

## color score labels
col_score_labs_old = c(
  "1 [4]", "2 [5]", "3 [6]", "4 [7]", "5 [8]", " [9]", 
  "1 [10]", "2 [11]", "3 [12]", "4 [13]", "5 [14]")
col_score_labs = c(
  "1","2","3","4","5","",
  "1","2","3","4","5"
  )

# 02 - load and format data -----------------------------------------------

## load invidiuals data
data_ind = fread("data/data_individuals_clean_env_masked.csv", na.strings = "", stringsAsFactors = T)
data_ind[, morph := factor(morph, levels=c("gray","brown"))]
data_ind[, ring := factor(ring)]
data_ind[, couple_id := factor(couple_id)]
data_ind[, nest_id := factor(nest_id)]
data_ind[, year_index := year - min(year) + 1]
data_ind = copy(data_ind[year >= 1980,])

## unique data - ring only appears on first occurence 
data_ind_uni = data_ind[!duplicated(ring),] 

## load recruits file
data_recr = fread("data/data_recruits_env_masked.csv", stringsAsFactors = T)
data_recr[, morph := factor(morph, levels=c("gray","brown"))]
data_recr[, parent_combination := factor(
  parent_combination, levels=c("gray_gray","gray_brown", "brown_gray","brown_brown"))]

## merge recruitment success into individuals
data_recr_m = data_recr[!str_detect(couple_id, "NA"), .(couple_id, year)]
data_recr_m = data_recr_m[,.N, by=.(couple_id, year)]
setnames(data_recr_m, "N", "n_recruits")
data_recr_m[, recr_success := 1]
data_ind = merge(data_ind, data_recr_m, all.x=T)
data_ind[is.na(recr_success), recr_success:=0]

# 03 - prep animal model --------------------------------------------------

## pedigree
data_pedigree = data.table(nadiv::prepPed(data.table(
  animal=factor(data_recr$ring),
  sire=factor(data_recr$ring_m),
  dam=factor(data_recr$ring_f)
)[, .(animal, sire, dam)]))

## additive & dominace matrix
A <- as.matrix(makeA(data_pedigree))
D <- as.matrix(makeD(data_pedigree)$D)

## animal model data (not all vars are used)
data_animal_mod = data_recr[, .(
  ring, morph, nest_id, sex, year, laying_date, breed_snow_depth_mean, breed_snow_days_n, breed_temp_air_mean, breed_temp_air_CV, 
  color_score_median, color_score_median_trans, parent_combination, colour_score_median_midparent, colour_score_median_midparent_trans)]
data_animal_mod[, nest_id := factor(nest_id)]
data_animal_mod[, ring := factor(ring)]
data_animal_mod[, dom := factor(ring, levels = rownames(D))]

## addive variance coding
data_animal_mod[, geno_add := NA_real_]
data_animal_mod[parent_combination == "gray_gray", geno_add := 0]
data_animal_mod[parent_combination == "brown_brown", geno_add := 2]
data_animal_mod[parent_combination %in% c("gray_brown", "brown_gray"), geno_add := 1]

## dominance coding 
data_animal_mod[, geno_dom := NA_real_]
data_animal_mod[parent_combination == "gray_gray", geno_dom := 0]
data_animal_mod[parent_combination == "brown_brown", geno_dom := 0]
data_animal_mod[parent_combination %in% c("gray_brown", "brown_gray"), geno_dom := 1]

## maternal effects (ultimately not included)
data_animal_mod[, mother_id := as.factor(data_pedigree[match(ring, animal), dam])]
maternal_geno <- data_animal_mod[, .(ring, geno_add, geno_dom)]
setnames(maternal_geno, c("ring", "geno_add", "geno_dom"), c("dam", "mat_geno_add", "mat_geno_dom"))
data_animal_mod <- merge(data_animal_mod, maternal_geno, by.x = "mother_id", by.y = "dam", all.x = TRUE)

## check data_animal_mod
table(data_animal_mod$parent_combination, data_animal_mod$geno_add)
table(data_animal_mod$parent_combination, data_animal_mod$geno_dom)
aggregate(cbind(geno_add, geno_dom) ~ parent_combination, data = data_animal_mod, mean, na.rm = TRUE)


# figure 1b ----------------------------------------------------------------

data = copy(data_ind)

labels = c("1980-1984", "1985-1989", "1990-1994","1995-1999",
           "2000-2004", "2005-2009","2010-2014", "2015-2019", "2020-2022")
data[, years := cut(
  year, 
  breaks=c(1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025),
  include.lowest=F,right=F,
  labels=labels),]

p1 =
ggplot(data) + xlab("Color score") +
  geom_density_ridges(
    aes(y=years, x=color_score_median, fill=morph, group=interaction(years, morph)), 
    bandwidth = 0.5) + #
  geom_vline(aes(xintercept=9), size=0.7, linetype=2) +
  scale_y_discrete(limits=rev, expand = c(0.025,0)) +
  scale_x_continuous(limits=c(1,17), breaks=c(4,9,14)) +
  scale_fill_manual("Morphs", values=c("gray","chocolate4"), labels=c("Gray morph", "Brown morph")) +
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        legend.justification = "right",
        legend.direction = "horizontal",
        legend.key = element_rect(fill = NA),
        axis.text.y = element_text(size=10),
        plot.margin = unit(c(15,5.5,5.5,10), "pt")
        )

ggsave(paste0("figures/figure1b.png"), 
       width=4, height=5, units = c("in"), bg="white")

 # figure 2 stats ----------------------------------------------------------------

data = copy(data_ind)

## color score
gam1 = gam(color_score_median_trans_ord ~ morph + 
              s(year_index, by=morph, k=5) +
              s(observer, by=morph, bs = 're'),
            method = "REML",
            data=data)
summary(gam1)
capture.output(anova(gam1), file = "tables/GAM1.txt")
lm1 = glm(color_score_median_trans_ord ~ 
              year_index * morph,
            data=data)
capture.output(anova(lm1), file = "tables/LM1.txt")
capture.output(summary(emtrends(lm1, specs = "morph", var = "year_index"), infer = c(TRUE, TRUE)), file = "tables/LM1_trend.txt")

## laying date
gamS2 = gam(laying_date ~ morph +
              s(year_index, by=morph, k=5),
            method = "REML",
            data=data)
capture.output(anova(gamS2), file = "tables/GAMS2.txt")
lm2 = glm(laying_date ~ year_index * morph, data=data[!is.na(laying_date)])
capture.output(anova(lm2), file = "tables/LM2.txt")
capture.output(summary(emtrends(lm2, specs = "morph", var = "year_index"), infer = c(TRUE, TRUE)), file = "tables/LM2_trend.txt")

## recruitment success
gamS3 = gam(recr_success  ~ morph +
              s(year_index, by=morph, k=5),
            family = "binomial",
            method = "REML",
            data=data)
capture.output(anova(gamS3), file = "tables/GAMS3.txt")
lm3 <- glm(recr_success ~ year_index * morph,
             family = binomial(link = "logit"),
             data = data)
capture.output(anova(lm3), file = "tables/LM3.txt")
capture.output(summary(emtrends(lm3, specs = "morph", var = "year_index"), infer = c(TRUE, TRUE)), file = "tables/LM3_trend.txt")


## color score
data_pred_gam1 <- data.table(predictions(gam1, exclude=list(
  "s(observer):morphgray", "s(observer):morphbrown"), newdata = data.table(expand.grid(
    year_index = 1:43, morph=c("gray", "brown"), observer=c("a","b","c")))))
data_pred_gam1 = data_pred_gam1[observer=="a"]
data_pred_gam1[morph == "gray",  c("estimate", "conf.low", "conf.high") := 
                 .(estimate + 3, conf.low + 3, conf.high + 3)]
data_pred_gam1[morph == "brown", c("estimate", "conf.low", "conf.high") := 
                 .(estimate + 9, conf.low + 9, conf.high + 9)]
data_pred_gam1[, year := year_index + 1979]

data_pred_lm1 <- data.table(predictions(lm1, newdata = data.table(expand.grid(
  year_index = 1:43, morph=c("gray", "brown")))))
data_pred_lm1[morph == "gray",  c("estimate", "conf.low", "conf.high") := 
                 .(estimate + 3, conf.low + 3, conf.high + 3)]
data_pred_lm1[morph == "brown", c("estimate", "conf.low", "conf.high") := 
                 .(estimate + 9, conf.low + 9, conf.high + 9)]
data_pred_lm1[, year := year_index + 1979]

## laying date
data_pred_gamS2 <- data.table(predictions(gamS2, newdata = data.table(expand.grid(
  year_index = 1:43, morph=c("gray", "brown")
))))
data_pred_gamS2[, year := year_index + 1979]
data_pred_gamS2_summ = data_pred_gamS2[, as.list(unlist(lapply(.SD, mean))), 
                                       by=c("year", "year_index"),
                                       .SDcols= c("estimate", "conf.low", "conf.high")] 

data_pred_lm2 <- data.table(predictions(lm2))
data_pred_lm2[, year := year_index + 1979]
data_pred_lm2_summ = data_pred_lm2[, as.list(unlist(lapply(.SD, mean))), 
                                       by=c("year", "year_index"),
                                       .SDcols= c("estimate", "conf.low", "conf.high")] 


## recruitment success
data_pred_gamS3 <- data.table(predictions(gamS3, newdata = data.table(expand.grid(
  year_index = 1:43, morph=c("gray", "brown")
))))
data_pred_gamS3[, year := year_index + 1979]
data_pred_gamS3_summ = data_pred_gamS3[, as.list(unlist(lapply(.SD, mean))), 
                                       by=c("year", "year_index"),
                                       .SDcols= c("estimate", "conf.low", "conf.high")] 

data_pred_lm3 <- data.table(predictions(lm3, newdata = data.table(expand.grid(
  year_index = 1:43, morph=c("gray", "brown")
))))
data_pred_lm3[, year := year_index + 1979]
data_pred_lm3_summ = data_pred_lm3[, as.list(unlist(lapply(.SD, mean))), 
                                       by=c("year", "year_index"),
                                       .SDcols= c("estimate", "conf.low", "conf.high")] 

# figure 2 ----------------------------------------------------------------

## create summary
data_summ1 = data[
  , as.list(unlist(lapply(.SD, m_se))), 
  by=c("morph", "year", 'year_index' ),
  .SDcols= c("color_score_median","color_score_median_trans_ord","laying_date", "recr_success")] 

## pull environmental data (averaged by year)
cols_env = c("year",
             "winter_snow_depth_mean", "winter_snow_days_n",
             "winter_temp_air_mean", "winter_temp_air_CV")
data_env = unique(data_ind[,..cols_env])

## calculate count numbers
data_n = data[, .N, by=.(year, morph)]
data_n[, n_year := sum(N), by=year]
data_n[, prop := N/n_year]
data_n[, morph := factor(morph, levels=c("gray","brown"))]
data_n[, linetype := "Proportion"]

## fake plot for legend
legend_data <- data.table(x = 2010, y = 0, xend = 2012, yend = 0, group = "Highlight")
p_leg<-
  ggplot(data[, .N, by=.(year, morph)]) +
  geom_errorbar(aes(x=year, ymin=N-1, ymax=N+1, group=morph)) +
  geom_point(aes(x=year, y=N, fill=morph), pch=21) +
  geom_line(data=data_pred_gam1, aes(x=year, y=estimate, group=morph, linetype=morph), size=1) +
  geom_ribbon(data=data_pred_gam1, aes(x=year, ymin=conf.low, ymax=conf.high, group=morph, linetype=morph), alpha=0.2) +
  geom_segment(data=legend_data, aes(x=x, y=y, xend=xend, yend=yend, color=group), size=5, alpha=0.3) +
  scale_fill_manual("Morph", values=morph_cols, labels=c("Gray", "Brown")) +
  scale_linetype_manual("Model fits", values=c(1,2), labels=c("Non-linear fits (see Table 1 / S2)", "Linear fits (see Table S2)")) +
  scale_color_manual("Extreme winter conditions", values=c("Highlight"="red"), labels="2010-2013") +
  guides(fill = guide_legend(order=1, title.position="top", override.aes=list(color="black", size=3)),
         linetype = guide_legend(order=2, keywidth=unit(2, 'cm'), title.position="top"),
         color = guide_legend(order=3, keywidth=unit(2, 'cm'), title.position="top")) +
  theme(legend.position="bottom",
        legend.text = element_text(size=12),
        legend.justification="center",
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.direction="vertical")
leg = cowplot::get_plot_component(p_leg, "guide-box", return_all = TRUE)[[3]]


p1 =
ggplot(data_summ1) + ggtitle("A") +
  coord_cartesian(ylim = c(4,14)) +  ylab("Color score") +
  annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
  geom_ribbon(data=data_pred_gam1, aes(x=year, ymin=conf.low, ymax=conf.high, group=morph), alpha=0.2) +
  geom_errorbar(aes(x=year, ymin=color_score_median.MEAN-color_score_median.SE, 
                    ymax=color_score_median.MEAN+color_score_median.SE, group=morph), 
                color="black", width=0) +
  geom_point(aes(x=year, y=color_score_median.MEAN, fill=morph), size=2, pch=21) +
  geom_hline(yintercept = 9, linetype=2) + 
  geom_line(data=data_pred_gam1, aes(x=year, y=estimate, color=morph), size=1, linetype=1) +
  geom_line(data=data_pred_lm1, aes(x=year, y=estimate, color=morph), 
            linewidth=1, linetype=2) +
  geom_line(data=data_pred_gam1, aes(x=year, y=estimate, group=morph), 
            linewidth=0.25, color="black") +
  geom_line(data=data_pred_lm1, aes(x=year, y=estimate, group=morph), 
            linewidth=0.25, color="black", linetype=2) +
  scale_y_continuous(
    breaks=c(4:14),
    labels = col_score_labs_old ) +
  scale_fill_manual("Morph", values=morph_cols) +
  scale_color_manual("Morph", values=morph_cols) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size=12),
    axis.title.x = element_blank()
  )

p2 =
  ggplot(data_n) + ggtitle("B") +
    ylab("Morph abundance") +
    annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
    geom_bar(aes(x=year, y=N, fill=morph), stat="identity", size=0.25, color="black") +
    geom_line(data=data_n[morph=="brown"], aes(x=year, y=prop*100, color=linetype),  linewidth=0.8) +
    scale_fill_manual("Morphs", values=morph_cols, guide="none") +
    scale_color_manual(values = "chocolate4", labels="Proportion of brown morph") + 
    scale_y_continuous(expand = c(0, 0), breaks=seq(0,100,20), limits = c(0,100)) +
    theme(legend.position.inside = c(0.35,1),
          legend.position = "inside",
          legend.title = element_blank(),
          legend.key.width=unit(1.5,"cm"),
          legend.text = element_text(size=10),
          legend.key = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=12))

p3 =
  ggplot(data_summ1) + ggtitle("C") +
  geom_abline(slope=0) + ylab("Laying date") + 
  annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
  coord_cartesian(ylim=c(-25,25)) + 
  scale_color_manual(values=morph_cols) +
  geom_ribbon(data=data_pred_gamS2, aes(x=year, ymin=conf.low, ymax=conf.high, group=morph), alpha=0.2) +
  geom_errorbar(aes(x=year, ymin=laying_date.MEAN-laying_date.SE, ymax=laying_date.MEAN+laying_date.SE, group=morph),  
                color="black", width=0) +
  geom_point(aes(x=year, y=laying_date.MEAN, fill=morph), size=2, pch=21) +
  geom_line(data=data_pred_gamS2, aes(x=year, y=estimate, color=morph), size=1, linetype=1) +
  geom_line(data=data_pred_lm2_summ, aes(x=year, y=estimate), size=1, linetype=2) +
  scale_fill_manual(values=morph_cols) + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size=12),
    legend.position = "none",
    plot.margin = unit(c(5.5, 15, 5.5, 5.5), "pt"))

p4 =
ggplot(data_summ1) + ggtitle("D") +
  ylab("Recruitment probability") + 
  annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
  geom_errorbar(aes(x=year, ymin=recr_success.MEAN-recr_success.SE, ymax=recr_success.MEAN+recr_success.SE, group=morph),  
                color="black", width=0) +
  geom_point(aes(x=year, y=recr_success.MEAN, fill=morph), size=2, pch=21) +
  geom_ribbon(data=data_pred_gamS3, aes(x=year, ymin=conf.low, ymax=conf.high, group=morph), alpha=0.2) +
  geom_line(data=data_pred_gamS3, aes(x=year, y=estimate, color=morph), size=1, linetype=1) +
  geom_line(data=data_pred_lm3_summ, aes(x=year, y=estimate), size=1, linetype=2) +
  scale_color_manual(values=morph_cols) +
  scale_fill_manual(values=morph_cols) +
  scale_linetype_manual("Breeding",values=c(1,2)) + 
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12))

p5 =
  ggplot(data_env) + ggtitle("E") +
    ylab("Air temperature\n(°C, mean)") + 
    annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
    geom_line(aes(x=year, y=winter_temp_air_mean), linewidth=0.5) +
    geom_smooth(aes(x=year, y=winter_temp_air_mean), method="lm", color="black", linetype=2) +
    scale_y_continuous(limits=c(-6,2)) +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_blank(),
          plot.margin = unit(c(5.5, 5.5, 5.5, 20), "pt"))

p6 =
  ggplot(data_env) + ggtitle("F") +
    ylab("Air temperature\n(°C, CV)") + 
    annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
    geom_line(aes(x=year, y=winter_temp_air_CV), linewidth=0.5) +
    geom_smooth(aes(x=year, y=winter_temp_air_CV), method="lm", color="black", linetype=2) +
    scale_y_continuous(limits=c(0.2,0.5)) +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_blank())

p7 =
  ggplot(data_env) + ggtitle("G") +
    annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
    ylab("Snow days\n(N)") + coord_cartesian(ylim=c(0,175)) + 
    geom_line(aes(x=year, y=winter_snow_days_n), linewidth=0.5) +
    geom_smooth(aes(x=year, y=winter_snow_days_n), method="lm", color="black", linetype=2) + 
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_blank()
      )

p8 =
  ggplot(data_env) + ggtitle("H") +
    ylab("Snow depth\n(cm, mean)") + 
    annotate("rect", xmin=2010, xmax=2013, ymin=-Inf, ymax=Inf, alpha=0.3, fill="red") + 
    geom_line(aes(x=year, y=winter_snow_depth_mean), linewidth=0.5) +
    geom_smooth(aes(x=year, y=winter_snow_depth_mean), method="lm", color="black", linetype=2) + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size=12)
          )

p_grid1 = plot_grid(p1, p2, ncol=1, rel_heights = c(0.75,0.25), align="hv")
p_grid2 = plot_grid(p3, p4, ncol=1, align="hv")
p_grid3 = plot_grid(p5, p6, p7, p8, ncol=1, align="hv")
p_grid_all = plot_grid(p_grid1, p_grid2, p_grid3, 
                       ncol=3, align="h", axis="b", rel_widths = c(0.365,0.365,0.275))
p = plot_grid(p_grid_all, leg, ncol=1, rel_heights = c(0.9, 0.1))

ggsave("figures/figure2.jpg", p,
       width = 32, height = 22, units = 'cm',bg="white")


# figure 3 stats 1 (parent-offspring regression)  -------------------------

## only complete cases for parent-offspring regression
data_recr_sub = data_recr[complete.cases(color_score_median, colour_score_median_midparent)]

## parent-offspring regression
por1 = with(data_recr_sub, gls(color_score_median ~ colour_score_median_midparent))
capture.output(summary(por1), file = "tables/POR1.txt")

## add residuals from parent-offspring reg and n per level
data_recr_sub[as.numeric(names(residuals(por1))), res_por1 := residuals(por1)]
data_recr_sub[, n := .N, by=c("parent_combination", "morph")]
data_recr_sub[, n_ypos := median(res_por1), by=c("parent_combination", "morph")]

## residuals over time
data_recr_sub[, year_index := year - min(year) + 1]
por1_gls = gls(res_por1 ~ year_index * morph, data = data_recr_sub[!is.na(res_por1)],
               correlation = corAR1())
capture.output(summary(por1_gls), file = "tables/GLS_POR1.txt")


# figure 3 stats 2 (animal model)  ---------------------------------------

## additive variance component as random effect
## -> baseline heritability under a purely additive-only model.
am1 <- mmes(
  color_score_median ~ 1,
  random=~vsm(ism(ring),Gu=A) + vsm(ism(nest_id)),
  rcov=~units,
  data=data_animal_mod
)
summary(am1)
VC <- summary(am1)$varcomp
Va <- VC["ring", "VarComp"]
Vec <- VC["nest_id", "VarComp"]
Ve <- VC["units", "VarComp"]
h2_1 <- Va / (Va + Vec + Ve)
h2_1
ce2_1 <- Vec / (Va + Vec + Ve)
ce2_1

capture.output(summary(am1), file = "tables/AM1.txt")

## additive variance component and nest id as random effects and fixed effects for additive and dominance,
## and four key environmental variables                                                                            
## genotype-based effect size of dominance (interpretable effect of heterozygotes).
am2 <- mmes(
  color_score_median ~ geno_add + geno_dom + breed_snow_depth_mean + breed_snow_days_n + breed_temp_air_mean + breed_temp_air_CV,
  random=~vsm(ism(ring),Gu=A) + vsm(ism(nest_id)),
  rcov=~units,
  data=data_animal_mod
)

summary(am2)
VC <- summary(am2)$varcomp
Va <- VC["ring", "VarComp"]
Vec <- VC["nest_id", "VarComp"]
Ve <- VC["units", "VarComp"]
h2_2 <- Va / (Va + Vec + Ve)
h2_2
ce2_2 <- Vec / (Va + Vec + Ve)
ce2_2
capture.output(summary(am2), file = "tables/AM2.txt")

## evaluation
AIC1    <- summary(am1)$logo["Value", "AIC"]
AIC2    <- summary(am2)$logo["Value", "AIC"]
AIC_table <- data.frame(
  Model = c("Additive (AM 1)", 
            "Additive + genotype add./dom. (AM 2)"),
  AIC = c(AIC1, AIC2),
  h2 = c(h2_1, h2_2)
  
)
AIC_table

## residuals over time
data_animal_mod[, year_index := year - min(year) + 1]
data_animal_mod[as.numeric(names(residuals(am2)[,1])), res_am2 := residuals(am2)[,1]]
am2_gls = gls(res_am2 ~ year_index * morph, data = data_animal_mod[complete.cases(res_am2)],
              correlation = corAR1())
summary(am2_gls)
capture.output(summary(am2_gls), file = "tables/GLS_AM2.txt")

# figure 3 ---------------------------------------------------------------

labels = c("gray-gray (82)", "gray-brown (41)","brown-gray (28)","brown-brown (19)")
title = "Parental morph combination (female-male, N)"
p1=
  ggplot(data_recr_sub, aes(x=colour_score_median_midparent, y=color_score_median, 
                            color=parent_combination, shape=parent_combination,
                            fill=parent_combination)) +  ggtitle("A") +
  coord_fixed(ylim=c(4,14), xlim=c(4,14)) +
  scale_x_continuous(breaks=c(seq(4,14,2))) + scale_y_continuous(breaks=c(seq(4,14,2))) +
  ylab("Recruit color score") + xlab("Midparent color score") +
  geom_jitter(size=3, width=0.25, height=0.25, stroke =2) +
  geom_hline(yintercept=9, linetype=2, linewidth=0.5) + 
  geom_vline(xintercept=9, linetype=2, linewidth=0.5) + 
  geom_abline(slope = 1, intercept = 0) +
  geom_abline(slope = coef(por1)[[2]], intercept = coef(por1)[[1]], size=1.2, linetype=2) +
  guides(color = guide_legend(title.position="top", ncol=2)) + 
  scale_color_manual(title,labels = labels,
                     values=c( "gray", "gray","chocolate4", "chocolate4")) +
  scale_fill_manual(title, labels = labels,
                    values=c("gray", "gray", "chocolate4", "chocolate4")) + 
  scale_shape_manual(title,labels = labels,
                     values=c(16,24,25,16)) + 
  theme(legend.position = "bottom", 
        legend.justification = "center")

p2=
  ggplot(data_recr_sub) + ggtitle("B") +
  geom_point(aes(x=year, y=res_por1, color=morph)) + ylab("Residuals parent-offspring regr.") +
  geom_hline(yintercept=0, linetype=2) +
  geom_smooth(aes(x=year, y=res_por1, color=morph),method="lm") +
  guides(color = guide_legend(title.position="top")) + 
  scale_y_continuous(breaks=seq(-8,6,2), limits=c(-8,6)) +
  scale_color_manual("Recruit morph",
                     values=morph_cols, 
                     labels=c("gray", "brown")) +
  theme(
    legend.position = c(0.5,0.1),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.margin = margin(t = 2, r =2 , b = 2, l = 2, unit = "pt"),
    legend.box.background = element_rect(colour = "black"),
    # plot.margin = margin(t = 200, r =2 , b = -200, l = 2, unit = "pt"),
    axis.title.x = element_blank()
  )

ggsave("figures/figure3a.png",p1,width = 15, height = 15, units = 'cm', bg="white")
ggsave("figures/figure3b.png",p2,width = 10, height = 15, units = 'cm', bg="white")


# figure 4a stats ----------------------------------------------------------

data = copy(data_ind)
data = data[!is.na(laying_date),]
data[, laying_date_outlier := outlier_mad(laying_date, threshold = 3), by=morph]
data = data[laying_date_outlier == F,]
data[, morph_num := ifelse(morph=="gray",0,1)]

scale_vals = c(
  "year",
  "winter_temp_air_mean",
  "winter_temp_air_CV",
  "winter_snow_depth_mean",
  "winter_snow_days_n",
  "laying_date",
  "color_score_median_trans_ord"
)
scale_factors = data[, lapply(.SD, function(x) list(mean = mean(x, na.rm=T), sd = sd(x, na.rm=T))), .SDcols = scale_vals]
data_scaled = copy(data)
data_scaled[, (scale_vals) := lapply(.SD, scale), .SDcols=scale_vals]

H1 <- "
  winter_temp_air_mean ~ year
  winter_temp_air_CV ~ year
  winter_snow_depth_mean ~ year
  winter_snow_days_n ~ year
  color_score_median_trans_ord ~ winter_temp_air_mean + winter_temp_air_CV + winter_snow_depth_mean + winter_snow_days_n + morph_num
  laying_date ~ winter_temp_air_mean + winter_temp_air_CV + winter_snow_depth_mean + winter_snow_days_n + color_score_median_trans_ord + morph_num
  recr_success ~ winter_temp_air_mean + winter_temp_air_CV + winter_snow_depth_mean + winter_snow_days_n + color_score_median_trans_ord + laying_date + morph_num
"

fitH1 <- sem(H1, data = data_scaled, estimator = "MLR")
summary(fitH1)

## bootstrapped version to confirm single run results
# fitH1 <- sem(H1, data = data_scaled, se = "bootstrap", test = "bootstrap", bootstrap = 5000,
# parallel ="snow", ncpus = 8, verbose=T)

data_results = data.table(standardizedsolution(fitH1), parameterestimates(fitH1)["est"])
data_results = data_results[op == "~", .(lhs, rhs, est, est.std, se, z, pvalue)]
fwrite(data_results, "tables/SEM1.txt")


# figure 4a ---------------------------------------------------------------

# Extract path estimates
mod_est <- data.table(standardizedSolution(fitH1))
mod_est <- mod_est[op == "~"]

# Define labels
labels <- c(
  "year" = "Year",
  "winter_temp_air_mean" = "Air temp.\n(Mean)",
  "winter_temp_air_CV" = "Air temp.\n(CV)",
  "winter_snow_depth_mean" = "Snow depth\n(Mean)",
  "winter_snow_days_n" = "Snow days\n(N)",
  "morph_num" = "Morph",
  "color_score_median_trans_ord" = "Color Score",
  "laying_date" = "Laying Date",
  "recr_success" = "Recruitment\nSuccess"
)

calculate_dimensions <- function(labels) {
  dimensions_list <- lapply(labels, function(label) {
    lines <- strsplit(label, "\n")[[1]]
    max_width <- (max(nchar(lines)) / 10) + 0.75
    height <- length(lines) * 0.45
    return(c(width = max_width, height = height))
  })
  dimensions_matrix <- do.call(rbind, dimensions_list)
  return(data.table(dimensions_matrix))
}
nodes = data.table(labels, names(labels))
setnames(nodes, c("V2", "labels"), c("node_name", "label_name"))
nodes[, c("width", "height") := calculate_dimensions(label_name) ]

# Prepare data for qgraph
sig_lvl = 0.05
mod_est[, independent:= labels[rhs]]
mod_est[, dependent := labels[lhs]]
mod_est[, colour := ifelse(pvalue > sig_lvl, "gray80", ifelse(est.std < 0, "red", "green"))]
mod_est[, edge.width := ifelse(pvalue < sig_lvl, abs(est.std) * 20, 1)]
mod_est[, edge.label := ifelse(pvalue < sig_lvl, round(est.std, 2), NA)]
mod_est[,"curve"] = 0
curve_val = 1.4
mod_est[rhs=="winter_temp_air_mean" & !lhs=="color_score_median_trans_ord", curve := curve_val * -1]
mod_est[rhs=="winter_temp_air_CV" & !lhs=="color_score_median_trans_ord", curve := curve_val * -1]
mod_est[rhs=="winter_snow_depth_mean" & !lhs=="color_score_median_trans_ord", curve := curve_val]
mod_est[rhs=="winter_snow_days_n" & !lhs=="color_score_median_trans_ord", curve := curve_val]
mod_est[rhs=="color_score_median_trans_ord" & lhs=="recr_success", curve := -2]
mod_est[, mar := ifelse(pvalue < sig_lvl, 0.01, 0)]
mod_est[, asize := 3]
mod_est[pvalue < sig_lvl, asize := pmax(abs(est.std) * 10, 5)]

png(paste0("figures/figure4a.png"), width = 30, height = 15, units = 'cm', res = 400)
qgraph(input=mod_est[,c("independent", "dependent")],
       shape="rectangle",
       layout=rbind(
         c(5, 4),   # year
         c(2, 3),   # winter_temp_air_mean
         c(4, 3),   # winter_temp_air_CV
         c(6, 3),   # winter_snow_depth_mean
         c(8, 3),   # winter_snow_days_n
         c(8, 0),   # morph
         c(5, 2),   # laying_date
         c(5, 1),   # color_score_median_trans_ord
         c(5, 0)    # recr_success
       ), 
       # title=morph_name,
       label.scale=F,
       label.cex=1.125, 
       label.scale.equal=TRUE,
       node.width=nodes$width,
       node.height=nodes$height,
       
       # node.width=c(1:6),
       edge.labels = mod_est$edge.label,
       edge.label.margin=mod_est$mar,
       # edge.label.position=mod_est$edge.lab_pos,
       # edgeConnectPoints=as.matrix(edges[,c("con_point1","con_point2")]),
       edge.label.cex = 0.75,
       edge.color=mod_est$colour,
       edge.width=1, 
       
       esize = mod_est$edge.width,
       
       directed = TRUE,
       arrows =TRUE, 
       asize=mod_est$asize,
       curve = mod_est$curve,
       
       curvePivot=0.4,
       color="gray80",
       mar=c(3,3,3,3)
       
)
dev.off()


# figure 4b-d stats --------------------------------------------------------------

data = copy(data_ind)
data = data[!is.na(laying_date),]
data[, laying_date_outlier := outlier_mad(laying_date, threshold = 3), by=morph]
data = data[laying_date_outlier == F,]

knots = 5
gam2 = gam(color_score_median_trans_ord ~
                   morph + 
                   ti(winter_temp_air_mean, by=morph, k=knots) +
                   ti(winter_snow_days_n, by=morph, k=knots) +
                   ti(winter_temp_air_mean, winter_snow_days_n, by=morph, k=knots),
                 method = "REML",
                 data=data)
summary(gam2)
anova(gam2)
# gratia::draw(gam2)
capture.output(anova(gam2), file = "tables/GAM2.txt")


knots = 5
gam3 = gam(laying_date ~
               morph + 
               ti(winter_temp_air_mean, by=morph, k=knots) +
               ti(color_score_median_trans_ord, by=morph, k=knots) +
               ti(winter_temp_air_mean, color_score_median_trans_ord, by=morph, k=knots),
             method = "REML",
             data=data)

summary(gam3)
anova(gam3)
# gratia::draw(gam3)
capture.output(anova(gam3), file = "tables/GAM3.txt")

knots = 5
gam4 = gam(recr_success ~ morph + 
               ti(laying_date, by=morph, k=knots) +
               ti(color_score_median_trans_ord, by=morph, k=knots) +
               ti(laying_date, color_score_median_trans_ord, by=morph, k=knots),
             family="binomial",
             method="REML",
             data=data)
summary(gam4)
anova(gam4)
# gratia::draw(gam4)
capture.output(anova(gam4), file = "tables/GAM4.txt")



# figure 4b ------------------------------------------------------------------

## model name
mod = copy(gam2)

## variables
dep_var = "color_score_median_trans_ord"
vary = "winter_temp_air_mean"
varx = "winter_snow_days_n"

## categorical variable
var_cat = "morph"

## give extra prediction-range, in percent
extra_range_x_left = 0
extra_range_x_right = 0
extra_range_y_top = 0
extra_range_y_bottom = 0

## 3d props
n_gridlines = 25
too_far_val = 0.333
theta_val = 315 - 270
phi_val = 35
distance = 1000
alpha_val = 1
point_size = 1.25

# generate range of mod data
data_mod = data.table(copy(mod$model[, c(vary, varx, dep_var, "morph")]))
setnames(data_mod, c(varx, vary), c("varx", "vary"))
vars = list(seq(min(data_mod$varx) - (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_left, 
                max(data_mod$varx) + (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_right, length.out = n_gridlines), 
            seq(min(data_mod$vary) - (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_bottom, 
                max(data_mod$vary) + (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_top, length.out = n_gridlines),
            c("gray","brown"))
data_new = data.table(expand.grid(vars))
setnames(data_new, c(varx, vary, var_cat))
data_new[, predicted := predict.gam(mod, data_new, type= "response")]

## reduce to data-range
for(i in unique(data_mod[[var_cat]])){
  too_far_res = exclude.too.far(
    data_new[c(data_new[,..var_cat]==i),][[varx]],
    data_new[c(data_new[,..var_cat]==i),][[vary]],
    data[c(data[,..var_cat]==i),][[varx]], 
    data[c(data[,..var_cat]==i),][[vary]],
    too_far_val)
  data_new[c(data_new[,..var_cat]==i), too_far := too_far_res
  ]
}
# 
data_new[too_far==FALSE & predicted <= 6, predicted_filter := predicted]


## figure panel
dims = 2500
png(paste0("figures/figure4b.png"), width=dims, height=dims, unit="px", res=400, bg="white")
par(mfrow=c(1,1), mar = c(2,2,2,1))
z_orig = data_mod[[dep_var]]
x_orig = data_mod$varx
y_orig = data_mod$vary 

xlow =  min(data_mod$varx)
xupp =  max(data_mod$varx)
ylow =  min(data_mod$vary)
yupp =  max(data_mod$vary)

bg_cols = morph_cols[data_mod$morph]

scatter3D(x=x_orig, y=y_orig, z=z_orig, 
          ticktype = "detailed", 
          cex.main=2, colkey=FALSE,font.main = 1,
          col="black",   pch=21, lwd=0, bg = bg_cols, cex=0, cex.lab=1.5, 
          bty = "b2", box=T, #type = "h",
          phi = phi_val, theta = theta_val, r=distance, 
          zlim=c(0,5), xlim=c(0,180), ylim=c(-4,2),
          xlab="Snow days (N)",
          ylab="Air temp. (Mean)",
          zlab="Color score")

data_sub = data_new[morph=="brown",]
z = data_sub$predicted_filter
x = data_sub[[varx]]
y = data_sub[[vary]]
surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols["brown"], NAcol=morph_cols["brown"],
       add=TRUE)

data_sub = data_new[morph=="gray",]
z = data_sub$predicted_filter
x = data_sub[[varx]]
y = data_sub[[vary]]
surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols["gray"], NAcol=morph_cols["gray"],
       add=TRUE)

dev.off()

# figure 4c ------------------------------------------------------------------

## model name
mod = copy(gam3)

## variables
dep_var = "laying_date"
varx = "winter_temp_air_mean"
vary = "color_score_median_trans_ord"

## categorical variable
var_cat = "morph"

## give extra prediction-range, in percent
extra_range_x_left = 0
extra_range_x_right = 0
extra_range_y_top = 0
extra_range_y_bottom = 0

## 3d props
n_gridlines = 25
too_far_val = 0.333
theta_val = 315 - 270
phi_val = 35
distance = 1000
alpha_val = 1
point_size = 1.25

# generate range of mod data
data_mod = data.table(copy(mod$model[, c(vary, varx, dep_var, "morph")]))
setnames(data_mod, c(varx, vary), c("varx", "vary"))
vars = list(seq(min(data_mod$varx) - (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_left, 
                max(data_mod$varx) + (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_right, length.out = n_gridlines), 
            seq(min(data_mod$vary) - (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_bottom, 
                max(data_mod$vary) + (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_top, length.out = n_gridlines),
            c("gray","brown"))
data_new = data.table(expand.grid(vars))
setnames(data_new, c(varx, vary, var_cat))
data_new[, predicted := predict.gam(mod, data_new, type= "response")]

## reduce to data-range
for(i in unique(data_mod[[var_cat]])){
  too_far_res = exclude.too.far(
    data_new[c(data_new[,..var_cat]==i),][[varx]],
    data_new[c(data_new[,..var_cat]==i),][[vary]],
    data[c(data[,..var_cat]==i),][[varx]], 
    data[c(data[,..var_cat]==i),][[vary]],
    too_far_val)
  data_new[c(data_new[,..var_cat]==i), too_far := too_far_res
  ]
}
# & predicted <= 5
data_new[too_far==FALSE, predicted_filter := predicted]


## figure panel
dims = 2500
png(paste0("figures/figure4c.png"), width=dims, height=dims, unit="px", res=400, bg="white")
par(mfrow=c(1,1), mar = c(2,2,2,1))
z_orig = data_mod[[dep_var]]
x_orig = data_mod$varx
y_orig = data_mod$vary 

xlow =  min(data_mod$varx)
xupp =  max(data_mod$varx)
ylow =  min(data_mod$vary)
yupp =  max(data_mod$vary)

bg_cols = morph_cols[data_mod$morph]

scatter3D(x=x_orig, y=y_orig, z=z_orig, 
          ticktype = "detailed", 
          # main = "All years", 
          cex.main=2, colkey=FALSE,font.main = 1,
          col="black",   pch=21, lwd=0, bg = bg_cols, cex=0, cex.lab=1.5, 
          bty = "b2", box=T, #type = "h",
          phi = phi_val, theta = theta_val, r=distance, 
          zlim=c(-20,20),
          zlab="Laying date",
          xlab="Air temp. (Mean)",
          ylab="Color score")

data_sub = data_new[morph=="brown",]
z = data_sub$predicted_filter
x = data_sub[[varx]]
y = data_sub[[vary]]
surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols["brown"], NAcol=morph_cols["brown"],
       add=TRUE)

data_sub = data_new[morph=="gray",]
z = data_sub$predicted_filter
x = data_sub[[varx]]
y = data_sub[[vary]]
surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols["gray"], NAcol=morph_cols["gray"],
       add=TRUE)

dev.off()

# figure 4d ------------------------------------------------------------------

## model name
mod = copy(gam4)

## variables
dep_var = "recr_success"
vary = "color_score_median_trans_ord"
varx = "laying_date"

## categorical variable
var_cat = "morph"

## give extra prediction-range, in percent
extra_range_x_left = 0
extra_range_x_right = 0
extra_range_y_top = 0
extra_range_y_bottom = 0

## 3d props
n_gridlines = 25
too_far_val = 0.333
theta_val = 315 - 270
phi_val = 35
distance = 1000
alpha_val = 1
point_size = 1.25

# generate range of mod data
data_mod = data.table(copy(mod$model[, c(vary, varx, dep_var, "morph")]))
setnames(data_mod, c(varx, vary), c("varx", "vary"))
vars = list(seq(min(data_mod$varx) - (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_left, 
                max(data_mod$varx) + (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_right, length.out = n_gridlines), 
            seq(min(data_mod$vary) - (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_bottom, 
                max(data_mod$vary) + (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_top, length.out = n_gridlines),
            c("gray","brown"))
data_new = data.table(expand.grid(vars))
setnames(data_new, c(varx, vary, var_cat))
data_new[, predicted := predict.gam(mod, data_new, type= "response")]

## reduce to data-range
for(i in unique(data_mod[[var_cat]])){
  too_far_res = exclude.too.far(
    data_new[c(data_new[,..var_cat]==i),][[varx]],
    data_new[c(data_new[,..var_cat]==i),][[vary]],
    data[c(data[,..var_cat]==i),][[varx]], 
    data[c(data[,..var_cat]==i),][[vary]],
    too_far_val)
  data_new[c(data_new[,..var_cat]==i), too_far := too_far_res
  ]
}
# & predicted <= 5
data_new[too_far==FALSE, predicted_filter := predicted]


## figure panel
dims = 2500
png(paste0("figures/figure4d.png"), width=dims, height=dims, unit="px", res=400, bg="white")
par(mfrow=c(1,1), mar = c(2,2,2,1))
z_orig = data_mod[[dep_var]]
x_orig = data_mod$varx
y_orig = data_mod$vary 

xlow =  min(data_mod$varx)
xupp =  max(data_mod$varx)
ylow =  min(data_mod$vary)
yupp =  max(data_mod$vary)

bg_cols = morph_cols[data_mod$morph]

scatter3D(x=x_orig, y=y_orig, z=z_orig, 
          ticktype = "detailed", 
          # main = "All years", 
          cex.main=2, colkey=FALSE,font.main = 1,
          col="black",   pch=21, lwd=0, bg = bg_cols, cex=0, cex.lab=1.5, 
          bty = "b2", box=T, #type = "h",
          phi = phi_val, theta = theta_val, r=distance, 
          zlim=c(0,0.5),
          zlab="\nRecruitment success (predicted)",
          xlab="Air temp. (Mean)",
          ylab="Color score")
# zlab="\nRecruitment success (predicted)",
# xlab="Laying date",
# ylab="Color score")

data_sub = data_new[morph=="brown",]
z = data_sub$predicted_filter
x = data_sub[[varx]]
y = data_sub[[vary]]
surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols["brown"], NAcol=morph_cols["brown"],
       add=TRUE)

data_sub = data_new[morph=="gray",]
z = data_sub$predicted_filter
x = data_sub[[varx]]
y = data_sub[[vary]]
surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols["gray"], NAcol=morph_cols["gray"],
       add=TRUE)

dev.off()

# figure 5a stats ---------------------------------------------------------


data = copy(data_recr)
data[, morph_m := ifelse(morph_m=="gray",0,1)]
data[, morph_f := ifelse(morph_f=="gray",0,1)]
data[, morph := ifelse(morph=="gray",0,1)]

scale_vals = c(
  "color_score_median_trans_ord_f",
  "color_score_median_trans_ord_m",
  "breed_snow_depth_mean", 
  "breed_snow_days_n", 
  "breed_temp_air_mean", 
  "breed_temp_air_CV",
  "color_score_median_trans_ord",
  "laying_date"
)


data[, (scale_vals) := lapply(.SD, scale), .SDcols=scale_vals]

H2 <- "
  color_score_median_trans_ord_m ~ morph_m 
  color_score_median_trans_ord_f ~ morph_f
  laying_date ~ color_score_median_trans_ord_m + color_score_median_trans_ord_f + morph_m + morph_f
  breed_temp_air_mean ~ laying_date
  breed_temp_air_CV ~ laying_date
  breed_snow_depth_mean ~ laying_date
  breed_snow_days_n ~ laying_date
  morph ~ morph_m + morph_f 
  color_score_median_trans_ord ~ morph + color_score_median_trans_ord_m + color_score_median_trans_ord_f + breed_temp_air_mean + breed_temp_air_CV + breed_snow_depth_mean + breed_snow_days_n
  "


fitH2 <- sem(H2, data = data)
# fitH2 <- sem(H2, data = data_scaled, se = "bootstrap", test = "bootstrap", bootstrap = 5000,
#              parallel ="snow", ncpus = 8, verbose=T)
data_results = data.table(standardizedsolution(fitH2), parameterestimates(fitH2)["est"])
data_results = data_results[op == "~", .(lhs, rhs, est, est.std, se, z, pvalue)]

# lavaanPlot(model = fitH2,
#            coefs = TRUE,
#            sig = .1,
#            # coef_labels = TRUE,
#            stars = c("regress"),
#            # edge_options = e_opts,
#            node_options = list( fontname = "Helvetica")
# )

fwrite(data_results, "tables/SEM2.txt")
# write(paste(utils::capture.output(data.frame(data_results)),
# collapse = "\n"), file = "tables/SEM1.txt")

# figure 5a -------------------------------------------------------------------

# Extract path estimates
mod_est <- data.table(parameterEstimates(fitH2, standardized = TRUE))
mod_est <- mod_est[op == "~"]

labels = c("morph_m" = "Paternal morph",
           "morph_f" = "Maternal morph",
           "color_score_median_trans_ord_f" = "Maternal color score",
           "color_score_median_trans_ord_m" = "Paternal color score",
           "laying_date" = 'Laying date',
           "morph" = "Recruit morph",
           "breed_snow_depth_mean" = "Snow depth\nbreeding (Mean)",
           "breed_snow_days_n" = "Snow days\nbreeding (N)",
           "breed_temp_air_mean" = "Air temp.\nbreeding (Mean)",
           "breed_temp_air_CV" = "Air temp.\nbreeding (CV)",
           "color_score_median_trans_ord" = "Recruit color score"
)

calculate_dimensions <- function(labels) {
  dimensions_list <- lapply(labels, function(label) {
    lines <- strsplit(label, "\n")[[1]]
    max_width <- (max(nchar(lines)) / 10) + 0.75
    height <- length(lines) * 0.45
    return(c(width = max_width, height = height))
  })
  dimensions_matrix <- do.call(rbind, dimensions_list)
  return(data.table(dimensions_matrix))
}

nodes = data.table(labels, names(labels))
setnames(nodes, c("V2", "labels"), c("node_name", "label_name"))
nodes[, c("width", "height") := calculate_dimensions(label_name) ]

# Prepare data for qgraph
sig_lvl = 0.05
mod_est[, independent:= labels[lhs]]
mod_est[, dependent := labels[rhs]]
mod_est[, colour := ifelse(pvalue > sig_lvl, "gray80", ifelse(std.all < 0, "red", "green"))]
mod_est[, edge.width := ifelse(pvalue < sig_lvl, abs(std.all) * 20, 1)]
mod_est[, edge.label := ifelse(pvalue < sig_lvl, round(std.all, 2), NA)]
mod_est[, curve := 0]
mod_est[rhs == "color_score_median_trans_ord_m" & lhs == "color_score_median_trans_ord", curve := -2.25]
mod_est[rhs=="color_score_median_trans_ord_f" & lhs=="color_score_median_trans_ord", curve := 2.25]
mod_est[rhs=="morph" & lhs=="color_score_median_trans_ord", curve := -1.65]
mod_est[, mar := ifelse(pvalue < sig_lvl, 0.01, 0)]
mod_est[, asize := 3]
mod_est[pvalue < sig_lvl, asize := pmax(abs(std.all) * 10, 5)]

png("figures/figure5a.png", width = 20, height = 12.5, units = 'cm', res = 400)
qgraph(input=mod_est[,c("dependent", "independent")],
       shape="rectangle",
       layout=rbind(c(3,8), #  paternal morph
                    c(6,8), #  maternal morph
                    c(1,6), #  paternal color score
                    c(8,6), #  maternal color score
                    c(5,4), # laying date 
                    c(3,6), # recruit morph 
                    c(2,1.5), # mean snow depth breeding
                    c(3.75,1.5), # mean air temp breeding
                    c(5.25,1.5), # cv air temp breeding
                    c(7,1.5), # n snow days
                    c(4.5,-1)  # recruit color score
       ) , 
       # title=morph_name,
       label.scale=F,
       label.cex=0.9, 
       label.scale.equal=TRUE,
       node.width=nodes$width,
       node.height=nodes$height,
       
       # node.width=c(1:6),
       edge.labels = mod_est$edge.label,
       edge.label.margin=mod_est$mar,
       # edge.label.position=mod_est$edge.lab_pos,
       # edgeConnectPoints=as.matrix(edges[,c("con_point1","con_point2")]),
       edge.label.cex = 0.75,
       edge.color=mod_est$colour,
       edge.width=1, 
       
       esize = mod_est$edge.width,
       
       directed = TRUE,
       arrows =TRUE, 
       asize=mod_est$asize,
       curve = mod_est$curve,
       
       curvePivot=0.4,
       color="gray80",
       mar=c(3,3,3,3)
       
)
dev.off()


# figure 5b stats -------------------------------------------------------

knots = 5
gam5 = gam(color_score_median_trans_ord ~  morph + 
              ti(color_score_median_trans_ord_m, by=morph, k=knots) +
              ti(color_score_median_trans_ord_f, by=morph, k=knots) +
              ti(color_score_median_trans_ord_m, color_score_median_trans_ord_f, by=morph, k=knots),
            method = "REML",
            data=data_recr)
summary(gam5)
anova(gam5)
capture.output(anova(gam5), file = "tables/GAM5.txt")


# figure 5b ---------------------------------------------------------------

## model
mod = gam5

## give extra prediction-range, in percent
extra_range_x_left = 0
extra_range_x_right = 0
extra_range_y_top = 0
extra_range_y_bottom = 0

## variables
dependent_var = "color_score_median_trans_ord"
varx = "color_score_median_trans_ord_m"
vary = "color_score_median_trans_ord_f"
var_cat = "morph"

## labels
ylab="\nMaternal color score"
xlab="\nPaternal color score"
zlab="\nRecrut color score\n(predicted)"

## 3d props
n_gridlines = 25
too_far_val = 0.333
theta_val = 315
phi_val = 35
distance = 1000
alpha_val = 1
point_size = 0

# generate range of mod data
data_mod = data.table(mod$model)
setnames(data_mod, c(varx, vary), c("varx", "vary"))
vars = list(seq(min(data_mod$varx) - (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_left, 
                max(data_mod$varx) + (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_right, length.out = n_gridlines), 
            seq(min(data_mod$vary) - (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_bottom, 
                max(data_mod$vary) + (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_top, length.out = n_gridlines),
            c("gray", "brown"))
data_new = data.table(expand.grid(vars))
setnames(data_new, c(varx, vary, var_cat))
data_new[, predicted := predict(mod, data_new, type= "response")]
for(i in unique(data_mod[[var_cat]])){
  too_far_res = exclude.too.far(
    data_new[c(data_new[,..var_cat]==i),][[varx]],
    data_new[c(data_new[,..var_cat]==i),][[vary]],
    data_mod[c(data_mod[,..var_cat]==i),]$varx, 
    data_mod[c(data_mod[,..var_cat]==i),]$vary,
    too_far_val)
  data_new[c(data_new[,..var_cat]==i), too_far := too_far_res
  ]
}
data_new[too_far==FALSE, predicted_filter := predicted]

## figure panel
z_orig = data_mod[[dependent_var]]
x_orig = data_mod$varx
y_orig = data_mod$vary

bg_cols = morph_cols[as.character(data_mod$morph)]
dims = 2500
png(paste0("figures/figure5b.png"), width=dims, height=dims, unit="px", res=400, bg="white")
par(mar=c(3,3,3,3))
scatter3D(x=x_orig , y=y_orig, z=z_orig,
          ticktype = "detailed", 
          cex.main=2, colkey=FALSE, font.main = 1,
          col="black",   pch=21, lwd=0, bg = bg_cols, cex=0, cex.lab=1.5, 
          bty = "b2", box=T, phi = phi_val, theta = theta_val, pch=19, r=distance, 
          zlim=c(min(z_orig), max(z_orig)),
          xlim=c(min(x_orig), max(x_orig)),
          ylim=c(min(y_orig), max(y_orig)),
          ylab=ylab,
          xlab=xlab,
          zlab=zlab
)

morph_name = "brown"
z = data_new[morph==morph_name,]$predicted_filter
x = data_new[morph==morph_name,][[varx]]
y = data_new[morph==morph_name,][[vary]]

surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols[morph_name], NAcol=morph_cols[morph_name],
       add=TRUE)

morph_name = "gray"
z = data_new[morph==morph_name,]$predicted_filter
x = data_new[morph==morph_name,][[varx]]
y = data_new[morph==morph_name,][[vary]]

surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols[morph_name], NAcol=morph_cols[morph_name],
       add=TRUE)

dev.off()

# figure 5c stats -------------------------------------------------------

knots = 5
gam6 = gam(color_score_median_trans_ord ~  
              morph + 
              ti(breed_temp_air_mean  , by=morph, k=knots) +
              ti(breed_temp_air_CV , by=morph, k=knots) +
              ti(breed_temp_air_mean  , breed_temp_air_CV , by=morph, k=knots),
            method = "REML",
            data=data_recr)
summary(gam6)
anova(gam6)
capture.output(anova(gam6), file = "tables/GAM6.txt")

# figure 5c -------------------------------------------------------------------

## model
mod = gam6 

## give extra prediction-range, in percent
extra_range_x_left = 0
extra_range_x_right = 0
extra_range_y_top = 0
extra_range_y_bottom = 0

## variables
dependent_var = "color_score_median_trans_ord"
varx = "breed_temp_air_mean"
vary = "breed_temp_air_CV"
var_cat = "morph"

## labels
ylab="\nAir temp.\nbreeding (CV)"
xlab="\nAir temp.\nbreeding (Mean)"
zlab="\nRecruit color score\n(predicted)"


## 3d props
n_gridlines = 25
too_far_val = 0.333
theta_val = 315 - 180
phi_val = 40
distance = 1000
alpha_val = 1

# generate range of mod data
data_mod = data.table(mod$model)
setnames(data_mod, c(varx, vary), c("varx", "vary"))
vars = list(seq(min(data_mod$varx) - (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_left, 
                max(data_mod$varx) + (max(abs(data_mod$varx))/n_gridlines)*extra_range_x_right, length.out = n_gridlines), 
            seq(min(data_mod$vary) - (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_bottom, 
                max(data_mod$vary) + (max(abs(data_mod$vary))/n_gridlines)*extra_range_y_top, length.out = n_gridlines),
            c("gray", "brown"))
data_new = data.table(expand.grid(vars))
setnames(data_new, c(varx, vary, var_cat))
data_new[, predicted := predict(mod, data_new, type= "response")]
for(i in unique(data_mod[[var_cat]])){
  too_far_res = exclude.too.far(
    data_new[c(data_new[,..var_cat]==i),][[varx]],
    data_new[c(data_new[,..var_cat]==i),][[vary]],
    data_mod[c(data_mod[,..var_cat]==i),]$varx, 
    data_mod[c(data_mod[,..var_cat]==i),]$vary,
    too_far_val)
  data_new[c(data_new[,..var_cat]==i), too_far := too_far_res
  ]
}
data_new[too_far==FALSE & predicted >= 1, predicted_filter := predicted]

## figure panel
z_orig = data_mod[[dependent_var]]
x_orig = data_mod$varx
y_orig = data_mod$vary

bg_cols = morph_cols[as.character(data_mod$morph)]
png(paste0("figures/figure5c.png"), width=2750, height=2750, unit="px", res=400, bg="white")
par(mar=c(3,3,3,3))
scatter3D(x=x_orig , y=y_orig, z=z_orig,
          ticktype = "detailed", 
          cex.main=2, colkey=FALSE, font.main = 1,
          col="black",   pch=21, lwd=0, bg = bg_cols, cex=0, cex.lab=1.5, 
          bty = "b2", box=T, phi = phi_val, theta = theta_val, pch=19, r=distance, 
          zlim=c(min(z_orig), max(z_orig)),
          xlim=c(min(x_orig), 17),
          # ylim=c(min(y_orig), 50),
          ylab=ylab,
          xlab=xlab,
          zlab=zlab
)

morph_name = "brown"
z = data_new[morph==morph_name,]$predicted_filter
x = data_new[morph==morph_name,][[varx]]
y = data_new[morph==morph_name,][[vary]]

surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols[morph_name], NAcol=morph_cols[morph_name],
       add=TRUE)

morph_name = "gray"
z = data_new[morph==morph_name,]$predicted_filter
x = data_new[morph==morph_name,][[varx]]
y = data_new[morph==morph_name,][[vary]]

surf3D(x=matrix(x, n_gridlines, n_gridlines, byrow=T),
       y=matrix(y, n_gridlines, n_gridlines, byrow=T),
       z=matrix(z, n_gridlines, n_gridlines, byrow=T),
       border="black", alpha=alpha_val,
       col=morph_cols[morph_name], NAcol=morph_cols[morph_name],
       add=TRUE)

dev.off()


# SUPPLEMENT --------------------------------------------------------------
# figure S1 ---------------------------------------------------------------

## was made prior to pre-processing: requires raw data - not included

# figure S2 ---------------------------------------------------------------

data = copy(data_ind)
data = data[order(year, laying_date)]
data[!is.na(laying_date), laying_date_pos := laying_date - min(laying_date) + 1]
data[, laying_date_seq := match(laying_date_pos, sort(unique(laying_date_pos))), by = year]
data[!is.na(laying_date), time := paste0(year,"-",laying_date_seq)]
data[, time := as.numeric(factor(time, levels=unique(time)))]
data[, year_index := year - min(year) + 1]
data[,ring:=factor(ring)]

## data summaries
data_summ1 = data[
  , as.list(unlist(lapply(.SD, m_se))), 
  by=c("morph", "year", 'year_index' ),
  .SDcols= c("color_score_median","color_score_median_trans_ord","laying_date", "recr_success")] 
data_summ1_obs = data_ind[, as.list(unlist(lapply(.SD, m_se))), 
                          by=c("morph", "year", "observer"),
                          .SDcols= c("color_score_median", "color_score_median_trans")] 
data_summ1_obs[, morph := factor(morph, levels=c("brown","gray"))]

## model predictions with observer effect
data_pred_gam1_obs <- data.table(predictions(gam1))
data_pred_gam1_obs[morph=="gray", estimate := estimate+3]
data_pred_gam1_obs[morph=="brown", estimate := estimate+9]
data[, year_index := year - min(year) + 1]
data_pred_gam1_obs[, mod_name := rep(c("GAM 1 (all observations)", "GAM S1a (immigrants only)", "GAM S1b (new cohort only)"), length.out = .N)]
data_pred_gam1_obs[, year := year_index + 1979]


plot_leg =
ggplot(data_summ1_obs) + ggtitle("A") +
  coord_cartesian(ylim = c(4,14)) +  ylab("Color score") +
  geom_errorbar(data=data_summ1, 
                aes(x=year, ymin=color_score_median.MEAN-color_score_median.SE, 
                    ymax=color_score_median.MEAN+color_score_median.SE, group=morph), 
                color="black", width=0) +
  geom_hline(yintercept=9, linetype=2) + 
  geom_point(data=data_summ1, aes(x=year, y=color_score_median.MEAN, fill=morph), size=4, pch=21) +
  geom_point(aes(x=year, y=color_score_median.MEAN, color=observer), size=2, pch=16) +
  geom_line(data=data_pred_gam1_obs, aes(x=year, y=estimate, color=observer), size=1) +
  geom_line(data=data_pred_gam1_obs, aes(x=year, y=estimate, linetype=mod_name), size=1) +
  scale_color_manual("Observer", values=scales::hue_pal()(3), labels=c("'a'", "'b'", "'c'")) +
  scale_fill_manual("Morph", values=morph_cols, labels=c("Gray", "Brown")) +
  scale_linetype_manual("Model fits", values=c(1,2,3)) +
  guides(
    fill=guide_legend(order=1),
    color=guide_legend(order=2),
    linetype=guide_legend(order=3, keywidth = unit(2, "cm"))) + 
  theme(legend.position="bottom",
        legend.text = element_text(size = 12),
        legend.justification = "center",
        legend.direction = "vertical",
        legend.key.width = unit(1.2,"cm"),
        legend.margin = margin(25,0,0,0),
        legend.spacing.x = unit(30, 'pt'),
        axis.title.x = element_blank())
leg1 = cowplot::get_plot_component(plot_leg, "guide-box", return_all = TRUE)[[3]]


p1 =
ggplot(data_summ1_obs) + ggtitle("A") +
  coord_cartesian(ylim = c(4,14)) +  ylab("Color score") +
  geom_errorbar(data=data_summ1, 
                aes(x=year, ymin=color_score_median.MEAN-color_score_median.SE, 
                    ymax=color_score_median.MEAN+color_score_median.SE, group=morph), 
                color="black", width=0) +
  geom_hline(yintercept=9, linetype=2) + 
  geom_point(data=data_summ1, aes(x=year, y=color_score_median.MEAN, fill=morph), size=3, pch=21) +
  geom_point(aes(x=year, y=color_score_median.MEAN, color=observer), size=2, pch=16) +
  geom_line(data=data_pred_gam1, aes(x=year, y=estimate, group=morph), size=1, linetype=1) +
  geom_ribbon(data=data_pred_gam1_obs, aes(x=year, ymin=estimate-std.error, ymax=estimate+std.error, group=interaction(observer, morph)),alpha=0.2) +
  geom_line(data=data_pred_gam1_obs, 
            aes(x=year, y=estimate, group=interaction(morph, observer), color=observer), size=1) +
  scale_y_continuous(breaks=c(4:14), labels= col_score_labs) +
  scale_x_continuous(expand = c(0.02, 0)) + 
  scale_fill_manual("Morphs", values=morph_cols, guide="none") +
  scale_color_manual("", values=scales::hue_pal()(3), labels=c("Observer 'a'", "Observer 'b'", "Observer 'c'")) +
  guides(color = guide_legend(title.position="top", nrow=1)) + 
  theme(legend.position="none",
        axis.title.x = element_blank())

data = copy(data_ind)
data[, migration := factor(ifelse(ring %in% data_recr$ring,T,F)),]
gamS1a = gam(color_score_median_trans_ord ~ morph  +
             s(year, by=morph, k=5) +
             s(observer, by=morph, bs = 're') ,
           data=data[migration==T])
summary(gamS1a)
anova(gamS1a)
capture.output(anova(gamS1a), file = "tables/GAMS1a.txt")

## get estimates without observer RE
data_pred_modS1a <- data.table(predictions(gamS1a, exclude=list(
  "s(observer):morphgray", "s(observer):morphbrown"), newdata=expand.grid(
    year=seq(1980, 2022), morph=c("gray","brown"), observer="a"
  )))
data_pred_modS1a[morph=="gray", estimate := estimate+3]
data_pred_modS1a[morph=="brown", estimate := estimate+9]

p2 =
ggplot(data_summ1) + ggtitle("B") +
  coord_cartesian(ylim = c(4,14)) + 
  geom_hline(yintercept=9, linetype=2) + 
  geom_errorbar(aes(x=year, ymin=color_score_median.MEAN-color_score_median.SE, 
                    ymax=color_score_median.MEAN+color_score_median.SE, group=morph), 
                color="black", width=0) +
  geom_point(aes(x=year, y=color_score_median.MEAN, fill=morph), size=2, pch=21) +
  geom_line(data=data_pred_gam1, aes(x=year, y=estimate, group=morph), size=1) +
  geom_ribbon(data=data_pred_modS1a, aes(x=year, ymin=estimate-std.error, ymax=estimate+std.error, group=morph), alpha=0.3) +
  geom_line(data=data_pred_modS1a, aes(x=year, y=estimate, group=morph), size=1, linetype=2) +
  scale_x_continuous(expand = c(0.02, 0)) + 
  scale_y_continuous(breaks=c(4:14), labels=col_score_labs) +
  scale_fill_manual("Morphs", values=morph_cols, guide="none") +
  scale_color_manual("", values="black") +
  scale_linetype_manual("", values=2, labels=c("immigrants")) +
  guides(color = guide_legend(order=1, title.position="top"), linetype = guide_legend(order=2, title.position="top")) + 
  theme(
    legend.position="none",
    axis.title = element_blank()
  )

data = copy(data_ind_uni)
gamS1b = gam(color_score_median_trans_ord ~ morph  +
              s(year, by=morph, k=5) +
              s(observer, by=morph, bs = 're'),
            method="REML",
            data=data)
summary(gamS1b)
anova(gamS1b)
# gratia::draw(modS7)
capture.output(anova(gamS1b), file = "tables/GAMS1b.txt")

## get estimates without observer RE
data_pred_modS1b <- data.table(predictions(gamS1b, exclude=list(
  "s(observer):morphgray", "s(observer):morphbrown"), newdata = expand.grid(
    year=seq(1980, 2022), morph=c("gray", "brown"), observer="c"
  )))
data_pred_modS1b[morph=="gray", estimate := estimate+3]
data_pred_modS1b[morph=="brown", estimate := estimate+9]

p3 =
  ggplot(data_summ1) + ggtitle("C") +
    coord_cartesian(ylim = c(4,14)) + 
    geom_hline(yintercept=9, linetype=2) + 
    geom_errorbar(aes(x=year, ymin=color_score_median.MEAN-color_score_median.SE, 
                      ymax=color_score_median.MEAN+color_score_median.SE, group=morph), 
                  color="black", width=0) +
    geom_point(aes(x=year, y=color_score_median.MEAN, fill=morph), size=2, pch=21) +
    geom_line(data=data_pred_gam1, aes(x=year, y=estimate, group=morph), size=1) +
    geom_ribbon(data=data_pred_modS1b, aes(x=year, ymin=estimate-std.error, ymax=estimate+std.error,   group=morph), alpha=0.3) +
    geom_line(data=data_pred_modS1b, aes(x=year, y=estimate, group=morph), size=1, linetype=3) +
    scale_x_continuous(expand = c(0.02, 0)) + 
    scale_y_continuous(breaks=c(4:14), labels=col_score_labs) +
    scale_fill_manual("Morphs", values=morph_cols, guide="none") +
    scale_color_manual("", values="black") +
    scale_linetype_manual("", values=2, labels=c("new cohort")) +
    guides(color = guide_legend(order=1, title.position="top"), linetype = guide_legend(order=2, title.position="top")) + 
    theme(
      legend.position="none",
      axis.title = element_blank()
    )

p_grid1 = plot_grid(p1, p2, p3, ncol=3, align="hv")
p_grid2 = plot_grid(p_grid1, leg1, ncol=1, rel_heights = c(0.85,0.15))

ggsave("figures/figureS2.jpg", p_grid2,
       width = 32, height = 22, units = 'cm',bg="white")


# figure S3 ---------------------------------------------------------------

## only complete cases
data_recr_sub = data_recr[complete.cases(color_score_median, colour_score_median_midparent)]

## pedigree file
data_pedigree_sub = data.table(
  animal=data_recr_sub$ring,
  sire=data_recr_sub$ring_m,
  dam=data_recr_sub$ring_f
)

## add founders
all_animals <- data_pedigree_sub$animal
all_parents <- unique(c(data_pedigree_sub$sire, data_pedigree_sub$dam))
all_parents <- all_parents[!is.na(all_parents)]  # remove NAs
missing_parents <- setdiff(all_parents, all_animals)
if (length(missing_parents) > 0) {
  missing_rows <- data.table(
    animal = missing_parents,
    sire = NA_integer_,
    dam = NA_integer_
  )
  data_pedigree_sub <- rbind(data_pedigree_sub, missing_rows, fill=T)
}

data_pedigree_sub[,animal := factor(animal)]

## fixing + proper ordering 
data_pedigree_sub_ext = merge(data_pedigree_sub, unique(data_ind[, .(ring, sex, morph)]), by.x = "animal", by.y="ring", all.x=T)
data_pedigree_sub = data.table(with(data_pedigree_sub_ext, kinship2::fixParents(animal, sire, dam, sex)))
data_pedigree_sub = MasterBayes::orderPed(data_pedigree_sub)
setnames(data_pedigree_sub, c("animal", "dam", "sire", "sex"))
data_pedigree_sub_ext = MasterBayes::orderPed(data_pedigree_sub_ext)
data_pedigree_sub_ext[, sex:= data_pedigree_sub$sex ]


## extended pedigree with sex and morph
data_pedigree_sub_ext$id = NA
data_pedigree_sub_ext$affected = 1
ped_plot <- with(data_pedigree_sub_ext, kinship2::pedigree(
  id = animal,
  dadid = sire,
  momid = dam,
  sex = sex,
  affected = affected,
))

## mean depth of pedigree
depths <- kinship2::kindepth(ped_plot)
summary(depths)

jpeg("figures/figureS3.jpg", width = 5500, height = 1400, res = 150)
plot(ped_plot,
     col = morph_cols[data_pedigree_sub_ext$morph], 
     id = data_pedigree_sub_ext$id,
     mar = c(5,5,7,5),
)
title("Pedigree of tawny owls (1980-2022, complete records only)",
      cex.main=3, adj = 0)
legend(x = 100, y = 8,                         
       legend = c("  Male (gray)", "  Male (brown)", 
                  "  Female (gray)", "  Female (brown)"), 
       pch = c(15, 15, 16, 16),                
       pt.bg = rep(morph_cols, each = 2),      
       col = c("gray", "chocolate4", "gray", "chocolate4"),                     
       pt.cex = 5,        
       cex = 2,
       bty = "n",
       y.intersp = 4)                        
dev.off()
