## @knitr importDiv
source("data.R")
library(plyr)
library(dplyr)
library(glmmADMB)
library(agricolae)
library(phyloseq)
library(ape)
library(ggsci)

## @knitr rarefactionStat
# Rarefaction curves and species accumulation curves
r <- ggRareCurve(round(t(otu)), step_size = 100)
r <- r$ggData
r <- data.frame(r, meta[match(as.character(r$variable), rownames(meta)),])

rare <- ggplot(r, aes(x=sample_size, y=value, group = variable)) +
  geom_line(alpha = .3) +
  facet_wrap(~ subject_id, ncol = 2, scales="free_x") +
  scale_x_continuous(labels = function(x) paste0(x/1000, "k"), 
                     breaks = c(0, 100000, 200000), limits = c(0, 250000)) +
  scale_y_continuous(limits = c(0, 60), expand = c(0,0), breaks = seq(0,70, 10)) +
  xlab("Sample size") +
  ylab("Number of OTUs") +
  coord_cartesian(ylim = c(20, 60)) +
  theme(plot.margin = unit(c(.2,.2,.2,.2), "lines"),
        strip.placement = "inside",
        strip.text = element_text(margin = margin(0,0,0,0, "lines")))

# compute species accumulation curve with random method
sp.accum <- lapply(split(data.frame(round(t(otu))), meta$subject_id), specaccum, method = "random")
# Building data frame to use with ggplot2
sp.accum <- do.call(rbind, lapply(names(sp.accum), function(n){
  s <- sp.accum[[n]]
  data.frame(sites = s$sites, 
             richness = s$richness, 
             sd = s$sd, sample = n)
}))
sp.accum <- transform(sp.accum, min=richness - sd, max=richness + sd)

accum <- ggplot(sp.accum, aes(x = sites, y = richness, color = sample)) +
  geom_line(alpha = .8, size = pt2ggSize(1)) +
  geom_linerange(aes(ymin=min, ymax=max), show.legend = F, size = pt2ggSize(.8)) +
  scale_y_continuous(limits = c(0, 70), expand = c(0,0), breaks = seq(0,70, 10)) +
  scale_x_continuous(limits=c(.5, 15), expand = c(0,0), breaks = seq(1,15,1)) +
  scale_color_npg() +
  theme_publication(grid = "y", text = 10) +
  # theme(legend.justification = c(0,0),
  #       legend.position = c(0,0),
  #       legend.title = element_blank(),
  #       legend.key.height = unit(.7, "lines"),
  #       plot.margin = unit(c(1,1,.2,.2), "lines")) +
  guides(color = guide_legend(override.aes = list(size=rel(1.5)), ncol = 3)) +
  xlab("Number of samples") +
  ylab("Number of OTUs") 

#############################################
### DIVERSITY ESTIMATION (ALPHA AND BETA) ###
#############################################
## @knitr diversityEstimation
# otu <- MRcounts(exp)
div <- getDiversityIndexes(otu, metadata = meta) # main diversity indices
div$inv.simpson <- diversity(otu, "invsimpson", 2)[rownames(div)] # Adding inverse simpson

# Tests... probably useless
if(F){
  div_subjects <- ddply(div, .(metadata.subject_id), summarize, mean = mean(Good_coverage_estimator), 
                        error = sd(Good_coverage_estimator)/sqrt(nrow(div)))
  goodsCov <- data.frame(Subject = div_subjects$metadata.subject_id,
                         `Good's coverage` = with(div_subjects, sprintf("%.2f +/- %.4f", floor(mean * 100)/100, error)),
                         check.names = F)
}

div <- data.frame(shannon=diversity(round(t(otu))), 
                  richness=specnumber(round(t(otu))),
                  good=div$Good_coverage_estimator,
                  inv.simpson=div$inv.simpson,
                  meta)
div <- transform(div, evenness = shannon/log(richness))

# Testing difference between samples

div.index <- c("richness", "evenness", "inv.simpson", "good")
k <- lapply(div.index, function(i){kruskal(y = as.vector(div[i]), trt = div$subject_id)})
names(k) <- div.index

k.res <- data.frame(t(sapply(k, "[[", "statistics")))
k.res <- setNames(k.res, c("X2", "DF", "P-value"))
k.res$`Q-value` <- p.adjust(k.res$`P-value`, "fdr")

k.groups <- lapply(k, "[[", "groups")

summFun <- function(x) {
  m <- floor(mean(x) * 100)/100
  e <- ceiling(sd(x)/sqrt(length(x)) * 100) / 100
  sprintf("%.2f ± %.3f", m, e)
}

summary <- subset(div, select=c(subject_id, richness, 
                                evenness, inv.simpson, good)) %>%
  group_by(subject_id) %>%
  summarise_all(funs(summFun)) %>%
  as.data.frame

summary[,-1] <- do.call(cbind, lapply(div.index, function(i){
  paste0(summary[,i], " (", 
         k.groups[[i]][as.character(summary$subject_id), "groups"], 
         ")")
}))

# pander::pander(k.res)
pander::pander(summary)

pdf("species_curves.pdf", width = 4, height = 4, useDingbats = F)
grid.arrange(accum, rare, ncol = 1, heights = c(1/3, 2/3))
dev.off()

## @knitr rareAccumPlot
grid.arrange(accum, rare, ncol = 1, heights = c(1/3, 2/3))


## @knitr diversityModels
# ALPHA
# Population level models
div.index <- grep("good", div.index, invert = T, value = T)

models <- lapply(div.index, function(i) {
  frm <- as.formula(paste0(i, "~ days"))
  lm(frm, data = div)
})

ablines <- data.frame(t(sapply(models, coef)), variable = div.index,
                      check.names = F)

# Linear mixed-effects models random intercept (subject)
models <- lapply(div.index, function(i) {
  frm <- as.formula(paste0(i, "~ days + (1|subject_id)"))
  glmmadmb(frm, data = div, family = "gaussian", link = "identity")
})
models <- lapply(models, summary)
models <- t(sapply(models, function(x) x$coefficient[-1,]))
rownames(models) <- div.index

pander::pander(models) # no significant effect

# Building data frame for plotting
div <- subset(div, select = -good)
div <- melt(div, measure.vars = div.index)

# Capitalization function for plotting purposes
capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

# BETA
# testing beta diversities against days
kegg <- readRDS("kegg_paths.rds")

# Chao distance
distances <- by(round(t(otu)), meta$days, vegdist, method = "chao")
distances.p <- by(round(t(kegg)), meta$days, vegdist, method = "chao")

# Sorensen distance (number of species)
sorensen <- by(round(t(otu)), meta$days, vegdist, binary = T)
sorensen.p <- by(round(t(kegg)), meta$days, vegdist, binary = T)

# Unifrac distance
tree <- read.tree("RAxML_bestTree.otu.tree.nwk")
phylo <- phyloseq(otu_table(otu, taxa_are_rows = T), 
                  tree, 
                  sample_data(meta),
                  tax_table(as.matrix(tax)))

uni.distances <- by(meta, meta$days, function(m){
  extr <- rownames(m)
  sub <- prune_samples(extr, phylo)
  UniFrac(sub, weighted = F)
})

# Building data frame
distances <- do.call(rbind, lapply(seq_along(distances), function(i){
  v <- as.vector(distances[[i]])
  u <- as.vector(uni.distances[[i]])
  s <- as.vector(sorensen[[i]])
  d <- names(distances)[[i]]
  data.frame(chao = v, unifrac = u, sorensen = s, days = d)
}))

distances <- transform(distances, days = as.numeric(as.character(days)))
distances <- data.frame(distances,
                        meta[match(distances$days, meta$days), c("food", "experiment")])

# Testing like alpha before
div.index <- c("chao", "unifrac", "sorensen")

# Linear modelling (no glm since subjects are now shaded)
models <- lapply(div.index, function(i) {
  frm <- as.formula(paste0(i, "~ days"))
  lm(frm, data = distances)
})

# Getting regression lines
ablines.beta <- data.frame(t(sapply(models, coef)), variable = div.index,
                      check.names = F)

# Tabling
models <- lapply(models, summary)
models <- t(sapply(models, function(x) x$coefficient[-1,]))
rownames(models) <- div.index

pander::pander(models) # Significant effect FOUND!

# Meltign data frame for plotting
distances <- melt(distances, measure.vars = div.index)

# APHA and BETA plots
alphadivs <- ggplot(div, aes(x=days, y=value)) +
  geom_point(alpha = .2, col="blue", position=position_jitter(width = 1, height = 0)) +
  stat_summary(fun.y = mean, geom="point") +
  stat_summary(fun.data = mean_cl_boot, geom="linerange") +
  geom_abline(data = ablines, 
              aes(slope = days, intercept = `(Intercept)`),
              linetype = 2) + 
  switch_annotation +
  diet_annotation +
  scale_y_continuous(expand = c(0.01,0.05)) +
  scale_x_continuous(expand = c(0.03, 0.1)) +
  facet_wrap( ~ variable, scales = "free", 
              labeller = labeller(variable = capitalize),
              strip.position = "left", ncol =1) +
  theme(strip.placement = "outside",
        strip.text = element_text(size=rel(1)),
        axis.title.y = element_blank(),
        plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
  xlab("Days") +
  ggtitle("Within-sample alpha-diversity")

betadivs <- ggplot(distances, aes(x=days, y=value)) +
  geom_point(alpha = .2, col="blue", position=position_jitter(width = 1, height = 0)) +
  stat_summary(fun.y = mean, geom="point") +
  stat_summary(fun.data = mean_cl_boot, geom="linerange") +
  geom_abline(data = ablines.beta, 
              aes(slope = days, intercept = `(Intercept)`),
              linetype = 2) + 
  switch_annotation +
  diet_annotation +
  scale_y_continuous(expand = c(0.01,0.05)) +
  scale_x_continuous(expand = c(0.03, 0.1)) +
  facet_wrap( ~ variable, scales = "free", 
              labeller = labeller(variable = capitalize),
              strip.position = "left", ncol =1) +
  theme(strip.placement = "outside",
        strip.text = element_text(size=rel(1)),
        axis.title.y = element_blank(),
        plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
  xlab("Days") +
  ggtitle("Between-sample beta-diversity")

g1 <- ggplotGrob(alphadivs)
g2 <- ggplotGrob(betadivs)

g <- cbind(g1, g2, size = "last")
grid.newpage()

## @knitr diversityPlot
grid.draw(g)

pdf("diversity_within_between.pdf", width = 7, height = 4, useDingbats = F)
grid.draw(g)
dev.off()

# Building boxplots (handmade) - USELESS
if(F){
  p.up <- ggplot(distances, aes(x = days, y = dist, group = food)) +
    geom_boxplot(outlier.shape = NA)
  p.up.data <- ggplot_build(p.up)$data
  p.up.data <- p.up.data[[1]]
  p.up.data$label <- c("a", "b", "b")
  p.up.data$days <- p.up.data$x
  p.up.data$food <- levels(distances$food)
  p.up.data$dist <- p.up.data$ymax
  
  p.up <- p.up + geom_text(data = p.up.data, aes(label = label), vjust=-1) +
    scale_x_continuous(expand=c(0.03, 0.1)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,0.05,0.1)) +
    coord_cartesian(ylim = c(0, .1)) +
    ylab("Chao index") +
    myTheme +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  # plotting diversity Vs time
  p <- ggplot(distances, aes(x=days, y=dist)) +
    switch_annotation +
    geom_point(alpha = .2, col="blue", position="jitter") +
    stat_summary(fun.y = mean, geom="point") +
    stat_summary(fun.data = mean_cl_boot, geom="linerange") +
    geom_abline(slope = model$coefficients[2], intercept = model$coefficients[1],
                linetype = 2) +
    diet_annotation +
    myTheme +
    scale_x_continuous(expand=c(0.03, 0.1)) +
    scale_y_continuous(limits = c(0,.25), expand = c(0,0), breaks = seq(0,.25,.05)) +
    xlab("Days") +
    ylab("Chao index")
  
  # Building plot
  g1 <- ggplotGrob(p.up)
  g2 <- ggplotGrob(p)
  
  gg <- rbind(g1, g2, size = "last")
  
  grid.newpage()
  pdf("diversity_time.pdf", width = 5, height = 3, useDingbats = F)
  grid.draw(gg)
  dev.off()
  
  grid.newpage()
  ## @knitr plotDiff
  grid.draw(gg)
}

#################################################################
# EFFECT ALONG TIME AND WITHIN SUBJETCS
#################################################################

## @knitr totalCommunity
# Effect of days, subject and food (all significant)
# Permanova adonis2 funciton and post-hoc test
perm <- how(nperm = 1000, blocks = with(meta, subject_id))
model <- adonis2(t(round(otu)) ~ days + subject_id + food, data = meta, 
                 permutations = perm)
variance.explained <- data.frame(Factor = rownames(model),
                                 Sum.of.Squares = model$SumOfSqs,
                                 Percentage = (model$SumOfSqs / sum(model$SumOfSqs)) * 100)

pander::pander(model)
pander::pander(variance.explained)

## @knitr skipthis
pairwise.adonis <- function(x,factors, sim.method, p.adjust.m) {
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem])),] ~
                  factors[factors %in% c(as.character(co[1,elem]),as.character(co[2,elem]))] , method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}

posthoc <- rbind(pairwise.adonis(t(round(otu)), meta$food, "bray", "fdr"), # all significant! Each community is different
                 pairwise.adonis(t(round(otu)), meta$subject_id, "bray", "fdr"))

library(ape)
library(ggrepel)
# Ordination analysis probably useless
# Unifrac Distance
# source("https://bioconductor.org/biocLite.R")
# biocLite("BiocGenerics")
# library(BiocGenerics)
# devtools::install_github("joey711/phyloseq")
library(phyloseq)
tree <- read.tree("RAxML_bestTree.otu.tree.nwk")
phylo <- phyloseq(otu_table(otu, taxa_are_rows = T), 
                  tree, 
                  sample_data(meta),
                  tax_table(as.matrix(tax)))

plot(phy_tree(phylo))

sample_data(phylo)$days.fct <- factor(sample_data(phylo)$days)
unifrac.dist <- UniFrac(phylo, weighted = T)

ord <- ordinate(phylo, unifrac.dist, method = "PCoA", distance = "unifrac")
plot_ordination(phylo, ord, color = "days.fct") +
  stat_ellipse(type = "t", level = .95) +
  facet_wrap(~ subject_id)

library(scales)
plot_heatmap(phylo, "PCoA", "wunifrac", "days", "g",  trans = log_trans(2))

meta
otu.pcoa <- pcoa(unifrac.dist, correction = "lingoes")
env.fit <- envfit(ord$vectors, env = meta, permutations = 1000)

sub.phylo <- subset_taxa(phylo, p == "Proteobacteria")
plot_tree(sub.phylo, color = "subject_id", label.tips = "g", 
          size = "abundance", plot.margin = 0.5, ladderize = TRUE,
          sizebase = 10)

otu.pcoa$values$Rel_corr_eig
d <- data.frame(otu.pcoa$vectors, meta)
centroids <- data.frame(env.fit$factors$centroids)
centroids$label <- rownames(centroids)
centroids <- droplevels(centroids[grep("subject", rownames(centroids)),])

chulls <- ddply(d, .(subject_id), function(df) df[chull(df$Axis.1, df$Axis.2), ])

ggplot(d, aes(x = Axis.1, y = Axis.2, fill = subject_id)) +
  geom_polygon(data=chulls, aes(x=Axis.1, y=Axis.2, color=subject_id), alpha = .2, fill= NA) +
  geom_point(data = centroids, fill = "gray80", shape = 21) +
  geom_text_repel(data = centroids, aes(x = Axis.1, y = Axis.2,
                                         label = label), inherit.aes = F) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(shape = 21) +
  myTheme +
  theme(panel.background = element_rect(color = "black", fill = "transparent"))


## @knitr foodEffect
# Turnover along time for each subject
library(codyn)
turn <- melt(data.frame(t(otu), meta), id.vars = names(meta))
turnover <- do.call(rbind, lapply(c("total","appearance", "disappearance"), function(i){
  turnover <- turnover(turn, time.var = "days", 
                       species.var = "variable", 
                       abundance.var = "value", 
                       replicate.var = "subject_id",
                       metric = i)
  melt(turnover, measure.vars = i)
}))

# Families abundance
ab <- getAbundaces(t(otu), tax = tax$f, metadata = meta, melt = F, relative = T, exclude = 0.5)
ab.melt <- melt(ab, id.vars = names(meta))
ab.melt$value <- as.numeric(ab.melt$value)

# Correcting labels
levels(ab.melt$variable)[grep("\\.D\\.", levels(ab.melt$variable))] <- "Bacteria [D]"
levels(ab.melt$variable)[grep("others", levels(ab.melt$variable))] <- "< 0.5%"

# plotting families with a smoothing curve (try change the span)
library(ggsci)
p1 <- ggplot(ab.melt, aes(x = days, y = value, fill = variable)) +
  stat_smooth(geom="area", position = "stack", method = "loess", 
              se = F, alpha = .8, span = .5) +
  switch_annotation +
  facet_grid(. ~ subject_id, scales = "free", space = "free") +
  # scale_fill_brewer(palette = "Set3") +
  scale_fill_d3("category20") +
  scale_y_continuous(expand = c(0,0), labels = function(x) round(x * 100)) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0,1)) +
  myTheme +
  theme(legend.position = "bottom",
        legend.margin = margin(-.5,0,0,0, "lines"),
        legend.key.height = unit(.6, "lines"),
        legend.title = element_blank(),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust=1)) +
  guides(fill = guide_legend(ncol = 3)) +
  xlab("Days") +
  ylab("OTU abundance (%)")

# Plot total turnover, appearance and disappearance for each subject
p2 <- ggplot(turnover, aes(x = days, y = value, color = variable)) +
  switch_annotation +
  diet_annotation +
  geom_line(size = rel(.7), alpha = .7) +
  scale_y_continuous(limits = c(0,.35), expand = c(0.05,0), breaks = seq(0,.4,.2)) +
  geom_point(data = subset(turnover, variable == "total"), shape = 21, 
             aes(y = .31), fill = "white", color = "black", size = rel(.6)) +
  scale_x_continuous(expand=c(0.03, 0)) +
  scale_color_manual(values = c("black", "#4daf4a", "#e41a1c")) +
  scale_fill_grey() +
  facet_grid(. ~ subject_id, scales = "free", space = "free") +
  myTheme +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.margin = margin(0,0,-1,0, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=rel(1)),
                              ncol = 3)) +
  ylab("OTU turnover")

# Food effect
patients_diets <- do.call(rbind, lapply(levels(meta$subject_id), function(f) {
  sub <- valid[,meta$subject_id == f]
  sub <- cumNorm(sub, p = cumNormStat(sub))
  
  food <- pData(sub)$food
  food <- factor(food, labels = gsub(" ", ".", levels(food)))
  settings = zigControl(maxit = 100, verbose = FALSE)
  mod = model.matrix( ~ 0 + food)
  colnames(mod) = levels(food)
  # fitting the ZIG model
  res = fitZig(obj = sub, mod = mod, control = settings)
  # The output of fitZig contains a list of various useful
  # items. hint: names(res). Probably the most useful is the
  # limma 'MLArrayLM' object called fit.
  zigFit = res$fit
  finalMod = res$fit$design
  
  contrast.matrix = makeContrasts(third.variant - first.variant, 
                                  normal - third.variant, 
                                  normal - first.variant, levels = finalMod)
  
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit2 = eBayes(fit2)
  diets <- subset(topTableF(fit2, number = nrow(valid), adjust.method = "fdr"), adj.P.Val <= 0.05)
  
  # Contribution to the mean of each diet
  diets.means <- data.frame(aggregate(t(MRcounts(sub, norm = T)), 
                                      by = list(food), function(x) sum(x)/ncol(otu)), row.names = 1)
  diets.means <- data.frame(t(sweep(diets.means, 2, rowMeans(otu), "/")))
  
  # Setting up graphical parameters
  diets.means$color <- NA
  diets.means$sizes <- log2(rowMeans(otu))
  diets.means[rownames(diets), "color"] <- as.character(tax[rownames(diets), "p"])
  diets.means$label <- rownames(diets.means)
  
  diets.means$subject_id <- f
  
  diets.means
}))

library(ggtern)
library(ggrepel)
p3 <- ggtern(patients_diets, aes(x = first.variant, y = third.variant, z = normal)) +
  geom_point(aes(color = color, size = sizes), alpha = .6) +
  scale_size_continuous(range = c(.1, 3)) +
  theme_minimal(base_size = 10, base_family = "Helvetica") +
  Tlab("TV") + Llab("FV") + Rlab("NR")  +
  scale_color_brewer(palette = "Dark2", na.value = "grey80") +
  facet_wrap(~ subject_id, ncol = 2) +
  theme(legend.position = "bottom",
        legend.margin = margin(.5,0,0,0, "lines"),
        legend.key.height = unit(.6, "lines"),
        legend.title = element_blank(),
        strip.text = element_blank()) +
  guides(color = guide_legend(ncol = 3), size = FALSE)

# Building a single plot
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

g <- rbind(g2, g1, size = "last")
grid.newpage()

pdf("turnover_families.pdf", width = 7, height = 4, useDingbats = F)
grid.draw(g)
dev.off()

grid.newpage()
pdf("turnover_families_diet.pdf", width = 10, height = 6, useDingbats = F)
grid.arrange(g, p3, ncol = 2, widths = c(0.6, 0.4))
dev.off()

grid.newpage()
## @knitr plotFoodEffect
grid.draw(g)

##############################
# Persistence and abundance: #
##############################
## @knitr pers_abund
# Persistence is defined as the percentage of subjects that share
# a given OTU at a given time point
# Abundance is the sum of all reads assigned to a given OTU
# at a given time point
persistence <- data.frame(t(otu)) %>%
  group_by(days = as.factor(meta$days)) %>%
  mutate_all(funs(sign)) %>%
  summarise_all(funs(sum(.)/n())) %>%
  as.data.frame
  
persistence <- melt(persistence)

abundance <- data.frame(t(otu)) %>% 
  group_by(days = as.factor(meta$days)) %>%
  summarise_all(funs(sum)) %>%
  as.data.frame

abundance <- melt(abundance)

# Merging persistence and abundace into a signle data frame
all <- merge(persistence, abundance, by = c("days", "variable"))
all <- subset(all, value.y > 0) # removing undetected OTUs

# One liner model for each time point
models <- by(all, all$days, function(d) lm(log10(value.y) ~ value.x, data = d))
names(models) <- levels(all$days) # Naming models according to time points
coefs <- data.frame(do.call(rbind,lapply(models, coef))) # getting coefficients
coefs$days <- factor(rownames(coefs)) # adding time points to coefficients

# Reporting results of each model into a data frame
models.res <- data.frame(t(sapply(models, function(m) c(confint(m)[2,], 
                                                        summary(m)$coefficients[2,],
                                                        summary(m)$adj.r.squared))))
models.res$Days <- rownames(models.res)
models.res$q.value <- p.adjust(models.res$Pr...t.., method = "fdr")
names(models.res) <- c("2.5%", "97.5%", "Estimate", "St.error", "t", "P", "R", "Days", "Q")
models.res <- models.res[,c(8, 3, 1, 2, 4, 5, 7, 6, 9)]

write.table(models.res, "persistence_abundace_models.csv", 
            row.names = F, col.names = T, sep = "\t",
            quote = F)

pander::pander(models.res)

# Plotting Persistence Vs Abundance for each time point
p1 <- ggplot(all, aes(x = value.x, y = value.y)) +
  geom_point(alpha = .2, col="blue", position="jitter") +
  scale_y_log10(label = function(x) paste0(x/1000, "k"), 
                breaks = c(100, 1000, 10000, 100000)) +
  geom_vline(xintercept = 0.8, linetype = 3, color = "gray20") +
  geom_abline(data = coefs, aes(slope = value.x, intercept = X.Intercept.),
              linetype = 2) +
  scale_x_continuous(labels = function(x) x * 100, breaks = c(.5, 1)) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  facet_wrap(~ days, ncol = 5) +
  theme(strip.background = element_blank()) +
  xlab("Persistence (%)") +
  ylab(expression(Abundance~(log[10])))


# Getting core and accessory
all$core <- "accessory"
all$core[which(all$value.x >= 0.8)] <- "core"

# Naming OTUs according to taxonomic assignments
all$taxa <- tax[as.character(all$variable),]$p

# adding relative abundance (probably useless)
all <- all %>%
  group_by(days, core) %>%
  mutate(rel = value.y / sum(value.y)) %>%
  as.data.frame

# Plotting box plot of core/accessory taxa
p2 <- ggplot(all, aes(x = core, y = value.y, group = interaction(core, taxa))) +
  geom_boxplot(outlier.shape = NA, aes(fill = taxa),lwd = rel(.3), 
               position = position_dodge(width = .9), show.legend = T) +
  facet_wrap(~ days, ncol = 5) +
  geom_vline(xintercept = 1.5, linetype = 3, color = "gray20") +
  scale_y_log10(label = function(x) paste0(x/1000, "k"),
                breaks = c(100, 1000, 10000, 100000)) +
  theme_bw(base_size = 10, base_family = "Helvetica") +
  scale_x_discrete(label = function(x) ifelse(x == "accessory", "A", "C")) +
  scale_fill_brewer(palette = "Set3") +
  theme(strip.background = element_blank(),
        axis.title.y = element_blank(),
        legend.key.height = unit(.6, "lines"),
        legend.key.width = unit(1.2, "lines"),
        legend.title = element_blank(),
        axis.ticks.x = element_blank()) +
  xlab("Accessory/Core")

# Combining plots
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

g <- cbind(g1, g2, size = "last")

grid.newpage()

## @knitr plot_pers_abund
grid.draw(g)

# Permutational Multivariate Analysis of Variance Using Distance Matrices
## @knitr adonis_pers_abund
# Define core and accessory OTUs
c <- unique(as.character(subset(all, core == "core")$variable))
a <- unique(as.character(subset(all, core == "accessory")$variable))

res <- lapply(names(tax)[-1], function(t){
  core <- melt(aggTax(valid[c,], lvl = t, norm = T, out = "matrix"))
  accessory <- melt(aggTax(valid[a,], lvl = t, norm = T, out = "matrix"))
  
  core$type <- "core"
  accessory$type <- "accessory"
  
  aggregated <- dcast(rbind(core, accessory), X2 + type ~ X1, 
                      value.var = "value", fun.aggregate = sum)
  
  meta.aggregated <- data.frame(type = aggregated[[2]], meta[as.character(aggregated$X2),])
  aggregated <- aggregated[,-c(1,2)]
  
  r <- do.call(rbind, lapply(unique(meta.aggregated$days), function(d){
    extr <- which(meta.aggregated$days == d)
    o <- aggregated[extr,]
    m <- meta.aggregated[extr,]
    r <- data.frame(adonis2(o ~ type, data = m, 
                            permutations = 1000))
    r$R <- r$SumOfSqs[1]/sum(r$SumOfSqs)
    r$days <- d
    r[1,]
  }))
  r$lvl <- ""
  r$lvl[1] <- t
  r
})
res <- data.frame(do.call(rbind, res), row.names = NULL)
res <- res[,c(7, 6, 2, 3, 5, 4)]
res$Q <- p.adjust(res$Pr..F., "fdr")

res <- setNames(res, c("Level", "Days", "Sum of Squares", "F-statistic", "R2", "P-value", "Q-value"))
res$Level <- factor(res$Level)
res$Level <- factor(res$Level, labels = c("", "Class", "Family", "Genus", "Order", "Phylum", "Species"))

pander::pander(res)
