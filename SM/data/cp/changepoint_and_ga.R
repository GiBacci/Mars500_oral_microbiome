library(tidyverse)
library(changepoint.np)
library(GA)
library(cowplot)
library(lmerTest)
library(MuMIn)
library(ggbeeswarm)
library(DESeq2)
library(adespatial)
library(phyloseq)
library(vegan)

source("../../utils.R")

## @knitr functionsCP

# Splitting data according to changing points.
# Observations correspondi to a cp are repeated so
# that they are icluded in both models.
# This function is meant to be used in combination
# with "geom_smooth" to create segmented regression
splitDF <- function(data, var, split){
  isQuo <- any(grepl("quosure", class(var)))
  if(isQuo){
    quo_var <- var
  }else{
    quo_var <- enquo(var)
  }
  spltted <- data %>% mutate(grp1=findInterval(!!quo_var, split, 
                                               left.open = T),
                             grp2=findInterval(!!quo_var, split, 
                                               left.open = F)) %>%
    mutate(grp = grp1 - grp2)
  
  spltted %>% filter(grp == -1) %>%
    mutate(grp1 = grp2) %>%
    bind_rows(spltted) %>%
    arrange(!!quo_var, grp1)%>%
    select(-grp, -grp2) %>%
    dplyr::rename(grp = grp1)
}

# Fitness for GA algorithm
lmFitness <- function(i, cps, data, formula,  correctByR = F,
                      returnModels = F){
  # get x variable
  iv <- formula[[3]]
  int.var <- data %>% pull(!!iv)
  
  # get y variable
  y <- formula[[2]]
  
  if(all(i == 0)){
    m <- lm(formula, data = data)
    if(returnModels)
      return(m)
    r <- summary(m)$r.squared
    obs <- m$model %>% pull(!!y)
    pred <- predict(m)
    rmse <- 1/sqrt(sum((obs - pred)^2)/length(obs))
    res <- ifelse(correctByR, rmse * r, rmse)
    return(res)
  }else{
    # Get cp selected by GA
    cps <- cps[i == 1]
    
    # Check singularities
    data <- data %>% mutate(split = findInterval(int.var, cps, 
                                                 left.open = T, 
                                                 rightmost.closed = T))
    
    singular <- data %>% group_by(split) %>%
      summarize(singular = length(unique(!!iv)) == 1) %>%
      pull(singular) %>% any()
    
    # If one model is singluar return 0
    if(singular)
      return(0)
    
    # Split data according to cps selected
    # fit a model for each segment
    # and return the reciplocal value of the RMSE
    # corrected by the r squared
    models <- data %>%
      split(.$split) %>%
      map(~lm(formula, data = .))
    
    if(returnModels)
      return(models)
    
    total <- models %>%
      map(function(m){
        r <- summary(m)$r.squared
        obs <- m$model %>% pull(!!y)
        pred <- predict(m)
        rmse <- 1/sqrt(sum((obs - pred)^2)/length(obs))
        ifelse(correctByR, rmse * r, rmse)
      }) %>%
      unlist() %>%
      mean()
    return(total)
  }
}


# Fitness for GA algorithm
# r suared is compuetd using MuMIn::r.squaredGLMM
# the first r was taken (see help for details)
lmerFitness <- function(i, cps, data, formula, correctByR = F,
                        returnModels = F){
  library(lme4)
  library(MuMIn)
  
  # get x variable
  iv <- formula[[3]][[2]]
  int.var <- data %>% pull(!!iv)
  
  # get y variable
  y <- formula[[2]]
  
  # if no changepoint selected return 0
  if(all(i == 0)){
    m <- lmer(formula, data = data, REML = F)
    
    if(returnModels)
      return(m)
    
    r <- r.squaredGLMM(m)[,1]
    obs <- m@frame %>% pull(!!y)
    pred <- predict(m)
    rmse <- 1/sqrt(sum((obs - pred)^2)/length(obs))
    res <- ifelse(correctByR, rmse * r, rmse)
    return(res)
  }else{
    
    # Get cp selected by GA
    cps <- cps[i == 1]
    
    # Check singularities
    data <- data %>% mutate(split = findInterval(int.var, cps, 
                                                 left.open = T, 
                                                 rightmost.closed = T))
    
    singular <- data %>% group_by(split) %>%
      summarize(singular = length(unique(!!iv)) == 1) %>%
      pull(singular) %>% any()
    
    # If one model is singluar return 0
    if(singular)
      return(0)
    
    # Split data according to cps selected
    # fit a model for each segment
    # and return the reciplocal value of the RMSE
    # corrected by the r squared
    models <- data %>%
      split(.$split) %>%
      map(~lmer(formula, data = ., REML = F))
    
    if(returnModels)
      return(models)
    
    if(models %>% map(isSingular) %>% any)
      return(0)
    
    total <- models %>%
      map(function(m){
        r <- r.squaredGLMM(m)[,1]
        obs <- m@frame %>% pull(!!y)
        pred <- predict(m)
        rmse <- 1/sqrt(sum((obs - pred)^2)/length(obs))
        ifelse(correctByR, rmse * r, rmse)
      }) %>%
      unlist() %>%
      mean()
    return(total)
  }
}

plotCPS <- function(data, mapping, cps, model = F, alpha = 1/3){
  # Plotting
  p <- ggplot(data, mapping) +
    geom_vline(xintercept = cps, color = "red",
               linetype = 3) +
    geom_beeswarm(alpha = alpha, priority = "density") +
    # geom_point(alpha = 1/3) +
    stat_summary(geom = "line", fun = mean, show.legend = F) +
    stat_summary(geom = "point", shape = 21, fill = "white",
                 fun = mean, show.legend = F) +
    theme_bw(base_size = 8, base_line_size = pt2ggSize(.5),
             base_family = "Helvetica") +
    theme_publication()
  
  if(model){
    lm.data <- splitDF(data, mapping$x, cps)
    p <- ggplot(data, mapping) + 
      geom_smooth(data = lm.data, 
                         method = "lm", se = F, 
                         color = "red",
                         aes(group = grp)) +
      geom_vline(xintercept = cps, color = "red",
                 linetype = 3) +
      geom_beeswarm(alpha = alpha, priority = "density") +
      # geom_point(alpha = 1/3) +
      stat_summary(geom = "line", fun = mean, show.legend = F) +
      stat_summary(geom = "point", shape = 21, fill = "white",
                   fun = mean, show.legend = F) +
      theme_bw(base_size = 8, base_line_size = pt2ggSize(.5),
               base_family = "Helvetica") +
      theme_publication()
  }
  p + scale_y_continuous(limits = c(0, 1), breaks = c(0, .5, 1))
}

## @knitr runCPanalysis

# Loading data
data <- readRDS("../../../main/data/phylo_obj.rds")

# Lower bound cutoffs.
# pers.cut : persistence cutoff (keep in more than)
# ab.cut : abundance cutoff (keep in more than)
pers.cut <- nsamples(data)*0.05
ab.cut <- 10

data <- filter_taxa(data, function(x){
  pers <- sum(sign(x))/length(x)
  ab <- sum(x)
  !(pers < pers.cut & ab < ab.cut)
}, T)

# Normalizing
X <- data.frame(t(otu_table(data)), check.names = F)
meta <- sample_data(data) %>%
  data.frame()
meta$subject_id <- factor(meta$subject_id)

dds <- DESeqDataSetFromMatrix(X, meta, design = ~ days + food + subject_id)
dds <- estimateSizeFactors(dds, type = "poscounts")

X <- t(counts(dds, norm = T))
X <- sqrt(wisconsin(X))

data.mod <- data
otu_table(data.mod) <- otu_table(X, taxa_are_rows = F)


# Setting up analysis
useMixedModels <- T
rsquareCorrect <- F

# Temporal Beta-diversity Index (Bray-Curtis)
# between multiple time points
days <- unique(sample_data(data)$days)

# Temporal Beta-diversity Index (Bray-Curtis)
# between multiple time points
tbi <- cbind(one=days, two=lead(days)) %>% 
  na.omit() %>%
  apply(., 1, function(d){
    print(paste(d[1], "Vs", d[2]))
    
    # Construct two matrices, one for time point 1
    # and the other for time point 2
    com.1 <- X[sample_data(data)$days == d[1],]
    com.2 <- X[sample_data(data)$days == d[2],]
    
    # Take subjects from each matrix
    s.1 <- sample_data(data)[rownames(com.1),]$subject_id
    s.2 <- sample_data(data)[rownames(com.2),]$subject_id
    
    if(nrow(com.1) != nrow(com.2)){
      # Ordering matrices following subjects.
      # If a matrix is bigger than the other
      # only subjects present in both matrices
      # will be retained.
      ord <- match(s.1, s.2) %>% na.omit()
      s.2 <- s.2[ord]
      com.2 <- com.2[ord,]
      
      extr <- s.1 %in% s.2
      s.1 <- s.1[extr]
      com.1 <- com.1[extr,]
    }
    
    # Computing TBI
    TBI.res <- TBI(com.1, com.2, method = "%difference",
                   save.BC = T, pa.tr = T, test.BC = F,
                   BCD = T)
    
    # Result
    data.frame(subject_id = s.1, one=d[1], two=d[2], 
               tbi=TBI.res$TBI, TBI.res$BC, row.names = NULL,
               TBI.res$BCD.mat, check.names = F)
  }) %>% bind_rows()

# saving table
write.table(tbi, file = "beta_within.csv", sep = "\t", 
            col.names = T, row.names = F,
            quote = F)

# Matrix of TBI (one line per subject)
grp <- paste0("t", dense_rank(tbi$one), "/t", dense_rank(tbi$one) + 1)
grp <- factor(grp, levels = unique(grp))
tbi.mod <- tbi %>% mutate(grp = grp)

mat.cp <- tbi.mod %>%
  select(subject_id, grp, tbi) %>%
  spread(grp, tbi) %>%
  na.omit() %>%
  select(-subject_id) %>%
  as.matrix()

# Detecting changepoint for each subject
cp <- cpt.np(mat.cp, method = "PELT", penalty = "MBIC")
cps <- sapply(cp, cpts) %>% unlist %>% unique %>% sort
cps <- tbi.mod$one[match(colnames(mat.cp)[cps], tbi.mod$grp)]

# Plotting
p1 <- plotCPS(tbi, aes(x = two, y = tbi), cps)

# Optimizing using genetic algorithm
if(!useMixedModels){
  res <- ga(type = "binary",
            fitness = lmFitness,
            cps = cps,
            data = tbi,
            formula = tbi ~ two,
            correctByR = rsquareCorrect,
            returnModels = F,
            popSize = 100,
            nBits = length(cps),
            monitor = plot, parallel = T,
            run = 20, keepBest = F,
            seed = 123)
}else{
  res <- ga(type = "binary",
            fitness = lmerFitness,
            cps = cps,
            data = tbi,
            formula = tbi ~ two + (1|subject_id),
            correctByR = rsquareCorrect,
            returnModels = F,
            popSize = 100,
            nBits = length(cps),
            monitor = plot, parallel = 2,
            run = 20, keepBest = F,
            seed = 123)
}

# Plotting after GA selection
p2 <- plotCPS(tbi, aes(x = two, y = tbi), cps[res@solution == 1])

# Model again removing non significant 
# changepoint (according to lme4)
models <- lmerFitness(res@solution, cps, data = tbi, 
                      formula = tbi ~ two + (1|subject_id),
                      correctByR = rsquareCorrect,
                      returnModels = T)

p.vals <- models %>% map(summary) %>%
  map("coefficients") %>%
  map(10) %>%
  unlist() %>%
  as.vector()

# Removing non significant change points
t.points <- unique(c(min(tbi$one), cps[res@solution == 1], max(tbi$two)))
t.points <- cbind(t.points, lead(t.points)) %>% na.omit() %>%
  data.frame(p.vals) %>% filter(p.vals <= 0.05)
t.points <- unique(as.vector(t.points[[1]]))

final.solution <- rep(0, length(res@solution))
final.solution[which(cps %in% t.points)] <- 1

# Corrected models
final.models <- lmerFitness(final.solution, cps, data = tbi, 
                            formula = tbi ~ two + (1|subject_id),
                            correctByR = F,
                            returnModels = T)
saveRDS(cps[final.solution == 1], "cps_within_crewmembers.rds")
saveRDS(final.models, "final_models_within.rds")

p3 <- plotCPS(tbi, aes(x = two, y = tbi), 
              cps[final.solution == 1],
              model = T)

# Trying with beta-diversity but nothign special
beta <- betaDiv(data.mod ~ subject_id + days, method = "t") %>%
  filter(days_1 == days_2)

b <- beta %>% 
  select(subject_id_1, subject_id_2, days_1, value) %>%
  spread(days_1, value) %>%
  na.omit() %>%
  select(-subject_id_1, -subject_id_2) %>%
  as.matrix()

cp <- cpt.np(b, method = "PELT", penalty = "MBIC")
cps <- sapply(cp, cpts) %>% unlist %>% unique %>% sort
cps <- as.numeric(colnames(b)[cps])  

beta.plot <- beta %>% group_by(days_1) %>%
  filter(!duplicated(value))

p4 <- plotCPS(beta.plot, aes(x = days_1, y = value), cps)

if(!useMixedModels){
  res <- ga(type = "binary",
            fitness = lmFitness,
            cps = cps,
            data = beta,
            formula = value ~ days_1,
            correctByR = rsquareCorrect,
            popSize = 100,
            nBits = length(cps),
            monitor = plot, parallel = T,
            run = 20, keepBest = F,
            seed = 123)
}else{
  res <- ga(type = "binary",
            fitness = lmerFitness,
            cps = cps,
            data = beta,
            formula = value ~ days_1 + (1|subject_id_1),
            correctByR = rsquareCorrect,
            popSize = 100,
            nBits = length(cps),
            monitor = plot, parallel = 2,
            run = 20, keepBest = F,
            seed = 123)
}

p5 <- plotCPS(beta.plot, aes(x = days_1, y = value), 
              cps[res@solution == 1])

# Model again removing non significant 
# changepoint (according to lme4)
models <- lmerFitness(res@solution, cps, 
                      data = beta, 
                      formula = value ~ days_1 + (1|subject_id_1),
                      correctByR = rsquareCorrect,
                      returnModels = T)

p.vals <- models %>% map(summary) %>%
  map("coefficients") %>%
  map(10) %>%
  unlist() %>%
  as.vector()

# Removing non significant change points
t.points <- unique(c(min(beta$days_1), cps[res@solution == 1], max(beta$days_2)))
t.points <- cbind(t.points, lead(t.points)) %>% na.omit() %>%
  data.frame(p.vals) %>% filter(p.vals <= 0.05)
t.points <- unique(as.vector(t.points[[1]]))

final.solution <- rep(0, length(res@solution))
final.solution[which(cps %in% t.points)] <- 1

# The final result produced two non significant models
# removing everithing
final.solution <- rep(0, length(res@solution))

# Corrected models
final.models <- lmerFitness(final.solution, cps, data = beta, 
                            formula = value ~ days_1 + (1|subject_id_1),
                            correctByR = F,
                            returnModels = T)
saveRDS(cps[final.solution == 1], "cps_between_crewmembers.rds")
saveRDS(final.models, "final_models_between.rds")


p6 <- plotCPS(beta.plot, aes(x = days_1, y = value), 
              cps[final.solution == 1],
              model = T)

t <- theme(axis.title = element_blank())
col1 <- plot_grid(p1 + xlab("") + ylab("") + 
                    theme(plot.margin = unit(c(.5,.5,-.5,.5), "lines")), 
                  p2 + xlab("") + ylab("Differences within crewmembers") + 
                    theme(plot.margin = unit(c(-1,.5,-.5,.5), "lines"),
                          plot.background = element_rect(fill = "transparent",
                                                         colour = NA)), 
                  p3 + xlab("") + ylab("") + 
                    theme(plot.margin = unit(c(-1,.5,-.5,.5), "lines"),
                          plot.background = element_rect(fill = "transparent",
                                                          colour = NA)), 
                  ncol = 1, align = "vh")

col2 <- plot_grid(p4 + xlab("") + ylab("") + 
                    theme(plot.margin = unit(c(.5,.5,-.5,.5), "lines")), 
                  p5 + xlab("") + ylab("Differences between crewmembers") + 
                    theme(plot.margin = unit(c(-1,.5,-.5,.5), "lines"),
                          plot.background = element_rect(fill = "transparent",
                                                         colour = NA)), 
                  p6 + xlab("") + ylab("") + 
                    theme(plot.margin = unit(c(-1,.5,-.5,.5), "lines"),
                          plot.background = element_rect(fill = "transparent",
                                                         colour = NA)), 
                  ncol = 1, align = "vh")

all <- plot_grid(col1, col2, ncol = 2, labels = letters[1:2])
p <- ggdraw(add_sub(all, "Days of experiment", 
               vpadding=grid::unit(0,"lines"), 
               y=1, x=0.5,
               size = 8, hjust = 0.5))
saveRDS(p, "changepoint_analysis.rds")