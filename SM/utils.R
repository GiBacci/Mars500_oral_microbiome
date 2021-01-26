##################################################################################
#######  UTILS SCRIPT ############################################################
##################################################################################

#' Format taxa abundance
#'
#' The function will collapse count of taxa
#' lower than the chosen cutoff in N sample.
#' Taxa will be collapsed into a single taxa
#' called "other".
#'
#' @param x a matrix, a phyloseq object or a data frame
#' @param tax.level the level of the taxa table to be used for collapsing
#'                  counts (it must be a column name of the tax_table)
#' @param rel if TRUE, the abundance will be converted into
#'            relative abundance
#' @param cutoff the abundance cutoff (other)
#' @param pal the color palette to be used
#' @param DF if TRUE only a molteni data frame will be returned
#' @param N the number of samples in which the filter function
#'          must be true
#' @param na.rm if TRUE species with no defined taxonomic
#'              classification at the desired taxonomic level
#'              will be removed BEFORE computing abundance
#' @param ... other params
#'
#' @return a data frame or a plot
#' @export
#'
#' @examples
barTaxa <- function (x, tax.level, rel = F, cutoff = 0.05,
                     pal = pal_npg()(10), DF=F, N = NULL,
                     na.rm = T, other.col = "gray70", ...) {
  UseMethod("barTaxa", x)
}

# taxa : a data frame with taxonomic annotations
# sample : a data frame with sample annotations
barTaxa.matrix <- function (x, tax.level,
                            taxa, sample = NULL, rel = F,
                            cutoff = 0.05, pal = pal_npg()(10),
                            DF=F, N = NULL, na.rm = T,
                            other.col = "gray70") {
  library(phyloseq)
  
  phylo <- phyloseq(otu_table(x, taxa_are_rows = F),
                    tax_table(taxa))
  
  if(!is.null(sample))
    sample_data(phylo) <- sample_data(sample)
  
  barTaxa(phylo, tax.level = tax.level,
          rel = rel, cutoff = cutoff,
          pal = pal, DF=DF, N = N,
          na.rm = na.rm, other.col = other.col)
}

# taxa : a data frame with taxonomic annotations
# sample : a data frame with sample annotations
barTaxa.data.frame <- function (x, tax.level,
                                taxa, sample = NULL, rel = F,
                                cutoff = 0.05, pal = pal_npg()(10),
                                DF=F, N = NULL, na.rm = T,
                                other.col = "gray70") {
  
  barTaxa(as.matrix(x),
          taxa = taxa,
          tax.level = tax.level,
          sample = sample,
          rel = rel, cutoff = cutoff,
          pal = pal, DF=DF, N = N,
          na.rm = na.rm)
}

barTaxa.phyloseq <- function (x, tax.level, rel = F, cutoff = 0.05,
                              pal = pal_npg()(10), DF=F, N = NULL,
                              na.rm = T, other.col = "gray70") {
  library(dplyr)
  library(forcats)
  library(phyloseq)
  
  # Formatting label according to cutoff
  other.lab <- paste0("Other (<", cutoff * 100, "%)")
  
  # Removing not classified taxa
  if(na.rm){
    extr <- !is.na(as.vector(tax_table(x)[,tax.level]))
    data.rel <- prune_taxa(extr, x) %>%
      tax_glom(tax.level)
  }else{
    data.rel <- tax_glom(x, tax.level)
  }
  
  # Collapsing taxa that are lower than cutoff
  # fraction of taxa in N samples
  f <- filterfun_sample(bottomr(cutoff))
  N <- ifelse(is.null(N), nsamples(data.rel), N)
  other <- genefilter_sample(data.rel, f, A = N)
  
  # Number of taxa collapsed (useless)
  ntax.pre <- unique(tax_table(data.rel)[,tax.level]) %>% length()
  tax_table(data.rel)[other, tax.level] <- other.lab
  ntax.post <- unique(tax_table(data.rel)[,tax.level]) %>% length()
  lost <- ntax.pre - ntax.post
  
  # Transforming to relative abundance if needed
  if(rel)
    data.rel <- transform_sample_counts(data.rel, function(x) (x/sum(x)) * 100)
  
  # Melting everithing and reordering based on abundace
  # (other level will be the last anyways)
  bars <- psmelt(data.rel) %>%
    mutate_at(tax.level, ~fct_reorder(., Abundance, .desc = T)) %>%
    mutate_at(tax.level, ~fct_relevel(., other.lab, after = nlevels(.)))
  
  # Formatting color palette
  lev <- get(tax.level, bars)
  if(nlevels(lev) - 1 > length(pal)){
    col <- colorRampPalette(pal)(nlevels(lev) - 1)
  }else{
    col <- pal[1:(nlevels(lev)) - 1]
  }
  col <- append(col, other.col)
  names(col) <- levels(lev)
  
  if(DF){
    # If DF is true the data frame is returned here
    return(bars)
  }
  
  # Plotting
  ylab <- ifelse(rel, "Abundance (%)", "Number of reads")
  ggplot(bars, aes(x = Sample, y = Abundance)) +
    geom_col(aes_string(fill = tax.level)) +
    theme_minimal(base_family = "Helvetica", base_size = 8,
                  base_line_size = pt2ggSize(.5)) +
    scale_fill_manual(values = col) +
    scale_y_continuous(expand = c(0,0)) +
    theme_publication(legend.position = "right",
                      legend.direction = "vertical",
                      grid = F) +
    xlab("Sample") +
    ylab(ylab)
}

bottomr <- function(r, na.rm = T){
  # Function to be used in combination with:
  # phyloseq::filterfun_sample()
  # Make filter fun. return taxa lower then the r fraction of the
  # total abundance
  function(x){
    if (na.rm) {
      x = x[!is.na(x)]
    }
    ord = order(x)
    x[ord] = cumsum(x[ord]/sum(x)) <= r
    x == 1
  }
}

pt2ggSize <- function(pt) {
  # The number of points corresponding to size 1 in ggplot
  magic_number <- 2.14
  pt/magic_number
}

theme_publication <- function(lines = .5, rects = .5, text = 8, ticks.len = 2, legend = 10,
                              legend.position = "bottom", legend.direction = "horizontal",
                              grid = T, ticks = T, axis.line = F) {
  tm <- theme(text = element_text(size = text),
              plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              line = element_line(size = pt2ggSize(lines), lineend = "square"),
              rect = element_rect(size = pt2ggSize(rects)),
              
              # panel.background = element_rect(colour = NA),
              # plot.background = element_rect(colour = NA),
              # panel.border = element_rect(colour = NA),
              
              axis.title = element_text(size = rel(1)),
              axis.title.y = element_text(angle=90, vjust = 2, size = rel(1)),
              axis.title.x = element_text(vjust = -0.2, size = rel(1)),
              
              axis.text = element_text(size = rel(1)), 
              
              # axis.line = element_line(colour="black"),
              axis.ticks = if(ticks) element_line() else element_blank(),
              axis.ticks.length = unit(ticks.len, "pt"),
              
              panel.grid.major.x = if(grid == "x" | grid == T) element_line(colour="#f0f0f0") else element_blank(),
              panel.grid.major.y = if(grid == "y" | grid == T) element_line(colour="#f0f0f0") else element_blank(),
              panel.grid.minor = element_blank(),
              
              legend.key = element_rect(colour = NA),
              legend.position = legend.position,
              legend.direction = legend.direction,
              legend.key.size= unit(legend, "pt"),
              legend.margin = margin(0, 0, 0, 0, "cm"),
              
              legend.title = element_text(face="bold", size = rel(.8)),
              legend.text = element_text(size = rel(.8)),
              
              plot.margin=unit(c(3,2,3,3),"mm"),
              
              strip.background=element_blank(),
              strip.text = element_text(size = rel(1))
  )
  
  if(axis.line){
    tm + theme(axis.line = element_line(size = pt2ggSize(lines), lineend = "square"))
  }else{
    tm
  }
}

#' Within and between group distance
#'
#' Spits a distance matrix (dist object) into
#' tow components:
#'
#' 1. within group distance: distance between
#'        observation in the same group
#' 2. between groups distance: distance between
#'        observation in different groups
#'
#' @param dist a dist object
#' @param group a grouping factor (it will be converted into
#'              a character vector)
#'
#' @return a data frame with four columns:
#'           1. group.1: the name of the first group
#'           2. group.2: the name of the second group
#'           3. type: a column reporting within or
#'                    between specification
#'           4. dist: diversity
#' @export
#'
#' @examples
groupDist <- function(dist, group){
  UseMethod("groupDist", dist)
}

groupDist.dist <- function(dist, group){
  group <- as.character(group)
  
  within <- outer(group, group, "==")
  i <- which(lower.tri(within), arr.ind = T)
  i <- t(apply(i, 1, function(x){
    sort(c(group[x[1]], group[x[2]]))
  }))
  colnames(i) <- c("group.1", "group.2")
  
  within <- ifelse(within[lower.tri(within)],
                   "within", "between")
  
  data.frame(i, type=within,
             dist=as.vector(dist),
             stringsAsFactors = F)
}

betaDiv <- function(formula, data, method = "sor", dist.fun = betadiver, ...){
  # This function computes beta diversity by using the
  # betadiver or the vegdist (or any other function with a method argument)
  # function of vegan package.
  #
  # Args:
  # formula: a formula with the left side being the community data matrix
  #          and the right side being the interest variables. If the left
  #          part of the formula is a phyloseq object, then the function
  #          phyloseq::distance is applyed and `dist.fun` will be ignored.
  # data: a data frame with variables of interest. If the object in the left
  #       part of the formula is a phyloseq object data will be automatically
  #       exported from the object.
  # method: one of the method available in vegan::betadiver(help = T)
  # dist.fun: the distance function to be used (one of vegdit or betadiver).
  #           Ignored if the left part of the formula is a phyloseq object.
  #
  # Returns:
  # a data frame to be used in combination with ggridges:
  # ex:
  # library(phyloseq)
  # data(GlobalPatterns)
  #
  # topk <- genefilter_sample(GlobalPatterns, filterfun_sample(topk(10)))
  # X <- t(otu_table(prune_taxa(topk, GlobalPatterns)))
  # data <- sample_data(GlobalPatterns)
  #
  # beta.div <- betaDiv(X ~ SampleType, data = data)
  # ggplot(beta.div, aes(x = value, y = data1, fill = ..x..)) +
  #     geom_density_ridges_gradient(scale = 3, panel_scaling = F,
  #                                  size = pt2ggSize(.5)) +
  #     facet_wrap(~ data2, nrow = 1)
  
  library(vegan)
  library(reshape)
  
  call <- match.call()
  formula <- call$formula
  
  dataMatrix <- formula[[2]]
  X <- eval(dataMatrix, environment(formula), enclos = .GlobalEnv)
  
  if(class(X) == "phyloseq"){
    library(phyloseq)
    sample_names(X) <- 1:nsamples(X)
    betaDiv <- phyloseq::distance(X, method = method, ...)
    data <- data.frame(sample_data(X))
    N <- nsamples(X)
  }else if(is.numeric(X)){
    stopifnot(!missing(data))
    X <- as.matrix(X)
    rownames(X) <- 1:nrow(X)
    betaDiv <- dist.fun(data.frame(X), method = method, ...)
    N <- nsamples(X)
  }
  
  betaDiv <- melt.matrix(as.matrix(betaDiv))
  betaDiv <- betaDiv[!betaDiv[,1] == betaDiv[,2],]
  
  stopifnot(N == nrow(data))
  
  data <- data.frame(get_all_vars(formula[-2], data))
  data1 <- setNames(data[betaDiv[,1],], paste0(names(data), "_1"))
  data2 <- setNames(data[betaDiv[,2],], paste0(names(data), "_2"))
  
  res <- data.frame(value = betaDiv$value, data1, data2)
  return(res)
}

#' This function take a data.frame as input and
#' computes all combinations of variables.
#' A molteni data.frame will be produced with
#' at least 5 columns:
#'
#' x:          the colum with the first values
#' y:          the column with the second values
#' group.x:    the name of the variable in x
#' group.y:    the name of the variable in y
#' group.both: the contrast itself
#'
#' If specified, the function 'fun' will be
#' applyed to all variables selected (or
#' not excluded).
#'
#' args:
#'
#' data:     the data.frame to transform
#' all:      if TRUE all contrast will be reported
#' diag:     if TRUE contrast between the same variable
#'           will be reported
#' fun:      the function that will be applied to each
#'           variable
#' join.var: if variables are defined in another data frame
#'           it can be joined to the results produced by
#'           this function.
#' join.by:  if specified the join.var data frame will be
#'           joined by the column selected.
#' ...:      variables to include/exclude from the contrasts.
#'           Variables excluded will be added as new column
#'           and can be used for additional grouping (see
#'           tidyr::gather function for a similar syntax)
#'
#' returns:
#' A molteni data.frame with al least 5 columns
#'
#' examples:
#'
#' 1 - All pairwise comparisions in mtcars data.frame excluding cyl
#' that will be used for coloring.
#' Variables are transformed into 01 range in order to plot
#' using facet_grid:
#'
#' range01 <- function(x, ...){(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
#' mtcars %>% pairwise(-cyl, fun = range01, diag = T, all = T) %>%
#' ggplot(aes(x = x, y = y)) +
#'  geom_point(aes(color = cyl)) +
#'  facet_grid(group.y ~ group.x)
#'
#'
#' 2 - Pearson correlation coefficient for all variables within each level
#' of cyl (plotting heatmap):
#'
#' mtcars %>% pairwise(-cyl, fun = range01, diag = F, all = F) %>%
#'  group_by(group.x, group.y, cyl) %>%
#'  summarise(p = cor(x, y)) %>%
#' ggplot(aes(x = group.x, y = group.y, fill = p)) +
#'  geom_tile() +
#'  facet_wrap(~ cyl) +
#'  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'
#'
#' 3 - Plotting all against all plot of OTU abundance useful if
#' we had replicates to plot each replicate separately
#'
#' library(phyloseq)
#' data("GlobalPatterns")
#'
#' taxa <- genefilter_sample(GlobalPatterns, filterfun_sample(topk(2)))
#' f <- prune_taxa(taxa, GlobalPatterns)
#' data.frame(t(otu_table(f))) %>%
#'  bind_cols(data.frame(sample_data(f))) %>%
#'  pairwise(-one_of(sample_variables(f))) %>%
#' ggplot(aes(x = x, y = y, color = X.SampleID)) +
#'  geom_point() +
#'  scale_y_continuous(trans = "log1p") +
#'  scale_x_continuous(trans = "log1p")
#'
#' 4 - Plotting OTU abundance from samples coming from the same
#' environment and samples coming from different environments
#'
#' Pairwise contrasts between all variables.
#' In this case variables are samples (see phyloseq
#' package for details) so we can use the "join.var"
#' option to add sample specifications to the resulting
#' table. Variables and their descriptions are then
#' joined using the "X.SampleID" column.
#'
#' library(phyloseq)
#' library(forcats)
#' data("GlobalPatterns")
#'
#' data.all <- data.frame(otu_table(GlobalPatterns)) %>%
#' pairwise(fun = function(x) log2(x + 1),
#'          join.var = sample_data(GlobalPatterns),
#'          join.by = "X.SampleID",
#'          all = T, diag = T)
#'
#' Correlations between samples coming
#' from the same/different environment.
#' Apparently, samples coming from the same
#' environment seem to be more correlated
#' than those coming from different ones.
#'
#' data.all %>%
#'   group_by(group.x, group.y,
#'            SampleType.x, SampleType.y) %>%
#'   summarize(r = cor(x, y)) %>% # Summarizing correlations
#'   ungroup() %>%
#'   mutate(group.x = factor(group.x),     # Adjusting factor order for
#'          group.y = factor(group.y)) %>% # plotting diagonals correctly
#'   mutate(group.x = fct_relevel(group.x, rev(levels(group.y)))) %>%
#'   ggplot(aes(x = group.x, y = group.y, fill = r)) + # Plotting
#'   geom_tile() +
#'   facet_grid(SampleType.y ~ SampleType.x,
#'              scales = "free", space = "free") +
#'   theme_minimal(base_size = 8, base_family = "Helvetica") +
#'   coord_cartesian(expand = F) +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#'         strip.text.y = element_text(angle = 0, hjust = 0,
#'                                     face = "bold"),
#'         strip.text.x = element_text(angle = 90, hjust = 0,
#'                                     face = "bold"),
#'         panel.spacing = unit(.1, "lines")) +
#'   scale_fill_distiller(palette = "Spectral")
#'
#'
#' Repeating correlations excluding diagonals
#' and repeated contrasts. Anova anlysis
#' to test if samples coming from the
#' same environment are more correlated
#' than those coming from different ones.
#'
#' data <- data.frame(otu_table(GlobalPatterns)) %>%
#'  pairwise(fun = function(x) log2(x + 1),
#'            join.var = sample_data(GlobalPatterns),
#'            join.by = "X.SampleID") %>%
#'   mutate(effect = ifelse(SampleType.x == SampleType.y,
#'                          "within", "between"))
#'
#' correlation <- data %>%
#'   group_by(group.x, group.y,
#'            SampleType.x, SampleType.y,
#'            effect) %>%
#'   summarize(r = cor(x, y)) %>%
#'   ungroup()
#'
#' aov(r ~ effect, data = correlation) %>%
#'   summary() # Significant effect!
#'
pairwise <- function(data, ..., all = F, diag = F, fun = NULL,
                     join.var = NULL, join.by = NULL){
  # Required libraries
  library(tidyr)
  library(dplyr)
  library(forcats)
  library(purrr)
  
  # Parsing args
  quo_args <- quos(...)
  if(length(quo_args) == 0){ # If no variables are
    m <- data                # provided keep the whole
    f <- data.frame()        # data frame
  }else{
    m <- data %>% select(!!!quo_args)
    f <- data %>%
      select(setdiff(colnames(data), colnames(m)))
  }
  
  # Make contrasts
  c <- expand.grid(colnames(m),
                   colnames(m))
  colnames(c) <- c("key1", "key2")
  
  # If fun is provided transform
  # all variable selected using fun
  if(!is.null(fun))
    m <- m %>% mutate_all(fun)
  
  # Filter contrasts based on options
  if(all){ # If all keep all contrast
    if(!diag) # If diag = FALSE remove diags
      c <- c %>% filter(key1 != key2)
  }else{ # If all = FALSE remove redundant contrasts
    d <- c %>% filter(key1 == key2)
    c <- data.frame(t(combn(colnames(m), 2)))
    colnames(c) <- c("key1", "key2")
    if(diag){ # Same as above
      c <- suppressWarnings(bind_rows(c, d))
    }
    rm(d)
  }
  
  # Order levels based on occurrence
  ord <- levels(fct_infreq(c$key1))
  
  # Make the final data
  # c <- c %>% mutate_all(as.character) %>%
  #   rowwise() %>%
  #   do( # repeat all variables as x and y
  #     data.frame(x = m[,.$key1],
  #                y = m[,.$key2],
  #                group.x = .$key1,
  #                group.y = .$key2,
  #                group.both = paste(.$key1, .$key2, sep = "_vs_"),
  #                stringsAsFactors = F) %>%
  #       bind_cols(f)
  #   ) %>% ungroup()
  
  c <- c %>% mutate_all(as.character) %>%
    pmap(function(key1, key2){
      res <- data.frame(x = m[,key1],
                        y = m[,key2],
                        group.x = key1,
                        group.y = key2,
                        group.both = paste(key1, key2, sep = "_vs_"),
                        stringsAsFactors = F)
      if(nrow(f) > 0){
        res <- cbind(res,f)
      }
      return(res)
    }) %>% bind_rows() %>% ungroup()
  
  # Joining variables if specified
  if(!is.null(join.var)){
    join.var <- as.data.frame(join.var)
    if(is.null(join.by)){
      library(tibble)
      join.var <- join.var %>%
        rownames_to_column(var = "joining.by")
      join.by <- "joining.by"
    }
    c <- c %>% left_join(join.var, by = c("group.x" = join.by)) %>%
      left_join(join.var, by = c("group.y" = join.by))
  }
  
  ord.x <- ord[ord %in% unique(c$group.x)]
  ord.y <- ord[ord %in% unique(c$group.y)]
  
  c <- c %>% mutate(group.x = fct_relevel(group.x, ord.x),
                    group.y = fct_relevel(group.y, ord.y))
  
  return(c)
}

ggRareCurve <- function(x, ...) UseMethod("ggRareCurve")

ggRareCurve.default <- function(otu_table, step_size = 100){
  require(vegan)
  require(reshape)
  
  curve <- rarecurve(otu_table, step=step_size)
  
  n <- max(sapply(curve, length))
  m <- do.call(cbind, lapply(curve, "[", seq_len(n)))
  rownames(m) <- 1:nrow(m)
  d <- data.frame(m)
  sample_size <- seq(1:nrow(d)) * step_size
  
  d <- cbind(d,sample_size)
  
  names(d) <- c(row.names(otu_table), "sample_size")
  
  melt_d <- melt(d, id="sample_size")
  melt_d <- melt_d[!is.na(melt_d$value),]  
  
  melt_d <- melt_d[!is.na(melt_d$value ),]
  
  res <- list(original = curve,
              ggData = melt_d)
  class(res) <- "ggRareCurve"
  res
}

getDiversityIndexes <- function(otutable, by = "row", metadata = NULL) {
  # Computes and table the principal biodiversity indexes
  #
  # Args:
  #   otutable: the otutable
  #   by: the position of the OTUs
  #   metadata: a data frame containing
  #             sample information (one for each row)
  #
  # Returns:
  #   A table containing the main biodiversity indexes
  
  library(vegan)
  
  if(by == "row")
    otutable <- data.frame(t(otutable))
  
  N <- rowSums(otutable)
  otus <- specnumber(otutable)
  shannon <- diversity(otutable, index = "shannon")
  invsimpson <- diversity(otutable, index = "invsimpson")
  evenness <- shannon/log(otus)
  sing <- rowSums(round(otutable) == 1)
  doub <- rowSums(round(otutable) == 2)
  chao <- estimateR(round(otutable))["S.chao1",]
  goods <- (1 - (sing/N)) * 100
  
  
  # Tabling everything
  table <- data.frame(Number_of_clones=N,
                      Number_of_OTUs=otus,
                      Number_of_singletons=sing,
                      Number_of_doubletons=doub,
                      Chao1_richness=chao,
                      Shannon_diversity=shannon,
                      InvSimpson=invsimpson,
                      Evenness=evenness,
                      Good_coverage_estimator=goods)
  
  if(!is.null(metadata))
    return(data.frame(table, metadata = metadata))
  return(table)
}