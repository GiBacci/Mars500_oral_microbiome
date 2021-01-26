# Lading libraries
library(tidyverse)
library(ggsci)
library(cowplot)
library(patchwork)

# Loading script
source("~/Dropbox/Scripts/R/ggBio.R")

# Reading taxonomic tables
silva <- readRDS("../tax_assignments.rds")
homd <- readRDS("tax_assignments_homd.rds")

seqTab <- readRDS("../all_seqtab_nochim.rds")

# Compare lineage function
compareLineage <- function(x, y, clean = T){
  cleanNames <- function(x){
    # Cleaning additional labels
    x <- gsub("[[:punct:]]+", " ", x)
    sapply(strsplit(x, " "), "[", 1)
  }
  if(clean){
    x <- cleanNames(x)
    y <- cleanNames(y)
  }
  if(all(is.na(x))|all(is.na(y)))
    return(NULL)
  highest <- max(which(x == y))
  1:highest
}

# Compare lineages
res <- sapply(1:nrow(silva), function(i, clean){
  index <- compareLineage(homd[i,], silva[i,], clean)
  colnames(silva)[index]
}, clean = T)


# Comparing assignments
compareAssignments <- function(x, y){
  # Total number of ASVs
  N <- length(x)
  # Total number of reads clustered
  N.reads <- sum(seqTab)
  
  # Unknown
  na.x <- is.na(x)
  na.y <- is.na(y)
  
  # Number of reads of known origin
  reads.retained.x <- sum(seqTab[,!na.x])
  reads.retained.y <- sum(seqTab[,!na.y])
  
  # Number of unknown ASVs
  not.assigned.x <- sum(na.x)
  not.assigned.y <- sum(na.y)
  
  # Number of known ASVs
  assigned.x <- sum(!na.x)
  assigned.y <- sum(!na.y)
  
  # Number of ASVs assigne donly in one
  # database and not in the other
  assigned.only.x <- sum(na.y & !na.x)
  assigned.only.y <- sum(na.x & !na.y)
  
  # Results
  c(N=N,
    N.reads=N.reads,
    not.assigned.x=not.assigned.x,
    not.assigned.y=not.assigned.y,
    assigned.only.x=assigned.only.x,
    assigned.only.y=assigned.only.y,
    assigned.x=assigned.x, 
    assigned.y=assigned.y,
    reads.retained.x=reads.retained.x,
    reads.retained.y=reads.retained.y)
}

# Get stats for all taxonomic
# levels
stats <- mapply(function(x, y){
  compareAssignments(x, y)
}, data.frame(homd), data.frame(silva))
# Adding ASVs equally assigned (Concordant)
stats <- rbind(stats, equally.assigned=table(unlist(res))[colnames(stats)])
stats[is.na(stats)] <- 0

# Adding ASVs not equally assigned (Discordant)
not.equally.assigned <- stats["assigned.x",] - (stats["equally.assigned",] + stats["assigned.only.x",])
stats <- rbind(stats, not.equally.assigned=not.equally.assigned)

# Remove species assignements
stats <- subset(stats, select=-species)

# Database information
meta <- gsub("\\.x", "_homd", rownames(stats))
meta <- gsub("\\.y", "_silva", meta)
meta <- lapply(strsplit(meta, "_", fixed = T), "length<-", 2)
meta <- do.call(rbind, meta)
meta[is.na(meta)] <- "all"
colnames(meta) <- c("value", "db")

# Formatting ststs
stats <- data.frame(meta, stats, row.names = NULL)

####### Plotting ########
stats.mod <- stats %>%
  filter(value %in% c("not.assigned", "assigned.only", 
                      "equally.assigned", "not.equally.assigned")) %>%
  gather("tax.lvl", "asv", -value, -db) %>%
  mutate_at("tax.lvl", fct_inorder) %>%
  mutate(perc = (asv/nrow(silva)) * 100)

# a. Number of ASVs exclusively assigned by a database
# b. Number of unknown ASVs
ab <- stats.mod %>%
  filter(db != "all") %>%
  select(-asv) %>%
  spread(db, perc) %>%
  split(.$value) %>%
  imap(function(x, name){
    db.name <- c("HOMD", "Silva SSU")
    lab <- if(name == "assigned.only"){
      paste0("ASVs exclusively\nassigned (", db.name, ", %)")
    }else{
      paste0("Unknown ASVs (", db.name, ", %)")
    }
    
    tag <- if(name == "assigned.only"){
      "a"
    }else{
      "b"
    }
    
    values <- c(x$homd, x$silva)
    limits <- c(floor(min(values)), ceiling(max(values)))
    ggplot(x, aes(x = homd, y = silva, color = tax.lvl)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = 2, size = .25) +
      scale_color_npg(labels = str_to_title) +
      theme_bw(base_size = 8, base_line_size = pt2ggSize(.5),
               base_family = "Helvetica") +
      theme_publication(grid = F) +
      coord_cartesian(xlim = limits, ylim = limits) +
      xlab(lab[1]) +
      ylab(lab[2]) +
      labs(tag = tag) +
      theme(plot.tag = element_text(size = 10, face = "bold"))
  }) %>%
  wrap_plots(ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.title = element_blank(),
        legend.position = "right",
        legend.direction = "vertical")


# Number of equally assigned and not
# equally assigned ASVs
c <- stats.mod %>%
  filter(db == "all") %>%
  mutate(value = case_when(
    value == "equally.assigned" ~ "Concordant",
    T ~ "Discordant"
  )) %>%
  ggplot(aes(x = tax.lvl, y = perc, color = value)) +
  geom_line(aes(group = value)) +
  geom_point() +
  theme_classic(base_size = 8, base_line_size = pt2ggSize(.5),
           base_family = "Helvetica") +
  theme_publication(grid = F) +
  scale_color_aaas() +
  scale_x_discrete(labels = str_to_title) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        plot.tag = element_text(size = 10, face = "bold")) +
  ylab("Number of ASVs (%)") +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  labs(tag = "c")

# Number of reads recovered with both
# databases
d <- stats %>%
  select(value, db, domain) %>%
  filter(value == "reads.retained") %>%
  mutate(perc = (domain/sum(seqTab)) * 100,
         db = case_when(
           db == "homd" ~ "HOMD",
           T ~ "Silva SSU"
         )) %>%
ggplot(aes(x = db, y = perc)) +
  geom_col() +
  theme_classic(base_size = 8, base_line_size = pt2ggSize(.5),
           base_family = "Helvetica") +
  theme_publication(grid = F) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
  ylab("Read recovery (%)") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.tag = element_text(size = 10, face = "bold")) +
  labs(tag = "d")

# Building final plot
second.row <- plot_grid(c, d, rel_widths = c(2,1), 
                        align = "hv", axis = "tb")
p <- plot_grid(ab, second.row, nrow = 2, ncol = 1)
saveRDS(p, "database_comparison.rds")