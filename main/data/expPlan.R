# Packages used directly:
#    png         : for loading png files (earth and mars picrures)
#    scales      : for viridis palette
#    wesanderson : for wes palette

## ---- expPlan --------------------------------------------------------------------
meta <- sample_data(data)

exp.design <- meta %>% data.frame() %>%
  rownames_to_column(var = "id") %>%
  select(id, subject_id, collection_date, days, food) %>%
  mutate_at("collection_date", as.Date) %>%
  mutate_at("subject_id", as.factor) %>%
  arrange(days, subject_id)

plot_exp_design <- function(rect.height = 50, # height of diet bar
                            rect.dist = 20, # distance between bar and points
                            arrow.dist = 3, # distanc between arrows and bars/points
                            arrow.length = 4, # length of arrow ticks
                            arrow.size = pt2ggSize(.4), # size of arrows
                            point.stroke = 2, # Stroke around points
                            point.size = 2, # size of points
                            planet.size = 50, # size of planets
                            planet.dist = 5, # distance between planets and bar
                            text.size = 2, # text size for days
                            text.dist = -0.5, # distance of text from points
                            text.repel = -1.5, # repel for overlapping days
                            subject.pal = pal_npg(), # fill palette for subjects
                            diet.pal = pal_d3(alpha = .8) # fill palette for diets
){
  
  # Sampling points arranged by subject_id to adjust overlap
  exp.points <- exp.design %>%
    mutate(y = as.numeric(subject_id)) %>%
    arrange(desc(subject_id))
  
  # Days of experiment (bottom part of the plot)
  days.text <- exp.design %>%
    group_by(collection_date, days) %>%
    summarise(y = text.dist) %>%
    ungroup()
  move <- c(2,5,10,12) # manual moving to avoid overlap
  days.text[move, "y"] <- days.text[move, "y"] + text.repel
  
  # Planet images
  earth.img <- png::readPNG("data/earth.png")
  mars.img <- png::readPNG("data/mars.png")
  
  earth <- rasterGrob(earth.img)
  mars <- rasterGrob(mars.img)
  
  # Diets (all variants)
  rect <- data.frame(from = c(min(exp.design$days), 250, 270),
                     to = c(250, 270, 520)) %>%
    mutate_all(funs(min(exp.design$collection_date) + .)) %>%
    mutate(ymin = max(as.numeric(exp.design$subject_id)) + rect.dist,
           ymax = ymin + rect.height) %>%
    mutate(diet = factor(c("FV", "SV", "TV"), levels = c("FV", "SV", "TV"))) 
  
  # All part of the experiment (same as diets plus the follow-up)
  all <- data.frame(from = c(min(exp.design$days), 250, 270, 520),
                    to = c(250, 270, 520, max(exp.design$days) + 5)) %>%
    mutate_all(funs(min(exp.design$collection_date) + .))
  
  # Segment between planets
  rect.seg <- all[-2,]
  rect.seg$from[2] <- all$from[2]
  rect.seg <- cbind(rect.seg, rect[,-c(1:2)])
  # Manual adjust to avoid planet skew in the last bar
  # (added the same adjustment that will be removed during plot)
  rect.seg$to[nrow(rect.seg)] <- rect.seg$to[nrow(rect.seg)] + (planet.size/2) + 10
  
  # Journey lengths in days
  to <- all$to # manual adjustment (same as above)
  to[length(to)] <- to[length(to)] + (planet.size/2) + 10
  journey <- data.frame(length = as.numeric(to - all$from)) %>%
    mutate(x = all$from + (length/2),
           y = rect$ymax[1] + (planet.dist/2),
           label = length)
  # Removing the simulated landing
  extr <- which(journey$label == 20)
  journey$label[extr + 1] <- paste0(journey$label[extr], " + ", journey$label[extr + 1])
  journey <- journey[-extr, ]
  # Adjusting again
  journey$label[3] <- all$to[nrow(all)] - all$from[nrow(all)]
  
  
  # y axis breaks
  y.breaks <- c( rect$ymax[1] + (planet.dist/2),
                 (rect.height/2) + rect$ymin[1],
                 max(exp.points$y)/2,
                 (max(days.text$y) + min(days.text$y))/2
  )
  names(y.breaks) <- c("Phases", "Diets", "Samples", "Days")
  
  # Fill colors
  myCol <- c(subject.pal(nlevels(exp.design$subject_id)),
             diet.pal(nlevels(rect$diet)))
  
  # x-axis range
  x.range <- range(exp.design$collection_date)
  x.range[2] <- x.range[2] + 5
  
  # Plot
  ggplot(rect, aes(xmin = from, ymin = ymin, fill = diet)) +
    # dashed lines on top and bottom of diets
    geom_hline(yintercept = c(rect$ymin[1], rect$ymax[1]), linetype=2,
               size = pt2ggSize(.5))+ 
    # Diets
    geom_rect(aes(xmax = to, ymax = ymax), size = pt2ggSize(.4)) +
    # Lines on top and bottom of diets
    geom_segment(aes(x = from, xend = to, y = ymin, yend = ymin),
                 size = pt2ggSize(.5)) +
    geom_segment(aes(x = from, xend = to, y = ymax, yend = ymax),
                 size = pt2ggSize(.7)) +
    # Bars describing journeys and follow-up
    geom_errorbarh(data = rect.seg,
                   aes(xmin = from + (planet.size/2) + 10, 
                       xmax = to - (planet.size/2) - 10, 
                       y = ymax + (planet.dist/2)),
                   size = pt2ggSize(1)) +
    # Time of flight in days
    geom_text(data = journey, aes(x = x, y = y, label = label),
              inherit.aes = F,
              vjust = -0.8, size = text.size) +
    # Planet annotations
    annotation_custom(earth, 
                      xmin = rect$from[1] - (planet.size/2), 
                      xmax = rect$from[1] + (planet.size/2),
                      ymin = rect$ymax[1], 
                      ymax = rect$ymax[1] + planet.dist) +
    annotation_custom(mars, 
                      xmin = rect$from[2] - (planet.size/2), 
                      xmax = rect$from[2] + (planet.size/2),
                      ymin = rect$ymax[2], 
                      ymax = rect$ymax[2] + planet.dist) +
    annotation_custom(earth, 
                      xmin = rect$to[3] - (planet.size/2), 
                      xmax = rect$to[3] + (planet.size/2),
                      ymin = rect$ymax[3], 
                      ymax = rect$ymax[3] + planet.dist) +
    # Sampling points
    geom_point(data = exp.points, 
               aes(collection_date, y = y, fill = subject_id),
               shape = 21,
               size = point.size,
               stroke = point.stroke,
               inherit.aes = F,
               show.legend = F) +
    # Days of mission
    geom_text(data = days.text,
              aes(x = collection_date, y = y, label = days + 1),
              inherit.aes = F,
              vjust = 1,
              size = text.size) +
    # Arrows
    geom_segment(data = exp.design %>%
                   mutate(y=unique(rect$ymin) - arrow.dist, 
                          yend=max(as.numeric(exp.design$subject_id)) + arrow.dist),
                 aes(x=collection_date, 
                     xend=collection_date, 
                     y=y, 
                     yend=yend + .5),
                 color = "black",
                 size = arrow.size,
                 arrow = arrow(length = unit(arrow.length, "pt")),
                 inherit.aes = F) +
    # Axis specification
    scale_y_continuous(limits = c(-4, max(rect$ymax) + planet.dist),
                       expand = c(0,0), breaks = y.breaks, 
                       labels = names(y.breaks)) +
    scale_fill_manual(values = myCol) +
    scale_x_date(breaks = c(rect$from, max(rect$to), 
                            max(x.range)),  
                 position = "top",
                 limits = x.range,
                 labels = function(x){
                   last <- x[length(x) - 1]
                   last <- as.Date(as.Date(last) - 1)
                   x[length(x) - 1] <- last
                   x
                 }) +
    # Plot customiization
    theme_classic(base_size = 8, base_family = "Helvetica") +
    theme_publication(grid = F) +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_text(face = "bold"),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 0))
}

wes <- function(name){
  function(n){
    wesanderson::wes_palette(name, n, type = "discrete")
  }
}

figure.1a <- plot_exp_design(rect.height = 4, rect.dist = 7, arrow.dist = .5, 
                             planet.size = 40, point.size = 3.5, point.stroke = pt2ggSize(1),
                             diet.pal = scales::viridis_pal(), planet.dist = 5, text.size = 2.5,
                             text.repel = -2) +
  theme(plot.margin=unit(c(.1,1,.1,.1),"cm"),
        legend.title = element_blank(),
        axis.title = element_text(size = -10),
        legend.position = "none")
