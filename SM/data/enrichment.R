enrich <- function(grp, categories){
  if(length(grp) != length(categories)){
    stop("Group length and category length must be equal")
  }
  pop.size <- length(grp)
  
  cat.all <- table(categories) %>%
    dplyr::as_tibble()
  
  dplyr::tibble(grp = grp,
                categories = categories) %>%
    # na.omit() %>%
    dplyr::group_by(grp, categories) %>%
    dplyr::summarise(q = dplyr::n()) %>%
    dplyr::mutate(k = sum(q)) %>%
    dplyr::left_join(cat.all, by = "categories") %>%
    na.omit() %>%
    dplyr::mutate(m = n,
                  n = pop.size - m,
                  fraction.pop = m / pop.size,
                  fraction.obs =  q / k,
                  exp.q = k * fraction.pop,
                  log2.fold.enrichment = log2(fraction.obs / fraction.pop),
                  p.value = ifelse(log2.fold.enrichment > 0,
                                   phyper(q = q - 1, m = m, n = n, k = k, lower.tail = F),
                                   phyper(q = q, m = m, n = n, k = k, lower.tail = T)),
                  adj.p.value = p.adjust(p.value, "BH"))
}
