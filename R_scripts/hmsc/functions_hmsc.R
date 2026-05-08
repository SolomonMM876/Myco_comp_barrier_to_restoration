
# view trace for Gamma parameter
viewTraceGamma <- function(coda.object){
  out.list <- lapply(coda.object$Gamma, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'parameter', 'species'), sep=', ') %>% 
    select(-species)
  
  ggplot(out.df, aes(x=sample, y=value, colour=chain)) + 
    geom_line() + 
    geom_hline(aes(yintercept=0)) + 
    facet_wrap(vars(parameter), scales='free_y')
}

# view trace for Beta parameter
viewTraceBeta <- function(coda.object, file=NULL, height=10, width=10){
  if(is.null(file)) stop('specify file name')
  out.list <- lapply(coda.object$Beta, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'parameter', 'species'), sep=', ')
  
  p <- lapply(split(out.df, out.df$species), function(x){
    sp <- x$species[1]
    ggplot(x, aes(x=sample, y=value, colour=chain)) + 
      geom_line() + 
      geom_hline(aes(yintercept=0)) + 
      labs(title=sp) + 
      facet_wrap(vars(parameter), scales='free_y')
  })
  ggsave(filename=file, plot=gridExtra::marrangeGrob(p, nrow=1, ncol=1), 
    width=width, height=height)
}

# summarise trace for Beta parameter
summaryTraceBeta <- function(coda.object){
  out.list <- lapply(coda.object$Beta, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), 
                 names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'parameter', 'species'), sep=', ')
}


# updated summarise trace for Beta parameter
summaryTraceBetaAlt <- function(coda.object){
  process_chain <- function(x) {
    # Grab the true names directly from the coda object
    true_names <- colnames(x) 
    #force a matrix
    mat_clean <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x))
    # Convert to data frame 
    df <- as.data.frame(mat_clean)
    # FORCE the names back onto the dataframe
    colnames(df) <- true_names
    # Add sample index
    df$sample <- 1:nrow(df)
    return(df)
  }
  out.list <- lapply(coda.object$Beta, process_chain)
  # 3. Bind (Uses rbindlist to avoid vctrs error)
  out.df <- data.table::rbindlist(out.list, idcol = "chain")
  out.df %>% 
    as_tibble() %>% 
    pivot_longer(cols = -c(chain, sample), 
                 names_to = 'parameter', 
                 values_to = 'value') %>% 
    # Turn "B[Param, Sp]" -> "B, Param, Sp]"
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'parameter', 'species'), sep=', ')
}


# view trace for Sigma parameter
viewTraceSigma <- function(coda.object){
  out.list <- lapply(coda.object$Sigma, as.data.frame)
  out.list <- lapply(out.list, function(x){
    x$sample <- 1:nrow(x)
    return(x)})
  out.df <- bind_rows(out.list, .id='chain') %>% 
    pivot_longer(cols=ends_with(']'), names_to='parameter', values_to='value') %>% 
    mutate(parameter = gsub('\\[', ', ', parameter), 
           parameter = gsub('\\]$', '', parameter)) %>% 
    separate(parameter, c('estimate', 'species'), sep=', ')
  
  ggplot(out.df, aes(x=sample, y=log(value), colour=chain)) + 
    geom_line() + 
    facet_wrap(vars(species), scales='free_y')
}

