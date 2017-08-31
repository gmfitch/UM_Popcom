lambda_gm <- function(reps = 1000, bSD = 0, pSD = 0) {
  Eseq <- seq(0, 1, 0.001)
  lambda <- rep(NA, length(Eseq))
  
  for(i in 1:length(Eseq)) {    
    pSto <- rnorm(reps, mean = 1, sd = pSD)
    p <- (1 - Eseq[i]) * pSto
    
    bSto <- rnorm(reps, mean = 1, sd = bSD)
    b <- Eseq[i] * bSto
    
    lambda[i] <-  exp(mean(log(p + b)))

  }
  
  dat <- data.frame(E = Eseq, lambda = lambda)
  return(dat)
}

lambdaInteractive = function() {
  require(manipulate, quietly = TRUE)
  require(ggplot2, quietly = TRUE)
  
  manipulate(
  {
    params <- c(bSD = bSD, pSD = pSD)
    out <- lambda_gm(bSD = bSD, pSD = pSD)
    ggplot(out, aes(x = E, y = lambda, color = E)) +
      geom_point() +
      ylab("Geometric Mean of Lambdas") +
      stat_smooth(color = 'black', method = 'loess') +
      ylim(values = c(0.94,1.01)) +
      theme(legend.position = 'none')
  },
  bSD = slider(0, 0.35, label = 'StDev. b', init = 0),
  pSD = slider(0, 0.35, label = 'StDev. p', init = 0)
  )
}

