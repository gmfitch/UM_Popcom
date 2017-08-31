# Script for Lab 2: Exponential and Logistic Growth Models
# EEB 485, Fall 2012
# Ed Baskerville
#
# Version History:
# 12 September 2012: rewritten version for RStudio's manipulate package (Ed Baskerville)
# 14 September 2011: original version based on RPanel (Doug Jackson)

### DISCRETE-TIME EXPONENTIAL MODEL
### N_{t+1} = R * N_t

# Plot a single run
discreteExponential = function(N0 = 10, lambda = 1.1, tEnd = 50,
  autoScale = FALSE,
  axisNMax = 1000,
  useLogScale = FALSE, filename = NULL)
{
  cat(sprintf(
    "discreteExponential(N0=%d, lambda=%.1f, tEnd=%d, axisNMax=%d, autoScale=%s, useLogScale=%s)\n",
    N0, lambda, tEnd, axisNMax, autoScale, useLogScale))
  
  # Initialize time, population size variables
  tVec = seq(from=0, to=tEnd, by=1)
  NVec = numeric(tEnd + 1)
  NVec[1] = N0
  
  # Iterate exponential growth for each time step
  for(i in 1:tEnd)
  {
    NVec[i + 1] = lambda * NVec[i]
  }
  
  plotResult(discrete=TRUE, tVec=tVec, NVec=NVec, autoScale=autoScale, axisNMax=axisNMax,
             useLogScale=useLogScale)
  saveResult(tVec=tVec, NVec=NVec, filename=filename)
}

# Interactive manipulator
discreteExponentialInteractive = function()
{
  require(manipulate, quietly=TRUE)
  manipulate(
    invisible(discreteExponential(N0 = N0, lambda = lambda, tEnd = tEnd,
      autoScale = autoScale, axisNMax = axisNMax, useLogScale = useLogScale
    )),
    N0 = makeN0Slider(),
    lambda = makeLambdaSlider("lambda"),
    tEnd = makeTEndSlider(),
    axisNMax = makeAxisNMaxSlider(),
    autoScale = makeAutoScaleCheckbox(),
    useLogScale = makeUseLogScaleCheckbox()
  )
}

### DISCRETE-TIME EXPONENTIAL MODEL WITH ENVIRONMENTAL STOCHASTICITY
### N_{t+1} = R_t * N_t
### where R_t is drawn from a zero-truncated Normal(RMean, lambdaSD)

# Single simulator
discreteExponentialStochastic = function(N0=10, lambdaMean=1.1, lambdaSD=0.1,
  tEnd = 50, autoScale = FALSE,
  axisNMax = 1000, useLogScale = FALSE, seed = NULL,
  filename = NULL)
{
  if(is.null(seed))
  {
    # This is a terrible way to generate random seeds
    # but useful for a more important goal: keep the seed < 10000.
    seed = sample.int(n=10000, size=1)
  }
  set.seed(seed)
  
  cat(sprintf(
    "discreteExponentialStochastic(N0=%d, RMean=%.1f, lambdaSD=%.1f, tEnd=%d, axisNMax=%d, useLogScale=%s, autoScale=%s, seed=%d)\n",
    N0, lambdaMean, lambdaSD, tEnd, axisNMax, useLogScale, autoScale, seed))
  
  # Initialize time, population size variables
  tVec = seq(from=0, to=tEnd, by=1)
  NVec = numeric(tEnd + 1)
  NVec[1] = N0
  
  # Iterate exponential growth for each time step
  for(i in 1:tEnd)
  {
  	lambda = 0
  	while(lambda <= 0)
  	{
  	  lambda = rnorm(n = 1, mean = lambdaMean, sd = lambdaSD)
  	}
    NVec[i + 1] = lambda * NVec[i]
  }
  
  plotResult(discrete=TRUE, tVec=tVec, NVec=NVec, autoScale=autoScale, axisNMax=axisNMax,
             useLogScale=useLogScale)
  saveResult(tVec=tVec, NVec=NVec, filename=filename)
}

# Interactive controller
discreteExponentialStochasticInteractive = function()
{
  require(manipulate, quietly=TRUE)
  manipulate(
    invisible(discreteExponentialStochastic(
      N0 = N0,
      lambdaMean = lambdaMean, lambdaSD = lambdaSD,
      tEnd = tEnd,
      axisNMax = axisNMax, 
      autoScale=autoScale,
      useLogScale = useLogScale
    )),
    N0 = makeN0Slider(),
    lambdaMean = makeLambdaSlider("lambda mean"),
    lambdaSD = makelambdaSDSlider(),
    tEnd = makeTEndSlider(),
    axisNMax = makeAxisNMaxSlider(),
    autoScale = makeAutoScaleCheckbox(),
    useLogScale = makeUseLogScaleCheckbox()
  )
}

### CONTINUOUS-TIME EXPONENTIAL MODEL

continuousExponential = function(N0=10, r=0.1,
  tEnd = 50, autoScale = FALSE, axisNMax = 10000, useLogScale = FALSE, filename = NULL)
{
  require(deSolve, quietly=TRUE)
  
  cat(sprintf(
    "continuousExponential(N0=%d, r=%.1f, tEnd=%d, axisNMax=%d, autoScale=%s, useLogScale=%s)\n",
    N0, r, tEnd, axisNMax, autoScale, useLogScale
  ))
  
  dNdtFunc = function(t, states, params) 
  {
    with(as.list(c(states, params)),
    {
      dNdt = r*N
      return(list(c(dNdt)))
    })
  }
  
  # Set up initial conditions and parameters for ODE simulator
  tVec = seq(from=0, to=tEnd, by=0.01)
  states = c(N = N0)
  params = c(r = r)
  
  # Run ODE simulator
  results = as.data.frame(
    ode(y = states, times = tVec, func = dNdtFunc, parms = params)
  )
  
  plotResult(discrete=FALSE, tVec=results$time, NVec=results$N, autoScale=autoScale, axisNMax=axisNMax,
             useLogScale=useLogScale)
  saveResult(tVec=tVec, NVec=NVec, filename=filename)
}

continuousExponentialInteractive = function()
{
  require(manipulate, quietly=TRUE)
  manipulate(
    invisible(continuousExponential(
      N0 = N0,
      r = r, 
      tEnd = tEnd,
      axisNMax = axisNMax, 
      autoScale=autoScale,
      useLogScale = useLogScale
    )),
    N0 = makeN0Slider(),
    r = makeLittleRSlider("r"),
    tEnd = makeTEndSlider(),
    axisNMax = makeAxisNMaxSlider(),
    autoScale = makeAutoScaleCheckbox(),
    useLogScale = makeUseLogScaleCheckbox()
  )
}

### CONTINUOUS-TIME LOGISTIC MODEL

logistic = function(
  N0=10, r0=0.0, K=750, tEnd = 50, axisNMax = 10000, 
  autoScale = FALSE, useLogScale = FALSE, filename=NULL)
{
  require(deSolve, quietly=TRUE)
  
  cat(sprintf(
    "logistic(N0=%d, r0=%.1f, K=%d, tEnd=%d, axisNMax=%d, autoScale=%s, useLogScale=%s)\n",
    N0, r0, K, tEnd, axisNMax, autoScale, useLogScale
  ))
  
  dNdtFunc = function(t, states, params) 
  {
    with(as.list(c(states, params)),
    {
      dNdt = r0 *(1 - N / K) * N
      return(list(c(dNdt)))
    })
  }
  
  # Set up initial conditions and parameters for ODE simulator
  tVec = seq(from=0, to=tEnd, by=0.01)
  states = c(N = N0)
  params = c(r0 = r0)
  
  # Run ODE simulator
  results = as.data.frame(
    ode(y = states, times = tVec, func = dNdtFunc, parms = params)
  )
  
  plotResult(discrete=FALSE, tVec=results$time, NVec=results$N, autoScale=autoScale, axisNMax=axisNMax,
             useLogScale=useLogScale)
  saveResult(tVec=tVec, NVec=NVec, filename=filename)
}

logisticInteractive = function()
{
  require(manipulate, quietly=TRUE)
  manipulate(
    invisible(logistic(
      N0 = N0,
      r0 = r0, K = K,
      tEnd = tEnd,
      axisNMax = axisNMax, 
      autoScale=autoScale,
      useLogScale = useLogScale
    )),
    N0 = makeN0Slider(),
    r0 = makeLittleRSlider("r0"),
    K = makeKSlider(),
    tEnd = makeTEndSlider(),
    axisNMax = makeAxisNMaxSlider(),
    autoScale = makeAutoScaleCheckbox(),
    useLogScale = makeUseLogScaleCheckbox()
  )
}

### CONTINUOUS-TIME LOGISTIC MODEL WITH EXTRA PLOTS

logisticThreePlots = function(
  N0=10, r0=0.0, K=750, tEnd = 50, axisNMax = 10000, 
  autoScale = FALSE, useLogScale = FALSE, filename=NULL)
{
  require(deSolve, quietly=TRUE)
  
  cat(sprintf(
    "logisticThreePlots(N0=%d, r0=%.1f, K=%d, tEnd=%d, axisNMax=%d, autoScale=%s, useLogScale=%s)\n",
    N0, r0, K, tEnd, axisNMax, autoScale, useLogScale
  ))
  
  dNdtFunc = function(t, states, params) 
  {
    with(as.list(c(states, params)),
    {
      dNdt = r0 *(1 - N / K) * N
      return(list(c(dNdt)))
    })
  }
  
  # Set up initial conditions and parameters for ODE simulator
  tVec = seq(from=0, to=tEnd, by=0.01)
  states = c(N = N0)
  params = c(r0 = r0)
  
  # Run ODE simulator
  results = as.data.frame(
    ode(y = states, times = tVec, func = dNdtFunc, parms = params)
  )
  
  tVec = results$time
  NVec = results$N
  
  perCapitaGrowthRate = r0 * (1 - NVec / K)
  totalGrowthRate = r0 * (1 - NVec / K) * NVec
  
  par(mfrow=c(3,1))
  plotResult(discrete=FALSE, tVec=tVec, NVec=NVec, autoScale=autoScale, axisNMax=axisNMax,
             useLogScale=useLogScale)
  plotResult(discrete=FALSE, tVec=tVec, NVec=perCapitaGrowthRate, autoScale=TRUE, axisNMax=axisNMax,
             useLogScale=useLogScale, yLabel="per-capita growth rate")
  plotResult(discrete=FALSE, tVec=tVec, NVec=totalGrowthRate, autoScale=TRUE, axisNMax=axisNMax,
             useLogScale=useLogScale, yLabel="total growth rate")
  saveResult(tVec=tVec, NVec=NVec, filename=filename)
  par(mfrow=c(1,1))
}

logisticThreePlotsInteractive = function()
{
  require(manipulate, quietly=TRUE)
  manipulate(
    invisible(logisticThreePlots(
      N0 = N0,
      r0 = r0, K = K,
      tEnd = tEnd,
      axisNMax = axisNMax, 
      autoScale=autoScale,
      useLogScale = useLogScale
    )),
    N0 = makeN0Slider(),
    r0 = makeLittleRSlider("r0"),
    K = makeKSlider(),
    tEnd = makeTEndSlider(),
    axisNMax = makeAxisNMaxSlider(),
    autoScale = makeAutoScaleCheckbox(),
    useLogScale = makeUseLogScaleCheckbox()
  )
}

### DISCRETE-TIME RICKER MAP
### N_{t+1} = e^(r0 * (1 - N_t / K)) * N_t

# Plot a single run
ricker = function(N0 = 10, r0 = 0.1, K = 750,
  tEnd = 50, autoScale = FALSE, axisNMax = 10000, useLogScale = FALSE, filename = NULL)
{
  cat(sprintf(
    "ricker(N0=%d, r0=%.1f, K=%d, tEnd=%d, axisNMax=%d, autoScale=%s, useLogScale=%s)\n",
    N0, r0, K, tEnd, axisNMax, autoScale, useLogScale
  ))
  
  # Initialize time, population size variables
  tVec = seq(from=0, to=tEnd, by=1)
  NVec = numeric(tEnd + 1)
  NVec[1] = N0
  
  # Iterate exponential growth with density-dependent growth rate
  # for each time step
  for(i in 1:tEnd)
  {
    lambda = exp(r0 * (1 - NVec[i] / K))
    NVec[i + 1] = lambda * NVec[i]
  }
  
  plotResult(discrete=TRUE, tVec=tVec, NVec=NVec, autoScale=autoScale, axisNMax=axisNMax,
             useLogScale=useLogScale)
  saveResult(tVec=tVec, NVec=NVec, filename=filename)
}

# Interactive manipulator
rickerInteractive = function()
{
  require(manipulate, quietly=TRUE)
  manipulate(
    invisible(ricker(
      N0 = N0, r0=r0, K=K,
      tEnd = tEnd,
      axisNMax = axisNMax, 
      autoScale=autoScale,
      useLogScale = useLogScale
    )),
    N0 = makeN0Slider(),
    r0 = makeLittleRSlider("r0"),
    K = makeKSlider(),
    tEnd = makeTEndSlider(),
    axisNMax = makeAxisNMaxSlider(),
    autoScale = makeAutoScaleCheckbox(),
    useLogScale = makeUseLogScaleCheckbox()
  )
}

### DISCRETE-TIME RICKER MAP WITH ENVIRONMENTAL STOCHASTICITY
### N_{t+1} = [e^(r_t * (1 - N_t / K)) + Normal(0, lambdaSD)] * N_t

# Plot a single run
rickerStochastic = function(N0 = 10, r0 = 0.1, K = 750, lambdaSD = 0.1,
  tEnd = 50, autoScale = FALSE, axisNMax = 10000, useLogScale = FALSE, seed=NULL, filename=NULL)
{
  if(is.null(seed))
  {
    # This is a terrible way to generate random seeds
    # but useful for a more important goal: keep the seed < 10000.
    seed = sample.int(n=10000, size=1)
  }
  set.seed(seed)
  
  cat(sprintf(
  	"rickerStochastic(N0=%d, r0=%.1f, K=%d, lambdaSD=%.1f, tEnd=%d, autoScale=%s, axisNMax=%d, useLogScale=%s, seed=%d)\n",
  	N0, r0, K, lambdaSD, tEnd, autoScale, axisNMax, useLogScale, seed
  ))
  
  # Initialize time, population size variables
  tVec = seq(from=0, to=tEnd, by=1)
  NVec = numeric(tEnd + 1)
  NVec[1] = N0
  
  # Iterate exponential growth with density-dependent growth rate
  # for each time step, adding stochastic term to R
  for(i in 1:tEnd)
  {
    lambda = exp(r0 * (1 - NVec[i] / K))
    if(lambda > 0)
    {
		lambdaDelta = -lambda
		while(lambdaDelta <= -lambda)
		{
			lambdaDelta = rnorm(n=1, mean=0, sd=lambdaSD)
		}
		lambda = lambda + lambdaDelta
    }
    NVec[i + 1] = lambda * NVec[i]
  }
  
  plotResult(discrete=TRUE, tVec=tVec, NVec=NVec, autoScale=autoScale, axisNMax=axisNMax,
             useLogScale=useLogScale)
  saveResult(tVec=tVec, NVec=NVec, filename=filename)
}

# Interactive manipulator
rickerStochasticInteractive = function()
{
  require(manipulate, quietly=TRUE)
  manipulate(
    invisible(rickerStochastic(N0 = N0,
      r0=r0,
      K=K,
      lambdaSD=lambdaSD,
      tEnd = tEnd,
      axisNMax = axisNMax, 
      autoScale=autoScale,
      useLogScale = useLogScale
    )),
    N0 = makeN0Slider(),
    r0 = makeLittleRSlider("r0"),
    K = makeKSlider(),
    lambdaSD = makelambdaSDSlider(),
    tEnd = makeTEndSlider(),
    axisNMax = makeAxisNMaxSlider(),
    autoScale = makeAutoScaleCheckbox(),
    useLogScale = makeUseLogScaleCheckbox()
  )
}

### SHARED PLOTTING FUNCTION

plotResult = function(discrete=FALSE, tVec, NVec, autoScale, axisNMax, useLogScale, yLabel="population size")
{
  if(useLogScale)
  {
    axisNMin = 1
    logPar = "y"
  }
  else
  {
    axisNMin = 0
    logPar = ""
  }
  
  if(autoScale)
  {
    ylim = NULL
  }
  else
  {
    ylim = c(axisNMin, axisNMax)
  }
  
  if(discrete)
  {
  	plotType = 'b'
  }
  else
  {
  	plotType = 'l'
  }
  
  # Plot result
  plot(x = tVec, y = NVec, type = plotType, xlab = "time", ylab = yLabel,
       ylim=ylim, log=logPar)
}

saveResult = function(tVec, NVec, filename)
{
  if(is.null(filename))
    return(invisible())
  
  results = data.frame(t = tVec, N = NVec)
  write.csv(x = results, file = filename, row.names = FALSE)
}

makeN0Slider = function()
{
  slider(
    min = 1, max = 100, initial = 10, step = 1, ticks = FALSE
  )
}

makeLambdaSlider = function(label)
{
  slider(
    min = 0.0, max = 20, initial = 1.1, label = label,
    step = 0.1, ticks = FALSE
  )
}

makeLittleRSlider = function(label)
{
  slider(
    min = -5.0, max = 5.0, initial = 0.1, label = label,
    step = 0.1, ticks = FALSE
  )
}

makelambdaSDSlider = function()
{
  slider(
    min = 0.0, max = 1.0, initial = 0.1, label = "lambda std. dev.",
    step = 0.1, ticks = FALSE
  )
}

makeTEndSlider = function()
{
  slider(
    min = 10, max = 200, initial = 50, step = 1, ticks = FALSE, label = "maximum t"
  )
}

makeAxisNMaxSlider = function()
{
  slider(
    min = 1, max = 10000, initial = 1000, label = "y-axis scale",
    step = 1, ticks = FALSE
  )
}

makeAutoScaleCheckbox = function()
{
  checkbox(label = "auto-scale?")
}

makeUseLogScaleCheckbox = function()
{
  checkbox(label = "plot on log scale?")
}

makeKSlider = function()
{
  slider(
    min = 10, max = 2000, initial = 750, label = "K", step = 1, ticks = FALSE
  )
}
