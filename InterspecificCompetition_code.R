# Script for Lab 3: Lotka-Volterra Competition Model
# EEB 485, Fall 2012
# Ed Baskerville
#
# Version History:
# 19 September 2012 (Ed Baskerville):
# * added interactivity via RStudio manipulate package 
# * added phase plotting
# * added ugly hard-coded HTML table-printing code for R Markdown compatibility
# 28 September 2011 (Doug Jackson):
# * original version for EEB 485 Fall 2011 (was Lab 4)

### STANDARD VERSION OF MODEL ###

# Simulate Lotka-Volterra competition without stochasticity
lvCompetition = function(
	plotStyle="time",
	N1_0=100,
	r1=0.1,
	K1=150,
	alpha12=0.5,
	N2_0=99,
	r2=0.1,
	K2=200,
	alpha21=0.5,
	endTime=100,
	printCode=FALSE
)
{
	require(deSolve, quietly=TRUE)
	
	if(printCode)
	{
		cat(sprintf(
			'lvCompetition(plotStyle=\"%s\", N1_0=%d, r1=%.1f, K1=%d, alpha12=%.1f, N2_0=%d, r2=%.1f, K2=%d, alpha21=%.1f, endTime=%d)\n',
			plotStyle, N1_0, r1, K1, alpha12, N2_0, r2, K2, alpha21, endTime
		))
	}
	
	dNdtFunc = function(t, states, param) 
	{
		with(as.list(c(states, param)), 
		{
			dN1dt = r1 * N1 * ((K1-N1-alpha12*N2) / K1)
			dN2dt = r2 * N2 * ((K2-N2-alpha21*N1) / K2)
			list(c(dN1dt, dN2dt))
		})
	}
	
	# Parameter, initial state, and times for deSolve
	params = c(r1=r1, K1=K1, alpha12=alpha12, r2=r2, K2=K2, alpha21=alpha21) 
	states = c(N1=N1_0, N2=N2_0) 
	times = seq(0, endTime, by = 0.01)

	# Execute ODE
	out = as.data.frame(ode(y = states, times = times, func = dNdtFunc, parms = params))
	
	mainTitle = sprintf(
		'N1_0=%d, r1=%.1f, K1=%d, alpha12=%.1f,\nN2_0=%d, r2=%.1f, K2=%d, alpha21=%.1f',
		N1_0, r1, K1, alpha12, N2_0, r2, K2, alpha21
	)
	
	# Plot results
	if(plotStyle == "time")
	{
		lvTimePlot(
			mainTitle=mainTitle,
			times=times, N1=out$N1, N2=out$N2
		)
	}
	else if(plotStyle == "phase")
	{
		lvPhasePlot(
			mainTitle=mainTitle,
			N1=out$N1, N2=out$N2,
			N1_0=N1_0, r1=r1, K1=K1, alpha12=alpha12,
			N2_0=N2_0, r2=r2, K2=K2, alpha21=alpha21
		)
	}
	
	return(invisible(out))
}
		
# Interactive controller
lvCompetitionInteractive = function(
	N1_0=NULL,
	r1=NULL,
	K1=NULL,
	alpha12=NULL,
	N2_0=NULL,
	r2=NULL,
	K2=NULL,
	alpha21=NULL,
	endTime=NULL
)
{
	require(manipulate, quietly=TRUE)
	
	manipExpr = expression(invisible(
		lvCompetition(
			plotStyle=plotStyle,
			N1_0=N1_0,
			r1=r1,
			K1=K1,
			alpha12=alpha12,
			N2_0=N2_0,
			r2=r2,
			K2=K2,
			alpha21=alpha21,
			endTime=endTime,
			printCode=TRUE
		)))
	
	controls = list()
	controls$`_expr` = manipExpr
	controls$plotStyle = picker(
		"Abundances Over Time" = "time",
		"Phase Space" = "phase",
		label = "Plot Style"
	)
	if(is.null(N1_0))
		controls$N1_0 = slider(min = 1, max = 200, initial = 100, step = 1, ticks = FALSE)
	if(is.null(r1))
		controls$r1 = slider(min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE)
	if(is.null(K1))
		controls$K1 = slider(min = 1, max = 200, initial = 150, step = 1, ticks = FALSE)
	if(is.null(alpha12))
		controls$alpha12 = slider(
			min = 0.0, max = 3.0, initial = 1.2, step = 0.1, ticks = FALSE
		)
	if(is.null(N2_0))
		controls$N2_0 = slider(min = 1, max = 200, initial = 100, step = 1, ticks = FALSE)
	if(is.null(r2))
		controls$r2 = slider(min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE)
	if(is.null(K2))
		controls$K2 = slider(min = 1, max = 200, initial = 150, step = 1, ticks = FALSE)
	if(is.null(alpha21))
		controls$alpha21 = slider(
			min = 0.0, max = 3.0, initial = 0.5, step = 0.1, ticks = FALSE
		)
	if(is.null(endTime))
		controls$endTime = slider(
			min = 0, max = 400, initial = 100, step = 1, ticks = FALSE,
			label = "End Time"
		)
	
	do.call(manipulate, controls)
}

### STOCHASTIC VERSION ###

# Simulate LV competition with (somewhat weird) stochasticity
lvCompetitionStochastic = function(
	plotStyle="time",
	N1_0=100,
	r1Mean=0.1,
	r1StdDev=0.1,
	K1=150,
	alpha12=0.5,
	N2_0=99,
	r2Mean=0.1,
	r2StdDev=0.1,
	K2=200,
	alpha21=0.5,
	endTime=100,
	seed=NULL,
	showPlot=TRUE,
	printCode=FALSE
)
{
	require(deSolve, quietly=TRUE)
	
	if(is.null(seed))
	{
	# This is a terrible way to generate random seeds
	# but useful for a more important goal: keep the seed < 10000.
	seed = sample.int(n=10000, size=1)
	}
	set.seed(seed)
	
	if(printCode)
	{
		cat(sprintf(
			'lvCompetitionStochastic(plotStyle=\"%s\", N1_0=%d, r1Mean=%.1f, r1StdDev=%.1f, K1=%d, alpha12=%.1f, N2_0=%d, r2Mean=%.1f, r2StdDev=%.1f, K2=%d, alpha21=%.1f, endTime=%d, seed=%d)\n',
			plotStyle, N1_0, r1Mean, r1StdDev, K1, alpha12, N2_0, r2Mean, r2StdDev, K2, alpha21, endTime, seed
		))
	}
	
	# Parameter, initial state, and times for deSolve
	params = c(K1=K1, alpha12=alpha12, K2=K2, alpha21=alpha21) 
	states = c(N1=N1_0, N2=N2_0) 
	times = seq(0, endTime, by = 0.01)
	
	r1Vec = rnorm(n = length(times), mean = r1Mean, sd = r1StdDev)
	r2Vec = rnorm(n = length(times), mean = r2Mean, sd = r2StdDev)
	
	dNdtFunc = function(t, states, param) 
	{
		with(as.list(c(states, param)), 
		{
			dN1dt = r1Vec[1 + (t / 0.01)] * N1 * ((K1-N1-alpha12*N2) / K1)
			dN2dt = r2Vec[1 + (t / 0.01)] * N2 * ((K2-N2-alpha21*N1) / K2)
			list(c(dN1dt, dN2dt))
		})
	}

	# Execute ODE
	out = as.data.frame(ode(
		method = "euler",
		y = states, times = times, func = dNdtFunc, parms = params))
	
	mainTitle = sprintf(
		'N1_0=%d, r1Mean=%.1f, r1StdDev=%.1f, K1=%d, alpha12=%.1f,\nN2_0=%d, r2Mean=%.1f, r2StdDev=%.1f, K2=%d, alpha21=%.1f',
		N1_0, r1Mean, r1StdDev, K1, alpha12, N2_0, r2Mean, r2StdDev, K2, alpha21
	)
	
	# Plot results
	if(plotStyle == "time")
	{
		lvTimePlot(
			mainTitle=mainTitle,
			times=times, N1=out$N1, N2=out$N2
		)
	}
	else if(plotStyle == "phase")
	{
		lvPhasePlot(
			mainTitle=mainTitle,
			N1=out$N1, N2=out$N2,
			N1_0=N1_0, r1=r1Mean, K1=K1, alpha12=alpha12,
			N2_0=N2_0, r2=r2Mean, K2=K2, alpha21=alpha21
		)
	}
	
	return(invisible(out))
}

# Interactive controller
lvCompetitionStochasticInteractive = function(
	N1_0=NULL,
	r1Mean=NULL,
	r1StdDev=NULL,
	K1=NULL,
	alpha12=NULL,
	N2_0=NULL,
	r2Mean=NULL,
	r2StdDev=NULL,
	K2=NULL,
	alpha21=NULL,
	endTime=NULL
)
{
	require(manipulate, quietly=TRUE)
	
	
	
	manipExpr = expression(invisible(
		lvCompetitionStochastic(
			plotStyle=plotStyle,
			N1_0=N1_0,
			r1Mean=r1Mean,
			r1StdDev=r1StdDev,
			K1=K1,
			alpha12=alpha12,
			N2_0=N2_0,
			r2Mean=r2Mean,
			r2StdDev=r2StdDev,
			K2=K2,
			alpha21=alpha21,
			endTime=endTime,
			printCode=TRUE
		)))
	
	controls = list()
	controls$`_expr` = manipExpr
	controls$plotStyle = picker(
		"Abundances Over Time" = "time",
		"Phase Space" = "phase",
		label = "Plot Style"
	)
	if(is.null(N1_0))
		controls$N1_0 = slider(min = 1, max = 200, initial = 100, step = 1, ticks = FALSE)
	if(is.null(r1Mean))
		controls$r1Mean = slider(
			min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE
		)
	if(is.null(r1StdDev))
		controls$r1StdDev = slider(
			min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE
		)
	if(is.null(K1))
		controls$K1 = slider(min = 1, max = 200, initial = 150, step = 1, ticks = FALSE)
	if(is.null(alpha12))
		controls$alpha12 = slider(
			min = 0.0, max = 3.0, initial = 1.2, step = 0.1, ticks = FALSE
		)
	if(is.null(N2_0))
		controls$N2_0 = slider(min = 1, max = 200, initial = 100, step = 1, ticks = FALSE)
	if(is.null(r2Mean))
		controls$r2Mean = slider(
			min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE
		)
	if(is.null(r2StdDev))
		controls$r2StdDev = slider(
			min = 0, max = 20, initial = 0.1, step = 0.1, ticks = FALSE
		)
	if(is.null(K2))
		controls$K2 = slider(min = 1, max = 200, initial = 150, step = 1, ticks = FALSE)
	if(is.null(alpha21))
		controls$alpha21 = slider(
			min = 0.0, max = 3.0, initial = 0.5, step = 0.1, ticks = FALSE
		)
	if(is.null(endTime))
		controls$endTime = slider(
			min = 0, max = 400, initial = 100, step = 1, ticks = FALSE,
			label = "End Time"
		)
	
	do.call(manipulate, controls)
}

### PARAMETER SWEEP ###

# Sweeps parameters to explore time-to-exclusion as a12 and a21 change

lvCompetitionSweep = function(
	N1_0 = 10, r1 = 20, K1 = 50, N2_0 = 10, r2 = 20, K2 = 50, endTime = 100,
	start=0.2, end=2.0, step=0.2,
	alpha12Start=start, alpha12End=end, alpha12Step=step,
	alpha21Start=alpha12Start, alpha21End=alpha12End, alpha21Step=alpha12Step
)
{
	require(multicore, quietly=TRUE)
	
	runOne = function(params)
	{
		output = lvCompetition(
			plotStyle = "none",
			N1_0 = N1_0,
			r1 = r1,
			K1 = K1,
			alpha12 = params$alpha12,
			N2_0 = N2_0,
			r2 = r2,
			K2 = K2,
			alpha21 = params$alpha21,
			endTime = endTime
		)
		
		# Did N1 survive?
		excludedTimesN1 = output$time[output$N1 < 1]
		if(length(excludedTimesN1) == 0)
		{
			exclusionTimeN1 = NA
		}
		else
		{
			exclusionTimeN1 = excludedTimesN1[1]
		}
		
		# Did N2 survive?
		excludedTimesN2 = output$time[output$N2 < 1]
		if(length(excludedTimesN2) == 0)
		{
			exclusionTimeN2 = NA
		}
		else
		{
			exclusionTimeN2 = excludedTimesN2[1]
		}
		
		return(list(exclusionTimeN1 = exclusionTimeN1, exclusionTimeN2 = exclusionTimeN2, alpha12 = params$alpha12, alpha21 = params$alpha21))
	}
	
	# Set up sweep for alpha12, alpha21
	alpha12Values = seq(from=alpha12Start, to=alpha12End, by=alpha12Step)
	alpha21Values = seq(from=alpha21Start, to=alpha21End, by=alpha21Step)
	alphaIndexes = expand.grid(
		alpha12Index=seq(length(alpha12Values)),
		alpha21Index=seq(length(alpha21Values))
	)
	
	runParams = apply(X = alphaIndexes,
		MARGIN = 1,
		FUN = function(indexes)
		{
			list(alpha12 = alpha12Values[indexes[1]], alpha21 = alpha21Values[indexes[2]])
		}
	)
	
	# Run parameter sweep in parallel on multiple cores
	results = mclapply(X = runParams, FUN = runOne)
	
	exclusionTimes = matrix(data=0, nrow=length(alpha12Values), ncol=length(alpha21Values))
	
	for(i in seq(dim(alphaIndexes)[1]))
	{
		alpha12Index = alphaIndexes[i,1]
		alpha21Index = alphaIndexes[i,2]
		
		etN1 = results[[i]]$exclusionTimeN1
		etN2 = results[[i]]$exclusionTimeN2
		
		stopifnot(is.na(etN1) || is.na(etN2))
		
		if(!is.na(etN1))
		{
			exclusionTimes[alpha12Index, alpha21Index] = etN1
		}
		else
		{
			exclusionTimes[alpha12Index, alpha21Index] = etN2
		}
	}
	
	rownames(exclusionTimes)=lapply(alpha12Values, function(x) { sprintf('%.1f', x) })
	colnames(exclusionTimes)=lapply(alpha21Values, function(x) { sprintf('%.1f', x) })
	
	return(exclusionTimes)
}

lvCompetitionSweepWindows = function(
	N1_0 = 10, r1 = 20, K1 = 50, N2_0 = 10, r2 = 20, K2 = 50, endTime = 100,
	start=0.2, end=2.0, step=0.2,
	alpha12Start=start, alpha12End=end, alpha12Step=step,
	alpha21Start=alpha12Start, alpha21End=alpha12End, alpha21Step=alpha12Step
)
{
	runOne = function(params)
	{
		output = lvCompetition(
			plotStyle = "none",
			N1_0 = N1_0,
			r1 = r1,
			K1 = K1,
			alpha12 = params$alpha12,
			N2_0 = N2_0,
			r2 = r2,
			K2 = K2,
			alpha21 = params$alpha21,
			endTime = endTime
		)
		
		# Did N1 survive?
		excludedTimesN1 = output$time[output$N1 < 1]
		if(length(excludedTimesN1) == 0)
		{
			exclusionTimeN1 = NA
		}
		else
		{
			exclusionTimeN1 = excludedTimesN1[1]
		}
		
		# Did N2 survive?
		excludedTimesN2 = output$time[output$N2 < 1]
		if(length(excludedTimesN2) == 0)
		{
			exclusionTimeN2 = NA
		}
		else
		{
			exclusionTimeN2 = excludedTimesN2[1]
		}
		
		return(list(exclusionTimeN1 = exclusionTimeN1, exclusionTimeN2 = exclusionTimeN2, alpha12 = params$alpha12, alpha21 = params$alpha21))
	}
	
	# Set up sweep for alpha12, alpha21
	alpha12Values = seq(from=alpha12Start, to=alpha12End, by=alpha12Step)
	alpha21Values = seq(from=alpha21Start, to=alpha21End, by=alpha21Step)
	alphaIndexes = expand.grid(
		alpha12Index=seq(length(alpha12Values)),
		alpha21Index=seq(length(alpha21Values))
	)
	
	runParams = apply(X = alphaIndexes,
		MARGIN = 1,
		FUN = function(indexes)
		{
			list(alpha12 = alpha12Values[indexes[1]], alpha21 = alpha21Values[indexes[2]])
		}
	)
	
	# Run parameter sweep in parallel on multiple cores
	results = lapply(X = runParams, FUN = runOne)
	
	exclusionTimes = matrix(data=0, nrow=length(alpha12Values), ncol=length(alpha21Values))
	
	for(i in seq(dim(alphaIndexes)[1]))
	{
		alpha12Index = alphaIndexes[i,1]
		alpha21Index = alphaIndexes[i,2]
		
		etN1 = results[[i]]$exclusionTimeN1
		etN2 = results[[i]]$exclusionTimeN2
		
		stopifnot(is.na(etN1) || is.na(etN2))
		
		if(!is.na(etN1))
		{
			exclusionTimes[alpha12Index, alpha21Index] = etN1
		}
		else
		{
			exclusionTimes[alpha12Index, alpha21Index] = etN2
		}
	}
	
	rownames(exclusionTimes)=lapply(alpha12Values, function(x) { sprintf('%.1f', x) })
	colnames(exclusionTimes)=lapply(alpha21Values, function(x) { sprintf('%.1f', x) })
	
	return(exclusionTimes)
}

### PLOTTING ###

# Plot time course
lvTimePlot = function(
	mainTitle=NULL,
	times, N1, N2
)
{
	plot(times, N1, type = "l", col="red",
		lwd=2, xlab="", ylab="",
		ylim=c(0, 1.2*max(max(N1), max(N2)))
	)
	lines(times, N2, col="blue", lwd=2,)
	
	legend('topright', c('N1', 'N2'), col=c("red", "blue"), lty=c(1,1))
	title(main=mainTitle, xlab="time", ylab="abundance", cex.main=0.8)
}

# Plot LV vector field
lvPhasePlot = function(
	mainTitle=NULL,
	N1, N2,
	N1_0, r1=0.1, K1=150, alpha12=0.5, N2_0, r2=0.1, K2=200, alpha21=0.5
)
{
	options(warn=-1) # Prevents zero-length arrow warnings from being printed
	
	plot.new()
	plot.window(
	  xlim=c(0, 1.1 * max(N2_0, K2, K1/alpha12)),
	  ylim = c(0, 1.1 * max(N1_0, K1, K2/alpha21))
	)
	
	title(main=mainTitle, xlab=expression(N[2]), ylab=expression(N[1]), cex.main=0.8)
	
	axis(side=1, pos=0)
	axis(side=2, pos=0)
	
	# Plot N1 nullcline: dN1/dt = 0
	# intercept on N1 axis: N1 = K1, N2 = 0
	# intercept on N2 axis: N1 = 0, N2 = K1/alpha12
	lines(x = c(K2, 0), y = c(0, K2/alpha21), col='red')
	mtext(text=expression(K[2]), at = K2, side = 1, line=-2)
	mtext(text=expression(K[2] / alpha[21]), at = K2/alpha21, side = 2, line=-2)
	
	# Plot N2 nullcline: dN2/dt = 0
	# intercept on N2 axis: N2 = K2, N1 = 0
	# intercept on N1 axis: N2 = 0, N1 = K2/alpha21
	lines(x = c(0, K1/alpha12), y = c(K1, 0), col='blue')
	mtext(text=expression(K[1] / alpha[12]), at = K1/alpha12, side = 1, line=-2)
	mtext(text=expression(K[1]), at = K1, side = 2, line=-2)
	
	dN1dtFunc = function(N)
	{
		r1 * N[1] * (1 - (N[1] + alpha12 * N[2])/K1)
	}
	
	dN2dtFunc = function(N)
	{
		r2 * N[2] * (1 - (N[2] + alpha21 * N[1])/K2)
	}
	
	points = expand.grid(
		seq(0, 1.1 * max(K1, K2/alpha21), by=20),
		seq(0, 1.1 * max(K2, K1/alpha21), by=20)
	)
	
	dN1dt = apply(X=points, MARGIN=1, FUN=dN1dtFunc)
	dN2dt = apply(X=points, MARGIN=1, FUN=dN2dtFunc)
	
	magnitude = sqrt(mean(dN1dt*dN1dt + dN2dt*dN2dt)) / 20
	
	startN1 = points[,1] - dN1dt / magnitude
	endN1 = points[,1] + dN1dt / magnitude
	startN2 = points[,2] - dN2dt / magnitude
	endN2 = points[,2] + dN2dt / magnitude
	
	arrows(y0=startN1, x0=startN2, y1=endN1, x1=endN2, length=0.05)
			
	points(y=N1[1], x=N2[1], pch=19)
	lines(y=N1, x=N2, type='l')
	points(y=N1[length(N1)], x=N2[length(N2)], pch=1)
	
	legend('topright', c('N1 nullcline', 'N2 nullcline'),
		col=c("blue", "red"), lty=c(1,1))
}

### TABLE PRINTING ###

# Because spitting out hard-coded HTML formatting is sometimes the easiest way...
lvCompetitionSweepRmd = function(
	N1_0 = 10, r1 = 20, K1 = 50, N2_0 = 10, r2 = 20, K2 = 50, endTime = 100,
	start=0.2, end=2.0, step=0.2,
	alpha12Start=start, alpha12End=end, alpha12Step=step,
	alpha21Start=alpha12Start, alpha21End=alpha12End, alpha21Step=alpha12Step
)
{
	results = lvCompetitionSweep(N1_0=N1_0, r1=r1, K1=K1,
		N2_0=N2_0, r2=r2, K2=K2, endTime=endTime,
		start=start, end=end, step=step,
		alpha12Start=alpha12Start, alpha12End=alpha12End, alpha12Step=alpha12Step,
		alpha21Start=alpha21Start, alpha21End=alpha21End, alpha21Step=alpha21Step
	)
	
	cat('<table>\n')
	
	# Print header
	cat('<tr>\n')
	cat(sprintf('<td></td><td></td><td align="center" colspan="%d">$\\alpha_{21}$</td>\n', dim(results)[2]))
	cat('</tr>\n')
	
	cat('<tr>\n')
	cat('<td></td><td></td>')
	for(j in seq(dim(results)[2]))
	{
		cat(sprintf('<td><b>%s</b></td>\n', colnames(results)[j]))
	}
	cat('</tr>\n')
	
	# Print rows
	for(i in seq(dim(results)[1]))
	{
		cat('<tr>\n')
		
		if(i == 1)
		{
			cat(sprintf('<td rowspan="%d">$\\alpha_{12}$</td>', dim(results)[1]))
		}
		
		cat(sprintf('<td><b>%s</b></td>\n', rownames(results)[i]))
		for(j in seq(dim(results)[2]))
		{
			cat(sprintf('<td>%.2f</td>\n', results[i,j]))
		}
		cat('</tr>\n')
	}
	
	cat('</table>')
}

lvCompetitionSweepRmdWindows = function(
	N1_0 = 10, r1 = 20, K1 = 50, N2_0 = 10, r2 = 20, K2 = 50, endTime = 100,
	start=0.2, end=2.0, step=0.2,
	alpha12Start=start, alpha12End=end, alpha12Step=step,
	alpha21Start=alpha12Start, alpha21End=alpha12End, alpha21Step=alpha12Step
)
{
	results = lvCompetitionSweepWindows(N1_0=N1_0, r1=r1, K1=K1,
		N2_0=N2_0, r2=r2, K2=K2, endTime=endTime,
		start=start, end=end, step=step,
		alpha12Start=alpha12Start, alpha12End=alpha12End, alpha12Step=alpha12Step,
		alpha21Start=alpha21Start, alpha21End=alpha21End, alpha21Step=alpha21Step
	)
	
	cat('<table>\n')
	
	# Print header
	cat('<tr>\n')
	cat(sprintf('<td></td><td></td><td align="center" colspan="%d">$\\alpha_{21}$</td>\n', dim(results)[2]))
	cat('</tr>\n')
	
	cat('<tr>\n')
	cat('<td></td><td></td>')
	for(j in seq(dim(results)[2]))
	{
		cat(sprintf('<td><b>%s</b></td>\n', colnames(results)[j]))
	}
	cat('</tr>\n')
	
	# Print rows
	for(i in seq(dim(results)[1]))
	{
		cat('<tr>\n')
		
		if(i == 1)
		{
			cat(sprintf('<td rowspan="%d">$\\alpha_{12}$</td>', dim(results)[1]))
		}
		
		cat(sprintf('<td><b>%s</b></td>\n', rownames(results)[i]))
		for(j in seq(dim(results)[2]))
		{
			cat(sprintf('<td>%.2f</td>\n', results[i,j]))
		}
		cat('</tr>\n')
	}
	
	cat('</table>')
}

