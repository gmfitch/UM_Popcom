# Script for Lab 12: Metapopulations
# EEB 485, Fall 2012
# Ed Baskerville
# 
# Version History:
# 21 November 2012 (Ed Baskerville):
# * first version for EEB 485
# * added two-species dynamics
# * added habitat destruction
# 16 October 2012 (Aaron King):
# * original code for Bio 281

metapopInteractive = function() {
  require(manipulate)

### INITIALIZE METAPOPULATION MODEL ###
metapopSetup = function (p01, p02, npatch=50, nstep=50)
{
	# Randomly choose occupancy of initial patches
	# 0 = unoccupied
	# 1 = species 1
	# 2 = species 2
	u = runif(n=npatch)
	occ = numeric(npatch)
	occ[u < p01] = 1
	occ[u >= p01 & u < p02 + p01] = 2
	
	# Size of layout grid on one side
	nside <<- ceiling(sqrt(9*npatch))
	
	# Randomly choose where to actually put the patches
	# and calculate patch locations
	k = sample.int(n=nside*nside,size=npatch)
	xx = expand.grid(
		x=seq.int(from=2,to=3*nside+2,by=3),
		y=seq.int(from=2,to=3*nside+2,by=3)
	)
	pos <<- xx[k,]+runif(n=2*npatch,min=-0.25,max=0.25)
	n1OverTime <<- sum(occ == 1)
	n2OverTime <<- sum(occ == 2)
	
	return(occ)
}

### STEP OCCUPANCY VECTOR ###
metapopStep = function (occ, c1, e1, c2, e2)
{
	# Calculate patch occupancy and
	# discrete-time probability of colonization/extinction
	# based on continuous-time rates of Poisson process.
	# Colonization by 2 is conditional on not being colonized by 1.
	p1 = mean(occ == 1)
	p2 = mean(occ == 2)
	
	# Extinction probabilities
	pext1 = 1 - exp(-e1)
	pext2 = 1 - exp(-e2)
	
	# P(colonization by 1) = 1 - P(zero colonizations by 1)
	pcol1 = 1 - exp(-c1*p1)
	
	# P(colonization by 2) = 1 - P(zero colonizations by 2)
	pcol2 = 1 - exp(-c2*p2)
	
	ifelse(
		occ == 0,
		ifelse(
			runif(n=length(occ)) < pcol1,
			1,
			ifelse(
				runif(n=length(occ)) < pcol2,
				2,
				0
			)
		),
		ifelse(
			occ == 1,
			ifelse(
				runif(n=length(occ)) < pext1,
				0,
				1
			),
			ifelse(
				occ == 2,
				ifelse(
					runif(n=length(occ)) < pext2,
					0,
					ifelse(
						runif(n=length(occ)) < pcol1,
						1,
						2
					)
				),
				3 # Destroyed habitat stays destroyed
			)
		)
	)
}

# Onscreen location of patches
pos = array(NA,dim=c(0,2))

# Fraction of patches occupied over time
# for species 1 and species 2
n1OverTime = numeric(0)
n2OverTime = numeric(0)

# Number of patches on a side of the square
nside = numeric(0)

# Initial occupancy
occ = metapopSetup(p01=0.2, p02=0.0, npatch=50, nstep=50)

# Run simulation
manipulate(
	{
		if(reset)
		{
			occ = metapopSetup(p01=p01, p02=0, npatch=size, nstep=nstep)
		}
		else
		{
			if(run)
			{	c2=0
				for (k in seq_len(nstep))
				{
					occ = metapopStep(occ,c1=c1,e1=e,c2=0,e2=e)
					n1OverTime <<- append(n1OverTime, sum(occ == 1))
					n2OverTime <<- append(n2OverTime, sum(occ == 2))
				}
			}
		}
		
		op = par(mfrow=c(2,1), mar=(c(5, 4, -0.1, 2) + 0.1))
		plot(
			c(0,3*(nside+1)),c(0,3*(nside+1)),type='n',
			ann=F,xaxt='n',yaxt='n',bty='n'
		)
		symbols(
			x=pos[,1],y=pos[,2],
			 circles=rep((0.8+e)^(-2),nrow(pos)),
			 inches=F,
			 bg=ifelse(
				occ == 0,
				'white',
				ifelse(
					occ == 1,
					rgb(0.9, 0.6, 0),
					ifelse(
						occ == 2,
						rgb(0.35, 0.7, 0.9),
						'black'
					)
				)
			 ),
			 add=T
		)
		plot(
			n1OverTime,
			ylim=c(0,length(occ)),
			xlab="time",
			ylab="number occupied",type='o', pch=16, col=rgb(0.9, 0.6, 0), lwd=2
		)
		if(n2OverTime[1] > 0)
		{
			points(
				n2OverTime, pch=16, col=rgb(0.35, 0.7, 0.9)
			)
			lines(n2OverTime, lwd=2, col=rgb(0.35, 0.7, 0.9))
			legend(x="topleft", legend=c("species 1", "species 2"),
				fill=c(rgb(0.9, 0.6, 0), rgb(0.35, 0.7, 0.9))
			)
		}
		par(op)
	},
	p01 = slider(0, 1, initial=0.2, step=0.1, label='initial occupancy'),
	#p02 = slider(0, 1, initial=0.0, step=0.1, label='initial occupancy (sp. 2)'),
	c1 = slider(0, 1, initial=0.2, step=0.01, label='colonization rate'),
	#c2 = slider(0, 1, initial=0.2, step=0.01, label='colonization rate (sp. 2)'),
	e = slider(0, 1, initial=0.2, step=0.01, label='extinction rate'),
	size = slider(1, 500, initial=50, step=1, label='number of patches'),
	nstep = slider(1, 100, initial=50, step=1, label='number of steps'),
	reset = button("reset"),
	run = button("run")
	#d = slider(0, 1, initial=0.05, step=0.01, label='fraction to destroy')
	#destroy = button("destroy habitat!")
)

}

dpInteractive = function() {
  require(manipulate, quietly = TRUE)
  require(ggplot2, quietly = TRUE)
  
  dp <- function (m, e, p) {
    dp <- m * p * (1 - p) - e * p
    return(dp)
  }
    
  manipulate(
{
  params <- c(m = m, e = e)
  x <- data.frame(p = seq(0, 1, 0.01))
  x$dp <- sapply(x$p, function(p) dp(m = m, e = e, p = p))
  ints <- data.frame(xx = 0, yy = 0)
  int2 <- 1 - e / m
  if(int2 > 0) {
    ints <- rbind(ints, c(int2, 0))
  }

  ggplot() +
    geom_line(data = x, aes(x = p, y = dp, color = 'red'), lwd = 1.25) +
    xlim(0, 1) +
    ylim(-1, 1) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(data = ints, aes(x = xx, y = yy), size = 10, pch = 1) +
    theme(legend.position="none") +
    labs(x = 'p', y = 'dp/dt')
},
m = slider(0, 1, step = 0.01, label = 'm', init = 0.5),
e = slider(0, 1, step = 0.01, label = 'e', init = 0.5)
  )
}

