library(deSolve)
library(ggplot2)
library(Rmisc)
library(tidyr)
library(dplyr)

initial_step = 1
final_step = 500
initial_size = c(N = 20, P = 20)

times <- seq(from = initial_step, to = final_step, by = 1)

model_function = function(Time, State, Parameters){
		with(as.list(c(State, Parameters)), {
			dN = r*N - a*N*P
			dP = a*c*N*P - d*P
		return(list(c(dN, dP)))
		})
}

parameter_values = c(r = 0.2, a = 0.01, c = 0.1, d = 0.1)

output = ode(initial_size, times, model_function, parameter_values)

#par(mfrow=c(1, 2))
#The next lines of code create abundance over time graphs

output <- as.data.frame(output) %>%
    rename(Time = time, Prey = N, Predator = P)

tmp <- gather(output, Population, Density, Prey, Predator)
  
p1 <- ggplot(tmp, aes(x = Time, y = Density, color = Population, lty = Population)) +
  geom_line() +
  scale_color_manual(values = c('red','black')) +
  scale_linetype_manual(values = c(2:1)) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.key = element_blank()) +
  ggtitle("Time Series")

#matplot(output[ ,1], output[ , 2:3], type = "l", xlab = "Time", ylab = "Density", main = 'Time Series')
#legend("topright", legend = c("Prey", "Predator"), col=1:2, lty=2)

#The next lines of code create a phase space plot
tmp <- output
tmp$PreyEnd <- tmp$Prey[c(2:final_step,NA)]
tmp$PredatorEnd <- tmp$Predator[c(2:final_step,NA)]
tmp <- tmp[1:nrow(tmp)-1,]

p2 <- ggplot(tmp, aes(x = Prey, y = Predator)) +
  geom_segment(aes(xend = PreyEnd, yend = PredatorEnd), arrow = arrow(length = unit(0.03, "npc"))) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    legend.key = element_blank()) +
  ggtitle("Phase Plot")

multiplot(p1, p2, cols = 2)
