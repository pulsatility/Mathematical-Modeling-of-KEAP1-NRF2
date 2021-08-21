### {Model 1: One-step ETGE binding, equilibrium KEAP1-NRF2 model for Class I-V NRF2 activators}
#Unit: time = second, concentration = nM

# Load Package
library(deSolve)

## ---------------------------------
## Parameter values
## --------------------------------- 
parameters <- c(
  ClassI_V 	= 0,
  k0		= 0.1721,
  k1		= 0.0141,
  kprime1		= 0.0141,
  k2		= 0.282,
  kprime2		= 0.282,
  k3		= 0.197543,
  kprime3		= 0.197543,
  k4		= 0.196,
  kprime4		= 0.196,
  k5		= log(2)/(40*60),
  k6		= 0.00203,
  kprime6		= 1.1783e-4,
  k7		= 0.01,
  k8		= 0.1,
  k9		= log(2)/(40*60),
  kprime9		= log(2)/(40*60)
)

## ----------------------------------
## Initial conditions
## ----------------------------------
y0 <- c(
  NRF2free 			= 0,
  KEAP1free			= 530,
  KEAP1o_free		= 0,
  KEAP1_NRF2open 		= 0,
  KEAP1o_NRF2open 	= 0,
  KEAP1_NRF2closed 	= 0,
  KEAP1o_NRF2closed 	= 0
)

## ----------------------------------
## Differential equations to solve
## ----------------------------------
Model_1 <- function(times, y, parms)
{
  with(as.list(c(y, parms)),
       {
         #ODEs for the state variables
dNRF2freedt		= k0 - k5 * NRF2free - 2 * k1 * KEAP1free * NRF2free + k2 * KEAP1_NRF2open - 2 * kprime1 * KEAP1o_free* NRF2free + kprime2 * KEAP1o_NRF2open
dKEAP1freedt		= - 2 * k1 * KEAP1free * NRF2free + k2 * KEAP1_NRF2open -  k7 * KEAP1free * ClassI_V + k8 * KEAP1o_free + k6 * KEAP1_NRF2closed + k9 * KEAP1_NRF2open
dKEAP1o_freedt		= - 2 * kprime1 * KEAP1o_free * NRF2free + kprime2 * KEAP1o_NRF2open + k7 * KEAP1free * ClassI_V - k8 * KEAP1o_free + kprime9 * KEAP1o_NRF2open + kprime6 * KEAP1o_NRF2closed
dKEAP1_NRF2opendt	=  2 * k1 * KEAP1free * NRF2free - k2 * KEAP1_NRF2open - k3 * KEAP1_NRF2open + k4 * KEAP1_NRF2closed - k7 * ClassI_V * KEAP1_NRF2open + k8 * KEAP1o_NRF2open - k9 * KEAP1_NRF2open
dKEAP1o_NRF2opendt	= 2 * kprime1* KEAP1o_free * NRF2free - kprime2* KEAP1o_NRF2open - kprime3 * KEAP1o_NRF2open + kprime4 * KEAP1o_NRF2closed + k7 * ClassI_V * KEAP1_NRF2open - k8 * KEAP1o_NRF2open - kprime9 * KEAP1o_NRF2open
dKEAP1_NRF2closeddt	=  k3 * KEAP1_NRF2open - k4 * KEAP1_NRF2closed - k6 * KEAP1_NRF2closed - k7 * ClassI_V * KEAP1_NRF2closed + k8 * KEAP1o_NRF2closed
dKEAP1o_NRF2closeddt	= kprime3 * KEAP1o_NRF2open - kprime4 * KEAP1o_NRF2closed + k7 * ClassI_V * KEAP1_NRF2closed - k8 * KEAP1o_NRF2closed - kprime6 * KEAP1o_NRF2closed
         #Return simulated value of Y at each time step
         list(c(dNRF2freedt, dKEAP1freedt, dKEAP1o_freedt, dKEAP1_NRF2opendt, dKEAP1o_NRF2opendt, dKEAP1_NRF2closeddt, dKEAP1o_NRF2closeddt)); #They need to follow the same order as the ODEs
       })
}

## ----------------------------------
## Run simulation
## ----------------------------------
#Time span of simulation
tspan <- seq(0, 108000, by = 10)

#Call the lsoda function to numerically solve the ODEs
output <- lsoda(y = y0, times = tspan, func = Model_1, parms = parameters)

#Optional: show the result of the first 6 time steps
head(output)

## ----------------------------------
## Plotting simulation results
## ----------------------------------
par(mar = c(5, 5, 1, 1), mgp = c(2, 0.5, 0)) # mar(figure margins):bottom, left, top, right; mgp(label position): distance from label, scale, tick marks to plot 
plot(output, xlab = "Time (S)", ylab = "conc (nM)", las=1, cex.lab=1.4, cex.axis=1.3, col="blue")

