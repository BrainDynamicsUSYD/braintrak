MCMC master code

Task
- Fit arbitrary parameters and parameter sets
- Adapt display code accordingly


Functions

fit - return a single MCMC fit to a spectrum, takes in priors. Depending on options, hold some parameters constant and/or keep priors constant
parallel_chain - wrapper for the chain function - no knownledge of the system
chain - calculate the markov chain for the problem, depending on likelihood. no knowledge of the system required
objective - calculate the likelihood. full knowledge of the actual system required

fit_cluster, fit_several, fit_single - these are just wrappers
	- requires initializing the fit, thus knowledge of the system

params_from_p - map accordingly

Requirements
- Easily specify new parameters to fit and link


Model
- Provides parameters, objective function, constraints
	- Also dynamic constraints

Display
- Displays parameters

MCMC
- Wraps parameters around fitting



markov_chain needs to know how many points are being fitted
which can be different...