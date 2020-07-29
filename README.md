# metronome-phase-model
This repo has matlab/numpy/julia code for all the simulations run by me as part of a phase model - based project at BITS Goa.

Descriptions of the files:

Two Metronomes:	Holds all the function files for studying the two metronome system
	
	rk42d.m          : custom RK4 integrator
	findLimitCycle.m : finds limit cycle in the phase space for the system of equations passed to the function
	malkin.m         : finds Q(t) for the van Der Pol oscillator or Andronov Hopf oscillator using Malkin's approach
	metSolver2.m     : solves the full system of equations for the two coupled metronome system
	metMalkin.m      : finds the Q(t) of the two coupled metronomes system with a small angle approximation (using Malkin's approach)
	fullMetMalkin.m  : find the Q(t) of the two coupled metronomes system withOUT the small angle approximation
	metH.m           : finds the H(chi) and G(chi) for the system

Chimera finder: for finding chimera states in the swing-metronomes system
	
	chimera.jl, chimera.m, chimeraSolver.py:
		All three numerically solve the ODEs of the system
	twoChim.jl, swing_metro.m:
		Solve the single swing, one/multiple oscillator system ODEs
		
The MATLAB files will all run natively on a MATLAB IDE.
The Python files have dependencies on numpy and matplotlib libraries.
The Julia files have dependencies on PyPlot (the matplotlib port for Julia) and LinearAlgebra libraries.

The first folder has very bad coding and documentation practices. I apologize for not commenting enough in them.
These days even I'm not able to understand what's happening in them. In the future when I get time I may add comments.

The Julia files are the most recent scripts, but just like the matlab files they have virtually no documentation.
Please see the Python files for any sort of documentation.

The PDF file of my final report is available in this repo. Please check it for all the explanations involved.
Prerequisite reading:
	Izhikevich (2007)
	Pantaleone (2002)
	Martens et. al. (2013)

Please contact me for any clarifications.