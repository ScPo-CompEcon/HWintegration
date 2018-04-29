


include("src/HWintegration.jl")	# load our code
HWintegration.runall()

#if plot didn't run on its own
using Plots
default(size=(5000,5000))
HWintegration.plot_q1()
plot(HWintegration.plotall...)
