module Graphs

using Plots

function initialize_plot()
	#theme(:lime)
	plot(size=(800,800), titlefontsize=28)
end

function plot_positions(sys::ParticleSystem)
	N = sys.N

	initialize_plot()
	scatter!(
		component(1, positions(sys.state)), 
		component(2, positions(sys.state)), 
		markersize=5.0,
		)
	xlabel!("x")
	ylabel!("y")
	xlims!(0, sys.L)
	ylims!(0, sys.L)
end

end
