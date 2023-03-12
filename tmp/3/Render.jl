module Render

@static if @isdefined(graph)
	using Plots
end
@static if @isdefined(render)
	using GLMakie
end

function start_plot()
	theme(:lime)
	plot(size=(800,800), titlefontsize=28)
end

function plot_positions(sys::ParticleSystem)
	N = sys.numParticles

	start_plot()
	scatter!(sys.x[1:2:2N], sys.x[2:2:2N], markersize=5.0,)
	xlabel!("x")
	ylabel!("y")
	xlims!(0, sys.length)
	ylims!(0, sys.length)
end

function plot_trajectories(sys::ParticleSystem, number::Int64=1)
	N = sys.numParticles
	#plot()
end

function plot_temperature(sys::ParticleSystem)
	plot()
	plot!(sys.tPoints, sys.tempPoints)
	xlabel!("t")
	ylabel!("T")
end

function plot_energy(sys::ParticleSystem)
	plot()
	plot!(sys.sampleTimePoints, sys.energyPoints)
	xlabel!("t")
	ylabel!("energy")
	ylims!(minimum(sys.energyPoints) - 2, maximum(sys.energyPoints) + 2)
end

function velocity_histogram(sys::ParticleSystem)
	plot()
	histogram!(sys.vPoints, normalize=:pdf)
	xlabel!("velocity")
	ylabel!("probability")
end

end
