# molecular dynamics of a gas of atoms in 2D
# this is a barebones reimplementation of the old python script
# units: m = ε = σ = k_B = 1

push!(LOAD_PATH, "./")
using Statistics
using Vectors

using Plots
using PlotThemes

#@static if @isdefined(render)
#	using GLMakie
#end

dt = 0.001
sampleInterval = 100

mutable struct ParticleSystem2D
	numParticles::Int64		# number of particles
	length::Float64			# length of square side
	initTemp::Float64
	t::Float64			# initial time
	tPoints::Vector{Float64}	# array of time steps added to during integration
	steps::Int64
	energyPoints::Vector{Float64}	# list of energy, sampled every sampleInterval steps
	sampleTimePoints::Vector{Float64}
	tempPoints::Vector{Float64}
	tempAccumulator::Float64
	squareTempAccumulator::Float64
	virialAccumulator::Float64
	x::Vector{Float64}		# array of N (x,y) positions
	v::Vector{Float64}		# array of N (vx,vy) velocities
	xPoints
	vPoints
	forceType::String
end

function ParticleSystem2D(numParticles::Int64=64, 
		length::Float64=10.0, initTemp::Float64=1.0)
	t = 0.0
	tPoints = [ 0.0 ]
	steps = 0
	energyPoints = []
	sampleTimePoints = []
	tempPoints = [ initTemp ]
	tempAccumulator = 0.0
	squareTempAccumulator = 0.0
	virialAccumulator = 0.0
	x = zeros(2 * numParticles)
	v = zeros(2 * numParticles)
	xPoints = []
	vPoints = []
	forceType = "lennardJones"

	return ParticleSystem2D(
		numParticles, length, initTemp, t, tPoints, steps, 
		energyPoints, sampleTimePoints, tempPoints, 
		tempAccumulator, squareTempAccumulator, 
		virialAccumulator, x, v, xPoints, vPoints, forceType
		)
end

# FORCES
################################################################################

function minimum_image(x::Vector{Float64}, L::Float64)
	halfL = 0.5 * L
	for i in 1:length(x)
		if x[i] > halfL
			x[i] -= L
		elseif x[i] < -halfL
			x[i] += L
		end
	end
	#return (x .+ halfL) .% L .- halfL
	return x
end

function force(sys::ParticleSystem2D)
	if sys.forceType == "lennardJones"
		force, virial = lennard_jones_force2(sys)
	elseif sys.forceType == "powerLaw"
		force, virial = power_law_force(sys)
	end

	sys.virialAccumulator += virial
	return force
end

function minimum_separation_vector(p1, p2, L::Float64)
	r12 = p2 .- p1 # vector from point p1 to point p2

	for i in 1:length(r12)
		if r12[i] > (L/2)
			r12[i] -= L
		elseif r12[i] < -(L/2)
			r12[i] += L
		end
	end

	return r12
end

function lennard_jones_force2(sys::ParticleSystem2D)
	N = sys.numParticles
	L = sys.length
	tiny = 1.0e-40
	virial = 0.0

	x = sys.x[1:2:2N]
	y = sys.x[2:2:2N]
	force = zeros(2*N)

	#minImg = minimumImage(s, x)
	#println(minImg)
	for i in 1:N
		for j in (i+1):N
			dx = minimum_image(x .- x[i], L)
			dy = minimum_image(y .- y[i], L)

			r2inv = 1.0 ./ (dx.^2 .+ dy.^2 .+ tiny)
			f = 48.0 .* r2inv.^7 - 24.0 .* r2inv.^4
			fx = dx .* f
			fy = dy .* f

			fx[i] = 0.0	# no self force
			fy[i] = 0.0


#			println("i = ", i, ", j = ", j)
#			println(dx)
#			println(dy)
#			println(f)
#			println(sum(fx))
#			println(sum(fy))
			#force[2*i-1] = sum(fx)
			#force[2*i] = sum(fy)
			#println(force)

			virial += dot(fx,dx) + dot(fy,dy)
		end
	end

	return force, (0.5 * virial)
end

function lennard_jones_force(sys::ParticleSystem2D)
	N = sys.numParticles
	L = sys.length
	tiny = 1.0e-40
	virial = 0.0

	x = sys.x[1:2:2N]
	y = sys.x[2:2:2N]
	force = zeros(2N)

	#minImg = minimumImage(s, x)
	#println(minImg)
	for i in 1:N
		for j in (i+1):N
			dx = minimum_image(x[i] .- x, L)
			dy = minimum_image(y[i] .- y, L)

			r2inv = 1.0 ./ (dx.^2 .+ dy.^2 .+ tiny)
			c = 48.0 .* r2inv.^7 - 24.0 .* r2inv.^4
			fx = dx .* c
			fy = dy .* c

			fx[i] = 0.0	# no self force
			fy[i] = 0.0
			force[2i-1] = sum(fx)
			force[2i] = sum(fy)

#			println(fx)
#			println(dx)
#			println()
			virial += dot(fx,dx) + dot(fy,dy)
		end
	end

	return force, (0.5 * virial)
end

function power_law_force(sys::ParticleSystem2D)
	N = sys.numParticles
	L = sys.length
	tiny = 1.0e-40
	virial = 0.0

	x = sys.x[1:2:2N]
	y = sys.x[2:2:2N]
	force = zeros(2N)

	for i in 1:N
		for j in (i+1):N
			dx = minimum_image(x[i] .- x, L)
			dy = minimum_image(y[i] .- y, L)

			r2 = dx.^2 .+ dy.^2 .+ tiny
			r6inv = r2.^(-3)
			fx = dx .* r6inv
			fy = dy .* r6inv

			fx[i] = 0.0	# no self force
			fy[i] = 0.0
			force[2i-1] = sum(fx)
			force[2i] = sum(fy)

#			println(fx)
#			println(dx)
#			println()
			virial += dot(fx,dx) + dot(fy,dy)
		end
	end

	return force, (0.5 * virial)
end

# TIME EVOLUTION
################################################################################

function verlet_step!(sys::ParticleSystem2D)
	a = force(sys)
	sys.x .+= sys.v .* dt .+ 0.5 * (dt)^2 .* a
	sys.x .%= sys.length
	sys.v .+= 0.5 * dt .* (a .+ force(sys))
end

function evolve!(sys::ParticleSystem2D, time::Float64=10.0)
	steps = Int64(abs(time/dt))

	for step in 1:steps
		verlet_step!(sys)
		zero_total_momentum!(sys)

		sys.t += dt
		push!(sys.tPoints, sys.t)

		if (step % sampleInterval == 0)
			push!(sys.energyPoints, energy(sys))
			push!(sys.sampleTimePoints, sys.t)
			push!(sys.xPoints, sys.x)
			push!(sys.vPoints, sys.v)
		end

		T = temperature(sys)
		sys.steps += 1
		push!(sys.tempPoints, T)
		sys.tempAccumulator += T
		sys.squareTempAccumulator += T^2
	end
end

function zero_total_momentum!(sys::ParticleSystem2D)
	N = sys.numParticles
	vx = sys.v[1:2:2N]
	vy = sys.v[2:2:2N]

	vx .-= mean(vx)
	vy .-= mean(vy)

	sys.v[1:2:2N] = vx
	sys.v[2:2:2N] = vy
end

function reverse_time(sys::ParticleSystem2D)
	dt = -dt
end

function cool(sys::ParticleSystem2D, time::Float64=1.0)
	steps = Int64(time/dt)
	for step in 1:steps
		verlet_step!(sys)
		sys.v .*= (1.0 - dt)
	end

	reset_statistics!(sys)
end

# INITIAL CONDITIONS
################################################################################

function random_positions!(sys::ParticleSystem2D)
	sys.x = rand(2*sys.numParticles) .* sys.length

	sys.forceType = "lennardJones"
	cool(sys)
end

function triangular_lattice_positions!(sys::ParticleSystem2D)
	random_positions!(sys)
	sys.v .+= rand(2*sys.numParticles) .- 0.5

	cool(sys, 10.0)
	sys.forceType = "lennardJones"
end

function rectangular_lattice_positions!(sys::ParticleSystem2D)
	nx = Int64(sqrt(sys.numParticles))
	ny = nx
	dx = sys.length / nx
	dy = sys.length / ny

	for i in 0:(nx-1)
		for j in 0:(ny-1)
			sys.x[2*(i*ny+j)+1] = (i + 0.5) * dx
			sys.x[2*(i*ny+j)+2] = (j + 0.5) * dy
		end
	end
end

function random_velocities!(sys::ParticleSystem2D)
	sys.v = rand(2*sys.numParticles) .- 0.5
	zero_total_momentum!(sys)
	T = temperature(sys)
	sys.v .*= sqrt(sys.initTemp/T)
end

# MEASUREMENTS
################################################################################

function kinetic_energy(sys::ParticleSystem2D)
	return 0.5 * mag_sq(sys.v)
end

function potential_energy(sys::ParticleSystem2D)
	return lennard_jones_potential(sys)
end

function lennard_jones_potential(sys::ParticleSystem2D)
	tiny = 1.0e-40
	L = sys.length
	halfL = 0.5 * L
	N = sys.numParticles

	x = sys.x[1:2:2N]
	y = sys.x[2:2:2N]
	U = 0.0

	for i in 1:N
		dx = minimum_image(x[i] .- x, L)
		dy = minimum_image(y[i] .- y, L)

		r2inv = 1.0 / (dx.^2 .+ dy.^2 .+ tiny)
		dU = r2inv.^6 .- r2inv.^3
		dU[i] = 0.0	# no self interaction
		U += sum(dU)
	end

	return 2.0 * U
end

function energy(sys::ParticleSystem2D)
	return potential_energy(sys) + kinetic_energy(sys)
end

function temperature(sys::ParticleSystem2D)
	return kinetic_energy(sys) / sys.numParticles
end

# STATISTICS
################################################################################

function reset_statistics!(sys::ParticleSystem2D)
	sys.steps = 0
	sys.tempAccumulator = 0.0
	sys.squareTempAccumulator = 0.0
	sys.virialAccumulator = 0.0
	sys.xPoints = []
	sys.vPoints = []
end

function mean_temperature(sys::ParticleSystem2D)
	return sys.tempAccumulator / sys.steps
end

function mean_square_temperature(sys::ParticleSystem2D)
	return sys.squareTempAccumulator / sys.steps
end

function mean_pressure(sys::ParticleSystem2D)
	# factor of half because force is calculated twice each step
	meanVirial = 0.5 * sys.virialAccumulator / sys.steps
	return 1.0 + 0.5 * meanVirial / (sys.numParticles * mean_temperature(sys))
end

function heat_capacity(sys::ParticleSystem2D)
	meanTemperature = mean_temperature(sys)
	meanSquareTemperature = mean_square_temperature(sys)
	σ2 = meanSquareTemperature - meanTemperature^2
	denom = 1.0 - σ2 * sys.numParticles / meanTemperature^2
	return sys.numParticles / denom
end

function mean_energy(sys::ParticleSystem2D)
	return mean(sys.energyPoints)
end

function std_energy(sys::ParticleSystem2D)
	return std(sys.energyPoints)
end

# GRAPHS
################################################################################

function start_plot()
	theme(:lime)
	plot(size=(800,800), titlefontsize=28)
end

function plot_positions(sys::ParticleSystem2D)
	N = sys.numParticles

	start_plot()
	scatter!(sys.x[1:2:2N], sys.x[2:2:2N], markersize=5.0,)
	xlabel!("x")
	ylabel!("y")
	xlims!(0, sys.length)
	ylims!(0, sys.length)
end

function plot_trajectories(sys::ParticleSystem2D, number::Int64=1)
	N = sys.numParticles
	#plot()
end

function plot_temperature(sys::ParticleSystem2D)
	plot()
	plot!(sys.tPoints, sys.tempPoints)
	xlabel!("t")
	ylabel!("T")
end

function plot_energy(sys::ParticleSystem2D)
	plot()
	plot!(sys.sampleTimePoints, sys.energyPoints)
	xlabel!("t")
	ylabel!("energy")
	ylims!(minimum(sys.energyPoints) - 2, maximum(sys.energyPoints) + 2)
end

function velocity_histogram(sys::ParticleSystem2D)
	plot()
	histogram!(sys.vPoints, normalize=:pdf)
	xlabel!("velocity")
	ylabel!("probability")
end

function console_log(sys::ParticleSystem2D)
	println(
		  "  time:			", sys.t, 
		"\n  total energy:		", energy(sys), 
		"\n  temperature:		", temperature(sys),
		)

	if sys.steps > 0
		println(
			"\n  mean energy:		", mean_energy(sys), 
			"\n  standard dev:		", std_energy(sys),
			"\n  C_V:			", heat_capacity(sys),
			"\n  PV/NkT:			", mean_pressure(sys),
			)
	end

end

# RUN
################################################################################

sys = ParticleSystem2D(64, 8.0, 1.0)
#lennard_jones_force(sys)

# equilibriation and statistics
random_positions!(sys)
random_velocities!(sys)
console_log(sys)

evolve!(sys, 1.0)
console_log(sys)
plot_positions(sys)
plot_energy(sys)

#a, b = lennard_jones_force2(sys)

#reset_statistics!(sys)
#evolve!(sys, 20.0)
#plot_energy(sys)
#plot_temperature(sys)
