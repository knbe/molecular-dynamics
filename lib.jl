# molecular dynamics of a gas of atoms in 2D
# this is a barebones reimplementation of the old python script
# units: m = ε = σ = k_B = 1

push!(LOAD_PATH, "./")
using Statistics
using Plots
using Vectors

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

function old_minimum_image(dx::Float64, dy::Float64, L::Float64)
	#halfL = 0.5 * sys.length
	#return (x .+ halfL) .% L .- halfL
	if dx > halfL
		dx -= L
	elseif dx < -halfL
		dx += L
	end

	if dy > halfL
		dy -= L
	elseif dy < -halfL
		dy += L
	end

	return dx, dy
end

function minimum_image(x::Vector{Float64}, L::Float64)
	halfL = 0.5 * L
	return (x .+ halfL) .% L .- halfL
end

function force(sys::ParticleSystem2D)
	# fairly sure this is more efficient than having a big old if statement
#	@static if @isdefined(lennardJones)
#		force, virial = lennardJonesForce(s)
#		s.virialAccum += virial
#		return force
#	end
	# add functionality for other forces
	if sys.forceType == "lennardJones"
		force, virial = lennard_jones_force(sys)
	elseif sys.forceType == "powerLaw"
		force, virial = power_law_force(sys)
	end

	sys.virialAccumulator += virial
	return force
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
	sys.x .+= sys.v .* sys.dt .+ 0.5 * (sys.dt)^2 .* a
	sys.x .%= sys.length
	sys.v .+= 0.5 * sys.dt .* (a .+ force(sys))
end

function evolve!(sys::ParticleSystem2D, time::Float64=10.0)
	steps = Int64(abs(time/sys.dt))

	for step in 1:steps
		verlet_step!(sys)
		zero_total_momentum!(sys)

		sys.t += sys.dt
		push!(sys.tPoints, sys.t)

		if (step % sys.sampleInterval == 0)
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

function zero_total_momentum(sys::ParticleSystem2D)
end

function reverse_time(sys::ParticleSystem2D)
end

function cool(sys::ParticleSystem2D, time::Float64=1.0)
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

	for i in 1:nx
		for j in 1:ny
			sys.x[2*(i*ny+j)-1] = (i - 0.5) * dx
			sys.x[2*(i*ny+j)] = (j - 0.5) * dy
		end
	end
end

function random_velocities!(sys::ParticleSystem2D)
	sys.v = rand(2*sys.numParticles) - 0.5
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
	return 1.0 + 0.5 * meanVirial / (sys.N * mean_temperature(sys))
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

function plot_positions(sys::ParticleSystem2D)
	N = sys.numParticles
	plot()
	scatter!(sys.x[1:2:2N], sys.x[2:2:2N])
	xlabel!("x")
	ylabel!("y")
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
end

function velocity_histogram(sys::ParticleSystem2D)
	plot()
	histogram!(sys.vPoints, normalize=:pdf)
	xlabel!("velocity")
	ylabel!("probability")
end

# RUN
################################################################################

sys = ParticleSystem2D(4)
lennard_jones_force(sys)
