# molecular dynamics in 2d

push!(LOAD_PATH, "./")
using Statistics
using Plots

dt = 0.001
sampleInterval = 100

mutable struct ParticleSystem
	N::Int64			# number of particles
	L::Float64			# length of square side
	T₀::Float64			# initial temperature
	t::Float64			# system time
	state::Vector{Float64}		# state space array
	steps::Int64
	tPoints::Vector{Float64}	# array of time steps added to during integration
	sampleTimePoints::Vector{Float64}
	energyPoints::Vector{Float64}	# array of sampled energy points
	tempPoints::Vector{Float64}
	tempAccumulator::Float64
	squareTempAccumulator::Float64
	virialAccumulator::Float64
	xPoints
	vPoints
	forceType::String
end

function ParticleSystem(N::Int64=64, L::Float64=10.0, T₀::Float64=1.0)
	t = 0.0
	tPoints = [ 0.0 ]
	steps = 0
	energyPoints = []
	sampleTimePoints = []
	tempPoints = [ T₀ ]
	tempAccumulator = 0.0
	squareTempAccumulator = 0.0
	virialAccumulator = 0.0
	state = zeros(4N)	# state space array, [x1, y1, x2, y2, ..., vx1, vy1, ...]
	xPoints = []
	vPoints = []
	forceType = "lennardJones"

	return ParticleSystem(
		N, L, T₀, t, state, 
		steps, tPoints, sampleTimePoints, energyPoints, 
		tempPoints, tempAccumulator, squareTempAccumulator, 
		virialAccumulator, xPoints, vPoints, forceType
		)
end

# some useful "views" of the state array
@views positions(state) = state[1:Int64(length(state)/2)]
@views velocities(state) = state[(Int64(length(state)/2)+1):end]
@views xcomponent(vector) = vector[1:2:end]
@views ycomponent(vector) = vector[2:2:end]
@views particle(n, vector) = [ vector[2n-1], vector[2n] ]
@views numParticles(state) = Int64(length(state) / 4)

# INITIALIZATION
################################################################################

function random_positions!(sys::ParticleSystem)
	positions(sys.state) .= rand(2*sys.N) .* sys.L
	cool(sys)
end

function triangular_lattice_positions!(sys::ParticleSystem)
end

function rectangular_lattice_positions!(sys::ParticleSystem)
	nx = Int64(sqrt(sys.N))		# num lattice points per axis
	ny = nx				# square lattice
	ax = sys.L / nx		# lattice point separation
	ay = sys.L / ny

	for i in 0:(nx-1)
		for j in 0:(ny-1)
			sys.state[2*(i*ny+j)+1] = (i + 0.5) * ax
			sys.state[2*(i*ny+j)+2] = (j + 0.5) * ay
		end
	end
end

function random_velocities!(sys::ParticleSystem)
	velocities(sys.state) .= rand(2*sys.N) .- 0.5
	zero_total_momentum!(velocities(sys.state))
	T = temperature(sys.state)
	velocities(sys.state) .*= sqrt(sys.T₀/T)
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
	return x
end

function force(state::Vector{Float64}, L::Float64)
#	if sys.forceType == "lennardJones"
#		force, virial = lennard_jones_force(sys)
#	elseif sys.forceType == "powerLaw"
#		force, virial = power_law_force(sys)
#	end

	force = lennard_jones_force(state, L)
	return force
end

function lennard_jones_force(state::Vector{Float64}, L::Float64)
	N = numParticles(state)
	tiny = 1.0e-40
	virial = 0.0
	force = zeros(2N)

	x = xcomponent(positions(state))
	y = ycomponent(positions(state))

	Threads.@threads for i in 1:N
		for j in (i+1):N
			dx = minimum_image(x .- x[i], L)
			dy = minimum_image(y .- y[i], L)

			r2inv = 1.0 ./ (dx.^2 .+ dy.^2 .+ tiny)
			f = 48.0 .* r2inv.^6 - 24.0 .* r2inv.^4
			fx = dx .* f
			fy = dy .* f

			fx[i] = 0.0	# self force is zero
			fy[i] = 0.0

			#virial += dot(fx,dx) + dot(fy,dy)
		end
	end

	return force
end

function lennard_jones_potential(state::Vector{Float64}, L::Float64)

	N = numParticles(state)
	tiny = 1.0e-40

	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	U = 0.0

	for i in 1:N
		dx = minimum_image(x[i] .- x, L)
		dy = minimum_image(y[i] .- y, L)

		r2inv = 1.0 / (dx.^2 .+ dy.^2 .+ tiny)
		dU = r2inv.^6 .- r2inv.^3
		dU[i] = 0.0	# self interaction is zero
		U += sum(dU)
	end

	return 2.0 * U
end


# TIME EVOLUTION
################################################################################

function verlet_step!(state::Vector{Float64}, L::Float64)
	acceleration = force(state, L)
	positions(state) .+= velocities(state) .* dt .+ 0.5 * (dt)^2 .* acceleration
	positions(state) .%= L
	velocities(state) .+= 0.5 * dt .* (acceleration .+ force(state, L))
end

function evolve!(sys::ParticleSystem, time::Float64=10.0)
	steps = Int64(abs(time/dt))

	for step in 1:steps
		verlet_step!(sys.state, sys.L)
		zero_total_momentum!(velocities(sys.state))

		sys.t += dt
		push!(sys.tPoints, sys.t)

		if (step % sampleInterval == 1)
			push!(sys.energyPoints, energy(sys.state, sys.L))
			push!(sys.sampleTimePoints, sys.t)
			push!(sys.xPoints, positions(sys.state))
			push!(sys.vPoints, velocities(sys.state))

		end

		T = temperature(sys.state)
		sys.steps += 1
		push!(sys.tempPoints, T)
		sys.tempAccumulator += T
		sys.squareTempAccumulator += T^2
	end
end

function zero_total_momentum!(velocities)
	xcomponent(velocities) .-= mean(xcomponent(velocities))
	ycomponent(velocities) .-= mean(ycomponent(velocities))
end

function reverse_time()
	global dt *= -1
end

function cool(sys::ParticleSystem, time::Float64=1.0)
	steps = Int64(time/dt)
	for step in 1:steps
		verlet_step!(sys.state, sys.L)
		velocities(sys.state) .*= (1.0 - dt)
	end

	reset_statistics!(sys)
end



# MEASUREMENTS
################################################################################

function kinetic_energy(state::Vector{Float64})
	return 0.5 * sum(velocities(state) .* velocities(state))
end

function potential_energy(state::Vector{Float64}, L::Float64)
	return lennard_jones_potential(state, L)
end

function temperature(state::Vector{Float64})
	return kinetic_energy(state) / numParticles(state)
end

function energy(state::Vector{Float64}, L::Float64)
	return potential_energy(state, L) + kinetic_energy(state)
end

# STATISTICS
################################################################################

function reset_statistics!(sys::ParticleSystem)
	sys.steps = 0
	sys.tempAccumulator = 0.0
	sys.squareTempAccumulator = 0.0
	sys.virialAccumulator = 0.0
	sys.xPoints = []
	sys.vPoints = []
end

function mean_temperature(sys::ParticleSystem)
	return sys.tempAccumulator / sys.steps
end

function mean_square_temperature(sys::ParticleSystem)
	return sys.squareTempAccumulator / sys.steps
end

function mean_pressure(sys::ParticleSystem)
	# factor of half because force is calculated twice each step
	meanVirial = 0.5 * sys.virialAccumulator / sys.steps
	return 1.0 + 0.5 * meanVirial / (sys.N * mean_temperature(sys))
end

function heat_capacity(sys::ParticleSystem)
	meanTemperature = mean_temperature(sys)
	meanSquareTemperature = mean_square_temperature(sys)
	σ2 = meanSquareTemperature - meanTemperature^2
	denom = 1.0 - σ2 * sys.N / meanTemperature^2
	return sys.N / denom
end

function mean_energy(sys::ParticleSystem)
	return mean(sys.energyPoints)
end

function std_energy(sys::ParticleSystem)
	return std(sys.energyPoints)
end

# GRAPHS
################################################################################

function initialize_plot()
	#theme(:lime)
	plot(size=(800,800), titlefontsize=28)
end

function plot_positions(sys::ParticleSystem)
	N = sys.N

	initialize_plot()
	scatter!(
		xcomponent(positions(sys.state)), 
		ycomponent(positions(sys.state)), 
		markersize=5.0,
		)
	xlabel!("x")
	ylabel!("y")
	xlims!(0, sys.L)
	ylims!(0, sys.L)
end

function plot_trajectories(sys::ParticleSystem, number::Int64=1)
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
	histogram!(sys.vPoints, normalize=:pdf, bins = -2.0:0.2:2.0)
	xlabel!("velocity")
	ylabel!("probability")
end

function console_log(sys::ParticleSystem)
	println(  "  num threads:		", Threads.nthreads(), "\n"),
	println(
		  "  time:			", sys.t, 
		"\n  total energy:		", energy(sys.state, sys.L), 
		"\n  temperature:		", temperature(sys.state),
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


sys = ParticleSystem(64, 8.0, 1.0)

# time reversal
function demo1()
	rectangular_lattice_positions!(sys)
	random_velocities!(sys)
	evolve!(sys, 1.0)
	reverse_time()
	evolve!(sys, 1.0)
	console_log(sys)
	plot_positions(sys)
end

# equilibration
function demo2()
	random_positions!(sys)
	random_velocities!(sys)
	evolve!(sys, 1.0)
	console_log(sys)
	p1 = plot_positions(sys)
	p2 = plot_energy(sys)
	#p3 = plot_temperature(sys)
	p4 = velocity_histogram(sys)
	plot(p1, p2, p4)
end
