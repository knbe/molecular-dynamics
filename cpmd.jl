# molecular dynamics in 2d.
# usage: 

using Statistics
using Plots

const global sampleInterval::Int64 = 100

mutable struct ParticleSystem
	N::Int64			# number of particles
	L::Float64			# square box side length
	T₀::Float64			# initial temperature
	t::Float64			# system time
	dt::Float64			# time step
	state::Vector{Float64}		# state space array
	steps::Int64			# number of steps
	timePoints::Vector{Float64}	# array of sampled time points
	energyPoints::Vector{Float64}	# array of sampled energy values
	tempPoints::Vector{Float64}	# array of sampled temperature values
	tempAccumulator::Float64	# temperature accumulator
	squareTempAccumulator::Float64	# T^2 accumulator
	virialAccumulator::Float64	# virial accumulator
	xPoints::Vector{Vector{Float64}} # array of sampled position data
	vPoints::Vector{Vector{Float64}} # array of sampled velocity data
	forceType::String		# force
end

function ParticleSystem(N::Int64=64, L::Float64=10.0, T₀::Float64=1.0)
	t = 0.0
	dt = 0.001
	state = zeros(4N) # state space array, [x1, y1, x2, y2, ..., vx1, vy1, ...]
	steps = 0
	timePoints = Float64[]
	energyPoints = Float64[]
	tempPoints = Float64[ T₀ ]
	tempAccumulator = 0.0
	squareTempAccumulator = 0.0
	virialAccumulator = 0.0
	xPoints = Vector{Float64}[]
	vPoints = Vector{Float64}[]
	forceType = "lennardJones"

	return ParticleSystem(
		N, 
		L, 
		T₀, 
		t, 
		dt, 
		state, 
		steps, 
		timePoints, 
		energyPoints, 
		tempPoints, 
		tempAccumulator, 
		squareTempAccumulator, 
		virialAccumulator, 
		xPoints, 
		vPoints, 
		forceType
	)
end

# some useful "views" of the state array
@views positions(state) = state[1:Int64(length(state)/2)]
@views velocities(state) = state[(Int64(length(state)/2)+1):end]
@views xcomponent(vector) = vector[1:2:end]
@views ycomponent(vector) = vector[2:2:end]
@views particle(n, vector) = [ vector[2n-1], vector[2n] ]

# ADDITIONAL FUNCTIONS

function dot(v1::Vector{Float64}, v2::Vector{Float64})
	return sum(v1 .* v2)
end

# INITIALIZATION
################################################################################

function random_positions!(sys::ParticleSystem)
	positions(sys.state) .= rand(2*sys.N) .* sys.L
	cool!(sys)
end

function triangular_lattice_positions!(sys::ParticleSystem)
end

function rectangular_lattice_positions!(sys::ParticleSystem)
	n = Int64(sqrt(sys.N))	# num lattice points per axis
	ny = nx			# square lattice
	ax = sys.L / nx		# lattice spacing x
	ay = sys.L / ny		# lattice spacing y

	for i in 0:(nx-1)
		for j in 0:(ny-1)
			sys.state[2*(i*ny+j)+1] = (i + 0.5) * ax
			sys.state[2*(i*ny+j)+2] = (j + 0.5) * ay
		end
	end
end

function random_velocities!(sys::ParticleSystem)
	velocities(sys.state) .= rand(2*sys.N) .- 0.5
	zero_total_momentum!(sys)
	velocities(sys.state) .*= sqrt(sys.T₀/temperature(sys))
end

function zero_total_momentum!(sys::ParticleSystem)
	xcomponent(velocities(sys.state)) .-= 
		mean(xcomponent(velocities(sys.state)))
	ycomponent(velocities(sys.state)) .-= 
		mean(ycomponent(velocities(sys.state)))
end


# FORCES
################################################################################

#function minimum_image(x::Vector{Float64}, L::Float64)
#	halfL = 0.5 * L
#	for i in 1:length(x)
#		if x[i] > halfL
#			x[i] -= L
#		elseif x[i] < -halfL
#			x[i] += L
#		end
#	end
#	return x
#end

function force(sys::ParticleSystem)
#	if sys.forceType == "lennardJones"
#		force, virial = lennard_jones_force(sys)
#	elseif sys.forceType == "powerLaw"
#		force, virial = power_law_force(sys)
#	end

	force, virial = lennard_jones_force(sys)
	sys.virialAccumulator += virial

	return force
end

function minimum_image(xij::Float64, L::Float64)
	if xij > (L/2)
		xij -= L
	elseif xij < -(L/2)
		xij += L
	end
	return xij
end

function lennard_jones_force(sys::ParticleSystem)
	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	virial = 0.0
	force = zeros(2*sys.N)

	Threads.@threads for i = 1:(sys.N-1)
		for j = (i+1):sys.N
			dx = minimum_image(x[i] - x[j], sys.L)
			dy = minimum_image(y[i] - y[j], sys.L)

			r2inv = 1.0 / (dx^2 + dy^2)
			f = 48.0 * r2inv^7 - 24.0 * r2inv^4
			fx = dx * f
			fy = dy * f

			force[2*i-1] += fx
			force[2*i] += fy
			force[2*j-1] -= fx
			force[2*j] -= fy

			virial += fx * dx + fy * dy
		end
	end

	return force, 0.5 * virial
end

function lennard_jones_potential(sys::ParticleSystem)
	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	U = 0.0

	Threads.@threads for i in 1:(sys.N-1)
		for j in (i+1):sys.N
			dx = minimum_image(x[i] - x[j], sys.L)
			dy = minimum_image(y[i] - y[j], sys.L)

			r2inv = 1.0 / (dx^2 + dy^2)
			U += r2inv^6 - r2inv^3
		end
	end
	return 4.0 * U
end

# TIME EVOLUTION
################################################################################

function keep_particles_in_box!(sys::ParticleSystem)
	for i in 1:(2*sys.N)
		if positions(sys.state)[i] > sys.L
			positions(sys.state)[i] -= sys.L
		elseif positions(sys.state)[i] < 0.0
			positions(sys.state)[i] += sys.L
		end
	end
end

function verlet_step!(sys::ParticleSystem)
	acceleration = force(sys)

	positions(sys.state) .+= 
		velocities(sys.state) .* sys.dt .+ 0.5 .* acceleration .* (sys.dt)^2

	keep_particles_in_box!(sys)

	velocities(sys.state) .+= 
		0.5 * sys.dt .* (acceleration + force(sys))
end

function evolve!(sys::ParticleSystem, runtime::Float64=10.0)
	numsteps = Int64(abs(runtime/sys.dt) + 1)

	print_evolution_message(runtime, numsteps)

	@time for step in 1:numsteps
		verlet_step!(sys)
		zero_total_momentum!(sys)

		if (step % sampleInterval == 1)
			push!(sys.timePoints, sys.t)
			push!(sys.energyPoints, energy(sys))
			push!(sys.xPoints, positions(sys.state))
			push!(sys.vPoints, velocities(sys.state))

			T = temperature(sys)
			push!(sys.tempPoints, T)
			sys.tempAccumulator += T
			sys.squareTempAccumulator += T^2
		end

		sys.t += sys.dt
		sys.steps += 1
	end
	println("done.") 
	print_data(sys)
end

function reverse_time!(sys)
	sys.dt *= -1
	println("\ntime reversed! dt = $(sys.dt)")
end

function cool!(sys::ParticleSystem, cooltime::Float64=1.0)
	numsteps = Int64(cooltime/sys.dt)
	for step in 1:numsteps
		verlet_step!(sys)
		velocities(sys.state) .*= (1.0 - sys.dt)
	end
	reset_statistics!(sys)
end

# MEASUREMENTS
################################################################################

function kinetic_energy(sys::ParticleSystem)
	return 0.5 * sum(velocities(sys.state) .* velocities(sys.state))
end

function potential_energy(sys::ParticleSystem)
	return lennard_jones_potential(sys)
end

function temperature(sys::ParticleSystem)
	return kinetic_energy(sys) / sys.N
end

function energy(sys::ParticleSystem)
	return potential_energy(sys) + kinetic_energy(sys)
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
	plot(
		size=(800,800), 
		titlefontsize=16, 
		guidefontsize=16,
	)
end

function plot_positions(sys::ParticleSystem)
	initialize_plot()
	scatter!(
		xcomponent(positions(sys.state)), 
		ycomponent(positions(sys.state)), 
		markersize = 5.0,
		grid = true,
		widtn = true,
		)
	xlims!(0, sys.L)
	ylims!(0, sys.L)
	xlabel!("x")
	ylabel!("y")
	title!("particle positions at time t=$(sys.t)")
end

function plot_trajectories(sys::ParticleSystem, number::Int64=1)
	#N = sys.N
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
	histogram!(sys.vPoints, normalize=:pdf, bins = -2.0:0.2:2.0)
	xlabel!("velocity")
	ylabel!("probability")
end

# CONSOLE PRINT/DEBUGGING
################################################################################

function print_init_message(sys::ParticleSystem)
	println("\nmolecular dynamics.")
	println("number of threads: ", Threads.nthreads()),

	println("\nsystem parameters:")
	println("\tN =  $(sys.N)   (number of particles)")
	println("\tL =  $(sys.L)   (side length of square box)")
	println("\tDT = $(sys.dt)  (time step)")

	println("\ninitial conditions:")
	println("\ttime:            $(sys.t)")
	println("\tenergy:          $(energy(sys))")
	println("\tkinetic energy:  $(kinetic_energy(sys))")
	println("\ttemperature:     $(sys.T₀)")
end

function print_evolution_message(runtime, numsteps)
	println("\nevolving...")
	#println("\trun time: $runtime")
	#println("\tnum steps: $numsteps")
end

function print_data(sys::ParticleSystem)
	println("\nmeasurements:")
	println("\ttime:            $(sys.t)")
	println("\tsteps evolved:   $(sys.steps)")
	println("\tenergy:          $(energy(sys))")
	println("\tkinetic energy:  $(kinetic_energy(sys))")
	println("\tmean energy:     $(mean_energy(sys))")
	println("\tstd energy:      $(std_energy(sys))")
	println("\ttemperature:     $(temperature(sys))")
	println("\theat capacity:   $(heat_capacity(sys))")
	println("\tPV/NkT:          $(mean_pressure(sys))")
end

# RUN
################################################################################


#sys = ParticleSystem(64, 8.0, 1.0)

#function plot_lj_potential(ϵ=1.0, σ=1.0)
#	r = Vector{Float64}(0.1:0.01:10.0)
#	U = ljU.(r, ϵ=ϵ, σ=σ)
#	initialize_plot()
#	plot!(r,U)
#	xlims!(0,10)
#	ylims!(-ϵ-1, 30)
#	title!("LJ potential with ϵ=$ϵ, σ=$σ")
#	xlabel!("r")
#	ylabel!("U(r)")
#	plot!(
#		[2.5σ, 2.5σ], [-ϵ-1,60], 
#		linestyle = :dot, 
#		label = "2.5σ", 
#	)
#end

#function demo_0()
#	sys = ParticleSystem(10000, 100.0, 1.0)
#	rectangular_lattice_positions!(sys)
#	random_velocities!(sys)
#
#	println("ljp1")
#	@time begin
#	U = lennard_jones_potential(sys.state, sys.L)
#	println(U)
#	end
#
#	println("\nljp2")
#	@time begin
#	U2 = ljp2(sys)
#	end
#	println(U2)
#
#	println("\nljp3")
#	@time begin
#	U3 = ljp3(sys)
#	end
#	println(U3)
#end

#sys = ParticleSystem(16, 100.0, 1.0)
#rectangular_lattice_positions!(sys)
#random_velocities!(sys)


function testlj()
	rectangular_lattice_positions!(sys)
	random_velocities!(sys)
	lennard_jones_potential2(sys)
	evolve!(sys, 1.0)
end

sys = ParticleSystem(64, 130.0, 1.0)

# DEMO 0: EQUILIBRIATE
function demo_0()
	rectangular_lattice_positions!(sys)
	random_velocities!(sys)
	print_init_message(sys)

	evolve!(sys, 10.0)
	plot_positions(sys)
end

# DEMO 1: TIME REVERSAL TEST
function demo_1()
	rectangular_lattice_positions!(sys)
	random_velocities!(sys)
	print_init_message(sys)

	evolve!(sys, 10.0)
	reverse_time(sys)
	evolve!(sys, 10.0)
	plot_positions(sys)
end

# EQUILIBRIATION AND STATISTICS
function demo2()
	random_positions!(sys)
	random_velocities!(sys)
	evolve!(sys, 1.0)
	console_log(sys)
	p1 = plot_positions(sys)
	p2 = plot_energy(sys)
	p3 = velocity_histogram(sys)
	plot(p1, p2, p3)
end
