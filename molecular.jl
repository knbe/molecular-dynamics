
using Statistics
using StatsPlots

mutable struct ParticleSystem
	N::Int64			# number of particles
	L::Float64			# square box side length
	T₀::Float64			# initial temperature
	t::Float64			# system time
	dt::Float64			# time step
	state::Vector{Float64}		# state space array
	steps::Int64			# number of steps
	sampleInterval::Int64		# interval for sampling data
	timeData::Vector{Float64}	# array of sampled time points
	energyData::Vector{Float64}	# array of sampled energy values
	tempData::Vector{Float64}	# array of sampled temperature values
	tempAccumulator::Float64	# temperature accumulator
	squareTempAccumulator::Float64	# T^2 accumulator
	virialAccumulator::Float64	# virial accumulator
	xData::Vector{Vector{Float64}} # array of sampled position data
	vData::Vector{Vector{Float64}} # array of sampled velocity data
	forceType::String		# force
end

function ParticleSystem(N::Int64=64, L::Float64=10.0, T₀::Float64=1.0)
	t = 0.0
	dt = 0.001
	state = zeros(4N) # state space array, [x1,y1,x2,y2,...,vx1,vy1,...]
	steps = 0
	timeData = Float64[]
	energyData = Float64[]
	sampleInterval = 100
	tempData = Float64[]
	tempAccumulator = 0.0
	squareTempAccumulator = 0.0
	virialAccumulator = 0.0
	xData = Vector{Float64}[]
	vData = Vector{Float64}[]
	forceType = "lennardJones"

	return ParticleSystem(
		N, 
		L, 
		T₀, 
		t, 
		dt, 
		state, 
		steps, 
		sampleInterval,
		timeData, 
		energyData, 
		tempData, 
		tempAccumulator, 
		squareTempAccumulator, 
		virialAccumulator, 
		xData, 
		vData, 
		forceType
	)
end

# some useful "views" of the state array
# (read the performance tips chapter of the julia manual)
@views positions(state) = state[ 1:Int64(length(state)/2) ]
@views velocities(state) = state[ (Int64(length(state)/2)+1):end ]
@views xcomponent(vector) = vector[ 1:2:end ]
@views ycomponent(vector) = vector[ 2:2:end ]
@views particle(n, vector) = [ vector[2n-1], vector[2n] ]



function keep_particles_in_box!(sys::ParticleSystem)
	for i in 1:(2*sys.N)
		if positions(sys.state)[i] > sys.L
			positions(sys.state)[i] -= sys.L
		elseif positions(sys.state)[i] < 0.0
			positions(sys.state)[i] += sys.L
		end
	end

#	# another way of doing this: use the ternary operator
#	for i in 1:(2 * sys.N)
#		positions(sys.state)[i] < 0.0 ?  
#		positions(sys.state)[i] % sys.L + sys.L : 
#		positions(sys.state)[i] % sys.L
#	end
end

function verlet_step!(sys::ParticleSystem)
	# compute acceleration at current time
	acceleration = force(sys)

	# compute positions at t + dt
	positions(sys.state) .+= 
		velocities(sys.state) .* sys.dt .+ 
		0.5 .* acceleration .* (sys.dt)^2

	# enforce boundary conditions
	# (basically check if any particles left the box and put them back)
	# see function implementation for deets
	keep_particles_in_box!(sys)

	# compute velocities at t + dt
	velocities(sys.state) .+= 
		0.5 * sys.dt .* (acceleration + force(sys))
end

function evolve!(sys::ParticleSystem, runtime::Float64=10.0)
	numsteps = Int64(abs(runtime/sys.dt) + 1)

	print_evolution_message(runtime, numsteps)

	@time for step in 1:numsteps
		verlet_step!(sys)
		zero_total_momentum!(sys)

		if (step % sys.sampleInterval == 1)
			push!(sys.timeData, sys.t)
			push!(sys.energyData, energy(sys))
			push!(sys.xData, positions(sys.state))
			push!(sys.vData, velocities(sys.state))

			T = temperature(sys)
			push!(sys.tempData, T)
			sys.tempAccumulator += T
			sys.squareTempAccumulator += T^2
		end

		sys.t += sys.dt
		sys.steps += 1
	end
	println("done.") 
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

function reset_statistics!(sys::ParticleSystem)
	sys.steps = 0
	sys.tempAccumulator = 0.0
	sys.squareTempAccumulator = 0.0
	sys.virialAccumulator = 0.0
	sys.xData = []
	sys.vData = []
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
	return mean(sys.energyData)
end

function std_energy(sys::ParticleSystem)
	return std(sys.energyData)
end

