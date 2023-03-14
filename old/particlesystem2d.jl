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

