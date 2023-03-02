# molecular dynamics of a gas of atoms
# units where m = ε = σ = k_B = 1
# bog-standard reimplementation of old python thing

dt = 0.001
sampleInterval = 100

mutable struct System
	N::Int64		# number of particles
	L::Float64		# length of square side
	initTemp::Float64
	t::Float64			# initial time
	tPoints::Vector{Float64}	# array of time steps added to during integration
	steps::Int64
	energyPoints::Vector{Float64}	# list of energy, sampled every sampleInterval steps
	sampleTimePoints::Vector{Float64}
	tempPoints::Vector{Float64}
	tempAccum::Float64
	sqTempAccum::Float64
	virialAccum::Float64
	x::Vector{Float64}		# array of N (x,y) positions
	v::Vector{Float64}		# array of N (vx,vy) velocities
	xPoints::Vector{Float64}
	vPoints::Vector{Float64}
	forceType::String
end

function System(N::Int64, L::Float64, initTemp::Float64)
	N = N
	L = L
	initTemp = initTemp
	t = 0.0
	tPoints = [ t ]
	steps = 0
	energyPoints = []
	sampleTimePoints = []
	tempPoints = [ initTemp ]
	tempAccum = 0.0
	sqTempAccum = 0.0
	virialAccum = 0.0
	x = zeros(2 * N)
	v = zeros(2 * N)
	xPoints = []
	vPoints = []
	forceType = "lennardjones"

	return System(
		N, L, initTemp, t, tPoints, steps, energyPoints, 
		sampleTimePoints, tempPoints, tempAccum, sqTempAccum, 
		virialAccum, x, v, xPoints, vPoints, forceType
	)
end

function minimumImage(s::System, x::Float64)
end

function force(s::System)
end

function lennardJonesForce(s::System)
end

function powerLawForce(s::System)
end

# time evol methods

function verletStep(s::System)
end

function evolve(s::System, time::Float64=10.0)
end

function zeroTotalMomentum(s::System)
end

function reverseTime(s::System)
end

function cool(s::System, time::Float64=1.0)
end

# initial condition methods

function randomPositions(s::System)
end

function triangularLatticePositions(s::System)
end

function rectangularLatticePositions(s::System)
end

function randomVelocities(s::System)
end

# measurement methods

function kineticEnergy(s::System)
end

function potentialEnergy(s::System)
end

function lennardJonesPE(s::System)
end

function energy(s::System)
end

function temperature(s::System)
end

# statistics

function resetStatistics(s::System)
end

function meanTemperature(s::System)
end

function meanSqTemperature(s::System)
end

function meanPressure(s::System)
end

function heatCapacity(s::System)
end

function meanEnergy(s::System)
end

function stdEnergy(s::System)
end

# results

function results(s::System)
end
