# molecular dynamics of a gas of atoms in 2D
# this is a barebones reimplementation of the old python script
# units: m = ε = σ = k_B = 1

using Plots
using Statistics

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
	xPoints
	vPoints
	forceType::String
end

function System(N::Int64=16, L::Float64=10.0, initTemp::Float64=1.0)
	N = N
	L = L
	initTemp = initTemp
	t = 0.0
	tPoints = [ 0.0 ]
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

# VECTORS/MATH

function dotProduct(v1::Vector{Float64}, v2::Vector{Float64})
	dotProduct = 0
	for i in 1:length(v)
		dotProduct += v1[i] * v2[i]
	end
	return dotProduct
end

function magnitudeSquared(v::Vector{Float64})
	magSq = 0
	for i in 1:length(v)
		magSq += v[i] * v[i]
	end
	return magSq
end

# FORCES

function minimumImage(s::System, x)
	L = s.L
	halfL = 0.5 * L

	return (x .+ halfL) .% L .- halfL
end

function force(s::System)
	f, virial = lennardJonesForce(s)
	s.virialAccum += virial
	return f
end

function lennardJonesForce(s::System)

	N = s.N
	L = s.L
	halfL = 0.5 * L
	virial = 0.0
	tiny = 1.0e-40

	x = s.x[1:N]
	y = s.x[(N+1):2N]
	f = zeros(2N)


	minImg = minimumImage(s, x)

	#println(minImg)
	for i in 1:N
		for j in 2:N
			dx = x[i] - x[j]
			if dx > halfL
				dx -= L
			elseif dx < -halfL
				dx += L
			end

			dy = y[i] -y[j]
			if dy > halfL
				dy -= L
			elseif dy < -halfL
				dy += L
			end

			r2inv = 1.0 / (dx^2 + dy^2 + tiny)
			r6inv = r2inv^3
			r8inv = r2inv * r6inv
			#c = 48.0 * r8inv * r6inv - 24.0 * r8inv
			c = 48.0 * r2inv^7 - 24.0 * r2inv^4
			fx = dx * c
			fy = dy * c

			f[i] += fx
			f[i+N] += fy
			f[j] -= fx
			f[j+N] -= fy

			virial += fx * dx + fy * dy
		end
	end

	return f, virial

end

function powerLawForce(s::System)
end

# TIME EVOLUTION

function verletStep!(s::System)
	a = force(s)
	s.x .+= s.v .* dt .+ 0.5 * dt^2 .* a
	s.x = s.x .% s.L
	s.v .+= 0.5 * dt .* (a + force(s))
end

function evolve!(s::System, time::Float64=10.0)
	numSteps = Int64(time/dt)

	for t in 1:numSteps
		verletStep!(s)
		zeroTotalMomentum!(s)

		s.t += dt
		push!(s.tPoints, s.t)

		if (t % sampleInterval == 0)
			push!(s.sampleTimePoints, s.t)
			#push!(s.energyPoints, energy(s))
			push!(s.xPoints, s.x)
			push!(s.vPoints, s.v)
		end

		T = temperature(s)
		s.steps += 1
		push!(s.tempPoints, T)
		s.tempAccum += T
		s.sqTempAccum += T^2
	end
end

function reverseTime(s::System)
end

function cool(s::System, time::Float64=1.0)
end

# INITIALISATION

function randomPositions(s::System)
end

function triangularLatticePositions(s::System)
end

function rectangularLatticePositions!(s::System)
	nx = Int64(sqrt(s.N))
	ny = nx
	dx = s.L / nx
	dy = s.L / ny

	n = 1
	for i in 1:nx
		for j in 1:ny
			s.x[n] = (i - 0.5) * dx
			s.x[n+s.N] = (j - 0.5) * dy
			n += 1
		end
	end
end

function zeroTotalMomentum!(s::System)
	vx = s.v[1:s.N]
	vy = s.v[(s.N + 1):(2 * s.N)]

	vx .-= mean(vx)
	vy .-= mean(vy)

	s.v[1:s.N] = vx
	s.v[(s.N + 1):(2 * s.N)] = vy

end

function randomVelocities!(s::System)
	s.v = rand(2 * s.N) .- 0.5
	zeroTotalMomentum!(s)
	T = temperature(s)
	rescale = sqrt(s.initTemp / T)
	s.v .*= rescale
end

# MEASUREMENTS

function kineticEnergy(s::System)
	return 0.5 * magnitudeSquared(s.v)
end

function potentialEnergy(s::System)
end

function lennardJonesPE(s::System)
end

function energy(s::System)
end

function temperature(s::System)
	return kineticEnergy(s) / s.N
end

# STATISTICS

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

# PRINT RESULTS

function results(s::System)
	println(
		  "  time:			", s.t, 
		#"\n  total energy:		", energy(s), 
		"\n  temperature:		", temperature(s),
		)
end

# time reversal test

#rectangularLatticePositions
#randomVelocities
#results
#evolve(time=20.0)
#reverseTime
#evolve
#results
#plotPositions
#showPlots

# equilibration and statistics

#triangularLatticePositions
#randomVelocities
#plotPositions
#results
#
#evolve(time=10.0) # initial time evol
#resetStatistics	# remove transient behaviour
#evolve(time=20.0) # accumulate statistics
#
#results
#plotEnergy
#plotTrajectories
#plotTemperature
#velocityHistogram
#showPlots

s = System(64, 10.0, 1.0)

rectangularLatticePositions!(s)
randomVelocities!(s)
results(s)
scatter(s.x[1:s.N], s.x[(s.N+1):(2*s.N)])

evolve!(s, 10.0)
results(s)
scatter(s.x[1:s.N], s.x[(s.N+1):(2*s.N)])
