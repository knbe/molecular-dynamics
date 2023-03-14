function powerLawForce(s::System)
end


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


function lennard_jones_force3(sys::ParticleSystem2D)
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

			dx = x[i] - x[j]
			if dx > (L/2)
				dx -= L
			elseif dx < -(L/2)
				dx += L
			end

			dy = y[i] - y[j]
			if dy > (L/2)
				dy -= L
			elseif dy < -(L/2)
				dy += L
			end

			r2inv = 1.0 / (dx^2 + dy^2)
			r6inv = r2inv^3
			r8inv = r2inv * r6inv
			c = 48.0 * r8inv * r6inv - 24.0 * r8inv
			fx = dx * c
			fy = dy * c

			force[2*i-1] += fx
			force[2*i] += fy
			force[2*j-1] -= fx
			force[2*j] -= fy

			virial += fx*dx + fy*dy
		end
	end

	return force, (0.5 * virial)
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
	return 0.5 * mag_sq(s.v)
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


function plot_positions(s::System)
	N = s.N
	scatter(s.x[1:N], s.x[(N+1):2N])
end

lennardJones = 1

s = System(64, 10.0, 1.0)

rectangularLatticePositions!(s)
randomVelocities!(s)
results(s)

evolve!(s, 10.0)
results(s)
plot_positions(s)
