push!(LOAD_PATH, "./")

using Statistics
using Vectors
using Plots

import LinearAlgebra: norm

mutable struct ParticleSystem
	N::Int64
	D::Int64
	L::Float64
	T₀::Float64
	t::Float64
	state::Matrix{Float64}
end

function ParticleSystem(N::Int64, D::Int64, L::Float64, T₀::Float64)
	state = zeros(N, 2D)
	t = 0.0
	
	return ParticleSystem(N, D, L, T₀, t, state)
end

@views position(state) = state[:,1:(Int64((size(state)[2])/2))]
@views velocity(state) = state[:,(Int64((size(state)[2])/2)+1):end]
@views particle(n, matrix) = matrix[n,:]
@views coordinate(i, matrix) = matrix[:,i]

function initialize_positions_rectangular!(sys::ParticleSystem)
	NL = Int64(floor((sys.N)^(1/sys.D)))
	if sys.N != NL^(sys.D)
		sys.N = NL^(sys.D)
	end

	a = sys.L / NL	# lattice spacing

	n = 1
	for i in 1:NL
		for j in 1:NL
			particle(n, position(sys.state))[1] = (i - 0.5) * a
			particle(n, position(sys.state))[2] = (j - 0.5) * a
			n += 1
		end
	end
end

function initialize_velocities!(sys::ParticleSystem)
	for n in 1:length(velocity(sys.state))
		velocity(sys.state)[n] = rand() - 0.5
	end

	T = temperature(sys)
	println(T)
	rescale = sqrt(sys.T₀ / T)
	velocity(sys.state) .*= rescale

end

function zero_total_momentum!(velocity)
	numCoordinates = size(velocity)[2]
	
	for i in 1:numCoordinates
		coordinate(i, velocity) .-= mean(coordinate(i, velocity))
	end
end

function kinetic_energy(velocity)
	N,D = size(velocity)
	total = 0.0
	for d in 1:D
		total += mag_sq(velocity[:,d])
	end
	return total
end

function temperature(sys::ParticleSystem)
	return kinetic_energy(velocity(sys.state)) / sys.N
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

function lennard_jones_force(sys::ParticleSystem)
	tiny = 1.0e-40

	accel = [ zeros(sys.d) for i in 1:sys.N ]

	n = 1
	for i in 1:sys.N
		accel_i = zeros(sys.d)
		for j in 1:sys.N
			if i != j
				rij = minimum_separation_vector(
					sys.state.pos[i], sys.state.pos[j], sys.L)

				r2inv = 1.0 / ((norm(rij))^2 + tiny)
				rhat = unit(rij)

				a = (1.0/mag(rij)) * 48.0 * r2inv^6 - 24.0 * r2inv^3
				accel_ij = rhat .* a
				accel_i += accel_ij

#				println("i = ", i, ", j = ", j)
#				println("rij = ", rij)
#				println("r2inv = ", r2inv)
#				println("rhat = ", rhat)
#				println("a = ", a)
#				println("accel_ij = ", accel_ij)
#				println("accel_i = ", accel_i)
#				println()
			end
		end
		n += 1
		accel[i] = accel_i
	end

	return accel
end

function ljforce(position, L::Float64)
	N,D = size(position)
	acceleration = zeros(N,D)
	tiny = 1.0e-40
	#println(accel)

	for i in 1:N
		accel_i = zeros(D)
		for j in 1:N
			if i != j
				sep = minimum_separation_vector(
					particle(i, position), particle(j, position), L
					)

	#			#println(typeof(sep)
				sepSquared = mag_sq(sep) + tiny

				magnitude = (1.0/mag(sep)) * 
					(48.0 * sepSquared^(-6) - 24.0 * sepSquared^(-3))
				direction = unit(sep)
				accel_ij = magnitude .* direction
				accel_i += accel_ij
	
#				println("i = ", i, "j = ", j)
#				println(sep, ", ", sepSquared)
#				println(magnitude)
#				println(direction)
#				println(accel_ij)
#				println(accel_i)
			end
		end
		particle(i, acceleration)[1] = accel_i[1]
		particle(i, acceleration)[2] = accel_i[2]
	end

	return acceleration
end

function verlet_step!(state::Matrix{Float64}, dt::Float64, L::Float64)
	acceleration = ljforce(position(state), L)

	position(state) .+= velocity(state) .* dt .+ 0.5 * dt^2 .* acceleration
	position(state) .% L
	velocity(state) .+= 0.5 * dt .* (acceleration + ljforce(position(state), L))
end

function evolve!()
end

function run!(sys::ParticleSystem, dt::Float64, tt::Float64=10.0)
	numSteps = Int64(floor(tt/dt))

#	statedata = [ sys.state ]
#	tdata = 0:dt:tt
#
	for t in 1:150
		verlet_step!(sys.state, dt, sys.L)
		zero_total_momentum!(velocity(sys.state))

		#s.t += dt
		#push!(statedata, sys.state)
	end

	#xdata = [ statedata ]
	#plot(sy
	#println(statedata[1][1])
end

function plot_positions(sys::ParticleSystem)
	N = sys.N
	xdata = coordinate(1, position(sys.state))
	ydata = coordinate(2, position(sys.state))
	plot(size=(800,800))
	scatter!(xdata, ydata)
end


sys = ParticleSystem(16, 2, 5.0, 1.0)

initialize_positions_rectangular!(sys)
initialize_velocities!(sys)
zero_total_momentum!(velocity(sys.state))
run!(sys, 0.001, 0.2)
plot_positions(sys)

#ljforce(position(sys.state), 5.0)
#verlet_step!(sys.state, 0.01, 5.0)


#accel = lennard_jones_force(sys)
