push!(LOAD_PATH, "./")

using Statistics
using Vectors
using Plots

#mutable struct State
#	pos
#	vel
#end

mutable struct ParticleSystem
	N::Int64
	D::Int64
	L::Float64
	T₀::Float64
	t::Float64
	#state::State
	state::Matrix{Float64}
end

function ParticleSystem(N::Int64, D::Int64, L::Float64, T₀::Float64)
	state = zeros(N, 2D)
	t = 0.0
	
	return ParticleSystem(N, D, L, T₀, t, state)
end

#@views position(state) = state[1:Int64((length(state)/2))]
#@views velocity(state) = state[Int64((length(state)/2)+1):end]

@views position(state) = state[:,1:(Int64((size(state)[2])/2))]
@views velocity(state) = state[:,(Int64((size(state)[2])/2)+1):end]

@views particle(n, matrix) = matrix[n,:]

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
end

function zero_total_momentum!(state::Matrix{Float64})
	vel = velocity(state)
	D = size(vel)[2]
	
	for d in 1:D
		vel[:,d] .-= mean(vel[:,d])
	end

	velocity(state) = vel
end

function minimum_separation_vector(p1::Vector{Float64}, p2::Vector{Float64}, L::Float64)
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

				r2inv = 1.0 / (mag_sq(rij) + tiny)
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

function verlet_step!(sys::ParticleSystem, dt::Float64)
	accel = lennard_jones_force(sys)
	for i in 1:sys.N
		sys.state.pos[i] .+= (sys.state.vel[i] .* dt) .+ (0.5 * dt^2 .* accel[i])
		sys.state.pos[i] = sys.state.pos[i] .% sys.L
		sys.state.vel[i] .+= 0.5 * dt .* (accel[i] + lennard_jones_force(sys)[i])
	end
end

function evolve!(sys::ParticleSystem, dt::Float64, tt::Float64=10.0)
	numSteps = Int64(floor(tt/dt))

	statedata = [ sys.state ]
	tdata = 0:dt:tt

	for t in 1:numSteps
		verlet_step!(sys, dt)
		#zero_total_momentum!(s)

		#s.t += dt
		#push!(statedata, sys.state)
	end

	#xdata = [ statedata ]
	#plot(sy
	println(statedata[1][1])
end

function plot_positions(sys::ParticleSystem)
	N = sys.N
	#scatter(s.x[1:N], s.x[(N+1):2N])
end


sys = ParticleSystem(16, 2, 5.0, 1.0)

initialize_positions_rectangular!(sys)
initialize_velocities!(sys)
zero_total_momentum!(sys.state)
#evolve!(sys, 0.01, 10.0)

#accel = lennard_jones_force(sys)
