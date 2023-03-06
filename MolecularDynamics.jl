push!(LOAD_PATH, "./")

using Statistics
using Vectors

mutable struct State
	pos
	vel
end

mutable struct ParticleSystem
	N::Int64
	d::Int64
	L::Float64
	T₀::Float64
	t::Float64
	state::State
end

function ParticleSystem(N::Int64, d::Int64, L::Float64, T₀::Float64)
	N_axis = Int64(floor(sqrt(N)))
	N = N_axis^d
	pos = [ zeros(d) for i in 1:N ]
	vel = [ zeros(d) for i in 1:N ]
	t = 0.0
	
	return ParticleSystem(N, d, L, T₀, t, State(pos, vel))
end

function initialize_positions_rectangular!(sys::ParticleSystem)
	N_axis = (sys.N)^(1/sys.d)
	a = sys.L / N_axis

	n = 1
	for i in 1:N_axis
		for j in 1:N_axis
			sys.state.pos[n][1] = (i - 0.5) * a
			sys.state.pos[n][2] = (j - 0.5) * a
			n += 1
		end
	end
end

function initialize_velocities!(sys::ParticleSystem)
	sys.state.vel = [ (rand(sys.d) .- 0.5) for i in 1:sys.N ]
end

function zero_total_momentum!(sys::ParticleSystem)
#	for coord in 1:sys.d
#		v = [ sys.state.vel[n][coord] for n in sys.N ]
#
#		v .-= mean(v)
#
#		for n in 1:sys.N
#			sys.state.vel[n][coord] = v[n]
#		end
#	end
	
	vx = [ sys.state.vel[n][1] for n in 1:sys.N ]
	vy = [ sys.state.vel[n][2] for n in 1:sys.N ]

	vx .-= mean(vx)
	vy .-= mean(vy)

	for n in 1:sys.N
		sys.state.vel[n][1] = vx[n]
		sys.state.vel[n][2] = vy[n]
	end
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


sys = ParticleSystem(16, 2, 5.0, 1.0)

initialize_positions_rectangular!(sys)
initialize_velocities!(sys)
zero_total_momentum!(sys)
evolve!(sys, 0.01, 10.0)

#accel = lennard_jones_force(sys)
