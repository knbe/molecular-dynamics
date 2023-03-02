
push!(LOAD_PATH, "./")
using VectorMath
using Statistics
using Plots

# this struct needs to be split up imo
mutable struct ParticleSystem
	N::Int64		# number of particles
	L::Float64		# side length
	D::Int64		# dimensions
	T₀::Float64		# initial temperature T₀
	T_accumulator::Float64
	T_sq_accumulator::Float64
	T_data::Vector{Float64}
	E_data::Vector{Float64}
	virial_accumulator::Float64
	t::Float64		# system time
	t_data::Vector{Float64}
	pos::Vector{Float64}	# position array (x1, ..., xN, y1, ..., yN)
	vel::Vector{Float64}	# velocity array (vx1, ..., vxN, vy1, ..., vyN)
	pos_data
	vel_data
end

function ParticleSystem(N::Int64=64, L::Float64=10.0, T₀::Float64=1.0)
	N = N
	L = L
	D = 2			# later generalise this to 3D
	T₀ = T₀
	T_accumulator = 0.0
	T_sq_accumulator = 0.0
	T_data = []
	E_data = []
	virial_accumulator = 0.0
	t = 0.0
	t_data = [ 0.0 ]
	pos = zeros(D * N)
	vel = zeros(D * N)
	pos_data = [ pos ]
	vel_data = [ vel ]
	return ParticleSystem(
		N, L, D, T₀, T_accumulator, T_sq_accumulator, T_data, 
		E_data, virial_accumulator, t, t_data, pos, vel, pos_data, 
		vel_data
		)
end

function minimum_image(x::Vector{Float64}, L::Float64)
	return (x .+ (L/2)) .% L .- (L/2)
end

function force(L::Float64, pos::Vector{Float64})
	# implement something here like an ifdef so can use different forces
	#ifdef lennard_jones
	force, virial = lennard_jones_force(L, pos)
	#endif

	return force, virial
end
	
function lennard_jones_force(L::Float64, pos::Vector{Float64})
	N = Int64(length(pos) / 2)
	x = pos[1:N]
	y = pos[(N+1):2N]

	force = zeros(2*N)
	virial = 0.0
	min = 1.0e-40

	for i in 1:N
		for j in 2:N
			dx = x[i] - x[j]
			if dx > (L/2)
				dx -= L
			elseif dx < -(L/2)
				dx += L
			end

			dy = y[i] -y[j]
			if dy > (L/2)
				dy -= L
			elseif dy < -(L/2)
				dy += L
			end

			r2inv = 1.0 / (dx^2 + dy^2 + min)
			c = 48.0 * r2inv^6 - 24.0 * r2inv^4
			fx = c * dx
			fy = c * dy

			force[i] += fx
			force[i+N] += fy
			force[j] -= fx
			force[j+N] -= fy

			virial += fx * dx + fy * dy
		end
	end

	return force, virial
end

function zero_total_momentum!(vel::Vector{Float64})
	N = Int64(length(vel) / 2)
	vx = vel[1:N]
	vy = vel[(N+1):2N]

	vx .-= mean(vx)
	vy .-= mean(vy)

	vel[1:N] = vx
	vel[(N+1):2N] = vy
	return vel
end

function verlet_step!(pos::Vector{Float64}, vel::Vector{Float64}, L::Float64, dt::Float64)
	acc = force(L, pos)
	println(acc)
	#pos .+= vel .* dt .+ (0.5 * dt^2) .* acc
#	pos = pos .% L	# ehh, I want to get this step out of here. remove dependency on L.
#	vel .+= 0.5 * dt .* (acc + force(L, pos))
	return pos, vel
end


function kinetic_energy(vel::Vector{Float64})
	return 0.5 * mag_sq(vel)
end


function temperature(vel::Vector{Float64})
	N = Int64(length(vel) / 2)
	return kinetic_energy(vel) / N
end

function rectangular_lattice!(sys::ParticleSystem)
	N = sys.N
	L = sys.L
	nx = Int64(sqrt(N))
	ny = nx
	dx = L / nx
	dy = L / ny

	n = 1
	for i in 1:nx
		for j in 1:ny
			sys.pos[n] = (i - 0.5) * dx
			sys.pos[n+N] = (j - 0.5) * dy
			n += 1
		end
	end
end

function random_velocities!(sys::ParticleSystem)
	sys.vel = rand(2 * sys.N) .- 0.5
	zero_total_momentum!(sys.vel)
	T = temperature(sys.vel)
	rescale_factor = sqrt(sys.T₀ / T)
	sys.vel .*= rescale_factor
end

function evolve!(sys::ParticleSystem, time::Float64=10.0, 
		dt::Float64=0.01, samplingInterval::Int64=100)
	numSteps = Int64(time/dt)

	for t in 1:1
		sys.pos, sys.vel = verlet_step!(sys.pos, sys.vel, sys.L, dt)
#		sys.vel = zero_total_momentum!(sys.vel)
#
#		sys.t += dt
#		sys.steps += 1
#
#		if (t % samplingInterval == 0)
#			push!(sys.t_data, sys.t)
#			#push!(s.energyPoints, energy(s))
#			push!(sys.pos_data, sys.pos)
#			push!(sys.vel_data, sys.vel)
#		end
#
#		T = temperature(sys.vel)
#		sys.T_accumulator += T
#		sys.T_sq_accumulator += T^2
#		push!(sys.T_data, T)
	end
end

function print_console(sys::ParticleSystem)
	println(
		  "  time:			", sys.t, 
		#"\n  total energy:		", energy(s), 
		"\n  temperature:		", temperature(sys.vel),
		)
end

function plot_positions(pos::Vector{Float64})
	N = Int64(length(pos) / 2)
	scatter(pos[1:N], pos[(N+1):2N])
end

sys = ParticleSystem(16, 10.0, 1.0)

rectangular_lattice!(sys)
#plot_positions(sys.pos)
random_velocities!(sys)
print_console(sys)

evolve!(sys, 10.0, 0.01, 100)
