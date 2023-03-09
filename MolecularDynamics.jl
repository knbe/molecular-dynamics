# molecular dynamics 2d.

push!(LOAD_PATH, "./")
using Statistics
using Vectors
using Plots

@static if @isdefined(graphs)
	using Graphs
end

@static if @isdefined(render)
	using Render
end

global DT::Float64 = 0.001
const global SAMPLEINTERVAL::Int64 = 100
const global D::Int64 = 2

mutable struct ParticleSystem
	N::Int64		# number of particles
	L::Float64		# side length of square box
	T₀::Float64	
	t::Float64		# time
	state::Vector{Float64}	# state space vector
end

function ParticleSystem(N, L, T₀)
	t = 0.0
	state = zeros(D*N*2) # [x1,y1,z1,x2,y2,z2,...,vx1,vy1,vz1,...]
	return ParticleSystem(N, L, T₀, t, state)
end

# some useful "views" of the state array
@views positions(state) = state[1:Int64(length(state)/2)]
@views velocities(state) = state[(Int64(length(state)/2)+1):end]
@views component(i, vector) = vector[i:D:end]
@views particle(n, vector) = [ vector[D*n-1], vector[D*n] ]
#@views particle(n, vector) = [ vector[i] for i in (D*n-(D-1)):(D*n) ]
@views numParticles(state) = Int64(length(state) / (2*D))

# initialization
################################################################################

function set_rectangular_lattice_positions!(sys::ParticleSystem)
	n = Int64(floor(sqrt(sys.N)))
	sys.N = n^2
	a = sys.L / n

	for i in 0:(n-1)
		for j in 0:(n-1)
			positions(sys.state)[2*(i*n+j)+1] = (i + 0.5) * a
			positions(sys.state)[2*(i*n+j)+2] = (j + 0.5) * a
		end
	end

	println("\tposition config: \trectangular lattice")
end

function set_random_positions!(sys::ParticleSystem)
	positions(sys.state) .= rand(2*sys.N) .* sys.L
	cool(sys)

	println("\tposition config: \trandom")
end

function set_random_velocities!(sys::ParticleSystem)
	velocities(sys.state) .= rand(D*sys.N) .- 0.5
	zero_total_momentum!(velocities(sys.state))
	velocities(sys.state) .*= sqrt(sys.T₀/temperature(sys))
	
	println("\tvelocity config: \trandom")
end

function zero_total_momentum!(velocities)
	for i = 1:D
		component(i, velocities) .-= mean(component(i, velocities))
	end
end

# forces
################################################################################

function minimum_image_vector(r12::Vector{Float64}, L::Float64)
	for i in 1:length(r12)
		if r12[i] > (L/2)
			r12[i] -= L
		elseif r12[i] < -(L/2)
			r12[i] += L
		end
	end
	return r12
end

#function force(state::Vector{Float64}, L::Float64)
#	force = lennard_jones_force(state, L)
#	return force
#end

function ljforce(sys::ParticleSystem)
	force = zeros(D * sys.N)
	for i in 1:sys.N
		Threads.@threads for j in (i+1):sys.N
			rij_min = minimum_image_vector( 
				particle(j, positions(sys.state)) .- 
				particle(i, positions(sys.state)), 
				sys.L
			)

			r2 = mag_sq(rij_min)

			f = (48.0 * r2^(1/6) - 24.0 * r2^(1/4)) .* unit(rij_min)

			force[(D*i-(D-1)):D*i] .+= f
			force[(D*j-(D-1)):D*j] .-= f

		end
	end

	return force
end

# time evolution
################################################################################

function verlet_step!(sys::ParticleSystem)
	# compute acceleration
	acceleration = ljforce(sys)

	# verlet update position
	positions(sys.state) .+= velocities(sys.state) .* DT .+ 
				 0.5 * (DT)^2 .* acceleration

	# adjust positions so particles remain in the box (periodic BCs)
	keep_particles_in_box!(sys)

	# verlet update velocity
	velocities(sys.state) .+= 0.5 * DT .* (acceleration .+ ljforce(sys))
end

function keep_particles_in_box!(sys::ParticleSystem)
	for i in 1:(D*sys.N)
		if positions(sys.state)[i] > sys.L
			positions(sys.state)[i] -= sys.L
		elseif positions(sys.state)[i] < 0.0
			positions(sys.state)[i] += sys.L
		end
	end
end

function evolve!(sys::ParticleSystem, time::Float64=10.0)
	print_evolve_message(sys, time)
	steps = Int64(abs(time/DT)+1)

	@static if @isdefined(render)
		start_render(sys)
	end

	for step in 1:steps
		verlet_step!(sys)
		zero_total_momentum!(velocities(sys.state))

		sys.t += DT

		@static if @isdefined(render)
			update_render(sys)
		end

		if (step % SAMPLEINTERVAL == 1)
#			push!(sys.energyPoints, energy(sys.state, sys.L))
#			push!(sys.sampleTimePoints, sys.t)
#			push!(sys.xPoints, positions(sys.state))
#			push!(sys.vPoints, velocities(sys.state))
		end

#		T = temperature(sys.state)
#		sys.steps += 1
#		push!(sys.tempPoints, T)
#		sys.tempAccumulator += T
#		sys.squareTempAccumulator += T^2
	end
end

function reverse_time()
	global DT *= -1
end

function cool(sys::ParticleSystem, time::Float64=1.0)
	steps = Int64(time/dt)
	for step in 1:steps
		verlet_step!(sys.state, sys.L)
		velocities(sys.state) .*= (1.0 - dt)
	end

	#reset_statistics!(sys)
end

# MEASUREMENTS

#function kinetic_energy(state::Vector{Float64})
function kinetic_energy(sys::ParticleSystem)
	return 0.5 * sum(velocities(sys.state).^2)
end

#function potential_energy(state::Vector{Float64}, L::Float64)
#	return lennard_jones_potential(state, L)
#end

#function temperature(state::Vector{Float64}, N::Int64)
function temperature(sys::ParticleSystem)
	return kinetic_energy(sys) / sys.N
end

#function energy(state::Vector{Float64}, L::Float64)
#	return potential_energy(state, L) + kinetic_energy(state)
#end

# console print/debug stuff
################################################################################

function print_system_params(sys::ParticleSystem)
	println("\nahoy there!\nbienvenue à molecular dynamics!")
	println("\nsystem parameters:")
	println("\tN =  $(sys.N)   (number of particles)")
	println("\tL =  $(sys.L)   (side length of square box)")
	println("\tT₀ = $(sys.T₀)  (initial temperature)")
	println("\tDT = $DT  (time step)")
end

function print_initialize()
	println("\ninitializing...")
end

function print_initial_measurements(sys::ParticleSystem)
	println("done.")
	println("\ninitial measurements:")
	println("\tenergy: \t\tfuck knows")
	println("\ttemperature: \t\tfuck knows")
end

function print_evolve_message(sys::ParticleSystem, time::Float64)
	println("\nevolving...")
end

# demos
################################################################################

function demo1()
	sys = ParticleSystem(64, 8.0, 1.0)
	set_rectangular_lattice_positions!(sys)
	set_random_velocities!(sys)

	@time begin

	evolve!(sys, 12.0)
	reverse_time()
	evolve!(sys, 12.0)

	end

	#plot_positions(sys)
end

function demo2()
	sys = ParticleSystem(16, 4.0, 1.0)
	print_system_params(sys)

	print_initialize()
	set_rectangular_lattice_positions!(sys)
	set_random_velocities!(sys)
	print_initial_measurements(sys)

	evolve!(sys, 12.0)
end
