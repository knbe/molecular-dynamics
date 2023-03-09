# molecular dynamics 2d.

push!(LOAD_PATH, "./")
using Statistics
using ParticleSystems
using Vectors

@static if @isdefined(graphs)
	#using Graphs
	using Plots
end

@static if @isdefined(render)
	#using Render
	print_GLMakie_message()
	@time using GLMakie
end

mutable struct Data
	steps::Int64
	tPoints::Vector{Float64}	# array of time steps added to during integration
	sampleTimePoints::Vector{Float64}
	energyPoints::Vector{Float64}	# array of sampled energy points
	tempPoints::Vector{Float64}
	tempAccumulator::Float64
	squareTempAccumulator::Float64
	xPoints
	vPoints
end

function Data(sys::ParticleSystem)
	steps = 0
	tPoints = [ 0.0 ]
	energyPoints = []
	sampleTimePoints = []
	tempPoints = [ sys.T₀ ]
	tempAccumulator = 0.0
	squareTempAccumulator = 0.0
	xPoints = []
	vPoints = []

	return Data(steps, tPoints, sampleTimePoints, 
		energyPoints, tempPoints, tempAccumulator, 
		squareTempAccumulator, xPoints, vPoints)
end

global DT::Float64 = 0.001
const global SAMPLEINTERVAL::Int64 = 100
const global D::Int64 = 2

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

function minimum_image(x::Vector{Float64}, L::Float64)
	halfL = 0.5 * L
	for i in 1:length(x)
		if x[i] > halfL
			x[i] -= L
		elseif x[i] < -halfL
			x[i] += L
		end
	end
	return x
end

function force(sys::ParticleSystem)
	force = lennard_jones_force(sys)
	return force
end

function lennard_jones_force(sys::ParticleSystem)
	N = sys.N
	L = sys.L
	tiny = 1.0e-40
	virial = 0.0
	force = zeros(2N)

	x = component(1, positions(sys.state))
	y = component(2, positions(sys.state))

	Threads.@threads for i in 1:N
		for j in (i+1):N
			dx = minimum_image(x .- x[i], L)
			dy = minimum_image(y .- y[i], L)

			r2inv = 1.0 ./ (dx.^2 .+ dy.^2 .+ tiny)
			f = 48.0 .* r2inv.^6 - 24.0 .* r2inv.^4
			fx = dx .* f
			fy = dy .* f

			fx[i] = 0.0	# self force is zero
			fy[i] = 0.0

			#virial += dot(fx,dx) + dot(fy,dy)
		end
	end

	return force
end


function lennard_jones_potential(sys::ParticleSystem)

	N = sys.N
	L = sys.L
	tiny = 1.0e-40

	x = component(1, positions(sys.state))
	y = component(2, positions(sys.state))
	U = 0.0

	for i in 1:N
		dx = minimum_image(x[i] .- x, L)
		dy = minimum_image(y[i] .- y, L)

		r2inv = 1.0 / (dx.^2 .+ dy.^2 .+ tiny)
		dU = r2inv.^6 .- r2inv.^3
		dU[i] = 0.0	# self interaction is zero
		U += sum(dU)
	end

	return 2.0 * U
end

function lennard_jones_potential(sys::ParticleSystem)
	U = 0.0
	for i in 1:sys.N
		for j in (i+1):sys.N
			rij_min = minimum_image_vector( 
				particle(j, positions(sys.state)) .- 
				particle(i, positions(sys.state)), 
				sys.L
			)

			r2 = mag_sq(rij_min)

			dU = r2.^(1/6) .- r2.^(1/4)
			U += dU
		end
	end

	return 2.0 * U
end

# time evolution
################################################################################

function verlet_step!(sys::ParticleSystem)
	# compute acceleration
	acceleration = force(sys)

	# verlet update position
	positions(sys.state) .+= velocities(sys.state) .* DT .+ 
				 0.5 * (DT)^2 .* acceleration

	# adjust positions so particles remain in the box (periodic BCs)
	keep_particles_in_box!(sys)

	# verlet update velocity
	velocities(sys.state) .+= 0.5 * DT .* (acceleration .+ force(sys))
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

function evolve!(sys::ParticleSystem, data::Data, time::Float64=10.0)
	print_evolve_message(sys, time)
	steps = Int64(abs(time/DT)+1)

	@static if @isdefined(render)
		#start_render(sys)
		ptx = Observable( sys.state[1:2:2*sys.N] )
		pty = Observable( sys.state[2:2:2*sys.N] )

		fontsize_theme = Theme(fontsize=40)
		set_theme!(fontsize_theme)
		fig = Figure(resolution=(600,400))
		ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y")

		@static if @isdefined(demo2)
			Label(fig[2,1],
				"DT = $DT",
				justification = :left, 
				lineheight = 1.1,
				tellwidth = false,
			)	
		end

		for n = 1:sys.N
			GLMakie.scatter!(
				ax1, 
				ptx, 
				pty, 
				markersize=40.0
			)
		end
		limits!(ax1, 0, sys.L, 0, sys.L)
		display(fig)
		sleep(2)
	end

	for step in 1:steps
		verlet_step!(sys)
		zero_total_momentum!(velocities(sys.state))

		sys.t += DT

		@static if @isdefined(render)
			#update_render(sys)
			ptx[] = sys.state[1:2:2*sys.N]
			pty[] = sys.state[2:2:2*sys.N]

			@static if @isdefined(demo2)
				sleep(0.0001)
			end

			yield()
		end

		@static if @isdefined(datacoll)
			if (step % SAMPLEINTERVAL == 1)
				push!(data.energyPoints, energy(sys))
				push!(data.sampleTimePoints, sys.t)
				push!(data.xPoints, positions(sys.state))
				push!(data.vPoints, velocities(sys.state))
			end
		end

		T = temperature(sys)
		data.steps += 1
		push!(data.tempPoints, T)
		data.tempAccumulator += T
		data.squareTempAccumulator += T^2
	end
end

function reverse_time()
	global DT *= -1

	println("time reversed! DT = $DT")
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

function kinetic_energy(sys::ParticleSystem)
	return 0.5 * sum(velocities(sys.state).^2)
end

function potential_energy(sys::ParticleSystem)
	return lennard_jones_potential(sys)
end

#function temperature(state::Vector{Float64}, N::Int64)
function temperature(sys::ParticleSystem)
	return kinetic_energy(sys) / sys.N
end

function energy(sys::ParticleSystem)
	return potential_energy(sys) + kinetic_energy(sys)
end

# console print/debug stuff
################################################################################

function print_demo_message(demo::String)
	println("\nahoy there!\nbienvenue à molecular dynamics!")
	println("\ndemo: $demo")
end

function print_system_params(sys::ParticleSystem)
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
	println("\tenergy: \t\t$(energy(sys))")
	println("\ttemperature: \t\t$(temperature(sys))")
end

function print_measurements(sys::ParticleSystem)
	println("done.")
	println("\nmeasurements:")
	println("\tenergy: \t\t$(energy(sys))")
	println("\ttemperature: \t\t$(temperature(sys))")
end

function print_evolve_message(sys::ParticleSystem, time::Float64)
	println("\nevolving...")
end

function print_GLMakie_message()
	println("\nloading GLMakie... go make coffee or something\n")
end

# plots
################################################################################

function initialize_plot()
	#theme(:lime)
	plot(size=(800,800), titlefontsize=28)
end

function plot_positions(sys::ParticleSystem)
	N = sys.N

	initialize_plot()
	scatter!(
		xcomponent(positions(sys.state)), 
		ycomponent(positions(sys.state)), 
		markersize=5.0,
		)
	xlabel!("x")
	ylabel!("y")
	xlims!(0, sys.L)
	ylims!(0, sys.L)
end

function plot_trajectories(sys::ParticleSystem, number::Int64=1)
end

function plot_temperature(sys::ParticleSystem)
	plot()
	plot!(sys.tPoints, sys.tempPoints)
	xlabel!("t")
	ylabel!("T")
end

function plot_energy(sys::ParticleSystem)
	plot()
	plot!(sys.sampleTimePoints, sys.energyPoints)
	xlabel!("t")
	ylabel!("energy")
	ylims!(minimum(sys.energyPoints) - 2, maximum(sys.energyPoints) + 2)
end

function velocity_histogram(sys::ParticleSystem)
	plot()
	histogram!(sys.vPoints, normalize=:pdf, bins = -2.0:0.2:2.0)
	xlabel!("velocity")
	ylabel!("probability")
end


# demos
################################################################################


function demo_0()
	print_demo_message("dummy run")
	sys = ParticleSystem(64, 8.0, 1.0, D)
	print_system_params(sys)
	print_initialize()
	set_rectangular_lattice_positions!(sys)
	set_random_velocities!(sys)

	print_initial_measurements(sys)

	data = Data(sys)
	evolve!(sys, data, 10.0)

	print_measurements(sys)
end

function demo_1()
	sys = ParticleSystem(64, 8.0, 1.0)
	set_rectangular_lattice_positions!(sys)
	set_random_velocities!(sys)

	@time begin

	data = Data(sys)
	evolve!(sys, data, 12.0)
	reverse_time()
	evolve!(sys, data, 12.0)

	end

	#plot_positions(sys)
end

function demo_2()
	print_demo_message("time reversal test")

	sys = ParticleSystem(16, 4.0, 1.0, D)

	print_system_params(sys)
	print_initialize()

	set_rectangular_lattice_positions!(sys)
	set_random_velocities!(sys)

	print_initial_measurements(sys)

	data = Data(sys)
	evolve!(sys, data, 5.0)
	reverse_time()
	evolve!(sys, data, 5.0)
end

function demo_3()
	print_demo_message("conservation of momentum")

	sys = ParticleSystem(64, 8.0, 1.0, D)

	print_system_params(sys)
	print_initialize()

	set_rectangular_lattice_positions!(sys)
	set_random_velocities!(sys)

	print_initial_measurements(sys)

	data = Data(sys)
	evolve!(sys, data, 5.0)
	print_measurements(sys)
end
