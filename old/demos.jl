
push!(LOAD_PATH, "./")
includet("particlesystem2d.jl")
includet("molecular.jl")
includet("initialize.jl")
includet("lennardjones.jl")
includet("graphs.jl")

# DEMO 0: APPROACH TO EQUILIBRIUM
function demo_0()
	println("\nDEMO 0: APPROACH TO EQUILIBRIUM")
	println("----------------------------------------") 

	sys = ParticleSystem(64, 120.0, 1.0)
	print_system_parameters(sys)

	set_square_lattice_positions!(sys)
	set_random_velocities!(sys)
	print_system_data(sys)
	p1 = plot_positions(sys)

	evolve!(sys, 1.0)
	print_system_data(sys)

	p2 = plot_trajectories(sys, collect(1:64))
	p3 = plot_energy(sys)
	p4 = plot_temperature(sys)

	plot(
		p1, p2, p3, p4,
		layout = grid(2,2, heights=[0.7,0.3]),
		size = (1280,720)
	)
end

# DEMO 1: TIME REVERSAL TEST
function demo_1()
	println("\nDEMO 1: TIME REVERSAL TEST")
	println("----------------------------------------") 

	sys = ParticleSystem(64, 120.0, 1.0)
	print_system_parameters(sys)

	set_square_lattice_positions!(sys)
	set_random_velocities!(sys)
	print_system_data(sys)
	p1 = plot_positions(sys)

	evolve!(sys, 50.0)
	#p2 = plot_trajectories(sys, collect(1:64))
	p2 = plot_positions(sys)

	reverse_time!(sys)
	evolve!(sys, 50.0)
	print_system_data(sys)
	#p3 = plot_trajectories(sys, collect(1:64))
	p3 = plot_positions(sys)

	plot(
		p1, p2, p3,
		layout = (1,3),
		size = (1200,400)
	)
end

# DEMO 2: SPEED DISTRIBUTION
function demo_2()
	println("\nDEMO 2: SPEED DISTRIBUTION")
	println("----------------------------------------") 

	sys = ParticleSystem[]

	# array for speed distribution plots
	ps = Plots.Plot{Plots.GRBackend}[]

	# array for trajectory plots
	pt = Plots.Plot{Plots.GRBackend}[]

	# initialize three systems with different initial conditions 
	# but same KE and PE, evolve, and save plots
	for i = 1:3
		push!(sys, ParticleSystem(64, 120.0, 1.0))

		println("\nSYSTEM $i")
		print_system_parameters(sys[i])

		set_square_lattice_positions!(sys[i])
		add_position_jitter!(sys[i])
		set_random_velocities!(sys[i])
		print_system_data(sys[i])

		evolve!(sys[i], 40.0)
		print_system_data(sys[i])
		push!(ps, plot_speed_distribution(sys[i], 5))
		push!(pt, plot_trajectories(sys[i], collect(1:64)) )
	end

	# plot speed distribution and trajectory plots
	plot(
		ps[1], ps[2], ps[3], 
		pt[1], pt[2], pt[3], 
		layout = (2,3),
		size = (1920,1080)
	)
end

# DEMO 3: MELTING TRANSITION
function demo_3()
	println("\nDEMO 3: MELTING TRANSITION")
	println("----------------------------------------")

	# initialize system of particles on square lattice with zero velocity
	sys = ParticleSystem(16, 4.0, 5.0)
	set_square_lattice_positions!(sys)
	print_system_data(sys)
	p1 = plot_positions(sys)

	# evolve the system and watch them "crystallize" 
	# into a triangular lattice formation
	evolve!(sys, 20.0)
	print_system_data(sys)
	p2 = plot_trajectories(sys, collect(1:16))

	# now, increase the temperature of the system by giving the particles
	# some velocity. evolve the system and plot the trajectories.
	set_random_velocities!(sys)
	evolve!(sys, 60.0)
	print_system_data(sys)
	p3 = plot_trajectories(sys, collect(1:16))

	# some more plots
	p4 = plot_energy(sys, 0.0)
	p5 = plot_temperature(sys)
	p6 = plot_speed_distribution(sys, 20)

	plot(
		p1, p2, p3, p4, p5, p6,
		layout = (2,3), 
		size = (1280,720)
	)
end
