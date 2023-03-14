module MolecularDynamics

using Statistics
using StatsPlots

mutable struct ParticleSystem
	N::Int64			# number of particles
	L::Float64			# square box side length
	T₀::Float64			# initial temperature
	t::Float64			# system time
	dt::Float64			# time step
	state::Vector{Float64}		# state space array
	steps::Int64			# number of steps
	sampleInterval::Int64		# interval for sampling data
	timeData::Vector{Float64}	# array of sampled time points
	energyData::Vector{Float64}	# array of sampled energy values
	tempData::Vector{Float64}	# array of sampled temperature values
	tempAccumulator::Float64	# temperature accumulator
	squareTempAccumulator::Float64	# T^2 accumulator
	virialAccumulator::Float64	# virial accumulator
	xData::Vector{Vector{Float64}} # array of sampled position data
	vData::Vector{Vector{Float64}} # array of sampled velocity data
	forceType::String		# force
end

function ParticleSystem(N::Int64=64, L::Float64=10.0, T₀::Float64=1.0)
	t = 0.0
	dt = 0.001
	state = zeros(4N) # state space array, [x1,y1,x2,y2,...,vx1,vy1,...]
	steps = 0
	timeData = Float64[]
	energyData = Float64[]
	sampleInterval = 100
	tempData = Float64[]
	tempAccumulator = 0.0
	squareTempAccumulator = 0.0
	virialAccumulator = 0.0
	xData = Vector{Float64}[]
	vData = Vector{Float64}[]
	forceType = "lennardJones"

	return ParticleSystem(
		N, 
		L, 
		T₀, 
		t, 
		dt, 
		state, 
		steps, 
		sampleInterval,
		timeData, 
		energyData, 
		tempData, 
		tempAccumulator, 
		squareTempAccumulator, 
		virialAccumulator, 
		xData, 
		vData, 
		forceType
	)
end

# some useful "views" of the state array
# (read the performance tips chapter of the julia manual)
@views positions(state) = state[ 1:Int64(length(state)/2) ]
@views velocities(state) = state[ (Int64(length(state)/2)+1):end ]
@views xcomponent(vector) = vector[ 1:2:end ]
@views ycomponent(vector) = vector[ 2:2:end ]
@views particle(n, vector) = [ vector[2n-1], vector[2n] ]

# INITIALIZATION
################################################################################

function set_random_positions!(sys::ParticleSystem)
	println("\tposition configuration: random")
	positions(sys.state) .= rand(2*sys.N) .* sys.L
	cool!(sys)
end

function set_square_lattice_positions!(sys::ParticleSystem)
	println("\tposition configuration: square lattice")

	n = Int64(floor(sqrt(sys.N))) # num lattice points per axis
	latticeSpacing = sys.L / n

	if sys.N != n^2
		println("\t\toops... your chosen N=$(sys.N) is not a square number") 
		println("\t\t-> resetting N to $(n^2).")
		sys.N = n^2
		sys.state = zeros(4 * sys.N)
	end

	for i in 0:(n-1)
		for j in 0:(n-1)
			sys.state[2*(i*n+j)+1] = (i + 0.5) * latticeSpacing
			sys.state[2*(i*n+j)+2] = (j + 0.5) * latticeSpacing
		end
	end
end

function set_triangular_lattice_positions!(sys::ParticleSystem)
end

function add_position_jitter!(sys::ParticleSystem, jitter::Float64=0.5)
	println("\tadding a wee bit of random jitter to particle positions...")

	for i = 1:length(positions(sys.state))
		sys.state[i] += rand() - jitter
	end
end

function set_random_velocities!(sys::ParticleSystem)
	println("\tvelocity configuration: random")

	velocities(sys.state) .= rand(2*sys.N) .- 0.5
	zero_total_momentum!(sys)
	velocities(sys.state) .*= sqrt(sys.T₀/temperature(sys))
end

function zero_total_momentum!(sys::ParticleSystem)
	xcomponent(velocities(sys.state)) .-= 
		mean(xcomponent(velocities(sys.state)))
	ycomponent(velocities(sys.state)) .-= 
		mean(ycomponent(velocities(sys.state)))
end


# FORCES / POTENTIALS
################################################################################

function force(sys::ParticleSystem)
	if sys.forceType == "lennardJones"
		force, virial = lennard_jones_force(sys)
	elseif sys.forceType == "powerLaw"
		force, virial = power_law_force(sys)
	end

	sys.virialAccumulator += virial

	return force
end

# the minimum image approximation
# (periodic boundary conditions)
function minimum_image(xij::Float64, L::Float64)
	if xij > (L/2)
		xij -= L
	elseif xij < -(L/2)
		xij += L
	end
	return xij
end

function lennard_jones_force(sys::ParticleSystem)
	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	virial = 0.0
	force = zeros(2*sys.N)

	Threads.@threads for i = 1:(sys.N-1)
		for j = (i+1):sys.N
			dx = minimum_image(x[i] - x[j], sys.L)
			dy = minimum_image(y[i] - y[j], sys.L)

			r2inv = 1.0 / (dx^2 + dy^2)
			f = 48.0 * r2inv^7 - 24.0 * r2inv^4
			fx = dx * f
			fy = dy * f

			force[2*i-1] += fx
			force[2*i] += fy
			force[2*j-1] -= fx
			force[2*j] -= fy

			virial += fx * dx + fy * dy
		end
	end

	return force, 0.5 * virial
end

function lennard_jones_potential(sys::ParticleSystem)
	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	U = 0.0

	Threads.@threads for i in 1:(sys.N-1)
		for j in (i+1):sys.N
			dx = minimum_image(x[i] - x[j], sys.L)
			dy = minimum_image(y[i] - y[j], sys.L)

			r2inv = 1.0 / (dx^2 + dy^2)
			U += r2inv^6 - r2inv^3
		end
	end
	return 4.0 * U
end

function power_law_force(sys::ParticleSystem)
end

function power_law_potential(sys::ParticleSystem)
end

# TIME EVOLUTION
################################################################################

function keep_particles_in_box!(sys::ParticleSystem)
	for i in 1:(2*sys.N)
		if positions(sys.state)[i] > sys.L
			positions(sys.state)[i] -= sys.L
		elseif positions(sys.state)[i] < 0.0
			positions(sys.state)[i] += sys.L
		end
	end

#	# another way of doing this: use the ternary operator
#	for i in 1:(2 * sys.N)
#		positions(sys.state)[i] < 0.0 ?  
#		positions(sys.state)[i] % sys.L + sys.L : 
#		positions(sys.state)[i] % sys.L
#	end
end

function verlet_step!(sys::ParticleSystem)
	# compute acceleration at current time
	acceleration = force(sys)

	# compute positions at t + dt
	positions(sys.state) .+= 
		velocities(sys.state) .* sys.dt .+ 
		0.5 .* acceleration .* (sys.dt)^2

	# enforce boundary conditions
	# (basically check if any particles left the box and put them back)
	# see function implementation for deets
	keep_particles_in_box!(sys)

	# compute velocities at t + dt
	velocities(sys.state) .+= 
		0.5 * sys.dt .* (acceleration + force(sys))
end

function evolve!(sys::ParticleSystem, runtime::Float64=10.0)
	numsteps = Int64(abs(runtime/sys.dt) + 1)

	print_evolution_message(runtime, numsteps)

	@time for step in 1:numsteps
		verlet_step!(sys)
		zero_total_momentum!(sys)

		if (step % sys.sampleInterval == 1)
			push!(sys.timeData, sys.t)
			push!(sys.energyData, energy(sys))
			push!(sys.xData, positions(sys.state))
			push!(sys.vData, velocities(sys.state))

			T = temperature(sys)
			push!(sys.tempData, T)
			sys.tempAccumulator += T
			sys.squareTempAccumulator += T^2
		end

		sys.t += sys.dt
		sys.steps += 1
	end
	println("done.") 
end

function reverse_time!(sys)
	sys.dt *= -1
	println("\ntime reversed! dt = $(sys.dt)")
end

function cool!(sys::ParticleSystem, cooltime::Float64=1.0)
	numsteps = Int64(cooltime/sys.dt)
	for step in 1:numsteps
		verlet_step!(sys)
		velocities(sys.state) .*= (1.0 - sys.dt)
	end
	reset_statistics!(sys)
end

# MEASUREMENTS
################################################################################

function kinetic_energy(sys::ParticleSystem)
	return 0.5 * sum(velocities(sys.state) .* velocities(sys.state))
end

function potential_energy(sys::ParticleSystem)
	return lennard_jones_potential(sys)
end

function temperature(sys::ParticleSystem)
	return kinetic_energy(sys) / sys.N
end

function energy(sys::ParticleSystem)
	return potential_energy(sys) + kinetic_energy(sys)
end

# STATISTICS
################################################################################

function reset_statistics!(sys::ParticleSystem)
	sys.steps = 0
	sys.tempAccumulator = 0.0
	sys.squareTempAccumulator = 0.0
	sys.virialAccumulator = 0.0
	sys.xData = []
	sys.vData = []
end

function mean_temperature(sys::ParticleSystem)
	return sys.tempAccumulator / sys.steps
end

function mean_square_temperature(sys::ParticleSystem)
	return sys.squareTempAccumulator / sys.steps
end

function mean_pressure(sys::ParticleSystem)
	# factor of half because force is calculated twice each step
	meanVirial = 0.5 * sys.virialAccumulator / sys.steps
	return 1.0 + 0.5 * meanVirial / (sys.N * mean_temperature(sys))
end

function heat_capacity(sys::ParticleSystem)
	meanTemperature = mean_temperature(sys)
	meanSquareTemperature = mean_square_temperature(sys)
	σ2 = meanSquareTemperature - meanTemperature^2
	denom = 1.0 - σ2 * sys.N / meanTemperature^2
	return sys.N / denom
end

function mean_energy(sys::ParticleSystem)
	return mean(sys.energyData)
end

function std_energy(sys::ParticleSystem)
	return std(sys.energyData)
end

# MATH / ADDITIONAL FUNCTIONS
################################################################################

function dot(v1::Vector{Float64}, v2::Vector{Float64})
	return sum(v1 .* v2)
end

# GRAPHS
################################################################################

function initialize_plot()
	plot(
		size=(800,800), 
		titlefontsize=12, 
		guidefontsize=12,
	)
end

function plot_positions(sys::ParticleSystem)
	initialize_plot()
	for n = 1:sys.N
		scatter!(
			[ xcomponent(positions(sys.state))[n] ], 
			[ ycomponent(positions(sys.state))[n] ], 
			markersize = 4.0,
			markercolor = n,
			markerstrokewidth = 0.4,
			grid = true,
			framestyle = :box,
			legend = false,
		)
	end
	xlims!(0, sys.L)
	ylims!(0, sys.L)
	xlabel!("x")
	ylabel!("y")
	title!("positions at t=$(round(sys.t, digits=4))")
end

function plot_trajectories(sys::ParticleSystem, particles::Vector{Int64}=[ 1 ])
	initialize_plot()
	for n = 1:sys.N
		scatter!(
			[ xcomponent(positions(sys.state))[n] ], 
			[ ycomponent(positions(sys.state))[n] ], 
			markersize = 4.0,
			markercolor = n,
			markerstrokewidth = 0.4,
			grid = true,
			framestyle = :box,
			legend = false,
		)
	end

	for n in collect(particles)
		xdata = [ sys.xData[i][2n-1] for i in 1:length(sys.xData) ]
		ydata = [ sys.xData[i][2n] for i in 1:length(sys.xData) ]

		# plot trajectory line for nth particle
		scatter!(
			xdata, 
			ydata,
			color = n,
			#markerstrokewidth = 0,
			markerstrokecolor = n,
			markersize = 0.7,
			markeralpha = 0.5,
			label = false,
			widen = false,
		)

		# plot initial position for nth particle
		scatter!(
			[ sys.xData[1][2n-1] ], 
			[ sys.xData[1][2n] ],
			markersize = 4.0, 
			markercolor = n,
			markerstrokewidth = 0.4,
			markeralpha = 0.3,
			#label = "pcl. $n @t=t₀",
			widen = false,
		)

		# plot final position for nth particle
		scatter!(
			[ sys.xData[end][2n-1] ], 
			[ sys.xData[end][2n] ],
			markersize = 4.0, 
			markercolor = n,
			markerstrokewidth = 0.4,
			markeralpha = 1.0,
			#label = "pcl $n @t=t",
			widen = false,
		)
	end
	title!("positions & trajectories at time t=$(round(sys.t, digits=2))")
	plot!()
end

function plot_temperature(sys::ParticleSystem)
	initialize_plot()
	plot!(
		sys.timeData, 
		sys.tempData,
		#widen = true,
	)
	ylims!(
		mean(sys.tempData) - std(sys.tempData) * 20, 
		mean(sys.tempData) + std(sys.tempData) * 20, 
	)
	xlabel!("t")
	ylabel!("T(t)")
	title!("temperature vs time")

end

function plot_energy(sys::ParticleSystem, ylimit::Float64=1.0)
	initialize_plot()
	plot!(
		sys.timeData, 
		sys.energyData,
		#widen = true,
	)
	ylims!(
		#ylimit * (mean(sys.energyData) - 1), 
		#ylimit * (mean(sys.energyData) + 1)
		mean(sys.energyData) - std(sys.energyData) * 10, 
		mean(sys.energyData) + std(sys.energyData) * 10, 
	)
	xlabel!("t")
	ylabel!("E(t)")
	title!("energy vs time")
end

function plot_speed_distribution(sys::ParticleSystem, numSamples::Int64=5)
	initialize_plot()

	numDataPoints = Int64(length(sys.vData))
	interval = Int64(floor(numDataPoints / numSamples))

	samples = collect(1:interval:numDataPoints)
	for s in samples
		speed = sqrt.(
			xcomponent(sys.vData[s]).^2 .* 
			ycomponent(sys.vData[s]).^2
		)
		density!(
			sys.vData[s],
			normalize = :pdf, 
			label = "t = $(round(sys.timeData[s], digits=2))",
		)
	end
	xlabel!("speed")
	title!("speed distribution")
end

# CONSOLE PRINT DATA
################################################################################

function print_hello()
	println("\nmolecular dynamics!")
	println("number of threads: ", Threads.nthreads())
end

function print_bonjour()
	println("\nbonjour")
end

function print_system_parameters(sys::ParticleSystem)
	println("\nsystem parameters:")
	println("\tN =  $(sys.N)   (number of particles)")
	println("\tL =  $(sys.L)   (side length of square box)")
	println("\tDT = $(sys.dt)  (time step)")
end

function print_system_data(sys::ParticleSystem)
	println("\nsystem data at time t=$(round(sys.t, digits=4))")

	if sys.steps == 0
		println("\ttemperature:     $(temperature(sys))")
		println("\tenergy:          $(energy(sys))")
	else
		println("\tsteps evolved:   $(sys.steps)")
		println("\ttemperature:     $(temperature(sys))")
		println("\tenergy:          $(energy(sys))")
		println("\tmean energy:     $(mean_energy(sys))")
		println("\tstd energy:      $(std_energy(sys))")
		println("\theat capacity:   $(heat_capacity(sys))")
		println("\tPV/NkT:          $(mean_pressure(sys))")
	end
end

function print_evolution_message(runtime, numsteps)
	println("\nevolving...")
end

# DEMOS
################################################################################


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

	evolve!(sys, 20.0)
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


end # module
