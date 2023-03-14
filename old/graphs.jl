
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
