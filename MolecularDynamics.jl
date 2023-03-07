# molecular dynamics in 2d

push!(LOAD_PATH, "./")
using Statistics
using Vectors
using Plots

mutable struct ParticleSystem
	N::Int64			# number of particles
	L::Float64			# length of square side
	T₀::Float64			# initial temperature
	t::Float64			# system time
	state::Vector{Float64}		# state space array
	steps::Int64
	tPoints::Vector{Float64}	# array of time steps added to during integration
	sampleTimePoints::Vector{Float64}
	energyPoints::Vector{Float64}	# list of energy, sampled every sampleInterval steps
	tempPoints::Vector{Float64}
	tempAccumulator::Float64
	squareTempAccumulator::Float64
	virialAccumulator::Float64
	xPoints
	vPoints
	forceType::String
end

function ParticleSystem(N::Int64=64, L::Float64=10.0, T₀::Float64=1.0)
	t = 0.0
	tPoints = [ 0.0 ]
	steps = 0
	energyPoints = []
	sampleTimePoints = []
	tempPoints = [ T₀ ]
	tempAccumulator = 0.0
	squareTempAccumulator = 0.0
	virialAccumulator = 0.0
	state = zeros(4N)	# state space array, [x1, y1, x2, y2, ..., vx1, vy1, ...]
	xPoints = []
	vPoints = []
	forceType = "lennardJones"

	return ParticleSystem(
		N, L, T₀, t, state, 
		steps, tPoints, sampleTimePoints, energyPoints, 
		tempPoints, tempAccumulator, squareTempAccumulator, 
		virialAccumulator, xPoints, vPoints, forceType
		)
end

# some useful "views" of the state array
@views positions(state) = state[1:Int64(length(state)/2)]
@views velocities(state) = state[(Int64(length(state)/2)+1):end]
@views xcomponent(vector) = vector[1:2:end]
@views ycomponent(vector) = vector[2:2:end]
@views particle(n, vector) = [ vector[2n-1], vector[2n] ]
@views numParticles(state) = length(state) / 4

# INITIALIZATION
################################################################################

function random_positions!(sys::ParticleSystem)
	positions(sys.state) .= rand(2*sys.N) .* sys.L
	cool(sys)
end

function triangular_lattice_positions!(sys::ParticleSystem)
end

function rectangular_lattice_positions!(sys::ParticleSystem)
	nx = Int64(sqrt(sys.N))		# num lattice points per axis
	ny = nx				# square lattice
	ax = sys.L / nx		# lattice point separation
	ay = sys.L / ny

	for i in 0:(nx-1)
		for j in 0:(ny-1)
			sys.state[2*(i*ny+j)+1] = (i + 0.5) * ax
			sys.state[2*(i*ny+j)+2] = (j + 0.5) * ay
		end
	end
end

function random_velocities!(sys::ParticleSystem)
	velocities(sys.state) .= rand(2*sys.N) .- 0.5
	zero_total_momentum!(velocities(sys.state))
	T = temperature(sys.state)
	velocities(sys.state) .*= sqrt(sys.T₀/T)
end

function zero_total_momentum!(velocities)
	xcomponent(velocities) .-= mean(xcomponent(velocities))
	ycomponent(velocities) .-= mean(ycomponent(velocities))
end


# FORCES
################################################################################

function minimum_image(x::Vector{Float64}, L::Float64)
	halfL = 0.5 * L
	for i in 1:length(x)
		if x[i] > halfL
			x[i] -= L
		elseif x[i] < -halfL
			x[i] += L
		end
	end
	#return (x .+ halfL) .% L .- halfL
	return x
end

function force(sys::ParticleSystem)
	if sys.forceType == "lennardJones"
		force, virial = lennard_jones_force(sys)
	elseif sys.forceType == "powerLaw"
		force, virial = power_law_force(sys)
	end

	sys.virialAccumulator += virial
	return force
end

function lennard_jones_force(sys::ParticleSystem)
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
		Threads.@threads for j in (i+1):N
			dx = minimum_image(x .- x[i], L)
			dy = minimum_image(y .- y[i], L)

			r2inv = 1.0 ./ (dx.^2 .+ dy.^2 .+ tiny)
			f = 48.0 .* r2inv.^7 - 24.0 .* r2inv.^4
			fx = dx .* f
			fy = dy .* f

			fx[i] = 0.0	# no self force
			fy[i] = 0.0


#			println("i = ", i, ", j = ", j)
#			println(dx)
#			println(dy)
#			println(f)
#			println(sum(fx))
#			println(sum(fy))
			#force[2*i-1] = sum(fx)
			#force[2*i] = sum(fy)
			#println(force)

			virial += dot(fx,dx) + dot(fy,dy)
		end
	end

	return force, (0.5 * virial)
end

# TIME EVOLUTION
################################################################################


# MEASUREMENTS
################################################################################

function kinetic_energy(state::Vector{Float64})
	return 0.5 * mag_sq(velocities(state))
end

function temperature(state::Vector{Float64})
	return kinetic_energy(state) / numParticles(state)
end

# GRAPHS
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



# RUN
################################################################################


function console_log(sys::ParticleSystem)
	println(  "  num threads:		", Threads.nthreads(), "\n"),
	println(
		  "  time:			", sys.t, 
		"\n  total energy:		", energy(sys), 
		"\n  temperature:		", temperature(sys),
		)

	if sys.steps > 0
		println(
			"\n  mean energy:		", mean_energy(sys), 
			"\n  standard dev:		", std_energy(sys),
			"\n  C_V:			", heat_capacity(sys),
			"\n  PV/NkT:			", mean_pressure(sys),
			)
	end

end

sys = ParticleSystem(16, 4.0, 1.0)

rectangular_lattice_positions!(sys)
random_velocities!(sys)
plot_positions(sys)
