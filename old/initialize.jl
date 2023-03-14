
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

