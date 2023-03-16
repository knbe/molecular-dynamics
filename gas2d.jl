include("statespace.jl")

const global D::UInt8 = 2

function set_square_lattice_positions!(sys::StateSpace, L::Float64)
	println("\tposition configuration: square lattice")

	n = Int(floor(sqrt(sys.N)))
	latticeSpacing = L / n

	if sys.N != n^2
		println("\t\toops... your chosen N=$(sys.N) is not a square number") 
		println("\t\t-> resetting N to $(n^2).")
		sys.N = n^2
		sys.position = zeros(D*sys.N)
		sys.velocity = zeros(D*sys.N)
	end

	for i in 0:(n-1)
		for j in 0:(n-1)
			#component(i, sys.position
			#sys.position[D*(i*n+j)] = (i+0.5) * latticeSpacing
			#sys.velocity[D*(i*n+j)] = (j+0.5) * latticeSpacing
		end
	end
end

