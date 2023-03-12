push!(LOAD_PATH, "./")
import Vectors: mag
using Plots

function square_lattice_2d(N::Int64, L::Float64)
	n = Int64(floor(sqrt(N))) # number of particles along each axis
	a = L / n # spacing between particles

	coordinates = [ zeros(2) for i in 1:N ]
	num = 1
	for i in 1:n
		for j in 1:n
			coordinates[num] = [a*(i-0.5), a*(j-0.5)]
			num += 1
		end
	end
	return coordinates
end

function minimum_separation_vector(p1::Vector{Float64}, p2::Vector{Float64}, L::Float64)
	r12 = p2 .- p1 # vector from point p1 to point p2
	println("r12 = ", r12)

	for i in 1:length(r12)
		if r12[i] > (L/2)
			r12[i] -= L
		elseif r12[i] < -(L/2)
			r12[i] += L
		end
	end

	println("r12 = ", r12)
	return r12
end

function demo(p1::Int64, p2::Int64; N::Int64=25, L::Float64=10.0)
	lattice = square_lattice_2d(N, L)

	xdata = [ lattice[i][1] for i in 1:N ]
	ydata = [ lattice[i][2] for i in 1:N ]

	r1 = lattice[p1]
	r2 = lattice[p2]
	r12 = minimum_separation_vector(r1, r2, L)

	println(r1)
	println(r2)
	println(r2 .- r1)
	println(r12)

	plot(size=(800,800), title="min distance between $p1 and $p2")
	scatter!(xdata, ydata, 
		markersize = 12,
		label = "lattice points",
		markercolor = :grey,
		)
	scatter!([ r1[1] ], [ r1[2] ], 
		markercolor = :black, 
		markersize = 16,
	        label = "point $p1"	
		)
	scatter!([ r2[1] ], [ r2[2] ], 
		markercolor = :green, 
		markersize = 16,
	        label = "point $p2"	
		)
	plot!([ r1[1], r1[1] + r12[1] ], 
		[ r1[2], r1[2] + r12[2] ], 
		linewidth=4, 
		arrow=true,
		color=:black,
		label = "MINIMUM DISTANCE",
		)

	if (r2 .- r1) != r12
		scatter!(  [ r1[1] + r12[1] ], [ r1[2] + r12[2] ] ,
			markercolor = :red,
			markersize = 16,
			label = "projected",
			)

		plot!( [ r1[1], r2[1] ], [ r1[2], r2[2] ] ,
			linewidth=5, 
			arrow=true,
			color=:grey,
			label = "'true' distance",
			)
	end
	plot!()
end
