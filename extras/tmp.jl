#macro integrator(type)
#	return :( println("integrator: ", $type) )
#end

# memory allocation

const global x::Vector{Float64} = rand(10000)

function sum_global()
	s = 0.0
	for i in x
		s += i
	end
	return s
end

function sum_arg(x)
	s = 0.0
	for i in x
		s += i
	end
	return s
end

# containers

a = Float64[]

# type declarations

mutable struct Thing{T<:AbstractFloat}
	a::T
end

# is better than
mutable struct OtherThing
	a::AbstractFloat
end

# annotate values from untyped locations
function foo(a::Array{Any,1})
	x = a[1]::Int64
end

# specialising
