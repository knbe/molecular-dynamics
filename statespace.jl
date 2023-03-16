mutable struct StateSpace{T<:Vector{Float64}}
	N::Int64
	position::T
	velocity::T
end

function StateSpace(N::Int64)
	return StateSpace(N, zeros(D*N), zeros(D*N))
end

@views component(i, vector) = vector[ i:D:end ]
@views particle(i, vector) = vector[ (D*i-(D-1)):(D*i) ]
