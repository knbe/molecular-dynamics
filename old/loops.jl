include("molecular.jl")

mutable struct ParticleSystem
	D::Int64	# dimensions
	N::Int64	# number of particles
	stateSpace::Vector{Float64}
end

function tmp(sys::ParticleSystem)
end
