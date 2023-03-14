module ParticleSystems


mutable struct ParticleSystem
	N::Int64		# number of particles
	L::Float64		# side length of square box
	T₀::Float64	
	t::Float64		# time
	state::Vector{Float64}	# state space vector
end

function ParticleSystem(N, L, T₀, D)
	t = 0.0
	state = zeros(D*N*2) # [x1,y1,z1,x2,y2,z2,...,vx1,vy1,vz1,...]
	return ParticleSystem(N, L, T₀, t, state)
end

export ParticleSystem

end
