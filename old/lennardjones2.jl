function minimum_image(x::Vector{Float64}, L::Float64)
	halfL = 0.5 * L
	for i in 1:length(x)
		if x[i] > halfL
			x[i] -= L
		elseif x[i] < -halfL
			x[i] += L
		end
	end
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
	N = sys.N
	L = sys.L
	tiny = 1.0e-40
	virial = 0.0
	force = zeros(2N)

	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))

	Threads.@threads for i in 1:N
		for j in (i+1):N
			dx = minimum_image(x .- x[i], L)
			dy = minimum_image(y .- y[i], L)

			r2inv = 1.0 ./ (dx.^2 .+ dy.^2 .+ tiny)
			f = 48.0 .* r2inv.^6 - 24.0 .* r2inv.^4
			fx = dx .* f
			fy = dy .* f

			fx[i] = 0.0	# self force is zero
			fy[i] = 0.0

			virial += dot(fx,dx) + dot(fy,dy)
		end
	end

	return force, 0.5 * virial
end

function lennard_jones_potential(sys::ParticleSystem)
	N = sys.N
	L = sys.L
	tiny = 1.0e-40

	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	U = 0.0

	Threads.@threads for i in 1:N
		dx = minimum_image(x[i] .- x, L)
		dy = minimum_image(y[i] .- y, L)

		r2inv = 1.0 / (dx.^2 .+ dy.^2 .+ tiny)
		dU = r2inv.^6 .- r2inv.^3
		dU[i] = 0.0	# self interaction is zero
		U += sum(dU)
	end

	return 2.0 * U
end
