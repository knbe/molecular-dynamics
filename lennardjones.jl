

function force(sys::ParticleSystem)
	if sys.forceType == "lennardJones"
		force, virial = lennard_jones_force(sys)
	elseif sys.forceType == "powerLaw"
		force, virial = power_law_force(sys)
	end

	sys.virialAccumulator += virial

	return force
end

# the minimum image approximation
# (periodic boundary conditions)
function minimum_image(xij::Float64, L::Float64)
	if xij > (L/2)
		xij -= L
	elseif xij < -(L/2)
		xij += L
	end
	return xij
end


function lennard_jones_force(sys::ParticleSystem)
	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	virial = 0.0
	force = zeros(2*sys.N)

	Threads.@threads for i = 1:(sys.N-1)
		for j = (i+1):sys.N
			dx = minimum_image(x[i] - x[j], sys.L)
			dy = minimum_image(y[i] - y[j], sys.L)

			r2inv = 1.0 / (dx^2 + dy^2)
			f = 48.0 * r2inv^7 - 24.0 * r2inv^4
			fx = dx * f
			fy = dy * f

			force[2*i-1] += fx
			force[2*i] += fy
			force[2*j-1] -= fx
			force[2*j] -= fy

			virial += fx * dx + fy * dy
		end
	end

	return force, 0.5 * virial
end

function lennard_jones_potential(sys::ParticleSystem)
	x = xcomponent(positions(sys.state))
	y = ycomponent(positions(sys.state))
	U = 0.0

	Threads.@threads for i in 1:(sys.N-1)
		for j in (i+1):sys.N
			dx = minimum_image(x[i] - x[j], sys.L)
			dy = minimum_image(y[i] - y[j], sys.L)

			r2inv = 1.0 / (dx^2 + dy^2)
			U += r2inv^6 - r2inv^3
		end
	end
	return 4.0 * U
end

