
push!(LOAD_PATH, "./")
using ParticleSystem


function demo_1()
	render = 1
	lennardjones = 1
	sys = ParticleSystem(64, 10.0, 1.0)

	rectangular_lattice!(sys)
	random_velocities!(sys)
	print_console(sys)

	evolve!(sys, 10.0, 0.01, 100)
	plot_positions(sys.pos)
end

function demo_2()
	lennardjones = 1
	sys = ParticleSystem(144, 10.0, 1.0)

	rectangular_lattice!(sys)
	random_velocities!(sys)
	print_console(sys)

	evolve!(sys, 10.0, 0.01, 100)
	plot_positions(sys.pos)
end
