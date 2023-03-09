module Render
	
using GLMakie
using ParticleSystems

#ptx = Observable
#pty = []

function start_render(sys::ParticleSystem)
	ptx = Observable( sys.state[1:2:2*sys.N] )
	pty = Observable( sys.state[2:2:2*sys.N] )

	fontsize_theme = Theme(fontsize=40)
	set_theme!(fontsize_theme)
	fig = Figure(resolution=(600,400))
	ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y")

	for n = 1:sys.N
		GLMakie.scatter!(
			ax1, 
			ptx, 
			pty, 
			markersize=40.0
		)
	end
	limits!(ax1, 0, sys.L, 0, sys.L)
	display(fig)
end

function update_render(sys::ParticleSystem)
	ptx[] = sys.state[1:2:2*sys.N]
	pty[] = sys.state[2:2:2*sys.N]
	yield()
end

export start_render, 
	update_render

end
