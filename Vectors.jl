module Vectors

function dot(v1::Vector{Float64}, v2::Vector{Float64})
	dotProduct = 0
	for i in 1:length(v1)
		dotProduct += v1[i] * v2[i]
	end
	return dotProduct
end

function cross(v1::Vector{Float64}, v2::Vector{Float64})
	v3 = zeros(3)
	v3[1] = v1[2] * v2[3] - v1[3] * v2[2]
	v3[2] = v1[3] * v2[1] - v1[1] * v2[3]
	v3[3] = v1[1] * v2[2] - v1[2] * v2[1]
	return v3
end

function mag_sq(v::Vector{Float64})
	magSq = 0
	for i in 1:length(v)
		magSq += v[i] * v[i]
	end
	return magSq
end

function mag(v::Vector{Float64})
	return sqrt(mag_sq(v))
end

function unit(v::Vector{Float64})
	return v ./ mag(v)
end

export dot, mag, mag_sq, unit

end
