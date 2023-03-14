use loops, not vectorised expressions
--------------------------------------------------------------------------------

don't create arrays in loops
--------------------------------------------------------------------------------

  i.e. don't write 

  for ...
    x = mean(x, params)
  end

  instead use

  x = Array(Float64, ...)
  for ...
    mean!(x, params)
  end

  ie preallocate the array and use the functions to update them.

profile code
--------------------------------------------------------------------------------

better to have a struct of arrays rather than an array of structs
--------------------------------------------------------------------------------

bad:

struct F
  a::Int
  b::Int
end

good:

struct FA
  a::Vector{Int}
  b::Vector{Int}
end


type annotations only when necessary
--------------------------------------------------------------------------------

- if a function _can_ run with any types, then leave it be. only add type 
  annotations you know the function can run with only one input type
