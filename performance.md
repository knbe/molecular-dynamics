# use loops rather than vectorised expressions

julia can run this

```
for i = 1:length(x)
    x[i] += sqrt(a * v[i]) 
end
```

faster than

```
x .+= sqrt.(a .* v)
```

# use @views for slicing a large array

rather than creating a bunch of subarrays in a loop. 

https://docs.julialang.org/en/v1/manual/performance-tips/

# try not to create arrays in loops

as this increases memory allocation.

  i.e. don't write 

```

```

  for ...
    x = mean(x, params)
  end

  instead use

  x = Array(Float64, ...)
  for ...
    mean!(x, params)
  end

  ie preallocate the array and use the functions to update them.

# use the Profile module

https://docs.julialang.org/en/v1/manual/profile/

to check how many times a method gets called in a loop.
