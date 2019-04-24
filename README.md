# alphaPredictableComponents

## Summary
A Julia package for working with spacetime field (`stfield`) objects.
An `stfield` object is a data field with data, longitude, latitude, time.
Currently, an `stfield` object (with a particular `extent`) can be created from a `NetCDF` file.
An `stfield` object can be manipulated with commands such as `convert` (to a `DataFrame` or `Matrix`), 
or `copy`. 

## Installation
```julia
Pkg.add("https://github.com/Mattriks/alphaPredictableComponents.jl")
```

## Documentation
Documentation with examples can be found [here](../master/docs/build/index.md)
