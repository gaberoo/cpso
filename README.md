 
# PARTICLE SWARM OPTIMIZATION IN C++

Copyright (c) 2011-2015 Gabriel E Leventhal

adapted from:
  SimPSOLib: A simple yet flexible PSO implementation in C++.
  Author: Tomas V. Arredondo


## Test run

The test function to optimize is

    f = -(x[0]-10)*(x[0]-10) - (x[1]-5)*(x[1]-5) + 10

The optimimum of this function is at `x* = (10,5)`.

Compile test:

    make libpso.a
    make test

Run test:

    ./test test.json
