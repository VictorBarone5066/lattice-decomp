Lattice Decomposition

Decomposes a (presumably large) set of lattices into unique chemical environments.
Have a look at 'inputexample.ldc' for an example of how to build an input file.

This program writes the matrix A in the equation Ax = b, where A_ij is the count of
chemical environment number j that configuration i contains.  
A can be written in dense:
A(1)(2) = 3 means that configuration 1 has three instances of chem env #2
...or sparse:
   row col val col val col val ... 
formats.

A chemical environment is a set of atoms within some radius defined by the distances 
between them and the element of each one.  The decomposition of a set of configurations
into chemical environments is useful since many physical quantities (energy especially)
can be written as a linear combo of local environments.  The number of unique envs is 
much less than the number of possible full-cell configurations (as long as your cutoff 
radius is reasonable!), and so such a decomposition allows for accurate estimation of 
physical quantities of "many" possible configurations by fitting only "a few" data points.  

This is based off of the "Cluster Expansion" idea.

Usage notes:

-I compile with the following command:
   "gcc -o decomp -Ofast -std=c11 -w -lm Constants.c IO.c Lattice.c Memory.c Structs.c Main.c extern/kdtree-master/kdtree.c"
which names the executable 'decomp'.  

-Assuming that you named your executable 'decomp', the command:
   "./decomp -i input.ldc -r 3.15"
will read input from the file input.ldc and OVERWRITE the radius within input.ldc with
the value 3.5 (the -r option is not required.  If left out, the input.ldc radius will
be used instead).

-To avoid slowdowns associated with console output, especially for long runs, do this:
   "./decomp -i input.ldc -r 3.15 >out.ldc"
this is especially important for remote connections (in which case you might also sandwich the above command
inbetween 'nohup' and '&').

-The files Constants.* include some important things that you may want to change.  Specifically:
--the "#define"s in Constants.h (see comments within actual file)
--the "const external"s in Constants.c (see comments within actual file)  

    
