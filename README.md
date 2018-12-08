
Readme
------

MicroC is a FE code used to model composite material microstructures.

Compilation
-----------

1. Download a PETSc Library from [www.mcs.anl.gov/petsc](www.mcs.anl.gov/petsc).
2. Set the environmental variables `PETSC_DIR` and `PETSC_ARCH` with

```bash
   export PETSC_DIR=<path>
   export PETSC_ARCH=<architecture>
```

3. Configure PETSc for not to run with MPI, compile and install

```bash
./configure --with-mpi=no
make all
```

Build MicroC with CMake:
-----------------------

1. Clone the repository [GitHub](https://github.com/GG1991/macroc)
2. `cd <cloned directory>`
3. `mkdir build` (can be also `build + anything`)
4. `cd build`
5. `cmake ..`
6. `make`

To build the optimized version:

```bash
cmake -DCMAKE_BUILD_TYPE=Release ..
```

and the debug version:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Testing
---------

In order follow DevOps a testing methodology is used. **MicroC** uses CTest toolt to perform tests. The recommended use of it consists in going first to a stable version of the code in the master branch and execute:

```bash
ctest --VV -O output_stable.dat
```

Then in the file `output_stable.dat` the stable reference output would have be written for all the tests cases. Finally after doing modifications if it is necessary run the tests cases and compare the solution running the same command but changing the output file name, for example:

```bash
ctest -VV -O output_develop.dat
```

You can compare the solution using Linux tools like `diff` or `vimdiff`.
