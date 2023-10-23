# Reduced model for the analysis of vacuum-enhanced air gap membrane distillation (V-AGMD) modules

This code implements the computations related to a reduced 0D model for steady state analysis of vacuum-enhanced air gap membrane
distillation (V-AGMD) modules. The code is implemented in C and uses the Portable, Extensible Toolkit for Scientific Computation
(PETSc) mainly to solve the system of nonlinear algebraic equations.

The code is implemented in a modular fashion, with folders encapsulating most of the logic for each needed component, that is
for thermophysical properties, entry data fetching, desalination module modeling, and system integration and numerical solution.
They communicate with each other through APIs defined in header files.

## Dependencies

We recommend GNU Make, C/C++ and Fortran compilers, and Git are installed before proceeding. In Ubuntu, these packages can be
installed with the following commands (for GCC and gfortran):

```bash
$ sudo apt install make gcc gfortran git
```

As already mentioned, this project depends on the Portable, Extensible Toolkit for Scientific Computation (PETSc). More
information on how to install and use PETSc can be found in https://petsc.org/release/

The code was developed using the following PETSc compile configuration:

```bash
$ ./configure --with-cc=gcc --with-cxx --with-fc=gfortran --with-fortran-bindings=0 --download-mpich --download-fblaslapack 
--download-superlu-dist --download-cmake --download-sundials2 --with-debugging=no
```

However, there is no need to include some extra modules and MPI (serial code) and a bare minimum configuration that may work
is the one that follows:

```bash
$ ./configure --with-cc=gcc --with-cxx --with-fc=gfortran --with-fortran-bindings=0 --with-mpi=0 --download-fblaslapack 
--with-debugging=no
```

## Downloading the code

To clone the repository, cd to a folder intended to hold the repository and type the following command in the terminal:

```bash
$ git clone https://github.com/labmems/vagmd0Dmodel.git vagmd0Dmodel
```

## Compiling the code

To compile the code, simply cd to the folder into which the files were cloned and type the following command in the terminal:

```bash
$ make build
```

## Getting help

To get help on how to run the binary for custom numerical configurations and entry data, simply cd to the folder into which
the files were cloned and type the following command in the terminal:

```bash
$ make help
```

## Running the binary

To run the binary, simply cd to the folder into which the files were cloned and type the following command in the terminal:

```bash
$ make run
```