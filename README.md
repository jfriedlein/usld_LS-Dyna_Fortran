# usld_LS-Dyna_Fortran
Basics to implement user-defined elements (usld, uel) in LS-Dyna with Fortran 

## What is this all about? 
LS-Dyna offers the interfaces and solvers to, among many other things, simulate mechanical systems and the related material behaviour. To obtain accuracte results we need to utilise adequate elemnet formulations. In case, standard available element formulations (ELFORM=1,-1,-2,...) cannot generate valid results, new user-defined (solid) elements can be implemented. The latter can be implemented as a standalone self-contained routine, in LS-Dyna called "resultant element", where the user implements the entire element formulation (shape functions, integrations, material model, ...) and is provided with the nodal coordinates and displacements, and has to compute the force vector, tangent matrix, and history update for an element. This guide introduces the basics to implement user-defined resultant solid elements in LS-Dyna using the standard Fortran interface. 

## Software requirements and suggestions 
The software requirements are similar to the [implementation of umats](https://github.com/jfriedlein/usrmat_LS-Dyna_Fortran), but here Linux is considered.

* An object version of the LS-Dyna version you wish to use. Everything outlined here refers to version R11.1, as this is one of the few versions that enable visualisation of history variables in LS-PrePost (as far as I know). Their might be slight differences compared to older version, e.g. where to find the files (see for instance [Older LS-Dyna versions/releases](https://github.com/jfriedlein/usrmat_LS-Dyna_Fortran?tab=readme-ov-file#older-ls-dyna-versionsreleases)). The object version is a compressed package (e.g. .zip) typically ending with _lib.zip. You can acquire this package from your LSTC support or, in case you possess the login credentials for the [https://files.dynamore.se/index.php/s/RPpD5rW5xmo65rX/authenticate/showShare](https://files.dynamore.se/index.php/s/RPpD5rW5xmo65rX/authenticate/showShare), you can download the versions where all available version (for Windows and Linux) are listed (e.g. the here used ‘ls-dyna_smp_d_R11_1_0_x64_redhat65_ifort160_sse2.zip’).
* For mpp versions of LS-Dyna you also need some MPP tools, such as MSMPI.

For coding:
* Under Windows you typically need Visual Studio and the Intel Fortran Compiler or OneAPI. Check the readme.txt in the _lib.zip, which states the required versions for both. If you want to avoid trouble, adhere to the tools and version given in there. Especially older versions of LS-Dyna like R9/R10 require rather old Visual Studio and Intel Fortran versions, so make sure you can get access to these dusty rusty things.
*  You could use any text editor for the typing the source code. Or you can use an IDE like Visual studio, Visual studio code (VS-code, e.g. with "Modern Fortran" extension), ... 
* For the compilation of the Fortran .f/.F files you need a Fortran compiler, e.g. Intel Parallel Studio XE 2017 or OneAPI. Be aware of the dependencies of the Fortran compiler and Visual Studio. As long as you stick with the versions in readme.txt you should be fine. 

## Test setup and compiler 
see [Test setup and compiler](https://github.com/jfriedlein/usrmat_LS-Dyna_Fortran?tab=readme-ov-file#test-setup-and-compiler)

## Implementation 
Please be aware that Fortran has some “features” that might (most certainly) be unknown or unexpected to programmers used to “more modern” languages, such as C++, Matlab, Python, ... For a quick starter, see [A few notes on FORTRAN](https://github.com/jfriedlein/usrmat_LS-Dyna_Fortran?tab=readme-ov-file#a-few-notes-on-fortran).

## Our first resultant user-defined solid element
1. Open your working directory (the folder with the unpacked object version, e.g. ls-dyna_smp_d_R11_1_0_139588_winx64_ifort2017vs2017_lib) in your IDE, e. g. VS-code.
2. Implement your element formulation code (computation of force vector, stiffness matrix, history variables ...) in the file dyn21usld.f, for instance a fully integrated linear elastic element. We code our model in the first unused user-solid routine usld_e101.

code and explanation

## Some notes on proper simulation and solver settings for testing umats 

## Outline of the interface for umat and utan 

## Generalised interface for use of separate element subroutines

## Additional (extra) degrees of freedom xdofs
