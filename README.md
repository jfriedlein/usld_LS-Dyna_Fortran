# usld_LS-Dyna_Fortran
Basics to implement user-defined elements (usld, uel) in LS-Dyna with Fortran 

## What is this all about? 
LS-Dyna offers the interfaces and solvers to, among many other things, simulate mechanical systems and the related material behaviour. To obtain accuracte results we need to utilise adequate elemnet formulations. In case, standard available element formulations (ELFORM=1,-1,-2,...) cannot generate valid results, new user-defined (solid) elements can be implemented. The latter can be implemented as a standalone self-contained routine, in LS-Dyna called "resultant element", where the user implements the entire element formulation (shape functions, integrations, material model, ...) and is provided with the nodal coordinates and displacements, and has to compute the force vector, tangent matrix, and history update for an element. This guide introduces the basics to implement user-defined resultant solid elements using the Total Lagrangian formulation in LS-Dyna with the standard Fortran interface. 

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
2. Implement your element formulation code (computation of force vector, stiffness matrix, history variables ...) in the file dyn21usld.f. We code our model in the first unused user-solid routine usld_e101. Below an example for this routine is presented for a fully integrated linear elastic 3D solid element with 8-nodes and linear shape functions. The code is quite lengthy and not optimised as it aims to show the steps one-by-one with minimal use of separate functions.

```fortran
      subroutine usld_e101(force,stiff,ndtot,istif,
     . x1,x2,x3,x4,x5,x6,x7,x8,
     . y1,y2,y3,y4,y5,y6,y7,y8,
     . z1,z2,z3,z4,z5,z6,z7,z8,
     . xdof,
     . dx1,dx2,dx3,dx4,dx5,dx6,dx7,dx8,
     . dy1,dy2,dy3,dy4,dy5,dy6,dy7,dy8,
     . dz1,dz2,dz3,dz4,dz5,dz6,dz7,dz8,
     . dxdof,
     . hsv,ihsv,nhsv,
     . cm,lmc,
     . cmtrx,lft,llt,dener)
c
c     compute force (and stiffness) for user defined solid 101
c
      ! Load the UEL_helper module to retrieve the undeformed nodal coordinates below via "get_initialNodalCoords_1element_R102(*)"
       use UEL_helper
c
      ! default LS-Dyna include and declarations
       include 'nlqparm'
c      dimension force(nlq,ndtot),stiff(nlq,ndtot,ndtot)
       dimension force(nlq,*),stiff(nlq,ndtot,*)
       dimension x1(nlq),x2(nlq),x3(nlq),x4(nlq)
       dimension x5(nlq),x6(nlq),x7(nlq),x8(nlq)
       dimension y1(nlq),y2(nlq),y3(nlq),y4(nlq)
       dimension y5(nlq),y6(nlq),y7(nlq),y8(nlq)
       dimension z1(nlq),z2(nlq),z3(nlq),z4(nlq)
       dimension z5(nlq),z6(nlq),z7(nlq),z8(nlq)
       dimension xdof(nlq,8,*)
       dimension dx1(nlq),dx2(nlq),dx3(nlq),dx4(nlq)
       dimension dx5(nlq),dx6(nlq),dx7(nlq),dx8(nlq)
       dimension dy1(nlq),dy2(nlq),dy3(nlq),dy4(nlq)
       dimension dy5(nlq),dy6(nlq),dy7(nlq),dy8(nlq)
       dimension dz1(nlq),dz2(nlq),dz3(nlq),dz4(nlq)
       dimension dz5(nlq),dz6(nlq),dz7(nlq),dz8(nlq)
       dimension dxdof(nlq,8,*)
       dimension hsv(nlq,nhsv),ihsv(nlq,nhsv),cm(lmc)
       dimension cmtrx(nlq,21)
       dimension dener(nlq)
c
      ! custom declarations for this user-element
       real*8 :: x_undeformed(3,8)   ! undeformed coordinates as matrix of dofs for each node
       real*8 :: x_deformed(3,8)     ! deformed coordinates as matrix of dofs for each node
       real*8 :: u_displacement(3,8) ! displacements as matrix of dofs for each node
       real*8 :: Jacobi(3,3)         ! element's Jacobian matrix
       real*8 :: inv_Jac(3,3)        ! inverse Jacobian
       real*8 :: m_matrix(6,6)       ! material matrix (relation of stress and strain)
       real*8 :: B_matrix(6,24)      ! B-Matrix (gradients of shape functions)
       real*8 :: strain(6)           ! strain components in LS-Dyna vector notation (xx,yy,zz,xy,yz,zx)
       real*8 :: stress(6)           ! stress components in LS-Dyna vector notation (xx,yy,zz,xy,yz,zx)
       real*8 :: dN_dXi(8,3)         ! gradient of ansatz functions in local coordinates
       real*8 :: dN_dX (8,3)         ! gradient in x,y,z coordinates
       real*8 :: qp_coord(3,8)       ! quadrature point coordinates for full integration 2x2x2
       real*8 :: detJ                ! determinant of the Jacobi
       real*8 :: E_modul             ! Young's modulus
       real*8 :: Poisson             ! Poisson ratio
       real*8 :: fct_c_matrix        ! factor used for material matrix
       real*8 :: c11                 ! factor used for (11,22,33) components of material matrix
       real*8 :: c44                 ! factor used for (44,55,66) components of material matrix
       real*8 :: sqrt13              ! square root of 1./3.
       real*8 :: qp_weight
c
       integer :: i_el               ! iterator for loop over elements
       integer :: qp                 ! iterator for loop over quadrature points
c
      ! Loop over each element in the list of elements (indices "lft" to "llt")
       do i_el = lft,llt
          ! Init the variables to zero for a fresh start of each element "i_el"
           x_undeformed(:,:) = 0.
           x_deformed(:,:) = 0.
           u_displacement(:,:) = 0.
           Jacobi(:,:) = 0.
           inv_Jac(:,:) = 0.
           m_matrix(:,:) = 0.
           B_matrix(:,:) = 0.
           strain(:) = 0.
           stress(:) = 0.
           dN_dXi(:,:) = 0.
           dN_dX(:,:) = 0.
           qp_coord(:,:) = 0.
c
          ! group the deformed coordinates of each node (1...8) for the current element "i_el" into a vector
           x_deformed(:,1) = (/ x1(i_el),y1(i_el),z1(i_el) /)
           x_deformed(:,2) = (/ x2(i_el),y2(i_el),z2(i_el) /)
           x_deformed(:,3) = (/ x3(i_el),y3(i_el),z3(i_el) /)
           x_deformed(:,4) = (/ x4(i_el),y4(i_el),z4(i_el) /)
           x_deformed(:,5) = (/ x5(i_el),y5(i_el),z5(i_el) /)
           x_deformed(:,6) = (/ x6(i_el),y6(i_el),z6(i_el) /)
           x_deformed(:,7) = (/ x7(i_el),y7(i_el),z7(i_el) /)
           x_deformed(:,8) = (/ x8(i_el),y8(i_el),z8(i_el) /)
c
          ! obtain the initial/undeformed nodal coordinates for the current element "i_el"                                                   
           X_undeformed = transpose(
     &                     get_initialNodalCoords_1element_R102( i_el )
     &                    )
c
          ! compute the displacements from the deformed and undeformed coordinates for each node and dof
           u_displacement(:,:) = x_deformed(:,:) - X_undeformed(:,:)
c
          ! Compute material matrix "m_matrix" for linear elasticity
           ! The material parameters are stored in the list of material parameters "cm"
            E_modul = cm(1)
            Poisson = cm(2)
           fct_c_matrix = E_modul / ( (1.+Poisson)*(1.-2.*Poisson) )
           c11 = 1. - Poisson
           c44 = ( 1.-2.*Poisson )/2.
c
           m_matrix(1,1) = fct_c_matrix*c11
           m_matrix(2,2) = fct_c_matrix*c11
           m_matrix(3,3) = fct_c_matrix*c11
           m_matrix(4,4) = fct_c_matrix*c44
           m_matrix(5,5) = fct_c_matrix*c44
           m_matrix(6,6) = fct_c_matrix*c44
           m_matrix(1,2) = fct_c_matrix*Poisson
           m_matrix(1,3) = fct_c_matrix*Poisson
           m_matrix(2,1) = fct_c_matrix*Poisson
           m_matrix(2,3) = fct_c_matrix*Poisson
           m_matrix(3,1) = fct_c_matrix*Poisson
           m_matrix(3,2) = fct_c_matrix*Poisson
c
          ! Set up a 2x2x2 quadrature rule for full integration of the linear element
           sqrt13 = sqrt(1./3.)   
           qp_coord(:,1) = (/ -sqrt13,-sqrt13,-sqrt13 /)
           qp_coord(:,2) = (/ +sqrt13,-sqrt13,-sqrt13 /)
           qp_coord(:,3) = (/ +sqrt13,+sqrt13,-sqrt13 /)
           qp_coord(:,4) = (/ -sqrt13,+sqrt13,-sqrt13 /)
           qp_coord(:,5) = (/ -sqrt13,-sqrt13,+sqrt13 /)
           qp_coord(:,6) = (/ +sqrt13,-sqrt13,+sqrt13 /)
           qp_coord(:,7) = (/ +sqrt13,+sqrt13,+sqrt13 /)
           qp_coord(:,8) = (/ -sqrt13,+sqrt13,+sqrt13 /)
c
           qp_weight = 1.
c
          ! Initialise the first 6 history variables used for storing the Cauchy stress
          ! for this element "i_el" to zero as we use "add-in" concept below
           hsv(i_el,1:6) = 0
c
          ! Loop over all 8 quadrature points of this element "i_el" to collect their contributions
           do qp = 1, 8
              ! Compute the gradient of the shape function (8x3) with respect
              ! to the local coordinates for the current quadrature point "qp"
               dN_dXi(1,:) =
     &           (/ (-1./8.)*(1.-qp_coord(2,qp))*(1.-qp_coord(3,qp)),
     &              (-1./8.)*(1.-qp_coord(1,qp))*(1.-qp_coord(3,qp)),
     &              (-1./8.)*(1.-qp_coord(1,qp))*(1.-qp_coord(2,qp)) /)
               dN_dXi(2,:) = 
     &           (/ (+1./8.)*(1.-qp_coord(2,qp))*(1.-qp_coord(3,qp)),
     &              (-1./8.)*(1.+qp_coord(1,qp))*(1.-qp_coord(3,qp)),
     &              (-1./8.)*(1.+qp_coord(1,qp))*(1.-qp_coord(2,qp)) /)
               dN_dXi(3,:) = 
     &           (/ (+1./8.)*(1.+qp_coord(2,qp))*(1.-qp_coord(3,qp)),
     &              (+1./8.)*(1.+qp_coord(1,qp))*(1.-qp_coord(3,qp)),
     &              (-1./8.)*(1.+qp_coord(1,qp))*(1.+qp_coord(2,qp)) /)
               dN_dXi(4,:) = 
     &           (/ (-1./8.)*(1.+qp_coord(2,qp))*(1.-qp_coord(3,qp)),
     &              (+1./8.)*(1.-qp_coord(1,qp))*(1.-qp_coord(3,qp)),
     &              (-1./8.)*(1.-qp_coord(1,qp))*(1.+qp_coord(2,qp)) /)
               dN_dXi(5,:) = 
     &           (/ (-1./8.)*(1.-qp_coord(2,qp))*(1.+qp_coord(3,qp)),
     &              (-1./8.)*(1.-qp_coord(1,qp))*(1.+qp_coord(3,qp)),
     &              (+1./8.)*(1.-qp_coord(1,qp))*(1.-qp_coord(2,qp)) /)
               dN_dXi(6,:) = 
     &           (/ (+1./8.)*(1.-qp_coord(2,qp))*(1.+qp_coord(3,qp)),
     &              (-1./8.)*(1.+qp_coord(1,qp))*(1.+qp_coord(3,qp)),
     &              (+1./8.)*(1.+qp_coord(1,qp))*(1.-qp_coord(2,qp)) /)
               dN_dXi(7,:) = 
     &           (/ (+1./8.)*(1.+qp_coord(2,qp))*(1.+qp_coord(3,qp)),
     &              (+1./8.)*(1.+qp_coord(1,qp))*(1.+qp_coord(3,qp)),
     &              (+1./8.)*(1.+qp_coord(1,qp))*(1.+qp_coord(2,qp)) /)
               dN_dXi(8,:) = 
     &           (/ (-1./8.)*(1.+qp_coord(2,qp))*(1.+qp_coord(3,qp)), 
     &              (+1./8.)*(1.-qp_coord(1,qp))*(1.+qp_coord(3,qp)),
     &              (+1./8.)*(1.-qp_coord(1,qp))*(1.+qp_coord(2,qp)) /)
     
              ! Compute the Jacobian (3x3) for this element and its inverse (3x3)
               Jacobi = matmul( X_undeformed, dN_dXi )
                  
               detJ =   Jacobi(1,1)*Jacobi(2,2)*Jacobi(3,3)
     &                + Jacobi(1,2)*Jacobi(2,3)*Jacobi(3,1)
     &                + Jacobi(1,3)*Jacobi(2,1)*Jacobi(3,2)
     &                - Jacobi(3,1)*Jacobi(2,2)*Jacobi(1,3)
     &                - Jacobi(3,2)*Jacobi(2,3)*Jacobi(1,1)
     &                - Jacobi(3,3)*Jacobi(2,1)*Jacobi(1,2)
   
               inv_Jac(1,1) = (+1./detJ)*( Jacobi(2,2)*Jacobi(3,3)
     &                                    - Jacobi(3,2)*Jacobi(2,3) )
               inv_Jac(2,1) = (-1./detJ)*( Jacobi(2,1)*Jacobi(3,3)
     &                                     - Jacobi(3,1)*Jacobi(2,3) )
               inv_Jac(3,1) = (+1./detJ)*( Jacobi(2,1)*Jacobi(3,2)
     &                                    - Jacobi(3,1)*Jacobi(2,2) )
               inv_Jac(1,2) = (-1./detJ)*( Jacobi(1,2)*Jacobi(3,3)
     &                                     - Jacobi(3,2)*Jacobi(1,3) )
               inv_Jac(2,2) = (+1./detJ)*( Jacobi(1,1)*Jacobi(3,3)
     &                                    - Jacobi(3,1)*Jacobi(1,3) )
               inv_Jac(3,2) = (-1./detJ)*( Jacobi(1,1)*Jacobi(3,2)
     &                                     - Jacobi(3,1)*Jacobi(1,2) )
               inv_Jac(1,3) = (+1./detJ)*( Jacobi(1,2)*Jacobi(2,3)
     &                                    - Jacobi(2,2)*Jacobi(1,3) )
               inv_Jac(2,3) = (-1./detJ)*( Jacobi(1,1)*Jacobi(2,3)
     &                                     - Jacobi(2,1)*Jacobi(1,3) )
               inv_Jac(3,3) = (+1./detJ)*( Jacobi(1,1)*Jacobi(2,2)
     &                                    - Jacobi(2,1)*Jacobi(1,2) )
c
              ! Transform the shape gradient from the local coordinates Xi to the undeformed X
               dN_dX = matmul( dN_dXi, inv_Jac )
c
              ! Set up the B-matrix (6x24)
	            B_matrix(1,:) = (/ dN_dX(1,1), 0., 0.,
     &                            dN_dX(2,1), 0., 0.,
     &                            dN_dX(3,1), 0., 0., 
     &                            dN_dX(4,1), 0., 0., 
     &                            dN_dX(5,1), 0., 0., 
     &                            dN_dX(6,1), 0., 0., 
     &                            dN_dX(7,1), 0., 0., 
     &                            dN_dX(8,1), 0., 0. /)
               B_matrix(2,:) = (/ 0., dN_dX(1,2), 0., 
     &                            0., dN_dX(2,2), 0.,
     &                            0., dN_dX(3,2), 0., 
     &                            0., dN_dX(4,2), 0.,
     &                            0., dN_dX(5,2), 0., 
     &                            0., dN_dX(6,2), 0.,  
     &                            0., dN_dX(7,2), 0., 
     &                            0., dN_dX(8,2), 0. /) 
               B_matrix(3,:) = (/ 0., 0., dN_dX(1,3), 
     &                            0., 0., dN_dX(2,3),
     &                            0., 0., dN_dX(3,3), 
     &                            0., 0., dN_dX(4,3), 
     &                            0., 0., dN_dX(5,3), 
     &                            0., 0., dN_dX(6,3), 
     &                            0., 0., dN_dX(7,3), 
     &                            0., 0., dN_dX(8,3) /)
               B_matrix(4,:) = (/ dN_dX(1,2), dN_dX(1,1), 0., 
     &                            dN_dX(2,2), dN_dX(2,1), 0.,
     &                            dN_dX(3,2), dN_dX(3,1), 0.,
     &                            dN_dX(4,2), dN_dX(4,1), 0.,
     &                            dN_dX(5,2), dN_dX(5,1), 0.,
     &                            dN_dX(6,2), dN_dX(6,1), 0.,
     &                            dN_dX(7,2), dN_dX(7,1), 0.,
     &                            dN_dX(8,2), dN_dX(8,1), 0. /)
               B_matrix(5,:) = (/ 0., dN_dX(1,3), dN_dX(1,2), 
     &                            0., dN_dX(2,3), dN_dX(2,2), 
     &                            0., dN_dX(3,3), dN_dX(3,2), 
     &                            0., dN_dX(4,3), dN_dX(4,2), 
     &                            0., dN_dX(5,3), dN_dX(5,2), 
     &                            0., dN_dX(6,3), dN_dX(6,2), 
     &                            0., dN_dX(7,3), dN_dX(7,2), 
     &                            0., dN_dX(8,3), dN_dX(8,2) /)
               B_matrix(6,:) = (/ dN_dX(1,3), 0., dN_dX(1,1),
     &                            dN_dX(2,3), 0., dN_dX(2,1),
     &                            dN_dX(3,3), 0., dN_dX(3,1),
     &                            dN_dX(4,3), 0., dN_dX(4,1),
     &                            dN_dX(5,3), 0., dN_dX(5,1),
     &                            dN_dX(6,3), 0., dN_dX(6,1),
     &                            dN_dX(7,3), 0., dN_dX(7,1),
     &                            dN_dX(8,3), 0., dN_dX(8,1) /)
c
              ! Compute the strains (6x1) from the B-matrix and the displacements
               strain = matmul( B_matrix,
     &                          reshape (u_displacement, (/ 24 /) )
     &                        )
c
              ! Compute the stresses (6x1) from the material matrix and the strain
               stress = matmul( m_matrix, strain )
c
              ! Update the contribution to the force vector
               force(i_el,1:24) = force(i_el,1:24) 
     &                         + matmul( transpose(B_matrix), stress )
     &                           * detJ * qp_weight
c         
              ! Update the contribution to the stiffnes matrix if requested (istif=1)
               if (istif == 1) then 
                  stiff(i_el,1:24,1:24) = stiff(i_el,1:24,1:24)
     &             + matmul(
     &                        matmul( transpose(B_matrix), m_matrix ),
     &                        B_matrix
     &                     )
     &               * detJ * qp_weight
               endif
c           
              ! Store the computed stress into the history variables
              ! The element-averaged stress is stored, hence each of the 8 quadrature points contributes one-eighth to the average
               hsv(i_el,1:6) = hsv(i_el,1:6) + stress/8.
c
           enddo ! end loop over quadrature points
c
      enddo ! end loop over elements
      return
      end
```

## Some notes on proper simulation and solver settings for testing umats 

## Outline of the interface for umat and utan 

## Generalised interface for use of separate element subroutines

## Additional (extra) degrees of freedom xdofs

## Limitations
currently:
- no element deletion
- no remeshing
- no database cross section for section force: use nodfor instead
