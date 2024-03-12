!DEC$ IF .NOT. DEFINED (get_initialNodalCoords_1element_R102_f)
!DEC$ DEFINE get_initialNodalCoords_1element_R102_f
      function get_initialNodalCoords_1element_R102(i) result (X_coords)
c Return the initial/undeformed nodal coordinates
c for the current element "i" from the list of initial
c nodal coordinates "dm_rots" obtained from a pointer
c INPUT:
c - "i": index of current element (between lft and llt)
c OUTPUT:
c - "X_coords": array with initial nodal coordinates for element "i"
c               size(8,3): first index: node number 1,2,...,8; second index: x,y,z
c USAGE:
c - "X_coords = get_initialNodalCoords_1element_R102(i)"
c   integer :: i ! loop index "i" from loop over elements "do i=lft,llt"
c   real*8, dimension(8,3) :: X_coords
c NOTE:
c - compatible with Linux and Windows
c - tested under:
c       * ls-dyna_smp_d_R11_1_0_139588_winx64_ifort2017vs2017
c       * ls-dyna_smp_d_R11_1_0_x64_redhat65_ifort160_sse2
c - implemented for Release R10.2 and newer
c   (Windows: older releaes require different accessor method for "dm_rots";
c    Linux: for older releases such as R9 this approach is currently not available)
c - hardcoded for dimension=3D and 8-noded hexaeder element
c   (UEL currently only supports 3D hex)
c ToDo:
c - measure computation time and improve performance
c - maybe faster to return an array of all coordinates for block lft:llt
c
      use userinterface ! additionally necessary for use of getreal8ptr
c
      ! Include "nlqparm" to get the size "nlq"
       include 'nlqparm'
c
c Retrieve the index map for the current element
C_TASKCOMMON (aux33loc)
      common/aux33loc/ix1(nlq),ix2(nlq),ix3(nlq),ix4(nlq),ix5(nlq),
     & ix6(nlq),ix7(nlq),ix8(nlq)
c
      ! define the dimension (3D: dim=3; 2D: dim=2)
       integer, parameter :: dim=3
      ! index of current element
       integer, intent(in) :: i
      ! array of nodal coordinates   
       real*8, dimension(8,dim) :: X_coords
c
      ! retrieve the initial nodal coordinates "dm_rots"
      ! from module "mod_nodes"
      ! "dm_rots" is a vector:
      !  shape: number of elements times dimension
      !  entries: for instance: ( x1, y1, z1, x2, y2, z2, ... ) for dimension=3D   
       real*8, CONTIGUOUS, POINTER :: dm_rots(:) ! [POINTER_CONTIGUOUS]
c$omp threadprivate (/aux33loc/)
       call getreal8ptr('dm_rots',dm_rots)
c
      ! We directly write the entries from "dm_rots" into "X_coords" to
      ! avoid working on the large list "dm_rots" with large reshapes
       X_coords(:,:) = reshape(
     &                   dm_rots( 
     &                           (/
     &                              ((ix1(i)-1)*dim+1):(ix1(i)*dim),
     &                              ((ix2(i)-1)*dim+1):(ix2(i)*dim),
     &                              ((ix3(i)-1)*dim+1):(ix3(i)*dim),
     &                              ((ix4(i)-1)*dim+1):(ix4(i)*dim),
     &                              ((ix5(i)-1)*dim+1):(ix5(i)*dim),
     &                              ((ix6(i)-1)*dim+1):(ix6(i)*dim),
     &                              ((ix7(i)-1)*dim+1):(ix7(i)*dim),
     &                              ((ix8(i)-1)*dim+1):(ix8(i)*dim)
     &                            /)
     &                           ), (/8,dim/), order=(/2,1/)
     &                         )
c
      end function
!DEC$ ENDIF
