c UEL-helper, e.g. get_initialNodalCoords      
      include "../UEL_helper_Fortran_LS-Dyna/UEL_helper.f"
c
!...
c
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
c     compute force and stiffness for user defined solid 101
c
      ! Load the UEL_helper module to retrieve the undeformed nodal coordinates below via "get_initialNodalCoords_1element_R102(*)"
       use UEL_helper
      ! default LS-Dyna include and declarations      
       include 'nlqparm'
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
       dimension hsv(nlq,nhsv),ihsv(nlq,nhsv)
       real*8 :: cm(lmc)
       dimension cmtrx(nlq,21)
       dimension dener(nlq)
c
      ! custom declarations for this user-element
       integer, parameter :: dimen = 3
       integer, parameter :: n_nodes = 8
c
       ! Current element ID
        integer :: ele
       ! Iterator for nodes
        integer :: iNode
       ! Undeformed coordinates of element's nodes
        real*8 :: xUndef(n_nodes, dimen)
       ! Deformed coordinates of element's nodes
        real*8 :: xDef(n_nodes, dimen)
       ! Total displacements of element's nodes
        real :: uDisp(n_nodes,dimen)
       ! Element force vector
        real*8 :: forceVector(ndtot)
       ! Element stiffness matrix
        real*8 :: stiffMatrix(ndtot,ndtot)
       ! List of all dofs (displacements and xdofs)
        real*8 :: pe(ndtot)
       ! Number of extra degrees of freedom
        integer :: nxdof
c
      ! Determine total number of dofs per node (integer-division intended)
       ndofs_per_node = ndtot/8  ! only valid for 3D hexahedral element
      ! Determine number of xdofs per node
       nxdof = ndofs_per_node - 3 ! only valid for 3D hexahedral element
c
c Loop over all elements in the list [lft,llt]
       ELEMENTS: do ele = lft, llt
         ! Fill matrix of deformed coordinates for current element 'ele'
          xDef = reshape(
     &              (/
     &                x1(ele), y1(ele), z1(ele),
     &                x2(ele), y2(ele), z2(ele),
     &                x3(ele), y3(ele), z3(ele),
     &                x4(ele), y4(ele), z4(ele),
     &                x5(ele), y5(ele), z5(ele),
     &                x6(ele), y6(ele), z6(ele),
     &                x7(ele), y7(ele), z7(ele),
     &                x8(ele), y8(ele), z8(ele)
     &               /), (/8,dimen/), order=(/2,1/)
     &           )
         ! Retrieve matrix of undeformed coordinates for current element 'ele'
          xUndef = get_initialNodalCoords_1element_R102( ele )
         ! Compute total displacements for current element 'ele'
          uDisp = xDef - xUndef
         ! Combine all dofs in the list 'pe' (displacements and possibly xdofs)
          do iNode=1,n_nodes
            if ( nxdof==0 ) then
                pe( 1+(iNode-1)*ndofs_per_node : iNode*ndofs_per_node )
     &            = (/ uDisp(iNode,:) /)
            else           
                pe( 1+(iNode-1)*ndofs_per_node : iNode*ndofs_per_node )
     &            = (/ uDisp(iNode,:) , xdof(ele,iNode,1:nxdof) /)
            endif
          enddo 
         ! Initialise output variables of UEL to zero, because of internal "AddIn"
          forceVector = 0.0;
          stiffMatrix = 0.0;
c          
         ! Call UEL-routine for current element 'ele'
          if ( nxdof==0 ) then
            call Q1X( xUndef, pe, cm,
     &                hsv(ele,1:nhsv), nhsv, ndtot, istif, ele,
     &                forceVector, stiffMatrix, hsv(ele,1:nhsv) )
          elseif ( nxdof==1 ) then
            call Q1X_GX( xUndef, pe, cm 
     &                   hsv(ele,1:nhsv), nhsv, ndtot, istif, ele,
     &                   forceVector, stiffMatrix, hsv(ele,1:nhsv) )
          else
            write(*,*) "usld_e101<<
     & Currently only nxdof={0,1} supported"
            call cstop("ERROR TERMINATION")
          endif
c
         ! Add internal force contribution of current element to 'force'
          force(ele,1:ndtot) = force(ele,1:ndtot)
     &                         + forceVector
c
         ! If stiffness matrix is requested (istiff==1), add it to 'stiff'
          if( istif == 1 ) then
            ! Check stiffMatrix for NaN
            ! and replace by 1e6 to avoid error "nan found in imasem_unsym"
            ! @todo-optimize Is 1e6 a good choice?
             do ii=1,ndtot
               do jj=1,ndtot
                  if (isNan( stiffMatrix(ii,jj) )) then
                     stiffMatrix(ii,jj) = 1e6
                  endif
               enddo
             enddo
            ! Add stiffness contribution of current element to 'stiff'
            ! @note Purposefully using "AddTo", to keep possible previous values of 'stiff' by LS-Dyna
             stiff(ele,1:ndtot,1:ndtot) = stiff(ele,1:ndtot,1:ndtot)
     &                                    + stiffMatrix
          endif   
c         
       end do ELEMENTS  
      return
      end
