!DEC$ IF .NOT. DEFINED (UEL_helper_F)
!DEC$ DEFINE UEL_helper_F
! -----------MODULE user-defined element helper---------------------------
      module UEL_helper
c 
      contains
c
!      ------BEGIN FUNCTIONS-------------------------------------
		include './get_initialNodalCoords_1element_R102.f'
c
      end module UEL_helper
!DEC$ ENDIF
