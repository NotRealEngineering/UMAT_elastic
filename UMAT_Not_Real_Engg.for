!Not Real Engineering - All Rights Reserved You may not use, 
!                       distribute and modify this code without 
!                       the written permission from Not Real Engineering.
!**********************************************************************************************************
!
      module NumKind
!
!**********************************************************************************************************
        implicit none
        integer (kind(1)), parameter :: ikind = kind(1), 
     &                                  rkind = kind(0.D0), 
     &                                  lkind = kind(.true.)
!       
      end module Numkind
!**********************************************************************************************************

!****************************************************************************************
      subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl,
     &                ddsddt, drplde, drpldt, stran, dstran, time, 
     &                dtime, temp, dtemp, predef, dpred, cmname, 
     &                ndi, nshr, ntens, nstatv, props, nprops, 
     &                coords, drot, pnewdt, celent, dfgrd0, dfgrd1, 
     &                noel, npt, layer, kspt, kstep, kinc)
!****************************************************************************************

        use NumKind
		INCLUDE 'ABA_PARAM.INC'

!**********************************************************************************************************
!	INFO:
      ! USER MATERIAL SUBROUTINE
	  ! Gets called for each integration point at each increment
	  ! Must define: "ddsdde" and "stress"
	  ! STRESS(NTENS): This array is passed in as the stress tensor 
	  !                at the beginning of the increment and must be 
	  !                updatedin this routine to be the stress tensor 
	  !                at the end of the increment.
	  ! DDSDDE(NTENS,NTENS): Jacobian matrix of the constitutive model
	  
	  ! NTENS: Size of the stress or strain component array (NDI + NSHR).
	  ! NDI: Number of direct stress components at this point. 
	  ! NSHR: Number of engineering shear stress components at this point. 
	  ! DSTRAN(NTENS): Array of strain increments.	  
!**********************************************************************************************************
          
        integer (ikind) :: ndi, nshr, ntens, nstatv, nprops
        integer (ikind) :: noel, npt, layer
        integer (ikind) :: kspt, kstep(4), kinc
        
        real (rkind) :: stress(ntens), statev(nstatv), 
     &                  ddsdde(ntens,ntens), ddsddt(ntens),
     &                  drplde(ntens), stran(ntens), dstran(ntens), 
     &                  props(nprops)
        real (rkind) :: coords(3), drot(3), time(2), dfgrd0(3,3), 
     &                  dfgrd1(3,3)
        real (rkind) :: predef, dpred, sse, spd, scd, rpl, drpldt, 
     &                  dtime, temp, dtemp, pnewdt, celent,
     &                  ELAM, EG2, EG, E, Nu 	
        character*8 cmname
		integer (ikind) :: k1, k2 

	   ! Material properties defined in input file        
		E=props(1) ! Young's modulus
		Nu=props(2)  ! Poisson's ratio
		
		EG=E/(2.0*(1.0+Nu))
		EG2=EG*2.0
		ELAM=EG2*Nu/(1.0-2.0*Nu)

!****************************************************************************************
!		Stiffness tensor ddsdde
!****************************************************************************************
		do k1=1, ntens
          do k2=1, ntens
			ddsdde(k2, k1)=0.0
          end do
		end do

		do k1=1, ndi
		  do k2=1, ndi
			ddsdde(k2, k1)=ELAM
		  end do
			ddsdde(k1, k1)=EG2+ELAM
		end do 
	
		do k1=ndi+1, ntens
			ddsdde(k1, k1)=EG
		end do
!****************************************************************************************
		
!****************************************************************************************
!	Calculate Stresses
!****************************************************************************************
       do k1=1, ntens
        do k2=1, ntens
         stress(k2)=stress(k2)+ddsdde(k2, k1)*dstran(k1)
        end do
       end do 
!****************************************************************************************
        return

      end subroutine umat	  
!****************************************************************************************
!Not Real Engineering - All Rights Reserved You may not use, 
!                       distribute and modify this code without 
!                       the written permission from Not Real Engineering.