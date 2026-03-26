module Module_lecture

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------

subroutine read_parameters(nfich)
!=====================================================================
!  Fonction : Lecture du fichier menu pour les parametres de calcul
!=====================================================================
! Parametres globaux
use Module_parametres

! Déclarations
  implicit none
  integer, intent(in) :: nfich ! Numero du fichier pour ouverture
! =========================================================================================================================
! Fin des declarations
  
! Main
WRITE(*,*) 'Debut lecture du fichier menu :'
open(nfich,file='menu', form='formatted')
read(nfich,'(17/)') ! skip header
read(nfich,*) N
read(nfich,'(/)')
read(nfich,*) L
read(nfich,'(/)')
read(nfich,*) gamma
read(nfich,'(/)')
read(nfich,*) cfl
read(nfich,'(/)')
read(nfich,*) Tf
read(nfich,'(/)')
read(nfich,*) i_init
read(nfich,'(/)')
read(nfich,*) i_BC
read(nfich,'(/)')
read(nfich,*) i_sc
read(nfich,'(/)')
read(nfich,*) i_ord
read(nfich,'(/)')
read(nfich,*) i_t
read(nfich,'(/)')
read(nfich,*) i_cor
read(nfich,'(/)')
read(nfich,*) delta_star
read(nfich,'(/)')
read(nfich,*) phi
if(phi==3) then
	phi = 1.0d0/3.0d0
endif
read(nfich,'(/)')
read(nfich,*) b
read(nfich,'(/)')
read(nfich,*) i_fct
read(nfich,'(/)')
read(nfich,*) beta

close(nfich)


end subroutine read_parameters


! -----------------------------


subroutine display_parameters()

!=====================================================================
!  Fonction : Affichage des parametres de calcul selectionnes
!=====================================================================
! Parametres globaux
use Module_parametres
 ! =========================================================================================================================
! Fin des declarations 
  
! Main
write(*,*) '--- Parametres du calcul  :  ---'
write(*,*) 'N = ', N
write(*, '(A, F8.3)')  'L     = ', L
write(*, '(A, F8.3)')  'gamma = ', gamma
write(*, '(A, F8.3)')  'cfl   = ', cfl
write(*, '(A, F8.3)')  'Tf    = ', Tf

select case(i_init)
	case(1)
	    write(*,*) 'Cdts initiales : Sod Subsonique'
	case(2)
		write(*,*) 'Cdts initiales : Sod Supersonique'
	case(3)
		write(*,*) 'Cdts initiales : Apparition de vide'
	case(4)
		write(*,*) 'Cdts initiales : Ligne de glissement stationnaire'
	case(5)
		write(*,*) 'Cdts initiales : Choc stationnaire à Mach 3 '
	case(6)
		write(*,*) 'Cdts initiales : Blast wave (forte pression) '
	case(7)
		write(*,*) 'Cdts initiales : Collision de chocs '
end select

select case(i_BC)
	case(1)
	    write(*,*) 'BC : Reflective wall'
	case(2)
		write(*,*) 'BC : Transmissive wall'
end select

write(*,*)'Choice of the flux scheme :'
select case(i_sc)
	case(1)
		write(*,*) 'Roe'
	case(2)
	    write(*,*) 'HLL (Davis estimation)'
	case(3)
		write(*,*) 'HLL (Mean estimation)'
	case(4)
		write(*,*) 'HLLE (Einfeldt estimation)'
	case(5)
		write(*,*) 'HLLC-ANRS'
	case(6)
		write(*,*) 'HLLC robuste'
end select


select case(i_ord)
	case(1)
	    write(*,*) 'ordre 1 en espace'
	case(2)
		if(dabs(phi) - 1.0d0/3.0d0 <= 1.0d-6) then
			write(*,*) 'ordre 3 en espace'
		else
			write(*,*) 'ordre 2 en espace'
		endif
		write(*, '(A, F8.3)')  'PHI   = ', phi
		write(*, '(A, F8.3)')  'parametre de compression b = ', b
		select case(i_fct)
			case(1)
				write(*,*) 'Choix du limiter : Minmod'
			case(2)
				write(*,*) 'Choix du limiter : VanLeer'
			case(3)
				write(*,*) 'Choix du limiter : VanAlbada'
			case(4)
				write(*,*) 'Choix du limiter : Superbee'
			case(5)
				write(*,*) 'Choix du limiter : Chakravarthy'
				write(*, '(A, F8.3)')  'beta     = ', beta
		end select
end select

select case(i_t)
	case(1)
	    write(*,*) 'Time integration : Euler ordre 1'
		write(*, *)
	case(2)
		write(*,*) 'Time integration : RK2'
		write(*, *)
end select

select case(i_cor)
	case(0)
	    write(*,*) 'Sans correction entropique'
		write(*, *)
	case(1)
		write(*,*) 'Correction entropique'
		write(*,'(A, F8.3)') 'Coefficient delta_star=',delta_star
		write(*, *)
	case(2)
		write(*,*) 'Correction entropique complete'
		write(*, '(A, F8.3)') 'Coefficient delta_star=',delta_star
		write(*, *)
end select


end subroutine display_parameters

end module Module_lecture