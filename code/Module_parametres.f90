module Module_parametres

implicit none

    ! --- Choix des schémas ---
	integer	:: i_init              ! choix de la cdt initiale voulue 
    integer :: i_t    ! Schéma en temps (1: Euler, 2: RK2)
    integer :: i_sc   ! Schéma spatial (1: Roe, 2: HLL, 3: HLLC)
    integer :: i_ord  ! Ordre spatial (1: Ordre 1, 2: MUSCL)
    integer :: i_fct  ! Limiteur de pente (1: Minmod, 2: Superbee, etc.)
    integer :: i_BC   ! Conditions aux limites (1: Transmissif, 2: Réflexion)
    integer :: i_cor  ! Correction d'entropie (0: Non, 1: Oui)
	
contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------
	
end module Module_parametres
