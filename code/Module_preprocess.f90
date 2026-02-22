module Module_preprocess

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------


subroutine Conditions_initiales(i_init,gamma,rho_L,u_L,P_L,rho_R,u_R,P_R)
! ================================================
! Calcul des conditions initiales (Riemann problem)
! différentes configurations possibles
! ================================================
  implicit none
  integer,intent(in) 				 :: i_init			 				    ! Entier pour choisir condition
  double precision,intent(in)		 :: gamma								! Constante adiabatique
  double precision,intent(out)		 :: rho_L,u_L,P_L,rho_R,u_R,P_R			! Valeurs initiales
	
! Choix de la cdt initiale voulue : 1 : Sod subsonique, 2 : Sod Superso, 3 : vide, 4 : glissement, 5 : choc a Mach 3

! Validate input
  if (i_init < 1 .or. i_init > 5) then
     print *, "Error: Invalid i_init value. Must be between 1 and 5."
     stop
  end if
	
select case(i_init)
	case(1)
		! Test de Sod Subsonique 
		rho_L = 1.0D0; u_L = 0.0D0; P_L = 1.0D0  					! À gauche
		rho_R = 0.125D0; u_R = 0.0D0; P_R = 0.1D0 					! À droite
	case(2)
		! Test de Sod Supersonique 
		rho_L = 1.0D0; u_L = 0.75D0; P_L = 1.0D0  					! À gauche
		rho_R = 0.125D0; u_R = 0.0D0; P_R = 0.1D0 					! À droite
	case(3)
		! Apparition de vide
		rho_L = 1.0D0; u_L = -2.0D0; P_L = 0.4D0  					! À gauche
		rho_R = 1.0D0; u_R = 2.0D0; P_R = 0.4D0 					! À droite
	case(4)
		! Ligne de glissement stationnaire
		rho_L = 1.0D1; u_L = 0.0D0; P_L = 1.0D1  					! À gauche
		rho_R = 0.1D0; u_R = 0.0D0; P_R = 1.0D1 					! À droite
	case(5)
		! Choc stationnaire à Mach 3 
		rho_L = 1.0d0;        u_L = 3.5496478d0;   p_L = 1.0d0 		   ! À gauche
		rho_R = 3.85714285d0; u_R = 0.920279072d0; p_R = 10.33333333d0 ! À droite
end select


 end subroutine Conditions_initiales
 

 ! -----------------------------
 
 subroutine mesh(N,dx,x_i)
! ================================================
! Creation du maillage VF cell-centered
! ================================================
  implicit none
  integer,intent(in) 							  :: N			 ! grille de points
  double precision,intent(in) 					  :: dx			 ! pas de discretisation
  integer 			 							  :: i 			 ! entier pour calcul de boucles
  double precision,dimension(:),intent(out)		  :: x_i		 ! maillage
  
  ! Mesh creation
	do i = 1, N
		x_i(i) = (i - 0.5D0) * dx		! position de la cell
	enddo

 end subroutine mesh
 
 ! -----------------------------
 
subroutine Initialisation_variables(N,gamma,L,dx,rho_L,u_L,P_L,rho_R,u_R,P_R,x_i,w)
! ================================================
! Initialise les composantes du vecteur w
! w1 = rho, w2 = rho*u et w3 = rho*E
! ================================================
  implicit none
  integer,intent(in) 							  :: N			 								! grille de points
  double precision,intent(in) 					  :: rho_L,u_L,P_L,rho_R,u_R,P_R,L,dx,gamma		! parametres d'entrée de la subroutine
  double precision,dimension(:),intent(in)		  :: x_i										! mailles
  integer 			 							  :: i 			 								! entier pour calcul de boucles
  double precision								  :: x_c									    ! réel local
  double precision, dimension(:,-1:),intent(out)  :: w 	 										! Variables conservatives au temps t{n}
  
	x_c = 0.5d0*L  ! position de la discontinuite intitiale
	
	! BC à gauche pour ghost cell 0
	w(1,0) = rho_L  
	w(2,0) = rho_L * u_L
	w(3,0) = P_L / (gamma - 1.0D0) + 0.5D0 * rho_L * u_L**2
	! BC à gauche pour ghost cell -1
	w(:,-1) = w(:,0)
	
	! Calcul pour cell internes
	do i = 1, N
		if (x_i(i) < x_c) then
			w(1,i) = rho_L
			w(2,i) = rho_L * u_L
			w(3,i) = P_L / (gamma - 1.0D0) + 0.5D0 * rho_L * u_L**2
		else
			w(1,i) = rho_R
			w(2,i) = rho_R * u_R
			w(3,i) = P_R / (gamma - 1.0D0) + 0.5D0 * rho_R * u_R**2
		end if
		
	enddo
	! BC à droite pour ghost cell N+1
	w(1,N+1) = rho_R 
	w(2,N+1) = rho_R * u_R
	w(3,N+1) = P_R / (gamma - 1.0D0) + 0.5D0 * rho_R * u_R**2
	! BC à droite pour ghost cell N+2
	w(:,N+2) = w(:,N+1)
	
	
 end subroutine Initialisation_variables


end module Module_preprocess