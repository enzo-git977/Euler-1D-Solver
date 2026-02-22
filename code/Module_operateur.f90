module Module_operateur

 ! modules utilisés
  use Module_sortie
  use Module_preprocess

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------

subroutine Solve_Eqn_Euler_MUSCL(N,i_init,i_fct,i_ord,i_cor,i_BC,delta_star,beta,L,dx,b,phi,gamma,CFL,Tf,rho,u,P,E)
! ======================================================================================= 
! Subroutine pour resoudre le systeme d'Euler fluide parfait 1D compressible
! avec MUSCL
! ======================================================================================= 
! Déclarations
  implicit none
  integer, intent(in) 			 			   	:: N         				! Nombre de points de grille
  integer, intent(in) 			 			   	:: i_init        		    ! Choix de la condition initiale voulue
  integer, intent(in) 			 			   	:: i_ord        		    ! Choix de l'ordre de precision spatial
  integer, intent(in) 			 			   	:: i_fct        		    ! Define what limiter_function to use
  integer, intent(in) 			 			   	:: i_BC					    ! Define what bc to use
  double precision, intent(in) 	 			   	:: L         				! Longueur du domaine
  double precision, intent(in) 	 			   	:: gamma      				! Constante adiabatique
  double precision, intent(in) 				   	:: CFL        				! Condition CFL
  double precision, intent(in) 	 			   	:: Tf   					! Temps final
  double precision, intent(in) 	 			   	:: beta  					! Coefficient pour limiter Chakravarthy
  double precision, intent(in)	 			   	:: dx						! pas de discrétisation spatiale
  integer,intent(in)							:: i_cor					! Activation de la correction
  double precision,intent(in)					:: delta_star				! Coefficient pour correction entropique de Harten
  double precision,intent(in)					:: b,phi				    ! Parametre de compression, Parametre ordre de reconstruction
  double precision, dimension(:),intent(inout) 	:: rho, u, P, E		    	! Caractéristiques physique du fluide
  double precision, dimension(:),allocatable	:: x_i		    			! Maillage
  double precision, dimension(:,:),allocatable 	:: w,w_L,w_R				! Vecteur w au temps t
  double precision, dimension(:,:),allocatable  :: w_new					! Vecteur w au temps t{n+1}
  double precision 				     		   	:: rho_L, u_L, P_L			! Cdt initiales a gauche
  double precision 					 		   	:: rho_R, u_R, p_R,dt,t		! Cdt initiales a droite
  double precision, dimension(:,:),allocatable 	:: F_L,F_R					! Vecteur Flux
  double precision, dimension(:,:),allocatable 	:: F_num					! Vecteur Flux numériques 
  double precision, dimension(:),allocatable 	:: Ma,s,Temp				! Mach number, entropy and Temperature 
  integer,parameter								:: N_sauv=10				! Nombre de sauvegardes
  integer 			 							:: sauv	    				! Compteur
  
  ! Fin des declarations
 !=======================================================================================
	! Allocations dynamiques des variables locales de la subroutine
	allocate(w(3,-1:N+2),w_new(3,-1:N+2))  ! intern cells : 1, N ; ghost cells : -1, 0 and N+1, N+2
	allocate(w_L(3,N+1),w_R(3,N+1))        ! Interpolation des etats gauches et droites du vecteur conservatif
	allocate(F_L(3,N+1), F_R(3,N+1))	   ! Flux 
	allocate(F_num(3,N+1))				   ! Flux numeriques calculés avec Roes scheme 
	allocate(x_i(N))					   ! Maillage 	
	allocate(Ma(N), s(N), Temp(N))		   ! Variables thermophysiques Mach number, entropy, temperature
	
	! Conditions initiales 
	call Conditions_initiales(i_init,gamma,rho_L,u_L,P_L,rho_R,u_R,P_R)
	
	t = 0.0D0 ! temps de la simulation [s]
	sauv = 0    ! Compteur pour sauvegarde
	
	! Mesh creation
	call mesh(N,dx,x_i)
	
	! Initialisation du vecteur des variables d'etats
	call Initialisation_variables(N,gamma,L,dx,rho_L,u_L,P_L,rho_R,u_R,P_R,x_i,w)
	
	! Write initial solution and plot it 
	call write_solution_initial_time(N,dx,gamma,w)
	
	! Boucle temporelle
	do while (t < Tf)
		call Compute_time_step(N,gamma,dx,CFL,w,dt) ! Calcul du pas de temps avec CFL
		if (t+dt.ge.Tf) then
			dt = Tf-t
		end if
		t =t + dt
		print*,'t = ', t, ' 	| dt =', dt
		
		! Interpolation (ordre 1 ou reconstruction MUSCL ordre 2 ou plus)
		call Interp(N,i_fct,i_ord,beta,b,phi,w,w_L,w_R)
		
		! Calcul des flux a gauche pour chaque cellule
		call Calcul_Flux(N,gamma,w_L,F_L) 
		
		! Calcul des flux a droite pour chaque cellule
		call Calcul_Flux(N,gamma,w_R,F_R) 
		
		! Calcul des Flux numérique aux interfaces avec le schéma de Roe
		call roe_flux(N,i_cor,delta_star, gamma,w_L,w_R,F_L,F_R,F_num)
		
		! Calcul des variables au temps t{n+1}
		call Update_variables(N,dx,dt,w ,F_num,w_new)
		
		! Appliquer les conditions aux bords
		call Boundary_conditions(N,gamma,i_BC,w_new)
		
		! Update du vecteur conservatif
		w = w_new
		
		!~! Save results at regular intervals 
		!~if ( floor(t*N_sauv/Tf)>floor((t-dt)*N_sauv/Tf) ) then
		!~	sauv=sauv+1
		!~	call Convert_conservative_to_primal(N,gamma,w_new,rho,u,P,E)  ! Conversion des variables d'etats conservatives en variables 'réels'
		!~	call save_results(N,"film",dx,t,rho,u,P,E)                  ! Sauvegarde des résultats
		!~	call create_gnuplot_script(t,"film")						! Creation des script
		!~	call system('gnuplot film/plot_script.plt')				    ! Plot automatique des images
		!~endif
		
	end do
	
	! sauvegarde
	call Convert_conservative_to_primal(N,gamma,w_new,rho,u,P,E)
	call Obtain_data(N, gamma, rho, P, u, Ma, s, Temp) ! Compute Mach number and entropy
	call save_results(N,"data_final",dx,Tf,rho,u,P,E, Ma,s,Temp) 
	
	! liberation de mémoire
	deallocate(w,w_L,w_R)
	deallocate(w_new)
	deallocate(F_L,F_R)
	deallocate(F_num)
	deallocate(x_i)
	deallocate(Ma,s,Temp)
    
end subroutine Solve_Eqn_Euler_MUSCL

! -----------------------------
 
 subroutine Compute_time_step(N,gamma,dx,CFL,w,dt)
! ==================================================
! Calcul le pas de temps pour assurer cdt stabilite
! ==================================================
	implicit none
    integer, intent(in) 						  :: N						! Nombre de points
    double precision,dimension(:,-1:),intent(in)  :: w						! Vecteur etat
	double precision, intent(in) 				  :: dx, gamma, CFL			! paramètres flottants d'entrées					
    double precision,dimension(-1:N+2)			  :: rho, u, E, P			! Variables locales
	double precision,dimension(N)			      :: c						! Vitesse du son locale
    integer 									  :: i						! Entier
    double precision 							  :: max_speed				! Calcul de la vit max
	double precision, intent(inout) 			  :: dt						! Pas de tps
! =========================================================================================================================
! Fin des declarations 

    ! Conversion des variables conservées en primitives
    call Convert_conservative_to_primal(N,gamma,w,rho,u,P,E)
	
	! Calcul de la vitesse du son
	do i = 1, N
        c(i) = dsqrt(gamma * P(i) / max(rho(i), 1.0d-12)) ! Célérité des ondes acoustiques [m/s]
    end do			
	
	!dt=CFL*dx/max(abs(u)+a);
	max_speed = maxval(dabs(u(1:N)) + c) ! on utilise que les mailles internes
	! Calcul du pas de temps
    dt = CFL * dx / max_speed
	
 end subroutine Compute_time_step
 
  ! -----------------------------
 
 subroutine Calcul_Flux(N,gamma,w_r,F)
! ================================================
! Calcul du vecteur flux F :
! F = {rho*u,rho*u² + P, (rho*E+P)*u}
! ================================================
  implicit none
  integer,intent(in) 							 :: N			 ! grille de points
  double precision, dimension(:,:),intent(in)    :: w_r 		 ! Composantes du vecteur w au temps t{n}
  double precision,intent(in) 					 :: gamma		 ! constante isentropique du fluide
  double precision,dimension(N+1)			     :: rho,u,P,E	 ! Variables locales
  integer 			 							 :: i 			 ! entier pour calcul de boucles
  double precision, dimension(:,:),intent(out)   :: F   		 ! Composantes du vecteur Flux
! =======================================================================================================
! Fin des declarations 
	
	! Initialisation
	rho = 0.0d0; u = 0.0d0; P = 0.0d0; E = 0.0d0
	
	! Conversion des variables conservées en primitives
	call Convert_conservative_to_primal_v2(N,gamma,w_r,rho,u,P,E)
	
	! Calcul des flux
	do i = 1, N+1 
		F(1,i) = rho(i)*u(i)
		F(2,i) = rho(i)*u(i)**2 + P(i)
		F(3,i) = u(i)*(rho(i)*E(i)+P(i)) 
	end do
	
 end subroutine Calcul_Flux
 
 ! -----------------------------
 
 subroutine Update_variables(N,dx,dt,w ,F_num,w_new)
! ================================================
! Calcul des variables au prochain temps 
! discrétisation avec Euler explicite ordre 1
! ================================================
  implicit none
  integer,intent(in) 							   :: N			! grille de points
  double precision,intent(in) 					   :: dx,dt		! paramètres d'entrées de la subroutine
  double precision, dimension(:,-1:),intent(in)    :: w 		! Composantes du vecteur w au temps t{n}
  double precision, dimension(:,:),intent(in)      :: F_num		! Flux numériques
  integer 			 							   :: i 		! entier pour calcul de boucles
  double precision, dimension(:,-1:),intent(inout) :: w_new 	! Composantes du vecteur w au temps t{n+1}
! =========================================================================================================================
! Fin des declarations 	

  ! Calcul des états au temps {t+1} (que les mailles internes)
	do i = 1, N 					
		w_new(:,i) = w(:,i) - (dt / dx) * (F_num(:,i+1) - F_num(:,i))
	end do

 end subroutine Update_variables

! -----------------------------
subroutine roe_flux(N,i_cor,delta_star, gamma,w_L,w_R,F_L,F_R,F_num)
! ======================================================================================= 
! Subroutine pour calculer les flux numériques a l'aide d'un schéma de Roe
! ======================================================================================= 
  ! Déclaration des variables
	implicit none
	integer,intent(in)								:: N									! Nombre de points
    double precision,intent(in) 					:: gamma								! Constante adiabatique
	double precision,dimension(:,:),intent(in) 		:: w_R,w_L								! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)		:: F_L,F_R								! Flux
	integer,intent(in)							    :: i_cor								! Activation de la correction
	double precision,intent(in)					    :: delta_star						    ! Coefficient pour correction entropique de Harten
	double precision, dimension(N+1)				:: rho_L,u_L,P_L,E_L			 		! Variables primitives a gauche
	double precision, dimension(N+1)				:: rho_R,u_R,P_R,E_R			 		! Variables primitives a droite
	double precision, dimension(N+1)				:: H_R,H_L								! Enthalpie totale
	integer 										:: j									! Entier pour boucle
    double precision 								:: u_moy, H_moy, c_moy,a				! Valeurs moyennes de Roe
    double precision,dimension(3)					:: wdif									! Vecteur différences du vecteur état
	double precision,dimension(3,3)					:: A_tilde								! Matrice A_tilde
	double precision,dimension(:,:), intent(out) 	:: F_num 								! Flux num à retourner
! =========================================================================================================================
! Fin des declarations 
	
    ! Initialisation du flux Roe
    F_num = 0.0d0
	
	! Calcul du left state
	call Convert_conservative_to_primal_v2(N,gamma,w_L,rho_L,u_L,P_L,E_L)
	
	! Calcul du right state
	call Convert_conservative_to_primal_v2(N,gamma,w_R,rho_R,u_R,P_R,E_R)
	
	! Calcul de l'enthalpie totale a gauche
	H_L = (gamma/(gamma-1))*(P_L/rho_L)+0.5d0*u_L**2
	
	! Calcul de l'enthalpie totale a droite
	H_R = (gamma/(gamma-1))*(P_R/rho_R)+0.5d0*u_R**2
	
    ! Boucle sur les interfaces
    do j = 1, N + 1
		
		! Calcul des moyennes de Roe
		a = dsqrt(rho_R(j)/ rho_L(j))									! rapport des masses volumiques
		u_moy = (u_R(j)+a*u_L(j))/(1.0d0+a)								! Moyenne de Roe de la vitesse matérielle
		H_moy = (a*H_R(j)+H_L(j))/(1.0d0+a)								! Moyenne de Roe de l'enthalpie	
		c_moy = dsqrt((gamma - 1.0d0) * (H_moy - 0.5d0 * u_moy**2))		! Moyenne de Roe de la vitesse de l'onde acoustique
		
		! Construction des matrices (voir excellent cours : P.S volpiani Roe's scheme)
		call Construct_matrix(i_cor, delta_star, gamma, u_moy, H_moy, c_moy, A_tilde)
		
		! Différence des états 
		wdif(:) = w_R(:,j) - w_L(:,j) 
		
        ! Roe : 0.5*[F(w_L)+F(w_R) - |A_tilde| * (W_R - W_L)]
		F_num(:,j) = 0.5d0*(F_L(:,j) + F_R(:,j)) - 0.5d0*matmul(A_tilde, wdif(:))
		
    end do
	
	
end subroutine roe_flux
 
 ! ----------------------------- 
 
  subroutine Boundary_conditions(N,gamma,i_BC,w_new)
! ================================================
! Calcul des BCs
! ================================================
  implicit none
  integer,intent(in) 							  :: N, i_BC		! grille de points, parametre de choix de bc
  double precision,intent(in) 					  :: gamma			! paramètre d'entrée de la subroutine
  double precision, dimension(:,-1:),intent(inout):: w_new 			! Variables des lois de conservations au temps t{n+1}
! =========================================================================================================================
! Fin des declarations  
	
  select case(i_BC)
	case(1)
		! Reflective BC : zero, gradient for density and pressure and dirichlet (0) for velocity -> reflective walls
		
		! Left BC
		! first ghost cell mirrors first internal cell
		w_new(1,0) = w_new(1,1)  ! zero gradient (neumann)
		w_new(2,0) = -w_new(2,1) ! so that at the interface u = 0     
		w_new(3,0) = w_new(3,1)  ! zero gradient (neumann)
		! Second ghost cell mirrors second internal cell
		w_new(1,-1)= w_new(1,2)
		w_new(2,-1)= - w_new(2,2)
		w_new(3,-1)= w_new(3,2)
		
		! Right BC
		! last ghost cell mirrors last internal cell
		w_new(1,N+1) = w_new(1,N)  ! zero gradient (neumann)
		w_new(2,N+1) = -w_new(2,N) ! so that at the interface u = 0     
		w_new(3,N+1) = w_new(3,N)  ! zero gradient (neumann)
		! Second to last cell mirrors second to last internal cell
		w_new(1,N+2)= w_new(1,N-1)
		w_new(2,N+2)= - w_new(2,N-1)
		w_new(3,N+2)= w_new(3,N-1)
		
	case(2)
		! Transmittive BC
		
		! Left BC
		w_new(:,0) = w_new(:,1)      ! zero gradient (neumann)
		w_new(:,-1) = w_new(:,2)     ! zero gradient (neumann)

	    ! Right BC
		w_new(:,N+1) = w_new(:,N)    ! zero gradient (neumann)
		w_new(:,N+2) = w_new(:,N-1)  ! zero gradient (neumann)
		
	 end select

 end subroutine Boundary_conditions
 
 ! ----------------------------- 
 
subroutine Convert_conservative_to_primal(N,gamma,w_new,rho,u,P,E)
! =========================================================
! Convertit les variables d'état en variables primaires
! =========================================================
  implicit none
  integer,intent(in)							:: N				! Number of cells
  double precision,intent(in) 					:: gamma			! paramètres d'entrées de la subroutine
  double precision, dimension(:,-1:),intent(in) :: w_new 			! Variables des lois de conservations au temps t{n+1}
  integer										:: j				! loop variables
  double precision, dimension(-1:),intent(out)  :: rho,u,P,E		! Variables primitives
! =========================================================================================================================
! Fin des declarations 

! Conversion des variables conservées en primitives
  do j = -1, N+2
     rho(j) = w_new(1, j)
     u(j)   = w_new(2, j) / max(rho(j), 1.0d-12) ! Sécurité contre la division par zéro ou densité négative
     ! E total par unité de masse
     E(j)   = w_new(3, j) / max(rho(j), 1.0d-12) 
     ! Calcul de la pression (Loi des gaz parfaits)
     P(j)   = (gamma - 1.0d0) * rho(j) * (E(j) - 0.5d0 * u(j)**2)
     
     ! Alerte si pression négative 
     if (P(j) < 0.0d0) then
		print*, "Erreur physique majeure : Pression negative j =", j, " P =", P(j)
		stop "Calcul interrompu"
	 endif
  end do

 end subroutine Convert_conservative_to_primal
 
 ! -----------------------------
 
 subroutine Convert_conservative_to_primal_v2(N,gamma,w_r,rho,u,P,E)
! =========================================================
! Convertit les variables d'état en variables primaires
! =========================================================
  implicit none
  integer,intent(in)							:: N				! Number of cells
  double precision,intent(in) 					:: gamma			! paramètres d'entrées de la subroutine
  double precision, dimension(:,:),intent(in)   :: w_r 			    ! Variables des lois de conservations au temps t{n+1}
  integer										:: j				! loop variables
  double precision, dimension(:),intent(out)    :: rho,u,P,E		! Variables primitives
! =========================================================================================================================
! Fin des declarations 
 
! Conversion des variables conservées en primitives
  do j = 1, N+1
     rho(j) = w_r(1, j)
     u(j)   = w_r(2, j) / max(rho(j), 1.0d-12) ! Sécurité contre la division par zéro ou densité négative
     ! E total par unité de masse
     E(j)   = w_r(3, j) / max(rho(j), 1.0d-12) 
     ! Calcul de la pression (Loi des gaz parfaits)
     P(j)   = (gamma - 1.0d0) * rho(j) * (E(j) - 0.5d0 * u(j)**2)
  end do

 end subroutine Convert_conservative_to_primal_v2
 
 ! -----------------------------
 
subroutine Construct_matrix(i_cor, delta_star, gamma, u_moy, H_moy, c_moy, A_tilde)
! =====================================================================================================
! Subroutine pour calculer les composantes des matrices A_tilde = P_tilde * Delta_tilde * P_tilde ^-1
! ======================================================================================================
	implicit none
    double precision,intent(in)	 				   :: u_moy, H_moy, c_moy, gamma        ! Valeurs moyennes de Roe
	integer,intent(in)							   :: i_cor								! Activation de la correction
	double precision,intent(in)					   :: delta_star						! Coefficient pour correction entropique de Harten
	double precision							   :: alph1, alph2, delta				! Variables locales
	double precision							   :: lambda_1,lambda_2,lambda_3		! Variables locales
	double precision							   :: psi_1,psi_2,psi_3					! Variables locales
    double precision,dimension(3,3) 			   :: Delta_tilde, Pinv, P_tilde		! Matrices pour Diagonalisation
	double precision,dimension(:,:),intent(inout)  :: A_tilde							! Matrice A_tilde a retourner
! =========================================================================================================================
! Fin des declarations   

    ! Variables auxiliaires pour P^{-1}
    alph1 = (gamma - 1.0d0) * u_moy**2 / (2.0d0 * c_moy**2)
    alph2 = (gamma - 1.0d0) / (c_moy**2)
	
    ! Matrice P^{-1}
    Pinv(1, 1) = 0.5d0 * (alph1 + u_moy / c_moy)
    Pinv(1, 2) = -0.5d0 * (alph2 * u_moy + 1.0d0 / c_moy)
    Pinv(1, 3) = alph2 / 2.0d0
    Pinv(2, 1) = 1.0d0 - alph1
    Pinv(2, 2) = alph2 * u_moy
    Pinv(2, 3) = -alph2
    Pinv(3, 1) = 0.5d0 * (alph1 - u_moy / c_moy)
    Pinv(3, 2) = -0.5d0 * (alph2 * u_moy - 1.0d0 / c_moy)
    Pinv(3, 3) = alph2 / 2.0d0

    ! Matrice P
    P_tilde(1, 1) = 1.0d0
	P_tilde(1, 2) = 1.0d0
    P_tilde(1, 3) = 1.0d0
    P_tilde(2, 1) = u_moy - c_moy
    P_tilde(2, 2) = u_moy
    P_tilde(2, 3) = u_moy + c_moy
    P_tilde(3, 1) = H_moy - c_moy * u_moy
    P_tilde(3, 2) = 0.5d0 * u_moy**2
    P_tilde(3, 3) = H_moy + c_moy * u_moy
	
	! Initialisation de la matrice Delta
	Delta_tilde = 0.0d0

	! Remplissage de la matrice Delta
	select case(i_cor)
		case(0)
			! Sans correction
			! Matrice Delta (valeurs propres absolues)
			Delta_tilde(1, 1) = dabs(u_moy - c_moy)
			Delta_tilde(2, 2) = dabs(u_moy)
			Delta_tilde(3, 3) = dabs(u_moy + c_moy)
		case(1)
			! Entropie seule
			delta = delta_star*(dabs(u_moy)+c_moy)
			lambda_1 = dabs(u_moy - c_moy) ! vitesse acoustique 1 
			lambda_3 = dabs(u_moy + c_moy) ! vitesse acoustique 2 
			if(lambda_1 <= delta) then
				psi_1 = 0.5d0*(lambda_1**2+delta**2)/delta
			else 
				psi_1 = lambda_1
			endif
			if(lambda_3 <= delta) then
				psi_3 = 0.5d0*(lambda_3**2+delta**2)/delta
			else 
				psi_3 = lambda_3
			endif
			! Matrice Delta 
			Delta_tilde(1, 1) = psi_1
			Delta_tilde(2, 2) = dabs(u_moy)
			Delta_tilde(3, 3) = psi_3
		case(2)
			! Complete
			delta = delta_star*(dabs(u_moy)+c_moy)
			lambda_1 = dabs(u_moy - c_moy) ! vitesse acoustique 1 
			lambda_2  = dabs(u_moy)		   ! vitesse materielle 
			lambda_3 = dabs(u_moy + c_moy) ! vitesse acoustique 2
			if(lambda_1 <= delta) then
				psi_1 = 0.5d0*(lambda_1**2+delta**2)/delta
			else 
				psi_1 = lambda_1
			endif
			if(lambda_2 <= delta) then
				psi_2 = 0.5d0*(lambda_2**2+delta**2)/delta
			else 
				psi_2 = lambda_2
			endif
			if(lambda_3 <= delta) then
				psi_3 = 0.5d0*(lambda_3**2+delta**2)/delta
			else 
				psi_3 = lambda_3
			endif
			! Matrice Delta 
			Delta_tilde(1, 1) = psi_1
			Delta_tilde(2, 2) = psi_2
			Delta_tilde(3, 3) = psi_3
	end select

    ! Calcul de la matrice |A_tilde|
	A_tilde = matmul(P_tilde, Delta_tilde)
    A_tilde = matmul(A_tilde, Pinv)
	
end subroutine Construct_matrix

 ! -----------------------------
 
 subroutine Interp(N,i_fct,i_ord,beta,b,phi,w,w_L,w_R)
! ======================================================================================
! Calcul les valeurs des variables aux interfaces
! avec reconstruction MUSCL : 
! ghost cell | Cell i-1 | Cell i | Cell i+1 | Cell i+2 | Cell i+3 | ghost cell
!     ...        ...       ...        ...        ...
!    w(i-2)     w(i-1)     w(i)     w(i+1)     w(i+2)
!
! Interfaces (flux calculés ici) : 
!       i-1/2     i+1/2      i+3/2      i+5/2
! Pour chaque interface on construit w_L et w_R qui sont ensuite donné au schéma de Roe
! ======================================================================================
  implicit none
  integer,intent(in) 							  :: N					 ! grille de points
  double precision, intent(in)					  :: beta				 ! Coefficient pour limiter Chakravarthy
  double precision, dimension(:,-1:),intent(in)   :: w 				     ! Composantes du vecteur w au temps t{n}
  integer,intent(in) 			 				  :: i_fct		         ! Define what limiter_function to use
  integer,intent(in) 			 				  :: i_ord		         ! Define order of spatial precision
  double precision,intent(in)					  :: b,phi				 ! Parametre de compression, Parametre ordre de reconstruction
  double precision							      :: eps,c1,c2			 ! Tolerance to avoid division by 0, coefficients
  integer										  :: i,j				 ! loop variables 
  double precision, dimension(:),allocatable	  :: psi1,psi2,psi3,psi4 ! Local variables to compute fct PSI(slope)
  double precision, dimension(:,:),allocatable	  :: r1,r2,r3,r4 	 	 ! Local slope of a cell
  double precision, dimension(:,:),allocatable	  :: dif	 			 ! Difference of neighbours cells
  double precision, dimension(:,:),allocatable	  :: l1,l2,l3,l4 	     ! final slope after limiter : X*PSI(Y/X) with Y/X : slope
  double precision, dimension(:,:),intent(inout)  :: w_L,w_R 		     ! Composantes du vecteur w au temps t{n+1}
! =========================================================================================================================
! Fin des declarations 
   
    select case(i_ord)
	case(1)
		! ordre 1 
		do i=1,N+1
			w_R(:,i) = w(:,i)
			w_L(:,i) = w(:,i-1)
		enddo
		
	case(2)
		! Reconstruction MUSCL (Ref : chap2 MNA ensma)
		! Allocation dynamique des variables locales de la subroutine
		allocate(psi1(3),psi2(3),psi3(3),psi4(3))
		allocate(dif(3,0:N+2), r1(3,N+1), r2(3,N+1), r3(3,N+1), r4(3,N+1))
		allocate(l1(3,N+1), l2(3,N+1), l3(3,N+1), l4(3,N+1))
		
		! Initialisation
		r1 = 0.0d0; r2 = 0.0d0; r3 = 0.0d0; r4 = 0.0d0
		l1 = 0.0d0; l2 = 0.0d0; l3 = 0.0d0; l4 = 0.0d0
		psi1 = 0.0d0; psi2 = 0.0d0; psi3 = 0.0d0; psi4 = 0.0d0
		w_L = 0.0d0
		w_R = 0.0d0
		eps = 1.0d-8
		
		! Compute difference of each cells
		do j = 0,N+2
			dif(:,j) = w(:,j)-w(:,j-1)
		enddo
		
		! Compute coefficients
		c1 = (1.0d0-phi)/4.0d0
		c2 = (1.0d0+phi)/4.0d0
		
		! Compute Left and right state for each interface
		do j = 1,N+1
			! Compute slopes
			r1(:,j) = dif(:,j-1)/(b*dif(:,j)+eps) ! eps : avoid 0 at denom
			r2(:,j) = dif(:,j)/(b*dif(:,j-1)+eps)
			r3(:,j) = dif(:,j)/(b*dif(:,j+1)+eps)
			r4(:,j) = dif(:,j+1)/(b*dif(:,j)+eps)
			! Compute limiter
			call limiter_function(i_fct,r1(:,j),beta,psi1)
			call limiter_function(i_fct,r2(:,j),beta,psi2)
			call limiter_function(i_fct,r3(:,j),beta,psi3)
			call limiter_function(i_fct,r4(:,j),beta,psi4)
			l1(:,j) = (b*dif(:,j))*psi1(:)
			l2(:,j) = (b*dif(:,j-1))*psi2(:)
			l3(:,j) = (b*dif(:,j+1))*psi3(:)
			l4(:,j) = (b*dif(:,j))*psi4(:)
			
			! Compute Left and right states of the interface
			w_L(:,j) = w(:,j-1) + c1*l1(:,j) + c2*l2(:,j)
			w_R(:,j) = w(:,j) - c2*l3(:,j) - c1*l4(:,j)
		enddo
		
		! Desallocation de memoire des variables propres a la subroutine
		deallocate(psi1,psi2,psi3,psi4)
		deallocate(dif,r1,r2,r3,r4,l1,l2,l3,l4)
		
    end select
  
 end subroutine Interp

! -----------------------------

subroutine limiter_function(i_fct,r,beta,psi)
! ================================================
! Limiteurs 
! ================================================
  implicit none
  integer, INTENT(IN) 		   		            :: i_fct         ! define what limiter_function to use
  double precision,INTENT(IN) 					:: beta    	     ! Parameter for Chakravarthy limiter
  double precision,dimension(3),INTENT(IN) 		:: r    		 ! Vector of the slope
  double precision 			   		            :: ii,jj,kk,eps	 ! Local values
  integer										:: i_boucle		 ! entier pour boucle
  double precision,dimension(3),intent(inout)   :: psi   		 ! Return value
! =========================================================================================================================
! Fin des declarations 

! Validate input
  if (i_fct < 1 .or. i_fct > 5) then
     print *, "Error: Invalid i_fct value. Must be between 1 and 5."
     stop
  end if
	
	
! Select the limiter function, 1: Minmod, 2 : VanLeer, 3 : VanAlbada, 4 : Superbee, 5 : Chakravarthy
select case(i_fct)
	case(1)
		! Minmod (most robust but most dissipative)
		do i_boucle = 1,3
			ii = min(1.0d0,r(i_boucle))
			psi(i_boucle) = max(0.0d0,ii)
		enddo
	case(2)
		! VanLeer
		eps = 1.0d-8
		do i_boucle = 1,3
			psi(i_boucle) = (r(i_boucle)+dabs(r(i_boucle)))/(1+r(i_boucle)+eps) ! avoid 0 division
		enddo
	case(3)
		! VanAlbada
		do i_boucle = 1,3
			ii = (r(i_boucle)+r(i_boucle)**2)/(1+r(i_boucle)**2)
			psi(i_boucle) = max(0.0d0,ii)
		enddo
	case(4)
		! Superbee (captures stiff front shock but less robust)
		do i_boucle = 1,3
			ii = min(2.0d0,r(i_boucle))
			jj = min(1.0d0,2.0d0*r(i_boucle))
			kk = max(ii,jj)
			psi(i_boucle) = max(0.0d0,kk)
		enddo
	case(5)
		! Chakravarthy
		do i_boucle = 1,3
			ii = min(beta,r(i_boucle))
			psi(i_boucle) = max(ii,0.0d0)
		enddo
end select

END subroutine limiter_function

! -----------------------------

subroutine Obtain_data(N, gamma, rho, P, u, Ma, s, Temp)
! ================================================
! Compute mach and entropy 
! ================================================
 implicit none
  integer, INTENT(IN) 		   		           :: N             ! Number of cells
  double precision, INTENT(IN) 		           :: gamma    	    ! Polytropic coeff
  double precision, dimension(-1:),intent(in)  :: rho,u,P  	    ! Variables primitives
  double precision							   :: c_sound	    ! Vitesse du son
  integer 			   		                   :: i			    ! Indice boucle
  double precision, dimension(:),intent(out)   :: Ma, s, Temp   ! Mach number, entropy and Temperature
! =========================================================================================================================
! Fin des declarations 
  
do i = 1, N
    ! Vitesse du son locale : c = sqrt(gamma * P / rho)
    c_sound = dsqrt(gamma * P(i) / rho(i))
    
    ! Nombre de Mach : M = |u| / c
    Ma(i) = dabs(u(i)) / c_sound
    
    ! Entropie (forme spécifique s = P / rho^gamma)
    s(i) = P(i) / (rho(i)**gamma)
	
	! Température adimensionnée : T = P / (rho * (gamma-1))
    Temp(i) = P(i) / (rho(i) * (gamma-1))
	
end do

END subroutine Obtain_data

end module Module_operateur
