module Module_solve

 ! modules utilisés
  use Module_sortie
  use Module_preprocess
  use Module_riemann

contains

! ----------------------------------------------------------------------
! METHODES
! ----------------------------------------------------------------------

subroutine Solve_Eqn_Euler_MUSCL(dx,rho,u,P,E)
! ======================================================================================= 
! Subroutine pour resoudre le systeme d'Euler fluide parfait 1D compressible
! avec MUSCL
! ======================================================================================= 
! Parametres globaux 
use Module_parametres, only : N, i_t, Tf

! Déclarations
  implicit none
  double precision, intent(in) 	 			   	:: dx								    ! Double from the main
  double precision, dimension(:),intent(inout) 	:: rho, u, P, E		    				! Caractéristiques physique du fluide
  double precision, dimension(:),allocatable	:: x_i		    						! Maillage
  double precision, dimension(:,:),allocatable 	:: w,w_L,w_R							! Vecteur w au temps t
  double precision, dimension(:,:),allocatable  :: w_mid								! Vecteur w au temps t{n+1/2}
  double precision, dimension(:,:),allocatable  :: w_new								! Vecteur w au temps t{n+1}
  double precision 				     		   	:: rho_L, u_L, P_L						! Cdt initiales a gauche
  double precision 					 		   	:: rho_R, u_R, p_R,dt,t					! Cdt initiales a droite
  double precision, dimension(:,:),allocatable 	:: F_L,F_R								! Vecteur Flux
  double precision, dimension(:,:),allocatable 	:: F_num								! Vecteur Flux numériques 
  double precision, dimension(:),allocatable 	:: Ma,s,Temp							! Mach number, entropy and Temperature 
  integer,parameter								:: N_sauv=10							! Nombre de sauvegardes
  integer 			 							:: sauv	    							! Compteur
  
  ! Fin des declarations
 !=======================================================================================
	! Allocations dynamiques des variables locales de la subroutine
	allocate(w(3,-1:N+2),w_new(3,-1:N+2))  ! intern cells : 1, N ; ghost cells : -1, 0 and N+1, N+2
	allocate(w_L(3,N+1),w_R(3,N+1))        ! Interpolation des etats gauches et droites du vecteur conservatif
	allocate(F_L(3,N+1), F_R(3,N+1))	   ! Flux 
	allocate(F_num(3,N+1))				   ! Flux numeriques calculés avec Roes scheme 
	allocate(x_i(N))					   ! Maillage 	
	allocate(Ma(N), s(N), Temp(N))		   ! Variables thermophysiques Mach number, entropy, temperature
	
	! Pour RK2
	if(i_t == 2) then
		allocate(w_mid(3,-1:N+2))
		w_mid = 0.0d0 ! initialisation du vecteur
	endif
	
	! Conditions initiales 
	call Conditions_initiales(rho_L,u_L,P_L,rho_R,u_R,P_R)
	
	t = 0.0D0 ! temps de la simulation [s]
	sauv = 0    ! Compteur pour sauvegarde
	
	! Mesh creation
	call mesh(dx,x_i)
	
	! Initialisation du vecteur des variables d'etats
	call Initialisation_variables(dx,rho_L,u_L,P_L,rho_R,u_R,P_R,x_i,w)
	
	! Write initial solution and plot it 
	call write_solution_initial_time(dx,w)
	
	! Boucle temporelle
	do while (t < Tf)
		call Compute_time_step(dx,rho,u,E,P,w,dt) ! Calcul du pas de temps avec CFL
		if (t+dt.ge.Tf) then
			dt = Tf-t
		end if
		t =t + dt
		print*,'t = ', t, ' 	| dt =', dt
		
		! Interpolation (ordre 1 ou reconstruction MUSCL ordre 2 ou plus)
		call Interp(w,w_L,w_R)
		
		! Calcul des flux a gauche pour chaque cellule
		call Calcul_Flux(w_L,F_L) 
		
		! Calcul des flux a droite pour chaque cellule
		call Calcul_Flux(w_R,F_R) 
		
		! Calcul des Flux numérique aux interfaces avec un schéma numerique (1 : Roe, 2:  HLL, 3 : HLLC)
		call Num_flux(w_L,w_R,F_L,F_R,F_num)
		
		! Calcul des variables au temps t{n+1}
		call Update_variables(dx,dt,w,w_L,w_R,w_mid,F_L,F_R,F_num,w_new)
		
		! Appliquer les conditions aux bords
		call Boundary_conditions(w_new)
		
		! Update du vecteur conservatif
		w = w_new
		
		!~! Save results at regular intervals 
		!~if ( floor(t*N_sauv/Tf)>floor((t-dt)*N_sauv/Tf) ) then
		!~	sauv=sauv+1
		!~	call Convert_conservative_to_primal(w_new,rho,u,P,E)  ! Conversion des variables d'etats conservatives en variables 'réels'
		!~	call save_results("film",dx,t,rho,u,P,EMa,s,Temp)           ! Sauvegarde des résultats
		!~	call create_gnuplot_script(t,"film")						! Creation des script
		!~	call system('gnuplot film/plot_script.plt')				    ! Plot automatique des images
		!~endif
		
	end do
	
	! sauvegarde
	call Convert_conservative_to_primal(w_new,rho,u,P,E)
	call Obtain_data(rho, P, u, Ma, s, Temp) ! Compute Mach number and entropy
	call save_results("data_final",dx,Tf,rho,u,P,E, Ma,s,Temp) 
	
	! liberation de mémoire
	deallocate(w,w_L,w_R)
	deallocate(w_new)
	deallocate(F_L,F_R)
	deallocate(F_num)
	deallocate(x_i)
	deallocate(Ma,s,Temp)
	! Pour RK2
	if(i_t == 2) then
		deallocate(w_mid)
	endif
    
end subroutine Solve_Eqn_Euler_MUSCL

! -----------------------------
 
 subroutine Compute_time_step(dx,rho,u,E,P,w,dt)
! ==================================================
! Calcul le pas de temps pour assurer cdt stabilite
! ==================================================
! Parametres globaux 
use Module_parametres, only : N, gamma, CFL

! Déclarations
	implicit none
    double precision,dimension(:,-1:),intent(in)  :: w						! Vecteur etat
	double precision, intent(in) 				  :: dx						! paramètre flottant d'entrée					
    double precision,dimension(-1:),intent(inout) :: rho, u, E, P			! Variables primaires
	double precision,dimension(N)			      :: c						! Vitesse du son locale
    integer 									  :: i						! Entier
    double precision 							  :: max_speed				! Calcul de la vit max
	double precision, intent(inout) 			  :: dt						! Pas de tps
! =========================================================================================================================
! Fin des declarations 

    ! Conversion des variables conservées en primitives
    call Convert_conservative_to_primal(w,rho,u,P,E)
	
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
 
 subroutine Calcul_Flux(w_r,F)
! ================================================
! Calcul du vecteur flux F :
! F = {rho*u,rho*u² + P, (rho*E+P)*u}
! ================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

! Déclarations
  implicit none
  double precision, dimension(:,:),intent(in)    :: w_r 		 ! Composantes du vecteur w au temps t{n}
  double precision 					 			 :: P			 ! Pression
  integer 			 							 :: i 			 ! entier pour calcul de boucles
  double precision, dimension(:,:),intent(out)   :: F   		 ! Composantes du vecteur Flux
! =======================================================================================================
! Fin des declarations 
	
	! Calcul des flux
	do i = 1, N+1 
		F(1,i) = w_r(2,i)
		F(2,i) = 0.5d0*(3.0d0-gamma)*(w_r(2,i)**2 / w_r(1,i)) + (gamma-1.0d0)*(w_r(3,i))
		P = (w_r(3,i)-0.5d0*(w_r(2,i)**2 / w_r(1,i)))*(gamma-1.0d0)
		F(3,i) = (w_r(2,i)*w_r(3,i)) / w_r(1,i) + P*w_r(2,i)/w_r(1,i)
	end do
	
 end subroutine Calcul_Flux
 
 ! -----------------------------
 
 subroutine Update_variables(dx,dt,w,w_L,w_R,w_mid,F_L,F_R,F_num,w_new)
! ================================================
! Calcul des variables au prochain temps 
! discrétisation avec Euler explicite ordre 1
! ================================================
! Parametres globaux 
use Module_parametres, only : i_t, N

! Déclarations
  implicit none
  double precision,intent(in) 					   :: dx,dt								 ! paramètres d'entrées de la subroutine
  double precision, dimension(:,-1:),intent(in)    :: w 								 ! Composantes du vecteur w au temps t{n}
  double precision, dimension(:,:),intent(inout)   :: F_num								 ! Flux numériques
  double precision, dimension(:,-1:),intent(inout) :: w_mid								 ! Vecteur w au temps t{n+1/2}
  double precision, dimension(:,:),intent(inout)   :: w_L,w_R   						 ! Interpolation du vecteur w
  double precision, dimension(:,:),intent(inout)   :: F_L,F_R   						 ! Composantes du vecteur Flux
  integer 			 							   :: i 								 ! entier pour calcul de boucles
  double precision, dimension(:,-1:),intent(inout) :: w_new 							 ! Composantes du vecteur w au temps t{n+1}
! =========================================================================================================================
! Fin des declarations

! Validate input
  if (i_t < 1 .or. i_t > 2) then
     print *, "Error: Invalid i_t value. Must be between 1 and 2."
     stop
  end if

  select case(i_t)
	case(1)
		! Euler ordre 1 
		! Calcul des états au temps {t+1} (que les mailles internes)
		do i = 1, N 					
			w_new(:,i) = w(:,i) - (dt / dx) * (F_num(:,i+1) - F_num(:,i))
		end do
		
	case(2)
		! RK2
		! Calcul des états au temps {tn+1/2} (que les mailles internes)
		do i = 1, N 					
			w_mid(:,i) = w(:,i) - (dt / (2.0d0*dx)) * (F_num(:,i+1) - F_num(:,i))
		end do
		
		! Appliquer les conditions aux bords
		call Boundary_conditions(w_mid)
		
		! Interpolation (ordre 1 ou reconstruction MUSCL ordre 2 ou plus)
		call Interp(w_mid,w_L,w_R)
			
		! Calcul des flux a gauche pour chaque cellule
		call Calcul_Flux(w_L,F_L) 
			
		! Calcul des flux a droite pour chaque cellule
		call Calcul_Flux(w_R,F_R) 
			
		! Calcul des Flux numérique aux interfaces avec un schéma numerique (1 : Roe, 2:  HLL, 3 : HLLC)
		call Num_flux(w_L,w_R,F_L,F_R,F_num)
		
		! Calcul des états au temps {t+1} (que les mailles internes)
		do i = 1, N 					
			w_new(:,i) = w(:,i) - (dt / dx) * (F_num(:,i+1) - F_num(:,i))
		end do 	
  end select

 end subroutine Update_variables

! -----------------------------

subroutine Num_flux(w_L,w_R,F_L,F_R,F_num)
! ======================================================================================= 
! Subroutine pour calculer les flux numériques 
! 1) a l'aide d'un schéma de Roe
! 2) Schema HLL Davis
! 3) Schema HLL Roe
! 4) Schema HLLE
! 5) Schema HLLC-ANRS (Adaptive Non-iterative Riemann Solver (ANRS))
! ===================================================================================================================================================== 
! Parametres globaux 
use Module_parametres, only : i_sc

  ! Déclaration des variables
	implicit none
	double precision,dimension(:,:),intent(in) 	:: w_R,w_L		! États conservés (rho, rho*u, rho*E) a doite et a gauche de l'interface
	double precision,dimension(:,:),intent(in)	:: F_L,F_R		! Flux
	double precision,dimension(:,:),intent(out) :: F_num 		! Flux num à retourner
! =======================================================================================================================================================
! Fin des declarations 

    ! Initialisation du flux numerique
    F_num = 0.0d0
	
	select case(i_sc)
		case(1)
			! Roe SCHEME
			call ROE_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(2)
			! HLL-Davis SCHEME
			call HLL_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(3)
			! HLLR SCHEME (Hll-Roe)
			call HLLR_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(4)
			! HLLE SCHEME (Hll-Einfield)
			call HLLE_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(5)
			! HLLC-ANRS SCHEME
			call HLLC_ANRS_flux(w_L,w_R,F_L,F_R,F_num)
			
		case(6)
			! HLLC SCHEME
			call HLLC_flux(w_L,w_R,F_L,F_R,F_num)
	end select
	
end subroutine Num_flux
 
 ! ----------------------------- 
 
  subroutine Boundary_conditions(w_new)
! ================================================
! Calcul des BCs
! ================================================
! Parametres globaux 
use Module_parametres, only : i_BC, N

 ! Déclaration des variables
  implicit none
  double precision, dimension(:,-1:),intent(inout):: w_new 			! Variables des lois de conservations au temps t{n+1}
! =========================================================================================================================
! Fin des declarations 

 ! Validate input
  if (i_BC < 1 .or. i_BC > 2) then
     print *, "Error: Invalid i_BC value. Must be between 1 and 2."
     stop
  end if 
	
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
 
subroutine Convert_conservative_to_primal(w_new,rho,u,P,E)
! =========================================================
! Convertit les variables d'état en variables primaires
! =========================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

 ! Déclaration des variables
  implicit none
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
 
 
 subroutine Interp(w,w_L,w_R)
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
! Parametres globaux 
use Module_parametres, only : i_ord, N, b, phi

 ! Déclaration des variables
  implicit none
  double precision, dimension(:,-1:),intent(in)   :: w 				     ! Composantes du vecteur w au temps t{n}
  double precision							      :: eps,c1,c2			 ! Tolerance to avoid division by 0, coefficients
  integer										  :: i,j				 ! loop variables 
  double precision, dimension(3)				  :: psi1,psi2,psi3,psi4 ! Local variables to compute fct PSI(slope)
  double precision, dimension(3,N+1)			  :: r1,r2,r3,r4 	 	 ! Local slope of a cell
  double precision, dimension(3,0:N+2)			  :: dif	 			 ! Difference of neighbours cells
  double precision, dimension(3,N+1)			  :: l1,l2,l3,l4 	     ! final slope after limiter : X*PSI(Y/X) with Y/X : slope
  double precision, dimension(:,:),intent(inout)  :: w_L,w_R 		     ! Composantes du vecteur w au temps t{n+1}
! =========================================================================================================================
! Fin des declarations 

   ! Validate input
  if (i_ord < 1 .or. i_ord > 2) then
     print *, "Error: Invalid i_ord value. Must be between 1 and 2."
     stop
  end if
   
    select case(i_ord)
	case(1)
		! ordre 1 
		do i=1,N+1
			w_R(:,i) = w(:,i)
			w_L(:,i) = w(:,i-1)
		enddo
		
	case(2)
		! Reconstruction MUSCL (Ref : chap2 MNA ensma)
		
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
			call limiter_function(r1(:,j),psi1)
			call limiter_function(r2(:,j),psi2)
			call limiter_function(r3(:,j),psi3)
			call limiter_function(r4(:,j),psi4)
			l1(:,j) = (b*dif(:,j))*psi1(:)
			l2(:,j) = (b*dif(:,j-1))*psi2(:)
			l3(:,j) = (b*dif(:,j+1))*psi3(:)
			l4(:,j) = (b*dif(:,j))*psi4(:)
			
			! Compute Left and right states of the interface
			w_L(:,j) = w(:,j-1) + c1*l1(:,j) + c2*l2(:,j)
			w_R(:,j) = w(:,j) - c2*l3(:,j) - c1*l4(:,j)
		enddo
		
    end select
  
 end subroutine Interp

! -----------------------------

subroutine limiter_function(r,psi)
! ================================================
! Limiteurs 
! ================================================
! Parametres globaux 
use Module_parametres, only : i_fct, beta

! Déclaration des variables
  implicit none
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

subroutine Obtain_data(rho, P, u, Ma, s, Temp)
! ================================================
! Compute mach and entropy 
! ================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

 ! Déclaration des variables
 implicit none
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
	Ma(i) = u(i) / c_sound
    
    ! Entropie (forme spécifique s = P / rho^gamma)
    s(i) = P(i) / (rho(i)**gamma)
	
	! Température adimensionnée : T = P / (rho * (gamma-1))
    Temp(i) = P(i) / (rho(i) * (gamma-1))
	
end do

END subroutine Obtain_data

end module Module_solve
