program main
! ========================================================================================
! Programme de résolution du système d'Euler constitué des lois de conservation : 
! conservation de la masse, de la quantité de mvt et de l'energie
! en 1D instationnaire : d_t w + d_x F(w) = 0
! ou w = {rho,rho u, rho E} et F = {rho u, rho u²+P, (rho E+P) u}
! Condition initiale : w(x,0) = {w_R si x< x_d w_L sinon}
! Discrétisation spatiale choisit : n pts; x_i = (i-0.5)*h; i=1,n (cell-centered mesh)
!
! Domain cells :
!
!              |           |   u(i)    |           |
!              |  u(i-1)   |___________|           |
!              |___________|           |   u(i+1)  |
!              |           |           |___________|
!           ...|-----x-----|-----x-----|-----x-----|...
!              |    i-1    |     i     |    i+1    |
!              |-         +|-         +|-         +|
!            i-3/2       i-1/2       i+1/2       i+3/2

! Avec Reconstruction MUSCL : les variables d états ne sont plus cstes mais 
! évoluent de façon linéaires dans chaque mailles
! ========================================================================================


!!!!!!!!!!!!!!!!!!!!!!
! to create executable : make 
! to execute prog 	   : ./solve
! to clean 			   : make clean
!!!!!!!!!!!!!!!!!!!!!!
	
! Modules utilisés
  use Module_operateur
  use Module_lecture
  use Module_parametres

! Déclarations
  implicit none
  integer			 		  	 			 :: N              		! Nombre de points de grille (nombre impair comme ça le noeud central est en L/2)
  double precision							 :: delta_star			! Coefficient pour correction entropique de Harten
  double precision				 			 :: L               	! Longueur du domaine
  double precision			 	 			 :: gamma   			! Constante adiabatique
  double precision			 	 			 :: CFL					! Condition CFL
  double precision			 	 			 :: Tf      			! Temps final
  double precision			 	 			 :: beta     			! Coefficient pour limiter Chakravarthy
  double precision					         :: b,phi				! Parametre de compression, Parametre ordre de reconstruction
  double precision 			  	 			 :: dx,t1,t2			! réels double precision
  double precision, dimension(:),allocatable :: rho, u, P, E		! Variables des lois de conservations du fluide 
! Fin des declarations 
! =========================================================================================================================
! Main
! =========================================================================================================================

! Lecture of menu file
	call read_parameters(10,N,L,gamma,cfl,Tf,delta_star,beta,b,phi)   ! reads variables 

! Computation of space step dx 
	dx = L /(N-1)		! dx : Taille de la cell ; (N-1) : nombre de cells
	
! Allocation dynamique de mémoire
	allocate(rho(-1:N+2))   ! Masse volumique           [kg/m3]
	allocate(P(-1:N+2))     ! Pression 	                [Pa]
	allocate(u(-1:N+2))	    ! Vitesse 		            [m/s]
	allocate(E(-1:N+2))	    ! Energie totale spécifique [J/kg]

! Résolution des equations d'Euler avec MUSCL pour un fluide parfait unidimensionnel compressible
	call CPU_time(t1)
	call Solve_Eqn_Euler_MUSCL(N,delta_star,beta,L,dx,b,phi,gamma,CFL,Tf,rho,u,P,E)
	call CPU_time(t2)
	write(*,*)'Temps de resolution :',t2-t1 ! tps de calcul de la subroutine
	
! Création d'un script GNUPLOT : plotdata_exact.plt
	call create_gnuplot_script(Tf,"data_final")
	
! Trace avec GNUPLOT les graphes au temps final exact et approchés
	call system('gnuplot data_final/plotdata_exact.plt') 
	
! Ecris les parametres calcul à l'ecran
	call display_parameters(N,L,gamma,cfl,Tf,delta_star,beta,b,phi)  ! display on screen the variables selected 
	
! Desallocation dynamique de mémoire
	deallocate(rho, u, P, E)

!--------------fin
write(*,*)
write(*,*)'Execution avec succes !'
!===============
end program main
!===============