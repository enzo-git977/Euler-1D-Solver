module Module_sortie

contains

! ----------------------------------------------------------------------
! ECRITURE VECTEURS ET MATRICES
! ----------------------------------------------------------------------

subroutine ecrit_vec(x,n) ! Ecrit le vecteur coeff par coeff
    implicit none
    real(kind=kind(0.d0)),dimension(:)			            :: x    
    integer                                    				:: i,n
    do i=1,n
        write(*,fmt = '(I4,1X,F7.3)')i,x(i)
    enddo
	write(*,*)
end subroutine ecrit_vec

! ----------------------------------------------------------------------

subroutine ecrit_mat(A,n) ! Ecrit la matrice à l'ecran sous un format 'reconnaissable'
    implicit none
    real(kind=kind(0.d0)),dimension(:,:)			     :: A    
    integer                                  			 :: i,j,n
    do i=1,3
        do j=1,n
            write(*,'(F8.4, "  ")', advance = 'no')A(i,j)
        enddo
      	write(*,*)
    enddo
	write(*,*)
end subroutine ecrit_mat

! ----------------------------------------------------------------------

subroutine write_solution_initial_time(dx,w)
!=====================================================================
! Ecris la solution au temps t=0 pour ensuite la plot avec GNUPLOT
!=====================================================================
! Parametres globaux 
use Module_parametres, only : N, gamma

 ! Déclaration des variables
	implicit none
	double precision,intent(in) 		          :: dx						! parametres d'entrée de la subroutine
	double precision, dimension(:,-1:),intent(in) :: w 					    ! Variables des lois de conservations 
    integer 			 				          :: i 			 	 	    ! entier pour calcul de boucles
	double precision							  :: c_sound			    ! Vitesse du son
    double precision, dimension(N)		          :: rho,u,P,E,Ma,s,Temp	! Caractéristiques physique du fluide
! =========================================================================================================================
! Fin des declarations 	
	
	! Conservative to primal conversion
	do i=1,N
		rho(i) = w(1, i)
		u(i) = w(2, i) / max(rho(i), 1.0d-12)
		E(i) = w(3, i) / max(rho(i), 1.0d-12)
		P(i) = (gamma - 1.0d0) * rho(i) * (E(i) - 0.5d0 * u(i)**2)
		c_sound = dsqrt(gamma * P(i) / rho(i))
		Ma(i) = dabs(u(i)) / c_sound
		s(i) = P(i) / (rho(i)**gamma)
        Temp(i) = P(i) / (rho(i) * (gamma-1))
	enddo
	
	! ecriture dans fichier pour plot
	open(11, file="data_final/results_at_0.00s.dat", status="unknown")
	do i = 1, N-1
		write(11, *) (i-0.5d0)*dx, rho(i), Ma(i), P(i), Temp(i), u(i), s(i)
	end do
	write(11,  *) 1.0d0, rho(N), Ma(N), P(N), Temp(N), u(N), s(N)
    close(11)
	
	! lance gnuplot
	call system('gnuplot data_final/plotdata_t=0.plt')
	
end subroutine write_solution_initial_time

! ----------------------------------------------------------------------


subroutine save_results(Save_choice,dx,Tf,rho,u,P,E,Ma,s,Temp)
! ================================================
! write results in a file in the data folder
! ================================================
! Parametres globaux 
use Module_parametres, only : N

 ! Déclaration des variables
  implicit none
  double precision, dimension(-1:),intent(in)   :: rho, u, P, E						! Caractéristiques physique du fluide
  double precision, dimension(:),intent(in)     :: Ma, s, Temp					   	! Caractéristiques physique du fluide
  double precision,intent(in) 					:: dx,Tf							! parametres d'entrée de la subroutine
  character(len=*),intent(in) 				    :: Save_choice						! Pick either data_final or film
  character(len=10) 							:: time_str							! char for time
  character(len=50) 							:: filename							! name of the file
  integer 			 							:: i,ierr 			 				! entier pour calcul de boucles
! =========================================================================================================================
! Fin des declarations 	  
	
	! Convert time into character
	write(time_str, '(F5.3)') Tf

	! Construct filename without spaces
	filename = trim(Save_choice) // '/results_at_t=' // trim(adjustl(time_str)) // 's.dat'
	
	! Verification 
	print *, "Save file : ", trim(filename)
    
    ! Open the file with the generated filename
    open(unit=10, file=trim(filename), status="unknown", form="formatted", iostat=ierr)
    if (ierr /= 0) then
        print*, 'Error opening file:', trim(filename)
        stop
    end if
	
	! writing
	do i = 1, N-1
		write(10, *) (i-0.5d0)*dx, rho(i), Ma(i), P(i), Temp(i), u(i), s(i)
	end do
	write(10, *) 1.0d0, rho(N), Ma(N), P(N), Temp(N), u(N), s(N)
	close(10)
  
  end subroutine save_results
  
  ! ----------------------------------------------------------------------
  
subroutine create_gnuplot_script(Tf,Save_choice)
! =============================================================
! Creates a gnuplot script in the data folder
! This script plots the results for t=Tf in a multiplot layout
! =============================================================
  implicit none
  character(len=100)          :: gnuplot_file        ! Name of the gnuplot script
  integer                     :: iunit               ! For handling errors
  double precision,intent(in) :: Tf                  ! Simulation time
  character(len=50)           :: filename            ! Name for the data file
  character(len=30)           :: title_plot          ! Title for the plots
  character(len=*),intent(in) :: Save_choice		 ! Pick either data_final or film
  character(len=100) 		  :: output_png		     ! Name of output
	
	! name of the script
    gnuplot_file = trim(Save_choice) // '/plotdata_exact.plt'

	! open file
    open(newunit=iunit, file=trim(gnuplot_file), status="replace")

    ! --- Configuration Terminal et Output ---
    write(iunit, '(A)') "set terminal pngcairo enhanced font 'arial,10' size 1200, 800"
    write(iunit, '(A,A,A)') "set output '", trim(Save_choice), "/Solution_Exact_et_approche.png'"
    write(iunit, *) ""

    ! --- Définition du temps final pour le sprintf ---
    write(iunit, '(A,F8.4)') "t_final = ", Tf
    write(iunit, '(A,A,A)') "filename_final = sprintf('", trim(Save_choice), "/results_at_t=%.3fs.dat', t_final)"
    write(iunit, *) ""

    ! --- Multiplot ---
    write(iunit, '(A)') 'set multiplot layout 3,2 title "Exact vs Approximate Riemann solver"'
    write(iunit, *) ""

    ! --- 1. Density ---
    write(iunit, '(A)') "set xlabel 'x (m)'"
    write(iunit, '(A)') "set ylabel 'Density '"
    write(iunit, '(A)') "set title 'Density Comparison'"
    write(iunit, '(A)') "plot 'data_final/Exact_Sod.dat' using 1:2 with lines lw 2 title 'Exact', \"
    write(iunit, '(A)') "     filename_final using 1:2 with lines lw 2 title sprintf('t=%.3fs', t_final)"

    ! --- 2. Mach ---
    write(iunit, '(A)') "set xlabel 'x (m)'"
    write(iunit, '(A)') "set ylabel 'Mach number (-)'"
    write(iunit, '(A)') "set title 'Mach Comparison'"
    write(iunit, '(A)') "plot 'data_final/Exact_Sod.dat' using 1:3 with lines lw 2 title 'Exact', \"
    write(iunit, '(A)') "     filename_final using 1:3 with lines lw 2 title sprintf('t=%.3fs', t_final)"

    ! --- 3. Pressure ---
    write(iunit, '(A)') "set xlabel 'x (m)'"
    write(iunit, '(A)') "set ylabel 'Pressure '"
    write(iunit, '(A)') "set title 'Pressure Comparison'"
    write(iunit, '(A)') "plot 'data_final/Exact_Sod.dat' using 1:4 with lines lw 2 title 'Exact', \"
    write(iunit, '(A)') "     filename_final using 1:4 with lines lw 2 title sprintf('t=%.3fs', t_final)"

    ! --- 4. Temperature ---
    write(iunit, '(A)') "set xlabel 'x (m)'"
    write(iunit, '(A)') "set ylabel 'Temperature (-)'"
    write(iunit, '(A)') "set title 'T adimensionne Comparison'"
    write(iunit, '(A)') "plot 'data_final/Exact_Sod.dat' using 1:5 with lines lw 2 title 'Exact', \"
    write(iunit, '(A)') "     filename_final using 1:5 with lines lw 2 title sprintf('t=%.3fs', t_final)"

    ! --- 5. Velocity ---
    write(iunit, '(A)') "set xlabel 'x (m)'"
    write(iunit, '(A)') "set ylabel 'velocity '"
    write(iunit, '(A)') "set title 'velocity Comparison'"
    write(iunit, '(A)') "plot 'data_final/Exact_Sod.dat' using 1:6 with lines lw 2 title 'Exact', \"
    write(iunit, '(A)') "     filename_final using 1:6 with lines lw 2 title sprintf('t=%.3fs', t_final)"

    ! --- 6. Entropy ---
    write(iunit, '(A)') "set xlabel 'x (m)'"
    write(iunit, '(A)') "set ylabel 's '"
    write(iunit, '(A)') "set title 'Entropy Comparison'"
    write(iunit, '(A)') "plot 'data_final/Exact_Sod.dat' using 1:7 with lines lw 2 title 'Exact', \"
    write(iunit, '(A)') "     filename_final using 1:7 with lines lw 2 title sprintf('t=%.3fs', t_final)"

    write(iunit, '(A)') "unset multiplot"

    close(iunit)

  print *, "Gnuplot script 'plotdata_exact.plt' has been generated"

end subroutine create_gnuplot_script


  ! ----------------------------------------------------------------------
end module Module_sortie


