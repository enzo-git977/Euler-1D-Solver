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

subroutine write_solution_initial_time(N,dx,gamma,w)
!=====================================================================
! Ecris la solution au temps t=0 pour ensuite la plot avec GNUPLOT
!=====================================================================
	implicit none
	integer,intent(in) 					          :: N			     		! grille de points
	double precision,intent(in) 		          :: dx,gamma				! parametres d'entrée de la subroutine
	double precision, dimension(:,-1:),intent(in) :: w 					    ! Variables des lois de conservations 
    integer 			 				          :: i 			 	 	    ! entier pour calcul de boucles
	double precision							  :: c_sound			    ! Vitesse du son
    double precision, dimension(N)		          :: rho,u,P,E,Ma,s,Temp	! Caractéristiques physique du fluide
	
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
		write(11, *) (i-0.5d0)*dx, rho(i), Ma(i), P(i), Temp(i), E(i), s(i)
	end do
	write(11,  *) 1.0d0, rho(N), Ma(N), P(N), Temp(N), E(N), s(N)
    close(11)
	
	! lance gnuplot
	call system('gnuplot data_final/plotdata_t=0.plt')
	
end subroutine write_solution_initial_time

! ----------------------------------------------------------------------


subroutine save_results(N,Save_choice,dx,Tf,rho,u,P,E,Ma,s,Temp)
! ================================================
! write results in a file in the data folder
! ================================================
  implicit none
  integer,intent(in) 							:: N			 					! grille de points
  double precision, dimension(-1:),intent(in)   :: rho, u, P, E						! Caractéristiques physique du fluide
  double precision, dimension(:),intent(in)     :: Ma, s, Temp					   	! Caractéristiques physique du fluide
  double precision,intent(in) 					:: dx,Tf							! parametres d'entrée de la subroutine
  character(len=*),intent(in) 				    :: Save_choice						! Pick either data_final or film
  character(len=50) 							:: filename							! name of the file
  integer 			 							:: i,ierr 			 				! entier pour calcul de boucles
  
	! Choose filename precision depending on Save_choice
    if (trim(Save_choice) == 'data_final') then
		write(filename,'(A, F4.2, A)') trim(Save_choice)//'/results_at_t=', Tf, 's.dat' ! Create filename with the time value included (rounded to 2 decimal places)
    else if (trim(Save_choice) == 'film') then
		write(filename,'(A, F0.2, A)') trim(Save_choice)//'/results_at_t=', Tf, 's.dat' ! Create filename with the time value included (rounded to 5 decimal places)
    else
        print *, "Unknown Save_choice: ", trim(Save_choice)
        stop
    end if
    
    ! Open the file with the generated filename
    open(unit=10, file=trim(filename), status="unknown", form="formatted", iostat=ierr)
    if (ierr /= 0) then
        print*, 'Error opening file:', trim(filename)
        stop
    end if
	
	! writing
	do i = 1, N-1
		write(10, *) (i-0.5d0)*dx, rho(i), Ma(i), P(i), Temp(i), E(i), s(i)
	end do
	write(10, *) 1.0d0, rho(N), Ma(N), P(N), Temp(N), E(N), s(N)
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
	
  write(output_png,'(A,F4.2,A)') trim(Save_choice)//'/results_at_t=', Tf, 's.png' ! maybe use G0 instaead of F6.2

  ! Gnuplot script file name
  gnuplot_file = trim(Save_choice)//'/plot_script.plt'

  ! Create filename with the time value included (rounded to 2 decimal places)
  write(filename,'(A, F4.2, A)') trim(Save_choice)//'/results_at_t=', Tf, 's.dat'

  ! Create title for the plots
  write(title_plot, '("t=", F4.2, "s")') Tf

  ! Open the Gnuplot script file for writing
  open(newunit=iunit, file=gnuplot_file, status="replace")

  ! Write Gnuplot commands to the file for multiplot layout
  write(iunit, *) "set terminal pngcairo enhanced font 'arial,10' size 1200, 800"
  write(iunit, *) "set output '" // trim(output_png) // "'"
  write(iunit, *) "set title 'Results ", trim(title_plot), "'"
  write(iunit, *) "set multiplot layout 2,2 title 'Variables at " // trim(title_plot) // "'"

  ! Plot Density (rho)
  write(iunit, *) "set xlabel 'x'"
  write(iunit, *) "set ylabel 'Density (rho)'"
  write(iunit, *) "plot '", trim(filename), "' using 1:2 with lines lw 2 title 'Density (rho)'"

  ! Plot Velocity (u)
  write(iunit, *) "set xlabel 'x'"
  write(iunit, *) "set ylabel 'Velocity (u)'"
  write(iunit, *) "plot '", trim(filename), "' using 1:3 with lines lw 2 title 'Velocity (u)'"

  ! Plot Pressure (P)
  write(iunit, *) "set xlabel 'x'"
  write(iunit, *) "set ylabel 'Pressure (P)'"
  write(iunit, *) "plot '", trim(filename), "' using 1:4 with lines lw 2 title 'Pressure (P)'"

  ! Plot Energy (E)
  write(iunit, *) "set xlabel 'x'"
  write(iunit, *) "set ylabel 'Energy (E)'"
  write(iunit, *) "plot '", trim(filename), "' using 1:5 with lines lw 2 title 'Energy (E)'"

  ! End multiplot mode
  write(iunit, *) "unset multiplot"

  ! Close the Gnuplot script file
  close(iunit)

  print *, "Gnuplot script 'plot_script.plt' has been generated"

end subroutine create_gnuplot_script


  ! ----------------------------------------------------------------------
end module Module_sortie


