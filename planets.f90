program planets
  use rk_force2   !We are going to use our module in this problem
  implicit none
  !-------------------------------------------------------------------
  !Declaration of variables:
  integer, parameter       :: dp=kind(0.0_8)  !We are going to have double precision to our problem
  real(dp),allocatable     :: t(:),x(:,:)     !We are going to have a matrix for x to store for every dimension the required data
  real(dp),dimension(NEQ)  :: x0              !We have initial value as many as the number of equations
  real(dp)                 :: ti,tf           !Initial and final time
  integer                  :: Nt              !The number of places between ti and tf
  integer                  :: i,j             !Useful integers for do loops
  real(dp)                 :: dum1,dum2       !Useful dummy indexes during our program for k1 and k2 (we have k1 and k2 for Runge-Kutta orders)
  !-------------------------------------------------------------------
  !User Interface:
  print *, '#We are going to have 4th order Runge-Kutta method for number of ODEs:'
  print *, '#NEQ=',NEQ
  print *, '#Enter k1 and k2 for our problem constants:'
  read  *, k1,k2
  print *, '#k1,k2=',k1,k2
  print *, '#Enter the number of spaces between ti and tf:'
  read  *, Nt
  print *, '#Enter the initial time ti and the final time tf:'
  read  *, ti,tf
  print *, '#Enter the initial values of our problem:'
  read  *, x0  !It asks for values as many as the size of the array
  print *, '#Nt,ti,tf,x0=',Nt,ti,tf,x0
  !------------------------------------------------------------------
  !Array Allocation:
  ALLOCATE(t(Nt))
  ALLOCATE(x(Nt,NEQ))
  !------------------------------------------------------------------
  !Calculations:
  dum1=k1
  dum2=k2 !We keep the values of k1 and k2 to dummy indexes
  call rk  !We don't need to give input to our subroutine because we are going to have contains statement in a few lines
  !-------------------------------------------------------------------
  !Output file:
  k1=dum1
  k2=dum2   !Now we have the correct k1 and k2 for the calculation of the energy
  open(unit=11,file='voli.dat')
  do i=1,Nt
     write(11,*) t(i),x(i,:),energy(t(i),x(i,:))  !We have the ith line and all the collumns over there
  end do
  !-------------------------------------------------------------------
  close(11)
  !-------------------------------------------------------------------
contains
  !-------------------------------------------------------------------
  subroutine rk
    !-----------------------------------------------------------------
    !Declaration of local values:
    real(dp)   :: dt,ts,xs(NEQ) !We having the dt about space, ts that holds the values of t and also the vector of size NEQ to hold the values of the quantities because of dimensions
    integer    :: i,j           !Useful integers for do loops
    !-----------------------------------------------------------------
    !Calculations:
    dt=(tf-ti)/(Nt-1)
    t(1)=ti
    x(1,:)=x0 !We have our initial values for our problem
    ts=ti
    xs=x0     !We keep these values to make the next calculation with Runge-Kutta
    do i=2,Nt
       call rkstep(ts,xs,dt)
       t(i)=ts
       x(i,:)=xs
    end do
    !----------------------------------------------------------------
  end subroutine rk
  !-------------------------------------------------------------------
end program planets
!=====================================================================
