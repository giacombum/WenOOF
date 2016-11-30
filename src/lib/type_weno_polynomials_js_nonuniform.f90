module type_weno_polynomials_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Lagrange polynomials for Jiang-Shu WENO schemes on non uniform grids.
!<
!< @note The provided polynomials implement the Lagrange polynomials defined in *Efficient implementation of WENO schemes to
!< nonuniform meshes*, Nelida Črnjarić-Žic, Senka Maćešić, Bojan Crnković, ANNALI DELL'UNIVERSITA' DI FERRARA, 2007, vol. 57,
!< issue 2, pp. 199--215, doi:10.1007/s11565-007-0013-1
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P
use type_weno_polynomials
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: weno_polynomials_js_nonuniform, associate_WENO_polynomials_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(weno_polynomials) :: weno_polynomials_js_nonuniform
  !< Lagrange polynomials for non Jiang-Shu WENO schemes on non uniform grids object.
  !<
  !< @note The provided polynomials implement the Lagrange polynomials defined in *Efficient implementation of WENO schemes to
  !< nonuniform meshes*, Nelida Črnjarić-Žic, Senka Maćešić, Bojan Crnković, ANNALI DELL'UNIVERSITA' DI FERRARA, 2007, vol. 57,
  !< issue 2, pp. 199--215, doi:10.1007/s11565-007-0013-1
  private
  contains
    procedure, pass(self), public :: destroy
    procedure, pass(self), public :: create
    procedure, nopass,     public :: description
    procedure, pass(self), public :: compute
endtype weno_polynomials_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public, non TBP
  function associate_WENO_polynomials_js_nonuniform(polyn_input) result(polyn_pointer)
    !< Check the type of the polynomials passed as input and return non uniform Jiang-Shu polynomials associated to the polynomials.
    class(weno_polynomials), intent(in),    target  :: polyn_input   !< Input optimal weights.
    class(weno_polynomials_js_nonuniform),  pointer :: polyn_pointer !< Non uniform Jiang Shu optimal weights.

    select type(polyn_input)
      type is(weno_polynomials_js_nonuniform)
        polyn_pointer => polyn_input
      class default
        write(stderr, '(A)')'error: wrong polynomials type chosen'
        stop
    end select
  end function associate_WENO_polynomials_js_nonuniform

  ! deferred public methods
  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy Jiang-Shu polynomial coefficients for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_polynomials_js_nonuniform), intent(inout) :: self   !< WENO polynomials.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%coef)) deallocate(self%coef)
  if (allocated(self%poly)) deallocate(self%poly)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure subroutine create(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO polynomials coefficients for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_polynomials_js_nonuniform), intent(inout) :: self   !< WENO polynomials.
  integer(I_P), intent(in) :: S                               !< Number of stencils used.
  real(R_P)                :: coord_l, coord_r, coord_tar     !< Abscissas of the reconstruction points, left and right interfaces.
  real(R_P)                :: stencil_coord(1:, 1 - S:)       !< Abscissas of the interpolation stencil, [1:2, 1-S:-1+S].
  real(R_P)                :: den, num_prod, num, frac, coeff !< Intermediate values for coefficients evaluation.
  integer(I_P)             :: s1, s2, m, l, q                 !< Counters.
  integer(I_P)             :: f, f1, f2                       !< Faces to be computed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  allocate(self%poly(1:2, 0:S - 1))
  allocate(self%coef(1:2, 0:S - 1, 0:S - 1))
  self%poly = 0._R_P
  associate(c => self%coef)
    do f = 1, 2  ! interfaces loop
      if (f==1) then
        coord_tar = coord_l
      else
        coord_tar = coord_r
      endif
      do s1 = 0, S - 1  ! stencils loop
        do s2 = 0, S - 1  ! values loop
          do m = s2 + 1, S
            den = 1._R_P
            num = 0._R_P
            do l = 0, S
              ! denominator
              if (l==m) cycle
              den = den * (stencil_coord(f,1 - S + s1 + m) - stencil_coord(f,1 - S + s1 + l))
              ! numerator
              ! numerator product
              num_prod = 1._R_P
              do q = 0, S
                if ((q==l).or.(q==m)) cycle
                num_prod = num_prod * (coord_tar - stencil_coord(f,1 - S + s1 + q))
              enddo
              ! numerator sum
              num = num + num_prod
            enddo
            frac = (num / den) * (stencil_coord(f,1 - S + s1 + s2) - stencil_coord(f,1 - S + s1 + s2 - 1))
            coeff = coeff + frac
          enddo
          c(f,s2,s1) = coeff
        enddo
      enddo
    enddo
  endassociate
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create

  pure subroutine description(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable, intent(out) :: string            !< String returned.
  character(len=1), parameter                :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO polynomial'//nl
  string = string//'  Based on the work by Črnjarić-Žic, Maćešić and Crnković "Efficient Implementation of WENO Schemes to non '// &
           'uniform meshes", see ANNALI DELL''UNIVERSITA'' DI FERRARA, 2007, vol. 57, issue 2, pp. 199--215, '// &
           'doi:10.1007/s11565-007-0013-1'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    poly(poly_coef,v)'//nl
  string = string//'  where:'//nl
  string = string//'    poly_coef: real(R_P), intent(IN), the polynomial coefficient of the value'//nl
  string = string//'    v: real(R_P), intent(IN), the selected value from the stencil'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure subroutine compute(self, S, stencil, f1, f2, ff)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the partial value of the interpolating polynomial.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(weno_polynomials_js_nonuniform), intent(inout) :: self                    !< WENO polynomial.
  integer(I_P),               intent(in)    :: S                       !< Number of stencils actually used.
  real(R_P),                  intent(in)    :: stencil(1:, 1 - S:)     !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P),               intent(in)    :: f1, f2, ff              !< Faces to be computed.
  integer(I_P)                              :: s1, s2, f               !< Counters
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  self%poly = 0._R_P
  do s1 = 0, S - 1 ! stencils loop
    do s2 = 0, S - 1 ! values loop
      do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
        self%poly(f, s1) = self%poly(f, s1) + self%coef(f, s2, s1) * stencil(f + ff, -s2 + s1)
      enddo
    enddo
  enddo
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute

!-----------------------------------------------------------------------------------------------------------------------------------
endmodule type_weno_polynomials_js_nonuniform
