module wenoof_smoothness_indicators_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing smoothness indicators for WENO schemes on non-uniform meshes.
!<
!< @note The provided polynomials implement the smoothness indicators defined in *Grid adaptation with WENO schemes for non-uniform
!< grids to solve convection-dominated partial differential equations*, J. Smit, M. van Sint Annaland, J. A. M. Kuipers,
!< Chemical Engineering Science, vol. 60, issue 10, pp. 2609--2619, doi:10.1016/j.ces.2004.12.017
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P
use wenoof_smoothness_indicators_abstract
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: IS_js_nonuniform, associate_IS_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(IS) :: IS_js_nonuniform
  !< Smoothness indicators for WENO schemes on non-uniform meshes object.
  !<
  !< @note The provided polynomials implement the smoothness indicators defined in *Grid adaptation with WENO schemes for non-uniform
  !< grids to solve convection-dominated partial differential equations*, J. Smit, M. van Sint Annaland, J. A. M. Kuipers,
  !< Chemical Engineering Science, vol. 60, issue 10, pp. 2609--2619, doi:10.1016/j.ces.2004.12.017
  private
  contains
    procedure, pass(self), public :: destroy
    procedure, pass(self), public :: create
    procedure, nopass,     public :: description
    procedure, pass(self), public :: compute
endtype IS_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public, non TBP
  function associate_IS_js_nonuniform(IS_input) result(IS_pointer)
    !< Check the type of smoothness indicator passed as input and return a Jiang-Shu IS associated to smoothness indicator.
    class(IS),            intent(in), target  :: IS_input   !< Input smoothness indicator.
    class(IS_js_nonuniform),          pointer :: IS_pointer !< Non uniform Jiang Shu smoothness indicators.

    select type(IS_input)
      type is(IS_js_nonuniform)
        IS_pointer => IS_input
      class default
        write(stderr, '(A)')'error: wrong smoothness indicator type chosen'
        stop
    end select
  end function associate_IS_js_nonuniform

  ! deferred public methods
  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy Jiang-Shu and Gerolymos-Sénéchal-Vallet WENO smoothness indicators coefficients for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(IS_js_nonuniform), intent(inout) :: self   !< WENO smoothenss indicators.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%coef)) deallocate(self%coef)
  if (allocated(self%si)) deallocate(self%si)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure subroutine create(self, S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO smoothness indicators coefficients for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(IS_js_nonuniform), intent(inout) :: self   !< WENO smoothness indicators.
  integer(I_P),            intent(in)    :: S      !< Number of stencils used.
  real(R_P)                   :: coord(- S:)       !< Abscissas of the whole interpolation stencil, [-S:-1+S].
  real(R_P)                   :: a, b, d, e        !< Intermediate values.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  allocate(self%si(1:2, 0:S - 1))
  self%si = 0._R_P
  allocate(self%coef(1:2, 0:S - 1, 0:S - 1))
  associate(c => self%coef)
    select case(S)
    case(2) ! 3rd order
      a = coord(2-S) - coord(1-S)

      ! stencil 0
      b = coord(S-2) - coord( -S)
      !              i*i                 ;             (i-1)*i
      c(0,0,0) = 4._R_P * (a / b)**2._R_P; c(1,0,0) = -8._R_P * (a / b)**2._R_P
      !               /                  ;             (i-1)*(i-1)
      c(0,1,0) = 0._R_P                  ; c(1,1,0) =  4._R_P * (a / b)**2._R_P

      ! stencil 1
      b = coord(S-1) - coord(1-S)
      !             (i+1)*(i+1)          ;                (i+1)*i
      c(0,0,1) = 4._R_P * (a / b)**2._R_P; c(1,0,1) = -8._R_P * (a / b)**2._R_P
      !               /                  ;                 i*i
      c(0,1,1) = 0._R_P                  ; c(1,1,1) =  4._R_P * (a / b)**2._R_P

    case(3) ! 5rd order
      ! stencil 0
      a =  4._R_P  * ((coord(3-S) - coord(2-S))/(coord(3-S) - coord( -S)))**2._R_P
      b = ((coord(3-S) - coord(1-S)) * (coord(2-S) - coord(1-S)) + &
           10._R_P * (coord(3-S) - coord(2-S))**2._R_P) / (coord(2-S) - coord( -S))**2._R_P
      d = (20._R_P * (coord(3-S) - coord(2-S))**2._R_P + &
            2._R_P * (coord(3-S) - coord(1-S)) * (coord(2-S) - coord(1-S)) + &
                     (coord(3-S) - coord( -S)) * (coord(2-S) + coord(2-S) - 2._R_P * (coord(1-S)))) / &
          ((coord(3-S) - coord(1-S)) * (coord(2-S) - coord( -S)))
      e = (10._R_P * (coord(3-S) - coord(2-S))**2._R_P + &
          ( 2._R_P *  coord(3-S) - coord( -S) - coord(1-S)) * &
          (coord(3-S) - coord(2-S) - coord(3-S) - coord( -S))) / (coord(3-S) - coord(1-S))**2._R_P
      !      i*i       ;             (i-1)*i             ;       (i-2)*i
      c(0,0,0) = a * e ; c(1,0,0) = -a * (d + 2._R_P * e); c(2,0,0) =  a * d
      !      /         ;             (i-1)*(i-1)         ;       (i-2)*(i-1)
      c(0,1,0) = 0._R_P; c(1,1,0) =  a * (b + d + e)     ; c(2,1,0) = -2._R_P * a * (b + d)
      !      /         ;              /                  ;       (i-2)*(i-2)
      c(0,2,0) = 0._R_P; c(1,2,0) =  0._R_P              ; c(2,2,0) =  a * b

      ! stencil 1
      a =  4._R_P  * ((coord(3-S) - coord(2-S))/(coord(4-S) - coord(1-S)))**2._R_P
      b = ((coord(4-S) - coord(2-S)) * (coord(4-S) - coord(3-S)) + &
           10._R_P * (coord(3-S) - coord(2-S))**2._R_P) / (coord(3-S) - coord(1-S))**2._R_P
      d = (20._R_P * (coord(3-S) - coord(2-S))**2._R_P - &
                     (coord(4-S) - coord(3-S)) * (coord(2-S) - coord(1-S)) - &
                     (coord(4-S) - coord(2-S)) * (coord(3-S) - coord(1-S))) / &
          ((coord(4-S) - coord(2-S)) * (coord(3-S) - coord(1-S)))
      e = (10._R_P * (coord(3-S) - coord(2-S))**2._R_P + &
          (coord(2-S) - coord(1-S)) * (coord(3-S) - coord(1-S))) / (coord(4-S) - coord(2-S))**2._R_P
      !     (i+1)*(i+1)  ;              i*(i+1)            ;       (i-1)*(i+1)
      c(0,0,1) =   a * e ; c(1,0,1) = -a * (d + 2._R_P * e); c(2,0,1) =  a * d
      !      /           ;              i*i                ;       (i-1)*i
      c(0,1,1) =   0._R_P; c(1,1,1) =  a * (b + d + e)     ; c(2,1,1) = -a * (2._R_P * b + d)
      !      /           ;              /                  ;       (i-1)*(i-1)
      c(0,2,1) =   0._R_P; c(1,2,1) =  0._R_P              ; c(2,2,1) =  a * b

      ! stencil 2
      a =  4._R_P  * ((coord(3-S) - coord(2-S))/(coord(5-S) - coord(2-S)))**2._R_P
      b = ((coord(4-S) - coord(2-S)) * (coord(4-S) - coord(3-S)) + &
           10._R_P * (coord(3-S) - coord(2-S))**2._R_P) / (coord(5-S) - coord(3-S))**2._R_P
      d = (20._R_P * (coord(3-S) - coord(2-S))**2._R_P + &
            2._R_P * (coord(4-S) - coord(2-S)) * (coord(4-S) - coord(3-S)) - &
                     (coord(5-S) - coord(2-S)) * (2._R_P * coord(4-S) - coord(3-S) - coord(2-S))) / &
          ((coord(5-S) - coord(3-S)) * (coord(4-S) - coord(2-S)))
      e = (10._R_P * (coord(3-S) - coord(2-S))**2._R_P + &
          (coord(5-S) - coord(4-S) - 2._R_P * coord(2-S)) * (coord(5-S) + coord(4-S) - coord(3-S) - coord(2-S))) / &
          (coord(4-S) - coord(2-S))**2._R_P
      !     (i+2)*(i+2)  ;             (i+1)*(i+2)         ;        i*(i+2)
      c(0,0,2) =   a * b ; c(1,0,2) = -a * (2._R_P * b + d); c(2,0,2) =  a * d
      !      /           ;             (i+1)*(i+1)         ;        i*(i+1)
      c(0,1,2) =   0._R_P; c(1,1,2) =  a * (b + d + e)     ; c(2,1,2) = -a * (d + 2._R_P * e)
      !      /           ;              /                  ;        i*i
      c(0,2,2) =   0._R_P; c(1,2,2) =  0._R_P              ; c(2,2,2) =  a * e
    endselect
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create

  pure subroutine description(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing WENO smoothness indicator for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable, intent(out) :: string !< String returned.
  character(len=1), parameter                :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO smoothness indicators'//nl
  string = string//'  Based on the work by J. Smit, M. van Sint Annaland and J. A. M. Kuipers "Grid adaptation with WENO '// &
           'schemes for non-uniform grids to solve convection-dominated partial differential equations", see Chemical '// &
           'Engineering Science, vol. 60, issue 10, pp. 2609--2619, doi:10.1016/j.ces.2004.12.017'//nl
  string = string//'  The "compute" method has the following public API'//nl
  string = string//'    IS(smooth_coef,v1,v2)'//nl
  string = string//'  where:'//nl
  string = string//'    smooth_coef: real(R_P), intent(IN), the smoothness indicator coefficient of the value'//nl
  string = string//'    v1: real(R_P), intent(IN), the pivotal value from the stencil'//nl
  string = string//'    v2: real(R_P), intent(IN), the second value from the stencil'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description

  pure subroutine compute(self, S, stencil, f1, f2, ff)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Compute the partial value of the smoothness indicator of a single WENO interpolating polynomial for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(IS_js_nonuniform), intent(inout) :: self                    !< WENO smoothness indicator.
  integer(I_P),            intent(in)    :: S                       !< Number of stencils actually used.
  real(R_P),               intent(in)    :: stencil(1:, 1 - S:)     !< Stencil used for the interpolation, [1:2, 1-S:-1+S].
  integer(I_P),            intent(in)    :: f1, f2, ff              !< Faces to be computed.
  real(R_P)                              :: coef(1:2, 0:S-1, 0:S-1) !< Coefficients of the smoothness indicators.
  integer(I_P)                           :: s1, s2, s3, f           !< Counters
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  do s1 = 0, S - 1 ! stencils loop
    do f = f1, f2 ! 1 => left interface (i-1/2), 2 => right interface (i+1/2)
      self%si(f, s1) = 0._R_P
      do s2 = 0, S - 1
        do s3 = 0, S - 1
          self%si(f, s1) = self%si(f, s1) + self%coef(s3, s2, s1) * stencil(f + ff, s1 - s3) * stencil(f + ff, s1 - s2)
        enddo
      enddo
    enddo
  enddo
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine compute
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule wenoof_smoothness_indicators_js_nonuniform
