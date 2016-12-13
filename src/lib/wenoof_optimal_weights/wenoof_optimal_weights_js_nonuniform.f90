module wenoof_optimal_weights_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------
!< Module providing Jiang-Shu and Gerolymos-Sénéchal-Vallet optimal weights for WENO schemes on non uniform grids.
!<
!< @note The provided polynomials implement the optimal weights defined in *Efficient implementation of WENO schemes to
!< nonuniform meshes*, Nelida Črnjarić-Žic, Senka Maćešić, Bojan Crnković, ANNALI DELL'UNIVERSITA' DI FERRARA, 2007, vol. 57,
!< issue 2, pp. 199--215, doi:10.1007/s11565-007-0013-1
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use penf, only : I_P, R_P
use wenoof_optimal_weights_abstract
use wenoof_polynomials_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
save
public :: optimal_weights_js_nonuniform, associate_weights_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------------------------------------------------
type, extends(optimal_weights):: optimal_weights_js_nonuniform
  !< Jiang-Shu and Gerolymos-Sénéchal-Vallet WENO optimal weights object.
  !<
  !< @note The provided polynomials implement the optimal weights defined in *Efficient implementation of WENO schemes to
  !< nonuniform meshes*, Nelida Črnjarić-Žic, Senka Maćešić, Bojan Crnković, ANNALI DELL'UNIVERSITA' DI FERRARA, 2007, vol. 57,
  !< issue 2, pp. 199--215, doi:10.1007/s11565-007-0013-1
  private
  contains
    ! deferred public methods
    procedure, pass(self), public :: destroy
    procedure, pass(self), public :: create
    procedure, nopass,     public :: description
endtype optimal_weights_js_nonuniform
!-----------------------------------------------------------------------------------------------------------------------------------
contains
  ! public, non TBP
  function associate_WENO_weights_js_nonuniform(weights_input) result(weights_pointer)
    !< Check the type of the optimal weights passed as input and return Jiang-Shu optimal weights associated to the optimal weights.
    class(optimal_weights),   intent(in), target  :: weights_input   !< Input optimal weights.
    class(optimal_weights_js_nonuniform), pointer :: weights_pointer !< Non uniform Jiang Shu optimal weights.

    select type(weights_input)
      type is(optimal_weights_js_nonuniform)
        weights_pointer => weights_input
      class default
        write(stderr, '(A)')'error: wrong optimal weights type chosen'
        stop
    end select
  end function associate_weights_js_nonuniform

  ! deferred public methods
  pure subroutine destroy(self)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Destroy Jiang-Shu and Gerolymos-Sénéchal-Vallet WENO optimal weights for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(optimal_weights_js_nonuniform), intent(inout)  :: self   !< WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  if (allocated(self%opt)) deallocate(self%opt)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine destroy

  pure subroutine create(self,S)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Create WENO optimal weights for non uniform grids.
  !---------------------------------------------------------------------------------------------------------------------------------
  class(optimal_weights_js_nonuniform), intent(inout) :: self !< WENO Jiang-Shu optimal weights.
  integer(I_P), intent(in) :: S                               !< Number of stencils used.
  real(R_P), allocatable   :: coef(:,:)                       !< Polynomial coefficients on the whole recontruction stencil.
  real(R_P)                :: coord_l, coord_r, coord_tar     !< Abscissas of the reconstruction points, left and right interfaces.
  real(R_P)                :: stencil_coord(1:, 1 - S:)       !< Abscissas of the whole interpolation stencil, [1:2, 1-S:-1+S].
  real(R_P)                :: den, num_prod, num, frac, coeff !< Intermediate values for coefficients evaluation.
  integer(I_P)             :: s1, s2, m, l, q                 !< Counters.
  integer(I_P)             :: f, f1, f2                       !< Faces to be computed.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  call self%destroy
  allocate(self%opt(1:2, 0:S - 1))
  allocate(coef(1:2, 0:2*S - 2))
  do f = 1, 2  ! interfaces loop
    if (f==1) then
      coord_tar = coord_l
    else
      coord_tar = coord_r
    endif
    do s1 = 0, 2*S - 2  ! values loop
      do m = s2 + 1, 2*S - 2
        den = 1._R_P
        num = 0._R_P
        do l = 0, 2*S - 2
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
      coef(f,s2) = coeff
    enddo
  enddo
  associate(opt => self%opt, poly_coef => weno_polynomials_js_nonuniform%coef)
    do f = 1, 2  ! interfaces loop
      do s1 = 0, S - 1
        coeff = 0._R_P
        do s2 = 0, s1 - 1
          coeff = coeff + opt(f,s2) * poly_coef(f,s2,S - s2)
        enddo
        opt(f,s1) = (coef(f,s1) - coeff) / poly_coef(f,s1,0)
      enddo
    enddo
  endassociate
  deallocate(coef)
  return
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine create

  pure subroutine description(string)
  !---------------------------------------------------------------------------------------------------------------------------------
  !< Return a string describing Jiang-Shu and Gerolymos-Sénéchal-Vallet WENO optimal weights.
  !---------------------------------------------------------------------------------------------------------------------------------
  character(len=:), allocatable,  intent(out) :: string !< String returned.
  character(len=1), parameter                 :: nl=new_line('a')  !< New line character.
  !---------------------------------------------------------------------------------------------------------------------------------

  !---------------------------------------------------------------------------------------------------------------------------------
  string = 'WENO optimal weights,'//nl
  string = string//'  Based on the work by Črnjarić-Žic, Maćešić and Crnković "Efficient Implementation of WENO Schemes to non '// &
           'uniform meshes", see ANNALI DELL''UNIVERSITA'' DI FERRARA, 2007, vol. 57, issue 2, pp. 199--215, '// &
           'doi:10.1007/s11565-007-0013-1'//nl
  string = string//'    The optimal weights are allocated in a two-dimensional array, in which the first index'//nl
  string = string//'    is the face selected (1 => i-1/2, 2 => i+1/2) and the second index is the number of the stencil '//nl
  string = string//'    (from 0 to S-1)'
  !---------------------------------------------------------------------------------------------------------------------------------
  endsubroutine description
!-----------------------------------------------------------------------------------------------------------------------------------
endmodule wenoof_optimal_weights_js_nonuniform
