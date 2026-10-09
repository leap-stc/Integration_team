module SAM_consts_mod
  !! Physical constants from the SAM Model required for variable conversions
  
  use netcdf

  implicit none
  public

  ! From params.f90 in SAM
  ! These are constants used in equations as part of the parameterisation.
  ! Since the net is trained with them they should be left as defined here, rather than
  ! forced to match the external model

  ! Physical Constants:

  != unit (m / s^2)
  real(8), parameter :: ggr = 9.81
      !! Gravity acceleration, m/s2

  != unit J / kg :: lcond, lsub, lfus
  real(8), parameter :: lfus = 0.3336e+06
      !! Latent heat of fusion
  real(8), parameter :: lcond = 2.5104e+06
      !! Latent heat of condensation
  real(8), parameter :: lsub = 2.8440e+06
      !! Latent heat of sublimation, J/kg

  != unit (J / kg) / K :: cp
  real(8), parameter :: cp = 1004.
      !! Specific heat of air
  
  ! unit K :: fac_cond, fac_fus, fac_sub
  real(8), parameter :: fac_cond = lcond/cp
      !!
  real(8), parameter :: fac_fus = lfus/cp
      !!
  real(8), parameter :: fac_sub = lsub/cp
      !!

  ! Temperatures limits for various hydrometeors
  != unit K :: tprmin
  real(8), parameter :: tprmin = 268.16
      !! Minimum temperature for rain
  
  != unit K :: tprmax
  real(8), parameter :: tprmax = 283.16
      !! Maximum temperature for snow+graupel, K
  
  != unit K :: tbgmin
  real(8), parameter :: tbgmin = 253.16
      !! Minimum temperature for cloud water.
  
  != unit K :: tbgmin
  real(8), parameter :: tbgmax = 273.16
      !! Maximum temperature for cloud ice, K

  != unit K :: tbgmin
  real(8), parameter :: tgrmin = 223.16
      !! Maximum temperature for snow, K

  != unit K :: tbgmin
  real(8), parameter :: tgrmax = 283.16
      !! Maximum temperature for graupel, K

  ! Misc. microphysics variables
  ! != unit 1 / K :: a_pr
  real(8), parameter :: a_pr = 1./(tprmax-tprmin)
      !! Misc. microphysics variables
  
  != unit 1 / K :: a_bg
  real(8), parameter :: a_bg = 1./(tbgmax-tbgmin)
      !! Misc. microphysics variables

  real(8), parameter :: an = 1./(tbgmax-tbgmin)
  real(8), parameter :: bn = tbgmin * an
  real(8), parameter :: ap = 1./(tprmax-tprmin)
  real(8), parameter :: bp = tprmin * ap
  real(8), parameter :: fac1 = fac_cond+(1+bp)*fac_fus
  real(8), parameter :: fac2 = fac_fus*ap
  real(8), parameter :: ag = 1./(tgrmax-tgrmin)

  ! ---------------------------
  ! SAM Grid Variables
  ! ---------------------------
  integer, parameter :: input_ver_dim = 49 ! XL XL XL
      !! The number of cells in a SAM atmospheric column on which the neural net was trained
  
  ! Outputs from NN are supplied at lowest 30 half-model levels for sedimentation fluxes,
  ! and at 29 levels for fluxes (as flux at bottom boundary is zero).
  integer, parameter :: nrf = 49 ! XL XL XL 
      !! number of vertical levels the NN uses
  integer, parameter :: nrfq = nrf - 1
      !! number of vertical levels the NN uses when boundary condition is set to 0

  real(8), parameter :: dt_sam = 24.0
      !! SAM timestep in seconds

!---------------------------------------------------------------------
! Functions and Subroutines

contains

  != unit 1 :: omegan
  real(8) function omegan(tabs)
      != unit K :: tabs
      real(8), intent(in) :: tabs
          !! Absolute temperature

      omegan = max(0., min(1., (tabs-tbgmin)*a_bg))

      return

  end function omegan


  subroutine check(err_status)
  !! Check error status after netcdf call and print message for
  !! error codes.

  integer, intent(in) :: err_status
    !! error status from nf90 function

  if(err_status /= nf90_noerr) then
     write(*, *) trim(nf90_strerror(err_status))
  end if

  end subroutine check

end module SAM_consts_mod
