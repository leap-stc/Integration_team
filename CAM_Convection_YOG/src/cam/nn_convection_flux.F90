module nn_convection_flux_mod
    !! Code to perform the Convection Flux parameterisation
    !! and interface to the nn_cf_net Neural Net
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
! Libraries to use
use nn_cf_net_mod, only: nn_cf_net_init, net_forward !nn_cf_net_finalize
use SAM_consts_mod, only: fac_cond, fac_fus, tprmin, a_pr, input_ver_dim, &
                          nrf, nrfq
use FTorch_cesm_interface, only: torch_model
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux, nn_convection_flux_init, nn_convection_flux_finalize, &
        esati, rsati, rsatw, dtrsatw, dtrsati


!---------------------------------------------------------------------
! local/private data

! Neural Net parameters used throughout module

integer :: n_inputs, n_outputs, n_h1, n_h2, n_h3, n_h4, n_features_out, n_lev
    !! Length of input/output vector to the NN
type(torch_model) :: model_ftorch
logical :: do_init=.true.
    !! model initialisation is yet to be performed
integer, allocatable :: feature_out_sizes(:)

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_init(nn_filename, metadata_filename)
        !! Initialise the NN module using FTorch

        character(len=136), intent(in) :: nn_filename
            !! PyTorch model filename (.pt file) from which to load model 
        character(len=136), intent(in) :: metadata_filename
        !! NetCDF file with dimensions and scaling parameters
        
        ! Local variables for error handling
        character(128) :: errstring
        integer :: iulog
    
        ! Set logging unit (use standard output or your logging system)
        iulog = 6  ! Standard output, or use your model's logging unit

        ! Initialise the Neural Net from PyTorch file
        call nn_cf_net_init(nn_filename, metadata_filename, n_inputs, n_outputs, nrf, &
                        model_ftorch, iulog, errstring)
        
        ! Check for errors
        if (len_trim(errstring) > 0) then
            write(iulog, *) 'ERROR in nn_convection_flux_init: ', trim(errstring)
            stop 'Failed to initialize neural network'
        end if

        write(iulog, *) 'Neural network initialized successfully'

        ! Set init flag as complete
        do_init = .false.

    end subroutine nn_convection_flux_init


    subroutine nn_convection_flux(tabs_i, q_i, &
                                  tabs, &
                                  t, q, &
                                  rho, adz, dz, dtn, &
                                  precsfc)
        !! Interface to the neural net that applies physical constraints and reshaping
        !! of variables.
        !! Operates on subcycle of timestep dtn to update t, q, precsfc, and prec_xy

        ! -----------------------------------
        ! Input Variables
        ! -----------------------------------

        ! ---------------------
        ! Fields from beginning of time step used as NN inputs
        ! ---------------------
        != unit s :: tabs_i
        real(8), intent(in) :: tabs_i(:, :)
            !! Temperature

        != unit 1 :: q_i
        real(8), intent(in) :: q_i(:, :)
            !! Non-precipitating water mixing ratio

        ! ---------------------
        ! Other fields from SAM
        ! ---------------------
        != unit K :: tabs
        real(8), intent(in) :: tabs(:, :)
            !! absolute temperature
        
        != unit (J / kg) :: t
        real(8), intent(inout) :: t(:, :)
            !! Liquid Ice static energy (cp*T + g*z − L(qliq + qice) − Lf*qice)

        != unit 1 :: q
        real(8), intent(inout) :: q(:, :)
            !! total water

        ! ---------------------
        ! reference vertical profiles:
        ! ---------------------
        != unit (kg / m**3) :: rho
        real(8), intent(in) :: rho(:)
            !! air density at pressure levels

        ! != unit mb :: pres
        ! real(8), intent(in) pres(nzm)
        !     !! pressure,mb at scalar levels

        != unit 1 :: adz
        real(8), intent(in) :: adz(:)
            !! ratio of the pressure level grid height spacing [m] to dz (lowest dz spacing)

        ! ---------------------
        ! Single value parameters from model/grid
        ! ---------------------
        != unit m :: dz
        real(8), intent(in) :: dz
            !! grid spacing in z direction for the lowest grid layer

        != unit s :: dtn
        real(8), intent(in) :: dtn
            !! current dynamical timestep (can be smaller than dt due to subcycling)

        ! -----------------------------------
        ! Output Variables
        ! -----------------------------------
        != unit (J / kg) :: t_delta_adv, t_delta_auto, t_delta_sed
        != unit 1 :: q_delta_adv, q_delta_auto, q_delta_sed
        real(8), dimension(size(tabs_i, 1), size(tabs_i, 2)) :: t_delta_adv, q_delta_adv, &
                                                             t_delta_auto, q_delta_auto, &
                                                             t_delta_sed, q_delta_sed
            !! delta values of t and q generated by the NN

        != unit kg / m**2 :: precsfc
        real(8), intent(out), dimension(:)   :: precsfc
            !! Surface precipitation due to autoconversion and sedimentation

        ! -----------------------------------
        ! Local Variables
        ! -----------------------------------
        integer  i, k, dim_counter, out_dim_counter
        integer ncol
            !! Number of columns in a subdomain
        integer nzm
            !! Number of z points in a subdomain - 1
        ! real(8) :: omn
        !     !! variable to store result of omegan function
        !     !! Note that if using you will need to import omegan() from SAM_consts_mod

        ! Other variables
        real(8),   dimension(nrf) :: omp, fac
        real(8),   dimension(size(tabs_i, 2)) :: rsat, irhoadz, irhoadzdz

        ! -----------------------------------
        ! variables for NN
        ! -----------------------------------
        real(4), dimension(n_inputs) :: features
            !! Vector of input features for the NN
        real(4), dimension(n_outputs) :: outputs
            !! vector of output features from the NN
        ! NN outputs
        real(8),   dimension(nrf) :: t_flux_adv, q_flux_adv, q_tend_auto, &
                                  q_sed_flux, t_rad_rest_tend

        ncol = size(tabs_i, 1)
        nzm = size(tabs_i, 2)

        ! Check that we have initialised all of the variables.
        if (do_init) call error_mesg('NN has not yet been initialised using nn_convection_flux_init.')

        ! Define useful variables relating to grid spacing to convert fluxes to tendencies
        do k=1,nzm
            irhoadz(k) = dtn/(rho(k)*adz(k)) !  Temporary factor for below
            irhoadzdz(k) = irhoadz(k)/dz ! 2.0 * dtn / (rho(k)*(z(k+1) - z(k-1))) [(kg.m/s)^-1]
        end do

        ! The NN operates on atmospheric columns which have been flattened into 2D
        do i=1,ncol
            ! Initialize variables
            features = 0.
            dim_counter = 0
            outputs = 0.
            t_rad_rest_tend = 0.
            t_flux_adv = 0.
            q_flux_adv = 0.
            t_delta_adv(i,:) = 0.
            q_delta_adv(i,:) = 0.
            q_tend_auto = 0.
            t_delta_auto(i,:) = 0.
            q_delta_auto(i,:) = 0.
            q_sed_flux = 0.
            t_delta_sed(i,:) = 0.
            q_delta_sed(i,:) = 0.
            precsfc(i) = 0.
            omp = 0.
            fac = 0.

            !-----------------------------------------------------
            ! Combine all features into one vector

            ! Add temperature as input feature
            features(dim_counter+1:dim_counter + input_ver_dim) = real(tabs_i(i,1:input_ver_dim),4)
            dim_counter = dim_counter + input_ver_dim

            ! Add non-precipitating water mixing ratio as input feature using random forest (rf) approach from earlier paper
            ! Currently we do not use rf_uses_rh option, but may add it back in later
            ! if (rf_uses_rh) then
            ! ! If using generalised relative humidity convert non-precip. water content to rel. hum
            !     do k=1,nzm
            !         omn = omegan(tabs(i,j,k))
            !         rsat(k) = omn * rsatw(tabs(i,j,k),pres(k)) + (1.-omn) * rsati(tabs(i,j,k),pres(k))
            !     end do
            !     features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,j,1:input_ver_dim)/rsat(1:input_ver_dim),4)
            !     dim_counter = dim_counter + input_ver_dim
            ! else
            ! ! if using non-precipitating water as water content
                features(dim_counter+1:dim_counter+input_ver_dim) = real(q_i(i,1:input_ver_dim),4)
                dim_counter =  dim_counter + input_ver_dim
            ! endif

            !-----------------------------------------------------
            ! Call the forward method of the NN on the input features
            ! Scaling and normalisation done as layers in NN

            call net_forward(features, model_ftorch,  outputs)

            !-----------------------------------------------------
            ! Separate physical outputs from NN output vector

            out_dim_counter = 0

            ! Moist Static Energy radiative + rest(microphysical changes) tendency
            t_rad_rest_tend(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
            out_dim_counter = out_dim_counter + nrf

            ! Moist Static Energy advective subgrid flux
            ! BC: advection surface flux is zero
            t_flux_adv(1) = 0.0
            t_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
            out_dim_counter = out_dim_counter + nrfq

            ! Total non-precip. water mix. ratio advective flux
            ! BC: advection surface flux is zero
            q_flux_adv(1) = 0.0
            q_flux_adv(2:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrfq)
            out_dim_counter = out_dim_counter + nrfq

            ! Total non-precip. water autoconversion tendency
            q_tend_auto(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)
            out_dim_counter = out_dim_counter + nrf

            ! total non-precip. water mix. ratio ice-sedimenting flux
            q_sed_flux(1:nrf) = outputs(out_dim_counter+1:out_dim_counter+nrf)

            !-----------------------------------------------------
            ! Apply physical constraints and update q and t

            ! Non-precip. water content must be >= 0, so ensure advective fluxes
            ! will not reduce it below 0 anywhere
            do k=2,nrf
                if (q_flux_adv(k) < 0) then
                    ! If flux is negative ensure we don't lose more than is already present
                    if ( q(i,k) < -q_flux_adv(k)* irhoadzdz(k)) then
                        q_flux_adv(k) = -q(i,k)/irhoadzdz(k)
                    end if
                else
                    ! If flux is positive ensure we don't gain more than is in the box below
                    if (q(i,k-1) < q_flux_adv(k)* irhoadzdz(k)) then
                        q_flux_adv(k) = q(i,k-1)/irhoadzdz(k)
                    end if
                end if
            end do

            ! Convert advective fluxes to deltas
            do k=1,nrf-1
                t_delta_adv(i,k) = - (t_flux_adv(k+1) - t_flux_adv(k)) * irhoadzdz(k)
                q_delta_adv(i,k) = - (q_flux_adv(k+1) - q_flux_adv(k)) * irhoadzdz(k)
            end do
            ! Enforce boundary condition at top of column
            t_delta_adv(i,nrf) = - (0.0 - t_flux_adv(nrf)) * irhoadzdz(nrf)
            q_delta_adv(i,nrf) = - (0.0 - q_flux_adv(nrf)) * irhoadzdz(nrf)
            ! q must be >= 0 so ensure delta won't reduce it below zero
            do k=1,nrf
                if (q(i,k) < -q_delta_adv(i,k)) then
                    q_delta_adv(i,k) = -q(i,k)
                end if
            end do

            ! Update q and t with delta values
            q(i,1:nrf) = q(i,1:nrf) + q_delta_adv(i,1:nrf)
            t(i,1:nrf) = t(i,1:nrf) + t_delta_adv(i,1:nrf)

            ! ensure autoconversion tendency won't reduce q below 0
            do k=1,nrf
                omp(k) = max(0.,min(1.,(tabs(i,k)-tprmin)*a_pr))
                fac(k) = (fac_cond + fac_fus * (1.0 - omp(k)))
                if (q_tend_auto(k) < 0) then
                    q_delta_auto(i,k) = - min(-q_tend_auto(k) * dtn, q(i,k))
                else
                    q_delta_auto(i,k) = q_tend_auto(k) * dtn
                endif
            end do

            ! Update with autoconversion q and t deltas (dt = -dq*(latent_heat/cp))
            q(i,1:nrf) = q(i,1:nrf) + q_delta_auto(i,1:nrf)
            t_delta_auto(i,1:nrf) = - q_delta_auto(i,1:nrf)*fac(1:nrf)
            t(i,1:nrf) = t(i,1:nrf) + t_delta_auto(i,1:nrf)

            ! Ensure sedimenting ice will not reduce q below zero anywhere
            do k=2,nrf
                if (q_sed_flux(k) < 0) then
                    ! If flux is negative ensure we don't lose more than is already present
                    if ( q(i,k) < -q_sed_flux(k)* irhoadzdz(k)) then
                        q_sed_flux(k) = -q(i,k)/irhoadzdz(k)
                    end if
                else
                    ! If flux is positive ensure we don't gain more than is in the box below
                    if (q(i,k-1) < q_sed_flux(k)* irhoadzdz(k)) then
                        q_sed_flux(k) = q(i,k-1)/irhoadzdz(k)
                    end if
                end if
            end do

            ! Convert sedimenting fluxes to deltas
            do k=1,nrf-1 ! One level less than I actually use
                q_delta_sed(i,k) = - (q_sed_flux(k+1) - q_sed_flux(k)) * irhoadzdz(k)
            end do
            ! Enforce boundary condition at top of column
            q_delta_sed(i,nrf) = - (0.0 - q_sed_flux(nrf)) * irhoadzdz(nrf)
            ! q must be >= 0 so ensure delta won't reduce it below zero
            do k=1,nrf
                if (q_delta_sed(i,k) < 0) then
                    q_delta_sed(i,k) = min(-q_delta_sed(i,k), q(i,k))
                    q_delta_sed(i,k) = -q_delta_sed(i,k)
                end if
            end do

            ! Update q and t with sed q and t deltas (dt = -dq*(latent_heat/cp))
            q(i,1:nrf) = q(i,1:nrf) + q_delta_sed(i,1:nrf)
            t(i,1:nrf) = t(i,1:nrf) - q_delta_sed(i,1:nrf)*(fac_fus+fac_cond)

            ! Apply radiation rest tendency (multiply by dtn to get delta value)
            t(i,1:nrf) = t(i,1:nrf) + t_rad_rest_tend(1:nrf)*dtn

            ! Calculate surface precipitation
            ! Combination of sedimentation at surface, and autoconversion in the column
            ! Apply sedimenting flux at surface to get rho*dq term
            precsfc(i) = precsfc(i) - q_sed_flux(1) * irhoadzdz(1) * rho(1) * adz(1) !!  *dtn/dz
            ! Loop up column for all autoconverted precipitation
            do k=1,nrf
                precsfc(i) = precsfc(i) - q_delta_auto(i,k) * rho(k) * adz(k)
            end do
            precsfc(i) = precsfc(i) * dz

            ! As a final check enforce q must be >= 0.0
            do k = 1,nrf
                q(i,k) = max(0.,q(i,k))
            end do
        end do
        ! End of loop over columns

    end subroutine nn_convection_flux


    subroutine nn_convection_flux_finalize()
        !! Finalize the NN module

        ! Finalize the Neural Net deallocating arrays
        ! call nn_cf_net_finalize()

    end subroutine nn_convection_flux_finalize


    !-----------------------------------------------------------------

    subroutine error_mesg (message)
      character(len=*), intent(in) :: message
          !! message to be written to output   (character string)

      ! Since masterproc is a SAM variable adjust to print from all proc for now
      ! if(masterproc) print*, 'Neural network  module: ', message
      print*, 'Neural network  module: ', message
      stop

    end subroutine error_mesg


    !-----------------------------------------------------------------
    ! Need rsatw functions to:
    !     - run with rf_uses_rh option (currently unused)
    !     - convert variables in CAM interface
    ! Ripped from SAM model:
    ! https://github.com/yaniyuval/Neural_nework_parameterization/blob/f81f5f695297888f0bd1e0e61524590b4566bf03/sam_code_NN/sat.f90
    ! Note r is used instead of q as q in CAM signifies a moist mixing ratio whilst SAM
    ! uses dry mixing ratios.
    !
    ! Saturation vapor pressure and mixing ratio.
    ! Based on Flatau et.al, (JAM, 1992:1507)

    != unit mb :: esatw
    real(8) function esatw(t)
      != unit K :: t
      real(8), intent(in) :: t  ! temperature (K)

      != unit :: a0
      != unit :: mb / k :: a1, a2, a3, a4, a5, a6, a7, a8
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
              6.105851, 0.4440316, 0.1430341e-1, &
              0.2641412e-3, 0.2995057e-5, 0.2031998e-7, &
              0.6936113e-10, 0.2564861e-13,-0.3704404e-15/
      !       6.11239921, 0.443987641, 0.142986287e-1, &
      !       0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
      !       0.640689451e-10, -0.952447341e-13,-0.976195544e-15/

      != unit K :: dt
      real(8) :: dt
      dt = max(-80.,t-273.16)
      esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function esatw


    != unit 1 :: rsatw
    real(8) function rsatw(t,p)
      != unit K :: t
      real(8), intent(in) :: t  ! temperature

      != unit mb :: p, esat
      real(8) :: p  ! pressure
      real(8) :: esat

      esat = esatw(t)
      rsatw = 0.622 * esat/max(esat, p-esat)
    end function rsatw


    real(8) function dtesatw(t)
      real(8), intent(in) :: t  ! temperature (K)
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
                0.443956472, 0.285976452e-1, 0.794747212e-3, &
                0.121167162e-4, 0.103167413e-6, 0.385208005e-9, &
               -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
      real(8) :: dt
      dt = max(-80.,t-273.16)
      dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function dtesatw


    real(8) function dtrsatw(t,p)
      real(8), intent(in) :: t  ! temperature (K)
      real(8), intent(in) :: p  ! pressure    (mb)
      dtrsatw=0.622*dtesatw(t)/p
    end function dtrsatw


    real(8) function esati(t)
      != unit K :: t
      real(8), intent(in) :: t  ! temperature
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
              6.11147274, 0.503160820, 0.188439774e-1, &
              0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
              0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
      real(8) :: dt
      dt = max(-80.0, t-273.16)
      esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    end function esati


    != unit 1 :: rsati
    real(8) function rsati(t,p)
      != unit t :: K
      real(8), intent(in) :: t  ! temperature

      != unit mb :: p
      real(8), intent(in) :: p  ! pressure

      != unit mb :: esat
      real(8) :: esat
      esat = esati(t)
      rsati = 0.622 * esat/max(esat,p-esat)
    end function rsati


    real(8) function dtesati(t)
      real(8), intent(in) :: t  ! temperature (K)
      real(8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 / &
              0.503223089, 0.377174432e-1,0.126710138e-2, &
              0.249065913e-4, 0.312668753e-6, 0.255653718e-8, &
              0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
      real(8) :: dt
      dt = max(-800. ,t-273.16)
      dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
      return
    end function dtesati


    real(8) function dtrsati(t,p)
      real(8), intent(in) :: t  ! temperature (K)
      real(8), intent(in) :: p  ! pressure    (mb)
      dtrsati = 0.622 * dtesati(t) / p
    end function dtrsati


end module nn_convection_flux_mod
