module nn_interface_CAM
    !! Interface to convection parameterisation for the CAM model
    !! Reference: https://doi.org/10.1029/2020GL091363
    !! Also see YOG20: https://doi.org/10.1038/s41467-020-17142-3

!---------------------------------------------------------------------
! Libraries to use
use netcdf
use nn_convection_flux_mod, only: nn_convection_flux, &
                                  nn_convection_flux_init, nn_convection_flux_finalize, &
                                  esati, rsati, rsatw, dtrsatw, dtrsati
use SAM_consts_mod, only: nrf, ggr, cp, tbgmax, tbgmin, tprmax, tprmin, &
                          fac_cond, fac_sub, fac_fus, &
                          a_bg, a_pr, an, bn, ap, bp, &
                          omegan, check
use cam_logfile,    only: iulog

implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  nn_convection_flux_CAM, &
        nn_convection_flux_CAM_init, nn_convection_flux_CAM_finalize

! Make these routines public for purposes of testing,
! otherwise not required outside this module.
public interp_to_sam, interp_to_cam
public SAM_var_conversion, CAM_var_conversion
public interp_flux_cam_midpoint

!---------------------------------------------------------------------
! local/private data

!= unit 1 :: nz_sam
integer :: nz_sam
    !! number of vertical values in the SAM sounding profiles

!= unit m :: z
!= unit hPa :: pres, presi
!= unit kg m-3 :: rho
!= unit 1 :: adz
!= unit K :: gamaz
real(8), allocatable, dimension(:) :: z, pres, presi, rho, adz, gamaz
    !! SAM sounding variables

!= unit m :: dz
real(8) :: dz
    !! grid spacing in z direction for the lowest grid layer

!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines

    subroutine nn_convection_flux_CAM_init(nn_filename, metadata_filename, sounding_filename)
        !! Initialise the NN module using FTorch

        character(len=136), intent(in) :: nn_filename
            !! PyTorch model filename (.pt file) from which to load model
        character(len=136), intent(in) :: metadata_filename
            !! NetCDF file with dimensions and scaling parameters
        character(len=136), intent(in) :: sounding_filename
            !! NetCDF filename from which to read SAM sounding and grid data

        ! Initialise the Neural Net from file and get info
        call nn_convection_flux_init(nn_filename, metadata_filename)

        ! Read in the SAM sounding and grid data required
        call sam_sounding_init(sounding_filename)

    end subroutine nn_convection_flux_CAM_init


    subroutine nn_convection_flux_CAM(pres_cam, pres_int_cam, pres_sfc_cam, &
                                      tabs_cam, qv_cam, qc_cam, qi_cam, phis_cam, zm_cam, &
                                      landfrac, ps_cam, &
                                      cp_cam, &
                                      dtn, &
                                      ncol, &
                                      precsfc, &
                                      dqi, dqv, dqc, ds, &
                                      q_delta_adv, q_delta_auto, q_delta_sed, &
                                      t_delta_adv, t_delta_auto, t_delta_sed, t_delta_rad, &
                                      t_flux_adv)

        !! Interface to the nn_convection parameterisation for the CAM model
      
        real(8), intent(in) :: dtn    ! Seconds
        integer, intent(in) :: ncol   ! Number of columns 'used' in each array
        
        real(8), intent(in) :: cp_cam
            !! specific heat capacity of dry air from CAM [J/kg/K]
        real(8), dimension(:,:), intent(in) :: pres_cam
            !! pressure [Pa] from the CAM model
        real(8), dimension(:,:), intent(in) :: pres_int_cam
            !! interface pressure [Pa] from the CAM model
        real(8), dimension(:), intent(in)   :: pres_sfc_cam
            !! surface pressure [Pa] from the CAM model
        real(8), dimension(:,:), intent(in) :: tabs_cam
        !! absolute temperature [K] from the CAM model        
        real(8), dimension(:,:), intent(in) :: qv_cam, qc_cam, qi_cam
            !! moisture content [kg/kg] from the CAM model
        real(8), dimension(:,:), intent(in) :: zm_cam
        real(8), dimension(:), intent(in) :: phis_cam
        real(8), dimension(:), intent(in)   :: landfrac, ps_cam
            !! surface pressure [Pa] from the CAM model

        
        real(8), dimension(:, :), intent(out) :: dqi, dqv, dqc, ds
            !! Output tendencies - CAM variables on the CAM grid
        real(8), dimension(:), intent(out) :: precsfc
            !! Output surface precipitation in CAM

        !! XL XL XL
        real(8), dimension(:,:), intent(out) :: q_delta_adv, q_delta_auto, q_delta_sed, t_delta_adv, t_delta_auto, t_delta_sed, t_delta_rad
        real(8), dimension(:,:), intent(out) :: t_flux_adv

        real(8), allocatable, dimension(:) :: ps0_cam
        real(8), allocatable, dimension(:,:) :: gamaz_cam
        real(8), allocatable, dimension(:,:) :: qi0_cam, qv0_cam, qc0_cam, tabs0_cam, q0_cam, t0_cam
        real(8), allocatable, dimension(:,:) :: qi1_cam, qv1_cam, qc1_cam, tabs1_cam, q1_cam, t1_cam
        real(8), allocatable, dimension(:,:) :: delta_t_cam, delta_q_cam

        real(8), allocatable, dimension(:, :) :: t_rad_rest_tend, q_flux_adv, q_tend_auto, q_sed_flux
        real(8), allocatable, dimension(:, :) :: t_flux_adv_f, q_flux_adv_f, q_sed_flux_f
        
        ! Local input variables on the SAM grid
        real(8), dimension(ncol, nrf) :: qi0_sam, qv0_sam, qc0_sam, tabs0_sam, q0_sam, t0_sam
        !! Variables at the start of the parameterisation on the SAM grid - used to calculate tendencies
        real(8), dimension(ncol) :: ps0_sam

        !! XL XL XL
        real(8), dimension(ncol, nrf) :: t_rad_rest_tend_sam, t_flux_adv_sam, q_flux_adv_sam, q_tend_auto_sam, q_sed_flux_sam
        real(8), dimension(ncol, nrf) :: t_flux_adv_sam_f, q_flux_adv_sam_f, q_sed_flux_sam_f

        real(8), dimension(ncol, nrf) :: q_delta_adv_sam, q_delta_sed_sam, t_delta_adv_sam, t_delta_sed_sam
        real(8), dimension(ncol, nrf) :: gamaz_sam
        
        !! CAM variable tendencies calculated on the SAM grid
        real(8), dimension(ncol)      :: qi0_surf, qv0_surf, qc0_surf, tabs0_surf, q0_surf
            !! Surface values
        real(8) :: y_in(ncol)
            !! Distance of column from equator (proxy for insolation and sfc albedo)

        real(8), dimension(:), allocatable :: landmask_1D
        real(8), dimension(:,:), allocatable :: landmask_2D  
        
        integer :: k, i
        integer, dimension(:), allocatable :: idx
        integer :: pver

        ! XL        
        pver = size(tabs_cam, 2)        
        allocate(gamaz_cam(ncol, pver))        
        
        allocate(q0_cam(ncol, pver))
        allocate(t0_cam(ncol, pver))        
        allocate(tabs0_cam(ncol, pver))
        allocate(qv0_cam(ncol, pver))        
        allocate(qc0_cam(ncol, pver))
        allocate(qi0_cam(ncol, pver))

        allocate(q1_cam(ncol, pver))
        allocate(t1_cam(ncol, pver))        
        allocate(tabs1_cam(ncol, pver))
        allocate(qv1_cam(ncol, pver))
        allocate(qc1_cam(ncol, pver))
        allocate(qi1_cam(ncol, pver))

        allocate(delta_t_cam(ncol, pver))
        allocate(delta_q_cam(ncol, pver))

        allocate(t_rad_rest_tend(ncol, pver))
        allocate(q_flux_adv(ncol, pver))
        allocate(q_sed_flux(ncol, pver))        
        allocate(q_tend_auto(ncol, pver))

        allocate(t_flux_adv_f(ncol, pver))
        allocate(q_flux_adv_f(ncol, pver))
        allocate(q_sed_flux_f(ncol, pver))
        
        allocate(landmask_1D(ncol))
        allocate(landmask_2D(ncol, pver))
        ! XL        
        
        ! Initialise precipitation to 0 if required and at start of cycle if subcycling
        precsfc(:)=0.

        do i = 1,ncol
           do k = 1,pver
              !              gamaz_cam(i,k) = zm_cam(i,k)*ggr/cp_cam 
              gamaz_cam(i,k) = zm_cam(i,k)*ggr/cp_cam + phis_cam(i)/cp_cam
           end do
        end do
        ! cpairv
        
        ! NOTE: tabs0_cam/qv0_cam/qc0_cam/qi0_cam are allocated to the actual
        ! ncol extent, while tabs_cam/qv_cam/qc_cam/qi_cam are assumed-shape
        ! dummy args bound to the caller's full-pcols-extent slices. Must be
        ! explicitly bounded to ncol on both sides to avoid a shape-mismatched
        ! whole-array assignment writing past the allocated extent (same class
        ! of bug as the dqv/dqc/dqi/ds fix below).
        tabs0_cam(1:ncol,:) = tabs_cam(1:ncol,:)
        qv0_cam(1:ncol,:)   = qv_cam(1:ncol,:)
        qc0_cam(1:ncol,:)   = qc_cam(1:ncol,:)
        qi0_cam(1:ncol,:)   = qi_cam(1:ncol,:)

        call  CAM_var_conversion(qv0_cam, qc0_cam, qi0_cam, gamaz_cam, q0_cam, tabs0_cam, t0_cam)
        ! Store variables on SAM grid at the start of the timestep
        ! for calculating tendencies later        
        !-----------------------------------------------------

        ! Set surface values
        call extrapolate_to_surface(pres_cam(1:ncol, :), tabs0_cam(1:ncol, :), pres_sfc_cam(1:ncol), tabs0_surf)
        call extrapolate_to_surface(pres_cam(1:ncol, :), q0_cam(1:ncol, :), pres_sfc_cam(1:ncol), q0_surf)        
        where (q0_surf>=0)
            q0_surf = q0_surf
        elsewhere
            q0_surf = 0
        end where
        call interp_to_sam(pres_cam(1:ncol, :), pres_sfc_cam(1:ncol), &
                           tabs0_cam(1:ncol, :), tabs0_sam, tabs0_surf)
        call interp_to_sam(pres_cam(1:ncol, :), pres_sfc_cam(1:ncol), &
                           q0_cam(1:ncol, :), q0_sam, q0_surf)

        ps0_sam=ps_cam(1:ncol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        ! Run the neural net parameterisation
        ! Updates q and t with delta values
        ! advective, autoconversion (dt = -dq*(latent_heat/cp)),
        ! sedimentation (dt = -dq*(latent_heat/cp)),
        ! radiation rest tendency (multiply by dtn to get dt)
        
        call nn_convection_flux(tabs0_sam, q0_sam, ps0_sam, landfrac, &
             rho, adz, dz, dtn, &
             t_rad_rest_tend_sam, t_flux_adv_sam_f, &
             q_flux_adv_sam_f, q_tend_auto_sam, q_sed_flux_sam_f)

        ! We interpolate fluxes at lower level interface from SAM to CAM grid
        call interp_flux_sam_to_cam(pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), t_flux_adv_sam_f, t_flux_adv_f)
        call interp_flux_sam_to_cam(pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), q_flux_adv_sam_f, q_flux_adv_f)
        call interp_flux_sam_to_cam(pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), q_sed_flux_sam_f, q_sed_flux_f)

        call interp_flux_cam_midpoint(pres_cam(1:ncol, :), pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), t_flux_adv_f, t_flux_adv)
        call interp_flux_cam_midpoint(pres_cam(1:ncol, :), pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), q_flux_adv_f, q_flux_adv)

        ! Interpolate SAM variables to the CAM pressure levels
        ! setting tendencies above the SAM grid to 0.0
        ! TODO Interpolate all variables in one call
        call interp_to_cam(pres_cam(1:ncol, :), pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), t_rad_rest_tend_sam, t_rad_rest_tend(1:ncol, :))
        call interp_to_cam(pres_cam(1:ncol, :), pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), q_tend_auto_sam, q_tend_auto(1:ncol, :))
        call update_tendencies(dtn, tabs0_cam, q0_cam, t0_cam, pres_int_cam, t_rad_rest_tend, t_flux_adv_f, q_flux_adv_f, q_tend_auto, q_sed_flux_f, t1_cam, t_delta_adv, t_delta_auto, t_delta_sed, t_delta_rad, q1_cam, q_delta_adv, q_delta_auto, q_delta_sed, precsfc)
                
        ! Convert precipitation from kg/m^2 to m by dividing by density (1000)
        precsfc = precsfc * 1.0D-3
        
        !-----------------------------------------------------
        ! Formulate the output variables to CAM as required. !
        call SAM_var_conversion(t0_cam, q0_cam, tabs0_cam, qv0_cam, qc0_cam, qi0_cam, gamaz_cam, pres_cam, t1_cam, q1_cam, tabs1_cam, qv1_cam, qc1_cam, qi1_cam)
        
        !-----------------------------------------------------

        ! Convert back into CAM variable tendencies (diff div by dtn) on SAM grid
        ! q into components and tabs to s (mult by cp)
        ! NOTE: dqv/dqc/dqi/ds are assumed-shape, bound to the caller's
        ! full-pcols-extent ptend%q(:,...) slices, while qv1_cam/qv0_cam etc.
        ! are allocated to the actual ncol extent. ncol < pcols for any
        ! partial/edge chunk, so these must be explicitly bounded to ncol on
        ! both sides to avoid a shape-mismatched whole-array assignment.
        dqv(1:ncol,:) = (qv1_cam(1:ncol,:) - qv0_cam(1:ncol,:)) / dtn
        dqc(1:ncol,:) = (qc1_cam(1:ncol,:) - qc0_cam(1:ncol,:)) / dtn
        dqi(1:ncol,:) = (qi1_cam(1:ncol,:) - qi0_cam(1:ncol,:)) / dtn
        ds(1:ncol,:)  = cp_cam * (tabs1_cam(1:ncol,:) - tabs0_cam(1:ncol,:)) / dtn
        ! Convert precipitation from total in m to date in m / s by div by dtn
        precsfc = precsfc / dtn
!        write(iulog,*)'precsfc',precsfc        
        !-----------------------------------------------------

        q_delta_adv(1:ncol, :)=q_delta_adv(1:ncol, :)/dtn
        q_delta_auto(1:ncol, :)=q_delta_auto(1:ncol, :)/dtn
        q_delta_sed(1:ncol, :)=q_delta_sed(1:ncol, :)/dtn
        t_delta_adv(1:ncol, :)=t_delta_adv(1:ncol, :)/dtn
        t_delta_auto(1:ncol, :)=t_delta_auto(1:ncol, :)/dtn
        t_delta_sed(1:ncol, :)=t_delta_sed(1:ncol, :)/dtn
        t_delta_rad(1:ncol, :)=t_delta_rad(1:ncol, :)/dtn        

!        !define landmask: 0 for land - 1 for ocean (convention reversed from landfrac) - any land fraction > 0 defines land 
!        landmask_1D(1:ncol)=1.0
!        landmask_2D(1:ncol,:)=1.0
!        do i=1,ncol      
!           if (landfrac(i) .gt. 0.0) then
!              landmask_1D(i)=0.0
!           end if
!        end do
!        do k=1,pver
!           landmask_2D(1:ncol,k)=landmask_1D(1:ncol)
!        end do

!        dqi(:ncol,:pver)=dqi(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        dqv(:ncol,:pver)=dqv(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        dqc(:ncol,:pver)=dqc(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        ds(:ncol,:pver)=ds(:ncol,:pver)*landmask_2D(:ncol,:pver)
        
!        precsfc(:ncol)=precsfc(:ncol)*landmask_1D(:ncol)   
   
!        q_delta_adv(:ncol,:pver)=q_delta_adv(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        q_delta_auto(:ncol,:pver)=q_delta_auto(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        q_delta_sed(:ncol,:pver)=q_delta_sed(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        t_delta_adv(:ncol,:pver)=t_delta_adv(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        t_delta_auto(:ncol,:pver)=t_delta_auto(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        t_delta_sed(:ncol,:pver)=t_delta_sed(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        t_delta_rad(:ncol,:pver)=t_delta_rad(:ncol,:pver)*landmask_2D(:ncol,:pver)
!        t_flux_adv(:ncol,:pver)=t_flux_adv(:ncol,:pver)*landmask_2D(:ncol,:pver)        
        
        ! deallocate
        deallocate(gamaz_cam)
        deallocate(q0_cam)
        deallocate(t0_cam)
        deallocate(tabs0_cam)        
        deallocate(qi0_cam)
        deallocate(qv0_cam)
        deallocate(qc0_cam)


        deallocate(q1_cam)
        deallocate(t1_cam)        
        deallocate(tabs1_cam)
        deallocate(qi1_cam)        
        deallocate(qv1_cam)
        deallocate(qc1_cam)

        deallocate(delta_t_cam)
        deallocate(delta_q_cam)

        deallocate(t_rad_rest_tend)
        deallocate(q_tend_auto)

        deallocate(landmask_1D)
        deallocate(landmask_2D)        
        
    end subroutine nn_convection_flux_CAM


    subroutine nn_convection_flux_CAM_finalize()
        !! Finalize the module deallocating arrays

        call nn_convection_flux_finalize()
        call sam_sounding_finalize()

    end subroutine nn_convection_flux_CAM_finalize


    !-----------------------------------------------------------------
    ! Private Subroutines

    ! ---------------------
    subroutine update_tendencies(dtn, tabs, q_in, t_in, pres_int, t_rad_rest_tend, &
         t_flux_adv, q_flux_adv, q_tend_auto, q_sed_flux, &
         t_out, t_delta_adv, t_delta_auto, t_delta_sed, t_delta_rad, &
         q_out, q_delta_adv, q_delta_auto, q_delta_sed, precsfc)
      
      ! all fluxes are defined on half levels (pres_int)
      
      integer :: ncol, nrf
      !! array sizes
      integer :: i, k

      real(8) :: fmaxfac
      
      real(8), intent(in) :: dtn    ! Seconds            
      real(8), dimension(:,:), intent(in) :: tabs, q_in, t_in
      real(8), dimension(:,:), intent(in) :: pres_int
      real(8), dimension(:), intent(inout) :: precsfc
            !! interface pressure [Pa] from the CAM model      
      real(8), dimension(:,:), intent(inout) :: t_rad_rest_tend, q_tend_auto, t_flux_adv, q_flux_adv, q_sed_flux
      real(8), dimension(:,:), intent(out) :: t_delta_adv, q_delta_adv, q_delta_sed
      real(8), dimension(:,:), intent(out) :: t_out, t_delta_auto, t_delta_sed, t_delta_rad
      real(8), dimension(:,:), intent(out) :: q_out, q_delta_auto     

      real(8), allocatable, dimension(:,:) :: irhoadzdz, dp, omp, fac, nn_mask

      ! NOTE: t_delta_adv is bound to the caller's full-pcols-extent slice,
      ! while tabs/q_in/t_in (and t_out/q_out) are correctly sized to the
      ! real ncol. Using size(t_delta_adv,1) here returned pcols instead of
      ! ncol, causing an out-of-bounds read on t_in/q_in/tabs for any
      ! partial/edge chunk (DEBUG=TRUE caught this: "Subscript #1 of T_IN
      ! has value 16 which is greater than the upper bound of 15").
      ncol = size(tabs, 1)
      nrf = size(tabs, 2)
      allocate(irhoadzdz(ncol,nrf))
      allocate(dp(ncol,nrf))
      allocate(omp(ncol,nrf))
      allocate(fac(ncol,nrf))
      allocate(nn_mask(ncol,nrf))
      
      do i=1,ncol
         ! Initialize variables
         t_out(i,:) = t_in(i,:)
         q_out(i,:) = q_in(i,:)
         t_delta_auto(i,:) = 0.
         q_delta_auto(i,:) = 0.
         t_delta_sed(i,:) = 0.
         t_delta_rad(i,:) = 0.
         t_delta_adv(i,:) = 0.
         q_delta_adv(i,:) = 0.
         q_delta_sed(i,:) = 0.
         
         omp(i,:) = 0.
         fac(i,:) = 0.
         precsfc(i) = 0.
         
         !-----------------------------------------------------
         ! Apply physical constraints and update q and t
         do k=1,nrf
            dp(i,k) = abs(pres_int(i, k+1) - pres_int(i, k))
            if (dp(i,k) > 0.0) then
               irhoadzdz(i,k) = ggr/dp(i,k)*dtn
            else
               irhoadzdz(i,k) = 0.0
            end if
         end do
         
         !-----------------------------------------------------
         ! Non-precip. water content must be >= 0, so ensure advective fluxes
         ! will not reduce it below 0 anywhere
         
         do k=1,nrf-1               
            t_delta_adv(i,k) = - (t_flux_adv(i,k+1) - t_flux_adv(i,k)) * irhoadzdz(i,k)
            q_delta_adv(i,k) = - (q_flux_adv(i,k+1) - q_flux_adv(i,k)) * irhoadzdz(i,k)
         end do
         ! Enforce boundary condition at top of column
         t_delta_adv(i,nrf) = - (0.0 - t_flux_adv(i,nrf)) * irhoadzdz(i,nrf)
         q_delta_adv(i,nrf) = - (0.0 - q_flux_adv(i,nrf)) * irhoadzdz(i,nrf)

         ! Update q and t with delta values                                                                                         
         q_out(i,1:nrf) = q_out(i,1:nrf) + q_delta_adv(i,1:nrf)
         t_out(i,1:nrf) = t_out(i,1:nrf) + t_delta_adv(i,1:nrf)        
                           
         ! Convert sedimenting fluxes to deltas
         do k=1,nrf-1
            q_delta_sed(i,k) = - (q_sed_flux(i,k+1) - q_sed_flux(i,k)) * irhoadzdz(i,k)
         end do
         ! Enforce boundary condition at top of column
         q_delta_sed(i,nrf) = - (0.0 - q_sed_flux(i,nrf)) * irhoadzdz(i,nrf)                    

         ! Compute sedimentation on dt (dt = -dq*(latent_heat/cp))
         t_delta_sed(i,1:nrf) = - q_delta_sed(i,1:nrf)*(fac_fus+fac_cond)
         ! Update q and t with delta values
         q_out(i,1:nrf)=q_out(i,1:nrf)+q_delta_sed(i,1:nrf)         
         t_out(i,1:nrf)=t_out(i,1:nrf)+t_delta_sed(i,1:nrf)

         ! ensure autoconversion tendency won't reduce q below 0
         do k=1,nrf
            omp(i,k) = max(0.,min(1.,(tabs(i,k)-tprmin)*a_pr))
            fac(i,k) = (fac_cond + fac_fus * (1.0 - omp(i,k)))
            q_delta_auto(i,k) = q_tend_auto(i,k) * dtn
         end do      
         ! Compute autoconversion on dt (dt = -dq*(latent_heat/cp))
         t_delta_auto(i,1:nrf) = - q_delta_auto(i,1:nrf)*fac(i,1:nrf)
         ! Update q and t with delta values 
         q_out(i,1:nrf) = q_out(i,1:nrf) + q_delta_auto(i,1:nrf)
         t_out(i,1:nrf) = t_out(i,1:nrf) + t_delta_auto(i,1:nrf)
         
         ! Apply radiation rest tendency (multiply by dtn to get delta value)
         t_delta_rad(i,1:nrf) = t_rad_rest_tend(i,1:nrf)*dtn
         ! Update t with delta values 
         t_out(i,1:nrf) = t_out(i,1:nrf) + t_delta_rad(i,1:nrf)                  

         ! As a final check enforce q must be >= 0.0 - reset q_out to q_in if q_out <0
         ! This reset is a source of non-precipitating water, which is assumed to be from autoconversion
         ! (e.g. droplet evaporation, ice sublimation, re-entrainment of falling condensate, breakup/collision droplets/crystals)
         ! This added autoconversion implies a decrease in specific energy ('cooling') locally.

         ! Note 1: (negative) autoconversion terms sumed over the column cannot be less than the positive source of moisture added by this reset - if so, we end up with net positive column autoconversion (negative precipitation) - if precipitation is negative, I reset the entire column by setting all tendencies to 0.
         ! Note 2: we do not consider moisture reset to affect either advection or sedimentation terms. Doing so would lead to column water not being conserved despite advection acting conservatively.
         do k = 1,nrf
            if ( q_out(i,k) < 0.0 ) then
!               write(iulog,*)'Negative q at level ',k, ' - total q tendency reset to 0 by autoconversion'
               q_delta_auto(i,k) = q_delta_auto(i,k) + (q_in(i,k) - q_out(i,k))
               t_out(i,k) = t_out(i,k) - fac(i,k) * (q_in(i,k) - q_out(i,k))
               q_out(i,k) = q_in(i,k)
            end if
         end do
         
         ! Compute net water loss in column (precipitation)
         do k=1,nrf
            precsfc(i) = precsfc(i) - (q_out(i,k) - q_in(i,k))* dp(i,k) / ggr                       
         end do

         if ( precsfc(i) < 0.0 ) then
!            write(iulog,*)'Negative Precipitation - All tendencies set to 0 in column'            
            q_out(i,1:nrf) = q_in(i,1:nrf)
            t_out(i,1:nrf) = t_in(i,1:nrf)
            q_delta_auto(i,1:nrf) = 0.0
            t_delta_auto(i,1:nrf) = 0.0            
            q_delta_adv(i,1:nrf) = 0.0
            t_delta_adv(i,1:nrf) = 0.0            
            q_delta_sed(i,1:nrf) = 0.0
            t_delta_sed(i,1:nrf) = 0.0            
            t_delta_rad(i,1:nrf) = 0.0
            precsfc(i) = 0.0
         end if
                  
      end do

      ! deallocate
      deallocate(irhoadzdz)
      deallocate(dp)
      deallocate(omp)
      deallocate(fac)
      
    end subroutine update_tendencies

    subroutine interp_flux_sam_to_cam(pres_int_cam, pres_sfc_cam, flux_adv_sam, flux_adv_cam)

      real(8), dimension(:,:), intent(in) :: pres_int_cam
      real(8), dimension(:), intent(in) :: pres_sfc_cam
      real(8), dimension(:,:), intent(in) :: flux_adv_sam
      real(8), dimension(:,:), intent(out) :: flux_adv_cam

      real(8), dimension(:,:), allocatable :: p_norm_cam
      real(8), dimension(:), allocatable :: p_norm_sam
      
      integer :: i, k, ncol, nz_sam, nz_cam, k_l, k_u

      ncol = size(flux_adv_sam, 1)
      nz_sam = size(flux_adv_sam, 2)
      ! NOTE: pres_int_cam is an interface-level array (pver+1 levels), but
      ! flux_adv_cam (the array actually being written below) is a mid-level
      ! array with only pver levels. Bounding nz_cam by pres_int_cam's extent
      ! (DEBUG=TRUE caught this: "Subscript #2 of FLUX_ADV_CAM has value 33
      ! which is greater than the upper bound of 32") wrote one level past
      ! flux_adv_cam's allocated extent. Use flux_adv_cam's own extent instead.
      nz_cam = size(flux_adv_cam, 2)

      allocate(p_norm_sam(nz_sam))
      allocate(p_norm_cam(ncol,nz_cam))
      
      p_norm_sam(:) = presi(:) / presi(1)
      do i = 1,ncol
         do k=1,nz_cam
            p_norm_cam(i,k) = pres_int_cam(i,k) / pres_sfc_cam(i)
         end do
      end do

      do i = 1,ncol
         flux_adv_cam(i,1)=flux_adv_sam(i,1)
         do k=2,nz_cam

            ! assess whether sam levels exists at cam level
            if (p_norm_sam(nz_sam) .gt. p_norm_cam(i,k)) then
               flux_adv_cam(i,k) = 0.0
            else               
               ! find sam level just below cam level: k_l
               k_l = 1
               do while (p_norm_sam(k_l) .ge. p_norm_cam(i,k) .and. k_l < nz_sam)
                  k_l = k_l + 1
               end do
               k_l = k_l - 1
               k_u = k_l + 1
               
               if (p_norm_sam(k_l)==p_norm_cam(i,k)) then
                  flux_adv_cam(i,k) = flux_adv_sam(i,k_l)
               else                                 
                  ! linearly interpolate
                  flux_adv_cam(i,k) = flux_adv_sam(i,k_l) + (flux_adv_sam(i,k_u)-flux_adv_sam(i,k_l))*(p_norm_cam(i,k)-p_norm_sam(k_l))/(p_norm_sam(k_u)-p_norm_sam(k_l))
               end if

            end if

         end do
      end do
      
      deallocate(p_norm_sam)
      deallocate(p_norm_cam)      
      
    end subroutine interp_flux_sam_to_cam
          
    subroutine interp_to_sam(p_cam, p_surf_cam, var_cam, var_sam, var_cam_surface)
        !! Interpolate from the CAM pressure grid to the SAM pressure grid.
        !! Uses linear interpolation between nearest grid points on domain (CAM) grid.

        !! The CAM and SAM grids have pressures measured at cell centres which is where
        !! the variables that we are interpolating are stored.

        !! It can be expected, but not guaranteed, that the domain (CAM) grid is
        !! coarser than the target (SAM) grid.

        !! NOTE: When reading this code remember that pressure DECREASES monotonically
        !!       with altitude, so comparisons are not intuitive: Pa > Pb => a is lower altitude than b

        != unit hPa :: p_cam, p_surf_cam
        real(8), dimension(:,:), intent(in) :: p_cam
        real(8), dimension(:), intent(in) :: p_surf_cam
            !! pressure [Pa] from the CAM model
        != unit 1 :: var_cam, var_sam
        real(8), dimension(:,:), intent(in) :: var_cam
            !! variable from CAM grid to be interpolated from
        real(8), dimension(:), intent(in) :: var_cam_surface
            !! Surface value of CAM variable for interpolation/boundary condition
        real(8), dimension(:,:), intent(out) :: var_sam
            !! variable from SAM grid to be interpolated to
        != unit 1 :: p_norm_sam, p_norm_cam
        real(8), dimension(nrf) :: p_mid_sam, p_norm_sam
        real(8), dimension(:,:), allocatable :: p_norm_cam
            !! Normalised pressures (using surface pressure) as a sigma-coordinate

        integer :: i, k, nz_cam, ncol_cam

        integer :: kc_u, kc_l
        real(8) :: pc_u, pc_l
            !! Pressures of the CAM cell above and below the target SAM cell

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(ncol_cam, nz_cam))

        ! Normalise pressures by surface value
        ! Note - We work in pressure space but use a pressure 'midpoint' based on altitude from
        !        SAM as this is concurrent to the variable value at this location.
        p_norm_sam(:) = pres(:) / presi(1)
        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / p_surf_cam(:)
        end do

        ! Loop over columns and SAM levels and interpolate each column
        do k = 1, nrf
        do i = 1, ncol_cam
            ! Check we are within array

            ! If SAM cell centre is below lowest CAM centre - interpolate to surface value
            if (p_norm_sam(k) > p_norm_cam(i, 1)) then
!                write(*,*) "Interpolating to surface."
                var_sam(i, k) = var_cam_surface(i) &
                                + (p_norm_sam(k)-1.0) &
                                * (var_cam(i, 1)-var_cam_surface(i))/(p_norm_cam(i, 1)-1.0)
            ! Check CAM grid top cell is above SAM grid top
            elseif (p_norm_cam(i, nz_cam) > p_norm_sam(nrf)) then
                ! This should not happen as CAM grid extends to higher altitudes than SAM
                ! TODO This check is run on every iteration - move outside
!                 write(*,*) "CAM upper pressure level is lower than that of SAM: Stopping."
                stop
            else
                ! Locate the neighbouring CAM indices to interpolate between
                ! TODO - this will be slow - speed up later
                kc_u = 1
                pc_u = p_norm_cam(i, kc_u)
                do while (pc_u > p_norm_sam(k))
                    kc_u = kc_u + 1
                    pc_u = p_norm_cam(i, kc_u)
                end do
                kc_l = kc_u - 1
                pc_l = p_norm_cam(i, kc_l)
                ! Redundant following the lines above
                ! do while (pc_l < p_norm_sam(k))
                !     kc_l = kc_l - 1
                !     pc_l = p_norm_cam(i, kc_l)
                ! end do

                ! interpolate variables - Repeat for as many variables as need interpolating
                var_sam(i, k) = var_cam(i, kc_l) &
                                + (p_norm_sam(k)-pc_l) &
                                * (var_cam(i, kc_u)-var_cam(i, kc_l))/(pc_u-pc_l)
            endif
        end do
        end do

        ! Clean up
        deallocate(p_norm_cam)

    end subroutine interp_to_sam

    subroutine interp_to_cam(p_cam, p_int_cam, p_surf_cam, var_sam, var_cam)
      ! TODO This would seriously benefit from a refactor to add the top interface to p_sam_int and p_cam_int
        !! Interpolate from the SAM NN pressure grid to the CAM pressure grid.
        !! Uses conservative regridding between grids.

        !! Requires dealing in interfaces rather than centres.
        !! These  start at the surface of both grids
        !! and extend to grid-top for CAM, but stop at the base of the top cell
        !! (i.e. no grid-top value) for SAM.

        != unit Pa :: p_cam, p_int_cam
        real(8), dimension(:,:), intent(in) :: p_cam
            !! pressure [Pa] from the CAM model
        real(8), dimension(:,:), intent(in) :: p_int_cam
            !! interface pressures [Pa] for each column in CAM
        real(8), dimension(:), intent(in) :: p_surf_cam
            !! surface pressures [Pa] for each column in CAM
        != unit 1 :: var_cam, var_sam
        real(8), dimension(:,:), intent(out) :: var_cam
            !! variable from CAM grid to be interpolated from
        real(8), dimension(:,:), intent(in) :: var_sam
            !! variable from SAM grid to be interpolated to
        != unit 1 :: p_norm_sam, p_norm_cam
        real(8), dimension(nrf) :: p_norm_sam
        real(8), dimension(nrf+1) :: p_int_norm_sam
        real(8), dimension(:,:), allocatable :: p_norm_cam
        real(8), dimension(:,:), allocatable :: p_int_norm_cam

        integer :: i, k, c, nz_cam, ncol_cam, ks, ks_u, ks_l
        real(8) :: ps_u, ps_l, pc_u, pc_l
        real(8) :: p_norm_sam_min

        ncol_cam = size(p_cam, 1)
        nz_cam = size(p_cam, 2)
        allocate(p_norm_cam(ncol_cam, nz_cam))
        allocate(p_int_norm_cam(ncol_cam, nz_cam+1))

        ! Normalise pressures by surface value (extending SAM to add top interface)
        p_norm_sam(:) = pres(1:nrf) / presi(1)
        p_int_norm_sam(1:nrf+1) = presi(1:nrf+1) / presi(1)
        p_int_norm_sam(1) = 1.0

        do k = 1,nz_cam
            p_norm_cam(:,k) = p_cam(:,k) / p_surf_cam(:)
            p_int_norm_cam(:,k) = p_int_cam(:,k) / p_surf_cam(:)
        end do
        p_int_norm_cam(:,nz_cam+1) = p_int_cam(:,nz_cam+1) / p_surf_cam(:)
        p_int_norm_cam(:,1) = 1.0

        ! Loop over columns and SAM levels and interpolate each column
        ! TODO Revert indices to i,k for speed
        do k = 1, nz_cam
        do i = 1, ncol_cam
            ! Check we are within array overlap:
            ! No need to compare bottom end as pressure is normalised to 1.0 for both.

            if (p_int_norm_cam(i, k) < p_int_norm_sam(nrf+1)) then
                ! Anything above the SAM NN grid should be set to 0.0
!                 write(*,*) "CAM pressure is lower than SAM top: set tend to 0.0."
                var_cam(i, k) = 0.0

            else
                ! Get the pressures at the top and bottom of the CAM cell
                pc_u = p_int_norm_cam(i, k+1)
                pc_l = p_int_norm_cam(i,k)

                ! Locate the SAM indices to interpolate between
                ! This will be slow - speed up later
                ! ks_u is the SAM block with an upper interface above the CAM block
                ks_u = 1
                ps_u = p_int_norm_sam(ks_u+1)
                do while ((ps_u > pc_u) .and. (ks_u < nrf))
                    ks_u = ks_u + 1
                    ps_u = p_int_norm_sam(ks_u+1)                    
                end do
                ! ks_l is the SAM block with a lower interface below the CAM block
                ks_l = ks_u
                ps_l = p_int_norm_sam(ks_l)
                do while ((ps_l < pc_l) .and. (ks_l > 1))
                    ks_l = ks_l - 1
                    ps_l = p_int_norm_sam(ks_l)
                end do

                ! Combine all SAM cells that this CAM_CELL(i, k) overlaps
                var_cam(i,k) = 0.0
                if (ks_u == ks_l) then
                    ! CAM cell fully enclosed by SAM cell => match 'density' of SAM cell
                    var_cam(i,k) = var_cam(i,k) + var_sam(i,ks_u)
                else
                  do c = ks_l, ks_u
                    if (c == ks_l) then
                      ! Bottom cell
                      var_cam(i,k) = var_cam(i,k) + (pc_l-p_int_norm_sam(ks_l+1))*var_sam(i,ks_l) / (pc_l - pc_u)
                    elseif (c == ks_u) then
                      ! Top cell
                      var_cam(i,k) = var_cam(i,k) + (p_int_norm_sam(ks_u)-pc_u)*var_sam(i, ks_u) / (pc_l - pc_u)
                    else
                      ! Intermediate cell (SAM Cell fully enclosed by CAM Cell => Absorb all)
                      var_cam(i,k) = var_cam(i,k) + (p_int_norm_sam(c)-p_int_norm_sam(c+1))*var_sam(i, c) / (pc_l - pc_u)

                    endif

                  enddo
                endif
            endif
        end do
        end do

        ! Clean up
        deallocate(p_norm_cam)
        deallocate(p_int_norm_cam)

      end subroutine interp_to_cam

    subroutine extrapolate_to_surface(p_cam, var_cam, p_surf_cam, var_surf_cam)
        !! Extrapolate variable profile to the surface, required for interpolation.

        != unit Pa :: p_cam
        real(8), dimension(:,:), intent(in) :: p_cam
            !! pressure [Pa] from the CAM model
        real(8), dimension(:), intent(in) :: p_surf_cam
            !! surface pressures [Pa] for each column in CAM
        != unit 1 :: var_cam
        real(8), dimension(:,:), intent(in) :: var_cam
            !! variable from CAM grid to be interpolated from
        real(8), dimension(:), intent(out) :: var_surf_cam
            !! variable from CAM grid to be interpolated from

        real(8) :: gdt(size(var_surf_cam))

        gdt(:) = (var_cam(:,2) - var_cam(:,1)) / (p_cam(:,2) - p_cam(:,1))
        var_surf_cam(:) = var_cam(:,1) + gdt(:) * (p_surf_cam(:) - p_cam(:,1))

    end subroutine extrapolate_to_surface

    subroutine sam_sounding_init(filename)
        !! Read various profiles in from SAM sounding file

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, k
        integer :: z_dimid, dz_dimid
        integer :: z_varid, dz_varid, pres_varid, presi_varid, rho_varid, adz_varid

        character(len=136), intent(in) :: filename
            !! NetCDF filename from which to read data


        !-------------allocate arrays and read data-------------------

        ! Open the file. NF90_NOWRITE tells netCDF we want read-only
        ! access
        ! Get the varid or dimid for each variable or dimension based
        ! on its name.

        call check( nf90_open(trim(filename),NF90_NOWRITE,ncid ))

        call check( nf90_inq_dimid(ncid, 'z', z_dimid))
        call check( nf90_inquire_dimension(ncid, z_dimid, len=nz_sam))
        ! call check( nf90_inq_dimid(ncid, 'dz', dz_dimid))
        ! call check( nf90_inquire_dimension(ncid, dz_dimid, len=ndz))

        ! Note that nz of sounding may be longer than nrf
        allocate(z(nz_sam))
        allocate(pres(nz_sam))
        allocate(presi(nz_sam))
        allocate(rho(nz_sam))
        allocate(adz(nz_sam))
        allocate(gamaz(nz_sam))

        ! Read data in from nc file - convert pressures to hPa
        call check( nf90_inq_varid(ncid, "z", z_varid))
        call check( nf90_get_var(ncid, z_varid, z))
        call check( nf90_inq_varid(ncid, "pressure", pres_varid))
        call check( nf90_get_var(ncid, pres_varid, pres))
        pres(:) = pres / 100.0
        call check( nf90_inq_varid(ncid, "interface_pressure", presi_varid))
        call check( nf90_get_var(ncid, presi_varid, presi))
        presi(:) = presi(:) / 100.0
        call check( nf90_inq_varid(ncid, "rho", rho_varid))
        call check( nf90_get_var(ncid, rho_varid, rho))
        call check( nf90_inq_varid(ncid, "adz", adz_varid))
        call check( nf90_get_var(ncid, adz_varid, adz))
        call check( nf90_inq_varid(ncid, "dz", dz_varid))
        call check( nf90_get_var(ncid, dz_varid, dz))

        ! Close the nc file
        call check( nf90_close(ncid))

        ! Calculate gamaz required elsewhere
        do k = 1, nz_sam
            gamaz(k) = ggr/cp*z(k)
        end do

        write(*,*) 'Finished reading SAM sounding file.'

    end subroutine sam_sounding_init


    subroutine sam_sounding_finalize()
        !! Deallocate module variables read from sounding

        deallocate(z, pres, presi, rho, adz, gamaz)

    end subroutine sam_sounding_finalize

    
    subroutine CAM_var_conversion(qv, qc, qi, gamaz, r, tabs, t)
        !! Convert CAM moist mixing ratios for species qv, qc, qi to
        !! dry mixing ratio r to used by SAM parameterisation
        !! q is total water qv/c/i is cloud vapor/liquid/ice
        !! Convert CAM absolute temperature to moist static energy t used by SAM
        !! by inverting lines 49-53 of diagnose.f90 from SAM where

        !! WARNING: This routine uses gamaz(k) which is defined on the SAM grid.
        !! OR
        !! The native geopotential height on CAM grid

        integer :: ncol, nz
            !! array sizes
        integer :: i, k
            !! Counters
        real(8) :: omn
            !! intermediate omn factor used in variable conversion
        real(8) :: rv, rc, ri
            !! intermediate r values to hold dry mixing ratios of vapour, cld, ice

        ! ---------------------
        ! Fields from CAM/SAM
        ! ---------------------
        != unit 1 :: r
        real(8), intent(out) :: r(:, :)
            !! Total non-precipitating water mixing ratio as required by SAM NN
        real(8), intent(in) :: qv(:, :)
            !! Cloud water vapour from CAM
        real(8), intent(in) :: qc(:, :)
            !! Cloud water from CAM
        real(8), intent(in) :: qi(:, :)
            !! Cloud ice from CAM
        real(8), intent(in) :: gamaz(:, :)
            !! Geopotential Height

        real(8), intent(out) :: t(:, :)
            !! Static energy as required by SAM NN
        real(8), intent(in) :: tabs(:, :)
            !! Absolute temperature from CAM


        ncol = size(r, 1)
        nz = size(r, 2)

        do k = 1, nz
          do i = 1, ncol
            ! Calculate dry mixing ratio as required by SAM from moist variables as provided by CAM
            rv = qv(i,k) / (1. - qv(i,k))
            rc = qc(i,k) * (1.+rv)
            ri = qi(i,k) * (1.+rv)
            r(i,k) = rv + rc + ri

            ! omp  = max(0.,min(1.,(tabs(i,k)-tprmin)*a_pr))
            omn  = max(0.,min(1.,(tabs(i,k)-tbgmin)*a_bg))
            t(i,k) = tabs(i,k) &
                   ! - (fac_cond+(1.-omp)*fac_fus)*qp(i,k) &  ! There is no qp in CAM
                     - (fac_cond+(1.-omn)*fac_fus)*(rc + ri) &
                     + gamaz(i,k)
          end do
        end do

    end subroutine CAM_var_conversion


    subroutine SAM_var_conversion(t0, r0, tabs0, qv0, qc0, qi0, gamaz, pres_cam, t, r, tabs, qv, qc, qi)
        !! Convert SAM t and r to tabs, qv, qc, qi used by CAM
        !! t is normalised liquid ice static energy, r is a dry mixing ratio
        !! tabs is absolute temperature, q is cloud vapor/liquid/ice,

        integer :: ncol, nz
            !! array sizes
        integer :: i, k, niter
            !! Counters

        ! ---------------------
        ! Fields from SAM/CAM
        ! ---------------------

        real(8), intent(in) :: t0(:, :)
            !! normalised liquid ice static energy (init) XL XL XL       
        real(8), intent(in) :: r0(:, :)
            !! Total non-precipitating water mixing ratio from SAM (init) XL XL XL
        real(8), intent(in) :: tabs0(:, :)
            !! Total non-precipitating water mixing ratio from SAM
        real(8), intent(in) :: qv0(:, :)
            !! Cloud water vapour in CAM
        real(8), intent(in) :: qc0(:, :)
            !! Cloud water (liquid) in CAM
        real(8), intent(in) :: qi0(:, :)
            !! Cloud ice in CAM       
        
        != unit kg/kg :: r, qv, qc, qi
        real(8), intent(in) :: t(:, :)
            !! normalised liquid ice static energy (init) XL XL XL               
        real(8), intent(in) :: r(:, :)
            !! Total non-precipitating water mixing ratio from SAM
        real(8), intent(out) :: qv(:, :)
            !! Cloud water vapour in CAM
        real(8), intent(out) :: qc(:, :)
            !! Cloud water (liquid) in CAM
        real(8), intent(out) :: qi(:, :)
            !! Cloud ice in CAM

        != unit K :: tabs, tabs1
        real(8), intent(out) :: tabs(:, :)
            !! absolute temperature
        real(8) :: tabs1
        !! Temporary variable for tabs
        
        != K :: t
        real(8), intent(in) :: gamaz(:, :)
            !! Geopotential Height

        real(8), intent(in) :: pres_cam(:, :)
        !! Mid-level pressure level

        real(8), allocatable, dimension(:,:) :: pres_m(:, :)
            !! Mid-level pressure level        
        
        ! Intermediate variables
        real(8) :: rsat, om, omn, dtabs, drsat, lstarn, dlstarn, fff, dfff, rn, rv, r_temp
        
        ncol = size(tabs, 1)
        nz = size(tabs, 2)

        allocate(pres_m(ncol,nz))
        
        ! Loop over all cells and iterate until equilibrium of condensed species is achieved.
        ! This code is adapted from cloud.f90 in SAM
        do k = 1, nz
        do i = 1, ncol

            if( (t0(i,k) .eq. t(i,k)) .and. (r0(i,k) .eq. r(i,k)) ) then
               tabs(i,k) = tabs0(i,k)
               qc(i,k) = qc0(i,k)
               qi(i,k) = qi0(i,k)
               qv(i,k) = qv0(i,k)

            else
           
               pres_m(i,k)=pres_cam(i,k)/100.0
               ! pres_m(i,k)=pres(k)
               ! Enforce r >= 0.0
               r_temp=max(0.,r(i,k))

               ! Initial guess for temperature assuming no cloud water/ice:
               tabs(i,k) = t(i,k)-gamaz(i,k)
               tabs1=tabs(i,k)

               ! Set saturation vapour based on cloud type
               ! Warm cloud:
               if(tabs1 >= tbgmax) then
                  rsat = rsatw(tabs1,pres_m(i,k))
                  ! Ice cloud:
               elseif(tabs1 <= tbgmin) then
                  rsat = rsati(tabs1,pres_m(i,k))
                  ! Mixed-phase cloud:
               else
                  om = an*tabs1-bn
                  rsat = om*rsatw(tabs1,pres_m(i,k))+(1.-om)*rsati(tabs1,pres_m(i,k))
               endif

               !  Test if condensation is possible (humidity is above saturation) and iterate:
               if(r_temp > rsat) then
                  niter=0
                  dtabs = 100.
                  !                do while(abs(dtabs) > 0.01 .and. niter < 10)
                  do while(abs(dtabs) > 0.001)                   
                     ! Warm cloud regime
                     if(tabs1 >= tbgmax) then
                        om=1.
                        lstarn=fac_cond
                        dlstarn=0.
                        rsat=rsatw(tabs1,pres_m(i,k))
                        drsat=dtrsatw(tabs1,pres_m(i,k))
                        ! Ice cloud regime
                     else if(tabs1 <= tbgmin) then
                        om=0.
                        lstarn=fac_sub
                        dlstarn=0.
                        rsat=rsati(tabs1,pres_m(i,k))
                        drsat=dtrsati(tabs1,pres_m(i,k))
                        ! Mixed cloud regime
                     else
                        om=an*tabs1-bn
                        lstarn=fac_cond+(1.-om)*fac_fus
                        dlstarn=an
                        rsat=om*rsatw(tabs1,pres_m(i,k))+(1.-om)*rsati(tabs1,pres_m(i,k))
                        drsat=om*dtrsatw(tabs1,pres_m(i,k))+(1.-om)*dtrsati(tabs1,pres_m(i,k))
                     endif

                     ! Update thermodynamics and check for convergence
                     fff = tabs(i,k)-tabs1+lstarn*(r_temp-rsat)
                     dfff=dlstarn*(r_temp-rsat)-lstarn*drsat-1.
                     dtabs=-fff/dfff
                     niter=niter+1
                     tabs1=tabs1+dtabs
                  end do

                  ! Update saturation point and then calculate water and residual vapour content
                  rsat = rsat + drsat * dtabs
                  rn = max(0., r_temp-rsat)
                  rv = max(0., r_temp-rn)

                  ! If condensation not possible rn is 0.0
               else
                  rn = 0.
                  rv = r_temp

               endif

               ! Set tabs to iterated tabs after convection
               tabs(i,k) = tabs1

               ! Code for calculating qc and qi from rn.
               ! Adapted from statistics.f90 in SAM assuming dokruegermicro=.false.
               ! Also implements conversion from SAM dry mixing ratios to CAM moist mixing ratios
               omn = omegan(tabs(i,k))
               qc(i,k) = rn*omn/(1.+rv)
               qi(i,k) = rn*(1.-omn)/(1.+rv)
               qv(i,k) = rv/(1.+rv)

            end if
            
         end do
      end do

      deallocate(pres_m)
        
    end subroutine SAM_var_conversion

    subroutine interp_flux_cam_midpoint(pres_cam, pres_int_cam, pres_sfc_cam, flux_adv_cam_int, flux_adv_cam)
      
      real(8), dimension(:,:), intent(in) :: pres_cam, pres_int_cam, flux_adv_cam_int
      real(8), dimension(:), intent(in) :: pres_sfc_cam
      real(8), dimension(:,:), intent(out) :: flux_adv_cam
      real(8), dimension(:,:), allocatable :: p_norm_cam, p_norm_cam_int
      integer :: i, k, nz, ncol
      
      ncol = size(flux_adv_cam_int, 1)
      nz = size(flux_adv_cam_int, 2)
      allocate(p_norm_cam(ncol, nz))
      allocate(p_norm_cam_int(ncol, nz+1))
      
      do i = 1,ncol
         do k=1, nz
            p_norm_cam(i,k) = pres_cam(i,k) / pres_sfc_cam(i)
            p_norm_cam_int(i,k) = pres_int_cam(i,k) / pres_sfc_cam(i)
         end do
         p_norm_cam_int(i,nz+1) = pres_int_cam(i,nz+1) / pres_sfc_cam(i)
      end do
      
      do i = 1,ncol
         do k=1, nz-1
            flux_adv_cam(i,k) = flux_adv_cam_int(i,k) + (flux_adv_cam_int(i,k+1)-flux_adv_cam_int(i,k))*(p_norm_cam(i,k)-p_norm_cam_int(i,k))/(p_norm_cam_int(i,k+1)-p_norm_cam_int(i,k))
         end do
         flux_adv_cam(i,nz) = flux_adv_cam_int(i,nz) + (0.0-flux_adv_cam_int(i,nz))*(p_norm_cam(i,nz)-p_norm_cam_int(i,nz))/(p_norm_cam_int(i,nz+1)-p_norm_cam_int(i,nz))     
      end do                  

      deallocate(p_norm_cam)
      deallocate(p_norm_cam_int)
      
    end subroutine interp_flux_cam_midpoint
!        ! We interpolate from lower level interface to mid level in SAM grid
!        call interp_sam_lower2mid(t_flux_adv_sam_f, t_flux_adv_sam)
!        call interp_sam_lower2mid(q_flux_adv_sam_f, q_flux_adv_sam)
!        call interp_sam_lower2mid(q_sed_flux_sam_f, q_sed_flux_sam)        
        ! We interpolate from mid level to lower level interface in CAM grid
!        call interp_cam_mid2lower(pres_cam(1:ncol, :), pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), t_flux_adv_sam_f(1:ncol,1), t_flux_adv, t_flux_adv_f)
!        call interp_cam_mid2lower(pres_cam(1:ncol, :), pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), q_flux_adv_sam_f(1:ncol,1), q_flux_adv, q_flux_adv_f)
!        call interp_cam_mid2lower(pres_cam(1:ncol, :), pres_int_cam(1:ncol, :), pres_sfc_cam(1:ncol), q_sed_flux_sam_f(1:ncol,1), q_sed_flux, q_sed_flux_f)       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
end module nn_interface_CAM
