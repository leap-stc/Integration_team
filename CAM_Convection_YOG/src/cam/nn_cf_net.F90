module nn_cf_net_mod
    !! neural_net convection emulator
    !! Module containing code pertaining only to the Neural Net functionalities of
    !! the M2LiNES Convection Flux parameterisation.
    ! Author: J.W. Atkinson
    ! January 2024
    ! Modified to use FTorch libraries: Addisu Semie
    ! November 2025
    !---------------------------------------------------------------------------------

!---------------------------------------------------------------------
! Libraries to use
use netcdf
use shr_kind_mod,   only: r8=>shr_kind_r8
use shr_kind_mod,   only: r4=>shr_kind_r4
use FTorch_cesm_interface, only: torch_kCPU, torch_tensor, torch_model,        &
                                 torch_tensor_from_array, torch_model_load,       &
                                 torch_model_forward, torch_delete
implicit none
private


!---------------------------------------------------------------------
! public interfaces
public  net_forward, nn_cf_net_init !nn_cf_net_finalize


!---------------------------------------------------------------------
! local/private data

! Neural Net Parameters
integer :: n_in
    !! Combined length of input features
integer :: n_features_out
    !! Number of output features (variables)
integer, allocatable, dimension(:) :: feature_out_sizes
    !! Vector storing the length of each of the input features
integer :: n_lev
    !! number of atmospheric layers (model levels) output is supplied for
! Dimension of each hidden layer
integer :: n_h1
integer :: n_h2
integer :: n_h3
integer :: n_h4
integer :: n_out


! Scale factors for inputs
real(4), allocatable, dimension(:)       :: xscale_mean
real(4), allocatable, dimension(:)       :: xscale_stnd


! Scale factors for outputs
real(4), allocatable, dimension(:)       :: yscale_mean
real(4), allocatable, dimension(:)       :: yscale_stnd


!---------------------------------------------------------------------
! Functions and Subroutines

contains

    !-----------------------------------------------------------------
    ! Public Subroutines
        

    subroutine net_forward(features, model_ftorch, logits)
        !! Run forward method of the Neural Net using FTorch.
        real(4), dimension(:), intent(inout) :: features
            !! Vector of input features
        real(4), dimension(:), intent(out)  :: logits
            !! Output vector
        type(torch_model), intent(in) :: model_ftorch        
        ! Local variables
        integer :: out_pos, feature_size, f
        integer :: num_features, num_outputs
        
        ! FTorch tensors
        type(torch_tensor), dimension(1) :: in_tensor, out_tensor
        integer :: in_layout(2), out_layout(2)
        real(r4), allocatable, target :: in_data(:,:)
        real(r4), allocatable, target :: out_data(:,:)
        
        ! Get dimensions
        num_features = size(features)
        num_outputs = size(logits)
        
        ! Allocate tensors with batch size of 1
        allocate(in_data(1, num_features))
        allocate(out_data(1, num_outputs))
        
        ! Normalize input features
        features = (features - xscale_mean) / xscale_stnd
        
        ! Copy normalized features to input tensor
        in_data(1, :) = features
        
        ! Set up tensor layouts (row-major: [batch, features])
        in_layout = [1, 2]
        out_layout = [1, 2]
        
        ! Create FTorch tensors from arrays
        call torch_tensor_from_array(in_tensor(1), in_data, in_layout, torch_kCPU)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        
        ! Forward pass through the model
        call torch_model_forward(model_ftorch, in_tensor, out_tensor)
        
        ! Copy output from tensor
        logits = out_data(1, :)
        
        ! Apply scaling and denormalization of each output feature in logits
        out_pos = 0
        do f = 1, n_features_out
        feature_size = feature_out_sizes(f)
        logits(out_pos+1:out_pos+feature_size) = &
            (logits(out_pos+1:out_pos+feature_size) * yscale_stnd(f)) + yscale_mean(f)
        out_pos = out_pos + feature_size
        end do
        
        
        
        ! Deallocate arrays
        deallocate(in_data)
        deallocate(out_data)

        ! Clean up tensors
        call torch_delete(in_tensor)
        call torch_delete(out_tensor)
    end subroutine net_forward

    subroutine nn_cf_net_init(nn_filename, metadata_filename, n_inputs, n_outputs, nrf, model_ftorch, iulog, errstring)
        !! Initialise the neural net using FTorch
        integer, intent(out)             :: n_inputs, n_outputs
        integer, intent(in)              :: nrf          ! number of atmospheric layers in each input
        type(torch_model), intent(out)   :: model_ftorch  
        integer, intent(in)              :: iulog
        character(128), intent(out)      :: errstring    ! output status (non-blank for error return)
        character(len=136), intent(in)   :: nn_filename   ! PyTorch model filename (.pt file)
        character(len=136), intent(in)   :: metadata_filename ! NetCDF file with dimensions and scaling parameters

        ! Local variables for reading NetCDF file
        integer :: ncid, varid, dimid
        integer :: temp_scalar(1)

        ! Initialize error string
        errstring = ''

        ! Open the metadata file
        call check(nf90_open(trim(metadata_filename), NF90_NOWRITE, ncid))

        ! Read dimension values (stored as scalar variables)
        call check(nf90_inq_varid(ncid, "n_inputs", varid))
        call check(nf90_get_var(ncid, varid, temp_scalar))
        n_inputs = temp_scalar(1)

         call check(nf90_inq_varid(ncid, "n_outputs", varid))
        call check(nf90_get_var(ncid, varid, temp_scalar))
        n_outputs = temp_scalar(1)
        
        call check(nf90_inq_varid(ncid, "n_h1", varid))
        call check(nf90_get_var(ncid, varid, temp_scalar))
        n_h1 = temp_scalar(1)

        call check(nf90_inq_varid(ncid, "n_h2", varid))
        call check(nf90_get_var(ncid, varid, temp_scalar))
        n_h2 = temp_scalar(1)
    
        call check(nf90_inq_varid(ncid, "n_h3", varid))
        call check(nf90_get_var(ncid, varid, temp_scalar))
        n_h3 = temp_scalar(1)
    
        call check(nf90_inq_varid(ncid, "n_h4", varid))
        call check(nf90_get_var(ncid, varid, temp_scalar))
        n_h4 = temp_scalar(1)
    
        call check(nf90_inq_varid(ncid, "n_features_out", varid))
        call check(nf90_get_var(ncid, varid, temp_scalar))
        n_features_out = temp_scalar(1)

        write(iulog, *) 'Read network dimensions from metadata file:'
        write(iulog, *) '  n_inputs = ', n_inputs
        write(iulog, *) '  n_outputs = ', n_outputs
        write(iulog, *) '  n_h1 = ', n_h1
        write(iulog, *) '  n_h2 = ', n_h2
        write(iulog, *) '  n_h3 = ', n_h3
        write(iulog, *) '  n_h4 = ', n_h4
        write(iulog, *) '  n_features_out = ', n_features_out

        ! Allocate scaling arrays
        if (allocated(xscale_mean)) deallocate(xscale_mean)
        if (allocated(xscale_stnd)) deallocate(xscale_stnd)
        if (allocated(yscale_mean)) deallocate(yscale_mean)
        if (allocated(yscale_stnd)) deallocate(yscale_stnd)
    
        allocate(xscale_mean(n_inputs))
        allocate(xscale_stnd(n_inputs))
        allocate(yscale_mean(n_features_out))
        allocate(yscale_stnd(n_features_out))

        ! Read input scaling parameters
        call check(nf90_inq_varid(ncid, "fscale_mean", varid))
        call check(nf90_get_var(ncid, varid, xscale_mean))
        call check(nf90_inq_varid(ncid, "fscale_stnd", varid))
        call check(nf90_get_var(ncid, varid, xscale_stnd))
    
        ! Read output scaling parameters
        call check(nf90_inq_varid(ncid, "oscale_mean", varid))
        call check(nf90_get_var(ncid, varid, yscale_mean))
        call check(nf90_inq_varid(ncid, "oscale_stnd", varid))
        call check(nf90_get_var(ncid, varid, yscale_stnd))

        call check(nf90_close(ncid))
    
        write(iulog, *) 'Loaded scaling parameters from: ', trim(metadata_filename)

        ! Load the model using FTorch interface
        call torch_model_load(model_ftorch, trim(nn_filename), torch_kCPU)

        write(iulog, *) 'Loaded PyTorch model from: ', trim(nn_filename)

        ! Set sizes of output features based on nrf
        n_lev = nrf

        ! Allocate and set sizes of the feature groups
        if (allocated(feature_out_sizes)) deallocate(feature_out_sizes)
        allocate(feature_out_sizes(n_features_out))
        feature_out_sizes(:)   = n_lev
        feature_out_sizes(2:3) = n_lev-1

    end subroutine nn_cf_net_init


    !subroutine nn_cf_net_finalize()
        !! Clean up NN space by deallocating arrays and destroying the FTorch model.
        
        ! Deallocate scaling parameter arrays
    !   if (allocated(xscale_mean)) deallocate(xscale_mean)
    !    if (allocated(xscale_stnd)) deallocate(xscale_stnd)
    !    if (allocated(yscale_mean)) deallocate(yscale_mean)
    !    if (allocated(yscale_stnd)) deallocate(yscale_stnd)
        
        ! Deallocate feature size array
    !    if (allocated(feature_out_sizes)) deallocate(feature_out_sizes)
        
        ! Delete the FTorch model to free GPU/CPU memory
    !    call torch_delete(model_ftorch)
        
        ! Reset initialization flag
    !   do_init = .true.
    
    !end subroutine nn_cf_net_finalize


    subroutine check(err_status)
        !! Check error status after netcdf call and print message for
        !! error codes.

        integer, intent(in) :: err_status
            !! error status from nf90 function

        if(err_status /= nf90_noerr) then
             write(*, *) trim(nf90_strerror(err_status))
        end if

    end subroutine check


end module nn_cf_net_mod
