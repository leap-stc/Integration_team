module FTorch_cesm_interface
#ifdef USE_FTORCH
  use ftorch,         only: torch_kCPU, torch_tensor, torch_model, torch_tensor_from_array
  use ftorch,         only: torch_model_load, torch_model_forward, torch_delete
#else
  use shr_abort_mod, only : shr_abort_abort
  use shr_kind_mod, only : r4 => shr_kind_r4
  type :: torch_model
  end type torch_model
  type :: torch_tensor
  end type torch_tensor
  integer :: torch_kcpu
  character(len=*), parameter :: message="ERROR: Using FTorch Interface without USE_FTORCH=TRUE"
contains
  subroutine torch_model_load(model, fname, i)
    type(torch_model) :: model
    character(len=*) :: fname
    integer :: i
    call shr_abort_abort(message)
  end subroutine torch_model_load
  subroutine torch_tensor_from_array(tensor, data, layout, i)
    type(torch_tensor) :: tensor
    real(r4) :: data(:,:)
    integer :: layout(:)
    integer :: i
    call shr_abort_abort(message)
  end subroutine torch_tensor_from_array
  subroutine torch_model_forward(model, tensor_in, tensor_out)
    type(torch_model) :: model
    type(torch_tensor) :: tensor_in(:), tensor_out(:)
    call shr_abort_abort(message)
  end subroutine torch_model_forward
  subroutine torch_delete(tensor)
    type(torch_tensor) :: tensor
    call shr_abort_abort(message)
  end subroutine torch_delete
#endif
end module FTorch_cesm_interface
