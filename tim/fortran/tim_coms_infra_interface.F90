module tim_coms_infra_interface

use iso_fortran_env, only : int32, int64
use iso_c_binding,   only : c_int64_t, c_double, c_size_t, c_ptr, c_null_ptr, c_loc

implicit none
private

public :: tim_chksum

interface tim_chksum_c
  function tim_chksum_c(field_ptr, field_size, mask_ptr) bind(c, name="tim_chksum_c")
    import c_ptr, c_int64_t, c_size_t
    integer(c_int64_t)                        :: tim_chksum_c
    type(c_ptr),            value, intent(in) :: field_ptr
    integer(c_size_t),      value, intent(in) :: field_size
    type(c_ptr),            value, intent(in) :: mask_ptr
  end function tim_chksum_c
end interface tim_chksum_c

interface tim_chksum
  module procedure tim_chksum_real_0d
  module procedure tim_chksum_real_1d
  module procedure tim_chksum_real_2d
  module procedure tim_chksum_real_3d
  module procedure tim_chksum_real_4d
end interface tim_chksum

contains

function tim_chksum_real_0d(field, pelist, mask_val) result(chksum)
  real,              target, intent(in) :: field               !< Input scalar
  integer, optional,         intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,    optional, target, intent(in) :: mask_val            !< FMS mask value
  type(c_ptr)                           :: field_loc, mask_loc !< c pointers to field and mask
  integer(kind=int64)                   :: chksum              !< checksum of array

  field_loc = c_loc(field)
  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  chksum = tim_chksum_c(field_loc, int(1, kind=c_size_t), mask_loc)
end function tim_chksum_real_0d

function tim_chksum_real_1d(field, pelist, mask_val) result(chksum)
  real, dimension(:), target, intent(in) :: field               !< Input array
  integer,  optional,         intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,     optional, target, intent(in) :: mask_val            !< FMS mask value
  type(c_ptr)                            :: field_loc, mask_loc !< c pointers to field and mask
  integer(kind=int64)                    :: chksum              !< checksum of array

  field_loc = c_loc(field(1))
  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  chksum = tim_chksum_c(field_loc, int(size(field), kind=c_size_t), mask_loc)
end function tim_chksum_real_1d

function tim_chksum_real_2d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:), target, intent(in) :: field               !< Unrotated input field
  integer,    optional,         intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,       optional, target, intent(in) :: mask_val            !< FMS mask value
  type(c_ptr)                              :: field_loc, mask_loc !< c pointers to field and mask
  integer(kind=int64)                      :: chksum              !< checksum of array

  field_loc = c_loc(field(1,1))
  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  chksum = tim_chksum_c(field_loc, int(size(field), kind=c_size_t), mask_loc)
end function tim_chksum_real_2d

function tim_chksum_real_3d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:), target, intent(in) :: field               !< Unrotated input field
  integer,      optional,         intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,         optional, target, intent(in) :: mask_val            !< FMS mask value
  type(c_ptr)                                :: field_loc, mask_loc !< c pointers to field and mask
  integer(kind=int64)                        :: chksum              !< checksum of array

  field_loc = c_loc(field(1,1,1))
  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  chksum = tim_chksum_c(field_loc, int(size(field), kind=c_size_t), mask_loc)
end function tim_chksum_real_3d

function tim_chksum_real_4d(field, pelist, mask_val) result(chksum)
  real, dimension(:,:,:,:), target, intent(in) :: field               !< Unrotated input field
  integer,        optional,         intent(in) :: pelist(:)           !< PE list of ranks to checksum
  real,           optional, target, intent(in) :: mask_val            !< FMS mask value
  type(c_ptr)                                  :: field_loc, mask_loc !< c pointers to field and mask
  integer(kind=int64)                          :: chksum              !< checksum of array

  field_loc = c_loc(field(1,1,1,1))
  if(present(mask_val)) then
    mask_loc = c_loc(mask_val)
  else
    mask_loc = c_null_ptr
  end if

  chksum = tim_chksum_c(field_loc, int(size(field), kind=c_size_t), mask_loc)
end function tim_chksum_real_4d

end module tim_coms_infra_interface

