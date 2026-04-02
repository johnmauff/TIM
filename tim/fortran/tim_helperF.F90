module tim_helperF

!#include <timH.h>
  logical, parameter :: use_AMREX = .TRUE.
  !integer, parameter :: TIMH_runAMREX   = _TIMH_RUNAMREX, &
  !                      TIMH_capture    = _TIMH_CAPTURE, &
  !                      TIMH_runFORTRAN = _TIMH_RUNFORTRAN
  !integer, parameter :: TIMH_CAPTURE_INPUT = _TIMH_CAPTURE_INPUT,  &
  !                      TIMH_CAPTURE_OUTPUT = _TIMH_CAPTURE_OUTPUT, &
  !                      TIMH_RUN = _TIMH_RUN
  integer, parameter :: TIMH_runAMREX   = 0, &
                        TIMH_capture    = 1, &
                        TIMH_runFORTRAN = 2
  integer, parameter :: TIMH_CAPTURE_INPUT = 101,  &
                        TIMH_CAPTURE_OUTPUT = 102, &
                        TIMH_RUN = 103

   public :: TIMH_runAMREX, TIMH_capture, TIMH_runFORTRAN
   public :: TIMH_CAPTURE_INPUT, TIMH_CAPTURE_OUTPUT, TIMH_RUN
   public :: get_mode_env

contains

  function get_mode_env(name, default) result(mode)
    character(len=*), intent(in) :: name
    integer, intent(in), optional :: default
    integer :: mode

    character(len=:), allocatable :: str
    integer :: length, status

    ! Get length first
    call get_environment_variable(name, length=length, status=status)

    if (status /= 0 .or. length == 0) then
      if (present(default)) then
        mode = default
      else
        mode = MODE_FORTRAN   ! sensible fallback
      end if
      return
    end if

    allocate(character(len=length) :: str)
    call get_environment_variable(name, str)

    ! Normalize (optional but recommended)
    call to_upper(str)

    select case (trim(str))
    case ("AMREX")
      mode = TIMH_runAMREX
    case ("FORTRAN")
      mode = TIMH_runFORTRAN
    case ("CAPTURE")
      mode = TIMH_capture
    case default
      print *, "ERROR: Invalid value for ", name, " = ", trim(str)
      print *, "Allowed values: AMREX, FORTRAN, CAPTURE"
      stop 1
    end select

  end function get_mode_env


  subroutine to_upper(s)
    character(len=*), intent(inout) :: s
    integer :: i
    do i = 1, len(s)
      select case (s(i:i))
      case ('a':'z')
        s(i:i) = achar(iachar(s(i:i)) - 32)
      end select
    end do
  end subroutine to_upper

end module tim_helperF
