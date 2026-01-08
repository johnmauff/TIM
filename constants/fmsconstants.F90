!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!> @defgroup fmsconstants FMSConstants
!> @ingroup libfms
!> @brief Defines useful constants for Earth. Constants are defined as real
!!
!>    FMSconstants have been declared as REAL(kind=sizeof(rvar)), PARAMETER.
!!
!!    The value of a constant defined and used from here cannot be changed
!!    in a users program. New constants can be defined in terms of values
!!    from the FMSconstants module and their includes using a parameter
!!    statement.<br><br>
!!
!!    The currently support contant systems are:
!!       GFDL constants (gfdl_constants.fh)
!!       GEOS constants (geos_constants.fh)
!!       GFS  constants (gfs_constants.fh)
!!       <br><br>
!!
!!    The name given to a particular constant may be changed.<br><br>
!!
!!    Constants can only be used on the right side on an assignment statement
!!    (their value can not be reassigned).
!!
!!    Example:
!!
!! @verbatim
!!    use FMSConstants, only:  TFREEZE, grav_new => GRAV
!!    real, parameter :: grav_inv = 1.0 / grav_new
!!    tempc(:,:,:) = tempk(:,:,:) - TFREEZE
!!    geopotential(:,:) = height(:,:) * grav_new
!! @endverbatim
!> @file
!> @brief File for @ref FMSconstants_mod

!> @addtogroup FMSconstants_mod
!> @{
module FMSconstants

  use platform_mod, only: r4_kind, r8_kind

#define CESM_CONSTANTS
#if defined(CESM_CONSTANTS)
use shr_kind_mod,   only : R8 => shr_kind_r8
use shr_const_mod,  only : &
     SHR_CONST_PSTD,                      & ! standard pressure ~ pascals
     PI              => SHR_CONST_PI,     &
     PI_8            => SHR_CONST_PI,     &
     SECONDS_PER_DAY => SHR_CONST_CDAY,   &
     OMEGA           => SHR_CONST_OMEGA,  &
     RADIUS          => SHR_CONST_REARTH, &
     GRAV            => SHR_CONST_G,      &
     STEFAN          => SHR_CONST_STEBOL, & ! Stefan-Boltzmann constant ~ W/m^2/K^4
     AVOGNO          => SHR_CONST_AVOGAD, & ! Avogadro's number ~ molecules/kmole
     WTMAIR          => SHR_CONST_MWDAIR, & ! molecular weight dry air ~ kg/kmole
     WTMH2O          => SHR_CONST_MWWV,   & ! molecular weight water vapor
     RDGAS           => SHR_CONST_RDAIR,  & ! Dry air gas constant     ~ J/K/kg
     RVGAS           => SHR_CONST_RWV,    & ! Water vapor gas constant ~ J/K/kg
     VONKARM         => SHR_CONST_KARMAN, & ! Von Karman constant
     TFREEZE         => SHR_CONST_TKFRZ,  & ! freezing T of fresh water          ~ K
     RHOAIR          => SHR_CONST_RHODAIR,& ! density of dry air at STP  ~ kg/m^3
     DENS_H2O        => SHR_CONST_RHOFW,  & ! density of fresh water     ~ kg/m^3
     RHO0            => SHR_CONST_RHOSW,  & ! density of sea water       ~ kg/m^3
     CP_AIR          => SHR_CONST_CPDAIR, & ! specific heat of dry air   ~ J/kg/K
     CP_VAPOR        => SHR_CONST_CPWV,   & ! specific heat of water vap ~ J/kg/K
     CP_OCEAN        => SHR_CONST_CPSW,   &
     HLF             => SHR_CONST_LATICE, & ! latent heat of fusion      ~ J/kg
     HLV             => SHR_CONST_LATVAP    ! latent heat of evaporation ~ J/kg
#endif

  !--- default scoping
  implicit none

  !--- needed with implicit none
  real :: dum  !< dummy real variable

#define RKIND sizeof(dum)

!--- set a default for the FMSConstants
#if defined(CESM_CONSTANTS)
#include <cesm_constants.fh>
#else

#if !defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#define GFDL_CONSTANTS
#endif

!--- perform error checking and include the correct system of constants
#if defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#warning "Using GFDL constants"
#include <gfdl_constants.fh>
#elif !defined(GFDL_CONSTANTS) && defined(GFS_CONSTANTS) && !defined(GEOS_CONSTANTS)
#warning "Using GFS constants"
#include <gfs_constants.fh>
#elif !defined(GFDL_CONSTANTS) && !defined(GFS_CONSTANTS) && defined(GEOS_CONSTANTS)
#warning "Using GEOS constants"
#include <geos_constants.fh>
#else
#error FATAL FMSConstants error -  multiple constants macros are defined for FMS
#endif

#endif ! ifdef CESM_CONSTANTS

  !--- public interfaces
  public :: FMSConstants_init

  contains

    !> @brief FMSconstants init routine
    subroutine FMSconstants_init
      use mpp_mod, only: stdlog
      integer :: logunit
      logunit = stdlog()

      write (logunit,'(/,80("="),/(a))') trim(constants_version)

    end subroutine FMSconstants_init

end module FMSconstants
!> @}
! close documentation grouping
