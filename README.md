# TIM: TURBO Infrastructure for MOM6

TIM (TURBO Infrastructure for MOM6) is a work-in-progress effort to build an
AMReX-based software infrastructure for the Modular Ocean Model version 6 (MOM6).

The overarching objective is to enable GPU acceleration, performance portability,
and ultimately ultra-high-resolution ocean simulations on current and emerging heterogeneous architectures.

TIM is being developed incrementally, with the initial phase focusing on extracting and refactoring
the subset of FMS functionality required by MOM6. Over time, this functionality will be re-implemented
in modern C++ and with AMReX to provide all necessary infrastructure features including but not limited to:

- Domain decomposition
- Parallelism via MPI and GPU offloading
- Device-portable kernels
- Memory management and tiling strategies for GPUs
- I/O
- Diagnostics
- etc.
