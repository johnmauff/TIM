#ifndef A4_HELPER_H_
#define A4_HELPER_H_
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_Array4.H>

namespace A4 {

struct A4Box
{
    amrex::Box bx;
    amrex::Real* data = nullptr;
    amrex::Array4<amrex::Real> arr;

    int nx, ny, nz, ncomp;
};

A4Box make_a4(int nx, int ny, int nz, int ncomp=1)
{
    using namespace amrex;

    A4Box a4;

    a4.nx = nx;
    a4.ny = ny;
    a4.nz = nz;
    a4.ncomp = ncomp;

    a4.bx = Box(IntVect(0,0,0), IntVect(nx-1, ny-1, nz-1));

    // std::size_t n = (std::size_t)nx * ny * nz * ncomp;

    const Long npts = a4.bx.numPts();

    // Allocate in AMReX arena (GPU-aware)
    a4.data = (Real*) The_Arena()->alloc(npts * sizeof(Real));

    a4.arr = Array4<Real>(a4.data, lbound(a4.bx), ubound(a4.bx),ncomp);

    return a4;
}

void free_a4(A4Box& a4)
{
    using namespace amrex;

    if (a4.data)
    {
        The_Arena()->free(a4.data);
        a4.data = nullptr;
    }
}

void copy_f_to_a4(const double* f, A4Box& a4)
{
    using namespace amrex;

    int nx = a4.nx;
    int ny = a4.ny;

    auto arr = a4.arr;
    Box bx = a4.bx;

    auto lo = lbound(bx);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        int ii = i - lo.x;
	int jj = j - lo.y;
	int kk = k - lo.z;
        int idx = ii + nx*(jj + ny*kk);  // Fortran layout
        arr(i,j,k) = f[idx];
    });
}

void copy_a4_to_f(const A4Box& a4, double* f)
{
    using namespace amrex;

    int nx = a4.nx;
    int ny = a4.ny;

    auto arr = a4.arr;
    Box bx = a4.bx;

    auto lo = lbound(bx);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        int ii = i - lo.x;
	int jj = j - lo.y;
	int kk = k - lo.z;
        int idx = ii + nx*(jj + ny*kk);  // Fortran layout
        f[idx] = arr(i,j,k);
    });
}
}
#endif
