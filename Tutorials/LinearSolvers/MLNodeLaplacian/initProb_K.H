#ifndef INIT_PROB_K_H_
#define INIT_PROB_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_INLINE
void actual_init_poisson (int i, int j, int k,
                          amrex::Array4<amrex::Real> const& rhs,
                          amrex::Array4<amrex::Real> const& exact,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                          amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    constexpr amrex::Real  pi = 3.1415926535897932;
    constexpr amrex::Real tpi = 2.0 * pi;
    constexpr amrex::Real fpi = 4.0 * pi;
    amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5) - 0.5;
    amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5) - 0.5;
    amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5) - 0.5;
    amrex::Real r2 = x*x + y*y + z*z;

    exact(i,j,k) = std::sin(tpi * r2);
    rhs(i,j,k) = fpi * (fpi * r2 * std::sin(tpi * r2) - 3. * std::cos(tpi * r2));
}

#endif
