#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

class MyTest
{
public:

    MyTest ();

    void solve ();
    void writePlotfile () const;


private:
    void readParameters ();
    void resizeArrays ();
    void initData ();
    void initProb ();

    int max_level = 1;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;
    
    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    
    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;
    amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>> ebfactory;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> sigma;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> exact_solution;
};

#endif
