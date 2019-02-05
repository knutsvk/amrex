#include "MyTest.H"
#include "makeEB.H"

#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();
    resizeArrays();
    makeEB();
    initData();
}

void
MyTest::solve ()
{
    LPInfo info;
    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    MLNodeLaplacian laplacian(geom, grids, dmap, info, GetVecOfConstPtrs(ebfactory));

    // This is a 3d problem with Dirichlet BC
    laplacian.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)},
                          {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)});

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        laplacian.setLevelBC(ilev, &solution[ilev]);
    }

    MLMG mlmg(laplacian);
    mlmg.setMaxIter(max_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
}


void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("ref_ratio", ref_ratio);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
}

void 
MyTest::resizeArrays ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);
    dmap.resize(nlevels);
    ebfactory.resize(nlevels);

    solution.resize(nlevels);
    sigma.resize(nlevels);
    rhs.resize(nlevels);
    exact_solution.resize(nlevels);
}

void 
MyTest::makeEB ()
{
    const int nlevels = geom.size();
    EBSupport m_eb_support_level = EBSupport::full;

	bool inside = true;
    Real radius = 0.5;
    Array<Real, 3> center = {0.5, 0.5, 0.5};
    EB2::SphereIF sphere(radius, center, inside);

    auto gshop = EB2::makeShop(sphere);
	EB2::Build(gshop, geom.back(), 0, 30);
    const EB2::IndexSpace& ebis = EB2::IndexSpace::top();

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        const EB2::Level& ebis_lev = ebis.getLevel(geom[ilev]);
        eb_level...?
        ebfactory[lev].reset(new EBFArrayBoxFactory(
    }

}


void
MyTest::initData ()
{
    const int nlevels = geom.size();

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(max_grid_size);
        domain.grow(-n_cell/4);   // fine level cover the middle of the coarse domain
        domain.refine(ref_ratio); 
        dmap[ilev].define(grids[ilev]);
    }

    initEB(geom);

    const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        const EB2::Level& ebis_lev = eb_is.getLevel(geom[ilev]);
        eb_level = & ebis_lev;
        ebfactory[lev].reset(new EBFArrayBoxFactory(*eb_level, 
                                                    geom[ilev], 
                                                    grids[ilev], 
                                                    dmap[ilev], 
                                                    {2, 2, 2}, 
                                                    EBSupport::full));
    }

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        solution      [ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *ebfactory[ilev]);
        sigma         [ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *ebfactory[ilev]);
        rhs           [ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *ebfactory[ilev]);
        exact_solution[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *ebfactory[ilev]);
    }

    // Initialize rhs and exact solution
    initProb();
}

