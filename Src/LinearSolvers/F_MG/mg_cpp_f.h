#ifndef _MG_CPP_F_H_
#define _MG_CPP_F_H_

#if defined(BL_FORT_USE_UPPERCASE)
#define mgt_alloc                 MGT_ALLOC
#define mgt_set_level             MGT_SET_LEVEL
#define mgt_finalize              MGT_FINALIZE
#define mgt_set_rh_2d             MGT_SET_RH_2D
#define mgt_get_rh_2d             MGT_GET_RH_2D
#define mgt_set_uu_2d             MGT_SET_UU_2D
#define mgt_get_uu_2d             MGT_GET_UU_2D
#define mgt_dealloc               MGT_DEALLOC
#define mgt_solve                 MGT_SOLVE
#define mgt_set_nu1               MGT_SET_NU1
#define mgt_set_nu2               MGT_SET_NU2
#define mgt_set_eps               MGT_SET_EPS
#define mgt_set_bottom_solver     MGT_SET_BOTTOM_SOLVER
#define mgt_set_bottom_max_iter   MGT_SET_BOTTOM_MAX_ITER
#define mgt_set_bottom_solver_eps MGT_SET_BOTTOM_SOLVER_EPS
#define mgt_set_max_nlevel        MGT_SET_MAX_NLEVEL
#define mgt_set_min_width         MGT_SET_MIN_WIDTH
#elif defined(BL_FORT_USE_UNDERSCORE)
#define mgt_alloc                 mgt_alloc_
#define mgt_set_level             mgt_set_level_
#define mgt_finalize              mgt_finalize_
#define mgt_set_rh_2d             mgt_set_rh_2d_
#define mgt_get_rh_2d             mgt_get_rh_2d_
#define mgt_set_uu_2d             mgt_set_uu_2d_
#define mgt_get_uu_2d             mgt_get_uu_2d_
#define mgt_dealloc               mgt_dealloc_
#define mgt_solve                 mgt_solve_
#define mgt_set_nu1               mgt_set_nu1_
#define mgt_set_nu2               mgt_set_nu2_
#define mgt_set_eps               mgt_set_eps_
#define mgt_set_bottom_solver     mgt_set_bottom_solver_
#define mgt_set_bottom_max_iter   mgt_set_bottom_max_iter_
#define mgt_set_bottom_solver_eps mgt_set_bottom_solver_eps_
#define mgt_set_max_nlevel        mgt_set_max_nlevel_
#define mgt_set_min_width         mgt_set_min_width_
#elif defined(BL_FORT_USE_DBL_UNDERSCORE)
#define mgt_alloc                 mgt_alloc__
#define mgt_set_level             mgt_set_level__
#define mgt_finalize              mgt_finalize__
#define mgt_set_rh_2d             mgt_set_rh_2d__
#define mgt_get_rh_2d             mgt_get_rh_2d__
#define mgt_set_uu_2d             mgt_set_uu_2d__
#define mgt_get_uu_2d             mgt_get_uu_2d__
#define mgt_dealloc               mgt_dealloc__
#define mgt_solve                 mgt_solve__
#define mgt_set_nu1               mgt_set_nu1__
#define mgt_set_nu2               mgt_set_nu2__
#define mgt_set_eps               mgt_set_eps__
#define mgt_set_bottom_solver     mgt_set_bottom_solver__
#define mgt_set_bottom_max_iter   mgt_set_bottom_max_iter__
#define mgt_set_bottom_solver_eps mgt_set_bottom_solver_eps__
#define mgt_set_max_nlevel        mgt_set_max_nlevel__
#define mgt_set_min_width         mgt_set_min_width__
#endif

#ifdef __cplusplus
extern "C" 
{
#endif

  void mgt_alloc(int* mgt, const int* dm, const int* nlevel);

  void mgt_set_level(const int* mgt, const int* lev, const int* nb, 
		     const int* dm, const int* lo, const int* hi, 
		     const int* pd_lo, const int* pd_hi, const int* bc, const int* pm);

  void mgt_finalize(const int* mgt);

  void mgt_set_rh_2d(const int* mgt, const int* lev, const int* n, const double* rh, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_rh_2d(const int* mgt, const int* lev, const int* n, double* rh, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_set_uu_2d(const int* mgt, const int* lev, const int* n, const double* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_get_uu_2d(const int* mgt, const int* lev, const int* n, double* uu, 
		     const int* plo, const int* phi, 
		     const int* lo, const int* hi);
  
  void mgt_dealloc(const int* mgt);
  
  void mgt_solve(const int* mgt);
  
  void mgt_set_nu1(const int* mgt, const int* nu1);
  
  void mgt_set_nu2(const int* mgt, const int* nu2);
  
  void mgt_set_eps(const int* mgt, const double* eps);
  
  void mgt_set_bottom_solver(const int* mgt, const int* bottom_solver);
  
  void mgt_set_bottom_max_iter(const int* mgt, const int* bottom_max_iter);
  
  void mgt_set_bottom_solver_eps(const int* mgt, const double* bottom_solver_eps);
  
  void mgt_set_max_nlevel(const int* mgt, const int* max_nlevel);
  
  void mgt_set_min_width(const int* mgt, const int* min_width);
  
#ifdef __cplusplus
}
#endif

#endif
