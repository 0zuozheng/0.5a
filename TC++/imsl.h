/*      Copyright:  2005 by IMSL, Inc.  All Rights Reserved. */

#ifndef IMSL_H
#define IMSL_H

#include "imslerr.h"

#ifndef TIME_H__
#  define TIME_H__
#  include <time.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

#ifdef WAVE_RENAME
#  include "rename.h"
#else
#  if defined(COMPUTER_VAX) || defined(COMPUTER_VAXG) || defined(COMPUTER_ALFAC_IEEE)
#    define imsl_d_random_normal_multivariate imsl_d_random_normal_multivaria
#    define imsl_f_random_normal_multivariate imsl_f_random_normal_multivaria
#  endif
#endif

/* IMSL C/Math/Library version 7.0.0 */
#define IMSL_CMATH_VERSION      700

#define IMSL_PROTO(P,Q)         P Q

#if defined(_WIN32) || defined(__WIN32__)
#  define IMSL_DECL     __cdecl
#  ifdef __BORLANDC__
#    ifdef _RTLDLL
#      ifndef IMSL_EF
#        ifdef _EXPFUNC
#          define IMSL_EF       _EXPFUNC
#        else
#          define IMSL_EF       __import
#        endif
#        ifdef _EXPDATA
#          define IMSL_ED       _EXPDATA
#        else
#          define IMSL_ED       __import
#        endif
#      endif
#    else
#      define IMSL_EF
#      define IMSL_ED
#    endif
#    define IMSL_CI
#  else
     /* Visual C++ */
#    if defined(_DLL) && !defined(IMSL_STATIC)
#      ifndef IMSL_CI
#        ifdef _CRTIMP
#          define IMSL_CI       _CRTIMP
#        else
#          define IMSL_CI       __declspec(dllimport)
#        endif
#      endif
#    else
#      define IMSL_CI
#    endif
#    define IMSL_EF
#    define IMSL_ED
#  endif
#else
#  define IMSL_DECL
#  define IMSL_CI
#  define IMSL_EF
#  define IMSL_ED
#endif

#if defined(__cplusplus)        /* C++ compiler is being used */
extern "C" {
#endif

#define IMSL_F                        float

/* Complex structure. */
/* Do not redeclare f_complex and d_complex if C/Stat has already done so. */
#ifndef IMSLS_H
typedef struct {
    float                             re;
    float                             im;
} f_complex;
typedef struct {
    double                            re;
    double                            im;
} d_complex;
#endif

/* Peicewise polynomial structures. */
typedef struct {
    int                               domain_dim;
    int                               target_dim;
    int                               *order;
    int                               *num_coef;
    int                               *num_breakpoints;
    float                             **breakpoints;
    float                             **coef;
} Imsl_f_ppoly;

typedef struct {
    int                               domain_dim;
    int                               target_dim;
    int                               *order;
    int                               *num_coef;
    int                               *num_breakpoints;
    double                            **breakpoints;
    double                            **coef;
} Imsl_d_ppoly;

/* B-spline structures. */
typedef struct {
    int                               domain_dim;
    int                               target_dim;
    int                               *order;
    int                               *num_coef;
    int                               *num_knots;
    float                             **knots;
    float                             **coef;
} Imsl_f_spline;

typedef struct {
    int                               domain_dim;
    int                               target_dim;
    int                               *order;
    int                               *num_coef;
    int                               *num_knots;
    double                            **knots;
    double                            **coef;
} Imsl_d_spline;

/* Constraint structures used for imsl_*_spline_lsq_constrained. */
typedef struct {
    float                             xval;
    int                               der;
    int                               type;
    float                             bl;
    float                             bu;
} f_constraint_struct;

typedef struct {
    double                            xval;
    int                               der;
    int                               type;
    double                            bl;
    double                            bu;
} d_constraint_struct;
 
/* Radial Basis structures. */
typedef IMSL_F IMSL_PROTO((IMSL_DECL *f_radial_fcn),(IMSL_F));
typedef double IMSL_PROTO((IMSL_DECL *d_radial_fcn),(double));
typedef struct {
    int                               dimension;
    int                               num_centers;
    int                               additional_terms;
    float                             *centers;
    float                             *coefficients;
    f_radial_fcn                      radial_function;
    float                             delta;
    void                              *fcn_2_data;
    float                             (*fcn_2) (float, void *);
} Imsl_f_radial_basis_fit;

typedef struct {
    int                               dimension;
    int                               num_centers;
    int                               additional_terms;
    double                            *centers;
    double                            *coefficients;
    d_radial_fcn                      radial_function;
    double                            delta;
    void                              *fcn_2_data;
    double                            (*fcn_2) (double, void *);
} Imsl_d_radial_basis_fit;

/* Sparse matrix data structures */
typedef struct {
    int                               row;
    int                               col;
    float                             val;
} Imsl_f_sparse_elem;

typedef struct {
    int                               row;
    int                               col;
    double                            val;
} Imsl_d_sparse_elem;

typedef struct {
    int                               row;
    int                               col;
    f_complex                         val;
} Imsl_c_sparse_elem;

typedef struct {
    int                               row;
    int                               col;
    d_complex                         val;
} Imsl_z_sparse_elem;

struct Imsl_f_sparse_list_element_struct {
    float                                    val;
    int                                      i;
    int                                      j;
    struct Imsl_f_sparse_list_element_struct *next_row;
    struct Imsl_f_sparse_list_element_struct *next_col;
};
typedef struct Imsl_f_sparse_list_element_struct Imsl_f_sparse_list_element;

struct Imsl_d_sparse_list_element_struct {
    double                                   val;
    int                                      i;
    int                                      j;
    struct Imsl_d_sparse_list_element_struct *next_row;
    struct Imsl_d_sparse_list_element_struct *next_col;
};
typedef struct Imsl_d_sparse_list_element_struct Imsl_d_sparse_list_element;
 
struct Imsl_c_sparse_list_element_struct {
    f_complex                                val;
    int                                      i;   
    int                                      j;
    struct Imsl_c_sparse_list_element_struct *next_row;
    struct Imsl_c_sparse_list_element_struct *prev_row;
    struct Imsl_c_sparse_list_element_struct *next_col;
    struct Imsl_c_sparse_list_element_struct *prev_col;
};
typedef struct Imsl_c_sparse_list_element_struct Imsl_c_sparse_list_element;
 
struct Imsl_z_sparse_list_element_struct {
    d_complex                                val;
    int                                      i;
    int                                      j;
    struct Imsl_z_sparse_list_element_struct *next_row;
    struct Imsl_z_sparse_list_element_struct *prev_row;
    struct Imsl_z_sparse_list_element_struct *next_col;
    struct Imsl_z_sparse_list_element_struct *prev_col;
};
typedef struct Imsl_z_sparse_list_element_struct Imsl_z_sparse_list_element;
 
typedef struct {
    int                               markowitz;
    Imsl_f_sparse_list_element        *ptr;
} Imsl_f_header_element;

typedef struct {
    int                               markowitz;
    Imsl_d_sparse_list_element        *ptr;
} Imsl_d_header_element;

typedef struct {
    int                               markowitz;
    Imsl_c_sparse_list_element        *ptr;
} Imsl_c_header_element;

typedef struct {
    int                               markowitz;
    Imsl_z_sparse_list_element        *ptr;
} Imsl_z_header_element;

struct Imsl_f_memory_list_struct {
    Imsl_f_sparse_list_element        *node;
    struct Imsl_f_memory_list_struct  *next;
};
typedef struct Imsl_f_memory_list_struct Imsl_f_memory_list;

struct Imsl_d_memory_list_struct {
    Imsl_d_sparse_list_element        *node;
    struct Imsl_d_memory_list_struct  *next;
};
typedef struct Imsl_d_memory_list_struct Imsl_d_memory_list;

struct Imsl_c_memory_list_struct {
    Imsl_c_sparse_list_element        *node;
    struct Imsl_c_memory_list_struct  *next;
};
typedef struct Imsl_c_memory_list_struct Imsl_c_memory_list;

struct Imsl_z_memory_list_struct {
    Imsl_z_sparse_list_element        *node;
    struct Imsl_z_memory_list_struct  *next;
};
typedef struct Imsl_z_memory_list_struct Imsl_z_memory_list;

typedef struct {
    int                               n;
    int                               nz;
    int                               *row_pivots;
    int                               *col_pivots;
    Imsl_f_header_element             *row_ptr;
    Imsl_f_header_element             *col_ptr;
} Imsl_f_sparse_lu_factor;

typedef struct {
    int                               n;
    int                               nz;
    int                               *row_pivots;
    int                               *col_pivots;
    Imsl_d_header_element             *row_ptr;
    Imsl_d_header_element             *col_ptr;
} Imsl_d_sparse_lu_factor;

typedef struct {
    int                               n;
    int                               nz;
    int                               *row_pivots;
    int                               *col_pivots;
    Imsl_c_header_element             *row_ptr;
    Imsl_c_header_element             *col_ptr;
} Imsl_c_sparse_lu_factor;

typedef struct {
    int                               n;
    int                               nz;
    int                               *row_pivots;
    int                               *col_pivots;
    Imsl_z_header_element             *row_ptr;
    Imsl_z_header_element             *col_ptr;
} Imsl_z_sparse_lu_factor;

typedef struct {
    int                               **nzsub;
    int                               **xnzsub;
    int                               maxsub;
    int                               **xlnz;
    int                               maxlnz;
    int                               **perm;
    int                               **invp;
    int                               multifrontal_space;
} Imsl_symbolic_factor;

typedef struct {
    int                               **nzsub;
    int                               **xnzsub;
    int                               **xlnz;
    float                             **alnz;
    int                               **perm;
    float                             **diag;
} Imsl_f_numeric_factor;

typedef struct { 
    int                               **nzsub; 
    int                               **xnzsub; 
    int                               **xlnz; 
    double                            **alnz; 
    int                               **perm; 
    double                            **diag; 
} Imsl_d_numeric_factor;

typedef struct {
    int                               **nzsub;
    int                               **xnzsub;
    int                               **xlnz;
    f_complex                         **alnz;
    int                               **perm;
    f_complex                         **diag;
} Imsl_c_numeric_factor;

typedef struct {
    int                               **nzsub;
    int                               **xnzsub;
    int                               **xlnz;
    d_complex                         **alnz;
    int                               **perm;
    d_complex                         **diag;
} Imsl_z_numeric_factor;

/*      Monte Carlo data structure */
typedef struct {
    int                               nSkip;
    int                               *y;
    int                               dim;
    int                               maxDigits;
    int                               base;
    int                               *digitsN;
    int                               *c;
    int                               *power;
    double                            scale;
} Imsl_faure;

/* Description of an LP or QP problem */
typedef struct {
    char                              *filename;
    char                              name[9];
    int                               nrows;
    int                               ncolumns;
    int                               nonzeros;
    int                               nhessian;
    int                               ninteger;
    int                               nbinary;
    float                             *objective;
    Imsl_f_sparse_elem                *constraint;
    Imsl_f_sparse_elem                *hessian;
    float                             *lower_range;
    float                             *upper_range;
    float                             *lower_bound;
    float                             *upper_bound;
    int                               *variable_type;
    char                              name_objective[9];
    char                              name_rhs[9];
    char                              name_ranges[9];
    char                              name_bounds[9];
    char                              **name_row;
    char                              **name_column;
    float                             positive_infinity;
    float                             negative_infinity;
} Imsl_f_mps;

typedef struct {
    char                              *filename;
    char                              name[9];
    int                               nrows;
    int                               ncolumns;
    int                               nonzeros;
    int                               nhessian;
    int                               ninteger;
    int                               nbinary;
    double                            *objective;
    Imsl_d_sparse_elem                *constraint;
    Imsl_d_sparse_elem                *hessian;
    double                            *lower_range;
    double                            *upper_range;
    double                            *lower_bound;
    double                            *upper_bound;
    int                               *variable_type;
    char                              name_objective[9];
    char                              name_rhs[9];
    char                              name_ranges[9];
    char                              name_bounds[9];
    char                              **name_row;
    char                              **name_column;
    double                            positive_infinity;
    double                            negative_infinity;
} Imsl_d_mps;

/* Various enumerated types. */
                                /* enums for error severity */
typedef enum {
    IMSL_NOTE                         = 1,
    IMSL_ALERT                        = 2,
    IMSL_WARNING                      = 3,
    IMSL_FATAL                        = 4,
    IMSL_TERMINAL                     = 5,
    IMSL_WARNING_IMMEDIATE            = 6,
    IMSL_FATAL_IMMEDIATE              = 7
} Imsl_error;

                                /* enums for ODEs */
typedef enum {
    IMSL_ODE_INITIALIZE               = 1,
    IMSL_ODE_CHANGE                   = 2,
    IMSL_ODE_RESET                    = 3
} Imsl_ode;

                                /* enums for DEAs */
typedef enum {
    IMSL_DEA_INITIALIZE               = 1,
    IMSL_DEA_RESET                    = 3
} Imsl_dea;

                                /* enums for PDEs */
typedef enum {
    IMSL_PDE_INITIALIZE               = 1,
    IMSL_PDE_CHANGE                   = 2,
    IMSL_PDE_RESET                    = 3
} Imsl_pde;

                                /* enums for fast_poisson_2d */
typedef enum {
    IMSL_DIRICHLET_BC                 = 1,
    IMSL_NEUMANN_BC                   = 2,
    IMSL_PERIODIC_BC                  = 3
} Imsl_bc_type;
typedef enum {
    IMSL_RIGHT_SIDE                   = 0,
    IMSL_BOTTOM_SIDE                  = 1,
    IMSL_LEFT_SIDE                    = 2,
    IMSL_TOP_SIDE                     = 3
} Imsl_pde_side;
                                /* enums for imsl_page */
typedef enum {
    IMSL_SET_PAGE_WIDTH               =-1,
    IMSL_GET_PAGE_WIDTH               = 1,
    IMSL_SET_PAGE_LENGTH              =-2,
    IMSL_GET_PAGE_LENGTH              = 2
} Imsl_page_options;
                                /* enums for imsl_write_options */
typedef enum {
    IMSL_SET_DEFAULTS                 = 0,
    IMSL_SET_CENTERING                =-1,
    IMSL_GET_CENTERING                = 1,
    IMSL_SET_ROW_WRAP                 =-2,
    IMSL_GET_ROW_WRAP                 = 2,
    IMSL_SET_PAGING                   =-3,
    IMSL_GET_PAGING                   = 3,
    IMSL_SET_NAN_CHAR                 =-4,
    IMSL_GET_NAN_CHAR                 = 4,
    IMSL_SET_TITLE_PAGE               =-5,
    IMSL_GET_TITLE_PAGE               = 5,
    IMSL_SET_FORMAT                   =-6,
    IMSL_GET_FORMAT                   = 6
} Imsl_write_options;
                                /* enums for quadrature routines */
typedef enum {
    IMSL_ALG                          = 1,
    IMSL_ALG_LEFT_LOG                 = 2,
    IMSL_ALG_RIGHT_LOG                = 3,
    IMSL_ALG_LOG                      = 4,
    IMSL_INF_BOUND                    = 5,
    IMSL_BOUND_INF                    = 6,
    IMSL_INF_INF                      = 7,
    IMSL_COS                          = 8,
    IMSL_SIN                          = 9
} Imsl_quad;

                                /* enums for inverse laplace transform */
typedef enum {
    IMSL_NORMAL_TERMINATION,
    IMSL_TOO_LARGE,
    IMSL_TOO_SMALL,
    IMSL_TOO_LARGE_BEFORE_EXPANSION,
    IMSL_TOO_SMALL_BEFORE_EXPANSION
} Imsl_laplace_flow;

                                /* enums for general sparse routines */
typedef enum {
    IMSL_ROW_MARKOWITZ                = 0,
    IMSL_COLUMN_MARKOWITZ             = 1,
    IMSL_SYMMETRIC_MARKOWITZ          = 2
} Imsl_pivot;

                                /* enums for financial functions */
typedef enum {
    IMSL_AT_END_OF_PERIOD             = 0,
    IMSL_AT_BEGINNING_OF_PERIOD       = 1,

    IMSL_BASISPART_30E360             = 0,
    IMSL_BASISPART_365                = 1,
    IMSL_BASISPART_ACTUAL             = 2,
    IMSL_BASISPART_NASD               = 3,

    IMSL_DAY_CNT_BASIS_ACTUALACTUAL   = 0,
    IMSL_DAY_CNT_BASIS_NASD           = 1,
    IMSL_DAY_CNT_BASIS_ACTUAL360      = 2,
    IMSL_DAY_CNT_BASIS_ACTUAL365      = 3,
    IMSL_DAY_CNT_BASIS_30E360         = 4,

    IMSL_ANNUAL                       = 1,
    IMSL_SEMIANNUAL                   = 2,
    IMSL_QUARTERLY                    = 4
} Imsl_finance;



#if (defined(__alpha) && defined(__osf__)) || defined(__ia64) || defined(__64BIT__) || defined(_FLOAT0) || defined(__x86_64)
typedef void    IMSL_PROTO((IMSL_DECL *Imsl_error_print_proc),
                           (Imsl_error, int, char*, char*));
#elif defined(COMPUTER_LOPT64) || defined(COMPUTER_INTEL64) || defined(COMPUTER_MACOSX) || defined(COMPUTER_LOPT64_PGI) || defined(COMPUTER_MACX64)
typedef void    IMSL_PROTO((IMSL_DECL *Imsl_error_print_proc),
                           (Imsl_error, int, char*, char*));
#else
typedef void    IMSL_PROTO((IMSL_DECL *Imsl_error_print_proc),
                           (Imsl_error, long, char*, char*));
#endif
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_gen,
                        (int, float*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_gen,
                        (int, double*, double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_posdef,
                        (int, float*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_posdef,
                        (int, double*, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_lin_sol_posdef,
                        (int, f_complex*, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_lin_sol_posdef,
                        (int, d_complex*, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_nonnegdef,
                        (int, float*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_nonnegdef,
                        (int, double*, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_lin_sol_gen,
                        (int, f_complex*, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_lin_sol_gen,
                        (int, d_complex*, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_least_squares_gen,
                        (int, int, float*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_least_squares_gen,
                        (int, int, double*, double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_mat_mul_rect,
                        (char*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_mat_mul_rect,
                        (char*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_mat_mul_rect,
                        (char*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_mat_mul_rect,
                        (char*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_svd_gen,
                        (int, int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_svd_gen,
                        (int, int, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_lin_svd_gen,
                        (int, int, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_lin_svd_gen,
                        (int, int, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_lsq_lin_constraints,
                        (int, int, int, float*, float*, float*,
                         float*, float*, int*, float*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_lsq_lin_constraints,
                        (int, int, int, double*, double*, double*,
                         double*, double*, int*, double*, double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_gen_coordinate,
                        (int, int, Imsl_f_sparse_elem*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_gen_coordinate,
                        (int, int, Imsl_d_sparse_elem*, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_lin_sol_gen_coordinate,
                        (int, int, Imsl_c_sparse_elem*, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_lin_sol_gen_coordinate,
                        (int, int, Imsl_z_sparse_elem*, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_posdef_coordinate,
                        (int, int, Imsl_f_sparse_elem*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_posdef_coordinate,
                        (int, int, Imsl_d_sparse_elem*, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_lin_sol_posdef_coordinate,
                        (int, int, Imsl_c_sparse_elem*, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_lin_sol_posdef_coordinate,
                        (int, int, Imsl_z_sparse_elem*, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_gen_band,
                        (int, float*, int, int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_gen_band,
                        (int, double*, int, int, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_lin_sol_gen_band,
                        (int, f_complex*, int, int, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_lin_sol_gen_band,
                        (int, d_complex*, int, int, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_posdef_band,
                        (int, float*, int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_posdef_band,
                        (int, double*, int, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_lin_sol_posdef_band,
                        (int, f_complex*, int, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_lin_sol_posdef_band,
                        (int, d_complex*, int, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_def_cg,
                        (int, void(IMSL_DECL *)(float*, float*),
                         float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_def_cg,
                        (int, void(IMSL_DECL *)(double*, double*),
                         double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_sol_gen_min_residual,
                        (int, void(IMSL_DECL *)(float*, float*),
                         float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_sol_gen_min_residual,
                        (int, void(IMSL_DECL *)(double*, double*),
                         double*, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_free_symbolic_factor,
                        (Imsl_symbolic_factor*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_free_numeric_factor,
                        (Imsl_f_numeric_factor*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_free_numeric_factor,
                        (Imsl_d_numeric_factor*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_free_numeric_factor,
                        (Imsl_c_numeric_factor*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_free_numeric_factor,
                        (Imsl_z_numeric_factor*));
IMSL_CI Imsl_f_sparse_elem    * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_generate_test_coordinate,
                        (int, int, int*, ...));
IMSL_CI Imsl_d_sparse_elem    * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_generate_test_coordinate,
                        (int, int, int*, ...));
IMSL_CI Imsl_c_sparse_elem    * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_generate_test_coordinate,
                        (int, int, int*, ...));
IMSL_CI Imsl_z_sparse_elem    * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_generate_test_coordinate,
                        (int, int, int*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_generate_test_band,
                        (int, int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_generate_test_band,
                        (int, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_generate_test_band,
                        (int, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_generate_test_band,
                        (int, int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_mat_add_band,
                        (int, int, int, IMSL_F, float*, int, int,
                         IMSL_F, float*, int*, int*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_mat_add_band,
                        (int, int, int, double, double*, int, int,
                         double, double*, int*, int*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_mat_add_band,
                        (int, int, int, f_complex, f_complex*, int, int,
                         f_complex, f_complex*, int*, int*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_mat_add_band,
                        (int, int, int, d_complex, d_complex*, int, int,
                         d_complex, d_complex*, int*, int*, ...));
IMSL_CI Imsl_f_sparse_elem      * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_mat_add_coordinate,
                        (int, int, IMSL_F, Imsl_f_sparse_elem*, int,
                         IMSL_F, Imsl_f_sparse_elem*, int*, ...));
IMSL_CI Imsl_d_sparse_elem      * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_mat_add_coordinate,
                        (int, int, double, Imsl_d_sparse_elem*, int,
                         double, Imsl_d_sparse_elem*, int*, ...));
IMSL_CI Imsl_c_sparse_elem      * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_mat_add_coordinate,
                        (int, int, f_complex, Imsl_c_sparse_elem*, int,
                         f_complex, Imsl_c_sparse_elem*, int*, ...));
IMSL_CI Imsl_z_sparse_elem      * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_mat_add_coordinate,
                        (int, int, d_complex, Imsl_z_sparse_elem*, int,
                         d_complex, Imsl_z_sparse_elem*, int*, ...));
IMSL_CI void          * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_mat_mul_rect_coordinate,
                        (char*, ...));
IMSL_CI void          * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_mat_mul_rect_coordinate,
                        (char*, ...));
IMSL_CI void          * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_mat_mul_rect_coordinate,
                        (char*, ...));
IMSL_CI void          * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_mat_mul_rect_coordinate,
                        (char*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_mat_mul_rect_band,
                        (char*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_mat_mul_rect_band,
                        (char*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_mat_mul_rect_band,
                        (char*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_mat_mul_rect_band,
                        (char*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_eig_sym,
                        (int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_eig_sym,
                        (int, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_eig_gen,
                        (int, float*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_eig_gen,
                        (int, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_eig_gen,
                        (int, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_eig_gen,
                        (int, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_eig_herm,
                        (int, f_complex*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_eig_herm,
                        (int, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_geneig_sym_posdef,
                        (int, float*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_geneig_sym_posdef,
                        (int, double*, double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_eig_symgen,
                        (int, float*, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_eig_symgen,
                        (int, double*, double*, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_geneig,
                        (int, float*, float*, f_complex*, float*,
                         ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_geneig,
                        (int, double*, double*, d_complex*, double*,
                         ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_geneig,
                        (int, f_complex*, f_complex*, f_complex*,
                         f_complex*, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_geneig,
                        (int, d_complex*, d_complex*, d_complex*,
                         d_complex*, ...));
IMSL_CI Imsl_f_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_ppoly_create,
                        (int, int, int*, int*, ...));
IMSL_CI Imsl_d_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_ppoly_create,
                        (int, int, int*, int*, ...));
IMSL_CI Imsl_f_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cub_spline_interp,
                        (int, float[], float[]));
IMSL_CI Imsl_d_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cub_spline_interp,
                        (int, double[], double[]));
IMSL_CI Imsl_f_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cub_spline_tcb,
                        (int, float[], float[], ...));
IMSL_CI Imsl_d_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cub_spline_tcb,
                        (int, double[], double[], ...));
IMSL_CI Imsl_f_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cub_spline_interp_e_cnd,
                        (int, float[], float[], ...));
IMSL_CI Imsl_d_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cub_spline_interp_e_cnd,
                        (int, double[], double[], ...));
IMSL_CI Imsl_f_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cub_spline_interp_shape,
                        (int, float[], float[], ...));
IMSL_CI Imsl_d_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cub_spline_interp_shape,
                        (int, double[], double[], ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cub_spline_value,
                        (IMSL_F, Imsl_f_ppoly*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cub_spline_value,
                        (double, Imsl_d_ppoly*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cub_spline_integral,
                        (IMSL_F, IMSL_F, Imsl_f_ppoly*));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cub_spline_integral,
                        (double, double, Imsl_d_ppoly*));
IMSL_CI Imsl_f_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_create,
                        (int, int, int*, int*, ...));
IMSL_CI Imsl_d_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_create,
                        (int, int, int*, int*, ...));
IMSL_CI Imsl_f_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_interp,
                        (int, float[], float[], ...));
IMSL_CI Imsl_d_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_interp,
                        (int, double[], double[], ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_knots,
                        (int, float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_knots,
                        (int, double[], ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_value,
                        (IMSL_F, Imsl_f_spline*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_value,
                        (double, Imsl_d_spline*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_integral,
                        (IMSL_F, IMSL_F, Imsl_f_spline*));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_integral,
                        (double, double, Imsl_d_spline*));
IMSL_CI Imsl_f_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_2d_interp,
                        (int, float[], int, float[], float[],...));
IMSL_CI Imsl_d_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_2d_interp,
                        (int, double[], int, double[], double[],...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_2d_value,
                        (IMSL_F, IMSL_F, Imsl_f_spline*,...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_2d_value,
                        (double, double, Imsl_d_spline*,...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_2d_integral,
                        (IMSL_F, IMSL_F, IMSL_F, IMSL_F, Imsl_f_spline*));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_2d_integral,
                        (double, double, double, double, Imsl_d_spline*));
IMSL_CI Imsl_f_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_2d_least_squares,
                        (int, float[], int, float[], float[], int,
                         int, ...));
IMSL_CI Imsl_d_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_2d_least_squares,
                        (int, double[], int, double[], double[], int,
                         int, ...));
IMSL_CI Imsl_f_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_least_squares,
                        (int, float[], float[], int,...));
IMSL_CI Imsl_d_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_least_squares,
                        (int, double[], double[], int,...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_poly_regression,
                        (int, float[], float[], int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_poly_regression,
                        (int, double[], double[], int, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_print,
                        (Imsl_f_spline*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_print,
                        (Imsl_d_spline*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_ppoly_print,
                        (Imsl_f_ppoly*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_ppoly_print,
                        (Imsl_d_ppoly*));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_user_fcn_least_squares,
                        (IMSL_F(IMSL_DECL *)(int, IMSL_F), int, int,
                         float[], float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_user_fcn_least_squares,
                        (double(IMSL_DECL *)(int, double), int, int,
                         double[], double[], ...));
IMSL_CI Imsl_f_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cub_spline_smooth,
                        (int, float[], float[], ...));
IMSL_CI Imsl_d_ppoly  * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cub_spline_smooth,
                        (int, double[], double[], ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_scattered_2d_interp,
                        (int, float[], float[], int, int, float[],
                         float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_scattered_2d_interp,
                        (int, double*, double[], int, int, double[], 
                         double[], ...));
IMSL_CI Imsl_f_radial_basis_fit   * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_radial_scattered_fit, 
                        (int, int, float*, float*, int, ...));
IMSL_CI Imsl_d_radial_basis_fit   * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_radial_scattered_fit, 
                        (int, int, double*, double*, int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_radial_evaluate, 
                        (int, float*, Imsl_f_radial_basis_fit*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_radial_evaluate, 
                        (int, double*, Imsl_d_radial_basis_fit*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_radial_function,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_radial_function,
                        (double));
IMSL_CI Imsl_f_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_spline_lsq_constrained,
                        (int, float[], float[], int, int,
                         f_constraint_struct[], ...));
IMSL_CI Imsl_d_spline * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_spline_lsq_constrained,
                        (int, double[], double[], int, int,
                         d_constraint_struct[], ...));
IMSL_CI float          * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_smooth_1d_data, (int, 
                        float[], float[],...));
IMSL_CI double         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_smooth_1d_data, (int, 
                        double[], double[],...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn,
                        (double(IMSL_DECL *)(double), double, double,
                         ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_sing_pts,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         int, float*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_sing_pts,
                        (double(IMSL_DECL *)(double), double, double,
                         int, double*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_sing,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_sing,
                        (double(IMSL_DECL *)(double), double, double,
                         ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_inf,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, Imsl_quad,
                         ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_inf,
                        (double(IMSL_DECL *)(double), double, Imsl_quad,
                         ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_trig,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         Imsl_quad, IMSL_F, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_trig,
                        (double(IMSL_DECL *)(double), double, double,
                         Imsl_quad, double, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_alg_log,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         Imsl_quad, IMSL_F, IMSL_F, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_alg_log,
                        (double(IMSL_DECL *)(double), double, double,
                         Imsl_quad, double, double, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_fourier,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, Imsl_quad,
                         IMSL_F, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_fourier,
                        (double(IMSL_DECL *)(double), double, Imsl_quad,
                         double, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_smooth,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_smooth,
                        (double(IMSL_DECL *)(double), double, double,
                         ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_cauchy,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         IMSL_F, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_cauchy,
                        (double(IMSL_DECL *)(double), double, double,
                         double, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_hyper_rect,
                        (IMSL_F(IMSL_DECL *)(int, float*), int,
                         float*, float*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_hyper_rect,
                        (double(IMSL_DECL *)(int, double*), int,
                         double*, double*, ...));
/* multivariate normal cdf version of hyper_rect */
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_hyper_rect2,
                        (IMSL_F(IMSL_DECL *)(int, float*), int,
                         float*, float*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_hyper_rect2,
                        (double(IMSL_DECL *)(int, double*), int,
                         double*, double*, ...));
/************************************************ */
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_2d,
                        (IMSL_F(IMSL_DECL *)(IMSL_F, IMSL_F), IMSL_F,
                         IMSL_F, IMSL_F(IMSL_DECL *)(IMSL_F),
                         IMSL_F(IMSL_DECL *)(IMSL_F), ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_2d,
                        (double(IMSL_DECL *)(double, double), double,
                         double, double(IMSL_DECL *)(double),
                         double(IMSL_DECL *)(double), ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_gauss_quad_rule,
                        (int, float[], float[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_gauss_quad_rule,
                        (int, double[], double[], ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fcn_derivative,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fcn_derivative,
                        (double(IMSL_DECL *)(double), double, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_ode_runge_kutta,
                        (int, float*, IMSL_F, float*, char*,
                         void(IMSL_DECL *)(int, IMSL_F, float*, float*)));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_ode_runge_kutta,
                        (int, double*, double, double*, char*,
                         void(IMSL_DECL *)(int, double, double*, double*)));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_ode_runge_kutta_mgr,
                        (int, char**, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_ode_runge_kutta_mgr,
                        (int, char**, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_ode_adams_gear,
                        (int, float*, IMSL_F, float*, char*,
                         void(IMSL_DECL *)(int, IMSL_F, float*, float*)));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_ode_adams_gear,
                        (int, double*, double, double*, char*,
                         void(IMSL_DECL *)(int, double, double*, double*)));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_ode_adams_gear_mgr,
                        (int, char**, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_ode_adams_gear_mgr,
                        (int, char**, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_pde_method_of_lines_mgr,
                        (int, void**, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_pde_method_of_lines_mgr,
                        (int, void**, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_pde_method_of_lines,
                        (int, float*, IMSL_F, int, float[], float[], void*,
                         void(IMSL_DECL *)(int, IMSL_F, IMSL_F, float*,
                         float*, float*, float*),
                         void(IMSL_DECL *)(int, IMSL_F, IMSL_F, float*,
                         float*, float*)));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_pde_method_of_lines,
                        (int, double*, double, int, double[], double[], void*,
                         void(IMSL_DECL *)(int, double, double, double*,
                         double*, double*, double*),
                         void(IMSL_DECL *)(int, double, double, double*,
                          double*, double*)));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fast_poisson_2d,
                        (IMSL_F(IMSL_DECL *)(IMSL_F, IMSL_F),
                         IMSL_F(IMSL_DECL *)(Imsl_pde_side, IMSL_F, IMSL_F),
                         IMSL_F, int, int, IMSL_F, IMSL_F, IMSL_F,
                         IMSL_F, Imsl_bc_type*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fast_poisson_2d,
                        (double(IMSL_DECL *)(double, double),
                         double(IMSL_DECL *)(Imsl_pde_side, double, double),
                         double, int, int, double, double, double,
                         double, Imsl_bc_type*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fft_real,
                        (int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fft_real,
                        (int, double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fft_real_init,
                        (int));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fft_real_init,
                        (int));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_fft_complex_init,
                        (int));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_fft_complex_init,
                        (int));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_fft_complex,
                        (int, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_fft_complex,
                        (int, d_complex*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_fft_2d_complex,
                        (int, int, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_fft_2d_complex,
                        (int, int, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fft_sine,
                        (int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fft_sine,
                        (int, double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fft_sine_init,
                        (int));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fft_sine_init,
                        (int));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fft_cosine,
                        (int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fft_cosine,
                        (int, double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fft_cosine_init,
                        (int));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fft_cosine_init,
                        (int));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_convolution,
                        (int, float*, int, float*, int*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_convolution,
                        (int, double*, int, double*, int*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_convolution,
                        (int, f_complex*, int, f_complex*, int*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_convolution,
                        (int, d_complex*, int, d_complex*, int*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_inverse_laplace,
                        (f_complex(IMSL_DECL *)(f_complex),
                         IMSL_F, int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_inverse_laplace,
                        (d_complex(IMSL_DECL *)(d_complex),
                         double, int, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_zeros_poly,
                        (int, float*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_zeros_poly,
                        (int, double*, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_zeros_poly,
                        (int, f_complex*, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_zeros_poly,
                        (int, d_complex*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_zeros_function,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_zeros_function,
                        (double(IMSL_DECL *)(double), ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_zeros_fcn,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_zeros_fcn,
                        (double(IMSL_DECL *)(double), ...));

IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_zeros_fcn,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_zeros_fcn,
                        (double(IMSL_DECL *)(double), ...));

IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_zeros_sys_eqn,
                        (void(IMSL_DECL *)(int, float[], float[]),
                         int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_zeros_sys_eqn,
                        (void(IMSL_DECL *)(int, double[], double[]),
                         int, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_min_uncon,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_min_uncon,
                        (double(IMSL_DECL *)(double), double, double,
                         ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_min_uncon_deriv,
                        (IMSL_F(IMSL_DECL *)(IMSL_F),
                         IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_min_uncon_deriv,
                        (double(IMSL_DECL *)(double),
                         double(IMSL_DECL *)(double), double, double,
                         ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_min_uncon_multivar,
                        (IMSL_F(IMSL_DECL *)(int, float[]), int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_min_uncon_multivar,
                        (double(IMSL_DECL *)(int, double[]), int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_nonlin_least_squares,
                        (void(IMSL_DECL *)(int, int, float[], float[]),
                         int, int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_nonlin_least_squares,
                        (void(IMSL_DECL *)(int, int, double[], double[]),
                         int, int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_min_con_nonlin,
                        (void(IMSL_DECL *)(int, int, int, float[],
                                           int[], float*, float[]),
                         int, int, int, int, float[], float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_min_con_nonlin,
                        (void(IMSL_DECL *)(int, int, int, double[],
                                           int[], double*, double[]), 
                         int, int, int, int, double[], double[], ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_lin_prog,
                        (int, int, float*, float[], float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_lin_prog,
                        (int, int, double*, double[], double[], ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_quadratic_prog,
                        (int, int, int, float*, float*, float*,
                         float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_quadratic_prog,
                        (int, int, int, double*, double*, double*,
                         double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bounded_least_squares,
                        (void(IMSL_DECL *)(int, int, float[], float[]),
                         int, int, int, float[], float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bounded_least_squares,
                        (void(IMSL_DECL *)(int, int, double[], double[]),
                         int, int, int, double[], double[], ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_min_con_gen_lin,
                        (void(IMSL_DECL *)(int, float*, float*),
                         int, int, int, float[], float[], float[],
                         float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_min_con_gen_lin,
                        (void(IMSL_DECL *)(int, double*, double*),
                         int, int, int, double[], double[], double[],
                         double[], ...));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_i_min,
                        (int, int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_min,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_min,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_vmin,
                        (int, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_vmin,
                        (int, ...));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_i_max,
                        (int, int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_max,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_max,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_vmax,
                        (int, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_vmax,
                        (int, ...));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_ii_power,
                        (int, int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_fi_power,
                        (IMSL_F, int));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_di_power,
                        (double, int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_ff_power,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_dd_power,
                        (double, double));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_ci_power,
                        (f_complex, int));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_zi_power,
                        (d_complex, int));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_cf_power,
                        (f_complex, IMSL_F));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_zd_power,
                        (d_complex, double));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_cc_power,
                        (f_complex, f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_zz_power,
                        (d_complex, d_complex));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_erf,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_erf,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_erfc,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_erfc,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_erf_inverse,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_erf_inverse,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_erfc_inverse,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_erfc_inverse,
                        (double));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_erfe,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_erfe,
                        (d_complex));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_erfce,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_erfce,
                        (double x));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_J0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_J0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_J1,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_J1,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_I0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_I0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_I1,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_I1,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_Y0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_Y0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_Y1,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_Y1,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_K0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_K0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_K1,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_K1,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_exp_I0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_exp_I0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_exp_I1,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_exp_I1,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_exp_K0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_exp_K0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bessel_exp_K1,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bessel_exp_K1,
                        (double));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Ix,
                        (IMSL_F, f_complex, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Ix,
                        (double, d_complex, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Kx,
                        (IMSL_F, f_complex, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Kx,
                        (double, d_complex, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Jx,
                        (IMSL_F, f_complex, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Jx,
                        (double, d_complex, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Yx,
                        (IMSL_F, f_complex, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Yx,
                        (double, d_complex, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Ix_adr,
                        (float*, f_complex*, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Ix_adr,
                        (double*, d_complex*, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Jx_adr,
                        (float*, f_complex*, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Jx_adr,
                        (double*, d_complex*, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Kx_adr,
                        (float*, f_complex*, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Kx_adr,
                        (double*, d_complex*, int, ...));
IMSL_CI f_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_bessel_Yx_adr,
                        (float*, f_complex*, int, ...));
IMSL_CI d_complex     * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_bessel_Yx_adr,
                        (double*, d_complex*, int, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_gamma,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_gamma,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_log_gamma,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_log_gamma,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_gamma_incomplete,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_gamma_incomplete,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_t_inverse_cdf,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_t_inverse_cdf,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_normal_inverse_cdf,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_normal_inverse_cdf,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_binomial_cdf,
                        (int, int, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_binomial_cdf,
                        (int, int, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_normal_cdf,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_normal_cdf,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_chi_squared_cdf,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_chi_squared_cdf,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_chi_squared_inverse_cdf,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_chi_squared_inverse_cdf,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_F_cdf,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_F_cdf,
                        (double, double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_gamma_cdf,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_gamma_cdf,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_t_cdf,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_t_cdf,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_F_inverse_cdf,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_F_inverse_cdf,
                        (double, double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_hypergeometric_cdf,
                        (int, int, int, int));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_hypergeometric_cdf,
                        (int, int, int, int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_poisson_cdf,
                        (int, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_poisson_cdf,
                        (int, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_beta,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_beta,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_beta_incomplete,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_beta_incomplete,
                        (double, double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_log_beta,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_log_beta,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_elliptic_integral_K,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_elliptic_integral_K,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_elliptic_integral_E,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_elliptic_integral_E,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_elliptic_integral_RC,
                        (IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_elliptic_integral_RC,
                        (double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_elliptic_integral_RD,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_elliptic_integral_RD,
                        (double, double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_elliptic_integral_RF,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_elliptic_integral_RF,
                        (double, double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_elliptic_integral_RJ,
                        (IMSL_F, IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_elliptic_integral_RJ,
                        (double, double, double, double)); 
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fresnel_integral_C,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fresnel_integral_C,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_fresnel_integral_S,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_fresnel_integral_S,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_airy_Ai,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_airy_Ai,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_airy_Bi,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_airy_Bi,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_airy_Ai_derivative,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_airy_Ai_derivative,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_airy_Bi_derivative,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_airy_Bi_derivative,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_ber0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_ber0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_bei0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_bei0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_ker0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_ker0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_kei0,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_kei0,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_ber0_derivative,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_ber0_derivative,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_bei0_derivative,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_bei0_derivative,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_ker0_derivative,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_ker0_derivative,
                        (double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_kelvin_kei0_derivative,
                        (IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_kelvin_kei0_derivative,
                        (double));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_i_write_matrix,
                        (char*, int, int, int[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_write_matrix,
                        (char*, int, int, float[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_write_matrix,
                        (char*, int, int, f_complex[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_write_matrix,
                        (char*, int, int, double[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_write_matrix,
                        (char*, int, int, d_complex[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_page,
                        (int, int*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_write_options,
                        (int, int*));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_table_oneway,
                        (int, float[], int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_table_oneway,
                        (int, double[], int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_regression,
                        (int, int, float*, float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_regression,
                        (int, int, double*, double[], ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_ranks,
                        (int, float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_ranks,
                        (int, double[], ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_simple_statistics,
                        (int, int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_simple_statistics,
                        (int, int, double*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_chi_squared_test,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), int, int,
                         float*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_chi_squared_test,
                        (double(IMSL_DECL *)(double), int, int,
                         double*, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_covariances,
                        (int, int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_covariances,
                        (int, int, double*, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_random_seed_set,
                        (int));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_random_seed_get,
                        (void));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_random_option,
                        (int));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_random_uniform,
                        (int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_random_uniform,
                        (int, ...));
IMSL_CI int           * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_random_poisson,
                        (int, IMSL_F, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_random_normal,
                        (int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_random_normal,
                        (int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_random_exponential,
                        (int, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_random_exponential,
                        (int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_random_gamma,
                        (int, IMSL_F, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_random_gamma,
                        (int, double, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_random_beta,
                        (int, IMSL_F, IMSL_F, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_random_beta,
                        (int, double, double, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_random_normal_multivariate,
                        (int, int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_random_normal_multivariate,
                        (int, int, double*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_beta_cdf,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_beta_cdf,
                        (double, double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bivariate_normal_cdf,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bivariate_normal_cdf,
                        (double, double, double));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_beta_inverse_cdf,
                        (IMSL_F, IMSL_F, IMSL_F));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_beta_inverse_cdf,
                        (double, double, double));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_output_file,
                        (int, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_ctime,
                        (void));
#if (defined(__alpha) && defined(__osf__)) || defined(__ia64) || defined(__64BIT__) || defined(_FLOAT0) || defined(__x86_64)
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_error_code,
                        (void));
#elif defined(COMPUTER_SG64XS) || defined(COMPUTER_FUJITSU) || defined(COMPUTER_LINUX64) || defined(COMPUTER_SL64XS) || defined(COMPUTER_HP64XS) || defined(COMPUTER_HPSI64) || defined(COMPUTER_FUJSOL64) || defined(COMPUTER_LOPT64) || defined(COMPUTER_INTEL64) || defined(COMPUTER_MACOSX) || defined(COMPUTER_LOPT64_PGI) || defined(COMPUTER_MACX64)
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_error_code,
                        (void));
#else
IMSL_CI long            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_error_code,
                        (void));
#endif
IMSL_CI Imsl_error      IMSL_DECL IMSL_EF IMSL_PROTO(imsl_error_type,
                        (void));
IMSL_CI char          * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_error_message,
                        (void));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_error_options,
                        (int, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_omp_options,
                        (int, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_skip_signal_handler,
                        ());
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_free,
                        (void *));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_free_allocated_memory,
                        (void *));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_i_machine,
                        (int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_machine,
                        (int));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_machine,
                        (int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_constant,
                        (char*, char*));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_constant,
                        (char*, char*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_days_to_date,
                        (int, int*, int*, int*));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_date_to_days,
                        (int, int, int));
IMSL_CI char          * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_version,
                        (int));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_vector_norm,
                        (int, float*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_vector_norm,
                        (int, double*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_vector_norm,
                        (int, f_complex*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_vector_norm,
                        (int, d_complex*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_matrix_norm,
                        (int, int, float*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_matrix_norm,
                        (int, int, double*, ...));
IMSL_CI float           IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_matrix_norm_coordinate,
                        (int, int, int, Imsl_f_sparse_elem*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_matrix_norm_coordinate,
                        (int, int, int, Imsl_d_sparse_elem*, ...));
IMSL_CI float           IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_matrix_norm_band,
                        (int, float*, int, int, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_matrix_norm_band,
                        (int, double*, int, int, ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_sort,
                        (int, float*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_sort,
                        (int, double*, ...));
IMSL_CI int           * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_i_sort,
                        (int, int*, ...));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_neg,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_neg,
                        (d_complex ));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_add,
                        (f_complex, f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_add,
                        (d_complex, d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_sub,
                        (f_complex, f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_sub,
                        (d_complex, d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_mul,
                        (f_complex, f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_mul,
                        (d_complex, d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_div,
                        (f_complex, f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_div,
                        (d_complex, d_complex));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_eq,
                        (f_complex, f_complex));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_eq,
                        (d_complex, d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_cz_convert,
                        (d_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_zc_convert,
                        (f_complex));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_aimag,
                        (f_complex));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_aimag,
                        (d_complex));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_fc_convert,
                        (f_complex));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_dz_convert,
                        (d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_cf_convert,
                        (IMSL_F, IMSL_F));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_zd_convert,
                        (double, double));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_conjg,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_conjg,
                        (d_complex));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_arg,
                        (f_complex));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_arg,
                        (d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_sqrt,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_sqrt,
                        (d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_log,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_log,
                        (d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_exp,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_exp,
                        (d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_sin,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_sin,
                        (d_complex));
IMSL_CI f_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_cos,
                        (f_complex));
IMSL_CI d_complex       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_cos,
                        (d_complex));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_c_abs,
                        (f_complex));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_z_abs,
                        (d_complex));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_future_value,
                        (IMSL_F rate, int nper, IMSL_F pmt, IMSL_F pv, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_future_value,
                        (double rate, int nper, double pmt, double pv, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_present_value,
                        (IMSL_F rate, int nper, IMSL_F pmt, IMSL_F fv, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_present_value,
                        (double rate, int nper, double pmt, double fv, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_payment,
                        (IMSL_F rate, int nper, IMSL_F pv, IMSL_F fv, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_payment,
                        (double rate, int nper, double pv, double fv, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_effective_rate,
                        (IMSL_F nominalRate, int nper));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_effective_rate,
                        (double nominalRate, int nper));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_nominal_rate,
                        (IMSL_F effectiveRate, int nper));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_nominal_rate,
                        (double effectiveRate, int nper));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_interest_payment,
                        (IMSL_F rate, int period, int nper, IMSL_F pv, IMSL_F fv, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_interest_payment,
                        (double rate, int period, int nper, double pv, double fv, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_dollar_decimal,
                        (IMSL_F fractionalDollar, int fraction));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_dollar_decimal,
                        (double fractionalDollar, int fraction));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_dollar_fraction,
                        (IMSL_F decimalDollar, int fraction));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_dollar_fraction,
                        (double decimalDollar, int fraction));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_number_of_periods,
                        (IMSL_F rate, IMSL_F pmt, IMSL_F pv, IMSL_F fv, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_number_of_periods,
                        (double rate, double pmt, double pv, double fv, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_principal_payment,
                        (IMSL_F rate, int period, int nper, IMSL_F pv, IMSL_F fv, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_principal_payment,
                        (double rate, int period, int nper, double pv, double fv, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_net_present_value,
                        (IMSL_F rate, int count, IMSL_F value[]));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_net_present_value,
                        (double rate, int count, double value[]));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_depreciation_ddb,
                        (IMSL_F cost, IMSL_F salvage, int life, int period, IMSL_F factor));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_depreciation_ddb,
                        (double cost, double salvage, int life, int period, double factor));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_depreciation_db,
                        (IMSL_F cost, IMSL_F salvage, int life, int period, int month));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_depreciation_db,
                        (double cost, double salvage, int life, int period, int month));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cumulative_interest,
                        (IMSL_F rate, int nper, IMSL_F pv,      int start, int end, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cumulative_interest,
                        (double rate, int nper, double pv,      int start, int end, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_cumulative_principal,
                        (IMSL_F rate, int nper, IMSL_F pv,      int start, int end, int when));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_cumulative_principal,
                        (double rate, int nper, double pv,      int start, int end, int when));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_future_value_schedule,
                        (IMSL_F principal, int count, IMSL_F schedule[]));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_future_value_schedule,
                        (double principal, int count, double schedule[]));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_depreciation_sln,
                        (IMSL_F cost, IMSL_F salvage, int life));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_depreciation_sln,
                        (double cost, double salvage, int life));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_depreciation_syd,
                        (IMSL_F cost, IMSL_F salvage, int life, int per));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_depreciation_syd,
                        (double cost, double salvage, int life, int per));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_modified_internal_rate,
                        (int count, IMSL_F value[], IMSL_F financeRate, IMSL_F reinvestRate));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_modified_internal_rate,
                        (int count, double value[],     double financeRate, double reinvestRate));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_depreciation_vdb,
                        (IMSL_F cost, IMSL_F salvage, int life,
                        int start, int end, IMSL_F factor, int no_sl));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_depreciation_vdb,
                        (double cost, double salvage, int life,
                        int start, int end, double factor, int no_sl));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_present_value_schedule,
                        (IMSL_F rate, int count, IMSL_F value[], struct tm dates[]));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_present_value_schedule,
                        (double rate, int count, double value[], struct tm dates[]));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_internal_rate_of_return,
                        (int count, IMSL_F pmt[], ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_internal_rate_of_return,
                        (int count, double pmt[], ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_interest_rate_annuity,
                        (int nper, IMSL_F pmt, IMSL_F pv, IMSL_F fv, int when, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_interest_rate_annuity,
                        (int nper, double pmt, double pv, double fv, int when, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_internal_rate_schedule,
                        (int count, IMSL_F pmt[], struct tm dates[], ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_internal_rate_schedule,
                        (int count, double pmt[], struct tm dates[], ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_year_fraction,
                        (struct tm start, struct tm end, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_year_fraction,
                        (struct tm start, struct tm end, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_coupon_days,
                        (struct tm settlement, struct tm maturity, int frequency, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_coupon_days,
                        (struct tm settlement, struct tm maturity, int frequency, int basis));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_coupon_number,
                        (struct tm settlement, struct tm maturity, int frequency, int basis));
IMSL_CI struct tm       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_next_coupon_date,
                        (struct tm settlement, struct tm maturity, int frequency, int basis));
IMSL_CI struct tm       IMSL_DECL IMSL_EF IMSL_PROTO(imsl_previous_coupon_date,
                        (struct tm settlement, struct tm maturity, int frequency, int basis));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_days_to_next_coupon,
                        (struct tm settlement, struct tm maturity, int frequency, int basis));
IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_days_before_settlement,
                        (struct tm settlement, struct tm maturity, int frequency, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_price,
                        (struct tm settlement, struct tm maturity, IMSL_F rate, IMSL_F yield, 
                        IMSL_F redemption, int frequency, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_price,
                        (struct tm settlement, struct tm maturity, double rate, double yield, 
                         double redemption, int frequency, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_discount_rate,
                        (struct tm settlement, struct tm maturity, IMSL_F price, 
                        IMSL_F redemption, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_discount_rate,
                        (struct tm settlement, struct tm maturity, double price, 
                        double redemption, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_discount_price,
                        (struct tm settlement, struct tm maturity, IMSL_F discount, 
                        IMSL_F redemption, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_discount_price,
                        (struct tm settlement, struct tm maturity, double discount, 
                        double redemption, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_duration,
                        (struct tm settlement, struct tm maturity, IMSL_F coupon, 
                        IMSL_F yield, int frequency, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_duration,
                        (struct tm settlement, struct tm maturity, double coupon, 
                        double yield, int frequency, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_modified_duration,
                        (struct tm settlement, struct tm maturity, IMSL_F coupon, 
                        IMSL_F yield, int frequency, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_modified_duration,
                        (struct tm settlement, struct tm maturity, double coupon, 
                        double yield, int frequency, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_price_maturity,
                        (struct tm settlement, struct tm maturity, struct tm issue, 
                        IMSL_F rate, IMSL_F yield, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_price_maturity,
                        (struct tm settlement, struct tm maturity, struct tm issue, 
                        double rate, double yield, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_received_maturity,
                        (struct tm settlement, struct tm maturity, IMSL_F investment, 
                        IMSL_F discount, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_received_maturity,
                        (struct tm settlement, struct tm maturity, double investment, 
                        double discount, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bond_equivalent_yield,
                        (struct tm settlement, struct tm maturity, IMSL_F rate));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bond_equivalent_yield,
                        (struct tm settlement, struct tm maturity, double rate));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_treasury_bill_price,
                        (struct tm settlement, struct tm maturity, IMSL_F rate));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_treasury_bill_price,
                        (struct tm settlement, struct tm maturity, double rate));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_treasury_bill_yield,
                        (struct tm settlement, struct tm maturity, IMSL_F price));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_treasury_bill_yield,
                        (struct tm settlement, struct tm maturity, double price));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_interest_rate_security,
                        (struct tm settlement, struct tm maturity, IMSL_F investment, 
                        IMSL_F redemption, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_interest_rate_security,
                        (struct tm settlement, struct tm maturity, double investment, 
                        double redemption, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_accr_interest_periodic,
                        (struct tm issue, struct tm first_coupon, struct tm settlement,
                        IMSL_F rate, IMSL_F par, int frequency, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_accr_interest_periodic,
                        (struct tm issue, struct tm first_coupon, struct tm settlement,
                        double rate, double par, int frequency, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_accr_interest_maturity,
                        (struct tm issue, struct tm maturity, IMSL_F rate, IMSL_F par, 
                        int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_accr_interest_maturity,
                        (struct tm issue, struct tm maturity, double rate, double par, 
                        int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_depreciation_amordegrc,
                        (IMSL_F cost, struct tm issue, struct tm first_period, IMSL_F salvage, 
                        int period, IMSL_F rate, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_depreciation_amordegrc,
                        (double cost, struct tm issue, struct tm first_period, 
                        double salvage, int period, double rate, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_depreciation_amorlinc,
                        (IMSL_F cost, struct tm issue, struct tm first_period, 
                        IMSL_F salvage, int period, IMSL_F rate, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_depreciation_amorlinc,
                        (double cost, struct tm issue, struct tm first_period, 
                        double salvage, int period, double rate, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_discount_yield,
                        (struct tm settlement, struct tm maturity, IMSL_F price, 
                        IMSL_F redemption, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_discount_yield,
                        (struct tm settlement, struct tm maturity, double price, 
                        double redemption, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_yield_maturity,
                        (struct tm settlement, struct tm maturity, struct tm issue, 
                        IMSL_F rate, IMSL_F price, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_yield_maturity,
                        (struct tm settlement, struct tm maturity, struct tm issue, 
                        double rate, double price, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_convexity,
                        (struct tm settlement, struct tm maturity, IMSL_F coupon, 
                        IMSL_F yield, int frequency, int basis));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_convexity,
                        (struct tm settlement, struct tm maturity, double coupon, 
                        double yield, int frequency, int basis));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_yield_periodic,
                        (struct tm settlement, struct tm maturity, IMSL_F rate, 
                        IMSL_F price, IMSL_F redemption, int frequency, int basis, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_yield_periodic,
                        (struct tm settlement, struct tm maturity, double rate, 
                        double price, double redemption, int frequency, int basis, ...));
IMSL_CI Imsl_faure   *  IMSL_DECL IMSL_EF IMSL_PROTO(imsl_faure_sequence_init,
                        (int, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_faure_sequence_free,
                        (Imsl_faure*));
IMSL_CI IMSL_F        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_faure_next_point,
                        (Imsl_faure*, ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_faure_next_point,
                        (Imsl_faure*, ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_qmc,
                        (IMSL_F(IMSL_DECL *)(int, IMSL_F*), int, IMSL_F*, IMSL_F*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_qmc,
                        (double(IMSL_DECL *)(int, double*), int, double*, double*, ...));
/* multivariate normal cdf version of qmc */
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_int_fcn_qmc2,
                        (IMSL_F(IMSL_DECL *)(int, IMSL_F*), int, IMSL_F*, IMSL_F*, ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_int_fcn_qmc2,
                        (double(IMSL_DECL *)(int, double*), int, double*, double*, ...));

IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_constrained_nlp,
                        (void (*fcn) (int, float[], int, float*, int*),
                        int m, int me, int n, int ibtype,
                        float xlb[], float xub[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_constrained_nlp,
                        (void (*fcn) (int, double[], int, double*, int*),
                        int m, int me, int n, int ibtype,
                        double xlb[], double xub[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_bvp_finite_difference,
                        (void (*fcneqn) (int, float, float *, float, float *),
                        void (*fcnjac) (int, float, float *, float, float *),
                        void (*fcnbc)  (int, float *, float *, float, float *),
                        int neqns, int nleft, int ncupbc, float xleft,
                        float xright, int linear, int *nfinal, float *xfinal,
                        float *yfinal, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_bvp_finite_difference,
                        (void (*fcneqn) (int, double, double *, double, double *),
                        void (*fcnjac) (int, double, double *, double, double *),
                        void (*fcnbc)  (int, double *, double *, double, double *),
                        int neqns, int nleft, int ncupbc, double xleft,
                        double xright, int linear, int *nfinal, double *xfinal,
                        double *yfinal, ...));

/* V6.0 */
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_dea_petzold_gear_mgr,
			(int ido, char ** state, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_dea_petzold_gear_mgr,
			(int ido, char ** state, ...));
IMSL_CI int            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_dea_petzold_gear,
			(int neq, float * t, float tout,
			 float * y, float * ypr, char * state,
			 int (*gcn) (int, float, float *, float *, float *), ...));
IMSL_CI int            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_dea_petzold_gear,
			(int neq, double * t, double tout,
			 double * y, double * ypr, char * state,
			 int (*gcn) (int, double, double *, double *, double *), ...));
IMSL_CI Imsl_f_mps    * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_read_mps,
                        (char*, ...));
IMSL_CI Imsl_d_mps    * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_read_mps,
                        (char*, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_mps_free,
                        (Imsl_f_mps*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_mps_free,
                        (Imsl_d_mps*));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_pde_1d_mg_mgr,
			(int ido, char ** state, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_pde_1d_mg_mgr,
			(int ido, char ** state, ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_pde_1d_mg,
                        (int npdes,  int ngrids,  float * t, float tout,
		        float * u,  float xl, float xr, char * state,
		        void (*spdef) (float t, float x, int npde, int ngrids, 
                        float *full_u, float *grid_u, float *dudx, float *c, float *q, 
                        float *r, int *ires), 
                        void (*bndr) (float t, float *beta, float *gamma, 
                        float *full_u, float *grid_u, float *dudx, int npde, 
                        int ngrids, int left, int *ires), ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_pde_1d_mg,
                        (int npdes,  int ngrids,  double * t, double tout,
		        double * u,  double xl, double xr, char * state,
		        void (*spdef) (double t, double x, int npde, int ngrids,
                        double *full_u, double *grid_u, double *dudx, double *c, 
                        double *q, double *r, int *ires), 
                        void (*bndr) (double t, double *beta, double *gamma, 
                        double *full_u, double *grid_u, double *dudx, int npde, 
                        int ngrids, int left, int *ires), ...));
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_linear_programming_NA,
                        (int, int, float*, float[], float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_linear_programming,
                        (int, int, double*, double[], double[], ...));
IMSL_CI IMSL_F          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_min_uncon_golden,
                        (IMSL_F(IMSL_DECL *)(IMSL_F), IMSL_F, IMSL_F,
                         ...));
IMSL_CI double          IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_min_uncon_golden,
                        (double(IMSL_DECL *)(double), double, double,
                         ...));

/* end V6.0. */

/* V7.0 */
IMSL_CI float         * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_feynman_kac_evaluate,
                        (int, int, float[], float[], float[], ...));
IMSL_CI double        * IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_feynman_kac_evaluate,
                        (int, int, double[], double[], double[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_f_feynman_kac,
                        (int nxgrid, int ntgrid, int nlbcd,
                        int nrbcd, float xgrid[], float tgrid[],
                        void (*fcn_fkcfiv)(float, float, int *, float *),
                        void (*fcn_fkbcp)(int, float, int *, float[]),
                        float y[], float y_prime[], ...));
IMSL_CI void            IMSL_DECL IMSL_EF IMSL_PROTO(imsl_d_feynman_kac,
                        (int nxgrid, int ntgrid, int nlbcd,
                        int nrbcd, double xgrid[], double tgrid[],
                        void (*fcn_fkcfiv)(double, double, int *, double *),
                        void (*fcn_fkbcp)(int, double, int *, double[]),
                        double y[], double y_prime[], ...));
/* end V7.0 */

IMSL_CI int             IMSL_DECL IMSL_EF IMSL_PROTO(imsl_initialize,(int));

#if defined(__cplusplus)
}                               /* end extern "C" for C++ */
#endif

        /* Keywords */

enum Imsl_keyword {
    IMSL_TRANSPOSE                    = 10001,
    IMSL_RESULT                       = 10002,
    IMSL_A_COL_DIM                    = 10003,
    IMSL_FACTOR                       = 10004,
    IMSL_FAC_COL_DIM                  = 10005,
    IMSL_FACTOR_ONLY                  = 10006,
    IMSL_SOLVE_ONLY                   = 10007,
    IMSL_BACKWARD                     = 10008,
    IMSL_PARAMS                       = 10009,
    IMSL_ERR_ABS                      = 10010,
    IMSL_ERR_REL                      = 10011,
    IMSL_ETA                          = 10012,
    IMSL_EPS                          = 10013,
    IMSL_GUESS                        = 10014,
    IMSL_ITMAX                        = 10016,
    IMSL_INFO                         = 10017,
    IMSL_NUM_ROOTS                    = 10018,
    IMSL_RULE                         = 10019,
    IMSL_ERR_EST                      = 10020,
    IMSL_MAX_SUBINTER                 = 10021,
    IMSL_N_SUBINTER                   = 10022,
    IMSL_N_EVALS                      = 10023,
    IMSL_ERR_LIST                     = 10024,
    IMSL_ERR_ORDER                    = 10025,
    IMSL_BREAKPOINTS                  = 10026,
    IMSL_COEFS                        = 10027,
    IMSL_DERIV                        = 10028,
    IMSL_CONCAVE                      = 10029,
    IMSL_PERIODIC                     = 10030,
    IMSL_LEFT                         = 10031,
    IMSL_RIGHT                        = 10032,
    IMSL_WEIGHT                       = 10033,
    IMSL_SMPAR                        = 10034,
    IMSL_KNOTS                        = 10035,
    IMSL_ORDER                        = 10036,
    IMSL_OPT                          = 10037,
    IMSL_MIN_PROJECTION               = 10038,
    IMSL_KNOTS_USER                   = 10039,
    IMSL_PRINT_ALL                    = 10040,
    IMSL_ROW_LABELS                   = 10042,
    IMSL_COL_LABELS                   = 10043,
    IMSL_WRITE_FORMAT                 = 10044,
    IMSL_X_COL_DIM                    = 10045,
    IMSL_COV_COL_DIM                  = 10046,
    IMSL_COEF_COVARIANCES             = 10047,
    IMSL_COEF_COVARIANCES_USER        = 10048,
    IMSL_RANK                         = 10049,
    IMSL_NO_INTERCEPT                 = 10050,
    IMSL_ANOVA_TABLE                  = 10051,
    IMSL_ANOVA_TABLE_USER             = 10052,
    IMSL_TOLERANCE                    = 10053,
    IMSL_FREQUENCIES                  = 10054,
    IMSL_FREQUENCIES_USER             = 10055,
    IMSL_BOUNDS                       = 10056,
    IMSL_N_PARAMETERS_ESTIMATED       = 10057,
    IMSL_CUTPOINTS                    = 10058,
    IMSL_CUTPOINTS_USER               = 10059,
    IMSL_CELL_COUNTS                  = 10060,
    IMSL_CELL_COUNTS_USER             = 10061,
    IMSL_CELL_EXPECTED                = 10062,
    IMSL_CELL_EXPECTED_USER           = 10063,
    IMSL_CELL_CHI_SQUARED             = 10064,
    IMSL_CELL_CHI_SQUARED_USER        = 10065,
    IMSL_DEGREES_OF_FREEDOM           = 10066,
    IMSL_TIES_OPTION                  = 10067,
    IMSL_FUZZ                         = 10068,
    IMSL_SCORE_OPTION                 = 10069,
    IMSL_RESULT_USER                  = 10070,
    IMSL_NORM                         = 10071,
    IMSL_TOL                          = 10072,
    IMSL_HINIT                        = 10073,
    IMSL_HMIN                         = 10074,
    IMSL_SCALE                        = 10075,
    IMSL_FLOOR                        = 10076,
    IMSL_MAX_NUMBER_STEPS             = 10077,
    IMSL_MAX_NUMBER_FCN_EVALS         = 10078,
    IMSL_INTERRUPT_1                  = 10079,
    IMSL_INTERRUPT_2                  = 10080,
    IMSL_NSTEP                        = 10081,
    IMSL_NFCN                         = 10082,
    IMSL_HTRIAL                       = 10083,
    IMSL_VNORM                        = 10084,
    IMSL_HMAX                         = 10085,
    IMSL_STAT_COL_DIM                 = 10086,
    IMSL_CONFIDENCE_MEANS             = 10087,
    IMSL_CONFIDENCE_VARIANCES         = 10088,
    IMSL_MEANS                        = 10090,
    IMSL_MEANS_USER                   = 10091,
    IMSL_COVARIANCE_COL_DIM           = 10092,
    IMSL_COMPUTE_OPTION               = 10093,
    IMSL_VECTORS                      = 10094,
    IMSL_VECTORS_USER                 = 10095,
    IMSL_EVECU_COL_DIM                = 10096,
    IMSL_RANGE                        = 10097,
    IMSL_INFO_USER                    = 10099,
    IMSL_XGUESS                       = 10100,
    IMSL_STEP                         = 10101,
    IMSL_BOUND                        = 10102,
    IMSL_MAX_FCN                      = 10103,
    IMSL_INIT_TRUST_REGION            = 10104,
    IMSL_GRAD                         = 10105,
    IMSL_XSCALE                       = 10106,
    IMSL_FSCALE                       = 10107,
    IMSL_GRAD_TOL                     = 10108,
    IMSL_STEP_TOL                     = 10109,
    IMSL_REL_FCN_TOL                  = 10110,
    IMSL_MAX_STEP                     = 10111,
    IMSL_GOOD_DIGIT                   = 10112,
    IMSL_MAX_ITN                      = 10113,
    IMSL_MAX_GRAD                     = 10114,
    IMSL_INIT_HESSIAN                 = 10115,
    IMSL_FVALUE                       = 10116,
    IMSL_GVALUE                       = 10117,
    IMSL_JACOBIAN                     = 10118,
    IMSL_FNORM                        = 10119,
    IMSL_RESULT_DUAL                  = 10120,
    IMSL_UPPER_LIMIT                  = 10121,
    IMSL_CONSTR_TYPE                  = 10122,
    IMSL_LOWER_BOUND                  = 10123,
    IMSL_UPPER_BOUND                  = 10124,
    IMSL_OBJ                          = 10125,
    IMSL_DUAL_USER                    = 10126,
    IMSL_DUAL                         = 10127,
    IMSL_JAC_TOL                      = 10128,
    IMSL_H_COL_DIM                    = 10130,
    IMSL_ADD_TO_DIAG_H                = 10131,
    IMSL_FDATA_COL_DIM                = 10140,
    IMSL_WEIGHTS                      = 10141,
    IMSL_ERROR_SSQ                    = 10142,
    IMSL_OPTIMIZE                     = 10143,
    IMSL_INTERCEPT                    = 10144,
    IMSL_SSE                          = 10145,
    IMSL_SMOOTHING_PAR                = 10146,
    IMSL_SUR_COL_DIM                  = 10147,
    IMSL_OPT_ITMAX                    = 10148,
    IMSL_CONCAVE_ITMAX                = 10149,
    IMSL_SOLUTION_USER                = 10150,
    IMSL_FACTOR_USER                  = 10151,
    IMSL_INVERSE                      = 10152,
    IMSL_INVERSE_USER                 = 10153,
    IMSL_INV_COL_DIM                  = 10154,
    IMSL_INVERSE_ONLY                 = 10155,
    IMSL_SSQ_POLY                     = 10156,
    IMSL_SSQ_LOF                      = 10157,
    IMSL_X_MEAN                       = 10158,
    IMSL_Y_MEAN                       = 10159,
    IMSL_X_VARIANCE                   = 10160,
    IMSL_Y_VARIANCE                   = 10161,
    IMSL_BASIS                        = 10170,
    IMSL_PIVOT                        = 10171,
    IMSL_RESIDUAL                     = 10172,
    IMSL_RESIDUAL_USER                = 10173,
    IMSL_Q                            = 10177,
    IMSL_Q_USER                       = 10178,
    IMSL_Q_COL_DIM                    = 10179,
    IMSL_P_COL_DIM                    = 10180,
    IMSL_A_MATRIX                     = 10181,
    IMSL_B_MATRIX                     = 10182,
    IMSL_X_VECTOR                     = 10183,
    IMSL_Y_VECTOR                     = 10184,
    IMSL_RETURN_COL_DIM               = 10185,
    IMSL_B_COL_DIM                    = 10186,
    IMSL_INVA_COL_DIM                 = 10187,
    IMSL_SET_PRINT                    = 10188,
    IMSL_GET_PRINT                    = 10189,
    IMSL_SET_STOP                     = 10190,
    IMSL_GET_STOP                     = 10191,
    IMSL_GET_TRACEBACK                = 10192,
    IMSL_SET_TRACEBACK                = 10193,
    IMSL_SET_ERROR_FILE               = 10194,
    IMSL_ERROR_PRINT_PROC             = 10195,
    IMSL_ERROR_MSG_PATH               = 10196,
    IMSL_ERROR_MSG_NAME               = 10197,
    IMSL_S_USER                       = 10198,
    IMSL_U                            = 10199,
    IMSL_U_USER                       = 10200,
    IMSL_U_COL_DIM                    = 10201,
    IMSL_V                            = 10202,
    IMSL_V_USER                       = 10203,
    IMSL_V_COL_DIM                    = 10204,
    IMSL_SET_OUTPUT_FILE              = 10208,
    IMSL_GET_OUTPUT_FILE              = 10209,
    IMSL_GET_ERROR_FILE               = 10210,
    IMSL_PRINT_LOWER                  = 10211,
    IMSL_PRINT_UPPER                  = 10212,
    IMSL_PRINT_LOWER_NO_DIAG          = 10213,
    IMSL_PRINT_UPPER_NO_DIAG          = 10214,
    IMSL_ROW_NUMBER_ZERO              = 10215,
    IMSL_NO_ROW_LABELS                = 10216,
    IMSL_COL_NUMBER_ZERO              = 10217,
    IMSL_NO_COL_LABELS                = 10218,
    IMSL_VARIANCE_COVARIANCE_MATRIX   = 10219,
    IMSL_CORRECTED_SSCP_MATRIX        = 10220,
    IMSL_CORRELATION_MATRIX           = 10221,
    IMSL_STDEV_CORRELATION_MATRIX     = 10222,
    IMSL_AVERAGE_TIE                  = 10223,
    IMSL_HIGHEST                      = 10224,
    IMSL_LOWEST                       = 10225,
    IMSL_RANDOM_SPLIT                 = 10226,
    IMSL_RANKS                        = 10227,
    IMSL_BLOM_SCORES                  = 10228,
    IMSL_TUKEY_SCORES                 = 10229,
    IMSL_VAN_DER_WAERDEN_SCORES       = 10230,
    IMSL_EXPECTED_NORMAL_SCORES       = 10231,
    IMSL_SAVAGE_SCORES                = 10232,
    IMSL_CHEBYSHEV_FIRST              = 10240,
    IMSL_CHEBYSHEV_SECOND             = 10241,
    IMSL_HERMITE                      = 10242,
    IMSL_COSH                         = 10243,
    IMSL_JACOBI                       = 10244,
    IMSL_GEN_LAGUERRE                 = 10245,
    IMSL_FIXED_POINT                  = 10246,
    IMSL_TWO_FIXED_POINTS             = 10247,
    IMSL_LEGENDRE                     = 10248,
    IMSL_GRADIENT                     = 10250,
    IMSL_PRINT                        = 10251,
    IMSL_FJAC_Q                       = 10253,
    IMSL_FJAC_Q_USER                  = 10254,
    IMSL_FJAC_R                       = 10255,
    IMSL_FJAC_R_USER                  = 10256,
    IMSL_RETURN_NUMBER                = 10259,
    IMSL_RETURN_USER                  = 10260,
    IMSL_X_MEAN_USER                  = 10261,
    IMSL_DF_PURE_ERROR                = 10262,
    IMSL_SSQ_PURE_ERROR               = 10263,
    IMSL_SSQ_POLY_USER                = 10264,
    IMSL_SSQ_POLY_COL_DIM             = 10265,
    IMSL_SSQ_LOF_USER                 = 10266,
    IMSL_SSQ_LOF_COL_DIM              = 10267,
    IMSL_ROW_NUMBER                   = 10268,
    IMSL_COL_NUMBER                   = 10269,
    IMSL_CONDITION                    = 10270,
    IMSL_MAX_MOMENTS                  = 10271,
    IMSL_OS_VERSION                   = 10272,
    IMSL_COMPILER_VERSION             = 10273,
    IMSL_LIBRARY_VERSION              = 10274,
    IMSL_MAX_CYCLES                   = 10275,
    IMSL_N_CYCLES                     = 10276,
    IMSL_MAX_EVALS                    = 10277,
    IMSL_CUTPOINTS_EQUAL              = 10280,
    IMSL_ABS_FCN_TOL                  = 10290,
    IMSL_MAX_JACOBIAN                 = 10291,
    IMSL_INTERN_SCALE                 = 10292,
    IMSL_FVEC                         = 10293,
    IMSL_FVEC_USER                    = 10294,
    IMSL_FJAC                         = 10295,
    IMSL_FJAC_USER                    = 10296,
    IMSL_FJAC_COL_DIM                 = 10297,
    IMSL_JTJ_INVERSE                  = 10298,
    IMSL_JTJ_INVERSE_USER             = 10299,
    IMSL_JTJ_INV_COL_DIM              = 10300,
    IMSL_FULL_TRACEBACK               = 10301,
    IMSL_MAX_ITER                     = 10305,
    IMSL_PRECOND                      = 10306,
    IMSL_REL_ERR                      = 10307,
    IMSL_LICENSE_NUMBER               = 10308,
    IMSL_METHOD                       = 10309,
    IMSL_MAXORD                       = 10310,
    IMSL_MITER                        = 10311,
    IMSL_NFCNJ                        = 10312,
    IMSL_MEDIAN                       = 10320,
    IMSL_MEDIAN_AND_SCALE             = 10321,
    IMSL_CHI_SQUARED                  = 10322,
    IMSL_ONE_NORM                     = 10323,
    IMSL_INF_NORM                     = 10324,
    IMSL_SECOND_VECTOR                = 10325,
    IMSL_ABSOLUTE                     = 10326,
    IMSL_PERMUTATION                  = 10327,
    IMSL_PERMUTATION_USER             = 10328,
    IMSL_MAXIMIZATION                 = 10329,
    IMSL_DATA_BOUNDS                  = 10330,
    IMSL_KNOWN_BOUNDS                 = 10331,
    IMSL_CLASS_MARKS                  = 10332,
    IMSL_KNOWN_BOUNDS_ADR             = 10333,
    IMSL_RANK_ADR                     = 11001,
    IMSL_BASIS_ADR                    = 11002,
    IMSL_TOLERANCE_ADR                = 11003,
    IMSL_BOUNDS_ADR                   = 11004,
    IMSL_FUZZ_ADR                     = 11005,
    IMSL_CONFIDENCE_MEANS_ADR         = 11006,
    IMSL_CONFIDENCE_VARIANCES_ADR     = 11007,
    IMSL_RANGE_ADR                    = 11008,
    IMSL_LEFT_ADR                     = 11009,
    IMSL_RIGHT_ADR                    = 11010,
    IMSL_SMOOTHING_PAR_ADR            = 11012,
    IMSL_JACOBI_ADR                   = 11013,
    IMSL_GEN_LAGUERRE_ADR             = 11014,
    IMSL_FIXED_POINT_ADR              = 11015,
    IMSL_TWO_FIXED_POINTS_ADR         = 11016,
    IMSL_TOL_ADR                      = 11017,
    IMSL_HINIT_ADR                    = 11018,
    IMSL_HMIN_ADR                     = 11019,
    IMSL_HMAX_ADR                     = 11020,
    IMSL_SCALE_ADR                    = 11021,
    IMSL_FLOOR_ADR                    = 11022,
    IMSL_ETA_ADR                      = 11023,
    IMSL_EPS_ADR                      = 11024,
    IMSL_FSCALE_ADR                   = 11025,
    IMSL_GRAD_TOL_ADR                 = 11026,
    IMSL_STEP_TOL_ADR                 = 11027,
    IMSL_REL_FCN_TOL_ADR              = 11028,
    IMSL_MAX_STEP_ADR                 = 11029,
    IMSL_XGUESS_ADR                   = 11030,
    IMSL_STEP_ADR                     = 11031,
    IMSL_ABS_FCN_TOL_ADR              = 11033,
    IMSL_INIT_TRUST_REGION_ADR        = 11034,
    IMSL_ERR_REL_ADR                  = 11035,
    IMSL_ERR_ABS_ADR                  = 11036,
    IMSL_GAP                          = 11037,
    IMSL_ZBOUND                       = 11038,
    IMSL_BB_ONLY                      = 11039,
    IMSL_NODE_SELECT                  = 11040,
    IMSL_GRID                         = 11050,
    IMSL_GRID_USER                    = 11051,
    IMSL_LINEAR_TERM                  = 11060,
    IMSL_CONSTANT_TERM                = 11061,
    IMSL_SUPPLY_BASIS                 = 11062,
    IMSL_SUPPLY_DELTA                 = 11063,
    IMSL_SUPPLY_DELTA_ADR             = 11064,
    IMSL_CENTERS                      = 11065,
    IMSL_CENTERS_RATIO                = 11066,
    IMSL_CENTERS_RATIO_ADR            = 11067,
    IMSL_RANDOM_SEED                  = 11068,
    IMSL_CORRELATION                  = 11070,
    IMSL_Z_TRANS                      = 11071,
    IMSL_Z_TRANS_USER                 = 11072,
    IMSL_FIRST_CALL                   = 11073,
    IMSL_CONTINUE_CALL                = 11074,
    IMSL_LAST_CALL                    = 11075,
    IMSL_NHARD                        = 11080,
    IMSL_NO_SVD                       = 11081,
    IMSL_COMPANION_METHOD             = 11082,
    IMSL_RETURN_SYMBOLIC_FACTOR       = 11083,
    IMSL_SYMBOLIC_FACTOR_ONLY         = 11084,
    IMSL_SUPPLY_SYMBOLIC_FACTOR       = 11085,
    IMSL_RETURN_NUMERIC_FACTOR        = 11086,
    IMSL_SUPPLY_NUMERIC_FACTOR        = 11087,
    IMSL_NUMERIC_FACTOR_ONLY          = 11088,
    IMSL_MULTIFRONTAL_FACTORIZATION   = 11089,
    IMSL_SMALLEST_DIAGONAL_ELEMENT    = 11090,
    IMSL_LARGEST_DIAGONAL_ELEMENT     = 11091,
    IMSL_NUM_NONZEROS_IN_FACTOR       = 11092,
    IMSL_D_MATRIX                     = 11093,
    IMSL_SYMMETRIC_STORAGE            = 11094,
    IMSL_BAND_FORMAT                  = 11095,
    IMSL_CHANGE_LOOP_MAXIMUM          = 11096,
    IMSL_RETURN_SPARSE_LU_FACTOR      = 11097,
    IMSL_SUPPLY_SPARSE_LU_FACTOR      = 11098,
    IMSL_NUMBER_OF_SEARCH_ROWS        = 11099,
    IMSL_ITERATIVE_REFINEMENT         = 11100,
    IMSL_DROP_TOLERANCE               = 11101,
    IMSL_STABILITY_FACTOR             = 11102,
    IMSL_GROWTH_FACTOR_LIMIT          = 11103,
    IMSL_GROWTH_FACTOR                = 11104,
    IMSL_SMALLEST_PIVOT               = 11105,
    IMSL_MEMORY_BLOCK_SIZE            = 11106,
    IMSL_CSC_FORMAT                   = 11107,
    IMSL_FREE_SPARSE_LU_FACTOR        = 11108,
    IMSL_RETURN_SPARSE_LU_IN_COORD    = 11109,
    IMSL_SUPPLY_SPARSE_LU_IN_COORD    = 11110,
    IMSL_PIVOTING_STRATEGY            = 11111,
    IMSL_C_MATRIX                     = 11112,
    IMSL_HYBRID_FACTORIZATION         = 11113,
    IMSL_RETURN_MATRIX_SIZE           = 11114,
    IMSL_RETURN_USER_VECTOR           = 11115,
    IMSL_DROP_TOLERANCE_ADR           = 11116,
    IMSL_HYBRID_FACTORIZATION_ADR     = 11117,
    IMSL_GROWTH_FACTOR_LIMIT_ADR      = 11118,
    IMSL_BLOCKING_FACTOR              = 11119,
    IMSL_MAX_KRYLOV_SUBSPACE_DIM      = 11120,
    IMSL_HOUSEHOLDER_REORTHOG         = 11121,
    IMSL_RESIDUAL_NORM                = 11122,
    IMSL_STABILITY_FACTOR_ADR         = 11123,
    IMSL_PSEUDO_ACCURACY              = 11124,
    IMSL_FIRST_LAGUERRE_PARAMETER     = 11125,
    IMSL_SECOND_LAGUERRE_PARAMETER    = 11126,
    IMSL_MAXIMUM_COEFFICIENTS         = 11127,
    IMSL_ERROR_EST                    = 11128,
    IMSL_DISCRETIZATION_ERROR_EST     = 11129,
    IMSL_TRUNCATION_ERROR_EST         = 11130,
    IMSL_CONDITION_ERROR_EST          = 11131,
    IMSL_DECAY_FUNCTION_COEFFICIENT   = 11132,
    IMSL_DECAY_FUNCTION_BASE          = 11133,
    IMSL_LOG_LARGEST_COEFFICIENTS     = 11134,
    IMSL_LOG_SMALLEST_COEFFICIENTS    = 11135,
    IMSL_UNDER_OVERFLOW_INDICATORS    = 11136,
    IMSL_INITIAL_STEPSIZE             = 11137,
    IMSL_RELATIVE_ERROR               = 11138,
    IMSL_PERFORMANCE_INDEX            = 11139,
    IMSL_ACTIVE_CONSTRAINTS           = 11140,
    IMSL_ACTIVE_CONSTRAINTS_USER      = 11141,
    IMSL_NUMBER_ACTIVE_CONSTRAINTS    = 11142,
    IMSL_LAGRANGE_MULTIPLIERS         = 11143,
    IMSL_LAGRANGE_MULTIPLIERS_USER    = 11144,
    IMSL_A_TRANSPOSE                  = 11145,
    IMSL_B_TRANSPOSE                  = 11146,
    IMSL_A_CONJUGATE_TRANSPOSE        = 11147,
    IMSL_B_CONJUGATE_TRANSPOSE        = 11148,
    IMSL_RETURN_MATRIX_CODIAGONALS    = 11149,
    IMSL_SYMMETRIC                    = 11150,
    IMSL_INITIAL_VALUE_DERIVATIVE     = 11151,
    IMSL_HERMITIAN                    = 11152,
    IMSL_REL_ERR_ADR                  = 11200,
    IMSL_DISTANCE                     = 12001,
    IMSL_STOPPING_CRITERION           = 12002,
    IMSL_BASE                         = 12100,
    IMSL_SKIP                         = 12101,
    IMSL_RETURN_SKIP                  = 12102,
    IMSL_TAU0                         = 13001,
    IMSL_DEL0                         = 13002,
    IMSL_SMALLW                       = 13003,
    IMSL_DELMIN                       = 13004,
    IMSL_SCFMAX                       = 13005,
    IMSL_EPSDIF                       = 13006,
    IMSL_EPSFCN                       = 13007,
    IMSL_TAUBND                       = 13008,
    IMSL_DIFFTYPE                     = 13009,
    IMSL_TAU0_ADR                     = 13010,
    IMSL_DEL0_ADR                     = 13011,
    IMSL_SMALLW_ADR                   = 13012,
    IMSL_DELMIN_ADR                   = 13013,
    IMSL_SCFMAX_ADR                   = 13014,
    IMSL_EPSDIF_ADR                   = 13015,
    IMSL_EPSFCN_ADR                   = 13016,
    IMSL_TAUBND_ADR                   = 13017,
    IMSL_GCN_W_DATA                   = 13100,
    IMSL_FCN_W_DATA                   = 13101,
    IMSL_GRADIENT_W_DATA              = 13102,
    IMSL_PRECOND_W_DATA               = 13103,
    IMSL_JACOBIAN_W_DATA              = 13104,
    IMSL_RHS_PDE_W_DATA               = 13105,
    IMSL_RHS_BC_W_DATA                = 13106,
    IMSL_FCN_UT_W_DATA                = 13107,
    IMSL_FCN_BC_W_DATA                = 13108,
    IMSL_HCN_W_DATA                   = 13109,
    IMSL_ERR_EST_USER                 = 13110,
    IMSL_PROBLEM_EMBEDDED             = 13111,
    IMSL_SUPPLY_BASIS_W_DATA          = 13112,
    IMSL_PROBLEM_EMBEDDED_W_DATA      = 13113,
    IMSL_CONSTRAINT_RESIDUALS         = 14001,
    IMSL_CONSTRAINT_RESIDUALS_USER    = 14002,
    IMSL_ATOL_RTOL_ARRAYS             = 14010,
    IMSL_RETURN_AT_INTERNAL_STEP      = 14011,
    IMSL_INITIAL_STEPSIZE_ADR         = 14012,
    IMSL_ALL_NONNEGATIVE              = 14013,
    IMSL_INITIAL_VALUES_INCONSISTENT  = 14014,
    IMSL_ATOL_RTOL_SCALARS            = 14017,
    IMSL_ATOL_RTOL_SCALARS_ADR        = 14018,
    IMSL_MAX_MAGNITUDE_CONSTRAINT     = 14019,
    IMSL_MAX_MAGNITUDE_CONSTRAINT_ADR = 14020,
    IMSL_MAX_BDF_ORDER                = 14022,
    IMSL_NERROR_TEST_FAILURES         = 14024,
    IMSL_NCONV_TEST_FAILURES          = 14025,
    IMSL_BDF_ORDER_NEXT_STEP          = 14026,
    IMSL_BDF_ORDER_PREVIOUS_STEP      = 14027,
    IMSL_NSTEPS_TAKEN                 = 14028,
    IMSL_T_BARRIER                    = 14029,
    IMSL_T_BARRIER_ADR                = 14030,
    IMSL_NORM_FCN                     = 14031,
    IMSL_NORM_FCN_W_DATA              = 14032,
    IMSL_USER_JAC_FACTOR_SOLVE        = 14033,
    IMSL_USER_JAC_FACTOR_SOLVE_W_DATA = 14034,
    IMSL_CART_COORDINATES             = 14035,
    IMSL_CYL_COORDINATES              = 14036,
    IMSL_SPH_COORDINATES              = 14037,
    IMSL_TIME_SMOOTHING               = 14038,
    IMSL_SPATIAL_SMOOTHING            = 14039,
    IMSL_MONITOR_REGULARIZING         = 14040,
    IMSL_RELATIVE_TOLERANCE           = 14041,
    IMSL_ABSOLUTE_TOLERANCE           = 14042,
    IMSL_PDE_SYS_W_DATA               = 14043,
    IMSL_BOUNDARY_COND_W_DATA         = 14044,
    IMSL_GRID_VALUES                  = 14045,
    IMSL_USER_FACTOR_SOLVE            = 14046,
    IMSL_USER_FACTOR_SOLVE_W_DATA     = 14047,
    IMSL_INITIAL_CONDITIONS           = 14048,
    IMSL_INITIAL_CONDITIONS_W_DATA    = 14049,
    IMSL_FILE                         = 15001,
    IMSL_NAME_OBJECTIVE               = 15002,
    IMSL_NAME_RHS                     = 15003,
    IMSL_NAME_RANGES                  = 15004,
    IMSL_NAME_BOUNDS                  = 15005,
    IMSL_POSITIVE_INFINITY            = 15006,
    IMSL_NEGATIVE_INFINITY            = 15007,

    IMSL_ITERATION_COUNT              = 15010,
    IMSL_FEASIBLE_PT_ONLY             = 15011,
    IMSL_REFINEMENT                   = 15012,
    IMSL_EXTENDED_REFINEMENT          = 15013,
    IMSL_PRODUCE_F77_TEST             = 15014,
    IMSL_USE_UPDATED_LP_ALGORITHM     = 15015,

    IMSL_TENSION                      = 15016,
    IMSL_CONTINUITY                   = 15017,
    IMSL_BIAS                         = 15018,

    IMSL_LOWER_ENDPOINT               = 15019,
    IMSL_UPPER_ENDPOINT               = 15020,
    IMSL_INTERVAL                     = 15021,

    IMSL_TOLERANCE_MULLER             = 15022,
    IMSL_MIN_SEPARATION               = 15023,
    IMSL_ERR_X                        = 15024,
    IMSL_NUM_ROOTS_FOUND              = 15025,

    IMSL_SET_FUNCTIONS_THREAD_SAFE    = 15026,
    IMSL_GET_FUNCTIONS_THREAD_SAFE    = 15027,

	IMSL_FCN_FKCFIV_W_DATA            = 15028,       
    IMSL_FCN_FKBCP_W_DATA             = 15029,
    IMSL_FCN_INIT                     = 15030,
    IMSL_FCN_INIT_W_DATA              = 15031,
    IMSL_FCN_FORCE                    = 15032,
    IMSL_FCN_FORCE_W_DATA             = 15033,
    IMSL_NDEGREE                      = 15034,
    IMSL_TDEPEND                      = 15035,
    IMSL_STEP_CONTROL                 = 15036,
    IMSL_ISTATE                       = 15037,
    IMSL_EVALS                        = 15038,
    IMSL_GET_OMP_MKL_SINGLE_THREAD    = 15039,
    IMSL_SET_OMP_MKL_SINGLE_THREAD    = 15040,

    IMSL_GET_PRINT_TYPE               = 50000,
    IMSL_SET_PRINT_TYPE               = 50001,
    IMSL_DEFAULT_PRINT_PROC           = 50002,
    IMSL_WINDOWS_PRINT_PROC           = 50003,
    IMSL_CONSOLE_PRINT_PROC           = 50004,
    IMSL_RETURN_STRING                = 50010,
    IMSL_WRITE_TO_CONSOLE             = 50011,
    IMSL_SET_SIGNAL_TRAPPING          = 50012
};
#endif /* IMSL_H */