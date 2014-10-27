/*      Copyright:  2005 by IMSL, Inc.  All Rights Reserved. */

#ifndef IMSLS_H
#define IMSLS_H


/*
 * uint32_t and uint64_t are 32-bit and 64-bit unsigned integers.
 * These types are defined in the C99 standard, but not yet in all
 * compilers.
 * UINT32_C and UINT64_C are macros which append the correct suffix to a
 * literal integer to make it a 32-bit or 64-bit unsigned literal integer.
 */
#if defined(COMPUTER_VMSIA64) || defined(COMPUTER_VAXG)
    typedef unsigned int uint32_t;
    typedef unsigned long long uint64_t;
#   define UINT32_C(VALUE)   VALUE##UL
#   define UINT64_C(VALUE)   VALUE##ULL
#elif defined(_WIN32)
    typedef unsigned int uint32_t;
    typedef unsigned long long uint64_t;
#   define UINT32_C(VALUE)   VALUE##ui32
#   define UINT64_C(VALUE)   VALUE##ui64
#else
   /* Needed uint32_t and uint64_t. */
#  include <inttypes.h>
#endif
/* Used in Mersenne twister code(s). */
#if defined(IMSL_WAVE_AXPOSF)
#   define UINT32_C(VALUE)   VALUE##LL
#   define UINT64_C(VALUE)   VALUE##LL
#endif

#include "imsls_e.h"

#ifdef _OPENMP
#  include <omp.h>
#endif

#ifdef WAVE_RENAME
#  include "renames.h"
#endif

#ifdef __sgi
#  include <sgidefs.h>
#endif

#ifdef vax
#    define imsls_d_random_normal_multivariate imsls_d_random_normal_multi
#  ifndef DOUBLE
#    define imsls_f_random_normal_multivariate imsls_f_random_normal_multi
#  endif /* !DOUBLE */
#endif /* vax */

/* IMSL C/Stat/Library version 7.0.0 */
#define IMSL_CSTAT_VERSION      700

#define IMSLS_PROTO(P,Q)        P Q

#if defined(_WIN32) || defined(__WIN32__)
#  define IMSLS_DECL    __cdecl
#  ifdef __BORLANDC__
#    ifdef _RTLDLL
#      ifndef IMSLS_EF
#        ifdef _EXPFUNC
#          define IMSLS_EF      _EXPFUNC
#        else
#          define IMSLS_EF      __import
#        endif
#        ifdef _EXPDATA
#          define IMSLS_ED      _EXPDATA
#        else
#          define IMSLS_ED      __import
#        endif
#      endif
#    else
#      define IMSLS_EF
#      define IMSLS_ED
#    endif
#    define IMSLS_CI
#  else
     /* Visual C++ */
#    if defined(_DLL) && !defined(IMSL_STATIC)
#      ifndef IMSLS_CI
#        ifdef _CRTIMP
#          define IMSLS_CI      _CRTIMP
#        else
#          define IMSLS_CI      __declspec(dllimport)
#        endif
#      endif
#    else
#      define IMSLS_CI
#    endif
#    define IMSLS_EF
#    define IMSLS_ED
#  endif
#else
#  define IMSLS_DECL
#  define IMSLS_CI
#  define IMSLS_EF
#  define IMSLS_ED
#endif

#ifdef __cplusplus              /* C++ compiler is being used */
extern "C" {
#endif

#define IMSLS_F                       float

/* Do not redeclare f_complex and d_complex if C/Math has already done so. */
#ifndef IMSL_H 
typedef struct {
    float                             re;
    float                             im;
} f_complex;
typedef struct {
    double                            re;
    double                            im;
} d_complex;
#endif


/* Piecewise polynomial structures. */
typedef struct {
    int                               domain_dim;
    int                               target_dim;
    int                               *order;
    int                               *num_coef;
    int                               *num_breakpoints;
    float                             **breakpoints;
    float                             **coef;
} Imsls_f_ppoly;

typedef struct {
    int                               domain_dim;
    int                               target_dim;
    int                               *order;
    int                               *num_coef;
    int                               *num_breakpoints;
    double                            **breakpoints;
    double                            **coef;
} Imsls_d_ppoly;


typedef struct {
    int                               nobs;
    int                               ncol;
    int                               model;
    int                               ifix;
    int                               intcep;
    int                               nclvar;
    int                               *indcl;
    int                               nef;
    int                               *nvef;
    int                               *indef;
    int                               ncoef;
    int                               *nclval;
    float                             *clval;
    float                             *coef;   
} Imsls_f_survival;
typedef struct {
    int                               nobs;
    int                               ncol;
    int                               model;
    int                               ifix;
    int                               intcep;
    int                               nclvar;
    int                               *indcl;
    int                               nef;
    int                               *nvef;
    int                               *indef;
    int                               ncoef;
    int                               *nclval;
    double                            *clval;
    double                            *coef;   
} Imsls_d_survival;

typedef struct {
    int                               n_observations;
    int                               n_cases_missing;
    int                               n_dependent;
    int                               n_parameters;
    int                               intercept;
    int                               rank;
    float                             df_error;
    float                             *beta_hat;
    float                             *sscp_error;
    float                             *rxx_transpose;
    float                             *d;
    float                             *rxy_transpose;
    float                             *x_mean;
    float                             *x_minimum;
    float                             *x_maximum;
    float                             *y_mean;
    float                             *y_minimum;
    float                             *y_maximum;
} Imsls_f_regression;
typedef struct {
    int                               n_observations;
    int                               n_cases_missing;
    int                               n_dependent;
    int                               n_parameters;
    int                               intercept;
    int                               rank;
    double                            df_error;
    double                            *beta_hat;
    double                            *sscp_error;
    double                            *rxx_transpose;
    double                            *d;
    double                            *rxy_transpose;
    double                            *x_mean;
    double                            *x_minimum;
    double                            *x_maximum;
    double                            *y_mean;
    double                            *y_minimum;
    double                            *y_maximum;
} Imsls_d_regression;

/* enums for arma options      */
typedef enum {
    IMSLS_CONST                       = 1,
    IMSLS_NO_CONST                    = 2
} Imsls_arma_method;

typedef struct {
    float                             a_variance;
    float                             constant;
    float                             *z;
    float                             *ar;
    float                             *ma;
    int                               *ar_lags;
    int                               *ma_lags;
    int                               n_observations;
    int                               p;
    int                               q;
    Imsls_arma_method                 constant_option;
} Imsls_f_arma;
typedef struct {
    double                            a_variance;
    double                            constant;
    double                            *z;
    double                            *ar;
    double                            *ma;
    int                               *ar_lags;
    int                               *ma_lags;
    int                               n_observations;
    int                               p;
    int                               q;
    Imsls_arma_method                 constant_option;
} Imsls_d_arma;

/* Structure for Wilcoxon Rank Sum Linked List */
struct Imsls_f_uList_struct{
	  int m;
	  int u;
	  float prob;
	  struct Imsls_f_uList_struct *next;
};
typedef struct Imsls_f_uList_struct Imsls_f_uList;
struct Imsls_d_uList_struct{
	  int m;
	  int u;
	  double prob;
	  struct Imsls_d_uList_struct *next;
};
typedef struct Imsls_d_uList_struct Imsls_d_uList;

/* structure for naive bayes classification */
typedef struct
{
	float *avg;        /* length=n_continuous*n_classes */
	float *s;          /* length=n_continuous*n_classes */
	float *logMeans;   /* length=n_continuous*n_classes */
	float *logStdev;   /* length=n_continuous*n_classes */
	float *a;          /* length=n_continuous*n_classes */
	float *b;          /* length=n_continuous*n_classes */
	float *theta;      /* length=n_continuous*n_classes */
	float *class_pdf;  /* length=n_classes */
	float *conditional_nominal_pdf; /* length=countTableLength-n_classes */
	int *i_pdf;        /* length=n_continuous */
	int *n_categories; /* length=n_nominal    */
	int n_classes;
	int n_nominal;
	int n_continuous;
	float dlambda;
	float clambda;
	float zero_correction;
} Imsls_f_nb_classifier;

typedef struct
{
	double *avg;
	double *s;
	double *logMeans;
	double *logStdev;
	double *a;
	double *b;
	double *theta;
	double *class_pdf;
	double *conditional_nominal_pdf;
	int *i_pdf;
	int *n_categories;
	int n_classes;
	int n_nominal;
	int n_continuous;
	double dlambda;
	double clambda;
	double zero_correction;
}Imsls_d_nb_classifier;

/* data structures for genetic algorithms */
typedef struct
{
	int binaryIndex;
	int nominalIndex;
	int integerIndex;
	int realIndex;
	int c_length;
	int total_length;
	int n_binary;
	int n_nominal;
	int n_integer;
	int n_intBits;
	int n_real;
	int n_realBits;
	int *n_categories;
	int *i_intervals;
	int *i_bits;
	int *i_bounds;
	int *r_intervals;
	int *r_bits;
	int *allele;
	float *r_bounds;
}Imsls_f_chromosome;

typedef struct
{
	int binaryIndex;
	int nominalIndex;
	int integerIndex;
	int realIndex;
	int c_length;
	int total_length;
	int n_binary;
	int n_nominal;
	int n_integer;
	int n_intBits;
	int n_real;
	int n_realBits;
	int *n_categories;
	int *i_intervals;
	int *i_bits;
	int *i_bounds;
	int *r_intervals;
	int *r_bits;
	int *allele;
	double *r_bounds;
}Imsls_d_chromosome;
typedef struct
{
	int n_parents;
	int encoding;
	int total_length;
	Imsls_f_chromosome* chromosome;
	int *parent;
	int *nominalPhenotype;
	int *binaryPhenotype;
	int *integerPhenotype;
	float *realPhenotype;
}Imsls_f_individual;
typedef struct
{
	int n_parents;
	int encoding;
	int total_length;
	Imsls_d_chromosome* chromosome;
	int *parent;
	int *nominalPhenotype;
	int *binaryPhenotype;
	int *integerPhenotype;
	double *realPhenotype;
}Imsls_d_individual;
typedef struct
{
	int n;
	int indexFittest;
	int indexWeakest;
	float avgFitness;
	float stdFitness;
	float maxFitness;
	float minFitness;
	float* fitness;
	Imsls_f_chromosome* chromosome;
	Imsls_f_individual** individual;
}Imsls_f_population;

typedef struct
{
	int n;
	int indexFittest;
	int indexWeakest;
	double avgFitness;
	double stdFitness;
	double maxFitness;
	double minFitness;
	double* fitness;
	Imsls_d_chromosome* chromosome;
	Imsls_d_individual** individual;
}Imsls_d_population;

typedef struct {
    float                             x_minimum;
    float                             x_maximum;
    float                             y_minimum;
    float                             y_maximum;
    float                             dfe;
    float                             sse;
    float                             smultc;
    float                             saddc;
    float                             *a;
    float                             *b;
    float                             *scoef;
    float                             *d;
    int                               output_degree;
} Imsls_f_poly_regression;
typedef struct {
    double                            x_minimum;
    double                            x_maximum;
    double                            y_minimum;
    double                            y_maximum;
    double                            dfe;
    double                            sse;
    double                            smultc;
    double                            saddc;
    double                            *a;
    double                            *b;
    double                            *scoef;
    double                            *d;
    int                               output_degree;
} Imsls_d_poly_regression;

/* Monte Carlo Functions */
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
} Imsls_faure;


/* enums for imsls_f_ts_arma  */
typedef enum {
    IMSLS_NO_ESTIMATION               = 0,
    IMSLS_MOMENTS                     = 1,
    IMSLS_DURBIN_LEVINSON             = 2,
    IMSLS_INNOVATIONS                 = 3,
    IMSLS_PRELIMINARY_ARMA            = 4,
    IMSLS_LSE                         = 5,
    IMSLS_MLE                         = 6
} Imsls_est_method;

typedef enum {
    IMSLS_NO_SE_CCF                   = 0,
    IMSLS_BARTLETT_SE_CCF             = 1,
    IMSLS_XY_NO_CORR_SE_CCF           = 2,
    IMSLS_X_NOISES_SE_CCF             = 3
} Imsls_ccf_se_option;
 
/* enums for permute functions */
typedef enum {
    IMSLS_FORWARD_PERMUTATION         = 1,
    IMSLS_BACKWARD_PERMUTATION        = 2,
    IMSLS_PERMUTE_ROWS                = 3,
    IMSLS_PERMUTE_COLUMNS             = 4
} Imsls_permute;

typedef enum {
    IMSLS_ALL                         = 1,
    IMSLS_LEAVE_OUT_LAST              = 2,
    IMSLS_SUM_TO_ZERO                 = 3
} Imsls_dummy_method;

typedef enum {
    IMSLS_NOTE                        = 1,
    IMSLS_ALERT                       = 2,
    IMSLS_WARNING                     = 3,
    IMSLS_FATAL                       = 4,
    IMSLS_TERMINAL                    = 5,
    IMSLS_WARNING_IMMEDIATE           = 6,
    IMSLS_FATAL_IMMEDIATE             = 7
} Imsls_error;

/* enums for imsls_page */
typedef enum {
    IMSLS_SET_PAGE_WIDTH              =-1,
    IMSLS_GET_PAGE_WIDTH              = 1,
    IMSLS_SET_PAGE_LENGTH             =-2,
    IMSLS_GET_PAGE_LENGTH             = 2
} Imsls_page_options;
/* enums for imsls_write_options */
typedef enum {
    IMSLS_SET_DEFAULTS                = 0,
    IMSLS_SET_CENTERING               =-1,
    IMSLS_GET_CENTERING               = 1,
    IMSLS_SET_ROW_WRAP                =-2,
    IMSLS_GET_ROW_WRAP                = 2,
    IMSLS_SET_PAGING                  =-3,
    IMSLS_GET_PAGING                  = 3,
    IMSLS_SET_NAN_CHAR                =-4,
    IMSLS_GET_NAN_CHAR                = 4,
    IMSLS_SET_TITLE_PAGE              =-5,
    IMSLS_GET_TITLE_PAGE              = 5,
    IMSLS_SET_FORMAT                  =-6,
    IMSLS_GET_FORMAT                  = 6
} Imsls_write_options;

/* enums for quadrature routines */
typedef enum {
    IMSLS_ALG                         = 1,
    IMSLS_ALG_LEFT_LOG                = 2,
    IMSLS_ALG_RIGHT_LOG               = 3,
    IMSLS_ALG_LOG                     = 4,
    IMSLS_INF_BOUND                   = 5,
    IMSLS_BOUND_INF                   = 6,
    IMSLS_INF_INF                     = 7,
    IMSLS_COS                         = 8,
    IMSLS_SIN                         = 9
} Imsls_quad;

/* enums for Neural Network routines */
typedef enum {
    IMSLS_LINEAR                      = 0,
    IMSLS_LOGISTIC                    = 1,
    IMSLS_TANH                        = 2,
    IMSLS_SQUASH                      = 3,
    IMSLS_LOGISTIC_TABLE_LOOKUP       = 4
} Imsls_mlff_activation_fcn; 

/* enums for Neural Network Weights Init */
typedef enum {
    IMSLS_NN_NETWORK                  = 0,
    IMSLS_EQUAL                       = 1,
    IMSLS_RANDOM                      = 2, 
    IMSLS_PRINCIPAL_COMPONENTS        = 3,
    IMSLS_DISCRIMINANT                = 4
} Imsls_mlff_weight_method;
/* enums for Naive Bayes routines */
typedef enum {
    IMSLS_GAUSSIAN                    = 0,
    IMSLS_LOG_NORMAL                  = 1,
    IMSLS_GAMMA                       = 2,
    IMSLS_POISSON                     = 3,
    IMSLS_CHI_SQUARE                  = 4,
    IMSLS_USER                        = 5
} Imsls_continuous_distribution; 

/* enums for genetric algorithms   */
typedef enum{
      IMSLS_ROULETTE_WITH           = 0,
      IMSLS_ROULETTE_WITHOUT        = 1,
      IMSLS_DETERMINISTIC           = 2,
      IMSLS_REMAINDER_WITH          = 3,
      IMSLS_REMAINDER_WITHOUT       = 4,
      IMSLS_SUS_SELECTION           = 5,
      IMSLS_RANK_SELECTION          = 6,
      IMSLS_TOURNAMENT_1            = 7,
      IMSLS_TOURNAMENT_2            = 8
}Imsls_selection_model;
typedef enum{
      IMSLS_NONE                    = 0,
      IMSLS_FINAL                   = 1,
      IMSLS_TRACE_GEN               = 2,
      IMSLS_TRACE_ALL               = 3,
      IMSLS_DATA_WARNINGS           = 4
}Imsls_print_level;

/* Datamining */
typedef struct
{
  int          layer_id;
  int          n_inLinks;
  int          n_outLinks;
  int          *inLinks;   /* index to Links array */
  int          *outLinks;  /* index to Links array */
  float        gradient;
  float        bias;
  int          ActivationFcn;  /* In future this should be a ptr to a fcn */
} Imsls_f_NN_Node;

typedef struct
{
  int          layer_id;
  int          n_inLinks;
  int          n_outLinks;
  int          *inLinks;   /* index to Links array */
  int          *outLinks;  /* index to Links array */
  double       gradient;
  double       bias;
  int          ActivationFcn;  /* In future this should be a ptr to a fcn */
} Imsls_d_NN_Node;

typedef struct
{
  int          n_nodes;
  int          *nodes;         /* An array containing the indices into the
                                  Node array that belong to this layer */
} Imsls_NN_Layer;

typedef struct
{
  float        weight;
  float        gradient;
  int          to_node;   /* index of to node */
  int          from_node;  /* index of from node */
} Imsls_f_NN_Link;

typedef struct
{
  double       weight;
  double       gradient;
  int          to_node;
  int          from_node;
} Imsls_d_NN_Link;

typedef struct
{
  int               n_inputs;
  int               n_outputs;
  int               n_layers;
  Imsls_NN_Layer    *layers;
  int               n_links;
  int               next_link;
  Imsls_f_NN_Link   *links;
  int               n_nodes;
  Imsls_f_NN_Node   *nodes;
} Imsls_f_NN_Network;

typedef struct
{
  int               n_inputs;
  int               n_outputs;
  int               n_layers;
  Imsls_NN_Layer    *layers;
  int               n_links;
  int               next_link;
  Imsls_d_NN_Link   *links;
  int               n_nodes;
  Imsls_d_NN_Node   *nodes;
} Imsls_d_NN_Network;


#if (defined(__alpha) && defined(__osf__)) || defined(__64BIT__) || defined(__ia64) || defined(__x86_64) || defined(_FLOAT0)
typedef void    IMSLS_PROTO((IMSLS_DECL *Imsls_error_print_proc),
                            (Imsls_error, int, char*, char*));
#elif defined(COMPUTER_LOPT64) || defined(COMPUTER_INTEL64) || defined(COMPUTER_MACOSX) || defined(COMPUTER_LOPT64_PGI) || defined(COMPUTER_MACX64)
typedef void    IMSLS_PROTO((IMSLS_DECL *Imsls_error_print_proc),
                            (Imsls_error, int, char*, char*));
#else
typedef void    IMSLS_PROTO((IMSLS_DECL *Imsls_error_print_proc),
                            (Imsls_error, long, char*, char*));
#endif

IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_neg,
                        (f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_neg,
                        (d_complex ));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_add,
                        (f_complex, f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_add,
                        (d_complex, d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_sub,
                        (f_complex, f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_sub,
                        (d_complex, d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_mul,
                        (f_complex, f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_mul,
                        (d_complex, d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_div,
                        (f_complex, f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_div,
                        (d_complex, d_complex));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_eq,
                        (f_complex, f_complex));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_eq,
                        (d_complex, d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_cz_convert,
                        (d_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_zc_convert,
                        (f_complex));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_aimag,
                        (f_complex));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_aimag,
                        (d_complex));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_fc_convert,
                        (f_complex));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_dz_convert,
                        (d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_cf_convert,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_zd_convert,
                        (double, double));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_conjg,
                        (f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_conjg,
                        (d_complex));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_arg,
                        (f_complex));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_arg,
                        (d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_sqrt,
                        (f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_sqrt,
                        (d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_log,
                        (f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_log,
                        (d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_exp,
                        (f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_exp,
                        (d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_sin,
                        (f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_sin,
                        (d_complex));
IMSLS_CI f_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_cos,
                        (f_complex));
IMSLS_CI d_complex       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_cos,
                        (d_complex));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_c_abs,
                        (f_complex));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_z_abs,
                        (d_complex));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_cholesky,
                        (int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_cholesky,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_table_twoway,
                        (int, float[], float[], int, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_table_twoway,
                        (int, double[], double[], int, int, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_hypothesis_partial,
                        (Imsls_f_regression*, int, float*, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_hypothesis_partial,
                        (Imsls_d_regression*, int, double*, ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_hypothesis_test,
                        (Imsls_f_regression*, IMSLS_F, float*, ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_hypothesis_test,
                        (Imsls_d_regression*, double, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_hypothesis_scph,
                        (Imsls_f_regression*, int, float*, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_hypothesis_scph,
                        (Imsls_d_regression*, int, double*, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_robust_covariances,
                        (int, int, float*, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_robust_covariances,
                        (int, int, double*, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_pooled_covariances,
                        (int, int, float*, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_pooled_covariances,
                        (int, int, double*, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_partial_covariances,
                        (int, int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_partial_covariances,
                        (int, int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_survival_estimates,
                        (Imsls_f_survival*, int, float*, IMSLS_F, int,
                         IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_survival_estimates,
                        (Imsls_d_survival*, int, double*, double, int,
                         double, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_survival_glm,
                        (int, int, int, int, float*, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_survival_glm,
                        (int, int, int, int, double*, ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_exact_network,
                        (int, int, float[], ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_exact_network,
                        (int, int, double[], ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_exact_enumeration,
                        (int, int, float[], ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_exact_enumeration,
                        (int, int, double[], ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_cochran_q_test,
                        (int, int, float*, ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_cochran_q_test,
                        (int, int, double*, ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_move_nan,
                        (int, int, float*, ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_move_nan,
                        (int, int, double*, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_test_nan,
                        (IMSLS_F));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_test_nan,
                        (double));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_factorial,
                        (int));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_factorial,
                        (int));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_binomial_coefficient,
                        (int, int));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_binomial_coefficient,
                        (int, int));
IMSLS_CI f_complex      *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_zeros_poly,
                        (int, float*, ...));
IMSLS_CI d_complex      *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_zeros_poly,
                        (int, double*, ...));
IMSLS_CI f_complex      *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_eig_gen,
                        (int, float*, ...));
IMSLS_CI d_complex      *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_eig_gen,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_nonlin_least_squares,
                        (void(IMSLS_DECL *)(int, int, float[], float[]),
                         int, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_nonlin_least_squares,
                        (void(IMSLS_DECL *)(int, int, double[], double[]),
                         int, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_min_con_nonlin,
                        (void(IMSLS_DECL *)(int, int, int, float[],
                         int[], float*, float[]),
                         int, int, int, int, float[], float[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_min_con_nonlin,
                        (void(IMSLS_DECL *)(int, int, int, double[],
                         int[], double*, double[]),
                         int, int, int, int, double[], double[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_min_uncon_multivar,
                        (IMSLS_F(IMSLS_DECL *)(int, float[]), int, ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_min_uncon,
                        (IMSLS_F(IMSLS_DECL *)(float), float a, float b, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_min_uncon,
                        (double(IMSLS_DECL *)(double), double a, double b, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_min_uncon_multivar,
                        (double(IMSLS_DECL *)(int, double[]), int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_zeros_sys_eqn,
                        (void(IMSLS_DECL *)(int, float[], float[]), int n, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_zeros_sys_eqn,
                        (void(IMSLS_DECL *)(int, double[], double[]), int n, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_nonlinear_optimization,
                        (IMSLS_F(IMSLS_DECL *)(int, float[], int, float[]),
                         int, int, int, float*, float[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_nonlinear_optimization,
                        (double(IMSLS_DECL *)(int, double[], int, double[]),
                         int, int, int, double*, double[], ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_discriminant_analysis,
                        (int, int, float*, int, ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_discriminant_analysis,
                        (int, int, double*, int, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_logarithmic,
                        (int, IMSLS_F, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_logarithmic,
                        (int, double, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_geometric,
                        (int, IMSLS_F, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_geometric,
                        (int, double, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_hypergeometric,
                        (int, int, int, int, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_hypergeometric,
                        (int, int, int, int, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_neg_binomial,
                        (int, IMSLS_F, IMSLS_F, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_neg_binomial,
                        (int, double, double, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_binomial,
                        (int, int, IMSLS_F, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_binomial,
                        (int, int, double, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_uniform_discrete,
                        (int, int, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_uniform_discrete,
                        (int, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_chi_squared,
                        (int, IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_chi_squared,
                        (int, double, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_lognormal,
                        (int, IMSLS_F, IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_lognormal,
                        (int, double, double, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_student_t,
                        (int, IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_student_t,
                        (int, double, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_triangular,
                        (int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_triangular,
                        (int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_von_mises,
                        (int, IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_von_mises,
                        (int, double, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_weibull,
                        (int, IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_weibull,
                        (int, double, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_cauchy,
                        (int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_cauchy,
                        (int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_aggregate,
                        (int, float[], int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_aggregate,
                        (int, double[], int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_arma,
                        (int, int, float[], int, float[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_arma,
                        (int, int, double[], int, double[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_box_cox_transform,
                        (int, float[], IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_box_cox_transform,
                        (int, double[], double, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_categorical_glm,
                        (int, int, int, int, float*, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_categorical_glm,
                        (int, int, int, int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_table_oneway,
                        (int, float[], int,...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_table_oneway,
                        (int, double[], int,...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_sort_data,
                        (int, int, float*, int, ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_sort_data,
                        (int, int, double*, int, ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_wilcoxon_rank_sum,
                        (int, float[], int, float[], ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_wilcoxon_rank_sum,
                        (int, double[], int, double[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_exponential,
                        (int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_exponential,
                        (int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_exponential_mix,
                        (int, IMSLS_F, IMSLS_F, IMSLS_F, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_exponential_mix,
                        (int, double, double, double, ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_beta_cdf,
                        (IMSLS_F, IMSLS_F, IMSLS_F));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_beta_cdf,
                        (double, double, double));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_bivariate_normal_cdf,
                        (IMSLS_F, IMSLS_F, IMSLS_F));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_bivariate_normal_cdf,
                        (double, double, double));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_beta_inverse_cdf,
                        (IMSLS_F, IMSLS_F, IMSLS_F));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_beta_inverse_cdf,
                        (double, double, double));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_anova_oneway,
                        (int, int*, float*, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_anova_oneway,
                        (int, int*, double*, ...));
IMSLS_CI float			*IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ancovar,
						(int ngroup, int ncov, int ni[], float y[], float x[], ...));
IMSLS_CI double			*IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ancovar,
						(int ngroup, int ncov, int ni[], double y[], double x[], ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_anova_factorial,
                        (int, int*, float*, ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_anova_factorial,
                        (int, int*, double*, ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_contingency_table,
                        (int, int, float*, ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_contingency_table,
                        (int, int, double*, ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_sign_test,
                        (int, float*, ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_sign_test,
                        (int, double*, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_multiple_comparisons,
                        (int, float[], int, IMSLS_F, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_multiple_comparisons,
                        (int, double[], int, double, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_poly_regression,
                        (int, float[], float[], int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_poly_regression,
                        (int, double[], double[], int, ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_regression_summary, 
                        (Imsls_f_regression*, ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_regression_summary, 
                        (Imsls_d_regression*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_regression_prediction,
                        (Imsls_f_regression*, int, float[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_regression_prediction,
                        (Imsls_d_regression*, int, double[], ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_regression_selection,
                        (int, int, float[], float[], ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_regression_selection,
                        (int, int, double[], double[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_data_sets,
                        (int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_data_sets,
                        (int, ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_regression_stepwise,
                        (int, int, float[], float[], ...));
IMSLS_CI void            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_regression_stepwise,
                        (int, int, double[], double[], ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_cluster_k_means,
                        (int, int, float[], int, float*, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_cluster_k_means,
                        (int, int, double[], int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_principal_components,
                        (int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_principal_components,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_factor_analysis, 
                        (int, float*, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_factor_analysis, 
                        (int, double*, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_nonlinear_regression,
                        (IMSLS_F(IMSLS_DECL *)(int, float[], int, float[]),
                         int, int, int, float[], float[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_nonlinear_regression,
                        (double(IMSLS_DECL *)(int, double[], int, double[]),
                         int, int, int, double[], double[], ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_regressors_for_glm,
                        (int, float[], int, int, ...));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_regressors_for_glm,
                        (int, double[], int, int, ...));
/*  No longer documented starting with CNL 7.0 */
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_normality_test,
                        (int, float[], ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_normality_test,
                        (int, double[], ...));
/*   Added for 7.0 normality test replaced with the following */
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_chi_squared_normality_test,
                        (int, int, float[], ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_chi_squared_normality_test,
                        (int, int, double[], ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_lilliefors_normality_test,
                        (int, float[], ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_lilliefors_normality_test,
                        (int, double[], ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_shapiro_wilk_normality_test,
                        (int, float[], ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_shapiro_wilk_normality_test,
                        (int, double[], ...));
/* ***  */
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_arma,
                        (int, float[], int, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_arma,
                        (int, double[], int, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_poly_prediction, 
                        (Imsls_f_poly_regression*, int, float[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_poly_prediction, 
                        (Imsls_d_poly_regression*, int, double[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_difference, 
                        (int, float[], int, int[], ...)); 
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_difference,
                        (int, double[], int, int[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_arma_forecast, 
                        (Imsls_f_arma*, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_arma_forecast, 
                        (Imsls_d_arma*, int, ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_normal_one_sample,
                        (int, float[], ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_normal_one_sample,
                        (int, double[], ...));
IMSLS_CI IMSLS_F         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_normal_two_sample,
                        (int, float[], int, float[], ...));
IMSLS_CI double          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_normal_two_sample,
                        (int, double[], int, double[], ...));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_i_min,
                        (int, int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_min,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_min,
                        (double, double));
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_i_vmin,
                        (int, ...));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_vmin,
                        (int, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_vmin,
                        (int, ...));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_i_max,
                        (int, int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_max,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_max,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_vmax,
                        (int, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_vmax,
                        (int, ...));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_ii_power,
                        (int, int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_fi_power,
                        (IMSLS_F, int));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_di_power,
                        (double, int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_ff_power,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_dd_power,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_erf,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_erf,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_erfc,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_erfc,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_erf_inverse,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_erf_inverse,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_erfc_inverse,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_erfc_inverse,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_gamma,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_gamma,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_log_gamma,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_log_gamma,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_gamma_incomplete,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_gamma_incomplete,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_t_inverse_cdf,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_t_inverse_cdf,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_normal_inverse_cdf,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_normal_inverse_cdf,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_binomial_cdf,
                        (int, int, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_binomial_cdf,
                        (int, int, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_normal_cdf,
                        (IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_normal_cdf,
                        (double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_chi_squared_cdf,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_chi_squared_cdf,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_complementary_chi_squared_cdf,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_complementary_chi_squared_cdf,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_chi_squared_inverse_cdf,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_chi_squared_inverse_cdf,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_F_cdf,
                        (IMSLS_F, IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_F_cdf,
                        (double, double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_complementary_F_cdf,
                        (IMSLS_F, IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_complementary_F_cdf,
                        (double, double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_gamma_cdf,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_gamma_cdf,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_complementary_t_cdf,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_complementary_t_cdf,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_t_cdf,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_t_cdf,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_F_inverse_cdf,
                        (IMSLS_F, IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_F_inverse_cdf,
                        (double, double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_hypergeometric_cdf,
                        (int, int, int, int));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_hypergeometric_cdf,
                        (int, int, int, int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_poisson_cdf,
                        (int, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_poisson_cdf,
                        (int, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_beta,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_beta,
                        (double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_beta_incomplete,
                        (IMSLS_F, IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_beta_incomplete,
                        (double, double, double));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_log_beta,
                        (IMSLS_F, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_log_beta,
                        (double, double));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_i_write_matrix,
                        (char*, int, int, int[], ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_write_matrix,
                        (char*, int, int, float[], ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_write_matrix,
                        (char*, int, int, double[], ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_page,
                        (int, int*));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_write_options,
                        (int, int*));
IMSLS_CI float         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_regression,
                        (int, int, float*, float[], ...));
IMSLS_CI double        *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_regression,
                        (int, int, double*, double[], ...));
IMSLS_CI float         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ranks,
                        (int, float[], ...));
IMSLS_CI double        *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ranks,
                        (int, double[], ...));
IMSLS_CI float         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_simple_statistics,
                        (int, int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_simple_statistics,
                        (int, int, double*, ...));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_chi_squared_test,
                        (IMSLS_F(IMSLS_DECL *)(IMSLS_F), int, int,
                         float*, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_chi_squared_test,
                        (double(IMSLS_DECL *)(double), int, int,
                         double*, ...));
IMSLS_CI float         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_covariances,
                        (int, int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_covariances,
                        (int, int, double*, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_seed_set,
                        (int));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_seed_get,
                        (void));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_option,
                        (int));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_uniform,
                        (int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_uniform,
                        (int, ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_poisson,
                        (int, IMSLS_F, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_normal,
                        (int, ...));
IMSLS_CI double          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_normal,
                        (int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_normal_multivariate,
                        (int, int, float*, ...));
IMSLS_CI double          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_normal_multivariate,
                        (int, int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_gamma,
                        (int, IMSLS_F, ...));
IMSLS_CI double          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_gamma,
                        (int, double, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_beta,
                        (int, IMSLS_F, IMSLS_F, ...));
IMSLS_CI double          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_beta,
                        (int, double, double, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_output_file,
                        (int, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_ctime,
                        (void));
#if (defined(__alpha) && defined(__osf__)) || defined(__64BIT__) || defined(__ia64) || defined(__x86_64) || defined(_FLOAT0)
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_error_code,
                        (void));
#elif defined(COMPUTER_SG64XS) || defined(COMPUTER_FUJITSU) || defined(COMPUTER_LINUX64) || defined(COMPUTER_SL64XS) || defined(COMPUTER_HP64XS) || defined(COMPUTER_HPSI64)   || defined(COMPUTER_FUJSOL64) || defined(COMPUTER_LOPT64) || defined(COMPUTER_INTEL64)  || defined(COMPUTER_MACOSX) || defined(COMPUTER_LOPT64_PGI) || defined(COMPUTER_MACX64)
IMSLS_CI int             IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_error_code,
                        (void));
#else
IMSLS_CI long           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_error_code,
                        (void));
#endif
IMSLS_CI Imsls_error    IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_error_type,
                         (void));
IMSLS_CI char         * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_error_message,
                        (void));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_error_options,
                        (int, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_omp_options,
                        (int, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_free,
                        (void *));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_free_allocated_memory,
                        (void *));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_skip_signal_handler,
                        ());
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_i_machine,
                        (int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_machine,
                        (int));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_machine,
                        (int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_constant,
                        (char*, char*));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_constant,
                        (char*, char*));
IMSLS_CI char           *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_version,
                        (int));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_vector_norm,
                        (int, float*, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_vector_norm,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mat_mul_rect,
                        (char*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mat_mul_rect,
                        (char*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_permute_vector,
                        (int, float*, int*, Imsls_permute, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_permute_vector,
                        (int, double*, int*, Imsls_permute, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_permute_matrix,
                        (int, int, float*, int*, Imsls_permute, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_permute_matrix,
                        (int, int, double*, int*, Imsls_permute, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_autocorrelation,
                        (int, float*, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_autocorrelation,
                        (int, double*, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_partial_autocorrelation,
                        (int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_partial_autocorrelation,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_lack_of_fit,
                        (int, float[], int, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_lack_of_fit,
                        (int, double[], int, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_Lnorm_regression,
                        (int, int,float*, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_Lnorm_regression,
                        (int, int,double*, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_kruskal_wallis_test,
                        (int, int*, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_kruskal_wallis_test,
                        (int, int*, double*, ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_friedmans_test,
                        (int, int, float*, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_friedmans_test,
                        (int, int, double*, ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_anova_nested,
                        (int, int, int*, float*, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_anova_nested,
                        (int, int, int*, double*, ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_randomness_test,
                        (int, float*, int, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_randomness_test,
                        (int, double*, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_multivar_normality_test,
                        (int, int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_multivar_normality_test,
                        (int, int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_kolmogorov_two,
                        (int, float*, int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_kolmogorov_two,
                        (int, double*, int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_kolmogorov_one,
                        (float (*fcn) (float), int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_kolmogorov_one,
                        (double (*fcn) (double), int, double*, ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_chi_sq_inv,
                        (float, float, float));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_chi_sq_inv,
                        (double, double, double));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_chi_sq,
                        (float, float, float));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_chi_sq,
                        (double, double, double));
/* Added in 7.0 */
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_chi_sq_pdf,
                        (float x, float df, float lambda));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_chi_sq_pdf,
                        (double x, double df, double lambda));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_t_pdf,
                        (float x, float df, float delta));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_t_pdf,
                        (double x, double df, double delta));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_F_pdf,
                        (float x, float df1, float df2, float lambda));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_F_pdf,
                        (double x, double df1, double df2, double lambda));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_F_cdf,
                        (float x, float df1, float df2, float lambda));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_F_cdf,
                        (double x, double df1, double df2, double lambda));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_F_inverse_cdf,
                        (float x, float df1, float df2, float lambda));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_F_inverse_cdf,
                        (double x, double df1, double df2, double lambda));

IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_t_inv_cdf,
                        (float, int, float));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_t_inv_cdf,
                        (double, int, double));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_non_central_t_cdf,
                        (float, int, float));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_non_central_t_cdf,
                        (double, int, double));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_wilcoxon_sign_rank,
                        (int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_wilcoxon_sign_rank,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_k_trends_test,
                        (int, int*, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_k_trends_test,
                        (int, int*, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_noether_cyclical_trend,
                        (int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_noether_cyclical_trend,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_tie_statistics,
                        (int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_tie_statistics,
                        (int, double*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_cox_stuart_trends_test,
                        (int, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_cox_stuart_trends_test,
                        (int, double*, ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_binomial_pdf,
                        (int, int, float));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_binomial_pdf,
                        (int, int, double));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_anova_balanced,
                        (int, int*, float*, int, int*, int, int*, int*, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_anova_balanced, 
                        (int, int*, double*, int, int*, int, int*, int*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_garch,
                        (int, int, int, float*, float*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_garch, 
                        (int, int, int, double*, double*, ...));
IMSLS_CI  void          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_lfcn,
                        (int*, int*, int*, float x[],
                        float y[], float *a, float sigma[]));
IMSLS_CI  void          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_lfcn,
                        (int*, int*, int*, double x[],
                        double y[], double *a, double sigma[])); 
IMSLS_CI int*           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_general_discrete, 
                        (int n_random, int imin, int nmass, float *probs, ...));
IMSLS_CI int*           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_general_discrete, 
                        (int n_random, int imin, int nmass, double *probs, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_drngda, 
                        (int *nr, int *iopt, int *imin,
                        int *nmass, double probs[], int iwk[], 
                        double wk[], int ir[]));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_discrete_table_setup, 
                        (float (*prf) (int), float del,
                        int nndx, int *imin, int *nmass, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_discrete_table_setup, 
                        (double (*prf) (int), double del,
                        int nndx, int *imin, int *nmass, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_rngdt, 
                        (int *nr, int *imin, int *nmass,
                        float cumpr[], int ir[]));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_drngdt, 
                        (int *nr, int *imin, int *nmass,
                        double cumpr[], int ir[]));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_continuous_table_setup,
                        (float (*cdf) (float), int iopt, int ndata, 
                        float *table, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_continuous_table_setup,
                        (double (*cdf) (double), int iopt, int ndata, 
                        double *table, ...));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_general_continuous,
                        (int n_random, int ndata, float *table, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_general_continuous,
                        (int n_random, int ndata, double *table, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_drngct, 
                        (int *nr, int *ndata, double *table,
                        int *ldtabl, double r[]));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_stable, 
                        (int nr, float alpha, float bprime, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_stable, 
                        (int nr, double alpha, double bprime, ...));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_orthogonal_matrix, 
                        (int n, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_orthogonal_matrix, 
                        (int n, ...));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_mvar_from_data,
                        (int n_random, int ndim, int nsamp,
                        float *x, int nn,  ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_mvar_from_data,
                        (int n_random, int ndim, int nsamp,
                        double *x, int nn,  ...));
IMSLS_CI int*           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_multinomial, 
                        (int nr, int n, int k, float p[], ...));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_npp, 
                        (float, float, float (*ftheta) (float),
                        float, float, int, int*, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_npp, 
                        (double, double, double (*ftheta) (double),
                        double, double, int, int*, ...));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_option_get, 
                        ());
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_substream_seed_get, 
                        (int iseed1));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_table_get, 
                        (float **table, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_table_get, 
                        (double **table, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_table_set,
                        (float[]));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_table_set,
                        (double[]));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_GFSR_table_get, 
                        (int **table, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_GFSR_table_set,
                        (int[]));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_order_uniform,
                        (int ifirst, int ilast, int n, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_order_uniform,
                        (int ifirst, int ilast, int n, ...));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_order_normal,
                        (int ifirst, int ilast, int n, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_order_normal,
                        (int ifirst, int ilast, int n, ...));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_sphere,
                        (int n_random, int k, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_sphere,
                        (int n_random, int k, ...));
IMSLS_CI int*           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_permutation,
                        (int k, ...));
IMSLS_CI int*           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_sample_indices,
                        (int nsamp, int npop, ...));
IMSLS_CI int*           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_table_twoway,
                        (int, int, int nrtot[], int nctot[], ...));
IMSLS_CI float*         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_random_sample,
                        (int nrow, int nvar, float *population, int nsamp, ...));
IMSLS_CI double*        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_random_sample,
                        (int nrow, int nvar, double *population, int nsamp, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_kalman, 
                        (int nb, float b[], float *covb,
                         int *n, float *ss, float *alndet,...)); 
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_kalman, 
                        (int nb, double b[], double *covb,
                        int *n, double *ss, double *alndet,...)); 
IMSLS_CI Imsls_faure    *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_faure_sequence_init,
                        (int, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_faure_sequence_free,
                        (Imsls_faure*));
IMSLS_CI IMSLS_F        *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_faure_next_point,
                        (Imsls_faure*, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_faure_next_point,
                        (Imsls_faure*, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_crosscorrelation,
                        (int, float*, float*, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_crosscorrelation,
                        (int, double*, double*, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_multi_crosscorrelation,
                        (int, int, float*, int, int, float*, int, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_multi_crosscorrelation,
                        (int, int, double*, int, int, double*, int, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_prop_hazards_gen_lin,
                        (int nrow, int ncol, float x[], int nef, int nvef[], 
                        int indef[], int maxcl, int * ncoef, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_prop_hazards_gen_lin,
                        (int nrow, int ncol, double x[], int nef, int nvef[], 
                        int indef[], int maxcl, int * ncoef, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_life_tables,
                        (int n, float age[], float a[], int ipop[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_life_tables,
                        (int n, double age[], double a[], int ipop[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_kaplan_meier_estimates,
                        (int nobs, int ncol, float * x, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_kaplan_meier_estimates,
                        (int nobs, int ncol, double * x, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_split_plot,
                        (int n, int n_locations, int n_whole, int n_split,
                        int rep[], int whole[], int split[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_split_plot,
                        (int n, int n_locations, int n_whole, int n_split,
                        int rep[], int whole[], int split[], double y[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_split_split_plot,
                        (int n, int n_locations, int n_whole, int n_split, int n_sub,
                        int rep[], int whole[], int split[], int sub[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_split_split_plot,
                        (int n, int n_locations, int n_whole, int n_split, int n_sub,
                        int rep[], int whole[], int split[], int sub[], double y[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_strip_plot,
                        (int n, int n_locations, int n_strip_a, int n_strip_b,
                        int block[], int strip_a[], int strip_b[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_strip_plot,
                        (int n, int n_locations, int n_strip_a, int n_strip_b,
                        int block[], int strip_a[], int strip_b[], double y[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_strip_split_plot,
                        (int n, int n_locations, int n_strip_a, int n_strip_b, int n_split,
                        int rep[], int strip_a[], int strip_b[], int split[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_strip_split_plot,
                        (int n, int n_locations, int n_strip_a, int n_strip_b, int n_split,
                        int rep[], int strip_a[], int strip_b[], int split[], double y[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_crd_factorial,
                        (int n_obs, int n_locations, int n_factors, 
                        int n_levels[], int model[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_crd_factorial,
                        (int n_obs, int n_locations, int n_factors, 
                        int n_levels[], int model[], double y[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_rcbd_factorial,
                        (int n_obs, int n_locations, int n_factors, 
                        int n_levels[], int model[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_rcbd_factorial,
                        (int n_obs, int n_locations, int n_factors, 
                        int n_levels[], int model[], double y[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_homogeneity,
                        (int n, int n_treatment, int treatment[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_homogeneity,
                        (int n, int n_treatment, int treatment[], double y[], ...));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_yates,
                        (int n, int n_independent, float x[], ...));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_yates,
                        (int n, int n_independent, double x[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_latin_square,
                        (int n, int n_locations, int n_treatments, int row[], int col[], 
                        int treatment[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_latin_square,
                        (int n, int n_locations, int n_treatments, int row[], int col[], 
                        int treatment[], double y[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_lattice,
                        (int n, int n_locations, int n_reps, int n_blocks, int n_treatments,
                        int rep[], int block[], int treatment[], float y[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_lattice,
                        (int n, int n_locations, int n_reps, int n_blocks, int n_treatments,
                        int rep[], int block[], int treatment[], double y[], ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_stud_range_cdf,
                        (float q, float v, float r));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_stud_range_cdf,
                        (double q, double v, double r));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_stud_range_inverse_cdf,
                        (float p, float v, float r));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_stud_range_inverse_cdf,
                        (double p, double v, double r));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_int_fcn,
                        (IMSLS_F(IMSLS_DECL *)(IMSLS_F), IMSLS_F, IMSLS_F, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_int_fcn,
                        (double(IMSLS_DECL *)(double), double, double, ...));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_int_fcn_inf,
                        (IMSLS_F(IMSLS_DECL *)(IMSLS_F), IMSLS_F, Imsls_quad, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_int_fcn_inf,
                        (double(IMSLS_DECL *)(double), double, Imsls_quad, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_dissimilarities,
                        (int nrow, int ncol, float *x, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_dissimilarities,
                        (int nrow, int ncol, double *x, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_cluster_hierarchical,
                        (int npt, float *dist, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_cluster_hierarchical,
                        (int npt, double *dist, ...));
IMSLS_CI int           *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_cluster_number,
                        (int npt, int *iclson, int *icrson, int k, ...));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_hypergeometric_pdf,
                        (int k, int n, int m, int l));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_hypergeometric_pdf,
                        (int k, int n, int m, int l));
IMSLS_CI float          IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_gamma_inverse_cdf,
                        (float p, float a));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_gamma_inverse_cdf,
                        (double p, double a));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_nonparam_hazard_rate,
                        (int nobs, float * x, int n_hazard, 
                        float hazard_min, float hazard_increment, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_nonparam_hazard_rate,
                        (int nobs, double * x, int n_hazard, 
                        double hazard_min, double hazard_increment, ...));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_poisson_pdf,
                        (int, IMSLS_F));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_poisson_pdf,
                        (int, double));
/* V6.0 */
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_scale_filter,
                        (int nobs, float x[], int method, ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_scale_filter,
                        (int nobs, double x[], int method, ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_time_series_filter,
                        (int n_obs, int n_var, int n_lags, float x[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_time_series_filter,
                        (int n_obs, int n_var, int n_lags, double x[], ...));
IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_time_series_class_filter,
                        (int n_obs, int n_lags, int n_classes, 
                        int i_class[], float x[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_time_series_class_filter,
                        (int n_obs, int n_lags, int n_classes, 
                        int i_class[], double x[], ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_unsupervised_nominal_filter,
                        (int n_obs, int *n_classes, int x[], ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_unsupervised_ordinal_filter, 
                        (int n_obs, int x[], float z[], ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_unsupervised_ordinal_filter, 
                        (int n_obs, int x[], double z[], ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_MT32_init,
                        (int key_length, uint32_t init_key[]));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_MT64_init,
                        (int key_length, uint64_t init_key[]));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_MT32_table_get, 
                        (uint32_t **table, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_MT32_table_set,
                        (uint32_t[]));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_MT64_table_get, 
                        (uint64_t **table, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_random_MT64_table_set,
                        (uint64_t[]));
IMSLS_CI Imsls_f_NN_Network *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_network_init,
                        (int, int));
IMSLS_CI Imsls_d_NN_Network *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_network_init,
                        (int, int));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_network,
                        (Imsls_f_NN_Network * ffnet, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_network, 
                        (Imsls_d_NN_Network * ffnet, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_network_free, 
                        (Imsls_f_NN_Network * ffnet));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_network_free, 
                        (Imsls_d_NN_Network * ffnet));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_network_trainer, 
                        (Imsls_f_NN_Network * ffnet,
                        int n_observations, int n_categorical,
                        int n_continuous, int categorial[],
                        float continuous[], float output[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_network_trainer, 
                        (Imsls_d_NN_Network * ffnet,
                        int n_observations, int n_categorical,
                        int n_continuous, int categorial[],
                        double continuous[], double output[], ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_network_forecast, 
                        (Imsls_f_NN_Network * ffnet,
                        int n_categorical, int n_continuous, 
                        int categorical[], float continuous[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_network_forecast, 
                        (Imsls_d_NN_Network * ffnet,
                        int n_categorical, int n_continuous, 
                        int categorical[], double continuous[], ...));
/* New Neural Net Classification Routines ******************************************************************/
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_classification_trainer, 
                        (Imsls_f_NN_Network * ffnet,
                        int n_observations, int n_categorical,
                        int n_continuous, int classification[], int categorial[],
                        float continuous[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_classification_trainer, 
                        (Imsls_d_NN_Network * ffnet,
                        int n_observations, int n_categorical,
                        int n_continuous, int classification[], int categorial[],
                        double continuous[], ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_pattern_classification, 
                        (Imsls_f_NN_Network * ffnet, int n_patterns,
                        int n_categorical, int n_continuous, 
                        int categorical[], float continuous[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_pattern_classification, 
                        (Imsls_d_NN_Network * ffnet, int n_patterns,
                        int n_categorical, int n_continuous, 
                        int categorical[], double continuous[], ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_initialize_weights, 
                        (Imsls_f_NN_Network *network,int n_patterns, int n_nominalAtt, 
			int n_n_continuous, int nominal[], float continuousAtt[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_initialize_weights, 
                        (Imsls_d_NN_Network *network, int n_patterns, int n_nominalAtt, 
                        int n_continuous, int nominal[], double continuousAtt[], ...));

/*   ***************************************************  
          For Internal Use only
     ***************************************************/ 
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_training_engine,
                        (int *, int *, int *, int *, int *, int *, float *,
                        float *, float *, float *, float *, float *,
                        float *, float *));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_training_engine,
                        (int *, int *, int *, int *, int *, int *, double *,
                        double *, double *, double *, double *,
                        double *, double *, double *));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_forecasting_engine, 
                        (int, int, int, int, int, int, int *, int *,
                        int *, int *, float *, float *, float *));
IMSLS_CI int            IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_forecasting_engine, 
                        (int, int, int, int, int, int, int *, int *,
                        int *, int *, double *, double *, double *));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_network_write, 
                        (Imsls_f_NN_Network* network, char* filename,  ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_network_write, 
                        (Imsls_d_NN_Network* network, char* filename,  ...));
IMSLS_CI Imsls_f_NN_Network *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_mlff_network_read, 
                        (char* filename, ...));
IMSLS_CI Imsls_d_NN_Network *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_mlff_network_read, 
                        (char* filename, ...));
/* ************************************************** */
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_f_max_arma,
                        (int nobservations, float w[], int p, int q, ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_d_max_arma,
                        (int nobservations, double w[], int p, int q, ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_f_auto_uni_ar,
                        (int n_observations, float z[], int maxlag, int *p,...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_d_auto_uni_ar,
                        (int n_observations, double z[], int maxlag, int *p,...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_f_seasonal_fit,
                        (int n_observations, float z[], int maxlag,
                            int n_differences, int n_s_initial,
                            int s_initial[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_d_seasonal_fit,
                        (int n_observations, double z[], int maxlag,
                            int n_differences, int n_s_initial,
                            int s_initial[], ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_f_estimate_missing,
                        (int n, int tpoints[], float x[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_d_estimate_missing,
                        (int n, int tpoints[], double x[], ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_f_ts_outlier_identification,
                        (int nobs, int model[], float w[], ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_d_ts_outlier_identification,
                        (int nobs, int model[], double w[], ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_f_ts_outlier_forecast,
                        (int n_obs, float series[], int num_outliers, int outlier_stat[],
                             float omega[], float delta, int model[], float parameters[],
                             int n_predict, ...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_d_ts_outlier_forecast,
                        (int n_obs, double series[], int num_outliers, int outlier_stat[],
                             double omega[], double delta, int model[], double parameters[],
                             int n_predict, ...));
IMSLS_CI float *        IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_f_auto_arima,
                        (int n_obs, int tpoints[], float x[],...));
IMSLS_CI double *       IMSLS_DECL IMSLS_EF IMSLS_PROTO (imsls_d_auto_arima,
                        (int n_obs, int tpoints[], double x[],...));
IMSLS_CI Imsls_f_ppoly *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ppoly_create,
                        (int, int, int*, int*, ...));
IMSLS_CI Imsls_d_ppoly *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ppoly_create,
                        (int, int, int*, int*, ...));
IMSLS_CI Imsls_f_ppoly *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_cub_spline_interp_e_cnd,
                        (int, float[], float[], ...));
IMSLS_CI Imsls_d_ppoly *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_cub_spline_interp_e_cnd,
                        (int, double[], double[], ...));
IMSLS_CI IMSLS_F        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_cub_spline_value,
                        (IMSLS_F, Imsls_f_ppoly*, ...));
IMSLS_CI double         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_cub_spline_value,
                        (double, Imsls_d_ppoly*, ...));


  /* end V6.0. */

/*   V6.0.1.1 */
/* Multivariate Normal CDF ******************************************************************************/

IMSLS_CI float         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_multivariate_normal_cdf,
                        (int k, float bounds[], float means[], float sigma[],  ...));
IMSLS_CI double        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_multivariate_normal_cdf,
                        (int k, double bounds[], double means[], double sigma[], ...));
IMSLS_CI float         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_int_fcn_qmc2,
                        (IMSLS_F(IMSLS_DECL *)(int, IMSLS_F*), int, IMSLS_F*, IMSLS_F*, ...));
IMSLS_CI double        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_int_fcn_qmc2,
                        (double(IMSLS_DECL *)(int, double*), int, double*, double*, ...));
IMSLS_CI float         IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_int_fcn_hyper_rect2,
                        (IMSLS_F(IMSLS_DECL *)(int, float*), int,
                         float*, float*, ...));
IMSLS_CI double        IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_int_fcn_hyper_rect2,
                        (double(IMSLS_DECL *)(int, double*), int,
                         double*, double*, ...));
IMSLS_CI float *       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_lin_sol_posdef,
                        (int, float*, float*, ...));
IMSLS_CI double *      IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_lin_sol_posdef,
                        (int, double*, double*, ...));
IMSLS_CI float *       IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_eig_sym,
                        (int, float*, ...));
IMSLS_CI double *      IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_eig_sym,
                        (int, double*, ...));

/* end V6.0.1.1 */


IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_unsupervised_discretization, 
                        (int n, float x[], ...));
IMSLS_CI int            *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_unsupervised_discretization, 
                        (int n, double x[], ...));

IMSLS_CI float          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_empirical_quantiles,
			(int n_observations, float x[], int n_qprop, float qprop[], ...));
IMSLS_CI double         *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_empirical_quantiles,
			(int n_observations, double x[], int n_qprop, double qprop[], ...));

/* New Naive Bayes Classification Routines **************************************************************/
IMSLS_CI int          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_naive_bayes_trainer,
                        (int n_patterns, int n_classes, int classification[], ...));
IMSLS_CI int          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_naive_bayes_trainer,
                        (int n_patterns, int n_classes, int classification[], ...));
IMSLS_CI int          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_naive_bayes_classification,
                        (Imsls_f_nb_classifier * nb_classifier, int n_patterns, ...));
IMSLS_CI int          *IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_naive_bayes_classification,
                        (Imsls_d_nb_classifier * nb_classifier, int n_patterns, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_nb_classifier_free,
                        (Imsls_f_nb_classifier * nb_classifier));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_nb_classifier_free,
                        (Imsls_d_nb_classifier * nb_classifier));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_nb_classifier_write,
                        (Imsls_f_nb_classifier * nb_classifier, char* filename, ...));
IMSLS_CI void           IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_nb_classifier_write,
                        (Imsls_d_nb_classifier * nb_classifier, char* filename, ...));
IMSLS_CI Imsls_f_nb_classifier * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_nb_classifier_read,
                        (char* filename, ...));
IMSLS_CI Imsls_d_nb_classifier * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_nb_classifier_read,
                        (char* filename, ...));
/* Genetic Algorithms ***************************************************************************** */

IMSLS_CI Imsls_f_chromosome * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_chromosome, (int, ...));
IMSLS_CI Imsls_d_chromosome * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_chromosome, (int, ...));
IMSLS_CI Imsls_f_individual * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_individual, (Imsls_f_chromosome *chromosome, ...));
IMSLS_CI Imsls_d_individual * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_individual, (Imsls_d_chromosome *chromosome, ...));
IMSLS_CI Imsls_f_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_random_population, (int, Imsls_f_chromosome*, ...));
IMSLS_CI Imsls_d_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_random_population,(int, Imsls_d_chromosome*, ...));
IMSLS_CI Imsls_f_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_population,(int, Imsls_f_chromosome*, Imsls_f_individual** individual, ...));
IMSLS_CI Imsls_d_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_population,(int, Imsls_d_chromosome*, Imsls_d_individual** individual, ...));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_free_population, (Imsls_f_population *population));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_free_population, (Imsls_d_population *population));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_copy_population, (Imsls_f_population *pop1, Imsls_f_population *pop2));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_copy_population, (Imsls_d_population *pop1, Imsls_d_population *pop2));
IMSLS_CI Imsls_f_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_clone_population, (Imsls_f_population *population, ...));
IMSLS_CI Imsls_d_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_clone_population, (Imsls_d_population *population, ...));
IMSLS_CI Imsls_f_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_merge_population, (Imsls_f_population *pop1, Imsls_f_population *pop2, ...));
IMSLS_CI Imsls_d_population * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_merge_population, (Imsls_d_population *pop1, Imsls_d_population *pop2, ...));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_grow_population, (int n, Imsls_f_individual **ind, Imsls_f_population *population, ...));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_grow_population, (int n, Imsls_d_individual **ind, Imsls_d_population *population, ...));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_free_individual, (Imsls_f_individual *individual));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_free_individual, (Imsls_d_individual *individual));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_copy_individual, (Imsls_f_individual *ind1, Imsls_f_individual *ind2));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_copy_individual, (Imsls_d_individual *ind1, Imsls_d_individual *ind2));
IMSLS_CI Imsls_f_individual * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_clone_individual,(Imsls_f_individual *individual, ...));
IMSLS_CI Imsls_d_individual * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_clone_individual,(Imsls_d_individual *individual, ...));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_copy_chromosome, (Imsls_f_chromosome *chrom1, Imsls_f_chromosome *chrom2));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_copy_chromosome, (Imsls_d_chromosome *chrom1, Imsls_d_chromosome *chrom2));
IMSLS_CI Imsls_f_chromosome * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_clone_chromosome,(Imsls_f_chromosome *chromosome, ...));
IMSLS_CI Imsls_d_chromosome * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_clone_chromosome,(Imsls_d_chromosome *chromosome, ...));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_decode,(Imsls_f_individual *ind));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_decode,(Imsls_d_individual *ind));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_encode,(Imsls_f_individual *ind));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_encode,(Imsls_d_individual *ind));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_ga_mutate,(float probability, Imsls_f_individual *ind, ...));
IMSLS_CI void                 IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_ga_mutate,(double probability, Imsls_d_individual *ind, ...));
IMSLS_CI Imsls_f_individual * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_f_genetic_algorithm, (float(IMSLS_DECL *) (Imsls_f_individual*),Imsls_f_population*, ...));
IMSLS_CI Imsls_d_individual * IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_d_genetic_algorithm, (double(IMSLS_DECL *)(Imsls_d_individual*),Imsls_d_population*, ...));

IMSLS_CI int                  IMSLS_DECL IMSLS_EF IMSLS_PROTO(imsls_initialize,(int reason));

#ifdef __cplusplus
}                               /* end extern "C" for C++ */
#endif

        /* Keywords */

enum Imsls_keyword {
    IMSLS_ABS_FCN_TOL                              = 10010,             
    IMSLS_ABS_FCN_TOL_ADR                          = 10020,             
    IMSLS_ACTIVE                                   = 10030,             
    IMSLS_ADD_TO_DIAG_H                            = 10040,             
    IMSLS_ADJ_R_SQUARED                            = 10050,             
    IMSLS_ALL_STEPS                                = 10060,             
    IMSLS_ALPHA                                    = 10070, 
    IMSLS_ALPHA_ADR                                = 10075,             
    IMSLS_ANOVA_TABLE                              = 10080,             
    IMSLS_ANOVA_TABLE_USER                         = 10090,             
    IMSLS_ARIMA_INFO                               = 10100,             
    IMSLS_AR_LAGS                                  = 10110,             
    IMSLS_ASCENDING                                = 10120,             
    IMSLS_AUTOCORRELATION_TEST                     = 10130,             
    IMSLS_AUTOCOVARIANCES                          = 10140,             
    IMSLS_AUXILIARY_EXACT_REPS                     = 10150,             
    IMSLS_AVERAGE_TIE                              = 10160,             
    IMSLS_A_COL_DIM                                = 10170,             
    IMSLS_A_MATRIX                                 = 10180,             
    IMSLS_BACKCASTING                              = 10190,             
    IMSLS_BACKWARD                                 = 10200,             
    IMSLS_BACKWARD_ORIGIN                          = 10210,
    /*   Added for V. 7.0 Arma_forecast */
    IMSLS_ONE_STEP_FORECAST                        = 10211,
    IMSLS_ONE_STEP_FORECAST_USER                   = 10212,
    /* */
    IMSLS_BALANCED_DATA_RESTRICTIONS               = 10220,             
    IMSLS_BARTLETT                                 = 10230,             
    IMSLS_BASIS                                    = 10240,             
    IMSLS_BASIS_ADR                                = 10250,             
    IMSLS_BLOM_SCORES                              = 10260,             
    IMSLS_BONFERRONI                               = 10270,             
    IMSLS_BOUND                                    = 10280,             
    IMSLS_BOUNDS                                   = 10290,             
    IMSLS_BOUNDS_ADR                               = 10300,             
    IMSLS_BREAKPOINTS                              = 10310,             
    IMSLS_B_COL_DIM                                = 10320,             
    IMSLS_B_MATRIX                                 = 10330,             
    IMSLS_CELL_CHI_SQUARED                         = 10340,             
    IMSLS_CELL_CHI_SQUARED_USER                    = 10350,             
    IMSLS_CELL_COUNTS                              = 10360,             
    IMSLS_CELL_COUNTS_USER                         = 10370,             
    IMSLS_CELL_EXPECTED                            = 10380,             
    IMSLS_CELL_EXPECTED_USER                       = 10390,             
    IMSLS_CELL_MEANS                               = 10400,             
    IMSLS_CELL_VARIANCES                           = 10410,             
    IMSLS_CENTRAL_COMPOSITE                        = 10420,             
    IMSLS_CHEBYSHEV_FIRST                          = 10430,             
    IMSLS_CHEBYSHEV_SECOND                         = 10440,             
    IMSLS_CHI_SQUARED                              = 10450,             
    IMSLS_CHI_SQUARED_ALTERNATIVE                  = 10460,             
    IMSLS_CONTRIBUTIONS                            = 10470,             
    IMSLS_CONTRIBUTIONS_USER                       = 10480,             
    IMSLS_CHI_SQUARED_STATS                        = 10490,             
    IMSLS_CHI_SQUARED_STATS_USER                   = 10500,             
    IMSLS_CHI_SQUARED_TEST                         = 10510,             
    IMSLS_CHI_SQUARED_TEST_NULL                    = 10520,             
    IMSLS_CI_COMMON_VARIANCE                       = 10530,              
    IMSLS_CI_DIFF_FOR_EQUAL_VARIANCES              = 10540,             
    IMSLS_CI_DIFF_FOR_UNEQUAL_VARIANCES            = 10550,             
    IMSLS_CI_MEAN                                  = 10560,             
    IMSLS_CI_NEW_SAMPLE_SUM_WEIGHTS                = 10570,             
    IMSLS_CI_RATIO_VARIANCE                        = 10580,             
    IMSLS_CI_VARIANCE                              = 10590,             
    IMSLS_CI_X0_GIVEN_NEW_SAMPLE_Y0                = 10600,             
    IMSLS_CI_X0_GIVEN_POP_MEAN_Y0                  = 10610,             
    IMSLS_CLASS                                    = 10620,             
    IMSLS_CLASS_MARKS                              = 10630,             
    IMSLS_CLUSTER_MEANS                            = 10640,             
    IMSLS_CLUSTER_SEEDS                            = 10650,             
    IMSLS_CLUSTER_SSQ                              = 10660,             
    IMSLS_CLUSTER_SSQ_USER                         = 10670,    
    IMSLS_COEFS                                    = 10680,             
    IMSLS_COEF_COVARIANCES                         = 10690,             
    IMSLS_COEF_COVARIANCES_USER                    = 10700,             
    IMSLS_COEF_T_TESTS                             = 10710,             
    IMSLS_COEF_T_TESTS_USER                        = 10720,             
    IMSLS_COEF_VIF                                 = 10730,             
    IMSLS_COEF_VIF_USER                            = 10740,             
    IMSLS_COL_LABELS                               = 10750,             
    IMSLS_COL_NUMBER                               = 10760,             
    IMSLS_COL_NUMBER_ZERO                          = 10770,             
    IMSLS_COMPILER_VERSION                         = 10780,             
    IMSLS_COMPLETELY_TESTABLE                      = 10790,             
    IMSLS_COMPLETELY_TESTABLE_USER                 = 10800,             
    IMSLS_CORRELATIONS                             = 10810,             
    IMSLS_COMPUTE_OPTION                           = 10820,             
    IMSLS_CONCAVE                                  = 10830,             
    IMSLS_CONCAVE_ITMAX                            = 10840,             
    IMSLS_CONDITION                                = 10850,             
    IMSLS_CONFIDENCE                               = 10860,             
    IMSLS_CONFIDENCE_MEAN                          = 10870,             
    IMSLS_CONFIDENCE_MEANS                         = 10880,             
    IMSLS_CONFIDENCE_MEANS_ADR                     = 10890,             
    IMSLS_CONFIDENCE_VARIANCE                      = 10900,             
    IMSLS_CONFIDENCE_VARIANCES                     = 10910,             
    IMSLS_CONFIDENCE_VARIANCES_ADR                 = 10920,             
    IMSLS_CONSTANT                                 = 10930,             
    IMSLS_CONSTR_TYPE                              = 10940,             
    IMSLS_CONTINUOUS                               = 10950,             
    IMSLS_CONTROL_POP_MEAN_Y0                      = 10960,             
    IMSLS_CONTROL_SAMPLE_MEAN_Y0                   = 10970,             
    IMSLS_CONVERGENCE_EPS                          = 10980,             
    IMSLS_CONVERGENCE_TOLERANCE                    = 10990,             
    IMSLS_COOKSD                                   = 11000,             
    IMSLS_COOKSD_USER                              = 11010,             
    IMSLS_CORRECTED_SSCP_MATRIX                    = 11020,             
    IMSLS_CORRELATION_MATRIX                       = 11030,             
    IMSLS_COSH                                     = 11040,             
    IMSLS_COVARIANCE_COL_DIM                       = 11050,
    IMSLS_COV_COL_DIM                              = 11055,             
    IMSLS_COEF_COV_COL_DIM                         = 11060,             
    IMSLS_CRITERIONS                               = 11070,             
    IMSLS_CUM_PERCENT                              = 11080,             
    IMSLS_CUTPOINTS                                = 11090,             
    IMSLS_CUTPOINTS_EQUAL                          = 11100,             
    IMSLS_CUTPOINTS_USER                           = 11110,             
    IMSLS_DATA_BOUNDS                              = 11120,             
    IMSLS_DATA_OPTION                              = 11130,             
    IMSLS_DEGREES_OF_FREEDOM                       = 11140,             
    IMSLS_DELETED_RESIDUAL                         = 11150,             
    IMSLS_DELETED_RESIDUAL_USER                    = 11160,             
    IMSLS_DERIV                                    = 11170,             
    IMSLS_DESCENDING                               = 11180,
    IMSLS_DF                                       = 11185,             
    IMSLS_DFFITS                                   = 11190,             
    IMSLS_DFFITS_USER                              = 11200,             
    IMSLS_DF_PURE_ERROR                            = 11210,             
    IMSLS_DF_S_SQUARED                             = 11220,             
    IMSLS_DIFFERENCE                               = 11230,             
    IMSLS_DIFFERENCE_ORDERS                        = 11240,             
    IMSLS_DUAL                                     = 11250,             
    IMSLS_DUAL_USER                                = 11260,             
    IMSLS_DUMMY                                    = 11270,             
    IMSLS_DUNN_SIDAK                               = 11280,             
    IMSLS_EFFECT_TESTS                             = 11290,             
    IMSLS_EIGENVALUES                              = 11300,             
    IMSLS_EIGENVALUES_USER                         = 11310,             
    IMSLS_EIGENVECTORS                             = 11320,             
    IMSLS_ELEMENTWISE                              = 11330,             
    IMSLS_ELEMENTWISE_USER                         = 11340,             
    IMSLS_EMS                                      = 11350,             
    IMSLS_EMS_USER                                 = 11360,             
    IMSLS_EPS                                      = 11370,             
    IMSLS_EPS_ADR                                  = 11380,             
    IMSLS_EQUAL_NUMBERS                            = 11390,             
    IMSLS_ERROR_MSG_NAME                           = 11400,             
    IMSLS_ERROR_MSG_PATH                           = 11410,             
    IMSLS_ERROR_PRINT_PROC                         = 11420,             
    IMSLS_ERROR_SSQ                                = 11430,             
    IMSLS_ERR_ABS                                  = 11440,             
    IMSLS_ERR_ABS_ADR                              = 11450,             
    IMSLS_ERR_EST                                  = 11460,             
    IMSLS_ERR_LIST                                 = 11470,             
    IMSLS_ERR_ORDER                                = 11480,             
    IMSLS_ERR_REL                                  = 11490,             
    IMSLS_ERR_REL_ADR                              = 11500,             
    IMSLS_ETA                                      = 11510,             
    IMSLS_ETA_ADR                                  = 11520,             
    IMSLS_EVECU_COL_DIM                            = 11530,             
    IMSLS_EXCLUDE_FIRST                            = 11540,             
    IMSLS_EXPECTED                                 = 11550,             
    IMSLS_EXPECTED_NORMAL_SCORES                   = 11560,             
    IMSLS_EXPECTED_USER                            = 11570,             
    IMSLS_FACTOR                                   = 11580,             
    IMSLS_FACTOR_ONLY                              = 11590,             
    IMSLS_FACTOR_USER                              = 11600,             
    IMSLS_FAC_COL_DIM                              = 11610,             
    IMSLS_FDATA_COL_DIM                            = 11620,             
    IMSLS_FIRST_STEP                               = 11630,             
    IMSLS_FIXED_POINT                              = 11640,             
    IMSLS_FIXED_POINT_ADR                          = 11650,             
    IMSLS_FJAC                                     = 11660,             
    IMSLS_FJAC_COL_DIM                             = 11670,             
    IMSLS_FJAC_Q                                   = 11680,             
    IMSLS_FJAC_Q_USER                              = 11690,             
    IMSLS_FJAC_R                                   = 11700,             
    IMSLS_FJAC_R_USER                              = 11710,             
    IMSLS_FJAC_USER                                = 11720,             
    IMSLS_FLOOR                                    = 11730,             
    IMSLS_FLOOR_ADR                                = 11740,             
    IMSLS_FNORM                                    = 11750,             
    IMSLS_FORCE                                    = 11760,             
    IMSLS_FORWARD                                  = 11770,             
    IMSLS_FRACTIONAL_FACTORIAL                     = 11780,             
    IMSLS_FREQUENCIES                              = 11790,             
    IMSLS_FREQUENCIES_USER                         = 11800,             
    IMSLS_FSCALE                                   = 11810,             
    IMSLS_FSCALE_ADR                               = 11820,             
    IMSLS_FULL                                     = 11830,             
    IMSLS_FULL_NORMAL                              = 11840,             
    IMSLS_FULL_TRACEBACK                           = 11850,             
    IMSLS_FUNCTION_MIN                             = 11860,             
    IMSLS_FUZZ                                     = 11870,             
    IMSLS_FUZZ_ADR                                 = 11880,             
    IMSLS_FVALUE                                   = 11890,             
    IMSLS_FVEC                                     = 11900,             
    IMSLS_FVEC_USER                                = 11910,             
    IMSLS_F_TEST                                   = 11920,             
    IMSLS_F_TEST_RHS                               = 11930,             
    IMSLS_GENERALIZED_LEAST_SQUARES                = 11940,             
    IMSLS_GEN_LAGUERRE                             = 11950,             
    IMSLS_GEN_LAGUERRE_ADR                         = 11960,             
    IMSLS_GET_ERROR_FILE                           = 11970,             
    IMSLS_GET_OUTPUT_FILE                          = 11980,             
    IMSLS_GET_PRINT                                = 11990,             
    IMSLS_GET_STOP                                 = 12000,             
    IMSLS_GET_TRACEBACK                            = 12010,             
    IMSLS_GLM_INFO                                 = 12020,             
    IMSLS_GOOD_DIGIT                               = 12030,             
    IMSLS_GRAD                                     = 12040,             
    IMSLS_GRADIENT                                 = 12050,             
    IMSLS_GRADIENT_USER                            = 12051,
    IMSLS_GRADIENT_EPS                             = 12060,             
    IMSLS_GRAD_TOL                                 = 12070,             
    IMSLS_GRAD_TOL_ADR                             = 12080,             
    IMSLS_GROUP_COUNTS                             = 12090,             
    IMSLS_GROUP_COUNTS_USER                        = 12100,             
    IMSLS_GROUP_MEANS                              = 12110,             
    IMSLS_GROUP_MEANS_USER                         = 12120,             
    IMSLS_GROUP_STD_DEVS                           = 12130,             
    IMSLS_GROUP_STD_DEVS_USER                      = 12140,             
    IMSLS_GUESS                                    = 12150,             
    IMSLS_GVALUE                                   = 12160,             
    IMSLS_HALF_NORMAL                              = 12170,             
    IMSLS_HERMITE                                  = 12180,             
    IMSLS_HIGHEST                                  = 12190,             
    IMSLS_HINIT                                    = 12200,             
    IMSLS_HINIT_ADR                                = 12210,             
    IMSLS_HMAX                                     = 12220,             
    IMSLS_HMAX_ADR                                 = 12230,             
    IMSLS_HMIN                                     = 12240,             
    IMSLS_HMIN_ADR                                 = 12250,             
    IMSLS_HTRIAL                                   = 12260,             
    IMSLS_H_COL_DIM                                = 12270,             
    IMSLS_IMAGE                                    = 12280,             
    IMSLS_INDEPENDENT_VARIABLES                    = 12290,             
    IMSLS_INDICES_EFFECTS                          = 12300,             
    IMSLS_INDICES_KEYS                             = 12310,             
    IMSLS_INFO                                     = 12320,             
    IMSLS_INFO_USER                                = 12330,             
    IMSLS_INF_NORM                                 = 12340,             
    IMSLS_INITIAL_ESTIMATES                        = 12350,             
    IMSLS_INITIAL_TRUST_REGION                     = 12360,             
    IMSLS_INIT_HESSIAN                             = 12370,             
    IMSLS_INIT_TRUST_REGION                        = 12380,             
    IMSLS_INIT_TRUST_REGION_ADR                    = 12390,             
    IMSLS_INTERCEPT                                = 12400,             
    IMSLS_INTERMEDIATE_STEP                        = 12410,             
    IMSLS_INTERN_SCALE                             = 12420,             
    IMSLS_INTERRUPT_1                              = 12430,             
    IMSLS_INTERRUPT_2                              = 12440,             
    IMSLS_INVA_COL_DIM                             = 12450,             
    IMSLS_INVERSE                                  = 12460,             
    IMSLS_INVERSE_ONLY                             = 12470,             
    IMSLS_INVERSE_USER                             = 12480,             
    IMSLS_INV_COL_DIM                              = 12490,             
    IMSLS_ITMAX                                    = 12500,             
    IMSLS_JACOBI                                   = 12510,             
    IMSLS_JACOBIAN                                 = 12520,             
    IMSLS_JACOBIAN_MATRIX                          = 12530,             
    IMSLS_JACOBIAN_VALUES                          = 12540,             
    IMSLS_JACOBI_ADR                               = 12550,             
    IMSLS_JAC_TOL                                  = 12560,             
    IMSLS_JTJ_INVERSE                              = 12570,             
    IMSLS_JTJ_INVERSE_USER                         = 12580,             
    IMSLS_JTJ_INV_COL_DIM                          = 12590,             
    IMSLS_KNOTS                                    = 12600,             
    IMSLS_KNOTS_USER                               = 12610,             
    IMSLS_KNOWN_BOUNDS                             = 12620,             
    IMSLS_LACK_OF_FIT_EXACT_REPS                   = 12630,             
    IMSLS_LACK_OF_FIT_EXACT_REPS_USER              = 12640,             
    IMSLS_LAMBDA_HAT                               = 12650,             
    IMSLS_LAST_GRADIENT                            = 12660,             
    IMSLS_LAST_GRADIENT_USER                       = 12670,             
    IMSLS_LAST_STEP                                = 12680,             
    IMSLS_LAST_STEP_USER                           = 12690,             
    IMSLS_LEAST_SQUARES                            = 12700,             
    IMSLS_LEFT                                     = 12710,             
    IMSLS_LEFT_ADR                                 = 12720,             
    IMSLS_LEGENDRE                                 = 12730,             
    IMSLS_LEVERAGE                                 = 12740,             
    IMSLS_LEVERAGE_USER                            = 12750,             
    IMSLS_LIBRARY_VERSION                          = 12760,             
    IMSLS_LICENSE_NUMBER                           = 12770,             
    IMSLS_LILLIEFORS                               = 12780,             
    /*  Added for lilliefors_normality_test in v7.0 */
    IMSLS_MAX_DIFFERENCE                           = 12781,
    /* ** */
    IMSLS_LISTWISE                                 = 12790,             
    IMSLS_LIST_CELLS                               = 12800,             
    IMSLS_LIST_CELLS_USER                          = 12810,             
    IMSLS_LOF_CRITERION                            = 12820,             
    IMSLS_LOG_PENALIZED_LIKELIHOOD                 = 12830,             
    IMSLS_LOST                                     = 12840,             
    IMSLS_LOWER_BOUND                              = 12850,             
    IMSLS_LOWEST                                   = 12860,             
    IMSLS_LP_NORM                                  = 12870,             
    IMSLS_LRT                                      = 12880,             
    IMSLS_MALLOWS_CP                               = 12890,             
    IMSLS_MAXIMUM_LIKELIHOOD                       = 12900,             
    IMSLS_MAXORD                                   = 12910,             
    IMSLS_MAX_CYCLES                               = 12920,             
    IMSLS_MAX_EVALS                                = 12930,             
    IMSLS_MAX_FCN                                  = 12940,             
    IMSLS_MAX_GRAD                                 = 12950,             
    IMSLS_MAX_ITER                                 = 12960,             
    IMSLS_MAX_ITERATIONS                           = 12970,             
    IMSLS_MAX_ITN                                  = 12980,             
    IMSLS_MAX_JACOBIAN                             = 12990,             
    IMSLS_MAX_JACOBIAN_EVALUATIONS                 = 13000,             
    IMSLS_MAX_MOMENTS                              = 13010,             
    IMSLS_MAX_NUMBER_FCN_EVALS                     = 13020,             
    IMSLS_MAX_NUMBER_STEPS                         = 13030,             
    IMSLS_MAX_N_BEST                               = 13040,             
    IMSLS_MAX_N_GOOD_SAVED                         = 13050,             
    IMSLS_MAX_SSE_EVALUATIONS                      = 13060,             
    IMSLS_MAX_STEP                                 = 13070,             
    IMSLS_MAX_STEPS_LINE_SEARCH                    = 13080,             
    IMSLS_MAX_STEP_ADR                             = 13090,             
    IMSLS_MAX_SUBINTER                             = 13100,             
    IMSLS_MAX_TERMS                                = 13101,
    IMSLS_MA_LAGS                                  = 13110,             
    IMSLS_MEANS                                    = 13120,             
    IMSLS_MEANS_USER                               = 13130,             
    IMSLS_MEAN_AND_VARIANCE                        = 13140,             
    IMSLS_MEDIAN                                   = 13150,             
    IMSLS_MEDIAN_AND_SCALE                         = 13160,             
    IMSLS_METHOD                                   = 13170,             
    IMSLS_METHOD_OF_MOMENTS                        = 13180,             
    IMSLS_MIN_PROJECTION                           = 13190,             
    IMSLS_MITER                                    = 13200,             
    IMSLS_MODEL_ORDER                              = 13210,             
    IMSLS_MORAN                                    = 13220,             
    IMSLS_MULTIVARIATE                             = 13230,             
    IMSLS_MULTIVARIATE_DEPENDENT                   = 13240,             
    IMSLS_MULTIVARIATE_TEST                        = 13250,             
    IMSLS_MULTIVARIATE_TEST_USER                   = 13260,             
    IMSLS_N                                        = 13270,             
    IMSLS_NEAR_REPS                                = 13280,             
    IMSLS_NFCN                                     = 13290,             
    IMSLS_NFCNJ                                    = 13300,             
    IMSLS_NORM                                     = 13310,             
    IMSLS_NO_BALANCED_DATA_RESTRICTIONS            = 13320,             
    IMSLS_NO_COL_LABELS                            = 13330,             
    IMSLS_NO_CONSTANT                              = 13340,             
    IMSLS_NO_INTERCEPT                             = 13350,             
    IMSLS_NO_ROW_LABELS                            = 13360,             
    IMSLS_NO_STANDARD_ERRORS                       = 13370,             
    IMSLS_NSTEP                                    = 13380,             
    IMSLS_NUM_ROOTS                                = 13390,             
    IMSLS_N_CASES_MISSING                          = 13400,             
    IMSLS_N_CYCLES                                 = 13410,             
    IMSLS_N_EVALS                                  = 13420,             
    IMSLS_N_ITERATIONS                             = 13430,             
    IMSLS_N_MISSING                                = 13440,             
    IMSLS_N_PARAMETERS_ESTIMATED                   = 13450,             
    IMSLS_N_POSITIVE_DEVIATIONS                    = 13460,             
    IMSLS_N_SAMPLE                                 = 13470,             
    IMSLS_N_SIMULATIONS                            = 13480,             
    IMSLS_N_STEPS                                  = 13490,             
    IMSLS_N_SUBINTER                               = 13500,             
    IMSLS_N_USER                                   = 13510,             
    IMSLS_N_X1_MISSING                             = 13520,             
    IMSLS_N_X2_MISSING                             = 13530,             
    IMSLS_N_ZERO_DEVIATION                         = 13540,             
    IMSLS_OBJ                                      = 13550,             
    IMSLS_ONE_AT_A_TIME                            = 13560,             
    IMSLS_ONE_NORM                                 = 13570,             
    IMSLS_OPT                                      = 13580,             
    IMSLS_OPTIMIZE                                 = 13590,             
    IMSLS_OPT_ITMAX                                = 13600,             
    IMSLS_ORDER                                    = 13610,             
    IMSLS_ORDERS                                   = 13620,             
    IMSLS_OS_VERSION                               = 13630,             
    IMSLS_OVERALL_TESTS                            = 13640,             
    IMSLS_OVERALL_TESTS_USER                       = 13650,             
    IMSLS_PAIRWISE                                 = 13660,             
    IMSLS_PAIRWISE_USER                            = 13670,             
    IMSLS_PARAMETER_EST_COVARIANCES                = 13680,             
    IMSLS_PARAMS                                   = 13690,             
    IMSLS_PARTIAL_ACF                              = 13700,             
    IMSLS_PASSIVE                                  = 13710,             
    IMSLS_PERCENTAGE                               = 13720,             
    IMSLS_PERCENTILE                               = 13730,             
    IMSLS_PERIODIC                                 = 13740,             
    IMSLS_PERMUTATION                              = 13750,             
    IMSLS_PERMUTATION_USER                         = 13760,             
    IMSLS_PIVOT                                    = 13770,             
    IMSLS_POINTWISE_CI                             = 13780,             
    IMSLS_POINTWISE_CI_USER                        = 13790,             
    IMSLS_POINTWISE_CI_NEW_SAMPLE                  = 13800,             
    IMSLS_POINTWISE_CI_NEW_SAMPLE_USER             = 13810,             
    IMSLS_POINTWISE_CI_POP_MEAN                    = 13820,             
    IMSLS_POINTWISE_CI_POP_MEAN_USER               = 13830,             
    IMSLS_POLY_REGRESSION_INFO                     = 13840,             
    IMSLS_POOLED_VARIANCE                          = 13850,             
    IMSLS_POOL_INTERACTIONS_INTO_ERROR             = 13860,             
    IMSLS_PRECOND                                  = 13870,             
    IMSLS_PRINCIPAL_FACTOR                         = 13880,             
    IMSLS_PRINCIPAL_COMPONENT                      = 13890,             
    IMSLS_PRINT                                    = 13900,             
    IMSLS_PRINT_ALL                                = 13910,             
    IMSLS_PRINT_BRIEF                              = 13920,             
    IMSLS_PRINT_LOWER                              = 13930,             
    IMSLS_PRINT_LOWER_NO_DIAG                      = 13940,             
    IMSLS_PRINT_NONE                               = 13950,             
    IMSLS_PRINT_UPPER                              = 13960,             
    IMSLS_PRINT_UPPER_NO_DIAG                      = 13970,             
    IMSLS_PURE_ERROR                               = 13980,             
    IMSLS_P_COL_DIM                                = 13990,             
    IMSLS_Q                                        = 14000,             
    IMSLS_QUANTILE_BOUNDS                          = 14010,             
    IMSLS_QUANTILE_BOUNDS_USER                     = 14020,             
    IMSLS_QUANTILE_PROPORTIONS                     = 14030,             
    IMSLS_Q_COL_DIM                                = 14040,             
    IMSLS_Q_USER                                   = 14050,             
    IMSLS_RANDOM_EFFECTS                           = 14060,             
    IMSLS_RANDOM_FACTORS                           = 14070,             
    IMSLS_RANDOM_SPLIT                             = 14080,             
    IMSLS_RANGE                                    = 14090,             
    IMSLS_RANGE_ADR                                = 14100,             
    IMSLS_RANK                                     = 14110,             
    IMSLS_RANKS                                    = 14120,             
    IMSLS_RANK_ADR                                 = 14130,             
    IMSLS_REGRESSION_INFO                          = 14140,             
    IMSLS_REGRESSORS_EXACT_REPS                    = 14150,             
    IMSLS_REL_ERR                                  = 14160,             
    IMSLS_REL_FCN_TOL                              = 14170,             
    IMSLS_REL_FCN_TOL_ADR                          = 14180,             
    IMSLS_RESIDUAL                                 = 14190,             
    IMSLS_RESIDUAL_INPUT                           = 14200,             
    IMSLS_RESIDUAL_SCALE                           = 14210,             
    IMSLS_RESIDUAL_USER                            = 14220,             
    IMSLS_RESULT                                   = 14230,             
    IMSLS_RESULT_DUAL                              = 14240,             
    IMSLS_RESULT_USER                              = 14250,             
    IMSLS_RETURN_COL_DIM                           = 14260,             
    IMSLS_RETURN_NUMBER                            = 14270,             
    IMSLS_RETURN_USER                              = 14280,             
    IMSLS_RIGHT                                    = 14290,             
    IMSLS_RIGHT_ADR                                = 14300,             
    IMSLS_ROWS_ARE_CASES                           = 14310,             
    IMSLS_ROWS_ARE_VARIABLES                       = 14320,             
    IMSLS_ROW_LABELS                               = 14330,             
    IMSLS_ROW_NUMBER                               = 14340,             
    IMSLS_ROW_NUMBER_ZERO                          = 14350,             
    IMSLS_RULE                                     = 14360,             
    IMSLS_R_SQUARED                                = 14370,             
    IMSLS_R_SQUARED_CRITERION                      = 14380,             
    IMSLS_SAVAGE_SCORES                            = 14390,             
    IMSLS_SCALE                                    = 14400,             
    IMSLS_SCALE_ADR                                = 14410,             
    IMSLS_SCHEFFE_CI                               = 14420,    
    IMSLS_SCHEFFE_CI_USER                          = 14425,         
    IMSLS_SCHEFFE_PARAMETERIZATION                 = 14430,             
    IMSLS_SCHEFFE_USER                             = 14440,             
    IMSLS_SCORE_OPTION                             = 14450,             
    IMSLS_SECOND_VECTOR                            = 14460,             
    IMSLS_SEQUENTIAL_TESTS                         = 14470,             
    IMSLS_SEQUENTIAL_TESTS_USER                    = 14480,             
    IMSLS_SET_ERROR_FILE                           = 14490,             
    IMSLS_SET_FIRST_TO_NAN                         = 14500,             
    IMSLS_SET_OUTPUT_FILE                          = 14510,             
    IMSLS_SET_PRINT                                = 14520,             
    IMSLS_SET_STOP                                 = 14530,             
    IMSLS_SET_TRACEBACK                            = 14540,             
    IMSLS_SHAPIRO_WILK_W                           = 14550,
    IMSLS_SHAPIRO_WILK                             = 14551,
    IMSLS_SHIFT                                    = 14560,  
    IMSLS_SHIFT_ADR                                = 14565,             
    IMSLS_SIMULATED_ENVELOPE                       = 14570,             
    IMSLS_SIMULATED_ENVELOPE_USER                  = 14580,             
    IMSLS_SMOOTHING_PAR                            = 14590,             
    IMSLS_SMOOTHING_PAR_ADR                        = 14600,             
    IMSLS_SMPAR                                    = 14610,             
    IMSLS_SOLUTION_USER                            = 14620,             
    IMSLS_SOLVE_ONLY                               = 14630,             
    IMSLS_SSE                                      = 14640,             
    IMSLS_SSE_ABS_EPS                              = 14650,             
    IMSLS_SSE_REL_EPS                              = 14660,             
    IMSLS_SSQ_LOF                                  = 14670,             
    IMSLS_SSQ_LOF_COL_DIM                          = 14680,             
    IMSLS_SSQ_LOF_USER                             = 14690,             
    IMSLS_SSQ_POLY                                 = 14700,             
    IMSLS_SSQ_POLY_COL_DIM                         = 14710,             
    IMSLS_SSQ_POLY_USER                            = 14720,             
    IMSLS_SSQ_PURE_ERROR                           = 14730,             
    IMSLS_SS_RESIDUAL                              = 14740,             
    IMSLS_STANDARDIZED_RESIDUAL                    = 14750,             
    IMSLS_STANDARDIZED_RESIDUAL_USER               = 14760,             
    IMSLS_STAT                                     = 14770,             
    IMSLS_STATISTICS                               = 14780,             
    IMSLS_STATISTICS_USER                          = 14790,             
    IMSLS_STAT_COL_DIM                             = 14800,             
    IMSLS_STAT_USER                                = 14810,             
    IMSLS_STDEV_CORRELATION_MATRIX                 = 14820,             
    IMSLS_STD_DEV                                  = 14830,             
    IMSLS_STD_DEVS                                 = 14840,             
    IMSLS_STEP                                     = 14850,             
    IMSLS_STEPWISE                                 = 14860,             
    IMSLS_STEP_ADR                                 = 14870,             
    IMSLS_STEP_EPS                                 = 14880,             
    IMSLS_STEP_TOL                                 = 14890,             
    IMSLS_STEP_TOL_ADR                             = 14900,             
    IMSLS_SUR_COL_DIM                              = 14910,             
    IMSLS_SWITCH_EXACT_HESSIAN                     = 14920,             
    IMSLS_S_SQUARED                                = 14930,             
    IMSLS_S_USER                                   = 14940,             
    IMSLS_TABLE                                    = 14950,             
    IMSLS_TABLE_USER                               = 14960,             
    IMSLS_TEST_EFFECTS                             = 14970,             
    IMSLS_TEST_EFFECTS_USER                        = 14980,             
    IMSLS_TEST_STATISTIC                           = 14990,             
    IMSLS_THETA_GUESS                              = 15000,             
    IMSLS_THETA_SCALE                              = 15010,             
    IMSLS_TIES_OPTION                              = 15020,             
    IMSLS_TOL                                      = 15030,             
    IMSLS_TOLERANCE                                = 15040,             
    IMSLS_TOLERANCE_ADR                            = 15050,             
    IMSLS_TOL_ADR                                  = 15060,             
    IMSLS_TRANSPOSE                                = 15070,             
    IMSLS_TUCKER_RELIABILITY_COEFFICIENT           = 15080,             
    IMSLS_TUKEY                                    = 15090,             
    IMSLS_TUKEY_SCORES                             = 15100,             
    IMSLS_TWO_FIXED_POINTS                         = 15110,             
    IMSLS_TWO_FIXED_POINTS_ADR                     = 15120,             
    IMSLS_T_TEST                                   = 15130,             
    IMSLS_T_TEST_ALTERNATIVE                       = 15140,             
    IMSLS_T_TEST_DIFF_FOR_EQUAL_VARIANCES          = 15150,             
    IMSLS_T_TEST_DIFF_FOR_UNEQUAL_VARIANCES        = 15160,             
    IMSLS_T_TEST_NULL                              = 15170,             
    IMSLS_U                                        = 15180,             
    IMSLS_UNBALANCED_PARAMETERIZATION              = 15190,             
    IMSLS_UNEQUAL_NUMBERS                          = 15200,             
    IMSLS_UNIQUE_VARIANCES_INPUT                   = 15210,             
    IMSLS_UNIQUE_VARIANCES_OUTPUT                  = 15220,             
    IMSLS_UNIVARIATE                               = 15230,             
    IMSLS_UNWEIGHTED_LEAST_SQUARES                 = 15240,             
    IMSLS_UPPER_BOUND                              = 15250,             
    IMSLS_UPPER_LIMIT                              = 15260,             
    IMSLS_U_COL_DIM                                = 15270,             
    IMSLS_U_USER                                   = 15280,             
    IMSLS_V                                        = 15290,             
    IMSLS_VAN_DER_WAERDEN_SCORES                   = 15300,             
    IMSLS_VARIANCE_COMPONENTS                      = 15310,             
    IMSLS_VARIANCE_COMPONENTS_USER                 = 15320,             
    IMSLS_VARIANCE_COVARIANCE_MATRIX               = 15330,             
    IMSLS_VECTORS                                  = 15340,             
    IMSLS_VECTORS_USER                             = 15350,             
    IMSLS_VNORM                                    = 15360,             
    IMSLS_V_COL_DIM                                = 15370,             
    IMSLS_V_USER                                   = 15380,             
    IMSLS_WEIGHT                                   = 15390,             
    IMSLS_WEIGHTS                                  = 15400,             
    IMSLS_WEIGHTS_USER                             = 15401,
    IMSLS_WRITE_FORMAT                             = 15410,             
    IMSLS_X                                        = 15420,             
    IMSLS_XGUESS                                   = 15430,             
    IMSLS_XGUESS_ADR                               = 15440,             
    IMSLS_XSCALE                                   = 15450,             
    IMSLS_CLASS_COL_DIM                            = 15460,             
    IMSLS_X_COL_DIM                                = 15470,             
    IMSLS_CONTINUOUS_COL_DIM                       = 15480,             
    IMSLS_X_MEAN                                   = 15490,             
    IMSLS_X_MEAN_USER                              = 15500,             
    IMSLS_X_USER                                   = 15510,             
    IMSLS_X_VARIANCE                               = 15520,             
    IMSLS_X_VECTOR                                 = 15530,
    IMSLS_Y                                        = 15535,             
    IMSLS_Y_COL_DIM                                = 15540,             
    IMSLS_Y_INPUT                                  = 15550,             
    IMSLS_Y_MEAN                                   = 15560,             
    IMSLS_Y_VARIANCE                               = 15570,             
    IMSLS_Y_VECTOR                                 = 15580,
    IMSLS_COEF_COL_DIM                             = 15600,
    IMSLS_COEF_STATISTICS                          = 15610,
    IMSLS_POOL_INTERACTIONS                        = 15620,
    IMSLS_X_DIMENSIONS                             = 15630,
    IMSLS_N_OBSERVATIONS                           = 15640,
    IMSLS_N_VARIABLES                              = 15650,
    IMSLS_NO_PRINT                                 = 15660,
    IMSLS_INPUT_COV                                = 15670,
    IMSLS_LEVEL                                    = 15680,
    IMSLS_P_VALUE_IN                               = 15690,
    IMSLS_P_VALUE_OUT                              = 15700,
    IMSLS_OUTPUT_INFO_USER                         = 15710,
    IMSLS_CLUSTER_SEEDS_COL_DIM                    = 15720,
    IMSLS_CLUSTER_MEANS_COL_DIM                    = 15730,
    IMSLS_CLUSTER_MEANS_USER                       = 15740,
    IMSLS_CLUSTER_VARIABLE_COLUMNS                 = 15750,
    IMSLS_CLUSTER_COUNTS                           = 15760,
    IMSLS_CLUSTER_COUNTS_USER                      = 15770,
    IMSLS_CUM_PERCENT_USER                         = 15780,
    IMSLS_EIGENVECTORS_USER                        = 15790,
    IMSLS_CORRELATIONS_USER                        = 15800,
    IMSLS_STD_DEV_USER                             = 15810,
    IMSLS_COVARIANCE_MATRIX                        = 15820,
    IMSLS_N_ZERO_DEVIATIONS                        = 15830,
    IMSLS_JACOBIAN_MATRIX_COL_DIM                  = 15840,
    IMSLS_JACOBIAN_MATRIX_RANK                     = 15850,
    IMSLS_JACOBIAN_MATRIX_USER                     = 15860,
    IMSLS_SCHEFFE                                  = 15870,
    IMSLS_IEND                                     = 15880,
    IMSLS_SWEPT_USER                               = 15890,
    IMSLS_HISTORY_USER                             = 15900,
    IMSLS_COV_SWEPT_USER                           = 15910,
    IMSLS_REGRESSORS                               = 16000,
    IMSLS_REGRESSORS_USER                          = 16010,
    IMSLS_REGRESSORS_COL_DIM                       = 16020,
    IMSLS_X_CLASS_COLUMNS                          = 16030,
    IMSLS_SIGNAL_TRAPPING_ON                       = 16040,
    IMSLS_SIGNAL_TRAPPING_OFF                      = 16041,
    IMSLS_SET_SIGNAL_TRAPPING                      = 16042,
    IMSLS_RELATIVE_ERROR                           = 16050,
    IMSLS_MEAN_ESTIMATE                            = 16051,
    IMSLS_PARAM_EST_COV                            = 16052,
    IMSLS_PARAM_EST_COV_USER                       = 16053,
    IMSLS_AUTOCOV                                  = 16054,
    IMSLS_AUTOCOV_USER                             = 16055,
    IMSLS_ARMA_INFO                                = 16056,
    IMSLS_TUKEY_USER                               = 16057,
    IMSLS_DUNN_SIDAK_USER                          = 16058,
    IMSLS_BONFERRONI_USER                          = 16059,
    IMSLS_ONE_AT_A_TIME_USER                       = 16061,
    IMSLS_CRITERIONS_USER                          = 16065,
    IMSLS_INDEPENDENT_VARIABLES_USER               = 16066,
    IMSLS_COEF_STATISTICS_USER                     = 16067,
    IMSLS_R                                        = 16068,
    IMSLS_R_USER                                   = 16069,
    IMSLS_R_RANK                                   = 16070,
    IMSLS_DFE                                      = 16071,
    IMSLS_R_COL_DIM                                = 16072,
    IMSLS_T_TEST_FOR_EQUAL_VARS                    = 16073,
    IMSLS_T_TEST_FOR_UNEQUAL_VARS                  = 16074,
    IMSLS_CI_RATIO_VARIANCES                       = 16075,
    IMSLS_CI_DIFF_FOR_EQUAL_VARS                   = 16076,
    IMSLS_CI_DIFF_FOR_UNEQUAL_VARS                 = 16077,
    IMSLS_MISSING_VALUE_METHOD                     = 16078,
    IMSLS_PREDICTED                                = 16079,
    IMSLS_PREDICTED_USER                           = 16080,
    IMSLS_INCIDENCE_MATRIX                         = 16081,
    IMSLS_INCIDENCE_MATRIX_USER                    = 16082,
    IMSLS_CONFIDENCE_ADR                           = 18000,
    IMSLS_MISSING_LISTWISE                         = 18010,
    IMSLS_MISSING_ELEMENTWISE                      = 18011,
    IMSLS_CONFIDENCE_MEAN_ADR                      = 20001,
    IMSLS_T_TEST_NULL_ADR                          = 20002,
    IMSLS_CONFIDENCE_VARIANCE_ADR                  = 20003,
    IMSLS_CHI_SQUARED_TEST_NULL_ADR                = 20004,
    IMSLS_PERCENTAGE_ADR                           = 20005,
    IMSLS_PERCENTILE_ADR                           = 20006,
    IMSLS_P_VALUE_IN_ADR                           = 20007,
    IMSLS_P_VALUE_OUT_ADR                          = 20008,
    IMSLS_GRADIENT_EPS_ADR                         = 20009,
    IMSLS_STEP_EPS_ADR                             = 20010,
    IMSLS_SSE_REL_EPS_ADR                          = 20011,
    IMSLS_SSE_ABS_EPS_ADR                          = 20012,
    IMSLS_INITIAL_TRUST_REGION_ADR                 = 20014,
    IMSLS_CONVERGENCE_EPS_ADR                      = 20015,
    IMSLS_SWITCH_EXACT_HESSIAN_ADR                 = 20016,
    IMSLS_BACKCASTING_ADR                          = 20017,
    IMSLS_CONVERGENCE_TOLERANCE_ADR                = 20018,
    IMSLS_RELATIVE_ERROR_ADR                       = 20019,
    IMSLS_MEAN_ESTIMATE_ADR                        = 20020,
    IMSLS_MODEL                                    = 20100,
    IMSLS_INTERVAL                                 = 20110,
    IMSLS_X_COL_FIXED_PARAMETER                    = 20120,
    IMSLS_X_COL_FREQUENCIES                        = 20130,
    IMSLS_RETURN_NEW_X                             = 20140,
    IMSLS_RETURN_NEW_X_USER                        = 20141,
    IMSLS_X_COL_DIST_PARAMETER                     = 20150,
    IMSLS_X_COL_RESPONSE                           = 20160,
    IMSLS_INFINITY_CHECK                           = 20170,
    IMSLS_NO_INFINITY_CHECK                        = 20180,
    IMSLS_CREATE_NEW_X                             = 20190,
    IMSLS_EFFECTS                                  = 20200,
    IMSLS_INITIAL_EST_INPUT                        = 20210,
    IMSLS_INITIAL_EST_INTERNAL                     = 20220,
    IMSLS_INITIAL_EST_OUTPUT                       = 20225, 
    IMSLS_MAX_CLASS                                = 20230,
    IMSLS_CLASS_INFO_USER                          = 20240,
    IMSLS_CLASS_INFO                               = 20250,
    IMSLS_COEF_STAT_USER                           = 20260,
    IMSLS_COEF_STAT                                = 20270,
    IMSLS_COEF_STAT_COL_DIM                        = 20280,
    IMSLS_CRITERION                                = 20290,
    IMSLS_COV_USER                                 = 20300,
    IMSLS_COV                                      = 20310,
    IMSLS_CASE_ANALYSIS                            = 20320,
    IMSLS_CASE_ANALYSIS_USER                       = 20330,
    IMSLS_CASE_ANALYSIS_COL_DIM                    = 20340,
    IMSLS_OBS_STATUS_USER                          = 20370,
    IMSLS_OBS_STATUS                               = 20380,
    IMSLS_LP_MAX                                   = 20390,
    IMSLS_N_ROWS_MISSING                           = 20400,
    IMSLS_MEAN                                     = 20410,
    IMSLS_MEAN_ADR                                 = 20415,
    IMSLS_VARIANCE                                 = 20420,
    IMSLS_VARIANCE_ADR                             = 20425,
    IMSLS_ACCEPT_REJECT_METHOD                     = 20430,
    IMSLS_IDO                                      = 20440,
    IMSLS_ROWS_ADD                                 = 20441,
    IMSLS_ROWS_DELETE                              = 20442,
    IMSLS_X_INDICES                                = 20443,
    IMSLS_PRIOR_PROPORTIONAL                       = 20444,
    IMSLS_PRIOR_INPUT                              = 20445,
    IMSLS_PRIOR_OUTPUT                             = 20446,
    IMSLS_PRIOR_OUTPUT_USER                        = 20447,
    IMSLS_COEF                                     = 20448,
    IMSLS_COEF_USER                                = 20449,
    IMSLS_CLASS_MEMBERSHIP                         = 20450,
    IMSLS_CLASS_MEMBERSHIP_USER                    = 20451,
    IMSLS_CLASS_TABLE                              = 20452,
    IMSLS_CLASS_TABLE_USER                         = 20453,
    IMSLS_MAHALANOBIS                              = 20454,
    IMSLS_MAHALANOBIS_USER                         = 20455,
    IMSLS_STATS                                    = 20456,
    IMSLS_STATS_USER                               = 20457,
    IMSLS_PROB                                     = 20458,
    IMSLS_PROB_USER                                = 20459,
    IMSLS_PRIOR_EQUAL                              = 20460,
    IMSLS_SIMPLE_LOWER_BOUNDS                      = 20470,
    IMSLS_SIMPLE_UPPER_BOUNDS                      = 20480,
    IMSLS_ACC                                      = 20490,
    IMSLS_ACC_ADR                                  = 20491,
    IMSLS_STOP_INFO                                = 20500,
    IMSLS_ACTIVE_CONSTRAINTS_INFO                  = 20510,
    IMSLS_ACTIVE_CONSTRAINTS_INFO_USER             = 20520,
    IMSLS_PRINT_LEVEL                              = 20530,
    IMSLS_LINEAR_CONSTRAINTS                       = 20540,
    IMSLS_ALL_COLUMNS                              = 20550,
    IMSLS_FIRST_N_COLUMNS                          = 20560,
    IMSLS_SPECIFY_COLUMNS                          = 20570,
    IMSLS_SWAP                                     = 20580,
    IMSLS_SWAP_USER                                = 20590,
    IMSLS_Q_STATISTIC                              = 20600,
    IMSLS_PROB_TABLE                               = 20610,
    IMSLS_P_VALUE                                  = 20620,
    IMSLS_CHECK_NUMERICAL_ERROR                    = 20630,
    IMSLS_WORKSPACE                                = 20640,
    IMSLS_APPROXIMATION_PARAMETERS                 = 20650, 
    IMSLS_NO_APPROXIMATION                         = 20660,
    IMSLS_B                                        = 20670,
    IMSLS_B_ADR                                    = 20675,
    IMSLS_X_COL_CENSORING                          = 20680,
    IMSLS_ITERATIONS                               = 20681,
    IMSLS_ITERATIONS_USER                          = 20682,
    IMSLS_X_COL_VARIABLES                          = 20683,
    IMSLS_SURVIVAL_INFO                            = 20684, 
    IMSLS_XBETA                                    = 20685,
    IMSLS_XBETA_USER                               = 20686,
    IMSLS_XPT_COLUMNS                              = 20689,
    IMSLS_PARTIAL_CORR                             = 20690,
    IMSLS_PARTIAL_COV                              = 20691,
    IMSLS_TEST                                     = 20692,
    IMSLS_TEST_USER                                = 20693,
    IMSLS_WEIGHT_FUNCTION                          = 20700,
    IMSLS_GROUPS_INPUT                             = 20710,
    IMSLS_GROUPS_OUTPUT                            = 20720,
    IMSLS_INIT_MEAN                                = 20730,
    IMSLS_INIT_MEDIAN                              = 20740,
    IMSLS_INIT_INPUT                               = 20750,
    IMSLS_HUBER                                    = 20760,
    IMSLS_STAHEL                                   = 20770,
    IMSLS_CONSTANTS                                = 20780,
    IMSLS_CONSTANTS_USER                           = 20790,
    IMSLS_SUM_WEIGHTS                              = 20800,
    IMSLS_GROUPS                                   = 20810,
    IMSLS_SUM_WEIGHTS_USER                         = 20820,
    IMSLS_INITIAL_EST_MEAN                         = 20830,
    IMSLS_INITIAL_EST_MEDIAN                       = 20840,
    IMSLS_MINIMAX_WEIGHTS                          = 20850,
    IMSLS_USER_WEIGHTS                             = 20860,
    IMSLS_BETA                                     = 20870,
    IMSLS_NO_CENTERING                             = 20880,
    IMSLS_CENTERING                                = 20890,
    IMSLS_N_DEPENDENT                              = 20900,
    IMSLS_BETA_MATRIX                              = 20910,
    IMSLS_BETA_MATRIX_USER                         = 20920,
    IMSLS_SCPE                                     = 20930,
    IMSLS_SCPE_USER                                = 20940,
    IMSLS_INDEX_REGRESSION                         = 20941,
    IMSLS_DFH                                      = 20942,
    IMSLS_G                                        = 20943,
    IMSLS_G_USER                                   = 20944,
    IMSLS_GP                                       = 20945,
    IMSLS_RANK_HP                                  = 20946,
    IMSLS_H_MATRIX_USER                            = 20947,
    IMSLS_H_MATRIX                                 = 20948,
    IMSLS_WILK_LAMBDA                              = 20950,
    IMSLS_ROY_MAX_ROOT                             = 20951,
    IMSLS_HOTELLING_TRACE                          = 20952,
    IMSLS_PILLAI_TRACE                             = 20953,
    IMSLS_KNOWN_BOUNDS_ADR                         = 20954,
    IMSLS_APPROXIMATION_PARAMETERS_ADR             = 20960,
    IMSLS_INVERSE_TRANSFORM                        = 25000,
    IMSLS_ARMA_CONSTANT                            = 25010,
    IMSLS_ARMA_CONSTANT_ADR                        = 25015,
    IMSLS_VAR_NOISE                                = 25020,
    IMSLS_VAR_NOISE_ADR                            = 25025,
    IMSLS_INPUT_NOISE                              = 25030,
    IMSLS_OUTPUT_NOISE                             = 25040,
    IMSLS_OUTPUT_NOISE_USER                        = 25050,
    IMSLS_NONZERO_ARLAGS                           = 25060,
    IMSLS_NONZERO_MALAGS                           = 25070,
    IMSLS_INITIAL_W                                = 25080,
    IMSLS_START_PERIOD                             = 25090,
    IMSLS_MEAN_METHOD                              = 25100,
    IMSLS_SEASONAL_LENGTH                          = 25110,
    IMSLS_SIMPLE_METHOD                            = 25115,
    IMSLS_SYMMETRIC_METHOD                         = 25120,
    IMSLS_NONSEASONAL_CONSTANT                     = 25125,
    IMSLS_NONSEASONAL_TREND                        = 25130,
    IMSLS_SEASONAL_CONSTANT                        = 25135,
    IMSLS_SEASONAL_TREND                           = 25140,
    IMSLS_USER_SUPPLIES_VALUE                      = 25150,
    IMSLS_MAX_LAG                                  = 25160,
    IMSLS_RETURN_MAX_LAG                           = 25170,
    IMSLS_USER_SUPPLIES_MEAN                       = 25180,
    IMSLS_RETURN_USER_MEAN                         = 25190,
    IMSLS_ACOV                                     = 25200, 
    IMSLS_ACOV_USER                                = 25210,
    IMSLS_SE_SACF                                  = 25220,
    IMSLS_SE_SACF_USER                             = 25230,
    IMSLS_PORTMANTEAU_TEST                         = 25240,
    IMSLS_PORTMANTEAU_TEST_USER                    = 25250,
    IMSLS_N_PARAMETERS                             = 25260,
    IMSLS_SPACF                                    = 25270,
    IMSLS_SPACF_USER                               = 25280,
    IMSLS_ESACF                                    = 25290,
    IMSLS_ESACF_USER                               = 25300,
    IMSLS_ESACF_SYMBOL_TABLE                       = 25310, 
    IMSLS_ESACF_SYMBOL_TABLE_USER                  = 25320,
    IMSLS_ESACF_AR_LAGS                            = 25330,
    IMSLS_ESACF_MA_LAGS                            = 25340,
    IMSLS_NO_SE                                    = 25345,
    IMSLS_MORAN_SE                                 = 25350,
    IMSLS_BARTLETT_SE                              = 25355,
    IMSLS_CHANGE_LOOP_MAXIMUM                      = 25360,
    IMSLS_COMPANION_METHOD                         = 25370,
    IMSLS_NO_CENTER                                = 25380,
    IMSLS_CENTER                                   = 25390,
    IMSLS_INPUT_SERIES_MEAN                        = 25400,
    IMSLS_N_DIFFERENCES                            = 25410,
    IMSLS_DIFF_PERIODS                             = 25420,
    IMSLS_DIFF_ORDERS                              = 25430,
    IMSLS_POWER_TRANSFORM                          = 25440,
    IMSLS_ESTIMATION_METHOD                        = 25450,
    IMSLS_NONZERO_AR_LAGS                          = 25460,
    IMSLS_NONZERO_MA_LAGS                          = 25470,
    IMSLS_INNOVATIONS_ORDER                        = 25480,
    IMSLS_FCN_RELATIVE_ERROR                       = 25490,
    IMSLS_FCN_ABSOLUTE_ERROR                       = 25500,
    IMSLS_EST_RELATIVE_ERROR                       = 25510,
    IMSLS_OUTPUT_SERIES_MEAN                       = 25520,
    IMSLS_OUTPUT_SERIES_CONSTANT                   = 25530,
    IMSLS_INPUT_NOISE_VARIANCE                     = 25540,
    IMSLS_OUTPUT_NOISE_VARIANCE                    = 25550,
    IMSLS_DEGREE_OF_FREEDOM                        = 25560,
    IMSLS_RESIDUALS                                = 25570,
    IMSLS_RESIDUALS_USER                           = 25580,
    IMSLS_LIKELIHOOD_VALUE                         = 25590,
    IMSLS_FITTED_VALUES                            = 25600,
    IMSLS_FITTED_VALUES_USER                       = 25610,
    IMSLS_FITS_MSE                                 = 25620,
    IMSLS_FITS_MSE_USER                            = 25630,
    IMSLS_MODEL_CRITERION                          = 25640, 
    IMSLS_MODEL_CRITERION_USER                     = 25650,
    IMSLS_RANDOMNESS_TEST                          = 25660,
    IMSLS_RANDOMNESS_TEST_USER                     = 25670,
    IMSLS_NORMALITY_TEST                           = 25680,
    IMSLS_NORMALITY_TEST_USER                      = 25690,
    IMSLS_MAX_LEAD_TIME                            = 25700,
    IMSLS_FORECAST_CONFIDENCE                      = 25710,
    IMSLS_FORECAST_ORIGIN                          = 25720,
    IMSLS_FORECASTS                                = 25730,
    IMSLS_FORECASTS_USER                           = 25740,
    IMSLS_FORECAST_MSE                             = 25750,
    IMSLS_FORECAST_MSE_USER                        = 25760,
    IMSLS_FUNCTION_WITH_PARAMS                     = 25770,
    IMSLS_JACOBIAN_WITH_PARAMS                     = 25780,
    IMSLS_GRAD_WITH_PARAMS                         = 25790,
    IMSLS_CHECK_ROOTS                              = 25800,
    IMSLS_NO_CHECK_ROOTS                           = 25810,
    IMSLS_EST_SCALE                                = 25820,
    IMSLS_HESSIAN                                  = 25830,
    IMSLS_HESSIAN_USER                             = 25840,
    IMSLS_KALMAN_INITIAL                           = 25850,
    IMSLS_TRANSITION_MATRIX                        = 25860,
    IMSLS_RELATION_MATRIX                          = 25870,
    IMSLS_OBS_ERROR_COVARIANCE                     = 25880,
    IMSLS_STATE_ERROR_COVARIANCE                   = 25890,
    IMSLS_CONTROL_EFFECT                           = 25900,
    IMSLS_NO_UPDATE                                = 25910,
    IMSLS_UPDATE                                   = 25920,
    IMSLS_PREDICTION_ERROR_USER                    = 25930,
    IMSLS_PREDICTION_ERROR_COV_USER                = 25940,
    IMSLS_CUMULATIVE_RANK                          = 25950,
    IMSLS_CUMULATIVE_SS                            = 25960,
    IMSLS_CUMULATIVE_LOG_DET                       = 25970,
    IMSLS_PREDICTION_STATE_USER                    = 25980,
    IMSLS_PREDICTION_STATE_COV_USER                = 25990,
    IMSLS_FILTER_USER                              = 26000,
    IMSLS_FILTER_COV_USER                          = 26010,
    IMSLS_OBS_STATE_ERROR_COV                      = 26020,
    IMSLS_SUM_METHOD                               = 26030,
    IMSLS_N_AGGREGATED_OBS                         = 26040,
    IMSLS_MAX_LAG_CCF                              = 26050,
    IMSLS_INPUT_MEANS                              = 26060,
    IMSLS_CCF_SE_OPTIONS                           = 26070,
    IMSLS_N_IMPULSE_WT                             = 26080,
    IMSLS_PREWHITE_AR_OPERATOR                     = 26090,
    IMSLS_PREWHITE_MA_OPERATOR                     = 26100,
    IMSLS_N_IMPULSE_WT_FOR_NOISE                   = 26110,
    IMSLS_OUTPUT_MEANS                             = 26120,
    IMSLS_VARIANCES                                = 26130,
    IMSLS_CROSS_COVARIANCES                        = 26140,
    IMSLS_CROSS_COVARIANCES_USER                   = 26150,
    IMSLS_SE_CCF                                   = 26160,
    IMSLS_SE_CCF_USER                              = 26170,
    IMSLS_IMPULSE_WT                               = 26180,
    IMSLS_IMPULSE_WT_USER                          = 26190,
    IMSLS_NOISE_SERIES                             = 26200,
    IMSLS_NOISE_SERIES_USER                        = 26210,
    IMSLS_PREWHITE_X                               = 26220,
    IMSLS_PREWHITE_X_USER                          = 26230,
    IMSLS_PREWHITE_Y                               = 26240,
    IMSLS_PREWHITE_Y_USER                          = 26250,
    IMSLS_CCF_PREWHITE_XY                          = 26300,
    IMSLS_CCF_PREWHITE_XY_USER                     = 26310,
    IMSLS_NUMBER_PARAMETERS                        = 26320,
    IMSLS_ACV                                      = 30001,
    IMSLS_ACV_USER                                 = 30002,
    IMSLS_SEAC                                     = 30003,
    IMSLS_SEAC_USER                                = 30004,
    IMSLS_X_MEAN_IN                                = 30005,
    IMSLS_X_MEAN_OUT                               = 30006,
    IMSLS_ISE_OPTION                               = 30007,
    IMSLS_LAGMIN                                   = 30008,
    IMSLS_METHOD_LAV                               = 30009,
    IMSLS_METHOD_LLP                               = 30010,
    IMSLS_METHOD_LMV                               = 30011,
    IMSLS_SEA                                      = 30012,
    IMSLS_MAX_RESIDUAL                             = 30013,
    IMSLS_RESIDUALS_LP_NORM                        = 30014,
    IMSLS_SUM_RANK                                 = 30015,
    IMSLS_SUM_RANK_USER                            = 30016,
    IMSLS_DISPERSION                               = 30017,
    IMSLS_DIFFERENCES                              = 30018,
    IMSLS_DIFFERENCES_USER                         = 30019,
    IMSLS_N_MISSING_X                              = 30020,
    IMSLS_N_MISSING_Y                              = 30021,
    IMSLS_DCUBE                                    = 30022,
    IMSLS_DCUBE_USER                               = 30023,
    IMSLS_DSQUARE                                  = 30024,
    IMSLS_DSQUARE_USER                             = 30025,
    IMSLS_PAIRS                                    = 30026,
    IMSLS_PAIRS_USER                               = 30027,
    IMSLS_RUNS                                     = 30028,
    IMSLS_RUNS_USER                                = 30029,
    IMSLS_EXPECT                                   = 30030,
    IMSLS_RUNS_EXPECT                              = 30031,
    IMSLS_RUNS_EXPECT_USER                         = 30032,
    IMSLS_Y_MEANS                                  = 30033,
    IMSLS_Y_MEANS_USER                             = 30034,
    IMSLS_SUM_FREQ                                 = 30035,
    IMSLS_EXPECTED_MEAN_SQUARE                     = 30036,
    IMSLS_EXPECTED_MEAN_SQUARE_USER                = 30037,
    IMSLS_VAR_COL_DIM                              = 30038,
    IMSLS_VAR_USER                                 = 30039,
    IMSLS_VAR                                      = 30040,
    IMSLS_A                                        = 30041,
    IMSLS_AIC                                      = 30042,
    IMSLS_MAX_SIGMA                                = 30043,
    IMSLS_MAX_SIGMA_ADR                            = 30044,
    IMSLS_X_MEAN_IN_ADR                            = 30045,
    IMSLS_GET_INDEX_VECTORS                        = 40001,
    IMSLS_GET_INDEX_VECTORS_USER                   = 40002,
    IMSLS_SET_INDEX_VECTORS                        = 40003,
    IMSLS_INDEX_ONLY                               = 40004,
    IMSLS_TABLE_COL_DIM                            = 40005,
    IMSLS_POPULATION_COL_DIM                       = 40006,
    IMSLS_FIRST_CALL                               = 40007,
    IMSLS_FIRST_CALL_USER                          = 40008,
    IMSLS_ADDITIONAL_CALL                          = 40009,
    IMSLS_Z_COL_DIM                                = 40010,
    IMSLS_T_COL_DIM                                = 40011,
    IMSLS_COVB_COL_DIM                             = 40012,
    IMSLS_COVV_USER                                = 40013,
    IMSLS_COVV                                     = 40014,
    IMSLS_T                                        = 40015,
    IMSLS_BASE                                     = 40016,
    IMSLS_SKIP                                     = 40017,
    IMSLS_RETURN_SKIP                              = 40018,
    IMSLS_FCN_W_DATA                               = 40020,
    IMSLS_GRADIENT_W_DATA                          = 40021,
    IMSLS_JACOBIAN_W_DATA                          = 40023,
    IMSLS_VARIANCES_USER                           = 40030,
    IMSLS_OUTPUT_MEANS_USER                        = 40031,
    IMSLS_RATIO                                    = 40032,
    IMSLS_CENSORING_CODE                           = 40033,
    IMSLS_STRATIFICATION_COL                       = 40034,
    IMSLS_FREQ_RESPONSE                            = 40035,
    IMSLS_VARIANCE_COVARIANCE_MATRIX_USER          = 40036,
    IMSLS_UPDATE_USER                              = 40037,
    IMSLS_STRATUM_NUMBER                           = 40038,
    IMSLS_STRATUM_NUMBER_USER                      = 40039,
    IMSLS_CLASS_VARIABLES                          = 40040,
    IMSLS_CLASS_VARIABLES_USER                     = 40041,
    IMSLS_POPULATION_SIZE                          = 40042,
    IMSLS_POPULATION_LIFE_TABLE                    = 40043,
    IMSLS_SORTED                                   = 40044,
    IMSLS_LOCATIONS                                = 40102,
    IMSLS_RCBD                                     = 40103,
    IMSLS_CRD                                      = 40104,
    IMSLS_LOC_FIXED                                = 40105,
    IMSLS_LOC_RANDOM                               = 40106,
    IMSLS_WHOLE_FIXED                              = 40107,
    IMSLS_WHOLE_RANDOM                             = 40108,
    IMSLS_SPLIT_FIXED                              = 40109,
    IMSLS_SPLIT_RANDOM                             = 40110,
    IMSLS_EFNCY                                    = 40111,
    IMSLS_GRAND_MEAN                               = 40112,
    IMSLS_WHOLE_PLOT_MEANS                         = 40113,
    IMSLS_WHOLE_PLOT_MEANS_USER                    = 40114,
    IMSLS_SPLIT_PLOT_MEANS                         = 40115,
    IMSLS_SPLIT_PLOT_MEANS_USER                    = 40116,
    IMSLS_TREATMENT_MEANS                          = 40117,
    IMSLS_TREATMENT_MEANS_USER                     = 40118,
    IMSLS_SQSS_MODEL_USER                          = 40119,
    IMSLS_STD_ERRORS                               = 40120,
    IMSLS_STD_ERRORS_USER                          = 40121,
    IMSLS_N_BLOCKS                                 = 40122,
    IMSLS_N_BLOCKS_USER                            = 40123,
    IMSLS_BLOCK_SS                                 = 40124,
    IMSLS_BLOCK_SS_USER                            = 40125,
    IMSLS_WHOLE_PLOT_SS                            = 40126,
    IMSLS_WHOLE_PLOT_SS_USER                       = 40127,
    IMSLS_SPLIT_PLOT_SS                            = 40128,
    IMSLS_SPLIT_PLOT_SS_USER                       = 40129,
    IMSLS_WHOLE_PLOT_ERROR_SS                      = 40130,
    IMSLS_WHOLE_PLOT_ERROR_SS_USER                 = 40131,
    IMSLS_SPLIT_PLOT_ERROR_SS                      = 40132,
    IMSLS_SPLIT_PLOT_ERROR_SS_USER                 = 40133,
    IMSLS_CV                                       = 40134,
    IMSLS_CV_USER                                  = 40135,
    IMSLS_TOTAL_SS                                 = 40136,
    IMSLS_TOTAL_SS_USER                            = 40137,
    IMSLS_WHOLEXSPLIT_PLOT_SS                      = 40138,
    IMSLS_WHOLEXSPLIT_PLOT_SS_USER                 = 40139,
    IMSLS_SUB_PLOT_MEANS                           = 40140,
    IMSLS_SUB_PLOT_MEANS_USER                      = 40141,
    IMSLS_LOCATION_ANOVA_TABLE                     = 40142,
    IMSLS_LOCATION_ANOVA_TABLE_USER                = 40143,
    IMSLS_SUB_PLOT_FIXED                           = 40144,
    IMSLS_SUB_PLOT_RANDOM                          = 40145,
    IMSLS_WHOLE_SPLIT_PLOT_MEANS                   = 40146,
    IMSLS_WHOLE_SPLIT_PLOT_MEANS_USER              = 40147,
    IMSLS_WHOLE_SUB_PLOT_MEANS                     = 40148,
    IMSLS_WHOLE_SUB_PLOT_MEANS_USER                = 40149,
    IMSLS_SPLIT_SUB_PLOT_MEANS                     = 40150,
    IMSLS_SPLIT_SUB_PLOT_MEANS_USER                = 40151,
    IMSLS_SUB_FIXED                                = 40152,
    IMSLS_SUB_RANDOM                               = 40153,
    IMSLS_LEVENES_MEAN                             = 40155,
    IMSLS_LEVENES_MEDIAN                           = 40156,
    IMSLS_STUDENTIZED_RESIDUALS                    = 40157,
    IMSLS_STUDENTIZED_RESIDUALS_USER               = 40158,
    IMSLS_STD_DEVS_USER                            = 40159,
    IMSLS_BARTLETTS                                = 40160,
    IMSLS_LEVENES                                  = 40161,
    IMSLS_N_MISSING_REPLACED                       = 40170,
    IMSLS_DESIGN                                   = 40171,
    IMSLS_INTITIAL_ESTIMATES                       = 40172,
    IMSLS_GET_SS                                   = 40173,
    IMSLS_ERROR_SS                                 = 40174,
    IMSLS_MISSING_INDEX                            = 40175,
    IMSLS_MISSING_INDEX_USER                       = 40176,
    IMSLS_SNK                                      = 40180,
    IMSLS_DUNCANS_MRT                              = 40181,
    IMSLS_LSD                                      = 40182,
    IMSLS_LOCATIONS_POOL                           = 40183,
    IMSLS_LOCATIONS_NOPOOL                         = 40184,
    IMSLS_FACTOR_MEANS                             = 40186,
    IMSLS_FACTOR_MEANS_USER                        = 40187,
    IMSLS_TWO_WAY_STD_ERRORS                       = 40188,
    IMSLS_TWO_WAY_STD_ERRORS_USER                  = 40189,
    IMSLS_TWO_WAY_MEANS                            = 40190,
    IMSLS_TWO_WAY_MEANS_USER                       = 40191,
    IMSLS_FACTOR_STD_ERRORS                        = 40192,
    IMSLS_FACTOR_STD_ERRORS_USER                   = 40193,
    IMSLS_TREATMENT_STD_ERROR                      = 40194,
    IMSLS_TREATMENT_STD_ERROR_USER                 = 40195,
    IMSLS_ANOVA_ROW_LABELS                         = 40200,
    IMSLS_ANOVA_ROW_LABELS_USER                    = 40201,
    IMSLS_STRIP_PLOT_A_FIXED                       = 40205,
    IMSLS_STRIP_PLOT_B_FIXED                       = 40206,
    IMSLS_STRIP_PLOT_A_RANDOM                      = 40207,
    IMSLS_STRIP_PLOT_B_RANDOM                      = 40208,
    IMSLS_STRIP_PLOT_A_MEANS                       = 40209,
    IMSLS_STRIP_PLOT_A_MEANS_USER                  = 40210,
    IMSLS_STRIP_PLOT_B_MEANS                       = 40211,
    IMSLS_STRIP_PLOT_B_MEANS_USER                  = 40212,
    IMSLS_STRIP_PLOT_AB_MEANS                      = 40213,
    IMSLS_STRIP_PLOT_AB_MEANS_USER                 = 40214,
    IMSLS_STRIP_PLOT_A_SUB_MEANS                   = 40215,
    IMSLS_STRIP_PLOT_A_SUB_MEANS_USER              = 40216,
    IMSLS_STRIP_PLOT_B_SUB_MEANS                   = 40217,
    IMSLS_STRIP_PLOT_B_SUB_MEANS_USER              = 40218,
    IMSLS_STRIP_PLOT_A_SPLIT_PLOT_MEANS            = 40219,
    IMSLS_STRIP_PLOT_A_SPLIT_PLOT_MEANS_USER       = 40220,
    IMSLS_STRIP_PLOT_B_SPLIT_PLOT_MEANS            = 40221,
    IMSLS_STRIP_PLOT_B_SPLIT_PLOT_MEANS_USER       = 40222,
    IMSLS_STRIP_A_FIXED                            = 40223,
    IMSLS_STRIP_B_FIXED                            = 40224,
    IMSLS_STRIP_A_RANDOM                           = 40225,
    IMSLS_STRIP_B_RANDOM                           = 40226,
    IMSLS_CLUSTERS                                 = 40300,
    IMSLS_CLUSTERS_USER                            = 40301,
    IMSLS_COLUMNS                                  = 40302,
    IMSLS_INDEX                                    = 40303,
    IMSLS_OBS_PER_CLUSTER                          = 40304,
    IMSLS_OBS_PER_CLUSTER_USER                     = 40305,
    IMSLS_ROWS                                     = 40306,
    IMSLS_TRANSFORMATION                           = 40307,
    IMSLS_ORTHOMAX_ROTATION                        = 40320,
    IMSLS_ORTHOMAX_ROTATION_USER                   = 40321,
    IMSLS_ORTHOGONAL_PROCRUSTES_ROTATION           = 40322,
    IMSLS_ORTHOGONAL_PROCRUSTES_ROTATION_USER      = 40323, 
    IMSLS_DIRECT_OBLIMIN_ROTATION                  = 40324,
    IMSLS_DIRECT_OBLIMIN_ROTATION_USER             = 40325,
    IMSLS_OBLIQUE_PROMAX_ROTATION                  = 40326,
    IMSLS_OBLIQUE_PROMAX_ROTATION_USER             = 40327,
    IMSLS_OBLIQUE_PIVOTAL_PROMAX_ROTATION          = 40328,
    IMSLS_OBLIQUE_PIVOTAL_PROMAX_ROTATION_USER     = 40329,
    IMSLS_OBLIQUE_PROCRUSTES_ROTATION              = 40330,
    IMSLS_OBLIQUE_PROCRUSTES_ROTATION_USER         = 40331,
    IMSLS_FACTOR_STRUCTURE                         = 40332,
    IMSLS_FACTOR_STRUCTURE_USER                    = 40333,
    IMSLS_CENSOR_CODES                             = 40400,
    IMSLS_SORT_OPTION                              = 40401,
    IMSLS_BETA_GRID                                = 40402,
    IMSLS_K                                        = 40403,
    IMSLS_SORTED_EVENT_TIMES                       = 40404,
    IMSLS_SORTED_EVENT_TIMES_USER                  = 40405,
    IMSLS_SORTED_CENSOR_CODES                      = 40406,
    IMSLS_SORTED_CENSOR_CODES_USER                 = 40407, 
    IMSLS_K_GRID                                   = 40408,
    IMSLS_CENSOR_CODES_COL                         = 40409,
    IMSLS_X_RESPONSE_COL                           = 40410,
    IMSLS_FREQ_RESPONSE_COL                        = 40411,
    IMSLS_STRATUM_NUMBER_COL                       = 40412,
    IMSLS_CONSTANT_COL                             = 40413,
    IMSLS_INPUT_N_CLASSES                          = 40500,
    IMSLS_EQUAL_FREQUENCY                          = 40501,
    IMSLS_OUTPUT_CUTPOINTS                         = 40502,
    IMSLS_OUTPUT_CUTPOINTS_USER                    = 40503,
    IMSLS_OUTPUT_N_CLASSES                         = 40504,
    IMSLS_N_VALID                                  = 40505,
    IMSLS_INPUT_CUTPOINTS                          = 40506,
    IMSLS_SUPPLY_CENTER_SPREAD                     = 40510,
    IMSLS_RETURN_CENTER_SPREAD                     = 40511,
    IMSLS_SCALE_LIMITS                             = 40512,
    IMSLS_LAGS                                     = 40513,
    IMSLS_ENCODE                                   = 40514,
    IMSLS_DECODE                                   = 40515,
    IMSLS_NO_TRANSFORM                             = 40516,
    IMSLS_SQUARE_ROOT                              = 40517,
    IMSLS_ARC_SIN                                  = 40518,
    IMSLS_SQSS                                     = 40519,
    IMSLS_SQSS_USER                                = 40520,
    IMSLS_N_CLASSES                                = 40521,

    IMSLS_CREATE_HIDDEN_LAYER                      = 40600,
    IMSLS_ACTIVATION_FCN                           = 40601,
    IMSLS_BIAS                                     = 40602,
    IMSLS_LINK_ALL                                 = 40603,
    IMSLS_LINK_NODE                                = 40604,
    IMSLS_LINK_LAYER                               = 40605,
    IMSLS_REMOVE_LINK                              = 40606,
    IMSLS_N_LINKS                                  = 40607,
    IMSLS_DISPLAY_NETWORK                          = 40608,
    IMSLS_STAGE_I                                  = 40609,
    IMSLS_NO_STAGE_II                              = 40610,
    IMSLS_CLASSIFICATION                           = 40611,
    IMSLS_PREDICTED_CLASS                          = 40612,
    IMSLS_PREDICTED_CLASS_USER                     = 40613,
    IMSLS_CLASS_ERROR                              = 40614,
    IMSLS_CLASS_ERROR_USER                         = 40615,
    IMSLS_WEIGHT_INITIALIZATION_METHOD             = 40616,
    IMSLS_N_WEIGHTS                                = 40617,
    IMSLS_PREDICTED_CLASS_PROB                     = 40618,
    IMSLS_PREDICTED_CLASS_PROB_USER                = 40619,
    IMSLS_LOGISTIC_TABLE                           = 40620,
	IMSLS_FILE                                     = 40621,
    IMSLS_GET_PRINT_TYPE                           = 50000,/* Added for WIN32 support */
    IMSLS_SET_PRINT_TYPE                           = 50010,/* Added for WIN32 support */
    IMSLS_DEFAULT_PRINT_PROC                       = 50020,/* Added for WIN32 support */
    IMSLS_WINDOWS_PRINT_PROC                       = 50030,/* Added for WIN32 support */
    IMSLS_CONSOLE_PRINT_PROC                       = 50040,/* Added for WIN32 support */
    IMSLS_RETURN_STRING                            = 50100,/* Added for WIN32 support */
    IMSLS_WRITE_TO_CONSOLE                         = 50110,/* Added for WIN32 support */

    /* New for CNL v.6.0 */
    IMSLS_BEST_PERIODS                             = 50120,
    IMSLS_BEST_PERIODS_USER                        = 50130,
    IMSLS_BEST_ORDERS                              = 50140,
    IMSLS_BEST_ORDERS_USER                         = 50150,
    IMSLS_AR_ORDER                                 = 50160,
    IMSLS_D_INITIAL                                = 50170,
    IMSLS_LOG_LIKELIHOOD                           = 50180,
    IMSLS_NTIMES                                   = 50190,
    IMSLS_GRID                                     = 50200,
    IMSLS_GRID_USER                                = 50210,
    IMSLS_DELTA                                    = 50220,
    IMSLS_CRITICAL                                 = 50230,
    IMSLS_EPSILON                                  = 50240,
    IMSLS_RESIDUAL_SIGMA                           = 50250,
    IMSLS_NUM_OUTLIERS                             = 50260,
    IMSLS_OUTLIER_STAT                             = 50270,
    IMSLS_OUTLIER_STAT_USER                        = 50280,
    IMSLS_OUT_FREE_FORECAST_USER                   = 50290,
    IMSLS_OUT_FREE_FORECAST                        = 50300,
    IMSLS_OUTLIER_STATISTICS                       = 50310,
    IMSLS_OUTLIER_STATISTICS_USER                  = 50320,
    IMSLS_TAU_STATISTICS                           = 50330,
    IMSLS_TAU_STATISTICS_USER                      = 50340,
    IMSLS_OMEGA_WEIGHTS                            = 50350,
    IMSLS_OMEGA_WEIGHTS_USER                       = 50360,
    IMSLS_ARMA_PARAM                               = 50370,
    IMSLS_ARMA_PARAM_USER                          = 50380,
    IMSLS_P_INITIAL                                = 50390,
    IMSLS_Q_INITIAL                                = 50400,
    IMSLS_S_INITIAL                                = 50410,
    IMSLS_OUT_FREE_SERIES                          = 50420,
    IMSLS_OUT_FREE_SERIES_USER                     = 50430,
    IMSLS_NUM_PREDICT                              = 50440,
    IMSLS_OUTLIER_FORECAST                         = 50450,
    IMSLS_OUTLIER_FORECAST_USER                    = 50460,
    IMSLS_TIMES_ARRAY                              = 50470,
    IMSLS_TIMES_ARRAY_USER                         = 50480,
    IMSLS_DELTA_ADR                                = 50481,
    IMSLS_CRITICAL_ADR                             = 50482,
    IMSLS_EPSILON_ADR                              = 50483,

    /* */
    IMSLS_XLO                                      = 50490,
    IMSLS_XLO_USER                                 = 50491,
    IMSLS_XHI                                      = 50492,
    IMSLS_XHI_USER                                 = 50493,

    /* Naive Bayes */
    IMSLS_DISCRETE_SMOOTHING_PARM                  = 50500,
    IMSLS_CONTINUOUS_SMOOTHING_PARM                = 50501,
    IMSLS_ZERO_CORRECTION                          = 50502,
    IMSLS_SELECTED_PDF                             = 50503,
    IMSLS_USER_PDF                                 = 50504,
	IMSLS_USER_PDF_WITH_PARMS                      = 50505,
    IMSLS_GAUSSIAN_PDF                             = 50506,
    IMSLS_GAMMA_PDF                                = 50507,
    IMSLS_LOG_NORMAL_PDF                           = 50508,
    IMSLS_POISSON_PDF                              = 50509,
    IMSLS_GAUSSIAN_PARAMETERS                      = 50510,
    IMSLS_GAMMA_PARAMETERS                         = 50511,
    IMSLS_LOG_NORMAL_PARAMETERS                    = 50512,
    IMSLS_POISSON_PARAMETERS                       = 50513,
    IMSLS_COUNT_TABLE                              = 50514,
    IMSLS_COUNT_TABLE_USER                         = 50515,
    IMSLS_USE_MISSING_VALUE_PATTERNS               = 50516,
    IMSLS_IGNORE_MISSING_VALUE_PATTERNS            = 50517,
    IMSLS_NB_CLASSIFIER                            = 50518,
    IMSLS_PRINT_WARNINGS                           = 50519,

    /* genetic algorithms  */
    IMSLS_BINARY                                   = 50700,
    IMSLS_NOMINAL                                  = 50701,
    IMSLS_INTEGER                                  = 50702,
    IMSLS_REAL                                     = 50703,
    IMSLS_PARENTS                                  = 50704,
    IMSLS_WITH_REPLACEMENT                         = 50705,
    IMSLS_WITHOUT_REPLACEMENT                      = 50706,
    IMSLS_BINARY_SELECTION_PROB                    = 50707,
    IMSLS_NOMINAL_SELECTION_PROB                   = 50708,
    IMSLS_INTEGER_SELECTION_MODEL                  = 50709,
    IMSLS_REAL_SELECTION_MODEL                     = 50710,
    IMSLS_GRAY_ENCODING                            = 50711,
    IMSLS_BASE2_ENCODING                           = 50712,
    IMSLS_DOMINANCE                                = 50713,
    IMSLS_NO_DOMINANCE                             = 50714,
    IMSLS_ELITISM                                  = 50715,
    IMSLS_NO_ELITISM                               = 50716,
    IMSLS_MAX_GENERATIONS                          = 50717,
    IMSLS_LINEAR_SCALING                           = 50718,
    IMSLS_SIGMA_SCALING                            = 50719,
    IMSLS_MUTATION_PROB                            = 50720,
    IMSLS_CROSSOVER_PROB                           = 50721,
    IMSLS_CROSSOVERS                               = 50722,
    IMSLS_SELECTION_MODEL                          = 50723,
    IMSLS_FITNESS                                  = 50724,
    IMSLS_FITNESS_FCN                              = 50725,
    IMSLS_FITNESS_FCN_WITH_PARMS                   = 50726,
    IMSLS_GENERATION_STATS                         = 50727,
    IMSLS_LAST_GENERATION                          = 50728,
    IMSLS_GENERATION_GAP                           = 50729,
    IMSLS_N_GENERATIONS                            = 50730,
    IMSLS_MAX_FITNESS                              = 50731,
    IMSLS_VELOCITY                                 = 50732,
    IMSLS_ON_LINE_PERFORMANCE                      = 50733,
    IMSLS_OFF_LINE_PERFORMANCE                     = 50734,
    IMSLS_N_PARENTS                                = 50735,
    IMSLS_NO_DECODE                                = 50736,
    IMSLS_PMX_CROSSOVER                            = 50737,
    IMSLS_INVERT_CROSSOVER                         = 50738,
    IMSLS_STANDARD_CROSSOVER                       = 50739,
    IMSLS_SWAP_MUTATION                            = 50740,

    /* OpenMP options */
    IMSLS_SET_FUNCTIONS_THREAD_SAFE                = 50741,
    IMSLS_GET_FUNCTIONS_THREAD_SAFE                = 50742,

    /* Multivariate Normal     */
    IMSLS_RANDOM_SEED                              = 50600,
	/* Wilcoxon Rank Sum Test  */             
	IMSLS_EXACT_P_VALUES                           = 50610,
	IMSLS_EXACT_P_VALUES_USER                      = 50611,
	IMSLS_MANN_WHITNEY                             = 50612,
	/* one-way analysis of covariance */
	IMSLS_ADJ_ANOVA                                = 50613,
	IMSLS_ADJ_ANOVA_USER                           = 50614,
	IMSLS_PARALLEL_TESTS                           = 50615,
	IMSLS_PARALLEL_TESTS_USER                      = 50616,
	IMSLS_R_MATRIX                                 = 50617,
	IMSLS_R_MATRIX_USER                            = 50618,
	IMSLS_XYMEAN                                   = 50619,
	IMSLS_XYMEAN_USER                              = 50620,
	IMSLS_COV_MEANS                                = 50621,
	IMSLS_COV_MEANS_USER                           = 50622,
	IMSLS_COV_COEF                                 = 50623,
	IMSLS_COV_COEF_USER                            = 50624,
	IMSLS_REG_ANOVA                                = 50625,
	IMSLS_REG_ANOVA_USER                           = 50626,
	IMSLS_COEF_TABLES                              = 50627,
	IMSLS_COEF_TABLES_USER                         = 50628,

	/* Normal Random Variables */
    IMSLS_ZIGGURAT_METHOD                          = 50629,

	/* Information Criteria */
	IMSLS_AICC                                     = 50630,
	IMSLS_BIC                                      = 50631,
	IMSLS_MODEL_SELECTION_CRITERION                = 50632,

        IMSLS_GET_OMP_MKL_SINGLE_THREAD                = 50633,
        IMSLS_SET_OMP_MKL_SINGLE_THREAD                = 50634,

    IMSLS_LAST_FLAG                                = 99999
};
#endif /* IMSLS_H */
