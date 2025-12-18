/* Copyright (C) 2025 Lennart Klebl; GPLv3+ (LICENSE.md) */

#include <diverge.h>
#include <diverge_misc.h>
LINALG_IMPLEMENT

#include <unistd.h>

#include <omp.h>
#include "surface_greens_function.h"
#include "progress.h"

typedef struct {
    double v, w, tp, U, mu;
    int nk, nkf;
    double tu_dist;

    int n_omega;
    double omega_min, omega_max;
    double omega_eta;

    char model_file[1024];
    char post_file[1024];
    char flow_file[1024];
} params_t;

static const params_t params_default = {
    .v = 0.0, .w = 0.0, .tp = 0.0, .U = 3.0, .mu = 0.0,
    .nk = 24, .nkf = 5,
    .tu_dist = 2.01,

    .n_omega = 100,
    .omega_min = -2., .omega_max = 2.,
    .omega_eta = 1.e-1,

    .model_file = "ssh_surface_mod.dvg",
    .post_file = "ssh_surface_post.dvg",
    .flow_file = "ssh_surface_flow.dat",
};

static params_t parse_args( int* pargc, char*** pargv );
static diverge_model_t* hubbard_ssh_surface( params_t p );

int main( int argc, char** argv ) {
    diverge_init( &argc, &argv );
    atexit( &diverge_finalize );

    params_t p = parse_args( &argc, &argv );
    diverge_model_t* m = hubbard_ssh_surface( p );
    diverge_model_output_set_npath(1);
    if (strcmp(p.model_file, "")) diverge_model_to_file(m, p.model_file);

    if (p.tu_dist >= 0.0) diverge_model_internals_tu(m, p.tu_dist);
    else                  diverge_model_internals_grid(m);
    diverge_flow_step_t* s = diverge_flow_step_init_any(m, "PCD");

    FILE* flow = NULL;
    if (diverge_mpi_comm_rank() == 0 && strcmp(p.flow_file, ""))
        flow = fopen(p.flow_file, "w");

    #define bprintf( ... ) { \
        if (flow) { \
            fprintf(flow, __VA_ARGS__); \
            fflush(flow); \
        } \
        mpi_printf(__VA_ARGS__); \
        fflush(stdout); \
    }

    diverge_euler_t eu = diverge_euler_defaults; eu.maxvert = 30;
    double maxvert = 0.0;
    do {
        diverge_flow_step_euler( s, eu.Lambda, eu.dLambda );
        double chmax[3] = {0};
        diverge_flow_step_chanmax( s, chmax );
        for (int i=0; i<3; ++i) chmax[i] = fabs(chmax[i]);
        bprintf( "%.5e %.5e %.5e %.5e\n", eu.Lambda, chmax[0], chmax[1], chmax[2] );
        maxvert = MAX(chmax[0],MAX(chmax[1],chmax[2]));
    } while (diverge_euler_next( &eu, maxvert ));

    if (flow) {
        fclose(flow);
        flow = NULL;
    }

    diverge_postprocess_conf_t pc = diverge_postprocess_conf_defaults;
    pc.tu_symmetry_maps = 1;
    pc.tu_storing_threshold = 0.9;
    pc.tu_storing_relative = 1;
    pc.tu_channel_calc_project = -1;
    pc.tu_which_solver_mode = 'B';
    pc.tu_susceptibilities_ff = 1;

    if (strcmp(p.post_file, "")) diverge_postprocess_and_write_fg(s, p.post_file, &pc);

    diverge_flow_step_free( s );
    diverge_model_free( m );

    return EXIT_SUCCESS;
}

typedef struct {
    double v;
    double w;
    double tp;
    double mu;

    index_t nthr;
    complex128_t *GLL;
    complex128_t *H00, *H01;
    complex128_t *workspace;
} hubbard_ssh_surface_params_t;

static hubbard_ssh_surface_params_t* diverge_model_get_hubbard_ssh_surface_params( diverge_model_t* m ) {
    index_t ndouble = sizeof(hubbard_ssh_surface_params_t) / sizeof(double) + 1;
    double* pos_ptr = m->positions[0];
    return (hubbard_ssh_surface_params_t*)(pos_ptr + 3*MAX_N_ORBS - ndouble - 1);
}

static void hubbard_ssh_surface_ham(const diverge_model_t* m, complex128_t* buf) {
    (void)m;
    (void)buf;
}

static greensfunc_op_t hubbard_ssh_surface_greens(const diverge_model_t* m,
        complex128_t Lambda, gf_complex_t* buf) {
    double* dLambda = (double*)&Lambda;
    const int skip_minus = qnan_isnan( *dLambda );
    complex64_t* xLambda = (complex64_t*)(dLambda+1);
    if (skip_minus) Lambda = *xLambda;

    hubbard_ssh_surface_params_t* p = diverge_model_get_hubbard_ssh_surface_params((diverge_model_t*)m);
    const index_t nktot = kdimtot(m->nk,m->nkf);
    const double* kfmesh = m->internals->kfmesh;
    #pragma omp parallel for num_threads(p->nthr)
    for (index_t k=0; k<nktot; ++k) {
        index_t t = omp_get_thread_num();
        complex128_t *GLL = p->GLL + 4*t,
                     *H00 = p->H00 + 4*t,
                     *H01 = p->H01 + 4*t,
                     *workspace = p->workspace + 4*9*t;
        double kx = kfmesh[3*k+0], ky = kfmesh[3*k+1];
        H00[0] = H00[3] = -2.*(cos(kx)+cos(ky))-4.*p->tp*cos(kx)*cos(ky);
        H00[1] = H00[2] = p->v;
        H01[2] = p->w;
        surface_gf( Lambda + p->mu, GLL, NULL, NULL, H00, H01, 2,
                20/*itermax*/, 1.e-6/*accuracy*/, workspace );
        buf[k] = GLL[0];
        if (!skip_minus) {
            surface_gf(-Lambda + p->mu, GLL, NULL, NULL, H00, H01, 2,
                20/*itermax*/, 1.e-6/*accuracy*/, workspace );
            buf[k+nktot] = GLL[0];
        }
    }
    return greensfunc_op_cpu;
}

static claM3d_t rotation( double theta ) {
    claM3d_t M = {0};
    M.m[0][0] = M.m[1][1] = cos(theta);
    M.m[0][1] = -(M.m[1][0] = sin(theta));
    M.m[2][2] = 1.0;
    return M;
}

static diverge_model_t* hubbard_ssh_surface( params_t p ) {
    diverge_model_t* m = diverge_model_init();
    m->lattice[0][0] = m->lattice[1][1] = m->lattice[2][2] = 1.0;
    m->n_spin = m->n_orb = m->SU2 = 1;

    m->nk[0] = m->nk[1] = p.nk;
    m->nkf[0] = m->nkf[1] = p.nkf;

    hubbard_ssh_surface_params_t* gdata = diverge_model_get_hubbard_ssh_surface_params(m);
    gdata->v = p.v;
    gdata->w = p.w;
    gdata->tp = p.tp;
    gdata->mu = p.mu;

    index_t nthr = gdata->nthr = diverge_omp_num_threads();

    p.n_omega = MAX(p.n_omega,0);
    index_t data_size = 2 * p.n_omega * sizeof(complex128_t) +
        nthr * POW2(2) * sizeof(complex128_t) * 14;
    char* data_ptr = calloc(data_size, sizeof*data_ptr);
    m->data = data_ptr;
    m->nbytes_data = 2 * p.n_omega * sizeof(complex128_t); // do not save the workspace
    if (m->nbytes_data == 0) m->nbytes_data = 1;

    { complex128_t* workspace_ptr = (complex128_t*)(data_ptr + 2 * p.n_omega * sizeof(complex128_t));

    gdata->GLL = workspace_ptr; workspace_ptr += nthr * POW2(2);
    gdata->H00 = workspace_ptr; workspace_ptr += nthr * POW2(2);
    gdata->H01 = workspace_ptr; workspace_ptr += nthr * POW2(2);
    gdata->workspace = workspace_ptr; } // scope workspace_ptr ends

    m->gfill = &hubbard_ssh_surface_greens;
    m->hfill = &hubbard_ssh_surface_ham;

    m->vert = calloc(1, sizeof*m->vert);
    m->vert[m->n_vert++] = (rs_vertex_t){ .chan='D', .V=p.U };

    m->orb_symmetries = calloc(8, sizeof*m->orb_symmetries);
    for (index_t i=0; i<8; ++i) {
        m->orb_symmetries[i] = 1.0;
        if (i < 4) {
            claM3d_t R = rotation(M_PI/2. * i);
            claM3d_rmap( R, m->rs_symmetries[i][0] );
        } else {
            claM3d_t R = rotation(M_PI/4. * (i-4));
            claM3d_t M = {{{-1,0,0},{0,1,0},{0,0,1}}};
            claM3d_rmap( claM3d_matmul(claM3d_matmul(R,M),claM3d_transpose(R)), m->rs_symmetries[i][0] );
        }
    }
    m->n_sym = 8;

    diverge_model_internals_common(m);

    complex128_t* omega = m->data;
    complex128_t* dos = omega + p.n_omega;
    double* omega_cpy = diverge_linspace( p.omega_min, p.omega_max, p.n_omega );

    progress_bar_t* pr = progress_bar_init( p.n_omega );
    progress_bar_set_prefix( pr, "𝒜(ω)" );
    progress_bar_set_width( pr, 50 );
    progress_bar_set_pchar_x( pr, "𝒜" );
    progress_bar_set_echar_x( pr, "ω" );
    // progress_bar_set_pchar( pr, 'A' );
    // progress_bar_set_echar( pr, 'w' );
    for (index_t i=0; i<p.n_omega; ++i) {
        omega[i] = omega_cpy[i] - I*p.omega_eta;
        index_t nktot = kdimtot(m->nk,m->nkf);
        gf_complex_t* gbuf = m->internals->greens;

        complex64_t Lambda = omega[i];
        double real = qnan_gen(0);
        complex128_t IN = 0;
        memcpy( &IN, &real, sizeof(real) );
        memcpy( (char*)&IN + sizeof(double), &Lambda, sizeof(Lambda) );
        (*m->gfill)( m, IN, gbuf );
        for (index_t k=0; k<nktot; ++k) dos[i] += gbuf[k];
        progress_bar_add( pr );
    }
    progress_bar_free( pr );
    return m;
}

static params_t parse_args( int* pargc, char*** pargv ) {
    params_t p = params_default;

    char pad[512] = "";
    sprintf( pad, "%*s", (int)strlen(*pargv[0]), " " );
    int opt;
    while ((opt = getopt(*pargc, *pargv, "K:k: U: v:w: m: t: o:O:n:e: M:F:P: T:h")) != -1) {
        switch (opt) {
            case 'K': p.nk = atoi(optarg); break;
            case 'k': p.nkf = atoi(optarg); break;

            case 'U': p.U = atof(optarg); break;

            case 'v': p.v = atof(optarg); break;
            case 'w': p.w = atof(optarg); break;
            case 'm': p.mu = atof(optarg); break;
            case 't': p.tp = atof(optarg); break;

            case 'o': p.omega_min = atof(optarg); break;
            case 'O': p.omega_max = atof(optarg); break;
            case 'n': p.n_omega = atoi(optarg); break;
            case 'e': p.omega_eta = atof(optarg); break;

            case 'M': strncpy( p.model_file, optarg, sizeof(p.model_file)-1 ); break;
            case 'P': strncpy( p.post_file, optarg, sizeof(p.post_file)-1 ); break;
            case 'F': strncpy( p.flow_file, optarg, sizeof(p.flow_file)-1 ); break;

            case 'T': p.tu_dist = atof(optarg); break;

            case 'h': fprintf( stderr,
                        "Usage: %s [-K:nk] [-k:nkf] [-U:Hubbard-U] [-v:hopping v] [-w:hopping w] [-m:mu]\n"
                        "       %s [-t:tprime] [-o:omega_min] [-O:omega_max] [-n:n_omega] [-e:omega_eta]\n"
                        "       %s [-M:model file] [-P:post file] [-F:flow file] [-T:tu_dist] [-h this help]\n"
                        , *pargv[0], pad, pad );
                      exit(EXIT_SUCCESS); break;
            case '?':
            default:
                      fprintf( stderr, "%s -h for help.\n", *pargv[0] );
                      exit(EXIT_FAILURE); break;
        }
    }
    *pargc -= optind;
    *pargv += optind;

    return p;
}

