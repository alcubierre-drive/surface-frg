#include "surface_greens_function.h"
#include <diverge_Eigen3.hpp>
#include <cassert>

/*
 * this subroutine is used to calculate surface state using
 * green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858
 * Using eq.(11) and eq.(13) (16), (17)
 * History:
 *         by Quan Sheng Wu on Oct/17/2012
 *         C++/C workspace version by Lennart Klebl on Dec/13/2025
 *+---------+---------+---------+---------+---------+---------+--------+
 */

static inline void invert_2x2( complex128_t* out, const complex128_t* in ) {
    const complex128_t invdet = 1./(in[IDX2(0,0,2)]*in[IDX2(1,1,2)] - in[IDX2(1,0,2)]*in[IDX2(0,1,2)]);
    out[IDX2(0,0,2)] = invdet * in[IDX2(1,1,2)];
    out[IDX2(0,1,2)] = -invdet * in[IDX2(0,1,2)];
    out[IDX2(1,0,2)] = -invdet * in[IDX2(1,0,2)];
    out[IDX2(1,1,2)] = invdet * in[IDX2(0,0,2)];
}
static inline void invert_nxn( complex128_t* out, const complex128_t* in, const uint64_t n ) {
    Map<cMatXcd>(out, n,n) = Map<const cMatXcd>(in, n,n).inverse();
}
static inline void invert_any( complex128_t* out, const complex128_t* in, const uint64_t n ) {
    switch (n) {
        case 0: break;
        case 1: out[0] = 1./in[0]; break;
        case 2: invert_2x2(out, in); break;
        default: invert_nxn(out, in, n); break;
    }
}

void surface_gf( const complex128_t omegac, complex128_t* GLL,
        complex128_t* GRR, complex128_t* GB, const complex128_t* H00_,
        const complex128_t* H01_, const index_t dim, const int itermax/*=100*/,
        const double accuracy/*=1.e-16*/, complex128_t* workspace/*size: POW2(dim)*9*/ ) {

    Map<const cMatXcd> H00(H00_, dim, dim);
    Map<const cMatXcd> H01(H01_, dim, dim);

    int idx = 0;
    Map<cMatXcd> ones(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> alphai(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> betai(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> epsiloni(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> epsilons(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> epsilons_t(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> mat1(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> mat2(workspace + (idx++) * POW2(dim), dim, dim);
    Map<cMatXcd> g0(workspace + (idx++) * POW2(dim), dim, dim);
    assert( idx == 9 );

    for (index_t i=0; i<dim; ++i)
    for (index_t j=0; j<dim; ++j)
        ones(i,j) = i==j ? 1.0 : 0.0;

    epsiloni = H00;
    epsilons = H00;
    epsilons_t = H00;
    alphai = H01;
    betai = H01.adjoint();

    for (int iter=0; iter<itermax; ++iter) {
        mat1 = omegac * ones - epsiloni;
        invert_any( g0.data(), mat1.data(), dim );
        mat1 = alphai * g0;
        mat2 = betai * g0;
        g0 = mat1 * betai;
        epsiloni = epsiloni + g0;
        epsilons = epsilons + g0;
        g0 = mat2 * alphai;
        epsiloni = epsiloni + g0;
        epsilons_t = epsilons_t + g0;
        g0 = mat1 * alphai;
        alphai= g0;
        g0 = mat2 * betai;
        betai = g0;
        if (alphai.array().abs().maxCoeff() < accuracy) break;
    }
    if (GLL) { mat1 = omegac*ones - epsilons;
               invert_any( GLL, mat1.data(), dim ); }
    if (GRR) { mat1 = omegac*ones - epsilons_t;
               invert_any( GRR, mat1.data(), dim ); }
    if (GB) { mat1 = omegac*ones - epsiloni;
              invert_any( GB, mat1.data(), dim ); }
}
