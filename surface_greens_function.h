/*
 * this subroutine is used to calculate surface state using
 * green's function method  ---  J.Phys.F.Met.Phys.15(1985)851-858
 * Using eq.(11) and eq.(13) (16), (17)
 * History:
 *         by Quan Sheng Wu on Oct/17/2012
 *         C++/C workspace version by Lennart Klebl on Dec/13/2025
 *+---------+---------+---------+---------+---------+---------+--------+
 */

#pragma once

#include <diverge.h>

#ifdef __cplusplus
extern "C" {
#endif

void surface_gf( const complex128_t omega, complex128_t* GLL_,
        complex128_t* GRR_, complex128_t* GB_, const complex128_t* H00_,
        const complex128_t* H01_, const index_t dim, const int itermax/*=100*/,
        const double accuracy/*=1.e-16*/, complex128_t* workspace/*size: POW2(dim)*9*/ );

#ifdef __cplusplus
}
#endif
