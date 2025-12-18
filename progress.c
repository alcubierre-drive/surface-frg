/* Copyright (C) 2025 Lennart Klebl; GPLv3+ (LICENSE.md) */

#include "progress.h"
#include <stdatomic.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MIN(x,y) ((x) < (y) ? (x) : (y))

struct progress_bar_t {
    atomic_int progress;
    atomic_int ntotal;

    atomic_int width;
    atomic_char pchar;
    atomic_char echar;
    atomic_int use_pchar_x;
    atomic_int use_echar_x;

    // unicode support
    char pchar_x[32];
    char echar_x[32];

    char prefix[32];
    atomic_flag has_finished;
};

progress_bar_t* progress_bar_init( int ntotal ) {
    progress_bar_t* p = calloc(1, sizeof *p);
    atomic_init( &p->progress, 0 );
    atomic_init( &p->ntotal, ntotal );
    atomic_init( &p->width, 60 );
    atomic_init( &p->pchar, '#' );
    atomic_init( &p->echar, '-' );
    atomic_init( &p->use_pchar_x, 0 );
    atomic_init( &p->use_echar_x, 0 );

    p->has_finished = (atomic_flag)ATOMIC_FLAG_INIT;
    atomic_flag_clear( &p->has_finished );
    return p;
}

void progress_bar_free( progress_bar_t* p ) {
    free( p );
}

void progress_bar_set_prefix( progress_bar_t* p, const char* prefix ) {
    strncpy( p->prefix, prefix, 31 );
}

void progress_bar_set_width( progress_bar_t* p, int width ) {
    atomic_store( &p->width, width );
}

void progress_bar_set_pchar( progress_bar_t* p, char pchar ) {
    atomic_store( &p->pchar, pchar );
}

void progress_bar_set_echar( progress_bar_t* p, char echar ) {
    atomic_store( &p->echar, echar );
}

void progress_bar_set_pchar_x( progress_bar_t* p, const char* pchar ) {
    if (!pchar) {
        atomic_store( &p->use_pchar_x, 0 );
    } else {
        size_t ncopy = MIN( strlen(pchar), sizeof(p->pchar_x)-1 );
        memcpy( p->pchar_x, pchar, ncopy );
        atomic_store( &p->use_pchar_x, 1 );
    }
}

void progress_bar_set_echar_x( progress_bar_t* p, const char* echar ) {
    if (!echar) {
        atomic_store( &p->use_echar_x, 0 );
    } else {
        size_t ncopy = MIN( strlen(echar), sizeof(p->echar_x)-1 );
        memcpy( p->echar_x, echar, ncopy );
        atomic_store( &p->use_echar_x, 1 );
    }
}

static void progress_bar_print( progress_bar_t* p ) {
    float percentage = (float)atomic_load( &p->progress ) / (float)( atomic_load( &p->ntotal ) );
    int np = (int)( percentage * (float)atomic_load(&p->width) + 0.5 );
    char pc = atomic_load( &p->pchar );
    char ec = atomic_load( &p->echar );

    int use_pc = !atomic_load( &p->use_pchar_x ),
        use_ec = !atomic_load( &p->use_echar_x );

    if (np < 0) np = 0;
    int w = atomic_load( &p->width );
    if (np > w) np = w;

    int has_finished = 0;
    if (np == w)
        has_finished = atomic_flag_test_and_set( &p->has_finished );

    char pstr[p->width*8 + 64];
    pstr[0] = '\r';
    pstr[1] = '\0';
    strcat(pstr, p->prefix);
    int c = strlen(pstr);
    if (c > 2) pstr[c++] = ' ';
    pstr[c++] = '[';
    int i = 0;
    if (use_pc) {
        for (; i<np; ++i) pstr[c++] = pc;
    } else {
        for (; i<np; ++i) {
            strcpy(pstr + c, p->pchar_x);
            c += strlen(p->pchar_x);
        }
    }

    if (use_ec) {
        for (; i<w; ++i) pstr[c++] = ec;
    } else {
        for (; i<w; ++i) {
            strcpy(pstr + c, p->echar_x);
            c += strlen(p->echar_x);
        }
    }
    char end[32] = {0};
    sprintf( end, "] %3.0f", percentage * 100.0 );
    for (const char* pe = end; *pe; ++pe)
        pstr[c++] = *pe;
    pstr[c++] = '%';
    pstr[c++] = '\0';

    if (!has_finished)
        fputs( pstr, stderr );
    if (!has_finished && (np == w))
        fputs( "\n", stderr );
    fflush( stderr );
}

void progress_bar_add( progress_bar_t* p ) {
    atomic_fetch_add( &p->progress, 1 );
    progress_bar_print( p );
}

void progress_bar_add_N( progress_bar_t* p, int N ) {
    atomic_fetch_add( &p->progress, N );
    progress_bar_print( p );
}

void progress_bar_resize( progress_bar_t* p, int ntotal ) {
    atomic_store( &p->ntotal, ntotal );
    progress_bar_print( p );
}

int progress_bar_idx( progress_bar_t* p ) {
    return atomic_load( &p->progress );
}

int progress_bar_size( progress_bar_t* p ) {
    return atomic_load( &p->ntotal );
}
