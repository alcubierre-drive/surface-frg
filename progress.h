/* Copyright (C) 2026 Lennart Klebl; GPLv3+ (LICENSE.md) */

#pragma once

typedef struct progress_bar_t progress_bar_t;

progress_bar_t* progress_bar_init( int ntotal );
void progress_bar_free( progress_bar_t* p );

void progress_bar_set_prefix( progress_bar_t* p, const char* prefix );
void progress_bar_set_width( progress_bar_t* p, int width );
void progress_bar_set_pchar( progress_bar_t* p, char pchar );
void progress_bar_set_echar( progress_bar_t* p, char echar );
void progress_bar_set_pchar_x( progress_bar_t* p, const char* pchar );
void progress_bar_set_echar_x( progress_bar_t* p, const char* echar );

void progress_bar_add( progress_bar_t* p );
void progress_bar_add_N( progress_bar_t* p, int N );

void progress_bar_resize( progress_bar_t* p, int ntotal );

int progress_bar_idx( progress_bar_t* p );
int progress_bar_size( progress_bar_t* p );

