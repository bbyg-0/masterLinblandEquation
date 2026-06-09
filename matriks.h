#ifndef MATRIKS_H
#define MATRIKS_H

#include <immintrin.h>
#include <stdio.h>

#define N 16

void matriks_perkalikan_avx2(const float *__restrict A, const float *__restrict B, float *__restrict Hasil);

void matriks_komutator_avx2(const float *__restrict A, const float *__restrict B, float *__restrict Hasil) {

void matriks_anticomutator_avx2(const float *__restrict A, const float *__restrict B, float *__restrict Hasil) {


#endif
