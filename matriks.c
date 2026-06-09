#include "matriks.h"

#include <immintrin.h>
#include <omp.h>
#include <stdio.h>

void matriks_perkalikan_avx2(const float *__restrict A, const float *__restrict B, float *__restrict Hasil) {
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j += 8) {
            
            __m256 reg_Hasil = _mm256_setzero_ps();

            for (int k = 0; k < N; k++) {
                __m256 reg_A = _mm256_set1_ps(A[i * N + k]);
                
                __m256 reg_B = _mm256_loadu_ps(&B[k * N + j]); 
                
                reg_Hasil = _mm256_fmadd_ps(reg_A, reg_B, reg_Hasil);
            }

            _mm256_storeu_ps(&Hasil[i * N + j], reg_Hasil);
        }
    }
}
