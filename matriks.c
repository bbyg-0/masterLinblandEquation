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


__attribute__((target("avx2,fma")))
void matriks_komutator_avx2(const float *__restrict A, const float *__restrict B, float *__restrict Hasil) {
    
    // Alokasi memori lokal untuk menyimpan hasil sementara AB (Matrix A x B)
    // Dipaksa selaras 32-byte agar load/store AVX2 berjalan super cepat
    float AB[N * N] __attribute__((aligned(32)));

    // LANGKAH 1: Hitung AB = A x B (Diparalelkan dengan OpenMP)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j += 8) {
            __m256 reg_Hasil = _mm256_setzero_ps();
            for (int k = 0; k < N; k++) {
                __m256 reg_A = _mm256_set1_ps(A[i * N + k]);
                __m256 reg_B = _mm256_loadu_ps(&B[k * N + j]);
                reg_Hasil = _mm256_fmadd_ps(reg_A, reg_B, reg_Hasil);
            }
            // Karena AB sudah ter-align 32-byte, kita bisa pakai store biasa (bukan loadu/storeu)
            _mm256_store_ps(&AB[i * N + j], reg_Hasil);
        }
    }

    // LANGKAH 2: Hitung BA = B x A dan LANGSUNG kurangkan dari AB (Loop Fusing)
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j += 8) {
            __m256 reg_BA = _mm256_setzero_ps();
            for (int k = 0; k < N; k++) {
                __m256 reg_B = _mm256_set1_ps(B[i * N + k]);
                __m256 reg_A = _mm256_loadu_ps(&A[k * N + j]);
                reg_BA = _mm256_fmadd_ps(reg_B, reg_A, reg_BA);
            }

            // OPTIMASI UTAMA: Load hasil AB dari RAM langsung ke register
            __m256 reg_AB = _mm256_load_ps(&AB[i * N + j]);

            // Lakukan pengurangan langsung di dalam register CPU: [A, B] = AB - BA
            __m256 reg_Komutator = _mm256_sub_ps(reg_AB, reg_BA);

            // Simpan hasil akhir ke matriks Hasil
            _mm256_storeu_ps(&Hasil[i * N + j], reg_Komutator);
        }
    }
}

__attribute__((target("avx2,fma")))
void matriks_anticomutator_avx2(const float *__restrict A, const float *__restrict B, float *__restrict Hasil) {
    
    float AB[N * N] __attribute__((aligned(32)));

    // LANGKAH 1: Hitung AB = A x B
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j += 8) {
            __m256 reg_Hasil = _mm256_setzero_ps();
            for (int k = 0; k < N; k++) {
                __m256 reg_A = _mm256_set1_ps(A[i * N + k]);
                __m256 reg_B = _mm256_loadu_ps(&B[k * N + j]);
                reg_Hasil = _mm256_fmadd_ps(reg_A, reg_B, reg_Hasil);
            }
            _mm256_store_ps(&AB[i * N + j], reg_Hasil);
        }
    }

    // LANGKAH 2: Hitung BA = B x A dan LANGSUNG jumlahkan dengan AB
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j += 8) {
            __m256 reg_BA = _mm256_setzero_ps();
            for (int k = 0; k < N; k++) {
                __m256 reg_B = _mm256_set1_ps(B[i * N + k]);
                __m256 reg_A = _mm256_loadu_ps(&A[k * N + j]);
                reg_BA = _mm256_fmadd_ps(reg_B, reg_A, reg_BA);
            }

            __m256 reg_AB = _mm256_load_ps(&AB[i * N + j]);

            // BEDANYA DI SINI: Anti-Komutator menggunakan operasi PENJUMLAHAN (AB + BA)
            __m256 reg_AntiKomutator = _mm256_add_ps(reg_AB, reg_BA);

            _mm256_storeu_ps(&Hasil[i * N + j], reg_AntiKomutator);
        }
    }
}
