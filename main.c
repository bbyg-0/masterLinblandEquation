#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include <immintrin.h>

#include <ctype.h>

#include "matriks.h"

#define MAX_GANGGUAN 100

// Struktur data untuk menyimpan konfigurasi setiap gangguan yang diparsing
typedef struct {
    int qubit_idx;      // Qubit ke berapa (0 to n)
    char jenis_gangguan;// Jenis operator ('a', 'b', 'c', dst)
    float nilai_gamma;  // Nilai float untuk gamma
} GangguanKuantum;

void cetak_panduan(char *nama_program) {
    printf("Penggunaan: %s [FLAGS]\n", nama_program);
    printf("Flags yang tersedia:\n");
    printf("  -q, --qubit     Jumlah qubit (Default: 1)\n");
    printf("  -w, --omega     Nilai kekuatan kopling Omega (Default: 1.0)\n");
    printf("  -e, --energy    Nilai energi Detuning E (Default: 0.0)\n");
    printf("  -d, --disturb   Detail gangguan dengan format [Qubit][Jenis][Gamma]\n");
    printf("                  Contoh: -d 0a0.15 -d 1b0.05\n");
    printf("                  (a: Emission, b: Dephasing, c: Pumping)\n");
    printf("  -h, --help      Menampilkan panduan ini\n");
}

int main(int argc, char *argv[]) {
    // 1. Inisialisasi Nilai Default
    int jumlah_qubit = 1;
    float omega = 1.0f;
    float detuning_E = 0.0f;

    // Array untuk menampung daftar gangguan dari user
    GangguanKuantum daftar_gangguan[MAX_GANGGUAN];
    int total_gangguan = 0;

    // 2. Definisi Flags Panjang
    static struct option long_options[] = {
        {"qubit",   required_argument, 0, 'q'},
        {"omega",   required_argument, 0, 'w'},
        {"energy",  required_argument, 0, 'e'},
        {"disturb", required_argument, 0, 'd'},
        {"help",    no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;

    // 3. Loop Membaca Flags
    // Menambahkan "d:" pada optstring untuk menerima argumen gangguan
    while ((opt = getopt_long(argc, argv, "q:w:e:d:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'q':
                jumlah_qubit = atoi(optarg);
                break;
            case 'w':
                omega = atof(optarg);
                break;
            case 'e':
                detuning_E = atof(optarg);
                break;
	    case 'd':
                if (total_gangguan < MAX_GANGGUAN) {
                    char *input = optarg;
                    
                    // TRANSLASI KRUSIAL:
                    // Ambil angka urutan manusia (1-9), lalu kurangi 1 
                    // agar menjadi indeks komputer (0-8)
                    int urutan_manusia = input[0] - '0'; 
                    daftar_gangguan[total_gangguan].qubit_idx = urutan_manusia - 1;
                    
                    // Karakter ke-1 dan seterusnya tetap sama
                    daftar_gangguan[total_gangguan].jenis_gangguan = input[1];
                    daftar_gangguan[total_gangguan].nilai_gamma = atof(&input[2]);
                    
                    total_gangguan++;
                }
                break;
            case 'h':
                cetak_panduan(argv[0]);
                return 0;
            default:
                cetak_panduan(argv[0]);
                return 1;
        }
    }

    // 4. Validasi Hasil Input Gangguan terhadap Jumlah Qubit
    for (int i = 0; i < total_gangguan; i++) {
        // Jika indeks komputer kurang dari 0 atau melebihi (jumlah_qubit - 1)
        if (daftar_gangguan[i].qubit_idx >= jumlah_qubit || daftar_gangguan[i].qubit_idx < 0) {
            fprintf(stderr, "Error: Anda ingin mengganggu Qubit ke-%d, tetapi Anda hanya mensimulasikan %d qubit.\n", 
                    daftar_gangguan[i].qubit_idx + 1, jumlah_qubit); // Ditambah 1 saat dicetak ke user
            return 1;
        }
    }

    // 5. Tampilkan Parameter yang Berhasil Diterima
    printf("=== PARAMETER SIMULASI KUANTUM ===\n");
    printf("Jumlah Qubit : %d\n", jumlah_qubit);
    printf("Nilai Omega  : %.4f\n", omega);
    printf("Nilai E      : %.4f\n", detuning_E);
    printf("Total Efek Lingkungan Aktif: %d\n", total_gangguan);
    
    for (int i = 0; i < total_gangguan; i++) {
        char *nama_efek = "Tidak Diketahui";
        if (daftar_gangguan[i].jenis_gangguan == 'a') nama_efek = "Spontaneous Emission (Sigma-)";
        if (daftar_gangguan[i].jenis_gangguan == 'b') nama_efek = "Dephasing (Sigma_Z)";
        if (daftar_gangguan[i].jenis_gangguan == 'c') nama_efek = "Incoherent Pumping (Sigma+)";

        printf("  -> Kebocoran pada Qubit [%d] | Jenis: %s | Gamma: %.4f\n", 
               daftar_gangguan[i].qubit_idx, nama_efek, daftar_gangguan[i].nilai_gamma);
    }

    float A[N * N];
    float B[N * N];
    float Hasil[N * N];

    // Mengisi data dummy untuk matriks A dan B
    for (int i = 0; i < N * N; i++) {
        A[i] = 1.0f; // Semua elemen A bernilai 1
        B[i] = 2.0f; // Semua elemen B bernilai 2
        Hasil[i] = 0.0f;
    }

    printf("Menghitung perkalian matriks 8x8 dengan AVX2 + FMA...\n\n");
    
    // Jalankan fungsi optimasi hardware
    matriks_perkalikan_avx2(A, B, Hasil);

    // Cetak sebagian hasil (Baris pertama saja untuk pembuktian)
    printf("Hasil Baris Pertama (Indeks 0-7):\n");
    for (int j = 0; j < N; j++) {
        printf("%.1f ", Hasil[0 * N + j]);
    }
    printf("\n\n(Ekspektasi jika semua A=1 dan B=2 dikali matriks 8x8 adalah 1*2*8 = 16.0)\n");
    return 0;
}
