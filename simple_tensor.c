#include "simple_tensor.h"

/* Pauli matrices - only keep I2 since X and Z are not used directly */
static const Complex I2[4] = {{1,0}, {0,0}, {0,0}, {1,0}};

/* Single qubit Hamiltonian: H = ω*Z + E*X */
Complex* single_qubit_hamiltonian(float omega, float E) {
    Complex* H = (Complex*)malloc(4 * sizeof(Complex));
    
    if (!H) return NULL;
    
    /* H = ω·σ_z + E·σ_x */
    H[0].real = omega;  H[0].imag = 0.0f;  /* ⟨0|H|0⟩ = ω */
    H[1].real = E;      H[1].imag = 0.0f;  /* ⟨0|H|1⟩ = E */
    H[2].real = E;      H[2].imag = 0.0f;  /* ⟨1|H|0⟩ = E */
    H[3].real = -omega; H[3].imag = 0.0f;  /* ⟨1|H|1⟩ = -ω */
    
    return H;
}

/* Rest of the functions remain the same */
/* tensor_product, create_hamiltonian, create_full_density_matrix, 
   print_matrix, free_matrix */

/* Tensor product implementation */
void tensor_product(Complex* A, int N_A, Complex* B, int N_B, Complex* C) {
    int N_C = N_A * N_B;
    
    /* Clear C first */
    memset(C, 0, N_C * N_C * sizeof(Complex));
    
    for (int i = 0; i < N_A; i++) {
        for (int j = 0; j < N_A; j++) {
            Complex a = A[i * N_A + j];
            
            if (a.real == 0 && a.imag == 0) continue;
            
            for (int k = 0; k < N_B; k++) {
                for (int l = 0; l < N_B; l++) {
                    Complex b = B[k * N_B + l];
                    int idx = (i * N_B + k) * N_C + (j * N_B + l);
                    
                    C[idx].real = a.real * b.real - a.imag * b.imag;
                    C[idx].imag = a.real * b.imag + a.imag * b.real;
                }
            }
        }
    }
}

/* Create N-qubit Hamiltonian */
Complex* create_hamiltonian(float* omegas, float* Es, int n_qubits) {
    int dim = 1 << n_qubits;
    int total_size = dim * dim;

    Complex* H_total = calloc(total_size, sizeof(Complex));
    if (!H_total) return NULL;

    // Inline identity — no dependency on global I2
    Complex ident[4] = {{1,0},{0,0},{0,0},{1,0}};

    for (int q = 0; q < n_qubits; q++) {
        Complex* H_q = single_qubit_hamiltonian(omegas[q], Es[q]);
        if (!H_q) { free(H_total); return NULL; }

        Complex* term = malloc(4 * sizeof(Complex));
        if (!term) { free(H_q); free(H_total); return NULL; }
        memcpy(term, (0 == q) ? H_q : ident, 4 * sizeof(Complex));

        for (int i = 1; i < n_qubits; i++) {
            int cur  = 1 << i;           // side length of `term`
            int new_total = cur * 2 * cur * 2;

            Complex* new_term = calloc(new_total, sizeof(Complex));
            if (!new_term) { free(term); free(H_q); free(H_total); return NULL; }

            Complex* mat2 = (i == q) ? H_q : ident;
            tensor_product(term, cur, mat2, 2, new_term);

            free(term);
            term = new_term;
        }

        for (int i = 0; i < total_size; i++) {
            H_total[i].real += term[i].real;
            H_total[i].imag += term[i].imag;
        }

        free(term);
        free(H_q);
    }

    return H_total;
}

/* Create full density matrix from individual qubit states */
Complex* create_full_density_matrix(LindbladConfig* config) {
    int n_qubits = config->num_qubits;
    if (n_qubits <= 0) return NULL;

    int dim = 1 << n_qubits;
    int total_size = dim * dim;

    // NOTE: assumes product state rho = rho_0 ⊗ rho_1 ⊗ ... ⊗ rho_{n-1}
    Complex* current = malloc(4 * sizeof(Complex));
    if (!current) return NULL;
    memcpy(current, config->qubits[0].rho, 4 * sizeof(Complex));
    int current_size = 2;

    for (int q = 1; q < n_qubits; q++) {
        int new_size  = current_size * 2;
        int new_total = new_size * new_size;

        Complex* new_rho = calloc(new_total, sizeof(Complex));
        if (!new_rho) { free(current); return NULL; }

        tensor_product(current, current_size, config->qubits[q].rho, 2, new_rho);

        free(current);
        current = new_rho;
        current_size = new_size;
    }

    // current is now the full dim x dim density matrix
    return current;
}

/* Print matrix */
void print_matrix(Complex* mat, int size, const char* name) {
    printf("\n=== %s (%dx%d) ===\n", name, size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            Complex c = mat[i * size + j];
            printf("%7.3f%+7.3fi ", c.real, c.imag);
        }
        printf("\n");
    }
    printf("==================\n");
}

/* Free matrix */
void free_matrix(Complex* mat) {
    if (mat) free(mat);
}
