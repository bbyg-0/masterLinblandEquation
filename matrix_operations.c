#include "matrix_operations.h"
#include <math.h>

/* Matrix multiplication: C = A * B */
void matrix_multiply(Complex* A, Complex* B, Complex* C, int size) {
    /* Clear C first */
    memset(C, 0, size * size * sizeof(Complex));
    
    /* Naive multiplication - optimized for small quantum matrices */
    for (int i = 0; i < size; i++) {
        for (int k = 0; k < size; k++) {
            Complex a = A[i * size + k];
            if (a.real == 0 && a.imag == 0) continue;
            
            for (int j = 0; j < size; j++) {
                Complex b = B[k * size + j];
                if (b.real == 0 && b.imag == 0) continue;
                
                int idx = i * size + j;
                
                /* C[i][j] += A[i][k] * B[k][j] */
                C[idx].real += a.real * b.real - a.imag * b.imag;
                C[idx].imag += a.real * b.imag + a.imag * b.real;
            }
        }
    }
}

/* Commutator: [A, B] = A*B - B*A */
void commutator(Complex* A, Complex* B, Complex* result, int size) {
    /* Temporary matrices for products */
    Complex* AB = (Complex*)malloc(size * size * sizeof(Complex));
    Complex* BA = (Complex*)malloc(size * size * sizeof(Complex));
    
    if (!AB || !BA) {
        printf("Error: Memory allocation failed in commutator\n");
        free(AB);
        free(BA);
        return;
    }
    
    /* Compute A*B and B*A */
    matrix_multiply(A, B, AB, size);
    matrix_multiply(B, A, BA, size);
    
    /* result = AB - BA */
    matrix_subtract(AB, BA, result, size);
    
    free(AB);
    free(BA);
}

/* Anti-commutator: {A, B} = A*B + B*A */
void anticommutator(Complex* A, Complex* B, Complex* result, int size) {
    Complex* AB = (Complex*)malloc(size * size * sizeof(Complex));
    Complex* BA = (Complex*)malloc(size * size * sizeof(Complex));
    
    if (!AB || !BA) {
        printf("Error: Memory allocation failed in anticommutator\n");
        free(AB);
        free(BA);
        return;
    }
    
    /* Compute A*B and B*A */
    matrix_multiply(A, B, AB, size);
    matrix_multiply(B, A, BA, size);
    
    /* result = AB + BA */
    matrix_add(AB, BA, result, size);
    
    free(AB);
    free(BA);
}

/* Copy matrix: dst = src */
void matrix_copy(Complex* src, Complex* dst, int size) {
    memcpy(dst, src, size * size * sizeof(Complex));
}

/* Add two matrices: C = A + B */
void matrix_add(Complex* A, Complex* B, Complex* C, int size) {
    for (int i = 0; i < size * size; i++) {
        C[i].real = A[i].real + B[i].real;
        C[i].imag = A[i].imag + B[i].imag;
    }
}

/* Subtract two matrices: C = A - B */
void matrix_subtract(Complex* A, Complex* B, Complex* C, int size) {
    for (int i = 0; i < size * size; i++) {
        C[i].real = A[i].real - B[i].real;
        C[i].imag = A[i].imag - B[i].imag;
    }
}

/* Scale matrix: C = scalar * A */
void matrix_scale(Complex* A, float scalar_real, float scalar_imag, Complex* C, int size) {
    for (int i = 0; i < size * size; i++) {
        /* (a+ib)*(c+id) = (ac - bd) + i(ad + bc) */
        C[i].real = A[i].real * scalar_real - A[i].imag * scalar_imag;
        C[i].imag = A[i].real * scalar_imag + A[i].imag * scalar_real;
    }
}

/* Check if matrix is Hermitian: A = A† */
int is_hermitian(Complex* A, int size, float tolerance) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            Complex a = A[i * size + j];
            Complex a_conj = A[j * size + i];
            a_conj.imag = -a_conj.imag;
            
            if (fabsf(a.real - a_conj.real) > tolerance ||
                fabsf(a.imag - a_conj.imag) > tolerance) {
                return 0; /* Not Hermitian */
            }
        }
    }
    return 1; /* Is Hermitian */
}

/* Create Lindblad operator based on disturbance type */
Complex* create_lindblad_operator(DisturbanceType type, int qubit_idx, int n_qubits) {
    int dim = 1 << n_qubits;
    Complex* L = (Complex*)calloc(dim * dim, sizeof(Complex));
    
    if (!L) return NULL;
    
    /* Create single-qubit operator first */
    Complex single_qubit_op[4] = {0};
    
    switch(type) {
        case DAMPING:  /* σ- = |0><1| */
            single_qubit_op[1].real = 1.0f;  /* |0><1| */
            break;
            
        case DEPHASING:  /* σz */
            single_qubit_op[0].real = 1.0f;
            single_qubit_op[3].real = -1.0f;
            break;
            
        case DEPOLARIZING:  /* σx for simplicity */
            single_qubit_op[1].real = 1.0f;
            single_qubit_op[2].real = 1.0f;
            break;
            
        default:
            free(L);
            return NULL;
    }
    
    /* Build full operator: I ⊗ I ⊗ ... ⊗ L_q ⊗ ... ⊗ I */
    Complex I[4] = {{1,0}, {0,0}, {0,0}, {1,0}};
    
    Complex* current = NULL;
    
    for (int i = 0; i < n_qubits; i++) {
        if (i == 0) {
            current = (Complex*)malloc(4 * sizeof(Complex));
            if (i == qubit_idx)
                memcpy(current, single_qubit_op, 4 * sizeof(Complex));
            else
                memcpy(current, I, 4 * sizeof(Complex));
        } else {
            int current_size = 1 << i;
            int new_size = current_size * 2;
            int new_total = new_size * new_size;
            
            Complex* new_current = (Complex*)calloc(new_total, sizeof(Complex));
            Complex* next_op = (i == qubit_idx) ? single_qubit_op : I;
            
            tensor_product(current, current_size, next_op, 2, new_current);
            
            free(current);
            current = new_current;
        }
    }
    
    if (current) {
        memcpy(L, current, dim * dim * sizeof(Complex));
        free(current);
    }
    
    return L;
}
