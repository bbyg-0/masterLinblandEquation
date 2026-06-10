#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include "lindblad_config.h"
#include "simple_tensor.h"

/* Matrix multiplication: C = A * B */
void matrix_multiply(Complex* A, Complex* B, Complex* C, int size);

/* Commutator: [A, B] = A*B - B*A */
void commutator(Complex* A, Complex* B, Complex* result, int size);

/* Anti-commutator: {A, B} = A*B + B*A */
void anticommutator(Complex* A, Complex* B, Complex* result, int size);

/* Helper: Copy matrix */
void matrix_copy(Complex* src, Complex* dst, int size);

/* Helper: Add two matrices: C = A + B */
void matrix_add(Complex* A, Complex* B, Complex* C, int size);

/* Helper: Subtract two matrices: C = A - B */
void matrix_subtract(Complex* A, Complex* B, Complex* C, int size);

/* Helper: Scale matrix: C = scalar * A */
void matrix_scale(Complex* A, float scalar_real, float scalar_imag, Complex* C, int size);

/* Helper: Check if matrix is Hermitian */
int is_hermitian(Complex* A, int size, float tolerance);

/* Create Lindblad operator based on disturbance type */
Complex* create_lindblad_operator(DisturbanceType type, int qubit_idx, int n_qubits);

#endif
