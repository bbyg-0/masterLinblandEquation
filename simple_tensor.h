#ifndef SIMPLE_TENSOR_H
#define SIMPLE_TENSOR_H

#include "lindblad_config.h"

/* Tensor product: C = A ⊗ B */
void tensor_product(Complex* A, int N_A, Complex* B, int N_B, Complex* C);

/* Single qubit Hamiltonian (used internally but needed for matrix_operations) */
Complex* single_qubit_hamiltonian(float omega, float E);

/* Create Hamiltonian for N qubits */
Complex* create_hamiltonian(float* omegas, float* Es, int n_qubits);

/* Create full density matrix from individual qubit states */
Complex* create_full_density_matrix(LindbladConfig* config);

/* Print matrix */
void print_matrix(Complex* mat, int size, const char* name);

/* Free matrix */
void free_matrix(Complex* mat);

#endif
