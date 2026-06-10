#include "matrix_operations.h"
#include "simple_tensor.h"

/* Test function for matrix operations */
void test_matrix_operations() {
    printf("\n========== TESTING MATRIX OPERATIONS ==========\n");
    
    /* Create 2x2 Pauli matrices for testing */
    Complex X[4] = {{0,0}, {1,0}, {1,0}, {0,0}};
    Complex Y[4] = {{0,0}, {0,-1}, {0,1}, {0,0}};
    Complex Z[4] = {{1,0}, {0,0}, {0,0}, {-1,0}};
    Complex I[4] = {{1,0}, {0,0}, {0,0}, {1,0}};
    
    /* Test 1: Matrix multiplication - X * X = I */
    printf("\n1. Matrix Multiplication: X * X\n");
    Complex X_squared[4];
    matrix_multiply(X, X, X_squared, 2);
    print_matrix(X_squared, 2, "X^2 (should be I)");
    
    /* Test 2: Commutator [X, Y] = 2iZ */
    printf("\n2. Commutator: [X, Y]\n");
    Complex XY_commutator[4];
    commutator(X, Y, XY_commutator, 2);
    print_matrix(XY_commutator, 2, "[X, Y] (should be 2iZ)");
    
    /* Test 3: Anti-commutator {X, Y} = 0 */
    printf("\n3. Anti-commutator: {X, Y}\n");
    Complex XY_anticommutator[4];
    anticommutator(X, Y, XY_anticommutator, 2);
    print_matrix(XY_anticommutator, 2, "{X, Y} (should be zero)");
    
    /* Test 4: Commutator [X, Z] = -2iY */
    printf("\n4. Commutator: [X, Z]\n");
    Complex XZ_commutator[4];
    commutator(X, Z, XZ_commutator, 2);
    print_matrix(XZ_commutator, 2, "[X, Z] (should be -2iY)");
    
    /* Test 5: Check if Hamiltonian is Hermitian */
    printf("\n5. Hermiticity Check\n");
    float omega = 1.5, E = 0.8;
    Complex* H = single_qubit_hamiltonian(omega, E);
    int hermitian = is_hermitian(H, 2, 1e-6);
    printf("Hamiltonian H = %.1f*Z + %.1f*X is %s\n", 
           omega, E, hermitian ? "Hermitian ✓" : "NOT Hermitian ✗");
    print_matrix(H, 2, "H");
    free(H);
    
    /* Test 6: For Lindblad equation: -i[H, ρ] */
    printf("\n6. Lindblad Commutator: -i[H, ρ]\n");
    Complex rho[4] = {{0.5,0}, {0.5,0}, {0.5,0}, {0.5,0}};
    Complex H_rho_comm[4];
    
    /* Compute [H, ρ] */
    commutator(H, rho, H_rho_comm, 2);
    
    /* Multiply by -i: -i * [H, ρ] */
    Complex lindblad_term[4];
    matrix_scale(H_rho_comm, 0, -1, lindblad_term, 2); /* -i * (a+ib) = -i*a + b */
    printf("Lindblad term -i[H, ρ]:\n");
    print_matrix(lindblad_term, 2, "-i[H, ρ]");
    
    free(H);
}

/* Test with full system (2 qubits) */
void test_full_system() {
    printf("\n========== FULL SYSTEM TEST (2 QUBITS) ==========\n");
    
    /* Create configuration */
    LindbladConfig config;
    config.num_qubits = 2;
    config.qubits[0].omega = 1.5; config.qubits[0].E = 0.8;
    config.qubits[1].omega = 2.3; config.qubits[1].E = 1.2;
    
    /* Create Hamiltonian */
    float omegas[] = {1.5, 2.3};
    float Es[] = {0.8, 1.2};
    Complex* H = create_hamiltonian(omegas, Es, 2);
    
    /* Create initial density matrix */
    Complex* rho = create_full_density_matrix(&config);
    
    /* Compute Lindblad term: -i[H, ρ] */
    int dim = 4; /* 2^2 = 4 */
    Complex* H_rho_comm = (Complex*)calloc(dim * dim, sizeof(Complex));
    Complex* lindblad_term = (Complex*)calloc(dim * dim, sizeof(Complex));
    
    printf("\nHamiltonian (4x4):\n");
    print_matrix(H, dim, "H");
    
    printf("\nInitial Density Matrix (4x4):\n");
    print_matrix(rho, dim, "ρ");
    
    /* Compute commutator [H, ρ] */
    commutator(H, rho, H_rho_comm, dim);
    
    /* Compute -i[H, ρ] */
    matrix_scale(H_rho_comm, 0, -1, lindblad_term, dim);
    
    printf("\nLindblad Evolution Term -i[H, ρ]:\n");
    print_matrix(lindblad_term, dim, "-i[H, ρ]");
    
    /* Verify trace is zero (for valid Lindblad term) */
    float trace = 0;
    for (int i = 0; i < dim; i++) {
        trace += lindblad_term[i * dim + i].real;
    }
    printf("\nTrace of -i[H, ρ] = %.6f (should be ~0)\n", trace);
    
    free(H);
    free(rho);
    free(H_rho_comm);
    free(lindblad_term);
}

/* Test Lindblad dissipator for single qubit */
void test_lindblad_dissipator() {
    printf("\n========== LINDBLAD DISSIPATOR TEST ==========\n");
    
    /* Single qubit density matrix */
    Complex rho[4] = {
        {0.8, 0}, {0.3, 0.1},
        {0.3, -0.1}, {0.2, 0}
    };
    
    /* Lindblad operator L (amplitude damping) */
    Complex L[4] = {
        {0, 0}, {1, 0},
        {0, 0}, {0, 0}
    };
    
    float gamma = 0.1;
    
    printf("Initial ρ:\n");
    print_matrix(rho, 2, "ρ");
    
    printf("\nLindblad operator L (damping):\n");
    print_matrix(L, 2, "L");
    
    /* Compute dissipator: D[ρ] = γ (LρL† - ½{L†L, ρ}) */
    
    /* Step 1: Compute L† */
    Complex L_dagger[4];
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            L_dagger[i * 2 + j].real = L[j * 2 + i].real;
            L_dagger[i * 2 + j].imag = -L[j * 2 + i].imag;
        }
    }
    
    /* Step 2: Compute LρL† */
    Complex Lrho[4], LrhoLdagger[4];
    matrix_multiply(L, rho, Lrho, 2);
    matrix_multiply(Lrho, L_dagger, LrhoLdagger, 2);
    
    /* Step 3: Compute L†L */
    Complex LdaggerL[4];
    matrix_multiply(L_dagger, L, LdaggerL, 2);
    
    /* Step 4: Compute {L†L, ρ} */
    Complex anticommutator_result[4];
    anticommutator(LdaggerL, rho, anticommutator_result, 2);
    
    /* Step 5: Compute ½{L†L, ρ} */
    Complex half_anticommutator[4];
    matrix_scale(anticommutator_result, 0.5, 0, half_anticommutator, 2);
    
    /* Step 6: Compute dissipator = γ (LρL† - ½{L†L, ρ}) */
    Complex dissipator[4];
    matrix_subtract(LrhoLdagger, half_anticommutator, dissipator, 2);
    matrix_scale(dissipator, gamma, 0, dissipator, 2);
    
    printf("\nDissipator D[ρ] = γ (LρL† - ½{L†L, ρ}):\n");
    print_matrix(dissipator, 2, "D[ρ]");
    
    /* Total evolution: dρ/dt = -i[H, ρ] + D[ρ] */
    float omega = 1.0, E = 0.5;
    Complex* H = single_qubit_hamiltonian(omega, E);
    
    Complex commutator_term[4];
    commutator(H, rho, commutator_term, 2);
    Complex lindblad_term[4];
    matrix_scale(commutator_term, 0, -1, lindblad_term, 2);
    
    Complex drho_dt[4];
    matrix_add(lindblad_term, dissipator, drho_dt, 2);
    
    printf("\nTotal evolution dρ/dt = -i[H,ρ] + D[ρ]:\n");
    print_matrix(drho_dt, 2, "dρ/dt");
    
    free(H);
}

int main() {
    /* Run all tests */
    test_matrix_operations();
    test_full_system();
    test_lindblad_dissipator();
    
    printf("\n========== ALL TESTS COMPLETE ==========\n");
    return 0;
}
