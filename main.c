#include "simple_tensor.h"
#include "matrix_operations.h"

int main(int argc, char* argv[]) {
    LindbladConfig config;
    const char* config_file = (argc > 1) ? argv[1] : "config.json";
    
    printf("\n========================================\n");
    printf("LINDBLAD MASTER EQUATION SIMULATION\n");
    printf("========================================\n");
    
    /* Load configuration */
    if (load_config_from_json(config_file, &config) != 0) {
        if (load_config_from_file(config_file, &config) != 0) {
            printf("WARNING: Could not load configuration from %s\n", config_file);
            printf("Using default configuration...\n\n");
            
            /* Create default configuration */
            config.num_qubits = 2;
            init_all_qubits_default(&config);
            config.qubits[0].omega = 1.5;
            config.qubits[0].E = 0.8;
            config.qubits[1].omega = 2.3;
            config.qubits[1].E = 1.2;
            config.num_disturbances = 0;
        }
    }
    
    /* Print configuration */
    print_config(&config);
    
    /* Create Hamiltonian */
    float* omegas = (float*)malloc(config.num_qubits * sizeof(float));
    float* Es = (float*)malloc(config.num_qubits * sizeof(float));
    
    for (int i = 0; i < config.num_qubits; i++) {
        omegas[i] = config.qubits[i].omega;
        Es[i] = config.qubits[i].E;
    }
    
    Complex* H = create_hamiltonian(omegas, Es, config.num_qubits);
    int dim = 1 << config.num_qubits;
    
    printf("\n=== SYSTEM HAMILTONIAN ===\n");
    print_matrix(H, dim, "Hamiltonian H");
    
    /* Check if Hamiltonian is Hermitian */
    if (is_hermitian(H, dim, 1e-6)) {
        printf("✓ Hamiltonian is Hermitian (valid)\n");
    } else {
        printf("✗ Hamiltonian is NOT Hermitian (invalid)\n");
    }
    
    /* Create full density matrix from initial states */
    Complex* rho = create_full_density_matrix(&config);
    printf("\n=== INITIAL DENSITY MATRIX ===\n");
    print_matrix(rho, dim, "ρ(0)");
    
    /* Verify trace of density matrix = 1 */
    float trace = 0;
    for (int i = 0; i < dim; i++) {
        trace += rho[i * dim + i].real;
    }
    printf("Trace of ρ(0) = %.6f (should be 1.0)\n", trace);
    
    /* Compute Lindblad term: -i[H, ρ] */
    Complex* H_rho_comm = (Complex*)calloc(dim * dim, sizeof(Complex));
    Complex* lindblad_term = (Complex*)calloc(dim * dim, sizeof(Complex));
    
    commutator(H, rho, H_rho_comm, dim);
    matrix_scale(H_rho_comm, 0, -1, lindblad_term, dim);
    
    printf("\n=== LINDBLAD EVOLUTION TERM ===\n");
    printf("dρ/dt = -i[H, ρ] + Σ γ_k (L_k ρ L_k† - ½{L_k†L_k, ρ})\n\n");
    print_matrix(lindblad_term, dim, "-i[H, ρ]");
    
    /* Process disturbances (Lindblad operators) */
    Complex* dissipator = NULL;
    
    if (config.num_disturbances > 0) {
        printf("\n=== LINDBLAD DISSIPATORS ===\n");
        
        /* Create full dissipator term */
        dissipator = (Complex*)calloc(dim * dim, sizeof(Complex));
        
        for (int d = 0; d < config.num_disturbances; d++) {
            Disturbance* dist = &config.disturbances[d];
            printf("\nDisturbance %d: Qubit %d, Type: %s, γ = %.4f\n", 
                   d, dist->qubit_index, disturbance_type_to_string(dist->type), dist->gamma);
            
            /* Create Lindblad operator L for this disturbance */
            Complex* L = create_lindblad_operator(dist->type, dist->qubit_index, config.num_qubits);
            
            if (L) {
                /* Compute L ρ L† */
                Complex* Lrho = (Complex*)calloc(dim * dim, sizeof(Complex));
                Complex* LrhoLdagger = (Complex*)calloc(dim * dim, sizeof(Complex));
                Complex* Ldagger = (Complex*)calloc(dim * dim, sizeof(Complex));
                
                /* Create L† */
                for (int i = 0; i < dim * dim; i++) {
                    Ldagger[i].real = L[i].real;
                    Ldagger[i].imag = -L[i].imag;
                }
                
                /* Compute L ρ */
                matrix_multiply(L, rho, Lrho, dim);
                
                /* Compute (L ρ) L† */
                matrix_multiply(Lrho, Ldagger, LrhoLdagger, dim);
                
                /* Compute L†L */
                Complex* LdaggerL = (Complex*)calloc(dim * dim, sizeof(Complex));
                matrix_multiply(Ldagger, L, LdaggerL, dim);
                
                /* Compute {L†L, ρ} */
                Complex* anticommutator_result = (Complex*)calloc(dim * dim, sizeof(Complex));
                anticommutator(LdaggerL, rho, anticommutator_result, dim);
                
                /* Compute ½{L†L, ρ} */
                Complex* half_anticommutator = (Complex*)calloc(dim * dim, sizeof(Complex));
                matrix_scale(anticommutator_result, 0.5, 0, half_anticommutator, dim);
                
                /* Compute D = γ (LρL† - ½{L†L, ρ}) */
                Complex* D = (Complex*)calloc(dim * dim, sizeof(Complex));
                matrix_subtract(LrhoLdagger, half_anticommutator, D, dim);
                matrix_scale(D, dist->gamma, 0, D, dim);
                
                /* Add to total dissipator */
                for (int i = 0; i < dim * dim; i++) {
                    dissipator[i].real += D[i].real;
                    dissipator[i].imag += D[i].imag;
                }
                
                printf("  Added dissipator term with γ = %.4f\n", dist->gamma);
                
                /* Clean up */
                free(L);
                free(Lrho);
                free(LrhoLdagger);
                free(Ldagger);
                free(LdaggerL);
                free(anticommutator_result);
                free(half_anticommutator);
                free(D);
            }
        }
        
        printf("\n=== TOTAL DISSIPATOR ===\n");
        print_matrix(dissipator, dim, "Σ γ_k (L_k ρ L_k† - ½{L_k†L_k, ρ})");
    } else {
        printf("\nNo disturbances specified.\n");
        printf("Evolution is unitary: dρ/dt = -i[H, ρ]\n");
    }
    
    /* Total evolution: dρ/dt = -i[H, ρ] + dissipator */
    Complex* drho_dt = (Complex*)calloc(dim * dim, sizeof(Complex));
    
    if (dissipator) {
        matrix_add(lindblad_term, dissipator, drho_dt, dim);
        printf("\n=== TOTAL EVOLUTION ===\n");
        print_matrix(drho_dt, dim, "dρ/dt");
    } else {
        matrix_copy(lindblad_term, drho_dt, dim);
    }
    
    /* Simple Euler step example */
    float dt = 0.01;
    Complex* rho_next = (Complex*)calloc(dim * dim, sizeof(Complex));
    
    /* rho_next = rho + dt * drho_dt */
    for (int i = 0; i < dim * dim; i++) {
        rho_next[i].real = rho[i].real + dt * drho_dt[i].real;
        rho_next[i].imag = rho[i].imag + dt * drho_dt[i].imag;
    }
    
    printf("\n=== TIME EVOLUTION (Euler step, dt = %.3f) ===\n", dt);
    printf("ρ(t+dt) = ρ(t) + dt * dρ/dt\n");
    print_matrix(rho_next, dim, "ρ(t+dt)");
    
    /* Verify trace preservation */
    float new_trace = 0;
    for (int i = 0; i < dim; i++) {
        new_trace += rho_next[i * dim + i].real;
    }
    printf("Trace after evolution = %.6f (should still be 1.0)\n", new_trace);
    
    /* Cleanup */
    free(omegas);
    free(Es);
    free_matrix(H);
    free_matrix(rho);
    free(H_rho_comm);
    free(lindblad_term);
    if (dissipator) free(dissipator);
    free(drho_dt);
    free(rho_next);
    
    printf("\n========================================\n");
    printf("SIMULATION SETUP COMPLETE\n");
    printf("========================================\n");
    
    return 0;
}
