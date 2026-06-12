#include "simple_tensor.h"
#include "matrix_operations.h"
#include <math.h>

// Compute time derivative for given density matrix
void compute_time_derivative(Complex* H, Complex* rho, LindbladConfig* config, 
                              Complex* drho_dt, int dim);


// Full Runge-Kutta 4th order for density matrix evolution
void runge_kutta_4(Complex* H, Complex* rho, LindbladConfig* config, 
                   Complex* rho_next, int dim, float dt);

// Save density matrix to file for visualization
void save_density_matrix_to_file(Complex* rho, int dim, float time, const char* filename);

// Save expectation values for each qubit (for Bloch sphere visualization)
void save_expectation_values(Complex* rho, LindbladConfig* config, float time, const char* filename);

// Helper: Create single qubit operator for N-qubit system
Complex* create_single_qubit_operator(Complex* op, int qubit_idx, int n_qubits);


int main(int argc, char* argv[]) {
    LindbladConfig config;
    const char* config_file = (argc > 1) ? argv[1] : "config.json";
    const char* output_file = (argc > 2) ? argv[2] : "density_matrix.dat";
    const char* expectation_file = (argc > 3) ? argv[3] : "expectations.dat";
    
    printf("\n========================================\n");
    printf("LINDBLAD MASTER EQUATION SIMULATION\n");
    printf("Runge-Kutta 4th Order Integration\n");
    printf("========================================\n");
    
    // Load configuration
    if (load_config_from_json(config_file, &config) != 0) {
        if (load_config_from_file(config_file, &config) != 0) {
            printf("WARNING: Could not load configuration from %s\n", config_file);
            printf("Using default configuration...\n\n");
            
            // Create default configuration
            config.num_qubits = 2;
            init_all_qubits_default(&config);
            config.qubits[0].omega = 1.5;
            config.qubits[0].E = 0.8;
            config.qubits[1].omega = 2.3;
            config.qubits[1].E = 1.2;
            config.num_disturbances = 3;
            
            // Add some default disturbances
            parse_disturbance_string("0d0.1", &config.disturbances[0]);
            parse_disturbance_string("1a0.05", &config.disturbances[1]);
            parse_disturbance_string("0s0.02", &config.disturbances[2]);
        }
    }
    
    // Print configuration
    print_config(&config);
    
    // Create Hamiltonian
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
    
    // Check if Hamiltonian is Hermitian
    if (is_hermitian(H, dim, 1e-6)) {
        printf("✓ Hamiltonian is Hermitian (valid)\n");
    } else {
        printf("✗ Hamiltonian is NOT Hermitian (invalid)\n");
    }
    
    // Create full density matrix from initial states
    Complex* rho = create_full_density_matrix(&config);
    printf("\n=== INITIAL DENSITY MATRIX ===\n");
    print_matrix(rho, dim, "ρ(0)");
    
    // Verify trace of density matrix = 1
    float trace = 0;
    for (int i = 0; i < dim; i++) {
        trace += rho[i * dim + i].real;
    }
    printf("Trace of ρ(0) = %.6f (should be 1.0)\n", trace);
    
    // Simulation parameters
    float t_total = 5.0f;      // Total simulation time
    float dt = 0.001f;           // Time step
    int n_steps = (int)(t_total / dt);
    
    printf("\n=== SIMULATION PARAMETERS ===\n");
    printf("Total time: %.2f\n", t_total);
    printf("Time step: %.4f (Runge-Kutta 4th order)\n", dt);
    printf("Number of steps: %d\n", n_steps);
    printf("Output file: %s\n", output_file);
    printf("Expectation file: %s\n", expectation_file);
    
    // Save initial state
    save_density_matrix_to_file(rho, dim, 0.0f, output_file);
    save_expectation_values(rho, &config, 0.0f, expectation_file);
    
    // Time evolution using Runge-Kutta 4
    printf("\n=== EVOLVING SYSTEM ===\n");
    
    Complex* rho_current = (Complex*)malloc(dim * dim * sizeof(Complex));
    Complex* rho_next = (Complex*)malloc(dim * dim * sizeof(Complex));
    memcpy(rho_current, rho, dim * dim * sizeof(Complex));
    
    for (int step = 1; step <= n_steps; step++) {
        float t = step * dt;
        
        // Perform Runge-Kutta 4 step
        runge_kutta_4(H, rho_current, &config, rho_next, dim, dt);
        
        // Verify trace preservation (optional, for debugging)
        float new_trace = 0;
        for (int i = 0; i < dim; i++) {
            new_trace += rho_next[i * dim + i].real;
        }
        
        // Save state every 10 steps (to keep file size reasonable)
        if (step % 10 == 0 || step == n_steps) {
            save_density_matrix_to_file(rho_next, dim, t, output_file);
            save_expectation_values(rho_next, &config, t, expectation_file);
            
            // Progress indicator
            if (step % 100 == 0 || step == n_steps) {
                printf("  Progress: %.1f%% (t = %.2f, trace = %.6f)\n", 
                       100.0f * step / n_steps, t, new_trace);
            }
        }
        
        // Swap for next iteration
        Complex* temp = rho_current;
        rho_current = rho_next;
        rho_next = temp;
    }
    
    // Final state
    printf("\n=== FINAL DENSITY MATRIX ===\n");
    print_matrix(rho_current, dim, "ρ(final)");
    
    float final_trace = 0;
    for (int i = 0; i < dim; i++) {
        final_trace += rho_current[i * dim + i].real;
    }
    printf("Final trace = %.10f (should be 1.0)\n", final_trace);
    
    // Cleanup
    free(omegas);
    free(Es);
    free_matrix(H);
    free_matrix(rho);
    free(rho_current);
    free(rho_next);
    printf("\n========================================\n");
    printf("SIMULATION COMPLETE\n");
    printf("Results saved to:\n");
    printf("  - %s (density matrix evolution)\n", output_file);
    printf("  - %s (expectation values for Bloch sphere)\n", expectation_file);
    printf("========================================\n");
    
    return 0;
}


// Helper: build full N-qubit operator from a 2x2 single-qubit matrix acting on qubit `qubit_idx`
static Complex* build_full_operator_from_single(const Complex single[4], int qubit_idx, int n_qubits) {
    int dim = 1 << n_qubits;
    Complex* full = (Complex*)calloc(dim * dim, sizeof(Complex));
    if (!full) return NULL;

    for (int row = 0; row < dim; row++) {
        for (int col = 0; col < dim; col++) {
            int row_bit = (row >> qubit_idx) & 1;
            int col_bit = (col >> qubit_idx) & 1;
            int match = 1;
            for (int q = 0; q < n_qubits; q++) {
                if (q == qubit_idx) continue;
                if (((row >> q) & 1) != ((col >> q) & 1)) {
                    match = 0;
                    break;
                }
            }
            if (match) {
                full[row * dim + col] = single[row_bit * 2 + col_bit];
            }
        }
    }
    return full;
}

// Compute time derivative for given density matrix
void compute_time_derivative(Complex* H, Complex* rho, LindbladConfig* config, 
                              Complex* drho_dt, int dim) {
    // Compute -i[H, ρ]
    Complex* H_rho_comm = (Complex*)calloc(dim * dim, sizeof(Complex));
    commutator(H, rho, H_rho_comm, dim);
    matrix_scale(H_rho_comm, 0, -1, drho_dt, dim);
    free(H_rho_comm);

    // Add dissipator contributions
    if (config->num_disturbances == 0) return;

    int n_qubits = config->num_qubits;

    for (int d = 0; d < config->num_disturbances; d++) {
        Disturbance* dist = &config->disturbances[d];
        float gamma = dist->gamma;

        // For depolarizing, we need three Lindblad operators: √(γ/3) X, √(γ/3) Y, √(γ/3) Z
        if (dist->type == DEPOLARIZING) {
            float gamma3 = gamma / 3.0f;
            // Pauli X
            Complex X[4] = {{0,0}, {1,0}, {1,0}, {0,0}};
            // Pauli Y
            Complex Y[4] = {{0,0}, {0,-1}, {0,1}, {0,0}};
            // Pauli Z
            Complex Z[4] = {{1,0}, {0,0}, {0,0}, {-1,0}};

            Complex* ops[3] = {X, Y, Z};
            for (int p = 0; p < 3; p++) {
                Complex* L = build_full_operator_from_single(ops[p], dist->qubit_index, n_qubits);
                if (!L) continue;

                // L ρ L†
                Complex* Lrho = calloc(dim * dim, sizeof(Complex));
                Complex* LrhoLdagger = calloc(dim * dim, sizeof(Complex));
                Complex* Ldagger = calloc(dim * dim, sizeof(Complex));
		// Correct: conjugate transpose
for (int row = 0; row < dim; row++)
    for (int col = 0; col < dim; col++) {
        Ldagger[col * dim + row].real =  L[row * dim + col].real;
        Ldagger[col * dim + row].imag = -L[row * dim + col].imag;
    }
                matrix_multiply(L, rho, Lrho, dim);
                matrix_multiply(Lrho, Ldagger, LrhoLdagger, dim);

                // L†L
                Complex* LdaggerL = calloc(dim * dim, sizeof(Complex));
                matrix_multiply(Ldagger, L, LdaggerL, dim);

                // {L†L, ρ}
                Complex* anticom = calloc(dim * dim, sizeof(Complex));
                anticommutator(LdaggerL, rho, anticom, dim);

                // D = γ3 * (LρL† - ½{L†L, ρ})
                Complex* D = calloc(dim * dim, sizeof(Complex));
                matrix_scale(anticom, 0.5f, 0, D, dim);
                matrix_subtract(LrhoLdagger, D, D, dim);
                matrix_scale(D, gamma3, 0, D, dim);

                // Add to drho_dt
                for (int i = 0; i < dim * dim; i++) {
                    drho_dt[i].real += D[i].real;
                    drho_dt[i].imag += D[i].imag;
                }

                free(L); free(Lrho); free(LrhoLdagger); free(Ldagger);
                free(LdaggerL); free(anticom); free(D);
            }
            continue; // depolarizing done
        }

        // For DAMPING or DEPHASING: build single‑qubit operator
        Complex single_op[4] = {0};
        switch (dist->type) {
            case DAMPING:
                single_op[1].real = 1.0f;   // σ⁻ = |0><1|
                break;
            case DEPHASING:
                single_op[0].real = 1.0f;
                single_op[3].real = -1.0f;  // σ_z
                break;
            default:
                continue;
        }

        Complex* L = build_full_operator_from_single(single_op, dist->qubit_index, n_qubits);
        if (!L) continue;

        // L ρ L†
        Complex* Lrho = calloc(dim * dim, sizeof(Complex));
        Complex* LrhoLdagger = calloc(dim * dim, sizeof(Complex));
        Complex* Ldagger = calloc(dim * dim, sizeof(Complex));
        for (int i = 0; i < dim * dim; i++) {
            Ldagger[i].real = L[i].real;
            Ldagger[i].imag = -L[i].imag;
        }
        matrix_multiply(L, rho, Lrho, dim);
        matrix_multiply(Lrho, Ldagger, LrhoLdagger, dim);

        // L†L
        Complex* LdaggerL = calloc(dim * dim, sizeof(Complex));
        matrix_multiply(Ldagger, L, LdaggerL, dim);

        // {L†L, ρ}
        Complex* anticom = calloc(dim * dim, sizeof(Complex));
        anticommutator(LdaggerL, rho, anticom, dim);

        // D = γ (LρL† - ½{L†L, ρ})
        Complex* D = calloc(dim * dim, sizeof(Complex));
        matrix_scale(anticom, 0.5f, 0, D, dim);
        matrix_subtract(LrhoLdagger, D, D, dim);
        matrix_scale(D, gamma, 0, D, dim);

        // Add to drho_dt
        for (int i = 0; i < dim * dim; i++) {
            drho_dt[i].real += D[i].real;
            drho_dt[i].imag += D[i].imag;
        }

        free(L); free(Lrho); free(LrhoLdagger); free(Ldagger);
        free(LdaggerL); free(anticom); free(D);
    }
}



void runge_kutta_4(Complex* H, Complex* rho, LindbladConfig* config, 
                   Complex* rho_next, int dim, float dt) {
    // Allocate temporary matrices for RK4 stages
    Complex* k1 = (Complex*)calloc(dim * dim, sizeof(Complex));
    Complex* k2 = (Complex*)calloc(dim * dim, sizeof(Complex));
    Complex* k3 = (Complex*)calloc(dim * dim, sizeof(Complex));
    Complex* k4 = (Complex*)calloc(dim * dim, sizeof(Complex));
    Complex* temp = (Complex*)calloc(dim * dim, sizeof(Complex));
    
    // k1 = f(rho, t)
    compute_time_derivative(H, rho, config, k1, dim);
    
    // k2 = f(rho + dt/2 * k1, t + dt/2)
    for (int i = 0; i < dim * dim; i++) {
        temp[i].real = rho[i].real + 0.5f * dt * k1[i].real;
        temp[i].imag = rho[i].imag + 0.5f * dt * k1[i].imag;
    }
    compute_time_derivative(H, temp, config, k2, dim);
    
    // k3 = f(rho + dt/2 * k2, t + dt/2)
    for (int i = 0; i < dim * dim; i++) {
        temp[i].real = rho[i].real + 0.5f * dt * k2[i].real;
        temp[i].imag = rho[i].imag + 0.5f * dt * k2[i].imag;
    }
    compute_time_derivative(H, temp, config, k3, dim);
    
    // k4 = f(rho + dt * k3, t + dt)
    for (int i = 0; i < dim * dim; i++) {
        temp[i].real = rho[i].real + dt * k3[i].real;
        temp[i].imag = rho[i].imag + dt * k3[i].imag;
    }
    compute_time_derivative(H, temp, config, k4, dim);
    
    // rho_next = rho + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    for (int i = 0; i < dim * dim; i++) {
        rho_next[i].real = rho[i].real + dt/6.0f * (k1[i].real + 2.0f*k2[i].real + 2.0f*k3[i].real + k4[i].real);
        rho_next[i].imag = rho[i].imag + dt/6.0f * (k1[i].imag + 2.0f*k2[i].imag + 2.0f*k3[i].imag + k4[i].imag);
    }
    
    // Clean up
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(temp);
}

void save_density_matrix_to_file(Complex* rho, int dim, float time, const char* filename) {
    static int first_write = 1;
    FILE* file;
    
    if (first_write) {
        file = fopen(filename, "w");
        if (!file) {
            printf("Error: Cannot create output file %s\n", filename);
            return;
        }
        // Write header
        fprintf(file, "# Lindblad Master Equation Simulation Results\n");
        fprintf(file, "# Format: time dim real_00 imag_00 real_01 imag_01 ... real_NN imag_NN\n");
        fprintf(file, "# Density matrix is stored row-major order\n\n");
        first_write = 0;
    } else {
        file = fopen(filename, "a");
        if (!file) {
            printf("Error: Cannot append to output file %s\n", filename);
            return;
        }
    }
    
    // Write time and matrix dimension
    fprintf(file, "%.8f %d", time, dim);
    
    // Write all density matrix elements (real and imaginary parts)
    for (int i = 0; i < dim * dim; i++) {
        fprintf(file, " %.8f %.8f", rho[i].real, rho[i].imag);
    }
    fprintf(file, "\n");
    
    fclose(file);
}

void save_expectation_values(Complex* rho, LindbladConfig* config, float time, const char* filename) {
    static int first_write_exp = 1;
    int n_qubits = config->num_qubits;
    int dim = 1 << n_qubits;

    // Pauli matrices — row-major: [0,0], [0,1], [1,0], [1,1]
    const Complex X[4] = {{0,0}, {1,0}, {1,0}, {0,0}};
    const Complex Y[4] = {{0,0}, {0,1}, {0,-1}, {0,0}};   // Y[0,1]=+i, Y[1,0]=-i
    const Complex Z[4] = {{1,0}, {0,0}, {0,0}, {-1,0}};

    FILE* file;
    if (first_write_exp) {
        file = fopen(filename, "w");
        if (!file) return;
        fprintf(file, "# Time");
        for (int q = 0; q < n_qubits; q++)
            fprintf(file, " X%d Y%d Z%d", q, q, q);
        fprintf(file, "\n");
        first_write_exp = 0;
    } else {
        file = fopen(filename, "a");
        if (!file) return;
    }

    fprintf(file, "%.8f", time);

    for (int q = 0; q < n_qubits; q++) {
        float exp_x = 0, exp_y = 0, exp_z = 0;

        // <σ> = Tr(σ ρ_reduced) = Σ_{a,b} σ_{ab} * ρ_reduced_{ba}
        // ρ_reduced_{ba} = Σ_{m: m without q-bit} ρ[ (b<<q)|m, (a<<q)|m ]
        for (int a = 0; a < 2; a++) {
            for (int b = 0; b < 2; b++) {
                // σ_{ab}: row=a, col=b
                Complex sx = X[a * 2 + b];
                Complex sy = Y[a * 2 + b];
                Complex sz = Z[a * 2 + b];

                if (sx.real == 0 && sx.imag == 0 &&
                    sy.real == 0 && sy.imag == 0 &&
                    sz.real == 0 && sz.imag == 0) continue;

                // Sum over all basis states of the other qubits
                for (int m = 0; m < (dim >> 1); m++) {
                    // Reconstruct full indices by inserting qubit q's bit
                    // m is the index over the (n-1) remaining qubits
                    int low  = m & ((1 << q) - 1);
                    int high = (m >> q) << (q + 1);

                    int row = high | (b << q) | low;   // full row index, q-bit = b
                    int col = high | (a << q) | low;   // full col index, q-bit = a

                    // ρ_reduced_{ba} contribution: ρ[row, col]
                    Complex rho_elem = rho[row * dim + col];

                    // Tr contribution: Re(σ_{ab} * ρ_{ba})
                    exp_x += sx.real * rho_elem.real - sx.imag * rho_elem.imag;
                    exp_y += sy.real * rho_elem.real - sy.imag * rho_elem.imag;
                    exp_z += sz.real * rho_elem.real - sz.imag * rho_elem.imag;
                }
            }
        }

        fprintf(file, " %.8f %.8f %.8f", exp_x, exp_y, exp_z);
    }

    fprintf(file, "\n");
    fclose(file);
}



Complex* create_single_qubit_operator(Complex* op, int qubit_idx, int n_qubits) {
    int dim = 1 << n_qubits;
    Complex* full_op = (Complex*)calloc(dim * dim, sizeof(Complex));
    
    Complex I[4] = {{1,0}, {0,0}, {0,0}, {1,0}};
    Complex* current = NULL;
    
    for (int i = 0; i < n_qubits; i++) {
        if (i == 0) {
            current = (Complex*)malloc(4 * sizeof(Complex));
            if (i == qubit_idx)
                memcpy(current, op, 4 * sizeof(Complex));
            else
                memcpy(current, I, 4 * sizeof(Complex));
        } else {
            int current_size = 1 << i;
            int new_size = current_size * 2;
            int new_total = new_size * new_size;
            
            Complex* new_current = (Complex*)calloc(new_total, sizeof(Complex));
            Complex* next_op = (i == qubit_idx) ? op : I;
            
            tensor_product(current, current_size, next_op, 2, new_current);
            
            free(current);
            current = new_current;
        }
    }
    
    if (current) {
        memcpy(full_op, current, dim * dim * sizeof(Complex));
        free(current);
    }
    
    return full_op;
}

