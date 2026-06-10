#ifndef LINDBLAD_CONFIG_H
#define LINDBLAD_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_QUBITS 100
#define MAX_DISTURBANCES 1000
#define MAX_STRING_LEN 100

/* Complex number structure for density matrix */
typedef struct {
    float real;
    float imag;
} Complex;

/* Disturbance type enumeration */
typedef enum {
    DAMPING,        /* 'd' - Amplitude damping */
    DEPHASING,      /* 'a' - Dephasing (phase damping) */
    DEPOLARIZING,   /* 's' - Depolarizing noise */
    UNKNOWN_DISTURBANCE
} DisturbanceType;

/* Structure to store disturbance information */
typedef struct {
    int qubit_index;
    DisturbanceType type;
    float gamma;
    char raw_string[MAX_STRING_LEN];
} Disturbance;

/* Structure to store qubit information */
typedef struct {
    float omega;        /* Qubit frequency */
    float E;           /* Energy/amplitude parameter */
    Complex rho[4];    /* 2x2 density matrix for initial state {{0.5,0.5},{0.5,0.5}} */
} QubitInfo;

/* Main configuration structure */
typedef struct {
    int num_qubits;
    QubitInfo qubits[MAX_QUBITS];
    int num_disturbances;
    Disturbance disturbances[MAX_DISTURBANCES];
} LindbladConfig;

/* Function prototypes */
int parse_disturbance_string(const char* str, Disturbance* disturbance);
int load_config_from_json(const char* filename, LindbladConfig* config);
int load_config_from_file(const char* filename, LindbladConfig* config);
void print_config(const LindbladConfig* config);
void init_default_qubit_state(QubitInfo* qubit);
void init_all_qubits_default(LindbladConfig* config);
const char* disturbance_type_to_string(DisturbanceType type);

#endif
