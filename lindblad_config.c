#include "lindblad_config.h"

/* Initialize single qubit to default state {{0.5, 0.5}, {0.5, 0.5}} */
void init_default_qubit_state(QubitInfo* qubit) {
    /* Default density matrix: all elements = 0.5 */
    qubit->rho[0].real = 0.5f; qubit->rho[0].imag = 0.0f;  /* ρ[0][0] */
    qubit->rho[1].real = 0.5f; qubit->rho[1].imag = 0.0f;  /* ρ[0][1] */
    qubit->rho[2].real = 0.5f; qubit->rho[2].imag = 0.0f;  /* ρ[1][0] */
    qubit->rho[3].real = 0.5f; qubit->rho[3].imag = 0.0f;  /* ρ[1][1] */
}

/* Initialize all qubits to default state */
void init_all_qubits_default(LindbladConfig* config) {
    for (int i = 0; i < config->num_qubits; i++) {
        init_default_qubit_state(&config->qubits[i]);
        /* Default omega and E if not set */
        config->qubits[i].omega = 1.0f;
        config->qubits[i].E = 1.0f;
    }
}

/* Convert disturbance type to string */
const char* disturbance_type_to_string(DisturbanceType type) {
    switch(type) {
        case DAMPING: return "damping";
        case DEPHASING: return "dephasing";
        case DEPOLARIZING: return "depolarizing";
        default: return "unknown";
    }
}

/* Parse disturbance string: "0d0.1" means qubit 0, damping, gamma=0.1 */
int parse_disturbance_string(const char* str, Disturbance* disturbance) {
    if (!str || !disturbance) return -1;
    
    memset(disturbance, 0, sizeof(Disturbance));
    strncpy(disturbance->raw_string, str, MAX_STRING_LEN - 1);
    
    /* Parse qubit index */
    int i = 0;
    int qubit_index = 0;
    while (isdigit(str[i])) {
        qubit_index = qubit_index * 10 + (str[i] - '0');
        i++;
    }
    disturbance->qubit_index = qubit_index;
    if (i == 0) return -1;
    
    /* Parse disturbance type */
    if (!isalpha(str[i])) return -1;
    switch(str[i]) {
        case 'd': disturbance->type = DAMPING; break;
        case 'a': disturbance->type = DEPHASING; break;
        case 's': disturbance->type = DEPOLARIZING; break;
        default: return -1;
    }
    i++;
    
    /* Parse gamma value */
    char* endptr;
    disturbance->gamma = (float)strtod(&str[i], &endptr);
    if (endptr == &str[i]) return -1;
    
    return 0;
}

/* Simple JSON parser helper */
static int parse_json_value(const char* json, const char* key, char* value, int value_size) {
    char search_key[MAX_STRING_LEN];
    snprintf(search_key, sizeof(search_key), "\"%s\"", key);
    
    char* key_pos = strstr(json, search_key);
    if (!key_pos) return -1;
    
    char* colon_pos = strchr(key_pos, ':');
    if (!colon_pos) return -1;
    colon_pos++;
    
    while (*colon_pos == ' ' || *colon_pos == '\t') colon_pos++;
    
    if (*colon_pos == '"') {
        colon_pos++;
        char* end_quote = strchr(colon_pos, '"');
        if (!end_quote) return -1;
        int len = end_quote - colon_pos;
        if (len >= value_size) len = value_size - 1;
        strncpy(value, colon_pos, len);
        value[len] = '\0';
    } else if (isdigit(*colon_pos) || *colon_pos == '-') {
        char* end_num = colon_pos;
        while (*end_num == '-' || *end_num == '.' || isdigit(*end_num)) end_num++;
        int len = end_num - colon_pos;
        if (len >= value_size) len = value_size - 1;
        strncpy(value, colon_pos, len);
        value[len] = '\0';
    } else {
        return -1;
    }
    
    return 0;
}

/* Load configuration from JSON file */
int load_config_from_json(const char* filename, LindbladConfig* config) {
    FILE* file = fopen(filename, "r");
    if (!file) return -1;
    
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    fseek(file, 0, SEEK_SET);
    
    char* json_buffer = (char*)malloc(file_size + 1);
    if (!json_buffer) {
        fclose(file);
        return -1;
    }
    
    fread(json_buffer, 1, file_size, file);
    json_buffer[file_size] = '\0';
    fclose(file);
    
    memset(config, 0, sizeof(LindbladConfig));
    
    /* Parse num_qubits */
    char num_str[32];
    if (parse_json_value(json_buffer, "num_qubits", num_str, sizeof(num_str)) == 0) {
        config->num_qubits = atoi(num_str);
    } else {
        free(json_buffer);
        return -1;
    }
    
    /* Initialize all qubits to default state */
    init_all_qubits_default(config);
    
    /* Parse qubits array if exists */
    char* qubits_start = strstr(json_buffer, "\"qubits\"");
    if (qubits_start) {
        char* array_start = strchr(qubits_start, '[');
        if (array_start) {
            int qubit_count = 0;
            char* current = array_start + 1;
            
            while (qubit_count < config->num_qubits && current) {
                char* object_start = strchr(current, '{');
                if (!object_start) break;
                
                char* object_end = strchr(object_start, '}');
                if (!object_end) break;
                
                int obj_len = object_end - object_start - 1;
                char obj_buffer[MAX_STRING_LEN];
                if (obj_len >= MAX_STRING_LEN) obj_len = MAX_STRING_LEN - 1;
                strncpy(obj_buffer, object_start + 1, obj_len);
                obj_buffer[obj_len] = '\0';
                
                char omega_str[32], E_str[32];
                if (parse_json_value(obj_buffer, "omega", omega_str, sizeof(omega_str)) == 0) {
                    config->qubits[qubit_count].omega = atof(omega_str);
                }
                if (parse_json_value(obj_buffer, "E", E_str, sizeof(E_str)) == 0) {
                    config->qubits[qubit_count].E = atof(E_str);
                }
                
                qubit_count++;
                current = object_end + 1;
            }
        }
    }
    
    /* Parse disturbances array */
    char* dist_start = strstr(json_buffer, "\"disturbances\"");
    if (dist_start) {
        char* dist_array_start = strchr(dist_start, '[');
        if (dist_array_start) {
            config->num_disturbances = 0;
            char* dist_current = dist_array_start + 1;
            
            while (config->num_disturbances < MAX_DISTURBANCES && dist_current) {
                while (*dist_current == ' ' || *dist_current == '\t' || 
                       *dist_current == ',' || *dist_current == '\n') 
                    dist_current++;
                
                if (*dist_current == ']') break;
                
                if (*dist_current == '"') {
                    dist_current++;
                    char* end_quote = strchr(dist_current, '"');
                    if (!end_quote) break;
                    
                    int len = end_quote - dist_current;
                    if (len >= MAX_STRING_LEN) len = MAX_STRING_LEN - 1;
                    
                    char dist_str[MAX_STRING_LEN];
                    strncpy(dist_str, dist_current, len);
                    dist_str[len] = '\0';
                    
                    if (parse_disturbance_string(dist_str, 
                        &config->disturbances[config->num_disturbances]) == 0) {
                        config->num_disturbances++;
                    }
                    
                    dist_current = end_quote + 1;
                } else {
                    dist_current++;
                }
            }
        }
    }
    
    free(json_buffer);
    return 0;
}

/* Load configuration from simple text file */
int load_config_from_file(const char* filename, LindbladConfig* config) {
    FILE* file = fopen(filename, "r");
    if (!file) return -1;
    
    memset(config, 0, sizeof(LindbladConfig));
    
    char line[MAX_STRING_LEN];
    int line_num = 0;
    
    while (fgets(line, sizeof(line), file)) {
        line_num++;
        
        char* newline = strchr(line, '\n');
        if (newline) *newline = '\0';
        if (line[0] == '\0' || line[0] == '#') continue;
        
        if (line_num == 1) {
            config->num_qubits = atoi(line);
            if (config->num_qubits <= 0 || config->num_qubits > MAX_QUBITS) {
                fclose(file);
                return -1;
            }
            /* Initialize all qubits to default state */
            init_all_qubits_default(config);
        } else if (line_num <= config->num_qubits + 1) {
            int qubit_idx = line_num - 2;
            if (sscanf(line, "%f %f", 
                       &config->qubits[qubit_idx].omega, 
                       &config->qubits[qubit_idx].E) == 2) {
                /* Keep default rho state */
            }
        } else {
            if (config->num_disturbances < MAX_DISTURBANCES) {
                if (parse_disturbance_string(line, 
                    &config->disturbances[config->num_disturbances]) == 0) {
                    config->num_disturbances++;
                }
            }
        }
    }
    
    fclose(file);
    return 0;
}

/* Print configuration */
/* Print configuration */
void print_config(const LindbladConfig* config) {
    printf("\n========== LINDBLAD CONFIGURATION ==========\n");
    printf("Number of qubits: %d\n\n", config->num_qubits);
    
    printf("Qubit Parameters and Initial States:\n");
    for (int i = 0; i < config->num_qubits; i++) {
        printf("\nQubit %d:\n", i);
        printf("  ω = %.4f, E = %.4f\n", config->qubits[i].omega, config->qubits[i].E);
        printf("  Initial ρ = [[%.1f, %.1f], [%.1f, %.1f]]\n",
               config->qubits[i].rho[0].real, config->qubits[i].rho[1].real,
               config->qubits[i].rho[2].real, config->qubits[i].rho[3].real);
    }
    
    printf("\nDisturbances (%d):\n", config->num_disturbances);
    for (int i = 0; i < config->num_disturbances; i++) {
        /* Remove const by using a copy or just access directly */
        printf("  %d: Qubit %d, Type: %s, γ = %.4f [%s]\n",
               i, config->disturbances[i].qubit_index, 
               disturbance_type_to_string(config->disturbances[i].type),
               config->disturbances[i].gamma, config->disturbances[i].raw_string);
    }
    printf("============================================\n");
}

/* Parse custom initial state from JSON */
static int parse_initial_state(const char* json, QubitInfo* qubit) {
    /* Look for "initial_state" or "rho" key */
    char* state_start = strstr(json, "\"initial_state\"");
    if (!state_start) {
        state_start = strstr(json, "\"rho\"");
    }
    if (!state_start) return -1;
    
    char* array_start = strchr(state_start, '[');
    if (!array_start) return -1;
    
    /* Parse [[a,b],[c,d]] format */
    float values[4] = {0};
    int count = 0;
    char* current = array_start;
    
    while (count < 4 && current) {
        current = strchr(current, '[');
        if (!current) break;
        current++;
        
        /* Parse two numbers */
        for (int i = 0; i < 2 && count < 4; i++) {
            char* endptr;
            values[count] = strtof(current, &endptr);
            if (endptr == current) break;
            current = endptr;
            count++;
            
            /* Skip comma */
            current = strchr(current, ',');
            if (current) current++;
        }
        
        /* Skip to next row */
        current = strchr(current, ']');
        if (current) current++;
    }
    
    if (count == 4) {
        qubit->rho[0].real = values[0]; qubit->rho[0].imag = 0;
        qubit->rho[1].real = values[1]; qubit->rho[1].imag = 0;
        qubit->rho[2].real = values[2]; qubit->rho[2].imag = 0;
        qubit->rho[3].real = values[3]; qubit->rho[3].imag = 0;
        return 0;
    }
    
    return -1;
}
