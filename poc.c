#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// ========== CONFIGURATION ==========
#define NUM_CONSTRAINTS 11

// Constants
#define M 25214903917ULL
#define A 11ULL
#define P 341873128712ULL
#define Q 132897987541ULL
#define MOD_48 (1ULL << 48)
#define MOD_35 (1ULL << 35)
#define MOD_31 (1ULL << 31)
#define MOD_18 (1ULL << 18)
#define MOD_17 (1ULL << 17)
#define MOD_13 (1ULL << 13)

const int32_t spacing = 34;
const int32_t separation = 12;
const int32_t bound = 22;
const uint64_t salt = 94251327;

// Expected values
const uint64_t EXPECTED_W = 203501504278641ULL;
const uint64_t EXPECTED_W_H = EXPECTED_W >> 35;
const uint64_t EXPECTED_W_L = EXPECTED_W % MOD_35;

// ========== DATA STRUCTURES ==========
typedef struct {
    int32_t x_r, z_r;
    int32_t x_c, z_c;
    uint64_t D;
    uint64_t D_h;
    int32_t dx;
    int32_t dz;
} Constraint;

// ========== HELPER FUNCTIONS ==========
void compute_constraint_values(Constraint* c) {
    uint64_t D = (P * c->x_r + Q * c->z_r + salt) % MOD_48;
    c->D = D;
    c->D_h = D >> 35;
    c->dx = c->x_c - c->x_r * spacing;
    c->dz = c->z_c - c->z_r * spacing;
}

// From theory: w_l = (e0 xor M) - (D mod 2^35) mod 2^35
uint64_t compute_w_l_from_e0(uint64_t e0, uint64_t D) {
    uint64_t D_low = D & (MOD_35 - 1);
    uint64_t w_l = (e0 ^ M) - D_low;
    w_l &= (MOD_35 - 1);
    return w_l;
}

// From theory: e0 = ((w_l + D) mod 2^35) xor M
uint64_t compute_e0_from_w_l(uint64_t w_l, uint64_t D) {
    uint64_t D_low = D & (MOD_35 - 1);
    uint64_t sum = (w_l + D_low) & (MOD_35 - 1);
    return sum ^ M;
}

int32_t gcd(int32_t a, int32_t b) {
    while (b) {
        int32_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// Weak condition 1: (dx - ((e0 * M + A) >> 17)) % g == 0
int check_weak1(uint64_t e0, int32_t dx, int32_t g) {
    if (g == 0) return 1;
    uint64_t term = (e0 * M + A) % MOD_48;
    uint64_t shifted = term >> 17;
    int32_t result = (dx - (int32_t)shifted) % g;
    if (result < 0) result += g;
    return result == 0;
}

// Weak condition 2: (dz - ((e0 * M^2 + A*(M+1)) >> 17)) % g == 0
int check_weak2(uint64_t e0, int32_t dz, int32_t g) {
    if (g == 0) return 1;
    uint64_t M_sq = (M * M) % MOD_48;
    uint64_t term = (e0 * M_sq + A * (M + 1)) % MOD_48;
    uint64_t shifted = term >> 17;
    int32_t result = (dz - (int32_t)shifted) % g;
    if (result < 0) result += g;
    return result == 0;
}

// Compute X from theory: X = (dx - (B >> 17) + c1*2^31) mod bound
int32_t compute_X(uint64_t e0, uint64_t D_h, int c0, int c1, int32_t dx) {
    uint64_t X_val = (MOD_35 * (D_h + c0) + e0) % MOD_48;
    uint64_t B = (X_val * M + A) % MOD_48;
    uint64_t shifted = B >> 17;
    
    int64_t X = (int64_t)dx - (int64_t)shifted + (int64_t)(c1 * MOD_31);
    X %= bound;
    if (X < 0) X += bound;
    
    return (int32_t)X;
}

// Check if w_h produces X: (2^18 * M * w_h) mod 2^31 mod bound == X
int check_w_h_for_X(uint32_t w_h, int32_t X) {
    uint64_t term = MOD_18 * M * w_h;
    uint64_t mod_31 = term % MOD_31;
    int32_t result = mod_31 % bound;
    if (result < 0) result += bound;
    return result == X;
}

// Simple verification
int verify_w_simple(uint64_t w, const Constraint* constraints, int num_constraints) {
    for (int i = 0; i < num_constraints; i++) {
        const Constraint* c = &constraints[i];
        
        uint64_t seed = (w + P * c->x_r + Q * c->z_r + salt) % MOD_48;
        seed = seed ^ M;
        
        // First nextInt(bound)
        uint64_t s1 = (seed * M + A) % MOD_48;
        uint32_t bits1 = (uint32_t)(s1 >> 17);
        uint32_t value1 = bits1 % bound;
        
        if (value1 != (uint32_t)c->dx) return 0;
        
        // Second nextInt(bound)
        uint64_t s2 = (s1 * M + A) % MOD_48;
        uint32_t bits2 = (uint32_t)(s2 >> 17);
        uint32_t value2 = bits2 % bound;
        
        if (value2 != (uint32_t)c->dz) return 0;
    }
    return 1;
}

int64_t mod_inverse(int64_t a, int64_t m) {
    int64_t m0 = m, y = 0, x = 1;
    if (m == 1) return 0;
    
    while (a > 1) {
        int64_t q = a / m;
        int64_t t = m;
        m = a % m;
        a = t;
        t = y;
        y = x - q * y;
        x = t;
    }
    
    if (x < 0) x += m0;
    return x;
}

// ========== LOGGING FUNCTIONS ==========
FILE* log_file = NULL;

void open_log_file(int num_constraints) {
    log_file = fopen("found_seeds_detailed.log", "w");
    if (!log_file) {
        printf("ERROR: Could not open log file\n");
        return;
    }
    fprintf(log_file, "=== Seed Recovery Log ===\n");
    fprintf(log_file, "Started: %s", ctime(&(time_t){time(NULL)}));
    fprintf(log_file, "Searching for seeds matching %d constraints\n", num_constraints);
    fprintf(log_file, "Expected seed: %llu\n", (unsigned long long)EXPECTED_W);
    fprintf(log_file, "---\n");
    fflush(log_file);
}

void log_seed(uint64_t w, uint64_t w_h, uint64_t w_l, int32_t common_X, 
              uint64_t total_candidates, uint64_t weak_passed, 
              uint64_t strong_checked, double elapsed, int is_expected) {
    if (!log_file) return;
    
    fprintf(log_file, "=== FOUND SEED ===\n");
    fprintf(log_file, "w = %llu (0x%llx)\n", (unsigned long long)w, (unsigned long long)w);
    fprintf(log_file, "w_h = %llu, w_l = %llu\n", 
            (unsigned long long)w_h, (unsigned long long)w_l);
    fprintf(log_file, "Common X = %d\n", common_X);
    fprintf(log_file, "Matches expected: %s\n", is_expected ? "YES" : "NO");
    fprintf(log_file, "Search stats at find time:\n");
    fprintf(log_file, "  Total candidates: %llu\n", (unsigned long long)total_candidates);
    fprintf(log_file, "  Weak passed: %llu\n", (unsigned long long)weak_passed);
    fprintf(log_file, "  Strong checked: %llu\n", (unsigned long long)strong_checked);
    fprintf(log_file, "  Time: %.2f seconds\n", elapsed);
    fprintf(log_file, "---\n");
    fflush(log_file);
}

void close_log_file() {
    if (log_file) {
        fprintf(log_file, "\n=== SEARCH COMPLETE ===\n");
        fprintf(log_file, "Ended: %s", ctime(&(time_t){time(NULL)}));
        fclose(log_file);
    }
}

// ========== MAIN ALGORITHM ==========
void find_seeds(const Constraint* constraints, int num_constraints) {
    printf("Number of constraints: %d\n", num_constraints);
    
    // Open log file
    open_log_file(num_constraints);
    
    int32_t g = gcd(MOD_18, bound);
    uint64_t m = MOD_17 * g;
    
    printf("g = %d, m = %llu\n", g, (unsigned long long)m);
    printf("Expected w: %llu (0x%llx)\n", 
           (unsigned long long)EXPECTED_W, (unsigned long long)EXPECTED_W);
    printf("Expected w_h: %llu, w_l: %llu\n\n",
           (unsigned long long)EXPECTED_W_H, (unsigned long long)EXPECTED_W_L);
    
    // Get modular inverse
    uint64_t M_mod_m = M % m;
    int64_t M_inv = mod_inverse(M_mod_m, m);
    if (M_inv == 0) {
        printf("ERROR: No modular inverse\n");
        close_log_file();
        return;
    }
    
    clock_t start = clock();
    uint64_t total_candidates = 0;
    uint64_t weak_passed = 0;
    uint64_t strong_checked = 0;
    uint64_t found_count = 0;
    int expected_seed_found = 0;
    
    const Constraint* c0 = &constraints[0];
    
    // Precompute w_h -> X mapping
    uint8_t* w_h_to_X = (uint8_t*)malloc(MOD_13 * sizeof(uint8_t));
    for (uint32_t w_h = 0; w_h < MOD_13; w_h++) {
        uint64_t term = MOD_18 * M * w_h;
        uint64_t mod_31 = term % MOD_31;
        int32_t X = mod_31 % bound;
        if (X < 0) X += bound;
        w_h_to_X[w_h] = (uint8_t)X;
    }
    
    // Allocate arrays dynamically based on num_constraints
    uint64_t* e0_list = (uint64_t*)malloc(num_constraints * sizeof(uint64_t));
    uint8_t** possible_X = (uint8_t**)malloc(num_constraints * sizeof(uint8_t*));
    for (int i = 0; i < num_constraints; i++) {
        possible_X[i] = (uint8_t*)malloc(bound * sizeof(uint8_t));
    }
    uint8_t* common_X = (uint8_t*)malloc(bound * sizeof(uint8_t));
    
    // For progress
    int last_percent = -1;
    uint64_t last_time = 0;
    
    // BRUTE FORCE r (0 to 2^17-1)
    for (uint32_t r = 0; r < MOD_17; r++) {
        int percent = (int)((double)r / MOD_17 * 100);
        clock_t now = clock();
        double elapsed = (double)(now - start) / CLOCKS_PER_SEC;
        
        if (percent != last_percent || elapsed - last_time >= 1.0) {
            printf("\rProgress: %d%% | Cand: %llu | Weak: %llu | Strong: %llu | Found: %llu | Time: %.1fs",
                   percent, (unsigned long long)total_candidates,
                   (unsigned long long)weak_passed,
                   (unsigned long long)strong_checked,
                   (unsigned long long)found_count,
                   elapsed);
            fflush(stdout);
            last_percent = percent;
            last_time = elapsed;
        }
        
        // Compute base e0
        int64_t base = (int64_t)(MOD_17 * c0->dx) - A + (int64_t)r;
        base %= (int64_t)m;
        if (base < 0) base += m;
        
        uint64_t e0_base = ((uint64_t)base * M_inv) % m;
        uint32_t q_max = (MOD_35 - e0_base + m - 1) / m;
        
        for (uint32_t q = 0; q < q_max; q++) {
            total_candidates++;
            
            uint64_t e0_0 = e0_base + q * m;
            if (e0_0 >= MOD_35) continue;
            
            // Check weak condition 2 for first constraint
            if (!check_weak2(e0_0, c0->dz, g)) {
                continue;
            }
            
            // Compute w_l from this e0
            uint64_t w_l = compute_w_l_from_e0(e0_0, c0->D);
            
            // DEBUG: Only for expected w_l
            int is_expected_w_l = (w_l == EXPECTED_W_L);
            if (is_expected_w_l) {
                printf("\n\n*** DEBUG: Found expected w_l = %llu (0x%llx) ***\n", 
                       (unsigned long long)w_l, (unsigned long long)w_l);
                printf("  e0_0 = %llu, r = %u, q = %u\n",
                       (unsigned long long)e0_0, r, q);
            }
            
            // Compute e0 for all constraints from this w_l
            e0_list[0] = e0_0;
            int all_weak_ok = 1;
            
            for (int i = 1; i < num_constraints; i++) {
                e0_list[i] = compute_e0_from_w_l(w_l, constraints[i].D);
                
                if (!check_weak1(e0_list[i], constraints[i].dx, g) ||
                    !check_weak2(e0_list[i], constraints[i].dz, g)) {
                    all_weak_ok = 0;
                    break;
                }
            }
            
            if (!all_weak_ok) {
                continue;
            }
            
            weak_passed++;
            
            // STRONG CONDITION: Find common X values across all constraints
            // For each constraint, compute all possible X values (for all c0,c1)
            // Then find X values common to ALL constraints
            
            for (int i = 0; i < num_constraints; i++) {
                memset(possible_X[i], 0, bound * sizeof(uint8_t));
                
                // Mark all possible X values for this constraint
                for (int c0_bit = 0; c0_bit < 2; c0_bit++) {
                    for (int c1_bit = 0; c1_bit < 2; c1_bit++) {
                        int32_t X = compute_X(e0_list[i], constraints[i].D_h,
                                             c0_bit, c1_bit, constraints[i].dx);
                        possible_X[i][X] = 1;
                    }
                }
            }
            
            // Find X values common to ALL constraints
            memset(common_X, 1, bound * sizeof(uint8_t));  // Start with all true
            int num_common_X = 0;
            
            for (int x = 0; x < bound; x++) {
                for (int i = 0; i < num_constraints; i++) {
                    if (!possible_X[i][x]) {
                        common_X[x] = 0;
                        break;
                    }
                }
                if (common_X[x]) {
                    num_common_X++;
                }
            }
            
            // If no common X values, skip
            if (num_common_X == 0) {
                continue;
            }
            
            // For each common X value, find w_h candidates
            for (int X = 0; X < bound; X++) {
                if (!common_X[X]) {
                    continue;
                }
                
                strong_checked++;  // We're checking this X value
                
                // Find all w_h that produce this X
                for (uint32_t w_h = 0; w_h < MOD_13; w_h++) {
                    if (w_h_to_X[w_h] != X) {
                        continue;
                    }
                    
                    // DEBUG: Only for expected w_h with expected w_l
                    if (is_expected_w_l && w_h == EXPECTED_W_H) {
                        printf("  *** DEBUG: Checking expected w_h = %u (0x%x) for common X=%d ***\n",
                               w_h, w_h, X);
                    }
                    
                    // Build and verify the seed
                    uint64_t w = ((uint64_t)w_h << 35) | w_l;
                    
                    if (verify_w_simple(w, constraints, num_constraints)) {
                        found_count++;
                        
                        double elapsed_total = (double)(clock() - start) / CLOCKS_PER_SEC;
                        
                        int is_expected = (w == EXPECTED_W);
                        
                        // Console output - detailed for expected seed only
                        if (is_expected) {
                            printf("\n\n*** EXPECTED SEED FOUND! ***\n");
                            printf("w = %llu (0x%llx)\n",
                                   (unsigned long long)w, (unsigned long long)w);
                            printf("w_h = %u, w_l = %llu\n",
                                   w_h, (unsigned long long)w_l);
                            printf("Common X = %d\n", X);
                            printf("Found at: %llu candidates, %.2f seconds\n",
                                   (unsigned long long)total_candidates, elapsed_total);
                            expected_seed_found = 1;
                        }
                        
                        // Log to file (all seeds)
                        log_seed(w, w_h, w_l, X, total_candidates, weak_passed, 
                                strong_checked, elapsed_total, is_expected);
                        
                        // Continue searching (don't return)
                    }
                }
            }
        }
    }
    
    // Cleanup
    free(w_h_to_X);
    free(e0_list);
    for (int i = 0; i < num_constraints; i++) {
        free(possible_X[i]);
    }
    free(possible_X);
    free(common_X);
    
    double elapsed_total = (double)(clock() - start) / CLOCKS_PER_SEC;
    printf("\n\n=== SEARCH COMPLETE ===\n");
    printf("Number of constraints: %d\n", num_constraints);
    printf("Total candidates: %llu\n", (unsigned long long)total_candidates);
    printf("Weak passed: %llu (%.4f%%)\n", 
           (unsigned long long)weak_passed,
           (double)weak_passed / total_candidates * 100);
    printf("Strong checks: %llu (%.4f%% of weak)\n",
           (unsigned long long)strong_checked,
           (double)strong_checked / weak_passed * 100);
    printf("Total seeds found: %llu\n", (unsigned long long)found_count);
    printf("Expected seed found: %s\n", expected_seed_found ? "YES" : "NO");
    printf("Total time: %.2f seconds\n", elapsed_total);
    
    if (found_count == 0) {
        printf("ERROR: No seeds found!\n");
    } else if (!expected_seed_found) {
        printf("WARNING: Expected seed not found among %llu seeds!\n",
               (unsigned long long)found_count);
    }
    
    close_log_file();
}

// ========== MAIN FUNCTION ==========
int main() {
    // Setup constraints - Update this array with your constraints
    // Array size must match NUM_CONSTRAINTS define
    Constraint constraints[NUM_CONSTRAINTS];
    
    constraints[0].x_r = 1; constraints[0].z_r = 1;
    constraints[0].x_c = 41; constraints[0].z_c = 50;
    
    constraints[1].x_r = 2; constraints[1].z_r = 1;
    constraints[1].x_c = 87; constraints[1].z_c = 44;
    
    constraints[2].x_r = 3; constraints[2].z_r = 1;
    constraints[2].x_c = 118; constraints[2].z_c = 50;
    
    constraints[3].x_r = 4; constraints[3].z_r = 1;
    constraints[3].x_c = 137; constraints[3].z_c = 53;
    
    constraints[4].x_r = 5; constraints[4].z_r = 1;
    constraints[4].x_c = 178; constraints[4].z_c = 48;
    
    constraints[5].x_r = 6; constraints[5].z_r = 1;
    constraints[5].x_c = 204; constraints[5].z_c = 46;
    
    constraints[6].x_r = 7; constraints[6].z_r = 1;
    constraints[6].x_c = 247; constraints[6].z_c = 48;
    
    constraints[7].x_r = 8; constraints[7].z_r = 1;
    constraints[7].x_c = 288; constraints[7].z_c = 50;
    
    constraints[8].x_r = 9; constraints[8].z_r = 1;
    constraints[8].x_c = 309; constraints[8].z_c = 53;
    
    constraints[9].x_r = 10; constraints[9].z_r = 1;
    constraints[9].x_c = 360; constraints[9].z_c = 45;
    
    constraints[10].x_r = 11; constraints[10].z_r = 1;
    constraints[10].x_c = 382; constraints[10].z_c = 37;
    
    // Initialize all constraints
    for (int i = 0; i < NUM_CONSTRAINTS; i++) {
        compute_constraint_values(&constraints[i]);
    }
    
    printf("Constraints (%d total):\n", NUM_CONSTRAINTS);
    for (int i = 0; i < NUM_CONSTRAINTS; i++) {
        printf("  [%d] (x_r=%d, z_r=%d) -> (x_c=%d, z_c=%d)\n",
               i, constraints[i].x_r, constraints[i].z_r,
               constraints[i].x_c, constraints[i].z_c);
        printf("       dx=%d, dz=%d, D_h=%llu\n",
               constraints[i].dx, constraints[i].dz,
               (unsigned long long)constraints[i].D_h);
    }
    printf("\n");
    
    // Run the search
    find_seeds(constraints, NUM_CONSTRAINTS);
    
    return 0;
}
