#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* The order-five optimum with negative cycle edges. */
static const int BASE[5][5] = {
    {0, -1, 1, 1, -1},
    {-1, 0, -1, 1, 1},
    {1, -1, 0, -1, 1},
    {1, 1, -1, 0, -1},
    {-1, 1, 1, -1, 0},
};

/* Representatives for the trace-zero and absolute-trace-four classes. */
static const int HADAMARD_REPRESENTATIVES[2][4][4] = {
    {
        {1, 1, 1, 1},
        {1, -1, 1, -1},
        {1, 1, -1, -1},
        {1, -1, -1, 1},
    },
    {
        {1, -1, -1, -1},
        {-1, 1, -1, -1},
        {-1, -1, 1, -1},
        {-1, -1, -1, 1},
    },
};

static const int EDGES[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
};

static int sign_from_bit(unsigned int mask, int bit) {
    return ((mask >> (unsigned int)bit) & 1U) != 0U ? -1 : 1;
}

static int is_hadamard(const int matrix[4][4]) {
    for (int row = 0; row < 4; ++row) {
        for (int other = 0; other < 4; ++other) {
            int inner_product = 0;
            for (int column = 0; column < 4; ++column) {
                inner_product += matrix[row][column] * matrix[other][column];
            }
            if (inner_product != (row == other ? 4 : 0)) {
                return 0;
            }
        }
    }
    return 1;
}

static int matrices_equal(const int left[4][4], const int right[4][4]) {
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            if (left[row][column] != right[row][column]) {
                return 0;
            }
        }
    }
    return 1;
}

static void verify_base_global_sign_reduction(void) {
    /* This permutation sends the signed five-cycle base to its negative.
       Therefore the global sign used in the Hadamard classification can be
       absorbed by relabelling the five clouds. */
    static const int permutation[5] = {0, 2, 4, 1, 3};
    for (int row = 0; row < 5; ++row) {
        for (int column = 0; column < 5; ++column) {
            if (BASE[permutation[row]][permutation[column]] !=
                -BASE[row][column]) {
                fprintf(stderr, "base anti-isomorphism mismatch\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    printf("base_global_sign_reduction=0,2,4,1,3\n");
}

static int next_permutation(int permutation[4]) {
    int pivot = 2;
    while (pivot >= 0 && permutation[pivot] > permutation[pivot + 1]) {
        --pivot;
    }
    if (pivot < 0) {
        return 0;
    }
    int successor = 3;
    while (permutation[successor] < permutation[pivot]) {
        --successor;
    }
    int temporary = permutation[pivot];
    permutation[pivot] = permutation[successor];
    permutation[successor] = temporary;
    for (int left = pivot + 1, right = 3; left < right; ++left, --right) {
        temporary = permutation[left];
        permutation[left] = permutation[right];
        permutation[right] = temporary;
    }
    return 1;
}

static void verify_hadamard_classification(void) {
    int total = 0;
    int trace_counts[3] = {0, 0, 0};

    for (unsigned int bits = 0; bits < (1U << 10U); ++bits) {
        int matrix[4][4];
        int bit = 0;
        for (int row = 0; row < 4; ++row) {
            for (int column = row; column < 4; ++column) {
                matrix[row][column] = sign_from_bit(bits, bit);
                matrix[column][row] = matrix[row][column];
                ++bit;
            }
        }
        if (!is_hadamard(matrix)) {
            continue;
        }

        ++total;
        int trace = 0;
        for (int row = 0; row < 4; ++row) {
            trace += matrix[row][row];
        }
        if (trace != -4 && trace != 0 && trace != 4) {
            fprintf(stderr, "unexpected symmetric-Hadamard trace\n");
            exit(EXIT_FAILURE);
        }
        ++trace_counts[(trace + 4) / 4];

        int covered = 0;
        int permutation[4] = {0, 1, 2, 3};
        do {
            for (unsigned int diagonal = 0; diagonal < 16U && !covered;
                 ++diagonal) {
                for (int representative = 0; representative < 2 && !covered;
                     ++representative) {
                    for (int global_sign = -1; global_sign <= 1;
                         global_sign += 2) {
                        int candidate[4][4];
                        for (int row = 0; row < 4; ++row) {
                            for (int column = 0; column < 4; ++column) {
                                candidate[row][column] =
                                    global_sign *
                                    sign_from_bit(diagonal, row) *
                                    sign_from_bit(diagonal, column) *
                                    HADAMARD_REPRESENTATIVES[representative]
                                                                [permutation[row]]
                                                                [permutation[column]];
                            }
                        }
                        if (matrices_equal(matrix, candidate)) {
                            covered = 1;
                            break;
                        }
                    }
                }
            }
        } while (!covered && next_permutation(permutation));

        if (!covered) {
            fprintf(stderr, "unclassified symmetric Hadamard matrix\n");
            exit(EXIT_FAILURE);
        }
    }

    if (total != 64 || trace_counts[0] != 8 || trace_counts[1] != 48 ||
        trace_counts[2] != 8) {
        fprintf(stderr, "symmetric-Hadamard classification count mismatch\n");
        exit(EXIT_FAILURE);
    }
    int corrupted[4][4];
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            corrupted[row][column] = HADAMARD_REPRESENTATIVES[0][row][column];
        }
    }
    corrupted[0][0] *= -1;
    if (is_hadamard(corrupted)) {
        fprintf(stderr, "Hadamard-entry corruption was not detected\n");
        exit(EXIT_FAILURE);
    }
    printf("symmetric_hadamards=64\n");
    printf("trace_counts=-4:8,0:48,4:8\n");
}

static int cloud_energy(const int hadamard[4][4], unsigned int diagonal_bits,
                        const int cloud_masks[5]) {
    int energy = 0;
    for (int cloud = 0; cloud < 5; ++cloud) {
        for (int edge = 0; edge < 6; ++edge) {
            energy += sign_from_bit(diagonal_bits, edge) *
                      sign_from_bit((unsigned int)cloud_masks[cloud],
                                    EDGES[edge][0]) *
                      sign_from_bit((unsigned int)cloud_masks[cloud],
                                    EDGES[edge][1]);
        }
    }
    for (int left_cloud = 0; left_cloud < 5; ++left_cloud) {
        for (int right_cloud = left_cloud + 1; right_cloud < 5;
             ++right_cloud) {
            int block_energy = 0;
            for (int row = 0; row < 4; ++row) {
                for (int column = 0; column < 4; ++column) {
                    block_energy +=
                        sign_from_bit((unsigned int)cloud_masks[left_cloud], row) *
                        hadamard[row][column] *
                        sign_from_bit((unsigned int)cloud_masks[right_cloud],
                                      column);
                }
            }
            energy += BASE[left_cloud][right_cloud] * block_energy;
        }
    }
    return energy;
}

static int full_maximum(const int hadamard[4][4], unsigned int diagonal_bits,
                        uint32_t *maximizing_index) {
    int best = 0;
    uint32_t index = 0U;
    uint32_t best_index = 0U;

    /* The fourth coordinate of cloud zero is fixed positive.  This removes
       only the global spin symmetry and leaves exactly 2^19 states. */
    for (int first = 0; first < 8; ++first) {
        for (int second = 0; second < 16; ++second) {
            for (int third = 0; third < 16; ++third) {
                for (int fourth = 0; fourth < 16; ++fourth) {
                    for (int fifth = 0; fifth < 16; ++fifth, ++index) {
                        const int masks[5] = {
                            first, second, third, fourth, fifth,
                        };
                        const int energy =
                            cloud_energy(hadamard, diagonal_bits, masks);
                        const int absolute_energy = energy < 0 ? -energy : energy;
                        if (absolute_energy > best) {
                            best = absolute_energy;
                            best_index = index;
                        }
                    }
                }
            }
        }
    }
    if (index != (1U << 19U)) {
        fprintf(stderr, "incorrect projective spin count\n");
        exit(EXIT_FAILURE);
    }
    if (maximizing_index != NULL) {
        *maximizing_index = best_index;
    }
    return best;
}

static uint64_t fnv1a_u32(uint64_t digest, uint32_t value) {
    for (int byte = 0; byte < 4; ++byte) {
        digest ^= (uint64_t)((value >> (unsigned int)(8 * byte)) & 0xffU);
        digest *= UINT64_C(1099511628211);
    }
    return digest;
}

static int product_spin_maximum(const int hadamard[4][4],
                                unsigned int diagonal_bits) {
    int best = 0;
    for (unsigned int base_mask = 0; base_mask < 32U; ++base_mask) {
        for (unsigned int cloud_mask = 0; cloud_mask < 16U; ++cloud_mask) {
            int masks[5];
            for (int cloud = 0; cloud < 5; ++cloud) {
                masks[cloud] =
                    sign_from_bit(base_mask, cloud) == 1
                        ? (int)cloud_mask
                        : (int)(cloud_mask ^ 15U);
            }
            const int energy = cloud_energy(hadamard, diagonal_bits, masks);
            const int absolute_energy = energy < 0 ? -energy : energy;
            if (absolute_energy > best) {
                best = absolute_energy;
            }
        }
    }
    return best;
}

static unsigned int transformed_diagonal(const int switching[4],
                                         const int diagonal_signs[6]) {
    unsigned int result = 0U;
    for (int edge = 0; edge < 6; ++edge) {
        const int sign = diagonal_signs[edge] * switching[EDGES[edge][0]] *
                         switching[EDGES[edge][1]];
        if (sign == -1) {
            result |= 1U << (unsigned int)edge;
        }
    }
    return result;
}

static void verify_fixed_half_attainment(void) {
    const int diagonal_signs[6] = {-1, 1, -1, -1, -1, 1};
    const int switching[4] = {-1, -1, -1, 1};
    int hadamard[4][4];
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            hadamard[row][column] = switching[row] * switching[column] *
                                    HADAMARD_REPRESENTATIVES[0][row][column];
        }
    }
    const unsigned int diagonal_bits =
        transformed_diagonal(switching, diagonal_signs);
    int diagonal_total = 0;
    for (int edge = 0; edge < 6; ++edge) {
        diagonal_total += sign_from_bit(diagonal_bits, edge);
    }
    int base_total = 0;
    for (int row = 0; row < 5; ++row) {
        for (int column = row + 1; column < 5; ++column) {
            base_total += BASE[row][column];
        }
    }
    int hadamard_total = 0;
    for (int row = 0; row < 4; ++row) {
        for (int column = 0; column < 4; ++column) {
            hadamard_total += hadamard[row][column];
        }
    }
    const int total_sign = base_total * hadamard_total + 5 * diagonal_total;
    if (!is_hadamard(hadamard) || diagonal_total != 0 || total_sign != 0) {
        fprintf(stderr, "fixed-half construction is corrupted\n");
        exit(EXIT_FAILURE);
    }

    uint32_t maximizing_index = 0U;
    const int maximum = full_maximum(hadamard, diagonal_bits, &maximizing_index);
    const int product_maximum =
        product_spin_maximum(hadamard, diagonal_bits);
    if (maximum != 44 || product_maximum != 32) {
        fprintf(stderr, "fixed-half attainment or fine-spin gap mismatch\n");
        exit(EXIT_FAILURE);
    }

    const int witness_masks[5] = {9, 10, 6, 4, 13};
    const int witness_energy =
        cloud_energy(hadamard, diagonal_bits, witness_masks);
    int corrupted_masks[5] = {9, 10, 6, 4, 13};
    corrupted_masks[0] ^= 1;
    const int corrupted_energy =
        cloud_energy(hadamard, diagonal_bits, corrupted_masks);
    if (witness_energy != 44 || corrupted_energy == witness_energy ||
        maximizing_index >= (1U << 19U)) {
        fprintf(stderr, "fine-spin witness corruption was not detected\n");
        exit(EXIT_FAILURE);
    }

    printf("fixed_half_total_sign=0\n");
    printf("fixed_half_minimum_attained=44\n");
    printf("product_spin_maximum=32\n");
    printf("fine_spin_witness=44\n");
}

int main(void) {
    verify_hadamard_classification();
    verify_base_global_sign_reduction();

    int minima[2] = {1000, 1000};
    int minimum_counts[2] = {0, 0};
    uint64_t digest = UINT64_C(14695981039346656037);
    for (int representative = 0; representative < 2; ++representative) {
        for (unsigned int diagonal_bits = 0; diagonal_bits < 64U;
             ++diagonal_bits) {
            const int maximum = full_maximum(
                HADAMARD_REPRESENTATIVES[representative], diagonal_bits, NULL);
            digest = fnv1a_u32(digest, (uint32_t)representative);
            digest = fnv1a_u32(digest, diagonal_bits);
            digest = fnv1a_u32(digest, (uint32_t)maximum);
            if (maximum < minima[representative]) {
                minima[representative] = maximum;
                minimum_counts[representative] = 1;
            } else if (maximum == minima[representative]) {
                ++minimum_counts[representative];
            }
        }
    }
    if (minima[0] != 44 || minimum_counts[0] != 4 || minima[1] != 48 ||
        minimum_counts[1] != 8) {
        fprintf(stderr, "Hadamard-cloud maximum table mismatch\n");
        return EXIT_FAILURE;
    }
    if (digest != UINT64_C(0x273ea01435d2c8a5)) {
        fprintf(stderr, "Hadamard-cloud maximum-table digest mismatch\n");
        return EXIT_FAILURE;
    }

    printf("representative_minima=trace0:44x4,trace4:48x8\n");
    printf("maxima_table_fnv64=%016llx\n", (unsigned long long)digest);
    verify_fixed_half_attainment();
    printf("corruption_controls=hadamard_entry,base_anti_isomorphism,fixed_half_total,fine_spin_witness\n");
    printf("hadamard_cloud_lift_verification=PASSED\n");
    return EXIT_SUCCESS;
}
