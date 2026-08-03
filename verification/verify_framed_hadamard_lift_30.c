#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

enum { CLOUDS = 4, CLOUD_SIZE = 4, ORDER = 16, STATES = 32768 };

static const int h[4][4] = {
    {1, 1, 1, 1},
    {1, 1, -1, -1},
    {1, -1, 1, -1},
    {-1, 1, 1, -1},
};
static const int pairs[6][2] = {
    {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3},
};
static const int base[6] = {1, 1, 1, -1, 1, -1};
static const int p_edges[6] = {1, 1, 1, -1, 1, -1};
static const int r_edges[6] = {-1, -1, 1, 1, 1, -1};
static const int *const internal[4] = {
    p_edges, r_edges, p_edges, r_edges,
};

static const int obstruction[6][4][4] = {
    {{1, 1, 1, 1}, {1, 1, 1, 1}, {-1, -1, -1, -1}, {1, 1, 1, -1}},
    {{1, 1, 1, 1}, {-1, -1, 1, 1}, {-1, -1, -1, -1}, {-1, -1, -1, 1}},
    {{1, 1, -1, -1}, {1, 1, -1, -1}, {-1, -1, 1, 1}, {1, 1, -1, 1}},
    {{1, 1, -1, -1}, {-1, -1, 1, -1}, {-1, -1, 1, 1}, {-1, -1, 1, -1}},
    {{1, 1, -1, -1}, {1, 1, -1, 1}, {-1, -1, 1, 1}, {1, 1, -1, 1}},
    {{1, 1, -1, -1}, {-1, -1, -1, -1}, {-1, -1, 1, 1}, {-1, -1, 1, -1}},
};

static int edge_energy(const int edges[6], const int spins[4]) {
    int result = 0;
    for (int edge = 0; edge < 6; ++edge) {
        const int left = pairs[edge][0];
        const int right = pairs[edge][1];
        result += edges[edge] * spins[left] * spins[right];
    }
    return result;
}

static int cross_energy(const int clouds[4][4]) {
    int result = 0;
    for (int edge = 0; edge < 6; ++edge) {
        const int left = pairs[edge][0];
        const int right = pairs[edge][1];
        int block = 0;
        for (int row = 0; row < 4; ++row) {
            for (int column = 0; column < 4; ++column) {
                block += clouds[left][row] * h[row][column]
                         * clouds[right][column];
            }
        }
        result += base[edge] * block;
    }
    return result;
}

static int block_energy(const int spins[16]) {
    int clouds[4][4];
    int result;
    for (int cloud = 0; cloud < 4; ++cloud) {
        for (int coordinate = 0; coordinate < 4; ++coordinate) {
            clouds[cloud][coordinate] = spins[4 * cloud + coordinate];
        }
    }
    result = cross_energy(clouds);
    for (int cloud = 0; cloud < 4; ++cloud) {
        result += edge_energy(internal[cloud], clouds[cloud]);
    }
    return result;
}

static void assemble(int matrix[16][16]) {
    for (int row = 0; row < 16; ++row) {
        for (int column = 0; column < 16; ++column) {
            matrix[row][column] = 0;
        }
    }
    for (int cloud = 0; cloud < 4; ++cloud) {
        for (int edge = 0; edge < 6; ++edge) {
            const int left = 4 * cloud + pairs[edge][0];
            const int right = 4 * cloud + pairs[edge][1];
            matrix[left][right] = internal[cloud][edge];
            matrix[right][left] = internal[cloud][edge];
        }
    }
    for (int edge = 0; edge < 6; ++edge) {
        const int left_cloud = pairs[edge][0];
        const int right_cloud = pairs[edge][1];
        for (int row = 0; row < 4; ++row) {
            for (int column = 0; column < 4; ++column) {
                const int left = 4 * left_cloud + row;
                const int right = 4 * right_cloud + column;
                const int value = base[edge] * h[row][column];
                matrix[left][right] = value;
                matrix[right][left] = value;
            }
        }
    }
}

static int direct_energy(const int matrix[16][16], const int spins[16]) {
    int result = 0;
    for (int left = 0; left < 16; ++left) {
        for (int right = left + 1; right < 16; ++right) {
            result += matrix[left][right] * spins[left] * spins[right];
        }
    }
    return result;
}

static void decode_projective(unsigned int mask, int spins[16]) {
    spins[0] = 1;
    for (int coordinate = 1; coordinate < 16; ++coordinate) {
        const unsigned int bit = 1U << (unsigned int)(coordinate - 1);
        spins[coordinate] = (mask & bit) != 0U ? -1 : 1;
    }
}

static int projective_index(const int spins[4]) {
    const int sign = spins[0];
    int index = 0;
    for (int coordinate = 1; coordinate < 4; ++coordinate) {
        index *= 2;
        if (sign * spins[coordinate] == -1) {
            index += 1;
        }
    }
    return index;
}

static int projectively_equal(const int left[4], const int right[4]) {
    const int sign = left[0] * right[0];
    for (int coordinate = 0; coordinate < 4; ++coordinate) {
        if (left[coordinate] != sign * right[coordinate]) {
            return 0;
        }
    }
    return 1;
}

static int verify_hadamard(void) {
    for (int left = 0; left < 4; ++left) {
        for (int right = 0; right < 4; ++right) {
            int product = 0;
            for (int coordinate = 0; coordinate < 4; ++coordinate) {
                product += h[left][coordinate] * h[right][coordinate];
            }
            if (product != (left == right ? 4 : 0)) {
                return 0;
            }
        }
    }
    return 1;
}

static void representative(int index, int spins[4]) {
    spins[0] = 1;
    for (int coordinate = 3; coordinate >= 1; --coordinate) {
        spins[coordinate] = (index & 1) != 0 ? -1 : 1;
        index >>= 1;
    }
}

static void edges_from_mask(unsigned int mask, int edges[6]) {
    for (int edge = 0; edge < 6; ++edge) {
        const unsigned int bit = 1U << (unsigned int)edge;
        edges[edge] = (mask & bit) != 0U ? -1 : 1;
    }
}

static void verify_obstruction(void) {
    const int expected_sources[3] = {0, 3, 2};
    const int expected_targets[3] = {3, 2, 0};
    for (int pair = 0; pair < 3; ++pair) {
        const int *positive_flat = &obstruction[2 * pair][0][0];
        const int *negative_flat = &obstruction[2 * pair + 1][0][0];
        const int (*positive)[4] = (const int (*)[4])positive_flat;
        const int (*negative)[4] = (const int (*)[4])negative_flat;
        if (cross_energy(positive) != 28 || cross_energy(negative) != -28) {
            fprintf(stderr, "obstruction cross energy mismatch\n");
            exit(EXIT_FAILURE);
        }
        for (int cloud = 0; cloud < 4; ++cloud) {
            if (cloud != 1
                && !projectively_equal(positive[cloud], negative[cloud])) {
                fprintf(stderr, "obstruction cancellation mismatch\n");
                exit(EXIT_FAILURE);
            }
        }
        if (projective_index(positive[1]) != expected_sources[pair]
            || projective_index(negative[1]) != expected_targets[pair]) {
            fprintf(stderr, "obstruction cycle mismatch\n");
            exit(EXIT_FAILURE);
        }
    }

    for (unsigned int mask = 0U; mask < 64U; ++mask) {
        int edges[6];
        int u0[4];
        int u2[4];
        int u3[4];
        int q0;
        int q2;
        int q3;
        edges_from_mask(mask, edges);
        representative(0, u0);
        representative(2, u2);
        representative(3, u3);
        q0 = edge_energy(edges, u0);
        q2 = edge_energy(edges, u2);
        q3 = edge_energy(edges, u3);
        if (q0 - q2 != 2 * (edges[1] + edges[3] + edges[5])) {
            fprintf(stderr, "odd-sum identity mismatch\n");
            exit(EXIT_FAILURE);
        }
        if (q0 == q3 && q3 == q2) {
            fprintf(stderr, "impossible response equality occurred\n");
            exit(EXIT_FAILURE);
        }
    }
}

static void verify_order_four(void) {
    int optimum = 100;
    int p_maximum = 0;
    int r_maximum = 0;
    for (unsigned int mask = 0U; mask < 64U; ++mask) {
        int edges[6];
        int maximum = 0;
        edges_from_mask(mask, edges);
        for (int index = 0; index < 8; ++index) {
            int spins[4];
            int energy;
            representative(index, spins);
            energy = abs(edge_energy(edges, spins));
            if (energy > maximum) {
                maximum = energy;
            }
        }
        if (maximum < optimum) {
            optimum = maximum;
        }
    }
    for (int index = 0; index < 8; ++index) {
        int spins[4];
        int p_energy;
        int r_energy;
        representative(index, spins);
        p_energy = abs(edge_energy(p_edges, spins));
        r_energy = abs(edge_energy(r_edges, spins));
        if (p_energy > p_maximum) {
            p_maximum = p_energy;
        }
        if (r_energy > r_maximum) {
            r_maximum = r_energy;
        }
    }
    if (optimum != 4 || p_maximum != 4 || r_maximum != 4) {
        fprintf(stderr, "order-four maximum mismatch\n");
        exit(EXIT_FAILURE);
    }
}

static void verify_common_internal_profile(void) {
    static int cross_values[STATES];
    static uint8_t response_indices[STATES][4];
    int profile[47] = {0};
    const int expected[47] = {
        [38] = 6, [40] = 18, [42] = 32, [46] = 8,
    };
    for (unsigned int mask = 0U; mask < (1U << 15U); ++mask) {
        int spins[16];
        int clouds[4][4];
        decode_projective(mask, spins);
        for (int cloud = 0; cloud < 4; ++cloud) {
            for (int coordinate = 0; coordinate < 4; ++coordinate) {
                clouds[cloud][coordinate] = spins[4 * cloud + coordinate];
            }
            response_indices[mask][cloud]
                = (uint8_t)projective_index(clouds[cloud]);
        }
        cross_values[mask] = cross_energy(clouds);
    }
    for (unsigned int edge_mask = 0U; edge_mask < 64U; ++edge_mask) {
        int edges[6];
        int responses[8];
        int maximum = 0;
        edges_from_mask(edge_mask, edges);
        for (int index = 0; index < 8; ++index) {
            int spins[4];
            representative(index, spins);
            responses[index] = edge_energy(edges, spins);
        }
        for (unsigned int state = 0U; state < (1U << 15U); ++state) {
            int energy = cross_values[state];
            for (int cloud = 0; cloud < 4; ++cloud) {
                energy += responses[response_indices[state][cloud]];
            }
            energy = abs(energy);
            if (energy > maximum) {
                maximum = energy;
            }
        }
        if (maximum >= 47) {
            fprintf(stderr, "common-internal maximum outside profile\n");
            exit(EXIT_FAILURE);
        }
        profile[maximum] += 1;
    }
    for (int maximum = 0; maximum < 47; ++maximum) {
        if (profile[maximum] != expected[maximum]) {
            fprintf(stderr, "common-internal profile mismatch\n");
            exit(EXIT_FAILURE);
        }
    }
}

int main(void) {
    int matrix[16][16];
    int histogram[241] = {0};
    int maximum = 0;
    int maximizers = 0;
    int cross_maximum = 0;
    const int expected_counts[31] = {
        38, 88, 186, 456, 670, 1056, 1634, 1904, 2518, 2960, 3050,
        3640, 3774, 4088, 4514, 4384, 4514, 4088, 3774, 3640, 3050,
        2960, 2518, 1904, 1634, 1056, 670, 456, 186, 88, 38,
    };

    if (!verify_hadamard()) {
        fprintf(stderr, "Hadamard check failed\n");
        return EXIT_FAILURE;
    }
    assemble(matrix);
    for (int row = 0; row < 16; ++row) {
        for (int column = 0; column < 16; ++column) {
            if (matrix[row][column] != matrix[column][row]
                || (row == column && matrix[row][column] != 0)
                || (row != column && abs(matrix[row][column]) != 1)) {
                fprintf(stderr, "admissibility check failed\n");
                return EXIT_FAILURE;
            }
        }
    }

    for (unsigned int mask = 0U; mask < (1U << 15U); ++mask) {
        int spins[16];
        int clouds[4][4];
        int direct;
        int blocked;
        int cross;
        decode_projective(mask, spins);
        direct = direct_energy(matrix, spins);
        blocked = block_energy(spins);
        if (direct != blocked) {
            fprintf(stderr, "direct/block energy mismatch\n");
            return EXIT_FAILURE;
        }
        histogram[direct + 120] += 1;
        if (abs(direct) > maximum) {
            maximum = abs(direct);
        }
        for (int cloud = 0; cloud < 4; ++cloud) {
            for (int coordinate = 0; coordinate < 4; ++coordinate) {
                clouds[cloud][coordinate] = spins[4 * cloud + coordinate];
            }
        }
        cross = abs(cross_energy(clouds));
        if (cross > cross_maximum) {
            cross_maximum = cross;
        }
    }
    for (int energy = -30; energy <= 30; energy += 2) {
        const int expected = expected_counts[(energy + 30) / 2];
        if (2 * histogram[energy + 120] != expected) {
            fprintf(stderr, "full histogram mismatch at energy %d\n", energy);
            return EXIT_FAILURE;
        }
        if (abs(energy) == maximum) {
            maximizers += histogram[energy + 120];
        }
    }
    if (maximum != 30 || maximizers != 38 || cross_maximum != 28) {
        fprintf(stderr, "maximum summary mismatch\n");
        return EXIT_FAILURE;
    }

    verify_order_four();
    verify_obstruction();
    verify_common_internal_profile();

    puts("projective_spins_checked=32768");
    puts("direct_block_histogram_match=TRUE");
    puts("lift_maximum=30");
    puts("projective_maximizers=38");
    puts("cross_only_maximum=28");
    puts("six_state_obstruction=PASSED");
    puts("common_literal_internal_minimum=38");
    puts("strict_c_framed_hadamard_lift_30_verification=PASSED");
    return EXIT_SUCCESS;
}
