/*
 * Exact Boolean-alignment checks for small prime-order Paley conference
 * matrices.  This program is self-contained: it constructs each matrix over
 * F_q, verifies C C^T = q I, and scans one representative of every antipodal
 * Boolean pair by Gray code.
 *
 * Build:
 *   cc -std=c11 -O3 -Wall -Wextra -Wpedantic \
 *      verification/research_paley_alignment.c -lm \
 *      -o /tmp/research_paley_alignment
 */

#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { MAX_ORDER = 30 };

typedef struct {
    int prime;
    int expected_maximum;
    uint64_t expected_projective_maximizers;
    int ratio_square_numerator;
    int ratio_square_denominator;
    const char *ratio_formula;
} TestCase;

typedef struct {
    int maximum;
    uint64_t projective_maximizers;
} ScanResult;

static void build_paley_conference(int prime,
                                   int8_t matrix[MAX_ORDER][MAX_ORDER]) {
    const int order = prime + 1;
    int is_nonzero_square[MAX_ORDER] = {0};

    memset(matrix, 0, sizeof(int8_t) * MAX_ORDER * MAX_ORDER);
    for (int value = 1; value < prime; ++value) {
        is_nonzero_square[(value * value) % prime] = 1;
    }

    for (int index = 1; index < order; ++index) {
        matrix[0][index] = 1;
        matrix[index][0] = 1;
    }
    for (int row = 0; row < prime; ++row) {
        for (int column = 0; column < prime; ++column) {
            if (row == column) {
                continue;
            }
            int difference = (row - column) % prime;
            if (difference < 0) {
                difference += prime;
            }
            matrix[row + 1][column + 1] =
                is_nonzero_square[difference] ? 1 : -1;
        }
    }
}

static int has_conference_identity(
    const int8_t matrix[MAX_ORDER][MAX_ORDER], int order, int prime) {
    for (int row = 0; row < order; ++row) {
        if (matrix[row][row] != 0) {
            return 0;
        }
        for (int column = 0; column < order; ++column) {
            if (row != column &&
                ((matrix[row][column] != 1 &&
                  matrix[row][column] != -1) ||
                 matrix[column][row] != matrix[row][column])) {
                return 0;
            }

            int product = 0;
            for (int index = 0; index < order; ++index) {
                product += matrix[row][index] * matrix[column][index];
            }
            const int expected = row == column ? prime : 0;
            if (product != expected) {
                return 0;
            }
        }
    }
    return 1;
}

static ScanResult scan_projective_boolean_states(
    const int8_t matrix[MAX_ORDER][MAX_ORDER], int order) {
    int spin[MAX_ORDER];
    int field[MAX_ORDER];
    int energy = 0;

    for (int row = 0; row < order; ++row) {
        spin[row] = 1;
        field[row] = 0;
        for (int column = 0; column < order; ++column) {
            field[row] += matrix[row][column];
        }
        energy += field[row];
    }
    energy /= 2;

    ScanResult result = {abs(energy), 1};
    const uint64_t state_count = UINT64_C(1) << (order - 1);
    for (uint64_t step = 1; step < state_count; ++step) {
        const uint64_t gray = step ^ (step >> 1);
        const uint64_t previous_step = step - 1;
        const uint64_t previous_gray =
            previous_step ^ (previous_step >> 1);
        const uint64_t changed_bit = gray ^ previous_gray;
        const int vertex = (int)__builtin_ctzll(changed_bit) + 1;

        const int old_spin = spin[vertex];
        energy -= 2 * old_spin * field[vertex];
        spin[vertex] = -old_spin;
        for (int row = 0; row < order; ++row) {
            field[row] -= 2 * matrix[row][vertex] * old_spin;
        }

        const int absolute_energy = abs(energy);
        if (absolute_energy > result.maximum) {
            result.maximum = absolute_energy;
            result.projective_maximizers = 1;
        } else if (absolute_energy == result.maximum) {
            ++result.projective_maximizers;
        }
    }
    return result;
}

static int verify_case(const TestCase *test_case) {
    const int order = test_case->prime + 1;
    int8_t matrix[MAX_ORDER][MAX_ORDER];
    build_paley_conference(test_case->prime, matrix);

    if (!has_conference_identity(matrix, order, test_case->prime)) {
        fprintf(stderr, "conference identity failed for q=%d\n",
                test_case->prime);
        return 0;
    }

    const ScanResult result =
        scan_projective_boolean_states(matrix, order);
    if (result.maximum != test_case->expected_maximum ||
        result.projective_maximizers !=
            test_case->expected_projective_maximizers) {
        fprintf(stderr,
                "Boolean scan failed for q=%d: got M=%d, count=%" PRIu64
                "; expected M=%d, count=%" PRIu64 "\n",
                test_case->prime, result.maximum,
                result.projective_maximizers,
                test_case->expected_maximum,
                test_case->expected_projective_maximizers);
        return 0;
    }

    const int64_t left =
        INT64_C(4) * result.maximum * result.maximum *
        test_case->ratio_square_denominator;
    const int64_t right =
        (int64_t)order * order * test_case->prime *
        test_case->ratio_square_numerator;
    if (left != right) {
        fprintf(stderr, "exact ratio-square check failed for q=%d\n",
                test_case->prime);
        return 0;
    }

    const double ratio =
        (2.0 * result.maximum) /
        ((double)order * sqrt((double)test_case->prime));
    printf("q=%d order=%d conference=PASS M=%d "
           "projective_maximizers=%" PRIu64 " "
           "ratio=%s ratio_squared=%d/%d decimal=%.15f\n",
           test_case->prime, order, result.maximum,
           result.projective_maximizers, test_case->ratio_formula,
           test_case->ratio_square_numerator,
           test_case->ratio_square_denominator, ratio);
    return 1;
}

static int verify_corruption_control(void) {
    const int prime = 5;
    const int order = prime + 1;
    int8_t matrix[MAX_ORDER][MAX_ORDER];
    build_paley_conference(prime, matrix);

    matrix[0][1] = (int8_t)-matrix[0][1];
    matrix[1][0] = (int8_t)-matrix[1][0];
    if (has_conference_identity(matrix, order, prime)) {
        fprintf(stderr, "symmetric-edge corruption was not detected\n");
        return 0;
    }
    printf("corruption_control=symmetric_edge_flip_detected\n");
    return 1;
}

int main(void) {
    static const TestCase test_cases[] = {
        {5, 5, UINT64_C(12), 5, 9, "sqrt(5)/3"},
        {13, 21, UINT64_C(156), 9, 13, "3/sqrt(13)"},
        {17, 33, UINT64_C(204), 121, 153, "11/(3*sqrt(17))"},
        {29, 75, UINT64_C(812), 25, 29, "5/sqrt(29)"},
    };

    for (size_t index = 0;
         index < sizeof(test_cases) / sizeof(test_cases[0]); ++index) {
        if (!verify_case(&test_cases[index])) {
            return EXIT_FAILURE;
        }
    }
    if (!verify_corruption_control()) {
        return EXIT_FAILURE;
    }
    printf("all_checks=PASS\n");
    return EXIT_SUCCESS;
}
