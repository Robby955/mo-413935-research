/*
 * Exact threshold scanner for the F(13) certificate.
 *
 * Input is a newline-delimited graph6 stream of graphs on 11 vertices.  A
 * graph is converted to a root-normalized signing of K_12: graph edges have
 * coefficient -1, nonedges have coefficient +1, and a twelfth root vertex
 * has coefficient +1 to every residual vertex.
 *
 * The program enumerates all 2^11 projective spin classes by a Gray-code
 * update.  It stops scanning a signing as soon as it witnesses |Q| >= 19.
 * Since every order-12 energy is even, this rejects exactly the signings with
 * M >= 20.  A record printed as a survivor has received a complete scan and
 * therefore has its exact M <= 18.
 *
 * This program deliberately does not claim a complete graph catalogue.  The
 * Python driver hashes and counts the exact nauty stream fed to it and checks
 * both producer and consumer exit statuses.
 */

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum {
    RESIDUAL_ORDER = 11,
    SIGNING_ORDER = 12,
    GRAPH6_EDGE_COUNT = 55,
    GRAPH6_DATA_CHARS = 10,
    GRAPH6_RECORD_CHARS = 11,
    THRESHOLD = 19,
};

static int decode_graph6(const char *record, int8_t matrix[12][12]) {
    if (strlen(record) != GRAPH6_RECORD_CHARS ||
        (unsigned char)record[0] != (unsigned char)(63 + RESIDUAL_ORDER)) {
        return 0;
    }

    memset(matrix, 0, SIGNING_ORDER * SIGNING_ORDER * sizeof(matrix[0][0]));
    int bit_index = 0;
    for (int data_index = 0; data_index < GRAPH6_DATA_CHARS; ++data_index) {
        int value = (int)(unsigned char)record[1 + data_index] - 63;
        if (value < 0 || value >= 64) {
            return 0;
        }
        for (int shift = 5; shift >= 0; --shift) {
            int bit = (value >> shift) & 1;
            if (bit_index < GRAPH6_EDGE_COUNT) {
                int remaining = bit_index;
                int column = 1;
                while (remaining >= column) {
                    remaining -= column;
                    ++column;
                }
                int row = remaining;
                int coefficient = bit ? -1 : 1;
                matrix[row][column] = (int8_t)coefficient;
                matrix[column][row] = (int8_t)coefficient;
            } else if (bit != 0) {
                return 0;
            }
            ++bit_index;
        }
    }

    for (int vertex = 0; vertex < RESIDUAL_ORDER; ++vertex) {
        matrix[vertex][RESIDUAL_ORDER] = 1;
        matrix[RESIDUAL_ORDER][vertex] = 1;
    }
    return 1;
}

/* Returns 1 for a complete M < THRESHOLD scan and 0 after a safe rejection. */
static int exact_survivor_maximum(const int8_t matrix[12][12], int *maximum) {
    int spins[SIGNING_ORDER];
    int fields[SIGNING_ORDER];
    int energy_twice = 0;
    for (int row = 0; row < SIGNING_ORDER; ++row) {
        spins[row] = 1;
        fields[row] = 0;
        for (int column = 0; column < SIGNING_ORDER; ++column) {
            fields[row] += matrix[row][column];
        }
        energy_twice += fields[row];
    }
    int energy = energy_twice / 2;
    int best = abs(energy);
    if (best >= THRESHOLD) {
        return 0;
    }

    const unsigned spin_classes = 1U << (SIGNING_ORDER - 1);
    for (unsigned step = 1; step < spin_classes; ++step) {
        unsigned gray = step ^ (step >> 1);
        unsigned previous = (step - 1) ^ ((step - 1) >> 1);
        unsigned changed = gray ^ previous;
        int vertex = __builtin_ctz(changed);

        int old_spin = spins[vertex];
        energy -= 2 * old_spin * fields[vertex];
        spins[vertex] = -old_spin;
        for (int row = 0; row < SIGNING_ORDER; ++row) {
            fields[row] -= 2 * matrix[row][vertex] * old_spin;
        }

        int absolute_energy = abs(energy);
        if (absolute_energy >= THRESHOLD) {
            return 0;
        }
        if (absolute_energy > best) {
            best = absolute_energy;
        }
    }

    *maximum = best;
    return 1;
}

int main(void) {
    char line[128];
    uint64_t records = 0;
    uint64_t survivors = 0;

    while (fgets(line, sizeof(line), stdin) != NULL) {
        size_t length = strlen(line);
        while (length > 0 && (line[length - 1] == '\n' || line[length - 1] == '\r')) {
            line[--length] = '\0';
        }
        if (length == 0) {
            continue;
        }
        if (length + 1 == sizeof(line) && line[length - 1] != '\n') {
            fprintf(stderr, "overlong graph6 record at input record %" PRIu64 "\n", records + 1);
            return 2;
        }

        int8_t matrix[SIGNING_ORDER][SIGNING_ORDER];
        if (!decode_graph6(line, matrix)) {
            fprintf(stderr, "invalid graph6 record at input record %" PRIu64 "\n", records + 1);
            return 3;
        }

        ++records;
        int maximum = 0;
        if (exact_survivor_maximum(matrix, &maximum)) {
            ++survivors;
            printf("SURVIVOR %s M=%d\n", line, maximum);
        }
    }

    if (ferror(stdin)) {
        fprintf(stderr, "input read failure\n");
        return 4;
    }
    printf("SUMMARY records=%" PRIu64 " survivors=%" PRIu64 "\n", records, survivors);
    return 0;
}
