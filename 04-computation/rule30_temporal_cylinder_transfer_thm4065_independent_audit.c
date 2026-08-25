#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* No-import audit for THM-4065.

   Unlike the primary reset-anchored solver, this implementation composes the
   affine map around the cylinder from residue zero, proves that the images of
   incoming states zero and one coincide, and then makes a second pass from
   the forced cyclic initial bit.  The finite graph audit materializes the
   complete successor/predecessor tables before decomposing components. */

typedef struct {
    uint32_t a;
    uint32_t b;
} State;

static const uint64_t EXPECTED_STEPS = UINT64_C(1420778968);

static void fail(const char *label) {
    fprintf(stderr, "requirement failed: %s\n", label);
    exit(2);
}

static void require_true(int condition, const char *label) {
    if (!condition) {
        fail(label);
    }
}

static uint32_t period_mask(unsigned period) {
    require_true(period >= 1 && period <= 32, "period range");
    return period == 32 ? UINT32_MAX : (UINT32_C(1) << period) - UINT32_C(1);
}

static unsigned bit_at(uint32_t word, unsigned time) {
    return (unsigned)((word >> time) & UINT32_C(1));
}

static unsigned step_bit(unsigned a, unsigned b, unsigned x) {
    return a ^ b ^ x ^ (b & x);
}

static uint32_t solve_by_affine_monodromy(
    uint32_t a, uint32_t b, unsigned period
) {
    const uint32_t mask = period_mask(period);
    a &= mask;
    b &= mask;
    require_true(b != 0, "nonzero driver");

    unsigned image_zero = 0;
    unsigned image_one = 1;
    for (unsigned time = 0; time < period; ++time) {
        const unsigned av = bit_at(a, time);
        const unsigned bv = bit_at(b, time);
        image_zero = step_bit(av, bv, image_zero);
        image_one = step_bit(av, bv, image_one);
    }
    require_true(image_zero == image_one, "reset monodromy is constant");

    const unsigned initial = image_zero;
    unsigned current = initial;
    uint32_t x = 0;
    for (unsigned time = 0; time < period; ++time) {
        x |= (uint32_t)current << time;
        current = step_bit(bit_at(a, time), bit_at(b, time), current);
    }
    require_true(current == initial, "cyclic wrap");
    return x & mask;
}

static uint32_t inverse_formula(uint32_t b, uint32_t x, unsigned period) {
    uint32_t a = 0;
    for (unsigned time = 0; time < period; ++time) {
        const unsigned next = time + 1 == period ? 0 : time + 1;
        const unsigned bv = bit_at(b, time);
        const unsigned xv = bit_at(x, time);
        const unsigned av = bit_at(x, next) ^ bv ^ xv ^ (bv & xv);
        a |= (uint32_t)av << time;
    }
    return a & period_mask(period);
}

static size_t encode(State state, uint32_t word_count) {
    return (size_t)state.a * (size_t)word_count + (size_t)state.b;
}

static State decode(size_t index, uint32_t word_count) {
    State state;
    state.a = (uint32_t)(index / word_count);
    state.b = (uint32_t)(index % word_count);
    return state;
}

static void audit_materialized_graph(unsigned period) {
    require_true(period <= 8, "finite graph audit range");
    const uint32_t word_count = period_mask(period) + UINT32_C(1);
    const size_t vertex_count = (size_t)word_count * (size_t)word_count;
    const size_t absent = SIZE_MAX;
    size_t *successor = (size_t *)malloc(vertex_count * sizeof(size_t));
    size_t *predecessor = (size_t *)malloc(vertex_count * sizeof(size_t));
    unsigned char *seen = (unsigned char *)calloc(vertex_count, sizeof(unsigned char));
    unsigned char *endpoints = (unsigned char *)calloc(word_count, sizeof(unsigned char));
    require_true(
        successor != NULL && predecessor != NULL && seen != NULL && endpoints != NULL,
        "graph allocation"
    );
    for (size_t index = 0; index < vertex_count; ++index) {
        successor[index] = absent;
        predecessor[index] = absent;
    }

    uint64_t inverse_checks = 0;
    for (uint32_t a = 0; a < word_count; ++a) {
        for (uint32_t b = 1; b < word_count; ++b) {
            const uint32_t x = solve_by_affine_monodromy(a, b, period);
            require_true(inverse_formula(b, x, period) == a, "symbolic inverse");
            const State source = {a, b};
            const State target = {b, x};
            const size_t source_index = encode(source, word_count);
            const size_t target_index = encode(target, word_count);
            require_true(successor[source_index] == absent, "one successor");
            require_true(predecessor[target_index] == absent, "one predecessor");
            successor[source_index] = target_index;
            predecessor[target_index] = source_index;
            ++inverse_checks;
        }
    }
    for (uint32_t first = 0; first < word_count; ++first) {
        for (uint32_t second = 0; second < word_count; ++second) {
            const size_t index = encode((State){first, second}, word_count);
            require_true(
                (successor[index] != absent) == (second != 0),
                "outdegree boundary"
            );
            require_true(
                (predecessor[index] != absent) == (first != 0),
                "indegree boundary"
            );
        }
    }

    uint64_t path_vertices = 0;
    uint64_t path_edges = 0;
    uint64_t max_path_edges = 0;
    for (uint32_t start_b = word_count - 1; start_b >= 1; --start_b) {
        size_t index = encode((State){0, start_b}, word_count);
        uint64_t edges = 0;
        while (successor[index] != absent) {
            require_true(!seen[index], "disjoint source paths");
            seen[index] = 1;
            ++path_vertices;
            index = successor[index];
            ++edges;
            require_true(edges <= vertex_count, "path termination");
        }
        require_true(!seen[index], "distinct sinks");
        seen[index] = 1;
        ++path_vertices;
        const State sink = decode(index, word_count);
        require_true(sink.a != 0 && sink.b == 0, "sink type");
        require_true(!endpoints[sink.a], "endpoint injective");
        endpoints[sink.a] = 1;
        path_edges += edges;
        if (edges > max_path_edges) {
            max_path_edges = edges;
        }
        if (start_b == 1) {
            break;
        }
    }
    for (uint32_t terminal = 1; terminal < word_count; ++terminal) {
        require_true(endpoints[terminal], "endpoint surjective");
    }

    uint64_t cycle_count = 0;
    uint64_t cycle_vertices = 0;
    uint64_t max_cycle = 0;
    for (size_t start = 0; start < vertex_count; ++start) {
        const State state = decode(start, word_count);
        if (state.a == 0 || state.b == 0 || seen[start]) {
            continue;
        }
        size_t index = start;
        uint64_t length = 0;
        do {
            require_true(!seen[index], "cycle disjointness");
            seen[index] = 1;
            ++cycle_vertices;
            ++length;
            index = successor[index];
            require_true(index != absent, "cycle successor");
        } while (index != start);
        ++cycle_count;
        if (length > max_cycle) {
            max_cycle = length;
        }
    }
    require_true(path_vertices + cycle_vertices + 1 == vertex_count, "component census");
    require_true(
        inverse_checks == (uint64_t)word_count * (uint64_t)(word_count - 1),
        "inverse check count"
    );

    printf(
        "graph_P=%u;words=%" PRIu32 ";inverse_checks=%" PRIu64
        ";boundary_paths=%" PRIu32 ";path_vertices=%" PRIu64
        ";path_edges=%" PRIu64 ";max_path_edges=%" PRIu64
        ";max_zero_gap=%" PRIu64 ";cycles=%" PRIu64
        ";cycle_vertices=%" PRIu64 ";max_cycle=%" PRIu64 "\n",
        period, word_count, inverse_checks, word_count - 1, path_vertices,
        path_edges, max_path_edges, max_path_edges + 1, cycle_count,
        cycle_vertices, max_cycle
    );
    free(endpoints);
    free(seen);
    free(predecessor);
    free(successor);
}

static unsigned least_period(uint32_t word) {
    const unsigned divisors[] = {1, 2, 4, 8, 16, 32};
    for (size_t item = 0; item < sizeof(divisors) / sizeof(divisors[0]); ++item) {
        const unsigned candidate = divisors[item];
        int agrees = 1;
        for (unsigned time = candidate; time < 32; ++time) {
            if (bit_at(word, time) != bit_at(word, time % candidate)) {
                agrees = 0;
                break;
            }
        }
        if (agrees) {
            return candidate;
        }
    }
    fail("least period");
    return 0;
}

static uint32_t integrate_difference(uint32_t difference, unsigned initial) {
    uint32_t word = 0;
    unsigned value = initial;
    for (unsigned time = 0; time < 32; ++time) {
        if (value) {
            word ^= UINT32_C(1) << time;
        }
        value ^= bit_at(difference, time);
    }
    require_true(value == initial, "difference integrates cyclically");
    return word;
}

static uint64_t rotate_left_64(uint64_t value, unsigned amount) {
    const unsigned shift = amount & 63U;
    if (shift == 0) {
        return value;
    }
    return (value << shift) | (value >> (64U - shift));
}

static void audit_boundary_prefix(void) {
    State state = {UINT32_C(0x00000000), UINT32_C(0xcf3030cf)};
    for (uint64_t step = 1; step <= UINT64_C(12133); ++step) {
        require_true(state.b != 0, "no prefix zero");
        const uint32_t next = solve_by_affine_monodromy(state.a, state.b, 32);
        state.a = state.b;
        state.b = next;
    }
    require_true(
        state.a == UINT32_C(0xbe79924b)
        && state.b == UINT32_C(0x90f58380),
        "prefix endpoint"
    );
    printf("boundary_prefix_check=start_depths:87866,87867;words:00000000,cf3030cf;"
           "edges:12133;end_depths:99999,100000;words:be79924b,90f58380\n");
}

static void audit_long_orbit(void) {
    State state = {UINT32_C(0xbe79924b), UINT32_C(0x90f58380)};
    uint64_t xor_digest = 0;
    uint64_t sum_digest = 0;
    uint64_t fnv_digest = UINT64_C(1469598103934665603);
    const uint64_t checkpoint_step[] = {
        UINT64_C(1), UINT64_C(10), UINT64_C(1000), UINT64_C(1000000),
        UINT64_C(100000000), UINT64_C(1000000000), UINT64_C(1420778968)
    };
    const State checkpoint_state[] = {
        {UINT32_C(0x90f58380), UINT32_C(0xcb08e372)},
        {UINT32_C(0xbc973f1e), UINT32_C(0x25ec4c68)},
        {UINT32_C(0x73703da9), UINT32_C(0x07f4bd13)},
        {UINT32_C(0xe7f1f52a), UINT32_C(0xd37751d3)},
        {UINT32_C(0xd4f2a1bc), UINT32_C(0x7e8f9310)},
        {UINT32_C(0x36363150), UINT32_C(0xc4144da0)},
        {UINT32_C(0xdae372fa), UINT32_C(0x00000000)}
    };
    size_t checkpoint = 0;

    for (uint64_t step = 1; step <= EXPECTED_STEPS; ++step) {
        require_true(state.b != 0, "no premature zero");
        const uint32_t next = solve_by_affine_monodromy(state.a, state.b, 32);
        state.a = state.b;
        state.b = next;

        const uint64_t packed = ((uint64_t)state.a << 32) | state.b;
        xor_digest ^= rotate_left_64(packed, (unsigned)step)
                      ^ step * UINT64_C(0x9e3779b97f4a7c15);
        sum_digest += packed + step * UINT64_C(0xd6e8feb86659fd93);
        fnv_digest ^= packed ^ step;
        fnv_digest *= UINT64_C(1099511628211);

        if (checkpoint < sizeof(checkpoint_step) / sizeof(checkpoint_step[0])
            && step == checkpoint_step[checkpoint]) {
            require_true(
                state.a == checkpoint_state[checkpoint].a
                && state.b == checkpoint_state[checkpoint].b,
                "independent checkpoint"
            );
            printf(
                "orbit_checkpoint=step:%" PRIu64 ";a:%08" PRIx32
                ";b:%08" PRIx32 "\n",
                step, state.a, state.b
            );
            ++checkpoint;
        }
        if (state.b == 0) {
            require_true(step == EXPECTED_STEPS, "terminal step");
            break;
        }
    }
    require_true(state.a == UINT32_C(0xdae372fa) && state.b == 0,
                 "terminal state");
    require_true(xor_digest == UINT64_C(0x25469bdaa1c33b27), "xor digest");
    require_true(sum_digest == UINT64_C(0x65b743c483016d28), "sum digest");
    require_true(fnv_digest == UINT64_C(0x46485ada15e2f5b9), "fnv digest");

    const unsigned weight = (unsigned)__builtin_popcount(state.a);
    const uint32_t candidate_zero = integrate_difference(state.a, 0);
    const uint32_t candidate_one = integrate_difference(state.a, 1);
    require_true(weight == 20, "terminal weight");
    require_true(least_period(state.a) == 32, "terminal period");
    require_true(candidate_zero == UINT32_C(0x93425cac), "first candidate");
    require_true(candidate_one == UINT32_C(0x6cbda353), "second candidate");
    require_true(least_period(candidate_zero) == 32, "first candidate period");
    require_true(least_period(candidate_one) == 32, "second candidate period");

    printf("parent_state=depths:99999,100000;words:be79924b,90f58380\n");
    printf("transfer_steps_to_zero=%" PRIu64 "\n", EXPECTED_STEPS);
    printf("next_zero_depth=1420878968\n");
    printf("previous_zero_depth=87866;exact_zero_gap=1420791102\n");
    printf("boundary_path_edges=1420791101;prefix_edges_to_parent_state=12133\n");
    printf("terminal_predecessor=%08" PRIx32 ";weight=%u;least_period=32\n",
           state.a, weight);
    printf("next_successor_candidates=%08" PRIx32 ",%08" PRIx32
           ";physical_phase=unresolved;least_period=32\n",
           candidate_zero, candidate_one);
    printf("orbit_xor_digest=%016" PRIx64 "\n", xor_digest);
    printf("orbit_sum_digest=%016" PRIx64 "\n", sum_digest);
    printf("orbit_fnv_digest=%016" PRIx64 "\n", fnv_digest);
    printf("period32_crude_next_zero_gap_bound=18446744065119617027\n");
}

int main(void) {
    printf("RULE30_TEMPORAL_CYLINDER_TRANSFER_THM4065_INDEPENDENT_AUDIT\n");
    printf("solver=two_pass_affine_monodromy;graph=materialized_successor_predecessor\n");
    printf("universal_map=T_P(A,B)=(B,X);domain:B_nonzero;codomain:first_nonzero\n");
    audit_materialized_graph(1);
    audit_materialized_graph(2);
    audit_materialized_graph(4);
    audit_materialized_graph(8);
    audit_boundary_prefix();
    audit_long_orbit();
    printf("parent=THM4047_fixed_query_packet;no_primary_source_import\n");
    printf("scope=PROVED universal transfer and infinite isolated zeros;FINITE-EXACT next zero;moving center and Rule30 prizes OPEN\n");
    return 0;
}
