#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* Exact finite-word certificate for THM-4065.

   Bit t of a word is the value at absolute time residue t.  For cyclic
   F_2-words A,B, with B nonzero, solve the wrapped Rule-30 left-front law

       X_(t+1) = A_t + B_t + (1+B_t) X_t.

   The primary solver starts immediately after the first reset bit of B.
   It also exhausts the complete transfer graphs through period eight and
   walks the exact period-32 physical state inherited from THM-4047 until
   the next zero word.  All unsigned overflow below is intentional modulo
   2^64 and is used only for reproducible orbit digests. */

typedef struct {
    uint32_t a;
    uint32_t b;
} State;

static const uint64_t EXPECTED_STEPS = UINT64_C(1420778968);
static const uint32_t EXPECTED_PREDECESSOR = UINT32_C(0xdae372fa);

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

static unsigned word_bit(uint32_t word, unsigned index) {
    return (unsigned)((word >> index) & UINT32_C(1));
}

static unsigned rotate_index(unsigned index, unsigned period) {
    return index + 1 == period ? 0 : index + 1;
}

static uint32_t solve_after_reset(uint32_t a, uint32_t b, unsigned period) {
    const uint32_t mask = period_mask(period);
    a &= mask;
    b &= mask;
    require_true(b != 0, "nonzero reset word");

    const unsigned reset = (unsigned)__builtin_ctz(b);
    unsigned current = word_bit(a, reset) ^ 1U;
    unsigned index = rotate_index(reset, period);
    uint32_t x = (uint32_t)current << index;

    for (unsigned offset = 1; offset < period; ++offset) {
        unsigned time = reset + offset;
        if (time >= period) {
            time -= period;
        }
        const unsigned av = word_bit(a, time);
        const unsigned bv = word_bit(b, time);
        current = bv ? av ^ 1U : av ^ current;
        index = rotate_index(time, period);
        x |= (uint32_t)current << index;
    }
    return x & mask;
}

static uint32_t recover_first_word(uint32_t b, uint32_t x, unsigned period) {
    const uint32_t mask = period_mask(period);
    uint32_t a = 0;
    b &= mask;
    x &= mask;
    for (unsigned time = 0; time < period; ++time) {
        const unsigned next = rotate_index(time, period);
        const unsigned bv = word_bit(b, time);
        const unsigned xv = word_bit(x, time);
        const unsigned av = word_bit(x, next) ^ bv ^ xv ^ (bv & xv);
        a |= (uint32_t)av << time;
    }
    return a & mask;
}

static void verify_wrapped_law(
    uint32_t a, uint32_t b, uint32_t x, unsigned period
) {
    for (unsigned time = 0; time < period; ++time) {
        const unsigned next = rotate_index(time, period);
        const unsigned av = word_bit(a, time);
        const unsigned bv = word_bit(b, time);
        const unsigned xv = word_bit(x, time);
        const unsigned rhs = av ^ bv ^ xv ^ (bv & xv);
        require_true(word_bit(x, next) == rhs, "wrapped recurrence");
    }
}

static size_t state_index(State state, uint32_t word_count) {
    return (size_t)state.a * (size_t)word_count + (size_t)state.b;
}

static State transfer(State state, unsigned period) {
    require_true(state.b != 0, "transfer domain");
    State answer;
    answer.a = state.b;
    answer.b = solve_after_reset(state.a, state.b, period);
    return answer;
}

static void audit_complete_graph(unsigned period) {
    const uint32_t mask = period_mask(period);
    const uint32_t word_count = period == 32 ? 0 : mask + UINT32_C(1);
    require_true(period <= 8, "finite graph audit range");
    const size_t vertex_count = (size_t)word_count * (size_t)word_count;
    unsigned char *seen = (unsigned char *)calloc(vertex_count, sizeof(unsigned char));
    unsigned char *endpoint = (unsigned char *)calloc(word_count, sizeof(unsigned char));
    require_true(seen != NULL && endpoint != NULL, "finite graph allocation");

    uint64_t inverse_checks = 0;
    for (uint32_t b = 1; b < word_count; ++b) {
        for (uint32_t a = 0; a < word_count; ++a) {
            const uint32_t x = solve_after_reset(a, b, period);
            verify_wrapped_law(a, b, x, period);
            require_true(recover_first_word(b, x, period) == a, "explicit inverse");
            ++inverse_checks;
        }
    }
    require_true(
        inverse_checks == (uint64_t)word_count * (uint64_t)(word_count - 1),
        "domain cardinality"
    );

    uint64_t path_vertices = 0;
    uint64_t path_edges = 0;
    uint64_t max_path_edges = 0;
    for (uint32_t start_b = 1; start_b < word_count; ++start_b) {
        State state = {0, start_b};
        uint64_t edges = 0;
        while (state.b != 0) {
            const size_t index = state_index(state, word_count);
            require_true(!seen[index], "boundary paths disjoint");
            seen[index] = 1;
            ++path_vertices;
            state = transfer(state, period);
            ++edges;
            require_true(edges <= vertex_count, "boundary path terminates");
        }
        const size_t sink_index = state_index(state, word_count);
        require_true(state.a != 0, "nonzero boundary sink");
        require_true(!seen[sink_index], "distinct boundary sinks");
        require_true(!endpoint[state.a], "endpoint permutation injective");
        seen[sink_index] = 1;
        endpoint[state.a] = 1;
        ++path_vertices;
        path_edges += edges;
        if (edges > max_path_edges) {
            max_path_edges = edges;
        }
    }
    for (uint32_t terminal_a = 1; terminal_a < word_count; ++terminal_a) {
        require_true(endpoint[terminal_a], "endpoint permutation surjective");
    }

    uint64_t cycle_count = 0;
    uint64_t cycle_vertices = 0;
    uint64_t max_cycle = 0;
    for (uint32_t a = 1; a < word_count; ++a) {
        for (uint32_t b = 1; b < word_count; ++b) {
            State start = {a, b};
            const size_t start_index = state_index(start, word_count);
            if (seen[start_index]) {
                continue;
            }
            State state = start;
            uint64_t length = 0;
            do {
                const size_t index = state_index(state, word_count);
                require_true(!seen[index], "cycle overlap");
                seen[index] = 1;
                ++cycle_vertices;
                ++length;
                state = transfer(state, period);
                require_true(state.a != 0 && state.b != 0, "interior cycle closure");
            } while (state.a != start.a || state.b != start.b);
            ++cycle_count;
            if (length > max_cycle) {
                max_cycle = length;
            }
        }
    }
    require_true(path_vertices + cycle_vertices + 1 == vertex_count, "graph partition");

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
    free(endpoint);
    free(seen);
}

static unsigned least_period(uint32_t word, unsigned ambient_period) {
    for (unsigned candidate = 1; candidate <= ambient_period; ++candidate) {
        if (ambient_period % candidate != 0) {
            continue;
        }
        int matches = 1;
        for (unsigned time = 0; time < ambient_period; ++time) {
            if (word_bit(word, time) != word_bit(word, time % candidate)) {
                matches = 0;
                break;
            }
        }
        if (matches) {
            return candidate;
        }
    }
    fail("least period exists");
    return 0;
}

static uint32_t zero_successor_candidate(uint32_t difference, unsigned initial) {
    unsigned current = initial;
    uint32_t word = 0;
    for (unsigned time = 0; time < 32; ++time) {
        word |= (uint32_t)current << time;
        current ^= word_bit(difference, time);
    }
    require_true(current == initial, "even predecessor weight");
    return word;
}

static uint64_t rotate_left_64(uint64_t value, unsigned amount) {
    amount &= 63U;
    return amount == 0 ? value : (value << amount) | (value >> (64U - amount));
}

static void audit_boundary_prefix(void) {
    State state = {UINT32_C(0x00000000), UINT32_C(0xcf3030cf)};
    for (uint64_t step = 1; step <= UINT64_C(12133); ++step) {
        require_true(state.b != 0, "no zero between boundary and parent state");
        state = transfer(state, 32);
    }
    require_true(
        state.a == UINT32_C(0xbe79924b)
        && state.b == UINT32_C(0x90f58380),
        "boundary prefix reaches parent state"
    );
    printf("boundary_prefix_check=start_depths:87866,87867;words:00000000,cf3030cf;"
           "edges:12133;end_depths:99999,100000;words:be79924b,90f58380\n");
}

static void exact_long_orbit(void) {
    uint32_t a = UINT32_C(0xbe79924b);
    uint32_t b = UINT32_C(0x90f58380);
    uint64_t xor_digest = 0;
    uint64_t sum_digest = 0;
    uint64_t fnv_digest = UINT64_C(1469598103934665603);
    const uint64_t checkpoints[] = {
        UINT64_C(1), UINT64_C(10), UINT64_C(1000), UINT64_C(1000000),
        UINT64_C(100000000), UINT64_C(1000000000), UINT64_C(1420778968)
    };
    const uint32_t checkpoint_a[] = {
        UINT32_C(0x90f58380), UINT32_C(0xbc973f1e), UINT32_C(0x73703da9),
        UINT32_C(0xe7f1f52a), UINT32_C(0xd4f2a1bc), UINT32_C(0x36363150),
        UINT32_C(0xdae372fa)
    };
    const uint32_t checkpoint_b[] = {
        UINT32_C(0xcb08e372), UINT32_C(0x25ec4c68), UINT32_C(0x07f4bd13),
        UINT32_C(0xd37751d3), UINT32_C(0x7e8f9310), UINT32_C(0xc4144da0),
        UINT32_C(0x00000000)
    };
    size_t checkpoint_index = 0;

    for (uint64_t step = 1; step <= EXPECTED_STEPS; ++step) {
        require_true(b != 0, "no earlier zero word");
        const uint32_t x = solve_after_reset(a, b, 32);
        a = b;
        b = x;

        const uint64_t packed = ((uint64_t)a << 32) | (uint64_t)b;
        xor_digest ^= rotate_left_64(packed, (unsigned)step)
                      ^ step * UINT64_C(0x9e3779b97f4a7c15);
        sum_digest += packed + step * UINT64_C(0xd6e8feb86659fd93);
        fnv_digest ^= packed ^ step;
        fnv_digest *= UINT64_C(1099511628211);

        if (checkpoint_index < sizeof(checkpoints) / sizeof(checkpoints[0])
            && step == checkpoints[checkpoint_index]) {
            require_true(
                a == checkpoint_a[checkpoint_index]
                && b == checkpoint_b[checkpoint_index],
                "frozen orbit checkpoint"
            );
            printf(
                "orbit_checkpoint=step:%" PRIu64 ";a:%08" PRIx32
                ";b:%08" PRIx32 "\n",
                step, a, b
            );
            ++checkpoint_index;
        }
        if (b == 0) {
            require_true(step == EXPECTED_STEPS, "zero step");
            break;
        }
    }
    require_true(checkpoint_index == sizeof(checkpoints) / sizeof(checkpoints[0]),
                 "all orbit checkpoints");
    require_true(b == 0, "expected terminal zero");
    require_true(a == EXPECTED_PREDECESSOR, "expected terminal predecessor");
    require_true(xor_digest == UINT64_C(0x25469bdaa1c33b27), "xor digest");
    require_true(sum_digest == UINT64_C(0x65b743c483016d28), "sum digest");
    require_true(fnv_digest == UINT64_C(0x46485ada15e2f5b9), "fnv digest");

    const unsigned predecessor_weight = (unsigned)__builtin_popcount(a);
    const unsigned predecessor_period = least_period(a, 32);
    const uint32_t candidate_zero = zero_successor_candidate(a, 0);
    const uint32_t candidate_one = zero_successor_candidate(a, 1);
    require_true(predecessor_weight == 20, "predecessor weight");
    require_true(predecessor_period == 32, "predecessor least period");
    require_true(candidate_one == ~candidate_zero, "candidate complement");
    require_true(least_period(candidate_zero, 32) == 32, "candidate zero period");
    require_true(least_period(candidate_one, 32) == 32, "candidate one period");

    printf("parent_state=depths:99999,100000;words:%08" PRIx32 ",%08" PRIx32 "\n",
           UINT32_C(0xbe79924b), UINT32_C(0x90f58380));
    printf("transfer_steps_to_zero=%" PRIu64 "\n", EXPECTED_STEPS);
    printf("next_zero_depth=%" PRIu64 "\n", UINT64_C(100000) + EXPECTED_STEPS);
    printf("previous_zero_depth=87866;exact_zero_gap=1420791102\n");
    printf("boundary_path_edges=1420791101;prefix_edges_to_parent_state=12133\n");
    printf("terminal_predecessor=%08" PRIx32 ";weight=%u;least_period=%u\n",
           a, predecessor_weight, predecessor_period);
    printf("next_successor_candidates=%08" PRIx32 ",%08" PRIx32
           ";physical_phase=unresolved;least_period=32\n",
           candidate_zero, candidate_one);
    printf("orbit_xor_digest=%016" PRIx64 "\n", xor_digest);
    printf("orbit_sum_digest=%016" PRIx64 "\n", sum_digest);
    printf("orbit_fnv_digest=%016" PRIx64 "\n", fnv_digest);
    printf("period32_crude_next_zero_gap_bound=18446744065119617027\n");
}

int main(void) {
    printf("RULE30_TEMPORAL_CYLINDER_TRANSFER_THM4065\n");
    printf("solver=first_reset_cyclic_walk;word_bits=absolute_time_residues\n");
    printf("universal_map=T_P(A,B)=(B,X);domain:B_nonzero;codomain:first_nonzero\n");
    audit_complete_graph(1);
    audit_complete_graph(2);
    audit_complete_graph(4);
    audit_complete_graph(8);
    audit_boundary_prefix();
    exact_long_orbit();
    printf("parent=THM4047_fixed_physical_bank;parent_phase_packet_already_audited\n");
    printf("scope=PROVED universal transfer and infinite isolated zeros;FINITE-EXACT next zero;moving center and Rule30 prizes OPEN\n");
    return 0;
}
