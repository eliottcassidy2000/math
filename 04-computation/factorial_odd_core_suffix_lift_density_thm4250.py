#!/usr/bin/env python3
"""Primary exact companion for THM-4250's odd-core suffix language.

No repository modules are imported. The checker uses direct modular products,
a separate digit/carry transducer, exact lift prediction, finite-state sequence
certificates, density hostiles, and a frozen odd-core atlas.
"""

from collections import Counter
from collections import deque
from hashlib import sha256


MAX_B = 63
LIFTS = 4


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def ell0(b):
    """Least ell with 2^(ell-1)>=b."""
    return (b - 1).bit_length() + 1


def overlap_masks(b, ell, r):
    modulus = 1 << ell
    return tuple(
        ((s * r) % modulus) & (((b - s) * r) % modulus)
        for s in range(1, (b + 1) // 2)
    )


def direct_good(b, ell, r):
    return all(overlap_masks(b, ell, r))


def closure_set(b, ell):
    return tuple(r for r in range(1, 1 << ell, 2) if direct_good(b, ell, r))


def transducer_good(b, ell, r):
    """Digit/carry path independent of modular multiplication."""
    carries = [0] * b
    seen = [False] * ((b + 1) // 2)
    for bit_index in range(ell):
        bit = (r >> bit_index) & 1
        digits = [0] * b
        for u in range(1, b):
            total = carries[u] + u * bit
            digits[u] = total & 1
            carries[u] = total >> 1
            require(0 <= carries[u] < u, "carry escaped its finite state set")
        for s in range(1, (b + 1) // 2):
            seen[s] = seen[s] or bool(digits[s] & digits[b - s])
    return all(seen[1:])


def lift_prediction(b, ell, r):
    """Return the exact set of closing epsilon in r+epsilon*2^ell."""
    modulus = 1 << ell
    failures = []
    for s in range(1, (b + 1) // 2):
        t = b - s
        x = (s * r) % modulus
        y = (t * r) % modulus
        if x & y:
            continue
        even, odd = (s, t) if s % 2 == 0 else (t, s)
        h_even = (even * r // modulus) & 1
        h_odd = (odd * r // modulus) & 1
        failures.append((h_even, h_odd))
    if not failures:
        return (0, 1)
    if any(h_even == 0 for h_even, _ in failures):
        return ()
    odd_bits = {h_odd for _, h_odd in failures}
    if len(odd_bits) != 1:
        return ()
    h_odd = next(iter(odd_bits))
    return (1 - h_odd,)


def witness_hitting_number(masks, ell):
    """Minimum number of bit positions hitting every nonzero overlap mask."""
    require(all(masks), "hitting number requested for a nonclosing row")
    for size in range(1, ell + 1):
        for bits in range(1, 1 << ell):
            if bits.bit_count() == size and all(mask & bits for mask in masks):
                return size
    raise RuntimeError("no hitting set")


def fibonacci(n):
    a, c = 0, 1
    for _ in range(n):
        a, c = c, a + c
    return a


def all_one_run_forces(u, length):
    for initial_carry in range(u):
        carry = initial_carry
        output = None
        for _ in range(length):
            total = carry + u
            output = total & 1
            carry = total >> 1
        if output != 1:
            return False
    return True


def audit_general_theorems():
    direct_transducer = 0
    lift_cells = 0
    interval_cells = 0
    run_states = 0
    for b in range(3, MAX_B + 1, 2):
        base = ell0(b)
        # Full independent direct/transducer agreement through three lifts.
        for ell in range(base, base + 4):
            modulus = 1 << ell
            for r in range(1, modulus, 2):
                require(
                    direct_good(b, ell, r) == transducer_good(b, ell, r),
                    f"transducer mismatch b={b} ell={ell} r={r}",
                )
                predicted = lift_prediction(b, ell, r)
                actual = tuple(
                    epsilon
                    for epsilon in (0, 1)
                    if direct_good(b, ell + 1, r + epsilon * modulus)
                )
                require(predicted == actual, f"lift mismatch b={b} ell={ell} r={r}")
                lift_cells += 1
                direct_transducer += 1

        # The negative interval M-u theorem, through six lifts.
        for ell in range(base, base + 7):
            modulus = 1 << ell
            cutoff = modulus // (2 * (b - 1))
            for u in range(1, cutoff + 1, 2):
                require(direct_good(b, ell, modulus - u), "negative interval failed")
                interval_cells += 1

        # Carry-reset lemma behind exponential density one.
        length = (b - 2).bit_length() + 1  # ceil(log2(b-1))+1
        for u in range(1, b):
            require(all_one_run_forces(u, length), f"run reset failed b={b} u={u}")
            run_states += u
        require(
            not all_one_run_forces(b - 1, length - 1),
            f"claimed sharp run length was not sharp b={b}",
        )

    return direct_transducer, lift_cells, interval_cells, run_states


def atlas():
    semantic = []
    total_minimal = 0
    max_tau = 0
    max_tau_witness = None
    print("ODD_CORE_SUFFIX_CLOSURE_ATLAS")
    print("definition=low-ell residues of sr and (b-s)r have nonzero bit overlap for every complementary pair")
    for b in range(3, MAX_B + 1, 2):
        base = ell0(b)
        sets = [closure_set(b, base + shift) for shift in range(LIFTS + 1)]
        counts = tuple(len(values) for values in sets)
        new = tuple(counts[i + 1] - 2 * counts[i] for i in range(LIFTS))
        require(all(x >= 0 for x in new), "density monotonicity failed")
        hit_hist = Counter()
        for r in sets[0]:
            tau = witness_hitting_number(overlap_masks(b, base, r), base)
            hit_hist[tau] += 1
            if tau > max_tau:
                max_tau = tau
                max_tau_witness = (b, base, r, overlap_masks(b, base, r))
        total_minimal += counts[0]
        row = (
            b,
            base,
            counts,
            new,
            tuple(sorted(hit_hist.items())),
            sets[0],
        )
        semantic.append(repr(row))
        print(
            f"b={b:02d} ell0={base} counts={counts} new_one_lifts={new} "
            f"hit={dict(sorted(hit_hist.items()))} R0={sets[0]}"
        )
    digest = sha256("\n".join(semantic).encode()).hexdigest()
    print(
        f"atlas_summary=odd_b_3..{MAX_B} minimal_closures={total_minimal} "
        f"max_minimum_hitting_number={max_tau} witness={max_tau_witness} digest={digest}"
    )
    return digest


def specializations_and_hostiles():
    # Exact b=3 Fibonacci law.
    for ell in range(1, 23):
        expected = (1 << (ell - 1)) - fibonacci(ell)
        require(len(closure_set(3, ell)) == expected, "b=3 Fibonacci law failed")
    print("b3_exact=|R_3(ell)|=2^(ell-1)-F_ell for F_1=F_2=1 verified_ell=1..22")

    for b in (3, 5, 7):
        base = ell0(b)
        rows = tuple(len(closure_set(b, base + j)) for j in range(9))
        print(f"special_b={b} ell0={base} nine_counts={rows}")

    # No monotonicity in b, no negation or inversion symmetry.
    require(direct_good(3, 4, 3) and not direct_good(5, 4, 3), "b hostile 3 failed")
    require(direct_good(5, 4, 5) and not direct_good(3, 4, 5), "b hostile 5 failed")
    require(direct_good(3, 3, 3) and not direct_good(3, 3, 5), "negation hostile failed")
    require(direct_good(3, 4, 13) and not direct_good(3, 4, pow(13, -1, 16)), "inverse hostile failed")

    # Low suffix failure is not all-height failure: 51 mod 64 fails for b=5,
    # while its extension 1331 closes once the decisive high bits are retained.
    require(not direct_good(5, 6, 51), "low representative unexpectedly closes")
    require(1331 % 64 == 51 and direct_good(5, 11, 1331), "high extension hostile failed")

    # A failed base can have zero or exactly one newly closing lift; a closed
    # base always has both.
    require(lift_prediction(3, 3, 5) == (1,), "one-lift hostile failed")
    require(lift_prediction(3, 3, 1) == (), "zero-lift hostile failed")
    require(lift_prediction(3, 3, 3) == (0, 1), "two-lift inheritance failed")

    print(
        "hostiles=R3(4)_and_R5(4)_incomparable; negation_and_inverse_fail; "
        "b5_suffix_51_mod64_fails_at_51_but_extension_1331_closes_at_ell11; "
        "b3_base_lifts r1->0 r5->1 r3->2"
    )


def transition_state(b, state, bit):
    carries, flags = state
    digits = []
    following = []
    for u, carry in enumerate(carries, 1):
        total = carry + u * bit
        digits.append(total & 1)
        following.append(total >> 1)
    for s in range(1, (b + 1) // 2):
        if digits[s - 1] & digits[b - s - 1]:
            flags |= 1 << (s - 1)
    return tuple(following), flags


def reachable_automaton(b):
    start = transition_state(b, (tuple([0] * (b - 1)), 0), 1)
    queue = deque([start])
    states = [start]
    index = {start: 0}
    transitions = []
    while queue:
        state = queue.popleft()
        row = []
        for bit in (0, 1):
            following = transition_state(b, state, bit)
            if following not in index:
                index[following] = len(states)
                states.append(following)
                queue.append(following)
            row.append(index[following])
        transitions.append(tuple(row))
    require(len(transitions) == len(states), "incomplete BFS")
    return states, transitions


def automaton_sequence(b, length):
    states, transitions = reachable_automaton(b)
    dp = [0] * len(states)
    dp[0] = 1
    full = (1 << ((b - 1) // 2)) - 1
    answer = []
    for _ in range(length):
        answer.append(
            sum(count for count, (_, flags) in zip(dp, states) if flags == full)
        )
        following = [0] * len(states)
        for source, count in enumerate(dp):
            for target in transitions[source]:
                following[target] += count
        dp = following
    # Entry zero is ell=1, because the mandatory low bit has already run.
    return states, answer


def exact_small_closed_forms():
    # Coefficients C mean a_n+C1*a_(n-1)+...+Cd*a_(n-d)=0.
    data = {
        3: (
            (1, -3, 1, 2),
            "(x-2)(x^2-x-1)",
            ((1, -2), (1, -1, -1)),
        ),
        5: (
            (1, -5, 8, -3, -6, 15, -17, 5, 5, -10, 8, 0, -1, 2),
            "(x-2)(x^2+1)(x^2-x-1)(x^4-x^3-1)(x^4-x^3-x^2+x-1)",
            (
                (1, -2),
                (1, 0, 1),
                (1, -1, -1),
                (1, -1, 0, 0, -1),
                (1, -1, -1, 1, -1),
            ),
        ),
        7: (
            (1, -3, 1, -1, 7, 2, -2, -11, -5, 1, 7, 5, 2),
            "(x-2)(x^3-x-1)(x^3-x^2-1)(x^5-2x^2-2x-1)",
            (
                (1, -2),
                (1, 0, -1, -1),
                (1, -1, 0, -1),
                (1, 0, 0, -2, -2, -1),
            ),
        ),
    }
    for b, (coefficients, factorization, factors) in data.items():
        product = [1]
        for factor in factors:
            following = [0] * (len(product) + len(factor) - 1)
            for i, x in enumerate(product):
                for j, y in enumerate(factor):
                    following[i + j] += x * y
            product = following
        require(tuple(product) == coefficients, f"factor expansion failed b={b}")
        states, sequence = automaton_sequence(b, 120)
        degree = len(coefficients) - 1
        discrepancies = []
        for n in range(degree, degree + len(states)):
            discrepancy = sum(
                coefficients[j] * sequence[n - j]
                for j in range(degree + 1)
            )
            discrepancies.append(discrepancy)
        require(not any(discrepancies), f"small recurrence failed b={b}")
        # Independent direct enumeration over the feasible prefix.
        for ell in range(1, 17):
            require(sequence[ell - 1] == len(closure_set(b, ell)), "sequence mismatch")
        print(
            f"cfinite_b={b} reachable_states={len(states)} recurrence_degree={degree} "
            f"characteristic={factorization} cayley_hamilton_zero_checks={len(states)}"
        )


def density_controls():
    rows = []
    for b in (3, 5, 7, 15, 31, 63):
        length = (b - 2).bit_length() + 1
        for ell in (ell0(b), ell0(b) + 4, ell0(b) + 8):
            total = 1 << (ell - 1)
            actual = len(closure_set(b, ell))
            blocks = (ell - 1) // length
            # Integer form of actual/total >= 1-(1-2^-L)^blocks.
            rhs_num = (1 << (length * blocks)) - ((1 << length) - 1) ** blocks
            rhs_den = 1 << (length * blocks)
            require(actual * rhs_den >= total * rhs_num, "density bound failed")
            rows.append((b, ell, actual, total, length, blocks, rhs_num, rhs_den))
    print(f"density_controls={rows}")


def main():
    print("THM-4250 PRIMARY ODD-CORE SUFFIX LIFT/DENSITY COMPANION")
    direct, lifts, intervals, run_states = audit_general_theorems()
    print(
        f"general_audit=PASS direct_transducer_cells={direct} lift_cells={lifts} "
        f"negative_interval_cells={intervals} carry_initial_states={run_states}"
    )
    atlas()
    specializations_and_hostiles()
    exact_small_closed_forms()
    density_controls()
    print("verdict=PASS")


if __name__ == "__main__":
    main()
