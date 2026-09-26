"""Exact probes for the user's authoritative 35-colour triangular prefix.

No external packages and no floating-point arithmetic in any decision.
Run from the repository root; stdout is the retained experiment artifact.
"""

from math import isqrt


LISTED = "R K B R K R K B R K B R K R K B R K R K B R K B R K R K B R K B R K B".split()


def require(test, message):
    if not test:
        raise RuntimeError(message)


def floor_phi(k):
    return (k + isqrt(5 * k * k)) // 2


def blue_position(k):
    return 2 * floor_phi(k) + k


def beatty_colour_prefix(limit):
    blues = set()
    k = 1
    while blue_position(k) <= limit:
        blues.add(blue_position(k))
        k += 1
    result = []
    nonblue = 0
    for n in range(1, limit + 1):
        if n in blues:
            result.append("B")
        else:
            nonblue += 1
            result.append("R" if nonblue % 2 else "K")
    return result


def substitution_prefix(limit):
    rule = {"R": "RK", "K": "B", "B": "RK"}
    word = "R"
    while len(word) < limit:
        word = "".join(rule[c] for c in word)
    return list(word[:limit])


def fib_word_prefix(limit):
    word = "A"
    while len(word) < limit:
        word = "".join("AB" if c == "A" else "A" for c in word)
    return word[:limit]


def charge_b(n):
    m = n + 1
    return (3 * m - isqrt(5 * m * m) - 1) // 2


def charge(n):
    b = charge_b(n)
    return n - 2 * b, b


def charge_colour(n):
    return tuple(x % 2 for x in charge(n))


def zeckendorf(n):
    weights = [1, 2]
    while weights[-1] <= n:
        weights.append(sum(weights[-2:]))
    result = []
    for i in range(len(weights) - 1, -1, -1):
        if weights[i] <= n:
            result.append(i)
            n -= weights[i]
    require(n == 0, "Zeckendorf remainder")
    require(all(a - b > 1 for a, b in zip(result, result[1:])), "adjacent digits")
    return result


def collision(feature, colours):
    seen = {}
    for n, col in enumerate(colours, start=1):
        key = feature(n)
        if key in seen and seen[key][1] != col:
            return seen[key][0], n, key, seen[key][1], col
        seen[key] = n, col
    return None


def transition_collision(colours, transform, modulus):
    seen = {}
    for n in range(1, len(colours) + 1):
        tn = transform(n)
        if not 1 <= tn <= len(colours):
            continue
        state = colours[n - 1], n % modulus
        target = colours[tn - 1]
        if state in seen and seen[state][1] != target:
            return seen[state][0], n, state, seen[state][1], target
        seen[state] = n, target
    return None


def main():
    require(len(LISTED) == 35, "authoritative prefix length")
    blue = [i + 1 for i, c in enumerate(LISTED) if c == "B"]
    print("authoritative_colours", "".join(LISTED))
    print("authoritative_blue_positions", blue)
    print("nonblue_spaces", [b - a - 1 for a, b in zip(blue, blue[1:])])
    print("spaces_including_initial_boundary", [blue[0] - 1] + [b - a - 1 for a, b in zip(blue, blue[1:])])
    tokens = []
    pos = 0
    while pos < len(LISTED):
        if LISTED[pos] == "B":
            tokens.append("B")
            pos += 1
        else:
            require(LISTED[pos : pos + 2] == ["R", "K"], "red-black pairing")
            tokens.append("A")
            pos += 2
    word = "".join(tokens)
    print("paired_tokens", word)
    expected = fib_word_prefix(len(tokens))
    print("fibonacci_token_mismatches", [(i + 1, a, b) for i, (a, b) in enumerate(zip(word, expected)) if a != b])
    count_factors = [(i + 1, word[i : i + 5], word[i : i + 5].count("B")) for i in range(len(word) - 4)]
    lo = min(count_factors, key=lambda row: row[2])
    hi = max(count_factors, key=lambda row: row[2])
    require(hi[2] - lo[2] == 2, "mechanical-word hostile")
    print("unbalanced_length5_factors", lo, hi)

    limit = 100_000
    candidate = beatty_colour_prefix(limit)
    require(candidate == substitution_prefix(limit), "independent substitution/Beatty constructions")
    mismatches = [(i + 1, a, b) for i, (a, b) in enumerate(zip(LISTED, candidate)) if a != b]
    require(mismatches == [(35, "B", "R")], "all given terms comparison")
    print("candidate_colour_mismatches", mismatches)
    print("candidate_blue_positions", [blue_position(k) for k in range(1, 13)])
    print("independent_candidate_agreement_through", limit)

    features = {
        "inherited_four_state_charge": charge_colour,
        "lowest_Fibonacci_index_mod3": lambda n: min(zeckendorf(n)) % 3,
        "highest_Fibonacci_index_mod3": lambda n: max(zeckendorf(n)) % 3,
        "number_of_Fibonacci_terms_mod3": lambda n: len(zeckendorf(n)) % 3,
        "sum_of_Fibonacci_indices_mod3": lambda n: sum(zeckendorf(n)) % 3,
    }
    for name, feature in features.items():
        witness = collision(feature, LISTED)
        require(witness is not None, "unexpected fitting simple rule " + name)
        print("rejected_rule", name, "same_feature_different_listed_colours", witness)
    for modulus in (1, 2, 3, 5, 30, 60, 210):
        for name, transform in (("successor", lambda n: n + 1), ("triple_plus_one", lambda n: 3 * n + 1), ("shortcut", lambda n: n // 2 if n % 2 == 0 else (3 * n + 1) // 2)):
            witness = transition_collision(candidate, transform, modulus)
            require(witness is not None, "unexpected deterministic finite quotient")
            if modulus in (1, 30, 60):
                print("candidate_state_not_closed", name, "modulus", modulus, witness)

    carry_witnesses = {}
    charge_state_witness = None
    seen = {}
    for n in range(1, limit + 1):
        b = charge_b(n)
        d = charge_b(3 * n + 1) - 3 * b
        require(d in (-1, 0, 1, 2), "triple carry range")
        q = charge(n)
        qt = charge(3 * n + 1)
        require(qt == (3 * q[0] + 1 - 2 * d, 3 * q[1] + d), "exact triple carry")
        carry_witnesses.setdefault(d, n)
        state = charge_colour(n), n % 30
        target = charge_colour(3 * n + 1)
        if charge_state_witness is None and state in seen and seen[state][1] != target:
            charge_state_witness = seen[state][0], n, state, seen[state][1], target
        seen[state] = n, target
    require(set(carry_witnesses) == {-1, 0, 1, 2}, "all triple carry values")
    print("inherited_charge_triple_carries_first_sources", sorted(carry_witnesses.items()))
    print("charge_plus_mod30_not_closed", charge_state_witness)
    require(charge_state_witness is not None, "charge quotient hostile")

    weights = [1, 2]
    for _ in range(99):
        weights.append(sum(weights[-2:]))
    boundary_examples = []
    for k, fk in enumerate(weights):
        lower_indices = list(range(k - 1, -1, -2))
        require(sum(weights[j] for j in lower_indices) == fk - 1, "alternating lower decomposition")
        q = charge(fk)
        qm = charge(fk - 1)
        boundary = tuple(a - b for a, b in zip(q, qm))
        require(boundary == ((1, 0) if k % 2 == 0 else (-1, 1)), "boundary seed parity")
        if k < 8:
            boundary_examples.append((k, fk, [weights[j] for j in lower_indices], boundary))
    print("two_boundary_units_exact_through_Fibonacci_index", len(weights) - 1)
    print("boundary_unit_examples", boundary_examples)

    phase_states = [(charge_colour(weights[k]), tuple((a - b) % 2 for a, b in zip(charge(weights[k]), charge(weights[k] - 1)))) for k in range(6)]
    require(len(set(phase_states)) == 6, "six-state phase decoder")
    require(all((charge_colour(weights[k]), tuple((a - b) % 2 for a, b in zip(charge(weights[k]), charge(weights[k] - 1)))) == phase_states[k % 6] for k in range(len(weights))), "six-state phase periodicity")
    print("six_fibonacci_phase_states", phase_states)

    # Independent integer-only check of the all-odd family used in the
    # bounded-potential obstruction.  These are not claims about any
    # unverified continuation of the user's colour word.
    for k in (1, 2, 3, 10, 100, 1000):
        start = 2**k - 1
        x = start
        for j in range(1, k + 1):
            require(x % 2 == 1, "odd family parity")
            x = (3 * x + 1) // 2
            require(x == 3**j * 2 ** (k - j) - 1, "odd family formula")
            if k > 1:
                require(x > start, "odd family growth")
        require(x == 3**k - 1, "odd family endpoint")
    print("all_odd_growth_exact_controls", [1, 2, 3, 10, 100, 1000])
    print("PASS")


if __name__ == "__main__":
    main()
