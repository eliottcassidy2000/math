"""Independent exact audit; checks remain enabled under python -O."""

import hashlib
import itertools
import json
import random
from math import gcd
from collections import Counter
from pathlib import Path

from cut_response import compile_response, endpoint_tables, extension, from_code


COUNTS = Counter()
DIGEST = hashlib.sha256()


def check(condition, label, data=None):
    if not condition:
        raise RuntimeError((label, data))
    COUNTS[label] += 1


def direct_h(out):
    """Independent forward-extension DP, without compiler endpoint tables."""
    n = len(out)
    if not n:
        return 1
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for u in range(n):
        dp[1 << u][u] = 1
    for mask in range(1, full + 1):
        for u, count in enumerate(dp[mask]):
            if count:
                available = out[u] & (full ^ mask)
                while available:
                    bit = available & -available
                    available ^= bit
                    dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[full])


def literal_data(out):
    """Permutations and every possible break, no subset endpoint DP."""
    n = len(out)
    h = 0
    starts, ends = [0] * n, [0] * n
    covers = [[0] * n for _ in range(n)]
    for path in itertools.permutations(range(n)):
        valid = [bool(out[u] >> v & 1) for u, v in zip(path, path[1:])]
        if all(valid):
            h += 1
            if n:
                starts[path[0]] += 1
                ends[path[-1]] += 1
        for split in range(1, n):
            if all(valid[:split - 1]) and all(valid[split:]):
                covers[path[split - 1]][path[split]] += 1
    return h, starts, ends, covers


def component_count(out):
    """SCC count from independent transitive closure."""
    n = len(out)
    reach = [out[u] | 1 << u for u in range(n)]
    for k in range(n):
        for u in range(n):
            if reach[u] >> k & 1:
                reach[u] |= reach[k]
    remaining = (1 << n) - 1
    count = 0
    while remaining:
        u = (remaining & -remaining).bit_length() - 1
        component = sum(1 << v for v in range(n)
                        if reach[u] >> v & 1 and reach[v] >> u & 1)
        remaining &= ~component
        count += 1
    return count


def audit(out, literal=False, all_constraints=False):
    n = len(out)
    full = (1 << n) - 1
    response = compile_response(out)
    COUNTS["parents"] += 1
    check(response.h == direct_h(out), "base_direct_dp")
    if literal:
        h, starts, ends, covers = literal_data(out)
        check((h, starts, ends, covers) ==
              (response.h, response.starts, response.ends, response.covers),
              "literal_permutation_and_cover_tables", out)
    end, _ = endpoint_tables(out)
    for u in range(n):
        lhs = response.starts[u] - response.ends[u]
        rhs = sum(response.covers[u][v] for v in range(n) if out[u] >> v & 1)
        rhs -= sum(response.covers[v][u] for v in range(n) if out[v] >> u & 1)
        check(lhs == rhs, "endpoint_divergence")
        check(2 * response.degrees[u] == lhs + sum(response.covers[v][u]
                                                  for v in range(n)),
              "canonical_outdegree")
        for v in range(u + 1, n):
            rest = full ^ (1 << u) ^ (1 << v)
            lower = sum(end[rest]) if rest else 1
            check(response.weights[u][v] >= lower, "pair_deletion_lower_bound")
            # W does not depend on the orientation of its own endpoint pair.
            if n <= 4:
                flipped = list(out)
                flipped[u] ^= 1 << v
                flipped[v] ^= 1 << u
                alt = compile_response(flipped)
                check((response.covers[u][v], response.covers[v][u]) ==
                      (alt.covers[u][v], alt.covers[v][u]), "own_arc_invariance")
    values = response.all_values()
    for cut in range(1 << n):
        actual = direct_h(extension(out, cut))
        check(values[cut] == actual == response.value(cut), "all_cut_direct_dp",
              (out, cut, values[cut], actual))
        if n <= 4:
            check(actual == literal_data(extension(out, cut))[0],
                  "extension_literal_permutations")
        cross = sum(1 for u in range(n) if cut >> u & 1 for v in range(n)
                    if out[u] >> v & 1 and not (cut >> v & 1))
        check(actual >= response.h + 2 * cross, "crossing_arc_lower_bound")
        check((actual == response.h) == (cross == 0), "zero_cut_equivalence")
    check(values.count(response.h) == component_count(out) + 1,
          "scc_minimizer_count")
    check(response.gain_content() == gcd(*(value - response.h for value in values)),
          "gain_gcd_content")
    expected_positive = min((value for value in values if value > response.h), default=None)
    positive = response.minimum_positive()
    check((positive is None and expected_positive is None) or
          (positive is not None and positive[0] == expected_positive and
           values[positive[1]] == expected_positive), "minimum_positive_cut")
    if n <= 5:
        for u in range(n):
            for v in range(u + 1, n):
                for mask in range(1 << n):
                    if not (mask & ((1 << u) | (1 << v))):
                        delta = (values[mask | 1 << u] + values[mask | 1 << v]
                                 - values[mask] - values[mask | 1 << u | 1 << v])
                        check(delta == 2 * response.weights[u][v],
                              "mixed_second_difference")
    constraints = itertools.product(range(3), repeat=n) if all_constraints else [
        tuple(1 if u == 0 else 2 if u == n - 1 else 0 for u in range(n)),
        tuple(2 if u % 3 == 0 else 1 if u % 3 == 1 else 0 for u in range(n))]
    for states in constraints:
        inside = sum(1 << u for u, state in enumerate(states) if state == 1)
        outside = sum(1 << u for u, state in enumerate(states) if state == 2)
        expected = min(values[cut] for cut in range(1 << n)
                       if (cut & inside) == inside and not (cut & outside))
        actual, cut = response.minimum(inside, outside)
        check(actual == expected and values[cut] == expected,
              "constrained_min_cut", (out, inside, outside))
    DIGEST.update(json.dumps([out, response.h, response.weights, values],
                            separators=(",", ":")).encode() + b"\n")


def block_insertion_audit():
    # This isolates the lower-bound proof from all tournament path enumeration.
    # Every possible 0/1 endpoint relation word, lengths zero through six.
    for length in range(7):
        for states in itertools.product(range(4), repeat=length):
            a = [1] + [state & 1 for state in states] + [0]
            b = [1] + [state >> 1 for state in states] + [0]
            forward = sum(a[i] * (1 - b[i + 1]) + b[i] * (1 - a[i + 1])
                          for i in range(length + 1))
            backward = sum(a[i + 1] * (1 - b[i]) + b[i + 1] * (1 - a[i])
                           for i in range(length + 1))
            check(forward == 2 + backward, "block_insertion_telescoping")


def main():
    for n in range(6):
        for code in range(1 << (n * (n - 1) // 2)):
            audit(from_code(n, code), literal=True, all_constraints=n <= 4)
    rng = random.Random(124031)
    random_counts = {}
    for n in range(6, 10):
        random_counts[n] = 12
        for _ in range(12):
            code = rng.randrange(1 << (n * (n - 1) // 2))
            audit(from_code(n, code))
    block_insertion_audit()
    # Hostile controls: wrong pair orientation/factor would already fail the
    # directed 3-cycle extension of 0->1; separately exercise bad inputs.
    for bad in [(1,), (0, 0), (2, 0, 0), (-1,)]:
        try:
            compile_response(bad)
        except ValueError:
            COUNTS["malformed_tournament_rejected"] += 1
        else:
            raise RuntimeError(("accepted malformed tournament", bad))
    response = compile_response((2, 0))
    check(response.all_values() == [1, 3, 1, 1], "orientation_and_factor_control")
    try:
        response.minimum(1, 1)
    except ValueError:
        COUNTS["contradictory_constraints_rejected"] += 1
    else:
        raise RuntimeError("accepted contradictory forced memberships")
    print(json.dumps({
        "status": "PASS", "universe": {
            "all_labelled_parent_orders": [0, 1, 2, 3, 4, 5],
            "random_seed": 124031, "random_parent_counts": random_counts,
            "all_cuts_for_every_parent": True,
            "literal_cover_tables_through_parent_order": 5,
            "literal_extension_paths_through_parent_order": 4,
            "all_forced_membership_assignments_through_order": 4},
        "checks": dict(sorted(COUNTS.items())),
        "data_sha256": DIGEST.hexdigest(),
        "source_sha256": {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                          for name in ["cut_response.py", "verify.py"]}},
        sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
