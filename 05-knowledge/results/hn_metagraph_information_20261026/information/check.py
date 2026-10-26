"""Exact information-boundary checks; standard library only, no asserts.

Universe: every labelled tournament on n=1,...,6. No inherited filters.
Independent H implementation: subset dynamic programming, checked against
literal permutation enumeration through n=5. All arithmetic is integral or
Fraction. The proof is in RESULT.md; this finite census does not replace it.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial
import json

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def arcs(n):
    return list(combinations(range(n), 2))


def matrix(n, mask):
    a = [[0] * n for _ in range(n)]
    for bit, (i, j) in enumerate(arcs(n)):
        a[i][j] = 1 if not (mask >> bit & 1) else -1
        a[j][i] = -a[i][j]
    return a


def encode(a):
    return sum((a[i][j] < 0) << bit
               for bit, (i, j) in enumerate(arcs(len(a))))


def h_dp(a, endpoints=False):
    n = len(a)
    out = [sum((a[i][j] == 1) << j for j in range(n))
           for i in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for subset in range(1, 1 << n):
        for last in range(n):
            value = dp[subset][last]
            if not value:
                continue
            rest = out[last] & ~subset
            while rest:
                bit = rest & -rest
                rest -= bit
                nxt = bit.bit_length() - 1
                dp[subset | bit][nxt] += value
    return tuple(dp[-1]) if endpoints else sum(dp[-1])


def fwht(values):
    result = list(values)
    width = 1
    while width < len(result):
        for start in range(0, len(result), 2 * width):
            for j in range(start, start + width):
                x, y = result[j], result[j + width]
                result[j], result[j + width] = x + y, x - y
        width *= 2
    return result


def curl(a):
    return tuple(a[0][i] + a[i][j] + a[j][0]
                 for i, j in combinations(range(1, len(a)), 2))


def scores(a):
    return tuple(sum(value == 1 for value in row) for row in a)


def wedge_sign(e, f):
    shared = set(e) & set(f)
    if not shared:
        return 0
    middle = next(iter(shared))
    left = next(v for v in e if v != middle)
    right = next(v for v in f if v != middle)
    return (1 if left < middle else -1) * (1 if middle < right else -1)


def line_tree(edges):
    if not edges:
        return []
    seen, queue, result = {0}, [0], []
    for u in queue:
        for v in range(len(edges)):
            if v not in seen and set(edges[u]) & set(edges[v]):
                seen.add(v)
                queue.append(v)
                result.append((u, v))
    require(len(seen) == len(edges), "line graph disconnected")
    return result


def signature(mask, tree):
    return tuple(((mask >> e) ^ (mask >> f)) & 1 for e, f in tree)


def decode(signature_bits, tree, m, root_bit=0):
    bits = [None] * m
    if not m:
        return 0
    bits[0] = root_bit
    for (u, v), value in zip(tree, signature_bits):
        require(bits[u] is not None and bits[v] is None, "tree order")
        bits[v] = bits[u] ^ value
    return sum(bit << i for i, bit in enumerate(bits))


def permuted(a, order):
    return [[a[order[i]][order[j]] for j in range(len(a))]
            for i in range(len(a))]


def rational(x):
    return str(x)


def main():
    records = []
    history = {}
    for n in range(1, 7):
        edges = arcs(n)
        m, size = len(edges), 1 << len(edges)
        full = size - 1
        matrices = [matrix(n, mask) for mask in range(size)]
        ending_counts = [h_dp(a, endpoints=True) for a in matrices]
        hs = [sum(counts) for counts in ending_counts]
        if n <= 5:
            orders = list(permutations(range(n)))
            for a, h in zip(matrices, hs):
                literal = sum(all(a[p[i]][p[i + 1]] == 1
                                  for i in range(n - 1)) for p in orders)
                require(h == literal, "DP/permutation H disagreement")
        spectrum = fwht(hs)
        require(spectrum[0] == factorial(n) * size // (1 << (n - 1)),
                "mean H")
        require(all(spectrum[s] == 0 for s in range(size)
                    if s.bit_count() % 2), "odd Walsh support")
        for e, f in combinations(range(m), 2):
            eps = wedge_sign(edges[e], edges[f])
            expected = Fraction(eps * 2 * factorial(n - 2) * size,
                                1 << (n - 1))
            require(spectrum[(1 << e) | (1 << f)] == expected,
                    "quadratic coefficient")
        tree = line_tree(edges)
        fibres, joint = defaultdict(list), set()
        curl_groups = defaultdict(list)
        for mask, a in enumerate(matrices):
            sig = signature(mask, tree)
            fibres[sig].append(mask)
            joint.add((curl(a), sig))
            curl_groups[curl(a)].append(mask)
            reconstructed = decode(sig, tree, m)
            require(mask in (reconstructed, reconstructed ^ full),
                    "tree reconstruction")
        if m:
            require(len(fibres) == size // 2, "feature fibre count")
            require(all(len(v) == 2 and v[0] ^ v[1] == full
                        for v in fibres.values()), "converse fibres")
        else:
            require(len(fibres) == 1, "singleton boundary")
        if n >= 3:
            require(len(joint) == size, "aligned joint invariant")
        # Direct full response equality, independent of the Fourier proof.
        if n <= 4:
            response_groups = defaultdict(list)
            for mask in range(size):
                response_groups[tuple(hs[mask ^ flip] for flip in range(size))].append(mask)
            require(sorted(response_groups.values()) == sorted(fibres.values()),
                    "literal response versus binary features")
        # Direct shifted-correlation check, with exhaustive shifts through n4.
        shift_samples = range(size) if n <= 4 else (0, 1, size // 3, full)
        for mask in shift_samples:
            for e, f in tree:
                support = (1 << e) | (1 << f)
                direct = sum((1 if (flip & support).bit_count() % 2 == 0 else -1)
                             * hs[mask ^ flip] for flip in range(size))
                sign = 1 if (mask & support).bit_count() % 2 == 0 else -1
                require(direct == sign * spectrum[support], "shift correlation")
        collisions = []
        for group in curl_groups.values():
            first = defaultdict(list)
            for mask in group:
                first[tuple(hs[mask ^ (1 << e)] for e in range(m))].append(mask)
            collisions.extend(v for v in first.values() if len(v) > 1)
            for root in range(n):
                # Path reversal turns starts into ends of the converse.
                candidates = [mask for mask in group
                              if ending_counts[mask ^ full][root] > 0]
                require(len(candidates) == 1, "valid path-start selects curl cut")
        record = {
            "n": n, "labelled_tournaments": size, "curl_fibres": len(curl_groups),
            "response_converse_fibres": len(fibres), "binary_features": len(tree),
            "aligned_joint_fibres": len(joint),
            "curl_plus_single_flip_collision_fibres": len(collisions),
        }
        if n <= 5:
            iso, curl_iso, separate, aligned = set(), set(), set(), set()
            iso_self_converse = set()
            for mask, a in enumerate(matrices):
                pa = [permuted(a, p) for p in orders]
                codes = [encode(b) for b in pa]
                iso_code = min(codes)
                opposite_code = min(code ^ full for code in codes)
                k_code = min(curl(b) for b in pa)
                q_code = min(iso_code, opposite_code)
                iso.add(iso_code)
                curl_iso.add(k_code)
                separate.add((k_code, q_code))
                aligned.add(min((curl(b), signature(code, tree))
                                for b, code in zip(pa, codes)))
                if iso_code == opposite_code:
                    iso_self_converse.add(iso_code)
            if n >= 3:
                require(len(aligned) == len(iso), "joint relabelling")
            record.update({"tournament_isotypes": len(iso),
                           "self_converse_isotypes": len(iso_self_converse),
                           "curl_isotypes": len(curl_iso),
                           "separate_quotient_pairs": len(separate),
                           "aligned_joint_isotypes": len(aligned)})
        if n in (3, 5):
            mean = Fraction(sum(hs), size)
            variance = sum((h - mean) ** 2 for h in hs) / size
            third = sum((h - mean) ** 3 for h in hs) / size
            grouped = defaultdict(list)
            for a, h in zip(matrices, hs):
                grouped[scores(a)].append(h)
            projected_variance = sum(len(v) * (Fraction(sum(v), len(v)) - mean) ** 2
                                     for v in grouped.values()) / size
            quadratic_energy = Fraction(sum(v * v for s, v in enumerate(spectrum)
                                            if s.bit_count() == 2), size * size)
            history[str(n)] = {
                "H_distribution": dict(sorted(Counter(hs).items())),
                "mean": rational(mean), "variance": rational(variance),
                "Fisher_derivative_at_zero": rational(third),
                "variance_of_score_conditional_mean": rational(projected_variance),
                "quadratic_Walsh_energy": rational(quadratic_energy),
                "OCR": rational(projected_variance / variance),
                "quadratic_fraction": rational(quadratic_energy / variance),
            }
            if n == 5:
                require(projected_variance != quadratic_energy, "historical OCR hostile")
                require(third == Fraction(165, 32), "n5 Fisher hostile")
                require(mean + Fraction(sum(v for s, v in enumerate(spectrum)
                                            if s.bit_count() == 2), size) == 0,
                        "transitive degree-zero-plus-two value")
                require(hs[0] == 1, "transitive H")
        records.append(record)
    # Explicit minimal labelled operation hostile.
    a = matrix(3, 0)  # 0 -> 1 -> 2, 0 -> 2.
    b = matrix(3, 3)  # 1 -> 2 -> 0, 1 -> 0.
    require(curl(a) == curl(b), "minimal same curl")
    require(h_dp(matrix(3, 1)) == 1 and h_dp(matrix(3, 2)) == 3,
            "same named arc flip H=1/3")
    # Explicit smallest separate-isomorphism quotient alignment failure.
    a = matrix(4, 2)  # 0->1->2->0, all three ->3.
    b = matrix(4, 2 ^ ((1 << 2) | (1 << 4) | (1 << 5)))
    require(curl(a) == curl(b), "C3/singleton curl rotation")
    require(sorted(scores(a)) == [0, 2, 2, 2] and
            sorted(scores(b)) == [1, 1, 1, 3], "nonisomorphism witness")
    require(any(encode(permuted(a, p)) ^ 63 == encode(b)
                for p in permutations(range(4))), "converse isomorphism")
    low = {mask for mask in range(8) if h_dp(matrix(3, mask)) == 1}
    low_edges = {(u, v) for u in low for v in low if u < v and (u ^ v).bit_count() == 1}
    require(len(low) == 6 and len(low_edges) == 6, "H=1 six-cycle cells")
    require(all(sum(v in edge for edge in low_edges) == 2 for v in low),
            "H=1 six-cycle degrees")
    reached, pending = {next(iter(low))}, [next(iter(low))]
    for u in pending:
        for v in low:
            if (u ^ v).bit_count() == 1 and v not in reached:
                reached.add(v)
                pending.append(v)
    require(reached == low, "H=1 sublevel connected")
    require(not any({base, base ^ (1 << i), base ^ (1 << j),
                     base ^ (1 << i) ^ (1 << j)} <= low
                    for base in range(8) for i, j in combinations(range(3), 2)),
            "H=1 sublevel no filled square")
    for dimension in range(1, 11):
        require(sum((-1) ** d * comb(dimension, d) * 2 ** (dimension - d)
                    for d in range(dimension + 1)) == 1, "filled cube Euler")
    result = {
        "status": "PASS", "checks": CHECKS,
        "universe": "Every labelled tournament n=1..6; no filters; isomorphism census n<=5",
        "records": records, "historical_controls": history,
        "small_hostiles": {"same_curl_masks_n3": [0, 3], "flip_edge": [0, 1],
                           "flipped_H": [1, 3], "separate_quotient_masks_n4": [2, 54],
                           "H1_sublevel_n3": {"vertices": 6, "edges": 6,
                                               "squares": 0, "betti0": 1, "betti1": 1}},
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
