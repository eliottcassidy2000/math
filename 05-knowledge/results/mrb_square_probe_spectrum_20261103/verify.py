#!/usr/bin/env python3
"""Exact source-preserving square probes; no external packages or producer imports."""
import hashlib
import itertools
import json
import math
from fractions import Fraction


def require(ok, context):
    if not ok:
        raise RuntimeError(context)


def square(n):
    return n >= 0 and math.isqrt(n) ** 2 == n


def birth(rows, incoming):
    """Old vertex i beats the newborn iff incoming contains i."""
    n = len(rows)
    return tuple(r | ((1 << n) if incoming >> i & 1 else 0)
                 for i, r in enumerate(rows)) + (((1 << n) - 1) ^ incoming,)


def h_dp(rows):
    n = len(rows)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, value in enumerate(dp[mask]):
            todo = rows[v] & ~mask if value else 0
            while todo:
                bit = todo & -todo
                todo -= bit
                dp[mask | bit][bit.bit_length() - 1] += value
    return sum(dp[-1])


def cut_coefficients(rows):
    """Boundary data from OLD permutations and their unique exposed bad step."""
    n = len(rows)
    start, end = [0] * n, [0] * n
    q = [[0] * n for _ in rows]
    h = 0
    for p in itertools.permutations(range(n)):
        bad = [i for i in range(n - 1) if not (rows[p[i]] >> p[i + 1] & 1)]
        if not bad:
            h += 1
            start[p[0]] += 1
            end[p[-1]] += 1
            for a, b in zip(p, p[1:]):
                q[a][b] += 1
        elif len(bad) == 1:
            i = bad[0]
            q[p[i]][p[i + 1]] += 1
    linear = [end[i] - start[i] + sum(q[i]) for i in range(n)]
    pair = {(i, j): q[i][j] + q[j][i]
            for i in range(n) for j in range(i + 1, n)}
    require(sum(linear) == sum(pair.values()), "constant-cut identity")
    return h, linear, pair


def evaluate(coeff, mask):
    h, linear, pair = coeff
    return h + sum(v for i, v in enumerate(linear) if mask >> i & 1) - sum(
        v for (i, j), v in pair.items() if mask >> i & 1 and mask >> j & 1)


def cube_by_new_permutations(rows):
    """Independent path count: enumerate NEW permutations, then legal cut masks."""
    n = len(rows)
    counts = [0] * (1 << n)
    for p in itertools.permutations(range(n + 1)):
        valid, ones, zeroes = True, 0, 0
        for a, b in zip(p, p[1:]):
            if b == n:
                ones |= 1 << a
            elif a == n:
                zeroes |= 1 << b
            elif not rows[a] >> b & 1:
                valid = False
                break
        if not valid:
            continue
        free = ((1 << n) - 1) ^ (ones | zeroes)
        sub = free
        while True:
            counts[ones | sub] += 1
            if not sub:
                break
            sub = (sub - 1) & free
    return counts


def rows_from_code(n, code):
    rows = [0] * n
    for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        a, b = (i, j) if code >> bit & 1 else (j, i)
        rows[a] |= 1 << b
    return tuple(rows)


def recover_from_probes(n, probes, omitted=None):
    h = probes[0]
    linear = [probes[1 << i] - h for i in range(n)]
    pair = {}
    for i, j in itertools.combinations(range(n), 2):
        if (i, j) != omitted:
            pair[i, j] = h + linear[i] + linear[j] - probes[(1 << i) | (1 << j)]
    if omitted is not None:
        pair[omitted] = sum(linear) - sum(pair.values())
    return h, linear, pair


def square_bank(labels):
    """One common scaled seed; independent candidates for EVERY size<=2 mask."""
    n, top = len(labels), max(labels)
    dummy = top + 1
    lam = 1
    while lam < 32 * (dummy + 1):
        lam *= 4
    scaled = [4 * lam * z for z in labels]
    candidates = {0: max(scaled) + 2}
    for size in (1, 2):
        for ids in itertools.combinations(range(n), size):
            a, b = sorted([labels[i] for i in ids] + ([dummy] if size == 1 else []))
            candidates[sum(1 << i for i in ids)] = (lam * (b - a) - 1) ** 2 - 4 * lam * a
    require(len(set(candidates.values())) == len(candidates), "candidate collision")
    for mask, x in candidates.items():
        require(x > max(scaled), "new label not last")
        require(sum(1 << i for i, a in enumerate(scaled) if square(a + x)) == mask,
                "unintended probe incidence")
    require(all(not square(x + y) for x, y in itertools.combinations(candidates.values(), 2)),
            "probe bank is not independent")
    require(all(square(labels[i] + labels[j]) == square(scaled[i] + scaled[j])
                for i, j in itertools.combinations(range(n), 2)), "old graph changed")
    return scaled, candidates, 4 * lam


def append_label(labels, incoming):
    if not labels:
        return [1], 1
    indices = [i for i in range(len(labels)) if incoming >> i & 1]
    require(len(indices) <= 2, "outside declared two-neighbor grammar")
    m = max(labels)
    if len(indices) < 2:
        x = (m + 1) ** 2 - (labels[indices[0]] if indices else 0)
        return labels + [x], 1
    old, bank, scale = square_bank(labels)
    return old + [bank[incoming]], scale


def decode(labels):
    """Gauge: earlier name -> later name exactly for a square sum."""
    rows = [0] * len(labels)
    for i, j in itertools.combinations(range(len(labels)), 2):
        a, b = (i, j) if square(labels[i] + labels[j]) else (j, i)
        rows[a] |= 1 << b
    return tuple(rows)


def compile_source(rows):
    labels, clock = [], []
    for j in range(len(rows)):
        incoming = sum(1 << i for i in range(j) if rows[i] >> j & 1)
        labels, scale = append_label(labels, incoming)
        clock.append(scale)
        require(decode(labels) == tuple(r & ((1 << (j + 1)) - 1) for r in rows[:j + 1]),
                "source compiler changed a named arc")
    return labels, clock


def rank_mod(rows, p=1000003):
    basis = {}
    selected = []
    for index, row in enumerate(rows):
        row = [v % p for v in row]
        for k in sorted(basis):
            c = row[k]
            if c:
                row = [(x - c * y) % p for x, y in zip(row, basis[k])]
        k = next((i for i, v in enumerate(row) if v), None)
        if k is not None:
            inv = pow(row[k], -1, p)
            basis[k] = [v * inv % p for v in row]
            selected.append(index)
    return len(basis), selected


def root_bracket(n, denominator=10**8):
    """Rational enclosure of n^(1/n), certified only by integer powers."""
    lo, hi = denominator, 2 * denominator
    target = n * denominator**n
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**n <= target:
            lo = mid
        else:
            hi = mid
    require(lo**n <= target < hi**n, "root interval")
    return Fraction(lo, denominator), Fraction(hi, denominator)


def main():
    result = {"status": "FINITE-EXACT controls; all-size claims require README proofs"}
    require(all(1000003 % d for d in range(2, math.isqrt(1000003) + 1)), "rank modulus prime")
    census = []
    for n in range(1, 6):
        vector_rows, checked = [], 0
        for code in range(1 << (n * (n - 1) // 2)):
            rows = rows_from_code(n, code)
            coeff = cut_coefficients(rows)
            direct = cube_by_new_permutations(rows)
            require([evaluate(coeff, s) for s in range(1 << n)] == direct,
                    ("boundary vs new permutations", n, code))
            probes = {s: direct[s] for s in range(1 << n) if s.bit_count() <= 2}
            omitted = (0, 1) if n >= 2 else None
            if omitted:
                del probes[3]
            recovered = recover_from_probes(n, probes, omitted)
            require(recovered == coeff, ("omitted probe recovery", n, code))
            require(all(v % 2 == 1 for v in direct), "raw count parity")
            vector_rows.append([coeff[0]] + coeff[1] + list(coeff[2].values()))
            checked += len(direct)
        rank, selected = rank_mod(vector_rows)
        census.append({"order": n, "sources": len(vector_rows), "cuts": checked,
                       "linear_rank_mod_1000003": rank, "basis_source_codes": selected})
    result["exhaustive_cut_census"] = census

    square_cases = 0
    for n in range(1, 6):
        for seed in itertools.combinations(range(1, 11), n):
            square_bank(list(seed))
            square_cases += 1
    result["square_seed_universe"] = {"subsets_of": [1, 10], "sizes": [1, 5], "cases": square_cases}

    b = (42, 4, 1, 22, 7, 94, 31)
    labels, clocks = compile_source(b)
    old, bank, scale = square_bank(labels)
    require(decode(old) == b, "common source drift")
    witnesses = []
    for target, masks, ports in [(93, [32, 128], [8, 3]), (95, [40], [5, 7]), (97, [8], [7, 3])]:
        rows, out = b, old[:]
        for j, mask in enumerate(masks):
            rows = birth(rows, mask)
            if j == 0:
                out.append(bank[mask])
            else:
                out, later_scale = append_label(out, mask)
                require(later_scale == 1, "unexpected later rescaling")
        require(out[:7] == old, "numeric shared ports changed")
        require(decode(out) == rows, "square decoder lost arcs")
        count = h_dp(rows)
        require(count == target, "incorrect branch target")
        witnesses.append({"target": target, "incoming_masks": masks, "ports": ports,
                          "rows": rows, "labels": [str(x) for x in out]})
    values = {str(s): h_dp(birth(b, s)) for s in (0, 32, 8, 40)}
    require(list(values.values()) == [31, 53, 97, 95], "four-port response")
    result["H31_shared_source"] = {
        "rows": b, "initial_labels": [str(x) for x in labels], "compiler_scales": clocks,
        "common_scale": str(scale), "shared_labels": [str(x) for x in old],
        "pair_response": values, "mixed_difference": 31 - 53 - 97 + 95,
        "witnesses": witnesses,
        "one_birth_all_cut_image": sorted(set(cube_by_new_permutations(b))),
    }
    # Entire one-birth response is still not a recurrence state.
    t, u = (2, 4, 1), (4, 1, 2)
    require(cube_by_new_permutations(t) == cube_by_new_permutations(u), "triangle probe twin")
    two_step = [h_dp(birth(birth(r, 1), 2)) for r in (t, u)]
    require(two_step == [13, 11], "two-step orientation separator")
    result["one_step_not_closed"] = {"sources": [t, u], "common_cube": cube_by_new_permutations(t),
                                     "incoming_word": [1, 2], "two_step_counts": two_step}
    # The exact MRB-side annihilator: each four-corner probe kills constants and singles.
    signs = [1, -1, -1, 1]
    require(sum(signs) == 0 and signs[1] + signs[3] == 0 and signs[2] + signs[3] == 0,
            "signed filter moments")
    result["finite_signed_filter"] = {"weights": signs, "B_response": list(values.values()),
                                      "interaction": -24, "MRB_closed_form_claimed": False}
    # Period-four MRB filters alpha=(-,+,-,+), beta=(-,-,+,+).
    # beta-alpha = 2 * sum_j [f(4j+3)-f(4j+2)]; the tail after j=2 is negative.
    bounds = {n: root_bracket(n) for n in (2, 3, 6, 7, 10, 11)}
    upper = 2 * sum(bounds[4*j+3][1] - bounds[4*j+2][0] for j in range(3))
    require(upper < Fraction(-29, 1000), "MRB signed orientation separation")
    result["MRB_orientation_separator"] = {
        "root_intervals": {str(n): [str(x) for x in v] for n, v in bounds.items()},
        "beta_minus_alpha_upper_bound": str(upper),
        "source_swap_difference_upper_bound": str(22 * upper),
        "four_cube_walsh_coefficients": [69, 5, 27, -6],
        "scope": "Exact source-profile functional; no new evaluation or arithmetic classification of MRB"}
    # Fixed anchors 1 and 3 admit the third anchor 2? No integer common neighbor to all three.
    require(not any(all(square(x + s * s * a) for a in (1, 2, 3))
                    for s in range(1, 13) for x in range(1, 10001)), "triple hostile")
    result["triple_anchor_probe"] = {"anchors": [1, 2, 3], "square_scales_s": [1, 12],
                                     "new_label_x": [1, 10000], "hits": 0,
                                     "universal_exclusion": "inherited Fermat area-one theorem"}
    result["script_sha256"] = hashlib.sha256(open(__file__, "rb").read()).hexdigest()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
