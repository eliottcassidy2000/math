#!/usr/bin/env python3
"""Exact controls for Boolean-affine multiplier-equivariant necklace maps.
No external dependencies. Assertions remain active under python -O.
"""
import argparse, hashlib, json, math
from pathlib import Path

checks = 0

def check(ok, *context):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(context)

def rot(w, m):
    return ((w << 1) & ((1 << m) - 1)) | (w >> (m - 1))

def anti(w, m):
    return ((w << 1) & ((1 << m) - 1)) | (1 ^ (w >> (m - 1)))

def derivative(w, m):
    return w ^ anti(w, m)

def canonical(w, m, twisted=False):
    op = anti if twisted else rot
    result = w
    for _ in range((2 if twisted else 1) * m - 1):
        w = op(w, m)
        result = min(result, w)
    return result

def period(w, m):
    z = rot(w, m)
    d = 1
    while z != w:
        z = rot(z, m)
        d += 1
    return d

def good(w, m):
    d = period(w, m)
    return ((w & ((1 << d) - 1)).bit_count() & 1) == 0

def pdiv(a, b):
    q = 0
    while a and a.bit_length() >= b.bit_length():
        k = a.bit_length() - b.bit_length()
        a ^= b << k
        q ^= 1 << k
    return q, a

def pgcd(a, b):
    while b:
        a, b = b, pdiv(a, b)[1]
    return a

def val1(a, cap):
    v = 0
    while v < cap:
        q, r = pdiv(a, 3)
        if r:
            break
        a = q
        v += 1
    return v

def products(h, m):
    columns = [h]
    for _ in range(1, m):
        columns.append(rot(columns[-1], m))
    table = [0] * (1 << m)
    for u in range(1, 1 << m):
        bit = u & -u
        table[u] = table[u ^ bit] ^ columns[bit.bit_length() - 1]
    return table

def totient(n):
    return sum(math.gcd(a, n) == 1 for a in range(1, n + 1))

def shell_orbits_formula(m, k):
    total = 0
    for d in range(1, m + 1):
        if m % d:
            continue
        E = d & -d
        if k <= E:
            total += totient(m // d) * (1 << (d - E)) * (1 if k == 0 else 1 << (k - 1))
    check(total % m == 0, "Burnside integrality", m, k)
    return total // m

def is_power_two(k):
    return k > 0 and k & (k - 1) == 0

def encode_affine(table, m):
    b = table[0]
    cols = tuple(table[1 << j] ^ b for j in range(m))
    return (b, cols)

def enumerate_raw_affine(m):
    """Independent full b,A enumeration, no circulant/derivative filtering."""
    n = 1 << m
    result = set()
    for packed in range(1 << (m * m)):
        cols = tuple((packed >> (m * j)) & (n - 1) for j in range(m))
        linear = [0] * n
        for f in range(1, n):
            bit = f & -f
            linear[f] = linear[f ^ bit] ^ cols[bit.bit_length() - 1]
        for b in range(n):
            table = [v ^ b for v in linear]
            if all(table[anti(f, m)] == rot(table[f], m) for f in range(n)):
                result.add((b, cols))
    return result

def run(max_m):
    rows = []
    for m in range(1, max_m + 1):
        n = 1 << m
        e = m & -m
        M = n | 1
        parity_odd = [u for u in range(n) if u.bit_count() & 1]
        B = {w for w in range(n) if good(w, m)}
        source_reps = [f for f in range(n) if canonical(f, m, True) == f]
        target_reps = {canonical(w, m) for w in B}
        check(len(source_reps) == len(target_reps), "type equality", m)
        for w in range(n):
            k = e - val1(w, e)
            d = period(w, m)
            check(good(w, m) == (k < (d & -d)), "height sign criterion", m, w, k, d)
        affine_encodings = set()
        safe = surjective = type_bijective = zero_safe = 0
        max_covered = max_covered_zero = 0
        max_rank = -1
        maps_by_height = {}
        for h in range(n):
            table = products(h, m)
            outputs = {table[u] for u in parity_odd}
            a = val1(h, e)
            k = e - a
            admissible = outputs <= B
            predicted = not is_power_two(k)
            check(admissible == predicted, "admissibility", m, h, k)
            check((0 in outputs) == (k == 0), "zero image", m, h)
            source_table = [table[derivative(f, m)] for f in range(n)]
            affine_encodings.add(encode_affine(source_table, m))
            check(all(source_table[anti(f, m)] == rot(source_table[f], m)
                      for f in range(n)), "equivariance", m, h)
            if admissible:
                safe += 1
                maps_by_height[k] = maps_by_height.get(k, 0) + 1
                rank = len(outputs).bit_length() - 1
                check(len(outputs) == 1 << rank, "affine image", m, h)
                max_rank = max(max_rank, rank)
                if 0 in outputs:
                    zero_safe += 1
                    check(all(val1(w, e) == e for w in outputs), "ideal socle", m, h)
                onto = outputs == B
                covered = len({canonical(source_table[f], m) for f in source_reps})
                max_covered = max(max_covered, covered)
                if 0 in outputs:
                    max_covered_zero = max(max_covered_zero, covered)
                bijective = covered == len(source_reps)
                expected = (pgcd(h, M) == 3 if m & 1 else m == 2 and h == 0)
                check(onto == expected, "surjective classification", m, h)
                check(bijective == expected, "type classification", m, h)
                surjective += onto
                type_bijective += bijective
                if onto and m & 1:
                    counts = {w: 0 for w in B}
                    for out in source_table:
                        counts[out] += 1
                    check(set(counts.values()) == {2}, "two labelled fibers", m, h)
        check(len(affine_encodings) == n, "unique normal forms", m)
        predicted_count = (1 << (m-e)) * ((1 << e) - sum(1 << ((1 << j)-1) for j in range(e.bit_length())))
        check(safe == predicted_count, "safe count", m)
        check(max_rank == (m-1 if m & 1 else m-2), "rank maximum", m)
        check(zero_safe == 1 << (m-e), "zero safe count", m)
        if m <= 3:
            raw = enumerate_raw_affine(m)
            check(raw == affine_encodings, "all arbitrary affine maps", m, len(raw))
        if m > 2 and not m & 1:
            check(3 in B and val1(3, e) < e, "adjacent pair hostile", m)
            # Disjoint adjacent pairs sum to the excluded all-one word.
            z = 0
            for j in range(0, m, 2):
                check((3 << j) in B, "disjoint pair target", m, j)
                z ^= 3 << j
            check(z == n-1 and z not in B, "target not affine", m)
        allowed_heights = [k for k in range(e + 1) if not is_power_two(k)]
        shell_counts = {k: shell_orbits_formula(m, k) for k in allowed_heights}
        check(max_covered_zero == shell_counts[0], "maximum zero coverage formula", m)
        check(max_covered == max(shell_counts.values()), "maximum shell coverage formula", m)
        def source_type_count(length):
            total = sum(totient(length // d) * (1 << (d - 1))
                        for d in range(1, length + 1)
                        if length % d == 0 and (length // d) % 2 == 1)
            check(total % length == 0, "odd-necklace Burnside integrality", length)
            return total // length
        A = source_type_count(m)
        check(A == len(source_reps), "source type divisor formula", m)
        optimum = A if e == 1 else ((A + source_type_count(m // 2)) // 2 if e == 2 else A // 2)
        check(max_covered == optimum, "closed sharp type coverage", m)

        for k in allowed_heights:
            exact_shell_orbits = {canonical(w, m) for w in range(n) if e-val1(w, e) == k}
            check(len(exact_shell_orbits) == shell_counts[k], "shell orbit formula", m, k)
        rows.append(dict(m=m, dyadic_depth=e, source_types=len(source_reps),
                         target_types=len(target_reps), target_labelled=len(B),
                         affine_maps=n, admissible_affine_maps=safe,
                         admissible_by_nilpotent_height=maps_by_height,
                         admissible_maps_hitting_zero=zero_safe,
                         maximum_affine_rank=max_rank,
                         maximum_target_types_covered=max_covered,
                         maximum_target_types_covered_hitting_zero=max_covered_zero,
                         affine_type_bijections=type_bijective,
                         affine_labelled_surjections=surjective))
        print(json.dumps(rows[-1], sort_keys=True), flush=True)
    return rows

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-m', type=int, default=10)
    parser.add_argument('--output', default='verification.json')
    args = parser.parse_args()
    rows = run(args.max_m)
    result = dict(status='PASS', universe=f'all m=1,...,{args.max_m}; every polynomial h and every input word; all arbitrary affine maps for m<=3',
                  checks=checks, source_map_evaluations=sum(4**m for m in range(1, args.max_m+1)),
                  odd_word_map_evaluations=sum(2**(2*m-1) for m in range(1, args.max_m+1)), rows=rows,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(f'PASS: {checks} checks')
