#!/usr/bin/env python3
"""Literal small certificates for the prime-kernel native switch and rebase."""
import argparse
import functools
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path
from check import require, hp, cycles

QMINUS = [22, 28, 8, 17, 4]
QPLUS = [22, 28, 8, 1, 12]
C3 = [2, 4, 1]


@functools.lru_cache(None)
def g(a):
    if a < 2:
        return 1
    q = a // 2
    return g(q) if a & 1 else g(q) + q + (q & 1)


def rev(a, bits):
    return int(f'{a:0{bits}b}'[::-1], 2) if bits else 0


def construct(core, cut, a, b):
    left = list(map(int, bin(a)[2:])) if a else []
    right = [1 - int(x) for x in bin(b)[2:][::-1]] if b else []
    l, n, r = len(left), len(core), len(right)
    root = l + n + r
    rows = [0] * (root + 1)
    for i, j in itertools.combinations(range(root), 2):
        forward = bool(core[i - l] >> (j - l) & 1) if l <= i < j < l + n else True
        x, y = (i, j) if forward else (j, i)
        rows[x] |= 1 << y
    bits = left + [cut >> v & 1 for v in range(n)] + right
    for v, bit in enumerate(bits):
        x, y = (root, v) if bit else (v, root)
        rows[x] |= 1 << y
    return rows


def dp(rows, start=None):
    n = len(rows)
    arr = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        if start is None or v == start:
            arr[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(arr[mask]):
            avail = rows[v] & ~mask
            while count and avail:
                bit = avail & -avail
                avail -= bit
                arr[mask | bit][bit.bit_length() - 1] += count
    return arr


def state(rows):
    root = len(rows) - 1
    full = (1 << len(rows)) - 1
    arr = dp(rows)
    h = [sum(x) for x in arr]
    h[0] = 1
    e = arr[-1][root]
    fixed = dp(rows, root)
    u = sum(sum(fixed[mask]) * h[full ^ mask] for mask in range(1, full + 1))
    delete = full ^ (1 << root)
    k = sum(h[mask] * h[delete ^ mask] for mask in range(delete + 1))
    return [h[-1], e, u, h[delete], k]


def permutation_state(rows):
    """Independent exhaustive directed paths, including short paths."""
    n = len(rows)
    full = (1 << n) - 1
    root = n - 1
    totals = [0] * (1 << n)
    rooted = [0] * (1 << n)
    totals[0] = 1
    ended = 0
    for length in range(1, n + 1):
        for word in itertools.permutations(range(n), length):
            if all(rows[a] >> b & 1 for a, b in zip(word, word[1:])):
                mask = sum(1 << v for v in word)
                totals[mask] += 1
                rooted[mask] += word[0] == root
                ended += length == n and word[-1] == root
    delete = full ^ (1 << root)
    return [totals[-1], ended,
            sum(rooted[m] * totals[full ^ m] for m in range(1 << n)),
            totals[delete], sum(totals[m] * totals[delete ^ m] for m in range(delete + 1))]


def signature(core, cut):
    full = (1 << len(core)) - 1
    counts = [sum(x) for x in dp(core)]
    counts[0] = 1
    sums = [0] * 8
    for length in range(1, len(core) + 1):
        for word in itertools.permutations(range(len(core)), length):
            if not all(core[a] >> b & 1 for a, b in zip(word, word[1:])):
                continue
            mask = sum(1 << v for v in word)
            w = counts[full ^ mask]
            starts = bool(cut >> word[0] & 1)
            ends = not bool(cut >> word[-1] & 1)
            for cat, yes in enumerate([True, ends, starts, starts and ends]):
                if yes:
                    sums[2 * cat] += w
                    sums[2 * cat + 1] += (-1) ** length * w
    return [counts[-1]] + sums


def formula(a, b, upper):
    lower = (9 * (g(a) + g(b) - 1) + 84 * a * b + 22 * a + 4 * (a % 2)
             + 41 * b + 7 * (b % 2) + 14)
    return lower + (3 * a + a % 2 + 2 if upper else 0)


def strong(rows):
    full = (1 << len(rows)) - 1
    for source in range(len(rows)):
        seen = 1 << source
        while True:
            new = seen
            for v in range(len(rows)):
                if seen >> v & 1:
                    new |= rows[v]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def n9(b):
    return 9 * g(b) + 41 * b + 7 * (b % 2) + 14


def n3(c):
    return 3 * g(c) + 12 * c + 10


def run():
    fixed = [(0, 1), (0, 2), (2, 3), (3, 0), (1, 2), (1, 4), (4, 2)]
    free = [(0, 4), (1, 3), (3, 4)]
    cube = []
    for code in range(8):
        rows = [0] * 5
        for i, j in fixed + [p if code >> k & 1 else p[::-1] for k, p in enumerate(free)]:
            rows[i] |= 1 << j
        require(strong(rows), "cube strongness")
        triangles = [c for c in cycles({i: {j for j in range(5) if rows[i] >> j & 1}
                                       for i in range(5)}, True) if len(c) == 3]
        arcs = [set(zip(c, c[1:] + c[:1])) for c in triangles]
        overlap = [[i, j] for i, j in itertools.combinations(range(len(arcs)), 2) if arcs[i] & arcs[j]]
        cube.append({"code": code, "rows": rows, "H": hp(rows), "triangles": triangles,
                     "triangle_arc_overlap_edges": overlap})
    require([x['H'] for x in cube] == [13, 11, 9, 9, 15, 9, 13, 9], "cube HP")
    sigs = [signature(q, 21) for q in [QMINUS, QPLUS]]
    require(sigs == [[9, 75, -9, 22, -4, 41, -7, 12, 2],
                     [9, 75, -9, 25, -5, 41, -7, 14, 2]], "core signatures")
    seed = []
    for upper, q in enumerate([QMINUS, QPLUS]):
        rows = construct(q, 21, 0, 0)
        actual = state(rows)
        require(actual == permutation_state(rows), "independent seed state")
        require(actual == [23 + 2 * upper, 2, 50, 9, 84], "seed")
        seed.append({"rows": rows, "HEUdK": actual})
    grid_checks = 0
    for a, b, upper in itertools.product(range(6), range(6), range(2)):
        rows = construct([QMINUS, QPLUS][upper], 21, a, b)
        expected = [formula(a, b, upper), 9 if b else 2,
                    (84 * a + 50) * (1 << b.bit_length()) - 9 * rev(b, b.bit_length()),
                    9, 84 * (1 << (a.bit_length() + b.bit_length()))]
        require(state(rows) == expected, "flank formula")
        require(strong(rows), "flank strongness")
        grid_checks += 1
    suffix_checks = 0
    for bits in range(1, 7):
        M = 1 << bits
        for s in range(M):
            for x in range(16):
                require(g(M * x + s) == g(x) + (M - 1 - rev(s, bits)) * x
                        + int(s < M // 2) * (x % 2) + g(s) - 1, "suffix recurrence")
                suffix_checks += 1
    defect_checks = 0
    for b, c in itertools.product(range(0, 32, 2), repeat=2):
        bp, cp = 4096 * b + 1824, 4096 * c + 3162
        L = 218 * b - 61 * c + 50
        D = n9(b) - n3(c)
        require(218 * bp - 61 * cp + 50 == 4096 * L, "L update")
        require(n9(bp) - n3(cp) == D + 936 * L, "D update")
        defect_checks += 1
    rebase = []
    b, c = 16, 58
    for j in range(8):
        require(n9(b) == n3(c), "paired target")
        h = n9(b)
        require(Fraction(39149312 * 4096 ** j + 169513, 47775) + Fraction(234 * j, 35) == h,
                "closed form")
        U9 = 50 * (1 << b.bit_length()) - 9 * rev(b, b.bit_length())
        U3 = 15 * (1 << c.bit_length()) - 3 * rev(c, c.bit_length())
        delta = 640 * 4096 ** j + 60 + Fraction(254592 * (4096 ** j - 1), 4095)
        require(U9 - U3 == delta and delta >= 700, "positive response reduction")
        rebase.append({"j": j, "B": b, "C": c, "lower_target": h, "upper_target": h + 2,
                       "H9_common_EUd": [9, U9, 9], "C3_common_EUd": [3, U3, 3]})
        b, c = 4096 * b + 1824, 4096 * c + 3162
    actual_rebase = []
    for name, q, cut, right in [('H9−', QMINUS, 21, 16), ('H9+', QPLUS, 21, 16),
                               ('C3−', C3, 3, 117), ('C3+', C3, 1, 235)]:
        rows = construct(q, cut, 0, right)
        actual = state(rows)
        expected_h = 825 if '+' in name else 823
        require(actual[0] == expected_h, "literal target rebase")
        require(strong(rows), "literal rebase strongness")
        actual_rebase.append({"name": name, "rows": rows, "HEUdK": actual})
    widths = []
    for k in range(1, 9):
        a = 2 * k - 1
        lower = 44 * k + 9 * g(k - 1) - 4
        require(formula(a, 0, False) == lower and formula(a, 0, True) == lower + 6 * k,
                "native multiple-six width")
        widths.append({"k": k, "A": a, "lower": lower, "upper": lower + 6 * k,
                       "EUd": [2, 168 * k - 34, 9]})
    return {"status": "FINITE-EXACT controls supporting unbounded proved formulas",
            "universe": {"completion_cube": 8, "flank_grid_A_B": "0..5, both cores",
                         "suffix_bits": "1..6, all suffixes, prefixes0..15",
                         "defect_inputs": "all even B,C in0..30", "rebase_repetitions": "0..7"},
            "completion_cube": cube, "core_signatures": sigs, "unit_seed": seed,
            "flank_grid_checks": grid_checks, "suffix_checks": suffix_checks,
            "defect_law_checks": defect_checks, "native_multiple_six_pairs": widths,
            "paired_rebase": rebase, "literal_823_825_witnesses": actual_rebase,
            "all_odd_target_coverage_claim": False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    output = json.dumps(run(), indent=2, sort_keys=True) + '\n'
    target = Path(__file__).with_name('switch_receipt.json')
    if args.write:
        target.write_text(output)
    else:
        require(target.read_text() == output, "switch receipt drift")
    print(json.dumps({"ok": True, "switch_receipt_sha256": hashlib.sha256(output.encode()).hexdigest(),
                      "unit_switch": [23, 25], "native_width_six": [49, 55],
                      "paired_rebase_seed": [823, 825]}, sort_keys=True))
