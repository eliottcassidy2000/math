"""Independent physical audit of the positive compiler and selector drift.

Uses this lane's own path counter, no producer/compiler/selector imports.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

from check_five_cores import canonical, check, from_bits, paths_dp


def g(n):
    total = 1
    while n > 1:
        q, bit = divmod(n, 2)
        if not bit:
            total += q + (q & 1)
        n = q
    return total


def construct(rows, out, a, b):
    left = [int(bit) for bit in bin(a)[2:]] if a else []
    right = [1-int(bit) for bit in reversed(bin(b)[2:])] if b else []
    ell, n = len(left), len(rows)
    cut = left + [(out >> i) & 1 for i in range(n)] + right
    apex = len(cut)
    result = [0] * (apex + 1)
    for u, v in combinations(range(apex), 2):
        if ell <= u < v < ell+n:
            winner, loser = (u, v) if rows[u-ell] >> (v-ell) & 1 else (v, u)
        else:
            winner, loser = u, v
        result[winner] |= 1 << loser
    for v, bit in enumerate(cut):
        winner, loser = (apex, v) if bit else (v, apex)
        result[winner] |= 1 << loser
    return tuple(result)


def path_data(rows):
    n = len(rows)
    full = (1 << n) - 1
    hp = {}
    for mask in range(full + 1):
        vertices = [v for v in range(n) if mask >> v & 1]
        hp[mask] = [p for p in permutations(vertices)
                    if all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))]
    starts = [sum(p[0] == v for p in hp[full]) for v in range(n)]
    ends = [sum(p[-1] == v for p in hp[full]) for v in range(n)]
    positive_start, positive_end = [0]*n, [0]*n
    signed_start, signed_end = [0]*n, [0]*n
    for mask in range(1, full+1):
        weight = len(hp[full ^ mask])
        sign = (-1)**mask.bit_count()
        for p in hp[mask]:
            positive_start[p[0]] += weight
            positive_end[p[-1]] += weight
            signed_start[p[0]] += sign*weight
            signed_end[p[-1]] += sign*weight
    check(signed_start == [-x for x in ends], "alternating start/end identity")
    check(signed_end == [-x for x in starts], "alternating end/start identity")
    check(sum(positive_start) == sum(positive_end), "positive path sum inconsistent")
    return len(hp[full]), starts, ends, positive_start, positive_end


def coefficients(rows, out, data):
    d, starts, ends, pstart, pend = data
    h = paths_dp(construct(rows, out, 0, 0))
    inside = [v for v in range(len(rows)) if not (out >> v & 1)]
    outside = [v for v in range(len(rows)) if out >> v & 1]
    return (h, d, d+sum(pend), sum(pend[v] for v in inside),
            sum(pstart[v] for v in outside), sum(starts[v] for v in inside),
            sum(ends[v] for v in outside))


def formula(z, a, b):
    h, d, ell, i, o, s, e = z
    return h+d*(g(a)+g(b)-2)+ell*a*b+i*a+o*b+s*(a&1)+e*(b&1)


def recover(values):
    h, p10, p30, p20, p01, p03, p11 = map(Fraction, values)
    i, o = (p30-p10)/2, (p03-p01)/2
    s, e = p10-h-i, p01-h-o
    return (h, (p20-h-2*i)/2, p11-h-i-o-s-e, i, o, s, e)


def rational_rank(rows):
    m = [list(map(Fraction, row)) for row in rows]
    rank = 0
    for column in range(7):
        pivot = next((i for i in range(rank, len(m)) if m[i][column]), None)
        if pivot is None:
            continue
        m[rank], m[pivot] = m[pivot], m[rank]
        leading = m[rank][column]
        m[rank] = [x/leading for x in m[rank]]
        for i in range(rank+1, len(m)):
            amount = m[i][column]
            if amount:
                m[i] = [x-amount*y for x, y in zip(m[i], m[rank])]
        rank += 1
    return rank


def main():
    contexts = [(0,0),(1,0),(3,0),(2,0),(0,1),(0,3),(1,1)]
    probes = list(dict.fromkeys(contexts + list(product(range(3), repeat=2)) + [(4,5),(5,4)]))
    classes = []
    for n in range(1, 6):
        classes.extend(sorted({canonical(from_bits(n, code))
                              for code in range(1 << (n*(n-1)//2))}))
    vectors, physical_checks = [], 0
    for rows in classes:
        data = path_data(rows)
        for out in range(1 << len(rows)):
            z = coefficients(rows, out, data)
            check(all(v >= 0 for v in z) and z[2] >= 2*z[1], "positivity failed")
            check(z[2] % 2 == 0 and (z[3]+z[5]) % 2 == 0 and (z[4]+z[6]) % 2 == 0,
                  "parity constraints failed")
            values = {pair: paths_dp(construct(rows, out, *pair)) for pair in probes}
            check(recover([values[pair] for pair in contexts]) == z, "seven contexts failed")
            for pair, value in values.items():
                check(formula(z, *pair) == value, "positive formula failed")
                physical_checks += 1
            vectors.append(z)
    check(len(classes) == 20 and len(vectors) == 470, "class/cut universe mismatch")
    check(rational_rank(vectors) == 7, "actual coefficient span is not rank seven")
    words = ["".join(bits) for k in range(1, 9) for bits in product("01", repeat=k)]
    drift_checks = 0
    for word in words:
        base = 1 << len(word)
        drifts = []
        for prefix, suffix in product(("", "1", "101"), ("", "0", "11")):
            values = [g(int(prefix+word*m+suffix or "0", 2)) for m in range(1, 8)]
            second = [values[i+2]-(1+base)*values[i+1]+base*values[i]
                      for i in range(len(values)-2)]
            drift = Fraction(second[0], 1-base)
            check(all(v == second[0] for v in second), "drift annihilator failed")
            check((drift > 0) == ("0" in word and "1" in word), "mixed drift sign failed")
            drifts.append(drift)
            drift_checks += 1
        check(len(set(drifts)) == 1, "fixed prefix/suffix changes linear drift")
    here = Path(__file__).resolve().parent
    report = {"status": "PASS; independent-agent audit controls, not an infinite selector census",
              "core_classes": len(classes), "cuts": len(vectors),
              "physical_child_DP_checks": physical_checks, "largest_child_order": 12,
              "coefficient_rank": 7, "loop_words_length_1_through_8": len(words),
              "drift_prefix_suffix_checks": drift_checks,
              "m_zero_hostile": {"g_1_actual": g(1), "formal_2_to_m_plus_1_at_zero": 2,
                                 "scope": "periodic formulas need m>=1; pumping proof already does"},
              "imports": "Only this lane's independent path/canonicalization code; no producer imports",
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (here / "compiler_selector_audit.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
