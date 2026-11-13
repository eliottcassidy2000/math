#!/usr/bin/env python3
"""Small exact controls for the HP/Selmer and F8 subgroup interfaces.

No imported research implementation. All load-bearing gates survive -O.
"""
import hashlib
import itertools as it
import json
from collections import defaultdict
from fractions import Fraction
from math import comb, factorial, isqrt, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def adjacency(n, code):
    a = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(it.combinations(range(n), 2)):
        a[i][j] = (code >> k) & 1
        a[j][i] = 1 - a[i][j]
    return a


def code_of(a):
    return sum(a[i][j] << k for k, (i, j) in enumerate(it.combinations(range(len(a)), 2)))


def rowzero(a):
    b = [r[:] for r in a]
    for i, r in enumerate(b):
        r[i] = sum(r) % 2
    return b


def phi(a):
    n = len(a)
    return [[(sum(a[i][k] * a[k][j] for k in range(n)) + a[i][j] + (i == j)) % 2
             for j in range(n)] for i in range(n)]


def rank(a):
    b = [r[:] for r in a]
    k = 0
    for j in range(len(b[0]) if b else 0):
        p = next((i for i in range(k, len(b)) if b[i][j]), None)
        if p is None:
            continue
        b[k], b[p] = b[p], b[k]
        for i in range(len(b)):
            if i != k and b[i][j]:
                b[i] = [x ^ y for x, y in zip(b[i], b[k])]
        k += 1
    return k


def kernel(a):
    n = len(a[0]) if a else 0
    rows = [sum(x << j for j, x in enumerate(r)) for r in a]
    return [v for v in range(1 << n) if all((v & r).bit_count() % 2 == 0 for r in rows)]


def apply(a, v):
    return sum((sum(x * ((v >> j) & 1) for j, x in enumerate(r)) % 2) << i for i, r in enumerate(a))


def nu(a):
    return len(a) - rank(phi(rowzero(a)))


def h_direct(a):
    return sum(all(a[p[i]][p[i + 1]] for i in range(len(a) - 1))
               for p in it.permutations(range(len(a))))


def two_slot_wedges(a):
    """Literal covers by two ordered optional slots, independently of local tests."""
    n = len(a)
    found = set()
    for p in it.permutations(range(n)):
        for cut in range(n + 1):
            if all(j == cut - 1 or a[p[j]][p[j + 1]] for j in range(n - 1)):
                for j in range(n - 2):
                    if cut not in (j + 1, j + 2):
                        u, v, w = p[j:j + 3]
                        found.add((v, min(u, w), max(u, w)))
    return found


def all_h(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for s in range(1, 1 << n):
        for j in range(n):
            if s >> j & 1:
                t = s ^ (1 << j)
                if t:
                    dp[s][j] = sum(dp[t][i] * a[i][j] for i in range(n))
    return [1] + [sum(dp[s]) for s in range(1, 1 << n)]


def covers_dp(a):
    n = len(a)
    hs = all_h(a)
    dp = [[0] * (1 << n) for _ in range(n + 1)]
    dp[0][0] = 1
    for j in range(1, n + 1):
        for s in range(1, 1 << n):
            t = s
            while t:
                dp[j][s] += hs[t] * dp[j - 1][s ^ t]
                t = (t - 1) & s
    return [dp[j][-1] for j in range(1, n + 1)]


def covers_words(a):
    n = len(a)
    counts = [0] * n
    for p in it.permutations(range(n)):
        bad = sum((not a[p[i]][p[i + 1]]) << i for i in range(n - 1))
        for cuts in range(1 << (n - 1)):
            if cuts & bad == bad:
                counts[cuts.bit_count()] += 1
    return counts


def q_coeff(fs):
    out = [Fraction(0)] * (len(fs) + 1)
    b = [Fraction(1)]
    for j, f in enumerate(fs, 1):
        c = [Fraction(0)] * (len(b) + 1)
        for k, x in enumerate(b):
            c[k] -= x * (j - 1) / j
            c[k + 1] += x / j
        b = c
        for k, x in enumerate(b):
            out[k] += f * x
    need(all(x.denominator == 1 for x in out), "integral Q controls")
    return [int(x) for x in out]


def join(a, b):
    na, nb = len(a), len(b)
    return [a[i] + [1] * nb for i in range(na)] + [[0] * na + r for r in b]


def convolution(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def walsh(v):
    a = v[:]
    h = 1
    while h < len(a):
        for i in range(0, len(a), 2 * h):
            for j in range(i, i + h):
                x, y = a[j], a[j + h]
                a[j], a[j + h] = x + y, x - y
        h *= 2
    return a


def isprime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def prime_model(ps):
    need(len(ps) == len(set(ps)) and all(isprime(p) and p % 8 == 3 for p in ps), "prime source")
    n = len(ps)
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                # Euler criterion and direct square-set membership independently.
                a[i][j] = int(pow(ps[j] % ps[i], (ps[i] - 1) // 2, ps[i]) == ps[i] - 1)
                squares = {x * x % ps[i] for x in range(1, ps[i])}
                need(a[i][j] == int(ps[j] % ps[i] not in squares), "Legendre two paths")
    m = rowzero(a)
    block = []
    for i in range(n):
        block.append([m[i][j] ^ (i == j) for j in range(n)] + [int(i == j) for j in range(n)])
    for i in range(n):
        block.append([int(i == j) for j in range(n)] + m[i])
    need(2 * n - rank(block) == nu(a), "full Monsky specialization")
    return {"primes": ps, "N": prod(ps), "code": code_of(a), "nu": nu(a), "Monsky_rank": rank(block)}


def tournament_controls():
    minima, tables = [], {}
    for n in range(1, 6):
        m = comb(n, 2)
        hs, groups, slot_observations = [], defaultdict(set), []
        for c in range(1 << m):
            a = adjacency(n, c)
            h = h_direct(a)
            need(h == all_h(a)[-1], "HP algorithms")
            d = nu(a)
            need(len(kernel(phi(rowzero(a)))) == 1 << d, "nullity algorithms")
            need(d % 2 == 0, "Phi3 even nullity")
            op = adjacency(n, c ^ ((1 << m) - 1))
            need(d == nu(op), "converse descent")
            need(phi(rowzero(op)) == phi(rowzero(a)), "full labelled Phi converse invariance")
            need(kernel(phi(rowzero(op))) == kernel(phi(rowzero(a))), "embedded Phi kernel converse invariance")
            ks = kernel(phi(rowzero(a)))
            need(all(v.bit_count() % 2 == 0 for v in ks), "Phi kernel even coordinate sum")
            need(all(apply(rowzero(op), v) ^ apply(rowzero(a), v) == ((n % 2) * v) for v in ks), "exact Monsky embedding parity change")
            observed = two_slot_wedges(a)
            local = {(v, u, w) for v in range(n)
                     for u, w in it.combinations([x for x in range(n) if x != v], 2)
                     if a[u][v] * a[v][w] + a[w][v] * a[v][u]}
            need(observed == local, "literal two-slot/local wedge equality")
            slot_observations.append(observed)
            groups[h].add(d)
            hs.append(h)
        collisions = {h: sorted(ds) for h, ds in groups.items() if len(ds) > 1}
        minima.append({"n": n, "labelled_tournaments": 1 << m, "H_collisions": collisions})
        if n >= 3:
            fw = walsh(hs)
            edges = list(it.combinations(range(n), 2))
            tree = []
            reached = {0}
            while len(reached) < m:
                e, f = next((e, f) for e in sorted(reached) for f in range(m)
                            if f not in reached and set(edges[e]) & set(edges[f]))
                tree.append((e, f))
                reached.add(f)
            for e, f in it.combinations(range(m), 2):
                common = set(edges[e]) & set(edges[f])
                if common:
                    b = next(iter(common))
                    aa = next(x for x in edges[e] if x != b)
                    cc = next(x for x in edges[f] if x != b)
                    epsilon = (1 if aa < b else -1) * (1 if b < cc else -1)
                    expected = epsilon * 2 * factorial(n - 2) * (1 << (m - n + 1))
                else:
                    expected = 0
                need(fw[(1 << e) | (1 << f)] == expected, "quadratic Walsh coefficient")
            for code in range(1 << m):
                signs = {0: 1}
                slot_signs = {0: 1}
                for e, f in tree:
                    signs[f] = signs[e] * (1 if ((code >> e) ^ (code >> f)) & 1 == 0 else -1)
                    v = next(iter(set(edges[e]) & set(edges[f])))
                    u = next(x for x in edges[e] if x != v)
                    w = next(x for x in edges[f] if x != v)
                    epsilon = (1 if u < v else -1) * (1 if v < w else -1)
                    beta = (v, min(u, w), max(u, w)) in slot_observations[code]
                    slot_signs[f] = slot_signs[e] * epsilon * (2 * beta - 1)
                recovered = sum((x > 0) << e for e, x in signs.items())
                need(recovered in (code, code ^ ((1 << m) - 1)), "response quotient decode")
                slot_recovered = sum((x > 0) << e for e, x in slot_signs.items())
                need(slot_recovered == recovered, "literal two-slot quotient decode")
                need(nu(adjacency(n, slot_recovered)) == nu(adjacency(n, code)), "actual wedge Selmer consumer")
                need(phi(rowzero(adjacency(n, slot_recovered))) == phi(rowzero(adjacency(n, code))), "actual wedge full Phi consumer")
                need(kernel(phi(rowzero(adjacency(n, slot_recovered)))) == kernel(phi(rowzero(adjacency(n, code)))), "actual wedge embedded kernel consumer")
        tables[n] = hs
    need(all(not row["H_collisions"] for row in minima[:4]), "first failure not below five")
    need(minima[-1]["H_collisions"] == {9: [0, 2], 15: [0, 4]}, "first scalar collisions")
    witnesses = []
    for c in (8, 10, 40, 76):
        a = adjacency(5, c)
        fs = covers_dp(a)
        need(fs == covers_words(a), "ordered covers two paths")
        witnesses.append({"code": c, "H": h_direct(a), "scores": sorted(map(sum, a)),
                          "ordered_covers": fs, "Q": q_coeff(fs), "M": rowzero(a),
                          "Phi_M": phi(rowzero(a)), "kernel": kernel(phi(rowzero(a))), "nu": nu(a)})
    x, y = witnesses[:2]
    need(x["Q"] == y["Q"] == [0, 2, 0, 6, 0, 1], "full Q collision")
    need(x["scores"] == y["scores"] == [1, 1, 2, 3, 3], "score collision")
    need((x["nu"], y["nu"]) == (2, 0), "arithmetic boundary")
    # All joins of labelled orders 1..3: independently verify both invariants.
    for na in range(1, 4):
        for nb in range(1, 4):
            for ca in range(1 << comb(na, 2)):
                for cb in range(1 << comb(nb, 2)):
                    a, b = adjacency(na, ca), adjacency(nb, cb)
                    c = join(a, b)
                    pc, pa, pb = phi(rowzero(c)), phi(rowzero(a)), phi(rowzero(b))
                    target = [r + [0] * nb for r in pa] + [[0] * na + r for r in pb]
                    need(pc == target, "literal Phi join direct sum")
                    need(q_coeff(covers_dp(c)) == convolution(q_coeff(covers_dp(a)), q_coeff(covers_dp(b))), "Q join product")
    families = []
    for n in range(5, 10):
        suffix = [[int(i < j) for j in range(n - 5)] for i in range(n - 5)]
        pair = [join(adjacency(5, c), suffix) for c in (8, 10)]
        need([nu(a) for a in pair] == [2, 0], "all-size family finite probe")
        need(sorted(map(sum, pair[0])) == sorted(map(sum, pair[1])), "family scores")
        qs = [q_coeff(covers_dp(a)) for a in pair]
        need(qs[0] == qs[1] == [0] * (n - 5) + [0, 2, 0, 6, 0, 1], "family Q")
        families.append({"n": n, "nu": [2, 0], "Q": qs[0]})
    arithmetic = [prime_model(ps) for ps in ([3, 19, 43, 139, 251], [3, 19, 83, 499, 1163])]
    need([x["code"] for x in arithmetic] == [8, 10], "named arithmetic witnesses")
    return {"minimality": minima, "witnesses": witnesses, "joined_controls": families, "arithmetic": arithmetic}


def mul(a, b):
    c = 0
    while b:
        if b & 1:
            c ^= a
        a <<= 1
        if a & 8:
            a ^= 11  # X^3+X+1
        b >>= 1
    return c


def power(a, n):
    x = 1
    for _ in range(n):
        x = mul(x, a)
    return x


def trace(a):
    return a ^ power(a, 2) ^ power(a, 4)


def fourier(s):
    u = v = 0
    for i in range(7):
        if s >> i & 1:
            u ^= power(2, i)
            v ^= power(2, (-i) % 7)
    return (u, v)


def act(w, b, k):
    u, v = w
    return mul(power(2, b), power(u, 1 << k)), mul(power(2, (-b) % 7), power(v, 1 << k))


def span(vs):
    s = {0}
    for v in vs:
        s |= {x ^ v for x in list(s)}
    return s


def f8_controls():
    even = [s for s in range(128) if s.bit_count() % 2 == 0]
    need(len({fourier(s) for s in even}) == 64, "augmentation Fourier isomorphism")
    for s in even:
        for b in range(7):
            for k in range(3):
                mapped = sum(1 << ((power(2, 0) * (1 << k) * i + b) % 7)
                             for i in range(7) if s >> i & 1)
                need(fourier(mapped) == act(fourier(s), b, k), "full G21 equivariance")
    z = [a for a in range(8) if trace(a) == 0]
    need(z == [0, 2, 4, 6], "trace-zero sector")
    need({power(2, i) for i in (1, 2, 4)} == set(z) - {0}, "one-two-four Fano line")
    need(all(trace(a) in (0, 1) for a in range(8)), "trace codomain")
    need(all(power(a, 8) == a for a in range(8)), "F8 control")
    need(all(power(a, 4) ^ power(a, 2) ^ a == 0 for a in z), "Phi3 restriction")
    need(trace(4) == 0 and trace(mul(2, 4)) == 1, "C7 destroys sector")
    need(mul(4, 6) == 5 and trace(5) == 1, "trace-zero not a subfield")
    # A cyclic C7 submodule generated by a mixed vector is all W; an axis gives that axis.
    module_dimensions = defaultdict(int)
    for u, v in it.product(range(8), repeat=2):
        orbit = [act((u, v), b, 0) for b in range(7)]
        s = span([a | (b << 3) for a, b in orbit])
        dim = (len(s)).bit_length() - 1
        need(dim == (0 if not (u or v) else 3 if not u or not v else 6), "C7 generated submodule")
        module_dimensions[dim] += 1
    # Two reciprocal cubics acting on the two axes.
    for a in range(8):
        need(mul(power(2, 3), a) ^ mul(2, a) ^ a == 0, "first cubic")
        beta = power(2, 6)
        need(mul(power(beta, 3), a) ^ mul(power(beta, 2), a) ^ a == 0, "reciprocal cubic")
    lines = {tuple(sorted(mul(power(2, b), x) for x in z if x)) for b in range(7)}
    need(len(lines) == 7, "seven Fano lines")
    need(all(sum(set(p) <= set(line) for line in lines) == 1 for p in it.combinations(range(1, 8), 2)), "Fano unique line")
    # All 30 labelled simplex codes: same pair marginals, third moments distinguish incidence.
    codes = {}
    for labels in it.permutations(range(1, 8)):
        words = tuple(sorted(sum(((a & v).bit_count() % 2) << i for i, v in enumerate(labels)) for a in range(8)))
        if words in codes:
            continue
        triples = []
        for i in range(7):
            need(sum((-1) ** ((w >> i) & 1) for w in words) == 0, "first Fano moment")
        for i, j in it.combinations(range(7), 2):
            counts = [sum((((w >> i) & 1) + 2 * ((w >> j) & 1)) == b for w in words) for b in range(4)]
            need(counts == [2] * 4, "pairwise Fano blindness")
        for t in it.combinations(range(7), 3):
            moment = sum((-1) ** sum((w >> i) & 1 for i in t) for w in words)
            expected = 8 if labels[t[0]] ^ labels[t[1]] ^ labels[t[2]] == 0 else 0
            need(moment == expected, "cubic Fano recovery")
            if moment:
                triples.append(t)
        codes[words] = triples
    need(len(codes) == len({tuple(x) for x in codes.values()}) == 30, "thirty faithful cubic tables")
    return {"polynomial": "X^3+X+1", "trace_zero": z, "trace_zero_exponents": [1, 2, 4],
            "Phi3_sigma_kernel_dimension": 4, "C7_submodule_dimensions": sorted(module_dimensions),
            "cyclic_generator_dimension_counts": dict(module_dimensions),
            "translation_hostile": {"input": [4, 0], "output": [3, 0], "traces": [0, 1]},
            "fano_lines": sorted(lines), "fano_codes": len(codes)}


def main():
    sources = [
        "05-knowledge/results/bsd_selmer_growth_20261013a.md",
        "05-knowledge/results/bsd_component_recursion_20261013c.md",
        "05-knowledge/results/bsd_tunnell_terminals_20261013c.md",
        "05-knowledge/reference/CORE-PAPERS-OPENAI-FIVE-MANUSCRIPTS-2026-10-13.md",
        "05-knowledge/results/hn_metagraph_information_20261026/information/RESULT.md",
        "05-knowledge/results/session320_second_moments_and_native_contexts_20261113/WEDGE_RECONSTRUCTION.md",
        "05-knowledge/results/session320_second_moments_and_native_contexts_20261113/F1_WEDGE_JET.md",
        "05-knowledge/results/session61_hadwiger_six_seven_20261025/palette-deck/RESULT.md",
        "05-knowledge/results/session54_reframed_numerology_20261025/sixtythree-bridge/README.md",
        "01-canon/theorems/THM-5125-residue-tournament-model-heegner-domination-and-cycle-gas-identities.md",
    ]
    result = {"status": "FINITE-EXACT controls; general candidate proof separate and RESERVED pending audit",
              "base_commit": "68d0d86bef480b6a36efc3d9864fd01ff3244bf5",
              "gauge": "code bit 1 is increasing arc in lexicographic pair order; M is row-zero binary adjacency",
              "universe": "all labelled tournaments orders1..5; all joins of orders1..3; specified family controls5..9; all64 augmentation states and21 group elements; all30 labelled Fano codes",
              "tournaments": tournament_controls(), "F8": f8_controls()}
    result["source_pins"] = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in sources}
    result["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result["checks"] = CHECKS
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
