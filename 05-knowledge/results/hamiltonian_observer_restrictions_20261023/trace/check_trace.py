"""Exact small controls for guarded return chains. Requires sympy; no floats."""
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path
import json
import sympy as sp


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def psd(a):
    """Exact rational symmetric elimination; a zero pivot requires a zero row."""
    a = sp.Matrix(a)
    need(a == a.T, "PSD test requires symmetry")
    for k in range(a.rows):
        d = a[k, k]
        if d < 0:
            return False
        if d == 0:
            if any(a[k, j] != 0 for j in range(k + 1, a.rows)):
                return False
            continue
        for i in range(k + 1, a.rows):
            for j in range(i, a.rows):
                a[i, j] -= a[i, k] * a[k, j] / d
                a[j, i] = a[i, j]
    return True


def projection(keys):
    counts = Counter(keys)
    return sp.Matrix(len(keys), len(keys),
                     lambda i, j: sp.Rational(1, counts[keys[i]])
                     if keys[i] == keys[j] else 0)


def words(n):
    edges = list(combinations(range(n), 2))
    index = {e: k for k, e in enumerate(edges)}
    orders = list(permutations(range(n)))

    def encode(code, order, tournament):
        ans = 0
        for k, (i, j) in enumerate(edges):
            u, v = order[i], order[j]
            bit = (code >> index[tuple(sorted((u, v)))]) & 1
            if tournament and u > v:
                bit ^= 1
            ans |= bit << k
        return ans

    gkeys = [min(encode(x, p, False) for p in orders)
             for x in range(1 << len(edges))]
    tkeys = [min(encode(x, p, True) for p in orders)
             for x in range(1 << len(edges))]
    pg, pt = projection(gkeys), projection(tkeys)
    need(pg * pg == pg and pt * pt == pt, "projections")
    return pg * pt * pg, tkeys


def check_gap(r, pi, gamma):
    d = sp.diag(*pi)
    need(r * sp.ones(r.rows, 1) == sp.ones(r.rows, 1), "row sums")
    need(all(x >= 0 for x in r), "nonnegative kernel")
    need(d * r == r.T * d, "reversibility")
    need(psd(d * r), "PSD operator")
    p = sp.Matrix(pi)
    need(psd((1 - gamma) * d + gamma * p * p.T - d * r), "gap inequality")


def trace(q, pi, selected, gamma):
    n, s = q.rows, list(selected)
    c = [i for i in range(n) if i not in s]
    mass = sum(pi[i] for i in s)
    pis = [pi[i] / mass for i in s]
    qss = q.extract(s, s)
    if not c:
        r, clock = qss, sp.ones(len(s), 1)
    else:
        qsc, qcs, qcc = q.extract(s, c), q.extract(c, s), q.extract(c, c)
        b = (sp.eye(len(c)) - qcc).inv()
        r = qss + qsc * b * qcs
        extension = sp.zeros(n, len(s))
        for i, old in enumerate(s):
            extension[old, i] = 1
        for i, old in enumerate(c):
            extension[old, :] = (b * qcs)[i, :]
        energy = extension.T * sp.diag(*pi) * (sp.eye(n) - q) * extension
        need(energy == mass * sp.diag(*pis) * (sp.eye(len(s)) - r), "harmonic energy")
        clock = sp.ones(len(s), 1) + qsc * b * sp.ones(len(c), 1)
        derivative = qss + 2 * qsc * b * qcs + qsc * b * qcc * b * qcs
        need(derivative * sp.ones(len(s), 1) == clock, "clock derivative")
        # First four duration coefficients and the remaining mass.
        prefix = qss.copy()
        for ell in range(2, 5):
            coeff = qsc * qcc ** (ell - 2) * qcs
            need(all(x >= 0 for x in coeff), "duration coefficient")
            prefix += coeff
        need(all(x >= 0 for x in r - prefix), "duration tail mass")
    check_gap(r, pis, gamma)
    need((sp.Matrix(1, len(s), pis) * clock)[0] == 1 / mass, "Kac mean clock")
    return r


def all_nonempty(n):
    for bits in range(1, 1 << n):
        yield [i for i in range(n) if (bits >> i) & 1]


def owner_fibres(q):
    """Order-four actual source, observer, pin and converse controls.

    Q_lift = Q x J_24 need not be materialized as a 192x192 matrix.
    For acceptance counts k_y among the 24 owner orders, D_yy=1-k_y/24.
    (I-QD)^-1 Q gives the sum of all rejected-word excursions.
    """
    n = 4
    edges = list(combinations(range(n), 2))
    edge_index = {e: i for i, e in enumerate(edges)}
    free = [(i, j) for i in range(n) for j in range(i + 2, n)]
    aux = list(combinations(range(n - 1), 2))
    free.sort(key=lambda e: aux.index((e[0], e[1] - 1)))
    orders = list(permutations(range(n)))
    _, tkeys = words(n)
    reps = sorted(set(tkeys))

    def backwards(code, u, v):
        bit = (code >> edge_index[tuple(sorted((u, v)))]) & 1
        return bit ^ (u > v)

    def frame(code, p):
        if any(backwards(code, p[i], p[i + 1]) for i in range(n - 1)):
            return None
        return sum(int(backwards(code, p[i], p[j])) << k
                   for k, (i, j) in enumerate(free))

    records, checks = [], 0
    for code in reps:
        fibre = [(frame(code, p), p) for p in orders if frame(code, p) is not None]
        # Independently reconstruct the full literal tournament from each frame.
        for x, p in fibre:
            rebuilt = 0
            for i, j in combinations(range(n), 2):
                b = 0 if j == i + 1 else ((x >> free.index((i, j))) & 1)
                u, v = p[i], p[j]
                rebuilt |= int(b ^ (u > v)) << edge_index[tuple(sorted((u, v)))]
            need(rebuilt == code, "actual source reconstruction")
            conv = frame(code ^ 63, p[::-1])
            reflected = sum(((x >> free.index((n - 1 - j, n - 1 - i))) & 1) << k
                            for k, (i, j) in enumerate(free))
            need(conv == reflected, "converse with actual reversed order")
        for ids in all_nonempty(len(fibre)):
            accepted = [fibre[i] for i in ids]
            counts = Counter(x for x, p in accepted)
            d = sp.diag(*[1 - sp.Rational(counts[i], 24) for i in range(8)])
            b = (sp.eye(8) - q * d).inv()
            f = b * q
            r = sp.Matrix(len(ids), len(ids), lambda i, j: f[accepted[i][0], accepted[j][0]] / 24)
            check_gap(r, [sp.Rational(1, len(ids))] * len(ids), sp.Rational(8, 9))
            clocks = b * sp.ones(8, 1)
            mean_clock = sum(clocks[x] for x, p in accepted) / len(ids)
            need(mean_clock == sp.Rational(192, len(ids)), "owner return clock")
            checks += 1
        # Every source/pin subset is among the exhaustive subsets above.
        pins = {f"{i}:{v}": sum(p[i] == v for x, p in fibre)
                for i in range(n) for v in range(n)}
        records.append({"tournament_code": code, "Hamiltonian_paths": len(fibre),
                        "full_mean_clock": str(sp.Rational(192, len(fibre))),
                        "indexed_pin_counts": pins})
    return checks, records


def main():
    q, _ = words(3)
    gamma = sp.Rational(8, 9)
    check_gap(q, [sp.Rational(1, 8)] * 8, gamma)
    need((q - sp.eye(8) / 9).det() == 0, "sharp ambient gap")
    restrictions = 0
    spectra = Counter()
    for selected in all_nonempty(8):
        r = trace(q, [sp.Rational(1, 8)] * 8, selected, gamma)
        spectra[str(r.trace())] += 1
        restrictions += 1
    # Unequal weights detect accidental use of uniform stationary laws.
    pi = [sp.Rational(1, 2), sp.Rational(1, 3), sp.Rational(1, 6)]
    nonuniform = (sp.eye(3) + sp.ones(3, 1) * sp.Matrix(1, 3, pi)) / 2
    for selected in all_nonempty(3):
        trace(nonuniform, pi, selected, sp.Rational(1, 2))
    refresh = sp.ones(8, 8) / 8
    watched = trace(refresh, [sp.Rational(1, 8)] * 8, [0, 1], sp.Integer(1))
    rejected = refresh.extract([0, 1], [0, 1]) + sp.eye(2) * sp.Rational(3, 4)
    need(watched == sp.ones(2, 2) / 2, "refresh trace")
    need(rejected * sp.Matrix([1, -1]) == sp.Rational(3, 4) * sp.Matrix([1, -1]),
         "rejection gap is only 1/4")
    owner_checks, owner_records = owner_fibres(q)
    receipt = {"status": "FINITE-EXACT", "arithmetic": "SymPy rational matrices; exact PSD elimination",
               "word_universe": 8, "all_nonempty_word_subsets": restrictions,
               "nonuniform_subsets": 7, "owner_subset_checks": owner_checks,
               "owner_fibres": owner_records,
               "guards": "all nonempty finite subsets, not an invented native mod8 formula",
               "return_gap_lower_bound": "8/9", "rejection_control_gap": "1/4",
               "trace_matrix_trace_histogram": dict(sorted(spectra.items()))}
    path = Path(__file__).with_name("trace_receipt.json")
    path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ["status", "all_nonempty_word_subsets",
                     "nonuniform_subsets", "owner_subset_checks", "return_gap_lower_bound"]}))
    print("ALL EXACT TRACE, CLOCK, OWNER AND CONVERSE CHECKS PASSED")


if __name__ == "__main__":
    main()
