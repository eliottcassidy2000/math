"""Exact controls for a uniform diagonal-density theorem; no census extrapolation.

Python 3 + SymPy.  Every gate remains active under -O.  The all-size proof is
in the companion report; this program checks its normalizations, local
expansion, integrals, and independent small labeled carriers.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial
import hashlib
import json
import sys
import sympy as sp

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(test, label):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(label)


def falling(n, k):
    ans = 1
    for j in range(k):
        ans *= n - j
    return ans


def skeleton(parts):
    edges, owner, off = [], {}, 0
    for block, size in enumerate(parts):
        for j in range(size):
            for col in (j, (j + 1) % size):
                edge = (off + j, off + col)
                edges.append(edge)
                owner[edge] = block
        off += size
    need(len(set(edges)) == 2 * off, "simple cycle carrier")
    return tuple(edges), owner


def parts_at_least_two(n, lo=2):
    if n == 0:
        yield ()
    for j in range(lo, n + 1):
        for tail in parts_at_least_two(n - j, j):
            yield (j,) + tail


def rook_product(parts, cap=6):
    ans = [1] + [0] * cap
    for size in parts:
        length = 2 * size
        factor = [Q(length, length - k) * comb(length - k, k)
                  if k <= size else 0 for k in range(cap + 1)]
        ans = [sum(ans[j] * factor[k - j] for j in range(k + 1))
               for k in range(cap + 1)]
    need(all(Q(x).denominator == 1 for x in ans), "integral rook coefficients")
    return list(map(int, ans))


def exact_matchings(edges, cap=6):
    ans = []
    for k in range(cap + 1):
        count = 0
        for selected in combinations(edges, k):
            if len({u for u, _ in selected}) == k and len({v for _, v in selected}) == k:
                count += 1
        ans.append(count)
    return ans


def full_collinear(board):
    return sum((b[0] - a[0]) * (c[1] - a[1]) ==
               (b[1] - a[1]) * (c[0] - a[0])
               for a, b, c in combinations(board, 3))


def graph_type(edges):
    """BFS selected-edge components; induced edges are deliberately irrelevant."""
    adjacency = {}
    for u, v in edges:
        u, v = (0, u), (1, v)
        adjacency.setdefault(u, set()).add(v)
        adjacency.setdefault(v, set()).add(u)
    remaining, signature = set(adjacency), []
    while remaining:
        seen, stack = set(), [next(iter(remaining))]
        while stack:
            vertex = stack.pop()
            if vertex not in seen:
                seen.add(vertex)
                stack.extend(adjacency[vertex] - seen)
        remaining -= seen
        edge_count = sum(len(adjacency[x]) for x in seen) // 2
        row_count = sum(x[0] == 0 for x in seen)
        col_count = len(seen) - row_count
        signature.append(("C" if edge_count == len(seen) else "P",
                          edge_count, row_count - col_count))
    return tuple(sorted(signature))


def type_normalization(signature):
    rows = cols = 0
    aut = 1
    for (kind, length, imbalance), multiplicity in Counter(signature).items():
        aut *= factorial(multiplicity)
        if kind == "C":
            r = c = length // 2
            component_aut = length
        elif length % 2:
            r = c = (length + 1) // 2
            component_aut = 1
        else:
            r = (length + 1 + imbalance) // 2
            c = (length + 1 - imbalance) // 2
            component_aut = 2
        rows += r * multiplicity
        cols += c * multiplicity
        aut *= component_aut ** multiplicity
    return rows, cols, aut


def diagonal_union_bank(n):
    events = []
    for d in range(1 - n, n):
        diagonal = [(r, r + d) for r in range(n) if 0 <= r + d < n]
        events.extend(map(frozenset, combinations(diagonal, 3)))
    bank = Counter()
    for i, event in enumerate(events):
        for j in range(i, len(events)):
            signature = graph_type(event | events[j])
            need(all(kind == "P" for kind, _, _ in signature),
                 "two parallel diagonals have no finite incidence cycle")
            bank[signature] += 1 if i == j else 2
    need(sum(bank.values()) == len(events) ** 2, "ordered event-pair multiplicity")
    return bank


def union_copy_variance(parts, target_bank):
    """Second path: selected skeleton edge subsets / all grid copies of the type."""
    n = sum(parts)
    edges, _ = skeleton(parts)
    source_bank = Counter()
    for k in range(3, 7):
        for selected in combinations(edges, k):
            source_bank[graph_type(selected)] += 1
    second = Q(0)
    for signature, count in target_bank.items():
        rows, cols, aut = type_normalization(signature)
        probability = Q(source_bank[signature] * aut, falling(n, rows) * falling(n, cols))
        need(0 <= probability <= 1, "non-induced grid-copy probability normalization")
        second += count * probability
    mean = Q(2 * n - 5, 3)
    return second - mean * mean


def label_census(parts):
    n = sum(parts)
    edges, _ = skeleton(parts)
    hist, joint, xhist = Counter(), Counter(), Counter()
    pset = list(permutations(range(n)))
    total = factorial(n) ** 2
    witness = None
    conditional_gates = 0
    for rho in pset:
        local = [0] * (n + 1)
        holder = [0, 0, 0]
        for sigma in pset:
            board = tuple((rho[u], sigma[v]) for u, v in edges)
            counts = Counter(v - u for u, v in board)
            score = sum(comb(y, 3) for y in counts.values() if y >= 3)
            hist[score] += 1
            joint[(counts[0], counts[1])] += 1
            for k in range(n + 1):
                local[k] += falling(counts[0], k)
            holder[0] += 2 ** (counts[0] + counts[1])
            holder[1] += 4 ** counts[0]
            holder[2] += 4 ** counts[1]
            if n == 4:
                x = full_collinear(board)
                xhist[x] += 1
                need(score <= x, "literal determinant contains slope-one triples")
                if score == 0 and x > 0 and witness is None:
                    witness = (tuple(sorted(board)), x)
        for k in range(n + 1):
            need(Q(local[k], factorial(n)) <= 2 ** k, "conditional single-line exponential bound")
            conditional_gates += 1
        need(holder[0] ** 2 <= holder[1] * holder[2], "conditional Holder, weighted union retained")
        conditional_gates += 1
    need(sum(hist.values()) == total, "complete shore-label universe")
    mean = sum(Q(x * count, total) for x, count in hist.items())
    second = sum(Q(x * x * count, total) for x, count in hist.items())
    need(mean == Q(2 * n - 5, 3), "exact diagonal mean")
    variance = second - mean * mean
    # An independent two-edge inclusion formula keeps shared row and column sites.
    L, M, R, C = n, n - 1, n - 1, n - 1
    p2 = Q(4 * n - 6, n * (n - 1) ** 2)
    q2 = Q(2, n * (n - 1))
    direct11 = sum(Q(y * z * count, total) for (y, z), count in joint.items())
    need(direct11 == (L * M - R - C) * p2 + (R + C) * q2,
         "mixed two-edge inclusion, overlaps retained")
    cov11 = direct11 - Q(2 * L, n) * Q(2 * M, n)
    need(cov11 != 0, "finite diagonal counts are dependent")
    return {"parts": parts, "labels": total, "mean_S": str(mean),
            "variance_S": str(variance), "P_S_zero": str(Q(hist[0], total)),
            "cov_Y0_Y1": str(cov11), "hist_S": dict(sorted(hist.items())),
            "hist_X_if_n4": dict(sorted(xhist.items())), "S_zero_X_positive": witness,
            "conditional_gates": conditional_gates}


def main():
    bank = []
    # Independent edge-subset and closed cycle-rook paths on every cycle type n<=8.
    for n in range(2, 9):
        for parts in parts_at_least_two(n):
            edges, owner = skeleton(parts)
            got, formula = exact_matchings(edges), rook_product(parts)
            need(got == formula, "edge subsets versus cycle rook product")
            if n >= 3:
                need(got[3] == Q(2 * n * (n - 2) * (2 * n - 5), 3), "cycle-blind m3")
            if n >= 4:
                need(got[4] == Q(n * (n - 3) * (2 * n - 7) * (2 * n - 5), 6)
                     + parts.count(2), "first cycle sidecar is c4")
            if all(x == 2 for x in parts) and n >= 4:
                occupancy = Counter()
                for chosen in combinations(edges, 3):
                    if len({u for u, _ in chosen}) == 3 and len({v for _, v in chosen}) == 3:
                        occupancy[tuple(sorted(Counter(owner[e] for e in chosen).values()))] += 1
                p21 = Q(6 * occupancy[(1, 2)], falling(n, 3) ** 2)
                p111 = Q(6 * occupancy[(1, 1, 1)], falling(n, 3) ** 2)
                need(p21 == Q(12, n * (n - 2) * (n - 1) ** 2), "all-C4 two-block probability")
                need(p111 == Q(8 * (n - 4), n * (n - 2) * (n - 1) ** 2),
                     "all-C4 three-block probability")
            bank.append((parts, got))

    n, z, theta, eta, overlap = sp.symbols("n z theta eta overlap")
    # Exact single-cycle rook polynomial independently checks the general local expansion.
    for k in range(1, 7):
        ordered = 2 * n * sp.prod(2 * n - k - j for j in range(1, k))
        exact = ordered / sp.prod(n - j for j in range(k)) ** 2
        relative = sp.cancel(exact / (2 / n) ** k).subs(n, 1 / z)
        series = sp.series(relative, z, 0, 2).removeO().expand()
        need(series == 1 + Q(k * (k - 1), 4) * z, "matching probability first correction")

    # The covariance expansion retains the one-path term.  Symbols overlap=(R+C)/n.
    for a in range(1, 4):
        for b in range(1, 4):
            k = a + b
            A = falling(theta / z, a) * falling(eta / z, b)
            D = a * b * overlap / z * (theta / z) ** (a - 1) * (eta / z) ** (b - 1)
            pk = (2 * z) ** k * (1 + Q(k * (k - 1), 4) * z)
            pa = (2 * z) ** a * (1 + Q(a * (a - 1), 4) * z)
            pb = (2 * z) ** b * (1 + Q(b * (b - 1), 4) * z)
            qk = (2 * z) ** k / 2
            cov = sp.expand((A - D) * pk + D * qk - A * pa * pb)
            want = a * b * 2 ** (k - 1) * (
                theta ** a * eta ** b - theta ** (a - 1) * eta ** (b - 1) * overlap)
            need(sp.expand(cov.coeff(z, 1) - want) == 0, "uniform mixed covariance coefficient")

    # The square of (Y)_3 uses genuine multiplication, including repeated events.
    y = sp.symbols("y")
    need(sp.expand(falling(y, 3) ** 2 - sum(
        comb(3, j) ** 2 * factorial(j) * falling(y, 6 - j) for j in range(4))) == 0,
        "falling factorial square identity")
    x, y = sp.symbols("x y", positive=True)
    same = 2 * sp.integrate(y ** 2 * sp.integrate(x ** 3, (x, 0, y)), (y, 0, 1))
    opposite = sp.integrate(sp.integrate(x ** 2 * y ** 2 * (x + y - 1),
                                        (y, 1 - x, 1)), (x, 0, 1))
    local_var = 2 * sp.integrate(8 * x ** 5 + 8 * x ** 4 + sp.Rational(4, 3) * x ** 3,
                                 (x, 0, 1))
    cross = 2 - 32 * (same + opposite)
    need(same == sp.Rational(1, 14), "same-side diagonal overlap integral")
    need(opposite == sp.Rational(71, 1260), "opposite-side diagonal overlap integral")
    need(local_var == sp.Rational(98, 15), "individual-line variance integral")
    need(cross == -sp.Rational(94, 45), "negative aggregate cross covariance")
    need(local_var + cross == sp.Rational(40, 9), "uniform variance coefficient")
    need((local_var + cross) / sp.Rational(2, 3) ** 2 == 10, "Chebyshev coefficient")

    # Exact finite diagonal geometry, including the two ends and the common main diagonal.
    for size in range(3, 81):
        lengths = [size - abs(d) for d in range(1 - size, size)]
        triples = sum(comb(L, 3) for L in lengths if L >= 3)
        need(triples == comb(size, 3) + 2 * comb(size, 4), "all slope-one triples")
        need(Q(4 * (2 * size - 5) * triples,
               size * (size - 1) ** 2 * (size - 2)) == Q(2 * size - 5, 3),
             "all-size exact diagonal mean")

    # Hostile to an invalid binary-permutation treatment of a union of target matchings.
    # Here rho is fixed, T0 and T1 are the two cyclic perfect matchings of G.
    edges7, _ = skeleton((7,))
    cyclic_target = set(edges7)
    weighted_second = 0
    for sigma in permutations(range(7)):
        hits = sum((u, sigma[v]) in cyclic_target for u, v in edges7)
        weighted_second += falling(hits, 2)
    weighted_second = Q(weighted_second, factorial(7))
    need(weighted_second == Q(49, 3) and weighted_second > 16,
         "weighted-union conditional factorial bound is false")
    census = [label_census(parts) for parts in ((3,), (4,), (2, 2), (5,), (2, 3))]
    n4 = {tuple(row["parts"]): row for row in census if sum(row["parts"]) == 4}
    need(n4[(4,)]["variance_S"] != n4[(2, 2)]["variance_S"],
         "universal asymptotic variance does not erase finite cycle effects")
    need(n4[(4,)]["S_zero_X_positive"] is not None, "diagonal gate is only sufficient for failure")
    need(Q(n4[(2, 2)]["hist_X_if_n4"].get(0, 0), 576) == Q(1, 2), "frozen 2C4 zero rate")
    need(Q(n4[(4,)]["hist_X_if_n4"].get(0, 0), 576) == Q(1, 36), "frozen C8 zero rate")
    union_banks = {n: diagonal_union_bank(n) for n in (3, 4, 5, 8)}
    copy_replay = []
    for row in census:
        parts = tuple(row["parts"])
        variance = union_copy_variance(parts, union_banks[sum(parts)])
        need(variance == Q(row["variance_S"]), "copy-union replay versus all labelings")
        copy_replay.append((parts, str(variance)))
    for parts in ((8,), (2, 2, 2, 2)):
        copy_replay.append((parts, str(union_copy_variance(parts, union_banks[8]))))
    payload = {"matching_carriers": bank, "census": census,
               "copy_replay": copy_replay,
               "integrals": list(map(str, (same, opposite, local_var, cross, local_var + cross)))}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    print("STATUS: PASS; exact controls plus companion all-size proof, no census extrapolation")
    print("MATCHING CARRIERS:", len(bank), "all cycle types, 2<=n<=8, degrees 0..6")
    for row in census:
        print(json.dumps(row, sort_keys=True))
    print("INDEPENDENT NON-INDUCED UNION-COPY VARIANCES:", copy_replay)
    print("VARIANCE COEFFICIENTS: single=98/15; cross=-94/45; total=40/9; Chebyshev=10")
    print("CONDITIONAL WEIGHTED-UNION HOSTILE: n=7, factorial second moment=49/3>16; Holder survives")
    print("SCOPE: uniform random-label density bound; no nonexistence, CLT, or exponential rate")
    print("SEMANTIC SHA256:", digest)
    print("ACTIVE GATES:", GATES)


if __name__ == "__main__":
    main()
