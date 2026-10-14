"""Exact marked-run quotient census for the oriented four-rim wheel.

No arithmetic interpolation: enumerate every decomposition of the eight
labelled wheel edges into ordered trails of prescribed lengths, then divide
by the wheel automorphism group (8). Coefficients are partition-Mobius weights.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations
from math import factorial
from functools import lru_cache
from pathlib import Path
import json

EDGES = [(0, i) for i in range(1, 5)] + [(1, 2), (2, 3), (3, 4), (4, 1)]
ADJ = [[] for _ in range(5)]
for ei, (u, v) in enumerate(EDGES):
    ADJ[u].append((v, ei, 1))
    ADJ[v].append((u, ei, -1))


def padd(a, b):
    out = [Fraction(0)] * max(len(a), len(b))
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    while len(out) > 1 and out[-1] == 0: out.pop()
    return tuple(out)


def pscale(a, c):
    return tuple(c * v for v in a)


def linear(a, c):
    return padd((Fraction(0),) + a, pscale(a, c))


def normalize(n, edges):
    """Injectivity: doubled edges square to1; oriented reversal gives -1."""
    sign, parity = 1, set()
    for a, b in edges:
        if a == b: return 0, n, ()
        if a > b: a, b, sign = b, a, -sign
        edge = (a, b)
        if edge in parity: parity.remove(edge)
        else: parity.add(edge)
    return sign, n, tuple(sorted(parity))


def evaluate(n, edges):
    """Internal evaluator for reductions of positive runs with <=8 edges.

    The public census does not admit extra isolated source vertices. Extending
    it to wheel-plus-isolate inputs needs a polynomial wheel coefficient.
    """
    sign, n, edges = normalize(n, edges)
    if sign == 0: return ((Fraction(0),), Fraction(0))
    p, c = reduced(n, edges)
    return pscale(p, sign), sign * c


@lru_cache(None)
def reduced(n, edges):
    """Exact injective sum=P(q)+c*unanchored-injective-wheel for q=3mod4."""
    if n == 0: return ((Fraction(1),), Fraction(0))
    neighbors = [[] for _ in range(n)]
    for a, b in edges:
        neighbors[a].append(b); neighbors[b].append(a)
    v = min(range(n), key=lambda a: (len(neighbors[a]), -a))
    degree = len(neighbors[v])
    if degree >= 3:
        if n == 4 and len(edges) == 6:
            # Transposing two vertices reverses the Vandermonde product.
            return ((Fraction(0),), Fraction(0))
        assert n == 5 and len(edges) == 8, (n, edges)
        hub = next(a for a in range(n) if len(neighbors[a]) == 4)
        rim0 = min(a for a in range(n) if a != hub)
        rim1 = min(a for a in neighbors[rim0] if a != hub)
        rim2 = next(a for a in neighbors[rim1] if a not in (hub, rim0))
        rim3 = next(a for a in range(n) if a not in (hub, rim0, rim1, rim2))
        model = [(hub, a) for a in (rim0, rim1, rim2, rim3)]
        model += [(rim0, rim1), (rim1, rim2), (rim2, rim3), (rim3, rim0)]
        sign = (-1) ** sum(a > b for a, b in model)
        assert set(tuple(sorted(e)) for e in model) == set(edges)
        return ((Fraction(0),), Fraction(sign))
    relabel = lambda a: a - (a > v)
    retained = [(relabel(a), relabel(b)) for a, b in edges if v not in (a, b)]
    if degree == 0:
        p, c = evaluate(n - 1, retained)
        assert c == 0  # Eight edges cannot leave an isolated vertex and a wheel here.
        return linear(p, -(n - 1)), Fraction(0)
    result, wheel = (Fraction(0),), Fraction(0)
    if degree == 2:
        # Orient both incident factors outward from v, costing (-1)^#incoming.
        orient = (-1) ** sum(a < v for a in neighbors[v])
        # Distinct neighbours: sum_x chi((a-x)(b-x))=-1.
        p, c = evaluate(n - 1, retained)
        result = padd(result, pscale(p, -orient))
        wheel -= orient * c
    # Full sum for degree1 is0; subtract every occupied non-neighbour value.
    for j in range(n):
        if j == v or j in neighbors[v]: continue
        extra = []
        for a, b in edges:
            if v in (a, b):
                extra.append((relabel(j if a == v else a),
                              relabel(j if b == v else b)))
        p, c = evaluate(n - 1, retained + extra)
        result = padd(result, pscale(p, -1))
        wheel -= c
    return result, wheel


def run_edges(lengths):
    n, edges = 0, []
    for length in lengths:
        edges.extend((n + i, n + i + 1) for i in range(length))
        n += length + 1
    return n, edges


def chi_table(q):
    squares = {x * x % q for x in range(1, q)}
    return [0] + [1 if x in squares else -1 for x in range(1, q)]


def poly_value(poly, q):
    ans = Fraction(0)
    for c in reversed(poly): ans = ans * q + c
    return ans


def brute_injective(q, lengths):
    """Independent affine-normalized tuple sum, no graph reductions."""
    n, edges = run_edges(lengths)
    if n > q: return 0
    chi = chi_table(q)
    ans = 0
    for tail in permutations(range(2, q), n - 2):
        values = (0, 1) + tail
        term = 1
        for a, b in edges: term *= chi[(values[b] - values[a]) % q]
        ans += term
    return q * (q - 1) * ans


def brute_wheel(q):
    chi = chi_table(q)
    ans = 0
    for tail in permutations(range(2, q), 3):
        values = (0, 1) + tail
        term = 1
        for a, b in EDGES: term *= chi[(values[b] - values[a]) % q]
        ans += term
    return q * (q - 1) * ans


def held_karp_coefficients(q):
    """Independent exact permutation product, fixing its initial vertex0."""
    chi = chi_table(q)
    states = {(1, 0): [1]}
    for length in range(1, q):
        nxt = {}
        for (mask, last), coeffs in states.items():
            for vertex in range(1, q):
                if mask & (1 << vertex): continue
                key = mask | (1 << vertex), vertex
                dst = nxt.setdefault(key, [0] * (length + 1))
                sign = chi[(vertex - last) % q]
                for i, value in enumerate(coeffs):
                    dst[i] += value; dst[i + 1] += sign * value
        states = nxt
    summed = [sum(coeffs[i] for coeffs in states.values()) for i in range(q)]
    return [Fraction(c, factorial(q - 1)) for c in summed]


def census(lengths):
    """Return coefficients for chi(-1)=+1 and -1, independently of q."""
    counts = Counter()
    terms = []
    def start(run, mask, word, sign):
        if run == len(lengths):
            if mask != 255:
                return
            sizes = Counter(word)
            assert len(sizes) == 5
            mu = 1
            for v in sizes.values():
                mu *= (-1) ** (v - 1) * factorial(v - 1)
            counts["labelled"] += 1
            counts["plus"] += mu
            counts["minus"] += mu * sign
            terms.append((word, sign, mu))
            return
        for vertex in range(5):
            walk(run, lengths[run], vertex, mask, word + (vertex,), sign)
    def walk(run, remaining, vertex, mask, word, sign):
        if remaining == 0:
            start(run + 1, mask, word, sign)
            return
        for nxt, edge, direction in ADJ[vertex]:
            if mask & (1 << edge) == 0:
                walk(run, remaining - 1, nxt, mask | (1 << edge),
                     word + (nxt,), sign * direction)
    start(0, 0, (), 1)
    # Independent quotient path: canonical restricted-growth partition keys.
    partitions = {}
    for word, sign, mu in terms:
        labels = {}
        key = tuple(labels.setdefault(v, len(labels)) for v in word)
        if key in partitions:
            assert partitions[key] == (sign, mu)
        partitions[key] = (sign, mu)
    assert len(terms) == 8 * len(partitions)
    assert sum(mu for _, mu in partitions.values()) * 8 == counts["plus"]
    assert sum(s * mu for s, mu in partitions.values()) * 8 == counts["minus"]
    return {"runs": lengths, "partitions": len(partitions),
            "coefficient_q1mod4": counts["plus"] // 8,
            "coefficient_q3mod4": counts["minus"] // 8}


def main():
    rows = [census(lengths) for lengths in
            [(8,), (7, 1), (6, 2), (5, 3), (4, 4),
             (4, 2, 2), (2, 2, 2, 2)]]
    out = {"status": "PROVED exact reduction identities; FINITE-EXACT symbolic coefficient census",
           "scope": "all odd prime powers q=3mod4, at most eight edges; prime modular substitution imported separately",
           "wheel_edges": EDGES, "rows": rows}
    polynomials = []
    total_p, total_w = (Fraction(0),), Fraction(0)
    weights = {(8,): Fraction(1), (6, 2): Fraction(1),
               (4, 4): Fraction(1, 2), (4, 2, 2): Fraction(1, 2),
               (2, 2, 2, 2): Fraction(1, 24)}
    for lengths, weight in weights.items():
        p, c = evaluate(*run_edges(lengths))
        expected = next(r["coefficient_q3mod4"] for r in rows if tuple(r["runs"]) == lengths)
        assert c == expected, (lengths, c, expected)
        polynomials.append({"runs": lengths, "polynomial_ascending": list(map(str, p)),
                            "injective_wheel_coefficient": str(c)})
        total_p = padd(total_p, pscale(p, weight)); total_w += weight * c
    assert total_w == 7
    out["exact_injective_run_decompositions"] = polynomials
    out["S8_polynomial_ascending"] = list(map(str, total_p))
    out["S8_injective_wheel_coefficient"] = str(total_w)
    smaller = {}
    for lengths in [(2,), (4,), (6,), (2, 2), (4, 2), (2, 2, 2)]:
        poly, wheel_coefficient = evaluate(*run_edges(lengths))
        assert wheel_coefficient == 0
        smaller[str(lengths)] = list(map(str, poly))
        for q in (3, 7, 11):
            assert poly_value(poly, q) == brute_injective(q, lengths)
    out["smaller_run_polynomials_ascending"] = smaller
    out["reducer_cache"] = str(reduced.cache_info())
    checks = []
    for q in (3, 7, 11):
        wheel = brute_wheel(q)
        for lengths in weights:
            p, c = evaluate(*run_edges(lengths))
            expected = poly_value(p, q) + c * wheel
            actual = brute_injective(q, lengths)
            assert expected == actual, (q, lengths, expected, actual)
            checks.append({"q": q, "runs": lengths, "sum": actual})
    for q, expected in [(19, 10243584), (23, 34327040), (31, 202487040)]:
        p, c = evaluate(*run_edges((8,)))
        assert c == 0 and poly_value(p, q) == expected
    hk = held_karp_coefficients(11)
    expected8 = (poly_value(total_p, 11) + total_w * brute_wheel(11)) / (factorial(11) // factorial(3))
    assert hk[8] == expected8
    assert sum(hk) * factorial(11) / 2 ** 10 == 95095
    assert all(hk[i] == 0 for i in range(1, 11, 2))
    out["independent_injective_checks"] = checks
    out["independent_Held_Karp_R11_coefficients"] = list(map(str, hk))
    out["R11_degree8"] = str(expected8)
    crt_modulus = 2147483647 * 2147483629 * 2147483587
    assert factorial(16) ** 30 < 10 ** 400 < crt_modulus ** 16
    n1 = 42745655163116955307561
    assert [n1 % p for p in (2147483647, 2147483629, 2147483587)] == [1586941875, 1066060341, 2021511562]
    assert n1 * 465 == 19876729650849384218015865
    out["H31_CRT_uniqueness"] = {"modulus": crt_modulus,
        "N1_upper_bound_strict": 10 ** 25,
        "integer_check": "(16!)^30 < (10^25)^16 < M^16",
        "N1_residue_reconstruction": n1,
        "scope": "uniqueness of incoming residues, not an independent rerun of the q31 DP"}
    path = Path(__file__).resolve().parents[1] / "05-knowledge/results/wheel_run_interface_20261014h.json"
    path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
