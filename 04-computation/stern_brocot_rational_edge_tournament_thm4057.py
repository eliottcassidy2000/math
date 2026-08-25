#!/usr/bin/env python3
"""Exact audit for THM-4057.

This companion types reduced rationals as directed coprime-graph arcs, checks
the reciprocal Stern--Brocot mirror and trace/norm/discriminant sidecar,
pulls Stern--Brocot depth parity through the three Berggren branches, and
tests the resulting depth-gauged tournament.  Empirical tournament profiles
are printed only as FINITE-EXACT data.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations
from math import gcd, prod


GATES = 0


def check(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def canonical_cf(p: int, q: int) -> tuple[int, ...]:
    check(p > 0 and q > 0, "CF positive")
    digits: list[int] = []
    while q:
        digits.append(p // q)
        p, q = q, p % q
    check(len(digits) == 1 or digits[-1] >= 2, "canonical CF tail")
    return tuple(digits)


def sb_path(p: int, q: int) -> str:
    check(p > 0 and q > 0 and gcd(p, q) == 1, "SB primitive")
    reverse_path: list[str] = []
    while p != q:
        if p < q:
            reverse_path.append("L")
            q -= p
        else:
            reverse_path.append("R")
            p -= q
    check((p, q) == (1, 1), "SB terminates at root")
    return "".join(reversed(reverse_path))


@lru_cache(maxsize=None)
def sb_depth(p: int, q: int) -> int:
    g = gcd(p, q)
    p //= g
    q //= g
    return sum(canonical_cf(p, q)) - 1


def swap_lr(path: str) -> str:
    return path.translate(str.maketrans("LR", "RL"))


Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[int, int]


def mat_vec(m: Matrix, v: Vector) -> Vector:
    return (m[0][0] * v[0] + m[0][1] * v[1],
            m[1][0] * v[0] + m[1][1] * v[1])


def mat_mul(a: Matrix, b: Matrix) -> Matrix:
    return (
        (a[0][0] * b[0][0] + a[0][1] * b[1][0],
         a[0][0] * b[0][1] + a[0][1] * b[1][1]),
        (a[1][0] * b[0][0] + a[1][1] * b[1][0],
         a[1][0] * b[0][1] + a[1][1] * b[1][1]),
    )


def mat_pow(m: Matrix, exponent: int) -> Matrix:
    ans: Matrix = ((1, 0), (0, 1))
    while exponent:
        if exponent & 1:
            ans = mat_mul(ans, m)
        m = mat_mul(m, m)
        exponent //= 2
    return ans


def determinant(m: Matrix) -> int:
    return m[0][0] * m[1][1] - m[0][1] * m[1][0]


A: Matrix = ((0, 1), (-1, 2))
B: Matrix = ((0, 1), (1, 2))
C: Matrix = ((1, 0), (2, 1))
BRANCHES: tuple[tuple[str, Matrix], ...] = (("A", A), ("B", B), ("C", C))
ROOT: Vector = (1, 2)


def gaussian(p: int, q: int) -> tuple[int, int, int]:
    return (q * q - p * p, 2 * p * q, p * p + q * q)


def sigma_depth(a: int, b: int) -> int:
    check(a > 0 and b > 0 and a != b, "tournament distinct positive")
    g = gcd(a, b)
    parity = -1 if sb_depth(a // g, b // g) % 2 else 1
    order = 1 if b > a else -1
    return order * parity


def arc(a: int, b: int) -> bool:
    return sigma_depth(a, b) > 0


def directed_triangles(n: int) -> list[tuple[int, int, int]]:
    cycles: list[tuple[int, int, int]] = []
    for a, b, c in combinations(range(1, n + 1), 3):
        if ((arc(a, b) and arc(b, c) and arc(c, a))
                or (arc(a, c) and arc(c, b) and arc(b, a))):
            cycles.append((a, b, c))
    return cycles


def score_profile(n: int) -> tuple[int, ...]:
    scores = []
    for a in range(1, n + 1):
        scores.append(sum(1 for b in range(1, n + 1) if b != a and arc(a, b)))
    return tuple(sorted(scores))


print("THM-4057 Stern--Brocot / rational-edge tournament audit")
print()

# 1. Reduced rationals are arcs of the coprimality graph, not a tournament.
bound = 48
rationals = {
    (p, q)
    for p in range(1, bound + 1)
    for q in range(1, bound + 1)
    if p != q and gcd(p, q) == 1
}
coprime_arcs = {
    (a, b)
    for a in range(1, bound + 1)
    for b in range(1, bound + 1)
    if a != b and gcd(a, b) == 1
}
check(rationals == coprime_arcs, "rational-edge bijection")
for a, b in rationals:
    check((b, a) in rationals, f"reciprocal arc {a}/{b}")
check((2, 4) not in rationals and (4, 2) not in rationals,
      "noncoprime pair has no raw arc")
natural_arcs = {(a, b) for a, b in rationals if a < b}
for a, b in natural_arcs:
    check(a < b, "natural orientation increasing")
check(all(not (b < a) for a, b in natural_arcs), "natural DAG")
pairwise_coprime = (2, 3, 5, 7)
check(all(gcd(a, b) == 1 for a, b in combinations(pairwise_coprime, 2)),
      "pairwise-coprime clique")
print(f"bounded rational arcs (endpoints <= {bound}): {len(rationals)}")
print("raw object: bidirected on coprime pairs, empty on noncoprime pairs")
print("p<q restriction: acyclic coprimality-graph orientation, not tournament")
print()

# 2. Stern--Brocot depth, reciprocal mirror, and Khinchin-type even data.
sb_cases = 0
for p in range(1, 121):
    for q in range(1, 121):
        if gcd(p, q) != 1:
            continue
        path = sb_path(p, q)
        reciprocal_path = sb_path(q, p)
        depth = sb_depth(p, q)
        check(depth == len(path), f"CF/path depth {p}/{q}")
        check(reciprocal_path == swap_lr(path), f"global mirror {p}/{q}")
        check(sb_depth(q, p) == depth, f"reciprocal depth {p}/{q}")
        digits = canonical_cf(p, q)
        reciprocal_digits = canonical_cf(q, p)
        positive_digits = digits[1:] if digits[0] == 0 else digits
        reciprocal_positive = (
            reciprocal_digits[1:] if reciprocal_digits[0] == 0
            else reciprocal_digits
        )
        check(positive_digits == reciprocal_positive,
              f"reciprocal positive digits {p}/{q}")
        check(prod(positive_digits) == prod(reciprocal_positive),
              f"product of full positive digit list {p}/{q}")
        if p < q:
            check(sb_depth(q - p, q) == depth,
                  f"subtree reflection depth {p}/{q}")
        sb_cases += 1
print(f"primitive mirror/depth cases: {sb_cases}")
print("reciprocal swaps L/R and preserves depth and positive CF digits: PASS")
print()

# 3. Trace/norm/discriminant sidecar and Gaussian-square orientation loss.
carrier_cases = 0
for p in range(1, 101):
    for q in range(1, 101):
        if p == q:
            continue
        trace = p + q
        norm = p * q
        delta = q - p
        check(trace * trace - 4 * norm == delta * delta,
              f"discriminant identity {p},{q}")
        check((trace - delta) // 2 == p and (trace + delta) // 2 == q,
              f"oriented reconstruction {p},{q}")
        check((gcd(p, q) == 1) == (gcd(trace, norm) == 1),
              f"coprime trace/norm {p},{q}")
        triple = gaussian(p, q)
        check(triple == (trace * delta, 2 * norm, trace * trace - 2 * norm),
              f"Gaussian carrier {p},{q}")
        reversed_triple = gaussian(q, p)
        check(reversed_triple == (-triple[0], triple[1], triple[2]),
              f"Gaussian reversal {p},{q}")
        if gcd(p, q) == 1:
            content = gcd(gcd(abs(triple[0]), triple[1]), triple[2])
            expected = 1 if (p - q) % 2 else 2
            check(content == expected, f"Pythagorean content {p},{q}")
        carrier_cases += 1
print(f"trace/norm/sheet and Gaussian cases: {carrier_cases}")
print("reversal fixes trace/norm and negates sheet/first Gaussian coordinate: PASS")
print()

# 4. Berggren branches pull Stern--Brocot depth into the A-Walsh character.
nodes: list[tuple[str, Vector]] = [("", ROOT)]
berggren_nodes = 0
print("h  nodes  depth-odd  depth-even  signed-sum")
for h in range(0, 9):
    odd = 0
    even = 0
    signed_sum = [0, 0]
    for word, u in nodes:
        d = sb_depth(*u)
        expected_depth = 1 + 2 * h - word.count("A")
        check(d == expected_depth, f"Berggren depth word={word}")
        sign = -1 if d % 2 else 1
        signed_sum[0] += sign * u[0]
        signed_sum[1] += sign * u[1]
        odd += d % 2
        even += 1 - d % 2
        check(gcd(*u) == 1 and 0 < u[0] < u[1],
              f"positive primitive Berggren word={word}")
        berggren_nodes += 1
    expected_vector = mat_vec(mat_pow(C, 2 * h), ROOT)
    check(tuple(signed_sum) == (-expected_vector[0], -expected_vector[1]),
          f"depth Walsh sum h={h}")
    check(odd == (3 ** h + 1) // 2, f"odd-depth count h={h}")
    check(even == (3 ** h - 1) // 2, f"even-depth count h={h}")
    print(f"{h:1d} {len(nodes):6d} {odd:10d} {even:11d} {tuple(signed_sum)}")
    next_nodes: list[tuple[str, Vector]] = []
    for word, u in nodes:
        parent_depth = sb_depth(*u)
        for name, matrix in BRANCHES:
            child = mat_vec(matrix, u)
            increment = 1 if name == "A" else 2
            check(sb_depth(*child) == parent_depth + increment,
                  f"branch increment {word}{name}")
            next_nodes.append((word + name, child))
    nodes = next_nodes
print(f"Berggren nodes checked across levels: {berggren_nodes}")
print("D(w(1/2))=1+2h-#A and pulled-back Walsh ray: PASS")
print()

# The depth cocycle holds for every reduced root in the positive chamber, not
# only for the standard 1/2 root.  The signed level sum acquires precisely the
# initial checkerboard sign.
general_roots: tuple[Vector, ...] = ((1, 3), (2, 5), (3, 7), (5, 8), (8, 13))
general_root_cases = 0
for root in general_roots:
    level: list[tuple[str, Vector]] = [("", root)]
    d0 = sb_depth(*root)
    for h in range(0, 6):
        signed = [0, 0]
        for word, u in level:
            d = sb_depth(*u)
            check(d == d0 + 2 * h - word.count("A"),
                  f"general depth root={root},word={word}")
            sign = -1 if d % 2 else 1
            signed[0] += sign * u[0]
            signed[1] += sign * u[1]
            general_root_cases += 1
        ray = mat_vec(mat_pow(C, 2 * h), root)
        initial_sign = -1 if d0 % 2 else 1
        check(tuple(signed) == (initial_sign * ray[0], initial_sign * ray[1]),
              f"general Walsh root={root},h={h}")
        level = [
            (word + name, mat_vec(matrix, u))
            for word, u in level
            for name, matrix in BRANCHES
        ]
print(f"all-root depth-cocycle nodes checked: {general_root_cases}")
print("D(wx)=D(x)+2h-#A and signed C^(2h) ray for five roots: PASS")
print()

# 5. Distinguish the three reflections used near this bridge.
H: Matrix = ((-1, 1), (1, 1))
TWO_I: Matrix = ((2, 0), (0, 2))
check(mat_mul(H, H) == TWO_I, "H square")
check(determinant(H) == -2, "H determinant")
check(mat_mul(H, A) == mat_mul(C, H), "H swaps A to C")
check(mat_mul(H, B) == mat_mul(B, H), "H fixes B")
check(mat_mul(H, C) == mat_mul(A, H), "H swaps C to A")
sample = (1, 3)
check((sample[1], sample[0]) != (sample[1] - sample[0], sample[1]),
      "global reciprocal differs from j")
check(abs(determinant(H)) != 1, "H not unimodular tree automorphism")
print("J:x->1/x, j:x->1-x, H:x->(1-x)/(1+x) are distinct: PASS")
print()

# 6. An explicit scale-invariant, nontransitive tournament gauge.
for n in range(2, 61):
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            if a == b:
                continue
            check(sigma_depth(a, b) == -sigma_depth(b, a),
                  f"antisymmetry {a},{b}")
            check(sigma_depth(a, b) in (-1, 1), f"totality {a},{b}")
            for scale in (2, 3, 5):
                check(sigma_depth(scale * a, scale * b) == sigma_depth(a, b),
                      f"scale invariance {scale}*({a},{b})")

first_triangle: tuple[int, int, int] | None = None
cycle_rows: list[tuple[int, int]] = []
for n in range(3, 21):
    cycles = directed_triangles(n)
    cycle_rows.append((n, len(cycles)))
    if cycles and first_triangle is None:
        first_triangle = cycles[0]
check(first_triangle == (1, 2, 5), "first directed triangle")
check(arc(2, 1) and arc(1, 5) and arc(5, 2), "triangle orientation")
check(not arc(1, 2), "depth gauge changes raw 1->2 direction")
print(f"first directed triangle: 2->1->5->2 on {first_triangle}")
print(f"directed C3 counts n=3..20: {[count for _, count in cycle_rows]}")
for n in (5, 10, 15, 20):
    print(f"score profile n={n}: {score_profile(n)}")
print("depth-gauged tournament antisymmetry/totality/scaling: PASS")
print()

# 7. Exact Pell and Fibonacci cycle families inside the depth tournament.
pell = [0, 1]
fib = [0, 1]
for _ in range(2, 31):
    pell.append(2 * pell[-1] + pell[-2])
    fib.append(fib[-1] + fib[-2])

pell_cycles: list[tuple[int, int, int]] = []
for j in range(3, 20, 2):
    check(sb_depth(pell[j], pell[j - 1]) == 2 * j - 3,
          f"Pell adjacent depth j={j}")
    check(arc(1, pell[j]) and arc(pell[j], pell[j - 1])
          and arc(pell[j - 1], 1), f"Pell cycle j={j}")
    pell_cycles.append((1, pell[j], pell[j - 1]))

fibonacci_cycles: list[tuple[int, int, int]] = []
for n in range(2, 25):
    check(sb_depth(fib[n + 1], fib[n]) == n - 1,
          f"Fibonacci adjacent depth n={n}")
    if n % 6 == 5:
        check(arc(1, fib[n]) and arc(fib[n], fib[n + 1])
              and arc(fib[n + 1], 1), f"Fibonacci cycle n={n}")
        fibonacci_cycles.append((1, fib[n], fib[n + 1]))
    elif n % 6 == 0:
        check(arc(1, fib[n + 1]) and arc(fib[n + 1], fib[n])
              and arc(fib[n], 1), f"Fibonacci cycle n={n}")
        fibonacci_cycles.append((1, fib[n + 1], fib[n]))

print(f"Pell directed-cycle family (first four): {pell_cycles[:4]}")
print(f"Fibonacci directed-cycle family (first four): {fibonacci_cycles[:4]}")
print("Pell all odd j>=3; Fibonacci n=0 or 5 mod 6 cycle laws: PASS")
print()

print(f"TOTAL EXACT GATES: {GATES}")
print("ALL CHECKS PASSED")
