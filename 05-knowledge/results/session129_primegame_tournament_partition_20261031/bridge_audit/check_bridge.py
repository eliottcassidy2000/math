#!/usr/bin/env python3
"""Exact, standard-library audit of the proposed PRIMEGAME block bridge.

No tournament census is repeated here. All fractions and assertions are local
to the user's fourteen fractions and the stated four-block partition.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import permutations, product
import json
import sys
from pathlib import Path

LETTERS = "ABCDEFGHIJKLMN"  # Conway program order, not user block names.
F = list(map(Q, ["17/91", "78/85", "19/51", "23/38", "29/33",
                  "77/29", "95/23", "77/19", "1/17", "11/13",
                  "13/11", "15/2", "1/7", "55"]))
BLOCKS = {"A": "JK", "B": "NMIA", "C": "CHFE", "D": "BDGL"}
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
CONTROL = [11, 13, 17, 19, 23, 29]
VERTICES = [0] + CONTROL

def require(condition, message):
    """Always-active verification gate, including under python -O."""
    if not condition:
        raise AssertionError(message)

if "--test-failing-gate" in sys.argv:
    require(False, "intentional negative control: gate must remain active")

def vp(n, p):
    z = 0
    while n % p == 0:
        n //= p
        z += 1
    return z

def valuation(q, p):
    return vp(q.numerator, p) - vp(q.denominator, p)

def mul(xs):
    z = Q(1)
    for x in xs:
        z *= x
    return z

def rank(rows):
    a = [[Q(x) for x in row] for row in rows]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        d = a[r][c]
        a[r] = [x / d for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                d = a[i][c]
                a[i] = [x - d*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]

def execute_word(n, word):
    trajectory = [n]
    for letter in word:
        x = n * F[LETTERS.index(letter)]
        if x.denominator != 1:
            return None
        n = x.numerator
        trajectory.append(n)
    return trajectory

def priority_steps(n, steps):
    word = ""
    trajectory = [n]
    for _ in range(steps):
        for letter, f in zip(LETTERS, F):
            x = n*f
            if x.denominator == 1:
                word += letter
                n = x.numerator
                trajectory.append(n)
                break
        else:
            raise AssertionError("The terminal 55 fraction always applies")
    return word, trajectory

V = [[valuation(f, p) for f in F] for p in PRIMES]
edges, gains = [], []
for f in F:
    neg = [p for p in CONTROL if valuation(f, p) == -1]
    pos = [p for p in CONTROL if valuation(f, p) == 1]
    require(all(abs(valuation(f, p)) <= 1 for p in CONTROL), "verification gate 91")
    require(len(neg) <= 1 and len(pos) <= 1, "verification gate 92")
    source, target = (neg or [0])[0], (pos or [0])[0]
    edges.append((source, target))
    gains.append(f * Q(source or 1, target or 1))

boundary = [[int(v == t) - int(v == s) for s, t in edges]
            for v in VERTICES]
require(boundary[1:] == V[4:], "verification gate 99")
require(rank(boundary) == 6 and rank(V) == 10, "verification gate 100")
simple_edges = {tuple(sorted(e)) for e in edges if e[0] != e[1]}
theta = {tuple(sorted((u, v))) for u in [11, 17] for v in [0, 13, 19]}
require(simple_edges == theta | {(11, 29), (19, 23)}, "verification gate 103")

# Integral chain map, with the zero-boundary A generator unbased.
P = [[int(letter in word) for word in BLOCKS.values()] for letter in LETTERS]
macro_vertices = [11, 13, 17]
d = [[0, 1, 1, 0], [0, -1, 0, 1], [0, 0, -1, -1]]
J = [[int(v == w) for w in macro_vertices] for v in VERTICES]
require(multiply(boundary, P) == multiply(J, d), "verification gate 110")
require(rank(P) == 4 and rank(d) == 2, "verification gate 111")
products = {g: mul(F[LETTERS.index(a)] for a in word) for g, word in BLOCKS.items()}
require(products == {"A": Q(1), "B": Q(55, 637), "C": Q(539, 153), "D": Q(585, 34)}, "verification gate 113")
macroV = multiply(V, P)
require(rank(macroV) == 3, "verification gate 115")
holonomy = products["B"] * products["D"] / products["C"]
require(holonomy == Q(2025, 4802) != 1, "verification gate 117")
def decode_macro(q):
    b = sum(valuation(q, p) for p in [2, 5, 23])
    return b, valuation(q, 11) - b, valuation(q, 13) + b

for b, c, dd in product(range(-4, 5), repeat=3):
    q = products["B"]**b * products["C"]**c * products["D"]**dd
    require(decode_macro(q) == (b, c, dd), "verification gate 124")
cycle_basis = [[1, 0, 0, 0], [0, 1, -1, 1]]
for c in cycle_basis:
    require(multiply(d, [[x] for x in c]) == [[0], [0], [0]], "verification gate 127")

# Priority semantics: two positive controls and the disconnected-cycle hostile.
positive = {}
for group, n, expected_word, expected_end in [
        ("A", 11, "KJ", 11), ("B", 637, "AIMN", 55),
        ("C", 153, "CHEF", 539)]:
    word, trajectory = priority_steps(n, len(expected_word))
    require(word == expected_word and trajectory[-1] == expected_end, "verification gate 135")
    require(Counter(word) == Counter(BLOCKS[group]), "verification gate 136")
    positive[group] = {"program_word": word, "trajectory": trajectory}
require(34*products["D"] == 585, "verification gate 138")
valid_D_permutations = ["".join(w) for w in permutations(BLOCKS["D"])
                        if execute_word(34, w) is not None]
require(valid_D_permutations == [], "verification gate 141")
require(priority_steps(34, 1)[0] == "I", "verification gate 142")

# Strong independent audit of the claimed nonnegative multiplicative kernel.
T = [V[0][i] + V[2][i] + V[8][i] for i in range(14)]
require(T == [0]*13 + [1], "verification gate 146")
kernel_words = ["JK", "ABDEFG", "AABCEFHKL", "ABEEFFLM"]
kernel_vectors = [[word.count(a) for a in LETTERS] for word in kernel_words]
for word in kernel_words:
    require(mul(F[LETTERS.index(a)] for a in word) == 1, "verification gate 150")
require(rank(kernel_vectors) == 4, "verification gate 151")
for c, dd, m, j in product(range(5), repeat=4):
    # Complete parametrization; converse and integrality are proved in report.
    z = [2*c+dd+m, c+dd+m, c, dd, c+dd+2*m, c+dd+2*m,
         dd, c, 0, j, c+j, c+m, m, 0]
    require(all(sum(V[r][i]*z[i] for i in range(14)) == 0 for r in range(10)), "verification gate 156")
    expected = [sum(q*vec[i] for q, vec in zip([j, dd, c, m], kernel_vectors))
                for i in range(14)]
    require(z == expected, "verification gate 159")

# Label-free graph symmetry check: no graph-automorphism excuse for arbitrary
# identification of arithmetic labels with tournament isomorphism classes.
directed_counts = Counter(edges)
directed_auts = []
simple_auts = []
for perm in permutations(VERTICES):
    f = dict(zip(VERTICES, perm))
    if Counter((f[s], f[t]) for s, t in edges) == directed_counts:
        directed_auts.append(list(perm))
    if {tuple(sorted((f[s], f[t]))) for s, t in simple_edges} == simple_edges:
        simple_auts.append(list(perm))
require(len(directed_auts) == 1 and len(simple_auts) == 2, "verification gate 172")

result = {
    "status": "PASS",
    "universe": "14 exact fractions; 10 prime valuations; four user blocks",
    "program_order": dict(zip(LETTERS, map(str, F))),
    "user_blocks_in_program_letters": BLOCKS,
    "control_edges": [{"letter": a, "source": s, "target": t, "gain": str(g)}
                      for a, (s, t), g in zip(LETTERS, edges, gains)],
    "group_products": {g: str(p) for g, p in products.items()},
    "macro_boundary_rows_11_13_17": d,
    "chain_map_verified": True,
    "ranks": {"control_incidence": 6, "control_cycle": 8,
              "full_valuation": 10, "multiplicative_kernel": 4,
              "macro_boundary": 2, "macro_cycles": 2,
              "macro_valuation": 3, "macro_multiplicative_kernel": 1},
    "macro_cycle_basis": cycle_basis,
    "triangle_holonomy": str(holonomy),
    "block_product_codec": {
        "encoding": "q = B^b C^c D^d",
        "decoding": "b=T(q), c=nu_11(q)-T(q), d=nu_13(q)+T(q)",
        "definition_T": "nu_2+nu_5+nu_23",
        "selected_valuation_matrix": [[1,0,0],[1,1,0],[-1,0,1]],
        "exact_roundtrips_checked": 729,
        "scope": "arbitrary integer exponent triples; proof uses integer left inverse"
    },
    "priority_positive_controls": positive,
    "D_hostile": {"input": 34, "formal_output": 585,
                  "executable_permutations": valid_D_permutations,
                  "full_program_first_step": priority_steps(34, 1)},
    "nonnegative_multiplicative_kernel_generators": kernel_words,
    "kernel_parametrization_checked_inputs": 625,
    "monotone_T_increment_by_program_letter": dict(zip(LETTERS, T)),
    "simple_graph_automorphisms": len(simple_auts),
    "directed_multigraph_automorphisms": len(directed_auts),
}
print(json.dumps(result, indent=2))
