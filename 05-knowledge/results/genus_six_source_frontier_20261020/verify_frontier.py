#!/usr/bin/env python3
"""Exact genus-six frontier census and explicit polynomial certificates.

The finite characteristic universe is justified in FRONTIER.md.
No producer or predecessor code is imported. No bare assertions are used.
"""
import itertools
import json
import math
from functools import lru_cache

import sympy as S

checks = []


def check(name, condition):
    if not condition:
        raise RuntimeError("FAILED: " + name)
    checks.append(name)


def zero(name, expression):
    check(name, S.expand(expression.doit()) == 0)


@lru_cache(None)
def monoid(gens, bound):
    have = {0}
    for n in range(1, bound + 1):
        if any(n >= g and n - g in have for g in gens):
            have.add(n)
    return frozenset(have)


def primitives(vals):
    return tuple(n for n in sorted(vals) if n > 0 and not any(
        a in vals and n - a in vals for a in range(1, n)))


expected_gens = [
    (2, 13), (3, 7), (4, 5), (4, 6, 9),
    (5, 7, 8, 9), (6, 7, 8, 9, 10),
]
expected_gap_sets = sorted([
    (1, 3, 5, 7, 9, 11), (1, 2, 4, 5, 8, 11),
    (1, 2, 3, 6, 7, 11), (1, 2, 3, 5, 7, 11),
    (1, 2, 3, 4, 6, 11), (1, 2, 3, 4, 5, 11),
])

# First census: all 4096 subsets, without assuming symmetric pair choices.
census = []
for mask in range(1 << 12):
    vals = {n for n in range(12) if mask & (1 << n)}
    if 0 not in vals or 11 in vals or len(vals) != 6:
        continue
    if any((n in vals) == (11 - n in vals) for n in range(12)):
        continue
    if any(a + b < 12 and a + b not in vals for a in vals for b in vals):
        continue
    census.append(tuple(n for n in range(12) if n not in vals))
check("all 4096 subsets give exactly six symmetric genus-six semigroups",
      sorted(census) == expected_gap_sets)

# Independent census: choose one member in each of the six complementary pairs.
second_census = []
for bits in itertools.product((0, 1), repeat=6):
    vals = {i if bit == 0 else 11 - i for i, bit in enumerate(bits)}
    if 0 not in vals:
        continue
    if all(a + b >= 12 or a + b in vals for a in vals for b in vals):
        second_census.append(tuple(n for n in range(12) if n not in vals))
check("independent complementary-pair census agrees",
      sorted(second_census) == expected_gap_sets)

semigroup_records = []
for gens in expected_gens:
    have = monoid(gens, 80)
    gaps = sorted(set(range(81)) - have)
    m = min(gens)
    apery = [min(n for n in have if n % m == r) for r in range(m)]
    check("primitive generators " + str(gens), primitives(have) == gens)
    check("six gaps and Frobenius eleven " + str(gens),
          len(gaps) == 6 and max(gaps) == 11)
    check("conductor twelve from Apery " + str(gens),
          max(apery) - m + 1 == 12)
    check("Apery genus formula " + str(gens),
          2 * sum(apery) == m * (12 + m - 1))
    check("conductor interval certificate " + str(gens),
          set(range(12, 12 + m)) <= have)
    semigroup_records.append({
        "minimal_generators": gens, "gaps": gaps, "apery_by_residue": apery,
    })
check("minimal generator census agrees with all gap sets",
      sorted(tuple(r["gaps"]) for r in semigroup_records) == expected_gap_sets)


def characteristic(seq):
    d = seq[0]
    previous_n = None
    conductor = 1 - seq[0]
    for i, value in enumerate(seq[1:], 1):
        new_d = math.gcd(d, value)
        n = d // new_d
        if i >= 2 and n <= 1:
            return None
        if n * value not in monoid(seq[:i], n * value):
            return None
        if i >= 2 and value >= previous_n * seq[i - 1]:
            return None
        conductor += (n - 1) * value
        d, previous_n = new_d, n
    return conductor if d == 1 else None


pairs = [(a, b) for a in range(1, 14) for b in range(1, 14)
         if a % b and b % a]
check("complete initial universe has gcds one through four",
      {math.gcd(*pair) for pair in pairs} == {1, 2, 3, 4})
check("length-four sequences cannot be omitted a priori",
      characteristic((8, 12, 10, 5)) == 20)
expected_sequences = sorted([
    (2, 13), (13, 2), (3, 7), (7, 3), (4, 5), (5, 4),
    (4, 10, 5), (10, 4, 5), (4, 6, 9), (6, 4, 9),
    (6, 9, 4), (9, 6, 4),
])

# First path is a direct Cartesian-product census of lengths two, three, four.
found = []
tested_by_length = {}
for length in (2, 3, 4):
    count = 0
    for a, b in pairs:
        for tail in itertools.product(range(1, 12), repeat=length - 2):
            seq = (a, b) + tail
            count += 1
            if characteristic(seq) == 12:
                found.append(seq)
    tested_by_length[length] = count
check("entire bounded characteristic universe at conductor twelve",
      sorted(found) == expected_sequences)
check("length-four sector contributes no conductor-twelve sequence",
      all(len(seq) <= 3 for seq in found))
for seq in found:
    have = monoid(seq, 80)
    check("sequence semigroup agrees with conductor formula " + str(seq),
          len(set(range(81)) - have) == 6 and 11 not in have
          and set(range(12, 81)) <= have)
    check("every primitive occurs in the sequence " + str(seq),
          set(primitives(have)) <= set(seq))

# Second path: target-filtered gcd recursion, testing monoids rather than
# using the characteristic conductor formula to select a candidate.
recursive = []
for gens in expected_gens:
    allowed = sorted(monoid(gens, 13) - {0})

    def extend(seq):
        d = math.gcd(*seq)
        if d == 1:
            have = monoid(seq, 80)
            if primitives(have) == gens and characteristic(seq) is not None:
                recursive.append(seq)
            return
        for value in allowed:
            if value <= 11 and math.gcd(d, value) < d:
                extend(seq + (value,))

    for a, b in itertools.permutations(allowed, 2):
        if a % b and b % a:
            extend((a, b))
check("independent target-filtered recursion agrees",
      sorted(recursive) == expected_sequences)
for nonplanar in expected_gens[-2:]:
    check("nonplanar symmetric semigroup " + str(nonplanar),
          all(primitives(monoid(seq, 80)) != nonplanar for seq in found))

t, u, v = S.symbols("t u v")
realizations = [
    {"sequence": (2, 13), "p": t**2, "q": t**13 + t, "r": None},
    {"sequence": (3, 7), "p": t**3, "q": t**7 + t, "r": None},
    {"sequence": (4, 5), "p": t**4, "q": t**5 + t, "r": None},
    {"sequence": (4, 6, 9), "p": t**4, "q": t**6 + t**3 + t,
     "r": v**2 - u**3},
    {"sequence": (4, 10, 5), "p": t**4 + t,
     "q": t**10 + S.Rational(5, 2)*t**7 + S.Rational(15, 8)*t**4
          + S.Rational(5, 16)*t,
     "r": v**2 - u**5 - S.Rational(5, 64)*u**2},
    {"sequence": (6, 9, 4), "p": t**6 + t**4,
     "q": t**9 + S.Rational(3, 2)*t**7 + S.Rational(3, 8)*t**5
          - S.Rational(5, 32)*t**3 - S.Rational(3, 128)*t,
     "r": v**2 - u**3 + S.Rational(3, 16)*u**2
          - S.Rational(7, 1024)*u},
]
realization_records = []
for item in realizations:
    seq, p, q, expression = [item[k] for k in ("sequence", "p", "q", "r")]
    n = int(S.degree(p, t))
    check("actual coordinate degrees " + str(seq),
          (n, int(S.degree(q, t))) == seq[:2])
    check("immersion derivative gcd " + str(seq),
          S.gcd(S.diff(p, t), S.diff(q, t)) == 1)
    if expression is None:
        basis_expr = [v**i for i in range(n)]
        primitive_degrees = seq
        remainder = None
    else:
        remainder = S.expand(expression.subs({u: p, v: q}))
        check("exact remainder degree " + str(seq),
              S.degree(remainder, t) == seq[2])
        primitive_degrees = (n, int(S.degree(q, t)),
                             int(S.degree(remainder, t)))
        basis_expr = [v**i * expression**j
                      for j in range(n // 2) for i in range(2)]
        transition = S.Matrix([
            [S.Poly(expr, v).nth(i) for expr in basis_expr]
            for i in range(n)])
        check("full monogenic module transition determinant " + str(seq),
              S.expand(transition.det()) == 1)
    check("primitive fraction field by degree gcd " + str(seq),
          math.gcd(*primitive_degrees) == 1)
    basis = [S.expand(expr.subs({u: p, v: q})) for expr in basis_expr]
    degrees = [int(S.degree(expr, t)) for expr in basis]
    check("full module distinct residues " + str(seq),
          sorted(d % n for d in degrees) == list(range(n)))
    check("module gap count exactly six " + str(seq),
          sum(d - d % n for d in degrees) // n == 6)
    observed = {d + n*j for d in degrees for j in range(81) if d+n*j <= 80}
    check("module proves entire claimed semigroup " + str(seq),
          observed == monoid(seq, 80))
    # Independent resultant: implicit relation has the required q-degree,
    # corroborating the field-degree + integrality proof of the full basis.
    resultant = S.Poly(S.resultant(p-u, q-v, t), v)
    check("resultant has full expected degree " + str(seq),
          resultant.degree() == n)
    realization_records.append({
        "sequence": seq, "p": str(p), "q": str(q),
        "remainder": None if remainder is None else str(remainder),
        "basis_degrees": degrees, "module_codimension": 6,
    })

r410 = realizations[4]["r"].subs(
    {u: realizations[4]["p"], v: realizations[4]["q"]})
zero("exact quartic-decic alternate remainder", r410-t**5/64-5*t**2/256)
check("quartic-decic derivative remainder",
      S.rem(S.diff(realizations[4]["q"], t),
            S.diff(realizations[4]["p"], t), t) == -S.Rational(5, 8))
r69 = realizations[5]["r"].subs(
    {u: realizations[5]["p"], v: realizations[5]["q"]})
zero("exact sextic-nonic alternate remainder",
     r69-t**4/2048-9*t**2/16384)
check("sextic-nonic derivative at zero",
      S.diff(realizations[5]["q"], t).subs(t, 0) == -S.Rational(3, 128))
check("sextic-nonic derivative at other critical points",
      S.rem(S.diff(realizations[5]["q"], t), 3*t*t+2, t)
      == -S.Rational(27, 128))

# Sharp genus-five mechanism boundary: derive the two odd-degree defects.
a, b, c, d, e, h, f, A5, B3 = S.symbols("a b c d e h f A5 B3")
p = t**4+a*t*t+b*t
q = t**6+A5*t**5+c*t**4+B3*t**3+e*t*t+h*t
r = S.Poly(q*q-p**3, t)
zero("first odd defect at degree eleven", r.nth(11)-2*A5)
q0 = q.subs(A5, 0)
r0 = S.Poly(q0*q0-p**3, t)
r0 = S.Poly(r0.as_expr()-r0.nth(10)*p*q0, t)
zero("second odd defect at degree nine", r0.nth(9)-(2*B3-3*b))
for degree, delta in [(11, 7), (9, 6), (7, 5)]:
    degrees = [0, 6, degree, 6+degree]
    check("odd-defect module codimension " + str(degree),
          sum(d-d % 4 for d in degrees)//4 == delta)
    check("odd-defect distinct module residues " + str(degree),
          sorted(d % 4 for d in degrees) == list(range(4)))

# Independent formal jet ring computes only the two needed top coefficients.
ad, bd, cd, dd, ed, hd, fd, Bd = S.symbols("ad bd cd dd ed hd fd Bd")
P = t**4+a*t*t+b*t+d
Q = t**6+c*t**4+B3*t**3+e*t*t+h*t+f
JP = ad*t*t+bd*t+dd
JQ = cd*t**4+Bd*t**3+ed*t*t+hd*t+fd
J = S.Poly(S.expand(JP*S.diff(Q, t)-S.diff(P, t)*JQ), t)
zero("upper Jacobian coefficient seven", J.nth(7)-(6*ad-4*cd))
zero("upper Jacobian coefficient six", J.nth(6)-(6*bd-4*Bd))
zero("constant-defect differential gate", J.nth(6)+2*(2*Bd-3*bd))
check("delta-six hostile escapes prior B=3b/2 condition",
      realizations[3]["p"] == t**4 and
      S.Poly(realizations[3]["q"], t).nth(3) == 1)

print(json.dumps({
    "status": "PASS",
    "scope": "genus-six classification and explicit realizations; no new Keller bound",
    "check_count": len(checks),
    "checks": checks,
    "symmetric_genus_six": semigroup_records,
    "characteristic_universe_counts_by_length": tested_by_length,
    "nonprincipal_sequences": expected_sequences,
    "realizations": realization_records,
    "keller_scope": {
        "excluded_by_inherited_quadratic_theorem": [[2, 13]],
        "remaining_unoriented_terminal_pairs": [[3, 7], [4, 5], [4, 6],
                                               [4, 10], [6, 9]],
        "new_conductor_bound_claimed": False,
    },
}, indent=2, sort_keys=True))
