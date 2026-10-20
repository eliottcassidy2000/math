#!/usr/bin/env python3
"""Independent audit: Kunz census, shortest-path Apery sets, exact identities.

Does not import either producer checker. No assertions are used.
All files written by the audit are confined to this directory.
"""
from functools import lru_cache
import hashlib
import heapq
import itertools
import json
from math import gcd
from pathlib import Path

import sympy as s

checks = []


def require(label, condition):
    if not bool(condition):
        raise RuntimeError("FAIL: " + label)
    checks.append(label)


def identity(label, expression):
    require(label, s.expand(expression) == 0)


def positive_compositions(total, length):
    if length == 1:
        if total > 0:
            yield (total,)
        return
    for first in range(1, total - length + 2):
        for tail in positive_compositions(total - first, length - 1):
            yield (first,) + tail


@lru_cache(None)
def apery(gens):
    """Dijkstra on residues, exact for gcd-one generators."""
    require_gcd = gcd(*gens)
    if require_gcd != 1:
        raise ValueError("Apery requested for a nonnumerical monoid")
    m = min(gens)
    dist = [None] * m
    dist[0] = 0
    queue = [(0, 0)]
    while queue:
        value, residue = heapq.heappop(queue)
        if value != dist[residue]:
            continue
        for g in gens:
            new = value + g
            r = new % m
            if dist[r] is None or new < dist[r]:
                dist[r] = new
                heapq.heappush(queue, (new, r))
    return tuple(dist)


def belongs(n, gens):
    d = gcd(*gens)
    if n < 0 or n % d:
        return False
    normalized = tuple(sorted(set(g // d for g in gens)))
    values = apery(normalized)
    n //= d
    return n >= values[n % len(values)]


def invariants(gens):
    ap = apery(tuple(sorted(set(gens))))
    m = len(ap)
    gaps = tuple(n for n in range(max(ap)) if n < ap[n % m])
    primitive = tuple(n for n in range(1, max(ap) + m + 1)
                      if n >= ap[n % m]
                      and not any(belongs(j, gens) and belongs(n-j, gens)
                                  for j in range(1, n)))
    return primitive, gaps, max(ap) - m + 1, ap


# Independent semigroup census. For genus six, m<=7. Positive Kunz
# coordinates sum to six; all Kunz inequalities are checked before symmetry.
kunz_records = []
composition_count = 0
for m in range(2, 8):
    for coords in positive_compositions(6, m - 1):
        composition_count += 1
        x = (0,) + coords
        if any(x[i] + x[j] + ((i+j)//m) < x[(i+j) % m]
               for i in range(1, m) for j in range(1, m)):
            continue
        ap = tuple(m*x[r] + r for r in range(m))
        conductor = max(ap) - m + 1
        if conductor != 12:
            continue
        gens = (m,) + ap[1:]
        primitive, gaps, c, ap2 = invariants(gens)
        require("Kunz conductor symmetry " + str(primitive),
                c == 12 and len(gaps) == 6
                and all((n in gaps) != (11-n in gaps)
                        for n in range(12)))
        kunz_records.append({"generators": primitive, "gaps": gaps, "apery": ap2})

expected_gens = [(2, 13), (3, 7), (4, 5), (4, 6, 9),
                 (5, 7, 8, 9), (6, 7, 8, 9, 10)]
require("Kunz census, all multiplicities two through seven",
        sorted(r["generators"] for r in kunz_records) == expected_gens)
require("Kunz universe has 32 positive compositions", composition_count == 32)


def characteristic(seq):
    ds = [seq[0]]
    ns = []
    for value in seq[1:]:
        ds.append(gcd(ds[-1], value))
        ns.append(ds[-2] // ds[-1])
    if ds[-1] != 1 or any(n <= 1 for n in ns[1:]):
        return False
    if any(not belongs(ns[i-1] * seq[i], seq[:i])
           for i in range(1, len(seq))):
        return False
    return all(seq[i] < ns[i-2] * seq[i-1]
               for i in range(2, len(seq)))


# Recursion never imposes a sequence-length cutoff. It stops precisely
# when the gcd becomes one, and otherwise adds only a strict gcd drop.
terminals = []
initial = [(a, b) for a in range(1, 14) for b in range(1, 14)
           if a % b and b % a]


def extend(seq):
    d = gcd(*seq)
    if d == 1:
        if characteristic(seq):
            terminals.append(seq)
        return
    for tail in range(1, 12):
        if gcd(d, tail) < d:
            extend(seq + (tail,))


for pair in initial:
    extend(pair)
sequences = sorted(seq for seq in terminals if invariants(seq)[2] == 12)
expected_seq = sorted([(2, 13), (13, 2), (3, 7), (7, 3), (4, 5), (5, 4),
                       (4, 10, 5), (10, 4, 5), (4, 6, 9), (6, 4, 9),
                       (6, 9, 4), (9, 6, 4)])
require("All genus-six characteristic sequences, exact Apery selection",
        sequences == expected_seq)
require("Initial universe and maximum gcd",
        len(initial) == 108 and max(gcd(*p) for p in initial) == 4)
require("No valid recursion chain exceeds length four",
        max(map(len, terminals)) == 4)
require("Length-four positive control", characteristic((8, 12, 10, 5))
        and invariants((8, 12, 10, 5))[2] == 20)
require("Length-four genus-six sector empty",
        not any(len(seq) == 4 for seq in sequences))
for seq in sequences:
    ds = [gcd(*seq[:i]) for i in range(1, len(seq)+1)]
    formula = 1-seq[0] + sum((ds[i-1]//ds[i]-1)*seq[i]
                            for i in range(1, len(seq)))
    primitive, gaps, c, ap = invariants(seq)
    require("Independent conductor formula " + str(seq), formula == c == 12)
    require("Primitive occurrence and genus " + str(seq),
            set(primitive) <= set(seq) and len(gaps) == 6)
require("Exactly four planar semigroups",
        sorted(set(invariants(seq)[0] for seq in sequences)) == expected_gens[:4])
require("Quadratic exclusion leaves exactly five degree pairs",
        sorted(set(tuple(sorted(seq[:2])) for seq in sequences
                   if invariants(seq)[0] != (2, 13)))
        == [(3, 7), (4, 5), (4, 6), (4, 10), (6, 9)])

t, a, b, c, e, h, A, B, X, Y, U, V, z = s.symbols(
    "t a b c e h A B X Y U V z")
R = s.Rational
examples = [
    ((2, 13), t**2, t**13+t, None),
    ((3, 7), t**3, t**7+t, None),
    ((4, 5), t**4, t**5+t, None),
    ((4, 6, 9), t**4, t**6+t**3+t, Y**2-X**3),
    ((4, 10, 5), t**4+t,
     t**10+R(5, 2)*t**7+R(15, 8)*t**4+R(5, 16)*t,
     Y**2-X**5-R(5, 64)*X**2),
    ((6, 9, 4), t**6+t**4,
     t**9+R(3, 2)*t**7+R(3, 8)*t**5-R(5, 32)*t**3-R(3, 128)*t,
     Y**2-X**3+R(3, 16)*X**2-R(7, 1024)*X),
]
for seq, p, q, rexpr in examples:
    n = seq[0]
    require("Actual degree pair " + str(seq),
            (s.degree(p, t), s.degree(q, t)) == seq[:2])
    require("Immersion over algebraic closure " + str(seq),
            s.gcd(s.diff(p, t), s.diff(q, t)) == 1)
    if rexpr is None:
        basis = [q**j for j in range(n)]
    else:
        r = s.expand(rexpr.subs({X: p, Y: q}))
        require("Cancellation degree " + str(seq), s.degree(r, t) == seq[2])
        basis = [r**j * q**i for j in range(n//2) for i in range(2)]
    degrees = tuple(int(s.degree(s.expand(f), t)) for f in basis)
    require("Field index degree-gcd certificate " + str(seq), gcd(*seq) == 1)
    require("Whole-module residue separation " + str(seq),
            sorted(d % n for d in degrees) == list(range(n)))
    # The residue minima below are the exact Apery set modulo n, found by
    # taking the smaller multiplicity Apery set and solving each residue.
    exact_min = tuple(next(k for k in range(r, max(degrees)+n+1, n)
                           if belongs(k, seq)) for r in range(n))
    require("All-degree semigroup equality " + str(seq),
            tuple(sorted(degrees, key=lambda x: x % n)) == exact_min)
    require("Exact module codimension " + str(seq),
            sum(d-d % n for d in degrees) == 6*n)

# Universal monic-quadratic transition, including the linear-q term needed
# by the beta-defect argument. U,V may be arbitrary polynomials in p.
for n in (4, 6):
    expressions = [Y**i * (Y**2-U*Y-V)**j
                   for j in range(n//2) for i in range(2)]
    matrix = s.Matrix([[s.expand(f).coeff(Y, k) for f in expressions]
                       for k in range(n)])
    identity("Universal full module transition " + str(n), matrix.det()-1)

p4 = t**4+b*t
q10 = t**10+R(5, 2)*b*t**7+R(15, 8)*b*b*t**4+R(5, 16)*b**3*t
identity("Quartic family cancellation",
         q10*q10-p4**5-R(5, 64)*b**4*p4*p4
         - b**5*t**5/64 - 5*b**6*t*t/256)
identity("Quartic family derivative Euclidean certificate",
         s.rem(s.diff(q10, t), s.diff(p4, t), t)+5*b**3/8)
require("Quartic exceptional composition index", gcd(4, 10) == 2
        and p4.subs(b, 0) == t**4 and q10.subs(b, 0) == t**10)

p6 = t**6+a*t**4+b*t*t
q9 = (t**9+R(3, 2)*a*t**7+(R(3, 2)*b+R(3, 8)*a*a)*t**5
      +a*(-5*a*a+36*b)*t**3/32+3*(16*b*b-a**4)*t/128)
D = a*a-4*b
AA = -3*a*D/16
BB = D*D*(7*a*a+8*b)/1024
remainder = s.expand(q9*q9-p6**3-AA*p6*p6-BB*p6)
identity("Sextic family cancellation",
         remainder-D**3*t*t*(9*a*a+8*a*t*t-4*b)/16384)
identity("Genus-six leading remainder coefficient",
         remainder.coeff(t, 4)-a*D**3/2048)
G = 3*z*z+2*a*z+b
qp = s.Poly(s.diff(q9, t), t)
Q1 = sum(coef*z**(monom[0]//2) for monom, coef in qp.terms())
identity("Parity derivative reconstruction", Q1.subs(z, t*t)-s.diff(q9, t))
identity("Sextic derivative at zero",
         Q1.subs(z, 0)+3*D*(a*a+4*b)/128)
identity("Sextic derivative remainder",
         s.rem(Q1, G, z)-3*D*(-a*a+12*a*z+12*b)/128)
identity("Sextic derivative resultant",
         s.resultant(G, Q1, z)-R(6561, 16384)*D**4)
require("Nonimmersed genus-six divisor has only the zero critical point",
        s.gcd(s.diff(p6, t).subs(b, -a*a/4),
              s.diff(q9, t).subs(b, -a*a/4)) == t)
identity("Genus-four degeneration coefficient",
         remainder.subs(a, 0)-b**4*t*t/64)
identity("Genus-four odd generator",
         q9.subs(a, 0)-(t**9+R(3, 2)*b*t**5+R(3, 8)*b*b*t))
hh = t**3+a*t/2
identity("Nonprimitive divisor p is cubic square", p6.subs(b, a*a/4)-hh*hh)
identity("Nonprimitive divisor q is cubic cube", q9.subs(b, a*a/4)-hh**3)
require("Genus-four degeneration exact semigroup",
        invariants((2, 9))[:3] == ((2, 9), (1, 3, 5, 7), 8))

# Derive the parity ansatz afresh from all high coefficients rather than
# assuming the claimed formulas for the t^7 and t^5 coefficients.
q7, q5, q3, q1, ca, cb = s.symbols("q7 q5 q3 q1 ca cb")
qans = t**9+q7*t**7+q5*t**5+q3*t**3+q1*t
raw = s.Poly(s.expand(qans*qans-p6**3-ca*p6*p6-cb*p6), t)
solved = {}
for degree, variable in [(16, q7), (14, q5), (12, ca), (10, q1)]:
    solutions = s.solve(raw.nth(degree).subs(solved), variable)
    require("Unique parity elimination at degree " + str(degree), len(solutions) == 1)
    solved[variable] = s.expand(solutions[0])
identity("Parity degree-eight factored obstruction",
         raw.nth(8).subs(solved)-D*(5*a**3-36*a*b+32*q3)/128)
solved[q3] = a*(-5*a*a+36*b)/32
identity("Parity generic family forced", qans.subs(solved)-q9)
identity("Parity generic A forced", solved[ca].subs(solved)-AA)
bs = s.solve(raw.nth(6).subs(solved), cb)
require("Parity degree-six B determined", len(bs) == 1)
identity("Parity generic B forced", bs[0]-BB)
exceptional = {key: val for key, val in solved.items() if key != q3}
cc = q3-a**3/8
identity("Parity exceptional composition q",
         qans.subs(exceptional).subs(b, a*a/4)-hh**3-cc*hh)
identity("Parity exceptional A",
         exceptional[ca].subs(b, a*a/4)-2*cc)
identity("Parity exceptional B and vanishing full remainder",
         (qans.subs(exceptional).subs(b, a*a/4))**2
         -hh**6-2*cc*hh**4-cc*cc*hh*hh)

# Derive the quartic/sextic boundary including actual scalar terms.
p = t**4+a*t*t+b*t
q = t**6+A*t**5+c*t**4+B*t**3+e*t*t+h*t
identity("Odd eleven defect", s.expand(q*q-p**3).coeff(t, 11)-2*A)
q = q.subs(A, 0)
r = s.expand(q*q-p**3-(2*c-3*a)*p*q)
identity("Odd nine defect", r.coeff(t, 9)-(2*B-3*b))
identity("Degree ten removed", r.coeff(t, 10))
for defect, delta in [(11, 7), (9, 6)]:
    degs = [0, 6, defect, defect+6]
    require("Odd defect exact codimension " + str(defect),
            sorted(x % 4 for x in degs) == list(range(4))
            and sum(x-x % 4 for x in degs) == 4*delta)
x = s.symbols("x")
af, bf, cf, Bf, df, ef, hf, ff = [
    s.Function(name)(x) for name in ("a", "b", "c", "B", "d", "e", "h", "f")]
P = t**4+af*t*t+bf*t+df
Q = t**6+cf*t**4+Bf*t**3+ef*t*t+hf*t+ff
J = s.Poly(s.expand(s.diff(P, x)*s.diff(Q, t)-s.diff(P, t)*s.diff(Q, x)), t)
identity("Full scalar-retaining Jacobian degree seven",
         J.nth(7)-6*s.diff(af, x)+4*s.diff(cf, x))
identity("Full scalar-retaining Jacobian degree six",
         J.nth(6)+2*s.diff(2*Bf-3*bf, x))

packet = Path(__file__).resolve().parents[1]
workspace = packet.parent
dependency_hashes = json.loads((packet/"DEPENDENCY_HASHES.json").read_text())
for relative, expected in dependency_hashes.items():
    require("Frozen dependency hash " + relative,
            hashlib.sha256((workspace/relative).read_bytes()).hexdigest() == expected)

print(json.dumps({
    "status": "PASS",
    "check_count": len(checks),
    "checks": checks,
    "independent_semigroup_method": "All 32 positive Kunz compositions; exact residue shortest paths",
    "semigroups": kunz_records,
    "characteristic_method": "Gcd-drop recursion with no length cutoff; exact Apery conductor selection",
    "nonprincipal_characteristic_sequences": sequences,
    "source_independence": "Neither supplied checker is imported",
    "scope": "Genus-six classification, realizations, exact beta boundary; no new Keller bound",
}, sort_keys=True, indent=2))
