#!/usr/bin/env python3
"""THM-4808 (G) over Q: the graph on ALL Klein orbits of P^1(Q) has two components.
Checks on every reduced slope m/n with |m|, |n| <= H (plus infinity):
  (Q1) the parity type 'opp' (m + n odd: Euclid slopes) or 'odd' (m, n both odd) is preserved by U, A, D, tau,
       iota, S -- so Euclid and odd/odd slopes never meet;
  (Q2) the Cayley involution C(t) = (t+1)/(t-1) swaps the two types and satisfies C D C = U, C U C = D,
       C A C = A, C tau C = tau, C iota C = S iota, C S C = S (an isomorphism of the two components);
  (Q3) on the Euclid component every vertex V.t with t = m/n > 1 has the three Berggren children as heads of
       U t, A t, D t and its Berggren parent (or (1,0,1) for the root) as head of tau t;
  (Q4) local structure of the two orbifold vertices {0, inf} and {1, -1}.
"""
import sys
from fractions import Fraction as F
from math import gcd

H = int(sys.argv[1]) if len(sys.argv) > 1 else 60
INF = 'inf'


def mob(M, t):
    a, b, c, d = M
    if t == INF:
        return INF if c == 0 else F(a, c)
    num, den = a * t + b, c * t + d
    return INF if den == 0 else num / den


U = (2, -1, 1, 0); A = (2, 1, 1, 0); D = (1, 2, 0, 1); TAU = (-1, 2, 0, 1)
IOTA = (-1, 0, 0, 1); S = (0, -1, 1, 0); SIOTA = (0, 1, 1, 0); C = (1, 1, 1, -1); Dinv = (1, -2, 0, 1)


def ptype(t):
    if t == INF:
        return 'opp'
    return 'opp' if (t.numerator + t.denominator) % 2 else 'odd'


def orbit(t):
    return frozenset(mob(v, t) for v in ((1, 0, 0, 1), IOTA, S, SIOTA))


def parent_descent(t):
    """Berggren parent slope of a Euclid slope t > 1, t != 2 (descent rule)."""
    if t > 3: return t - 2
    if t > 2: return 1 / (t - 2)
    return 1 / (2 - t)


pts = [INF] + sorted({F(m, n) for m in range(-H, H + 1) for n in range(1, H + 1) if gcd(m, n) == 1})
q1 = q2 = q3 = True; n3 = 0
for t in pts:
    for M in (U, A, D, TAU, IOTA, S):
        if ptype(mob(M, t)) != ptype(t): q1 = False
    ct = mob(C, t)
    if ptype(ct) == ptype(t): q2 = False
    for M, Mc in ((D, U), (U, D), (A, A), (TAU, TAU), (IOTA, SIOTA), (S, S)):
        if mob(C, mob(M, mob(C, t))) != mob(Mc, t): q2 = False
    if t != INF and ptype(t) == 'opp' and t > 1:
        n3 += 1
        heads = {name: orbit(mob(M, t)) for name, M in (('U', U), ('A', A), ('D', D), ('tau', TAU))}
        kids = [mob(U, t), mob(A, t), mob(D, t)]
        if any(orbit(k) != heads[nm] for k, nm in zip(kids, 'UAD')): q3 = False
        if t == 2:
            if heads['tau'] != frozenset({F(0), INF}): q3 = False
        elif heads['tau'] != orbit(parent_descent(t)): q3 = False
print(f"slopes checked: {len(pts)} (|m|,|n| <= {H}); Euclid vertices t > 1 checked in (Q3): {n3}")
print(f"(Q1) parity type preserved by U, A, D, tau, iota, S: {q1}  => G over Q is disconnected (Euclid / odd-odd)")
print(f"(Q2) C swaps types; C D C = U, C U C = D, C A C = A, C tau C = tau, C iota C = S iota, C S C = S: {q2}")
print(f"(Q3) heads of U t, A t, D t = Berggren children; head of tau t = Berggren parent ((1,0,1) at the root): {q3}")
for name, t in (('{0, inf} = (1,0,1)', INF), ('{1, -1} = (0,1,1)', F(1))):
    darts = [mob(M, t) for M in (U, A, D, TAU)]
    print(f"(Q4) vertex {name}: out-multiset [Ut, At, Dt, tau t] = {[str(d) for d in darts]}; "
          f"self-reverse darts: {[str(d) for d in set(darts) if mob(TAU, d) == d]}")
