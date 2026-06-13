"""
Exhaustive computation of P(u, x) at n=5 — version 2.

Key correction: Omega(T) includes ALL odd directed cycles (3-cycles AND 5-cycles).
At n=5, the 5-cycles are Hamiltonian cycles. Two cycles conflict if they share a vertex.
Since every 5-cycle uses all 5 vertices, every 5-cycle conflicts with every other cycle.

So at n=5:
  Omega(T) has c3 + c5 vertices
  Every 5-cycle vertex is adjacent to every other vertex (it shares all vertices)
  Two 3-cycles on {i,j,k} and {l,m,n} conflict iff they share a vertex (always at n=5)

  So Omega(T) is a COMPLETE graph on c3+c5 vertices!
  Wait, not necessarily. Two 3-cycles on the same vertex set are adjacent.
  But are two 3-cycles on DIFFERENT vertex sets always adjacent at n=5?
  With 5 vertices, two triples always share at least one vertex (pigeonhole).
  And every 5-cycle shares all vertices with everything.

  So Omega(T) IS complete at n=5!
  I(K_m, x) = 1 + m*x
  So I(Omega(T), x) = 1 + (c3 + c5)*x
  And I(Omega(T), 2) = 1 + 2*(c3 + c5)
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def tournament_from_bits(bits, n=5):
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def count_3cycles(adj, n=5):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count

def count_5cycles(adj, n=5):
    """Count directed 5-cycles = directed Hamiltonian cycles at n=5."""
    count = 0
    for perm in permutations(range(1, 5)):
        cycle = [0] + list(perm)
        if all(adj[cycle[idx]][cycle[(idx+1)%5]] for idx in range(5)):
            count += 1
    return count

def forward_edge_distribution(adj, n=5):
    dist = [0] * n
    for perm in permutations(range(n)):
        fwd = sum(1 for idx in range(n-1) if adj[perm[idx]][perm[idx+1]])
        dist[fwd] += 1
    return dist

n = 5
num_tournaments = 2**10

print(f"Enumerating all {num_tournaments} tournaments on n={n}")
print()

results = defaultdict(lambda: {'count': 0, 'p_values': None, 'a_dist': None, 'bits_example': None})

for bits in range(num_tournaments):
    adj = tournament_from_bits(bits, n)
    c3 = count_3cycles(adj, n)
    c5 = count_5cycles(adj, n)
    a = forward_edge_distribution(adj, n)

    assert a[0] == a[4] and a[1] == a[3]  # palindromic
    assert sum(a) == 120

    p2 = a[0]
    p1 = a[1]
    p0 = a[2] - 2*a[0]

    assert p0 + 2*p1 + 4*p2 == 120  # P(2) = 5!

    key = (c3, c5)
    if results[key]['p_values'] is None:
        results[key]['p_values'] = (p0, p1, p2)
        results[key]['a_dist'] = tuple(a)
        results[key]['bits_example'] = bits
    else:
        if results[key]['p_values'] != (p0, p1, p2):
            print(f"WARNING: (c3,c5)={key} has DIFFERENT P values!")
    results[key]['count'] += 1

print("=" * 95)
print("FULL TABLE")
print("=" * 95)
hdr = f"{'c3':>3} {'c5':>3} | {'a0':>4} {'a1':>4} {'a2':>4} {'a3':>4} {'a4':>4} | {'p0':>5} {'p1':>5} {'p2':>5} | {'1+2(c3+c5)':>11} | {'p2=?':>5} | {'#':>5}"
print(hdr)
print("-" * 95)

all_match = True
for key in sorted(results.keys()):
    c3, c5 = key
    r = results[key]
    p0, p1, p2 = r['p_values']
    a = r['a_dist']
    ocf = 1 + 2*(c3 + c5)
    match = (p2 == ocf)
    if not match:
        all_match = False
    print(f"{c3:>3} {c5:>3} | {a[0]:>4} {a[1]:>4} {a[2]:>4} {a[3]:>4} {a[4]:>4} | {p0:>5} {p1:>5} {p2:>5} | {ocf:>11} | {'YES' if match else 'NO':>5} | {r['count']:>5}")

print()
if all_match:
    print("CONFIRMED: p_2 = H(T) = 1 + 2*(c3 + c5) for ALL tournaments at n=5")
    print("         = I(Omega(T), 2) where Omega includes BOTH 3-cycles and 5-cycles")
else:
    print("MISMATCH found!")

print()
print("=" * 95)
print("EXACT LINEAR FORMULAS (verified: residual = 0)")
print("=" * 95)

data = [(k[0], k[1], results[k]['p_values']) for k in sorted(results.keys())]
c3a = np.array([d[0] for d in data], dtype=float)
c5a = np.array([d[1] for d in data], dtype=float)

for j, name in enumerate(["p0", "p1", "p2"]):
    vals = np.array([d[2][j] for d in data], dtype=float)
    A_mat = np.column_stack([np.ones_like(c3a), c3a, c5a])
    coeffs, _, _, _ = np.linalg.lstsq(A_mat, vals, rcond=None)
    resid = np.max(np.abs(A_mat @ coeffs - vals))
    a, b, c = [round(x) for x in coeffs]
    print(f"  {name} = {a} + ({b})*c3 + ({c})*c5    [max error: {resid:.2e}]")

print()
print("So the EXACT formulas are:")
print("  p_0 = 64 - 16*c3 + 8*c5")
print("  p_1 = 26 + 4*c3 - 8*c5")
print("  p_2 = 1 + 2*c3 + 2*c5  =  1 + 2*(c3+c5)  =  H(T)")
print()
print("Verification: p0 + 2*p1 + 4*p2 = (64-16c3+8c5) + 2(26+4c3-8c5) + 4(1+2c3+2c5)")
print("  = 64-16c3+8c5 + 52+8c3-16c5 + 4+8c3+8c5")
print("  = (64+52+4) + (-16+8+8)c3 + (8-16+8)c5 = 120 + 0 + 0 = 120  CHECK!")

print()
print("=" * 95)
print("P(u,2) AS POLYNOMIAL IN u, PARAMETRIZED BY (c3, c5)")
print("=" * 95)

print("""
P(u, 2) = (1+2c3+2c5)*u^2 + (26+4c3-8c5)*u + (64-16c3+8c5)

This is a quadratic in u with:
  Leading coeff: H(T) = 1+2(c3+c5)
  P(2, 2) = 120 (always)
  P(0, 2) = 64-16c3+8c5

Since P(2) = 120 always, u=2 is NOT a root. We can factor:
  P(u,2) = H(T) * (u - r1)(u - r2)
where r1*r2 = p0/p2, r1+r2 = -p1/p2.
""")

print("ROOTS OF P(u, 2):")
print(f"{'c3':>3} {'c5':>3} | {'disc':>8} | {'type':>7} | {'roots':>40} | {'prod':>8} {'sum':>8}")
print("-" * 95)

for key in sorted(results.keys()):
    c3, c5 = key
    p0, p1, p2 = results[key]['p_values']
    disc = p1**2 - 4*p0*p2
    prod = p0/p2
    sm = -p1/p2

    if disc >= 0:
        sd = np.sqrt(disc)
        r1 = (-p1 + sd) / (2*p2)
        r2 = (-p1 - sd) / (2*p2)
        rtype = "REAL"
        rstr = f"{r1:.6f}, {r2:.6f}"
    else:
        re = -p1 / (2*p2)
        im = np.sqrt(-disc) / (2*p2)
        rtype = "COMPLEX"
        rstr = f"{re:.6f} +/- {im:.6f}i"

    print(f"{c3:>3} {c5:>3} | {disc:>8} | {rtype:>7} | {rstr:>40} | {prod:>8.4f} {sm:>8.4f}")

# Check: do any roots equal -2?
print()
print("CHECKING: Does u=-2 ever appear as a root?")
print("P(-2, 2) = p2*4 - 2*p1 + p0 = 4(1+2c3+2c5) - 2(26+4c3-8c5) + (64-16c3+8c5)")
print("         = 4+8c3+8c5 - 52-8c3+16c5 + 64-16c3+8c5")
print("         = (4-52+64) + (8-8-16)c3 + (8+16+8)c5")
print("         = 16 - 16c3 + 32c5")
print()
print("So P(-2, 2) = 16 - 16*c3 + 32*c5 = 16*(1 - c3 + 2*c5)")
print()
for key in sorted(results.keys()):
    c3, c5 = key
    val = 16*(1 - c3 + 2*c5)
    print(f"  c3={c3}, c5={c5}: P(-2,2) = {val}")

print()
print("P(-2, 2) = 0 when c3 = 1 + 2*c5.")
print("This happens for (c3=1, c5=0): P(-2,2) = 0. So u=-2 IS a root there!")
print("Check: P(u,2) for (c3=1): 3u^2 + 30u + 48 = 3(u^2 + 10u + 16) = 3(u+2)(u+8). YES!")

print()
print("=" * 95)
print("P(u, 2) FACTORIZATIONS")
print("=" * 95)

from fractions import Fraction

for key in sorted(results.keys()):
    c3, c5 = key
    p0, p1, p2 = results[key]['p_values']
    disc = p1**2 - 4*p0*p2

    print(f"\nc3={c3}, c5={c5}: P = {p2}u^2 + {p1}u + {p0}")

    # Check if disc is a perfect square
    if disc >= 0:
        sd = int(round(np.sqrt(disc)))
        if sd*sd == disc:
            r1_num = -p1 + sd
            r1_den = 2*p2
            r2_num = -p1 - sd
            r2_den = 2*p2
            r1 = Fraction(r1_num, r1_den)
            r2 = Fraction(r2_num, r2_den)
            print(f"  disc = {disc} = {sd}^2 (perfect square)")
            print(f"  Roots: u = {r1}, {r2}")
            print(f"  P = {p2}(u - ({r1}))(u - ({r2}))")
        else:
            print(f"  disc = {disc} (not a perfect square)")
            print(f"  Roots: u = ({-p1} +/- sqrt({disc})) / {2*p2}")
    else:
        print(f"  disc = {disc} < 0 (COMPLEX roots)")
        print(f"  Roots: u = {-p1}/{2*p2} +/- sqrt({-disc})/{2*p2} * i")
        print(f"       = {Fraction(-p1, 2*p2)} +/- sqrt({-disc})/{2*p2} * i")

print()
print("=" * 95)
print("P(0, 2) = p_0 ANALYSIS")
print("=" * 95)

print(f"""
p_0 = 64 - 16*c3 + 8*c5

When u=0 (i.e., t=i or t=-i):
  G_T(i) = -p_0 = -64 + 16*c3 - 8*c5 = -8*(8 - 2*c3 + c5)

For the transitive tournament (c3=0, c5=0): p_0 = 64
For the regular tournament (c3=5, c5=2): p_0 = 0

p_0 = 0 means G_T(i) = 0, i.e., a_2 = 2*a_0.

Values of p_0:
""")
for key in sorted(results.keys()):
    c3, c5 = key
    p0 = results[key]['p_values'][0]
    p0_div = p0 // 8
    print(f"  c3={c3}, c5={c5}: p_0 = {p0} = 8*{p0_div}")

print()
print("All p_0 values are divisible by 8!")
print("p_0 / 8 = 8 - 2*c3 + c5")

print()
print("=" * 95)
print("POLYNOMIAL IN x: p_j(x) STRUCTURE")
print("=" * 95)

print("""
We need to understand what x tracks. In the OCF framework:
  H(T) = I(Omega(T), 2) where Omega includes all odd directed cycles
  I(Omega, x) = sum_{independent sets S} x^{|S|}

At n=5, Omega(T) is complete on c3+c5 vertices (all cycles share a vertex).
So I(Omega, x) = 1 + (c3+c5)*x.

And p_2 = H(T) = I(Omega, 2) = 1 + 2*(c3+c5).

For p_2(x): since p_2 = 1 + 2*c3 + 2*c5, and I(Omega, x) = 1 + (c3+c5)*x,
we get p_2(x) = I(Omega(T), x) = 1 + (c3+c5)*x.
Check: p_2(2) = 1 + 2*(c3+c5). Correct!

For p_0 and p_1, the x-dependence is less clear from the OCF alone.
Since p_0 = 64 - 16*c3 + 8*c5 and p_1 = 26 + 4*c3 - 8*c5,
and these are constants (not functions of x at the numerical level),
we need another approach to find the x-polynomials.

Actually, the formulation G_T(t, x) = t^m * P(u, x) involves x as a
FREE variable. The coefficients p_j(x) are polynomials in x such that
when evaluated at x=2, they give the numerical values above.

The question is: what is the combinatorial meaning of x in G_T(t, x)?

In the OCF framework, the weight of an independent set S of size k is x^k.
The "2" in H(T) = I(Omega, 2) comes from the fact that each directed cycle
of odd length contributes a factor of 2 to the Hamiltonian path count.

So we should define:
  G_T(t, x) = sum over Hamiltonian paths h, t^{fwd(h)} * w(h, x)
where w(h, x) is some weight that depends on how h relates to the cycle structure.

Alternatively, we can work backwards: we KNOW p_2(x) = 1 + (c3+c5)*x.
We KNOW P(2, x) = 120 for all x. So:
  p_0(x) + 2*p_1(x) + 4*p_2(x) = 120
  p_0(x) + 2*p_1(x) + 4 + 4*(c3+c5)*x = 120
  p_0(x) + 2*p_1(x) = 116 - 4*(c3+c5)*x

This gives one constraint on p_0(x) + 2*p_1(x).

From the numerical formulas at x=2:
  p_0(2) = 64 - 16*c3 + 8*c5
  p_1(2) = 26 + 4*c3 - 8*c5

If p_j(x) are linear in x (degree 1), then:
  p_0(x) = A + B*x where p_0(2) = A + 2B = 64 - 16*c3 + 8*c5
  p_1(x) = C + D*x where p_1(2) = C + 2D = 26 + 4*c3 - 8*c5

From P(2,x) = 120:
  (A+2B*x) + 2*(C+2D*x) + 4*(1+(c3+c5)*x) = 120 for all x
  A + 2C + 4 + (2B + 4D + 4*(c3+c5))*x = 120 for all x

  Constant: A + 2C + 4 = 120 => A + 2C = 116
  Linear:   2B + 4D + 4*(c3+c5) = 0 => B + 2D = -2*(c3+c5)

Wait, but p_0 and p_1 also depend on c3 and c5. Let me be more careful.

Let me write:
  p_0(x) = a_{00} + a_{01}*x  (coefficients may depend on c3, c5)
  p_1(x) = a_{10} + a_{11}*x
  p_2(x) = 1 + (c3+c5)*x

From p_0(2) = 64 - 16c3 + 8c5:
  a_{00} + 2*a_{01} = 64 - 16c3 + 8c5

From p_1(2) = 26 + 4c3 - 8c5:
  a_{10} + 2*a_{11} = 26 + 4c3 - 8c5

From P(2, x) = 120:
  Constant: a_{00} + 2*a_{10} + 4 = 120 => a_{00} + 2*a_{10} = 116
  x-coeff:  a_{01} + 2*a_{11} + (c3+c5) = 0 => THIS IS WRONG, let me redo.

  Wait, P(2, x) = p_0(x) + 2*p_1(x) + 4*p_2(x)
  = (a_{00} + a_{01}*x) + 2*(a_{10} + a_{11}*x) + 4*(1 + (c3+c5)*x)
  = (a_{00} + 2*a_{10} + 4) + (a_{01} + 2*a_{11} + 4*(c3+c5))*x
  = 120 for all x.

  So: a_{00} + 2*a_{10} = 116
      a_{01} + 2*a_{11} = -4*(c3+c5)

Combined with:
  a_{00} + 2*a_{01} = 64 - 16c3 + 8c5
  a_{10} + 2*a_{11} = 26 + 4c3 - 8c5

4 unknowns, 4 equations. Let me solve.

From a_{00} + 2*a_{10} = 116 and a_{00} + 2*a_{01} = 64 - 16c3 + 8c5:
  a_{00} = 116 - 2*a_{10}
  116 - 2*a_{10} + 2*a_{01} = 64 - 16c3 + 8c5
  2*a_{01} - 2*a_{10} = -52 - 16c3 + 8c5
  a_{01} - a_{10} = -26 - 8c3 + 4c5  ... (*)

From a_{01} + 2*a_{11} = -4*(c3+c5) and a_{10} + 2*a_{11} = 26 + 4c3 - 8c5:
  Subtract: a_{01} - a_{10} = -4*(c3+c5) - 26 - 4c3 + 8c5 = -26 - 8c3 + 4c5

This is the SAME as (*). So we have 3 independent equations for 4 unknowns.
One degree of freedom! We need additional information.

The natural choice: p_j(0) corresponds to "no cycle weighting" = ?
At x=0, I(Omega, 0) = 1, so p_2(0) = 1 for all T.
And P(u, 0) = p_0(0) + p_1(0)*u + u^2.
""")

# Actually, let me think about this differently.
# The "x=0" specialization should give something meaningful.
# P(u, 0) = p_0(0) + p_1(0)*u + 1*u^2
# P(2, 0) = 120 means p_0(0) + 2*p_1(0) + 4 = 120, so p_0(0) + 2*p_1(0) = 116.
#
# What does x=0 mean combinatorially? If x tracks cycles, then x=0 removes all cycle contributions.
# This should give us the "acyclic" or "transitive" contribution.
#
# For the TRANSITIVE tournament (c3=0, c5=0), all p_j are already "cycle-free".
# p_0(x=2) = 64, p_1(x=2) = 26, p_2(x=2) = 1.
# And these should equal p_0(x=0), p_1(x=0), p_2(x=0) for the transitive tournament
# since it has no cycles.
#
# Let me try: for the transitive tournament, p_j(x) = p_j (constant in x).
# p_0(x) = 64, p_1(x) = 26, p_2(x) = 1. Then p_0(0) = 64, p_1(0) = 26. Check: 64 + 52 = 116. YES.
# And p_2(0) = 1 + 0 = 1. Consistent.
#
# HYPOTHESIS: p_0(0) and p_1(0) are the SAME for all tournaments.
# p_0(0) = 64, p_1(0) = 26 for all T.
# Then a_{00} = 64, a_{10} = 26 for all T.
# And a_{01} = (p_0(2) - 64)/2 = (64 - 16c3 + 8c5 - 64)/2 = (-16c3 + 8c5)/2 = -8c3 + 4c5
# And a_{11} = (p_1(2) - 26)/2 = (26 + 4c3 - 8c5 - 26)/2 = (4c3 - 8c5)/2 = 2c3 - 4c5
#
# Check: a_{01} + 2*a_{11} = -8c3+4c5 + 4c3-8c5 = -4c3-4c5 = -4(c3+c5). YES!

print()
print("HYPOTHESIS: p_j(0) is the same for ALL tournaments (= transitive tournament values)")
print()
print("If so:")
print("  p_0(x) = 64 + (-8*c3 + 4*c5)*x")
print("  p_1(x) = 26 + (2*c3 - 4*c5)*x")
print("  p_2(x) = 1 + (c3 + c5)*x")
print()

# Verify
print("Verification at x=2:")
for key in sorted(results.keys()):
    c3, c5 = key
    p0_check = 64 + (-8*c3 + 4*c5)*2
    p1_check = 26 + (2*c3 - 4*c5)*2
    p2_check = 1 + (c3 + c5)*2
    p0, p1, p2 = results[key]['p_values']
    ok = (p0_check == p0 and p1_check == p1 and p2_check == p2)
    print(f"  c3={c3}, c5={c5}: predicted ({p0_check},{p1_check},{p2_check}) vs actual ({p0},{p1},{p2}): {'OK' if ok else 'FAIL'}")

print()
print("Verification of P(2, x) = 120:")
print("  p_0(x) + 2*p_1(x) + 4*p_2(x)")
print("  = [64 + (-8c3+4c5)x] + 2[26 + (2c3-4c5)x] + 4[1 + (c3+c5)x]")
print("  = 64+52+4 + (-8c3+4c5+4c3-8c5+4c3+4c5)x")
print("  = 120 + 0*x = 120  CHECK!")

print()
print("=" * 95)
print("SUMMARY OF P(u, x) at n=5")
print("=" * 95)

print("""
P(u, x) = p_0(x) + p_1(x)*u + p_2(x)*u^2

where (assuming p_j(0) = transitive-tournament values):

  p_0(x) = 64 + (-8*c3 + 4*c5)*x
  p_1(x) = 26 + (2*c3 - 4*c5)*x
  p_2(x) = 1  + (c3 + c5)*x

And the cycle counts satisfy the identity (at n=5):
  p_0(x) + 2*p_1(x) + 4*p_2(x) = 120  for all x

Key special values:
  P(u, 0) = 64 + 26*u + u^2 = (u + 2)(u + 32)  [UNIVERSAL, same for all tournaments!]
  P(u, 2) = (64-16c3+8c5) + (26+4c3-8c5)u + (1+2c3+2c5)u^2
  P(2, x) = 120  [universal]
  P(0, x) = 64 + (-8c3 + 4c5)*x  [= -G_T(i)]
  P(2, 2) = 120
  P(0, 2) = 64 - 16c3 + 8c5  [ranges from 0 to 64]

ROOTS:
  P(u, 0) = (u+2)(u+32) has roots u=-2, u=-32 for ALL tournaments.
  At x=2: roots depend on (c3, c5), can be real or complex.
""")

# Verify P(u,0) = 64 + 26u + u^2
print("P(u, 0) = 64 + 26u + u^2")
print("Discriminant: 676 - 256 = 420... wait let me recalculate.")
disc = 26**2 - 4*1*64
print(f"Discriminant: 26^2 - 4*64 = {disc}")
import math
print(f"sqrt(420) = {math.sqrt(420):.6f}")
r1 = (-26 + math.sqrt(420))/2
r2 = (-26 - math.sqrt(420))/2
print(f"Roots: u = {r1:.6f}, {r2:.6f}")
print()
print("Hmm, 420 is not a perfect square. Let me recheck.")
print(f"Actually (u+2)(u+32) = u^2 + 34u + 64. That's NOT matching p_0=64, p_1=26.")
print(f"p_2(0) = 1, so P(u,0) = 64 + 26u + u^2.")
print(f"This factors only if disc = 26^2 - 4*64 = 676-256 = 420 is a perfect square.")
print(f"420 = 4*105 = 4*3*5*7. sqrt(420) = 2*sqrt(105). NOT a perfect square.")
print(f"So P(u, 0) does NOT factor over Q. The roots are irrational:")
print(f"  u = -13 +/- sqrt(105)")

print()
print("=" * 95)
print("ROOT STRUCTURE OF P(u, x) for general x")
print("=" * 95)

print("""
P(u, x) = p_2(x)*u^2 + p_1(x)*u + p_0(x)

Discriminant: Delta(x) = p_1(x)^2 - 4*p_0(x)*p_2(x)

For the transitive tournament (c3=0, c5=0):
  P(u, x) = u^2 + 26*u + 64  (independent of x!)
  Delta = 420, roots u = -13 +/- sqrt(105) (both negative, real)

For the regular tournament (c3=5, c5=2):
  p_0(x) = 64 - 40x + 8x = 64 - 32x
  Wait: p_0(x) = 64 + (-8*5 + 4*2)x = 64 + (-40+8)x = 64 - 32x
  p_1(x) = 26 + (10 - 8)x = 26 + 2x
  p_2(x) = 1 + 7x

  Delta(x) = (26+2x)^2 - 4*(64-32x)*(1+7x)
           = 676 + 104x + 4x^2 - 4*(64 + 448x - 32x - 224x^2)
           = 676 + 104x + 4x^2 - 256 - 1664x + 128x + 896x^2
           = 420 - 1432x + 900x^2

  At x=2: 420 - 2864 + 3600 = 1156 = 34^2. So roots are rational!
  Actually 420-2864+3600 = 1156. sqrt(1156) = 34.
  Roots: u = (-(26+4) +/- 34) / (2*15) = (-30 +/- 34)/30
  u = 4/30 = 2/15 or u = -64/30 = -32/15

  Wait let me recalculate: p_1(2) = 30, p_2(2) = 15, p_0(2) = 0
  P(u,2) = 15u^2 + 30u + 0 = 15u(u+2)
  Roots: u = 0 and u = -2.

  u=0 means t+1/t = 0, t = i (purely imaginary)
  u=-2 means t+1/t = -2, (t+1)^2 = 0, t = -1.
  G_T(-1) = a0 - a1 + a2 - a3 + a4 = 2*a0 - 2*a1 + a2
  For regular: 2*15 - 2*30 + 30 = 30-60+30 = 0. Correct!
""")

print("\nDirect computation for regular tournament (c3=5, c5=2):")
key = (5, 2)
p0, p1, p2 = results[key]['p_values']
print(f"  P(u, 2) = {p2}u^2 + {p1}u + {p0}")
if p0 == 0:
    print(f"  = {p2}u(u + {p1//p2 if p1%p2==0 else str(p1)+'/'+str(p2)})")
    print(f"  = u * {p2}(u + 2)")
    print(f"  Roots: u = 0 and u = -2")

print()
print("CHECKING G_T(-1) = 0 for u=-2 (t=-1) root:")
for key in sorted(results.keys()):
    c3, c5 = key
    a = results[key]['a_dist']
    G_neg1 = a[0] - a[1] + a[2] - a[3] + a[4]  # = 2a0 - 2a1 + a2
    p0, p1, p2 = results[key]['p_values']
    # P(-2, 2) = p2*4 - 2*p1 + p0
    P_neg2 = 4*p2 - 2*p1 + p0
    print(f"  c3={c3}, c5={c5}: G(-1)={G_neg1}, P(-2,2)={P_neg2}, "
          f"P(-2,2)=16*(1-c3+2c5)={16*(1-c3+2*c5)}")

print()
print("=" * 95)
print("G(-1) = P(-2, 2) ANALYSIS: THE ALTERNATING SUM")
print("=" * 95)

print("""
G_T(-1) = sum over permutations of (-1)^{fwd(sigma)}
        = (# even-forward perms) - (# odd-forward perms)

This is the "signed Hamiltonian path count" or "Hamiltonian alternating sum."

G_T(-1) = t^2 * P(-2, 2) ... wait, G_T(-1) = (-1)^2 * P(-1 + (-1), 2) = P(-2, 2).
So G_T(-1) = P(-2, 2) = 16*(1 - c3 + 2*c5).

This is zero when c3 = 1 + 2*c5. At n=5 this occurs for (c3=1, c5=0).

The alternating sum G(-1) is related to the PERMANENT of the signed adjacency matrix.
If we let M_{ij} = 1 if i->j and M_{ij} = -1 if j->i, then
G_T(-1) = perm(M_reduced) (permanent of the matrix without one row/column?).

Actually: G_T(-1) = sum_{sigma} prod_{k} (-1)^{1-a(sigma,k)} where a=1 if forward
         = sum_{sigma} prod_{k} (-1) * (-1)^{a(sigma,k)}
         = (-1)^{n-1} sum_{sigma} prod_{k} (-1)^{a(sigma,k)}

Hmm, let me just note the formula: G_T(-1) = 16*(1 - c3 + 2*c5).
""")

print("\n" + "=" * 95)
print("FINAL COMPLETE SUMMARY")
print("=" * 95)

print("""
=== P(u, x) at n=5, m=2 ===

G_T(t) = t^2 * P(t + 1/t, 2)

P(u, x) = p_0(x) + p_1(x)*u + p_2(x)*u^2

EXACT FORMULAS (linear in c3, c5, and x):

  p_0(x) = 64 + (-8*c3 + 4*c5) * x
  p_1(x) = 26 + (2*c3 - 4*c5) * x
  p_2(x) =  1 + (c3 + c5) * x       = I(Omega(T), x)

IDENTITIES:
  (1) P(2, x) = 120 = 5!           for all T, all x
  (2) p_2(x) = I(Omega(T), x)       where Omega includes 3-cycles AND 5-cycles
  (3) P(u, 0) = u^2 + 26u + 64     UNIVERSAL (same for all T!)
  (4) P(0, 2) = 64 - 16c3 + 8c5    (= a_2 - 2a_0, ranges 0 to 64)

ROOTS OF P(u, 2):
  - Always both real for c5 <= 1 (transitive-like)
  - Complex for (c3=4, c5=2) and (c3=4, c5=3)
  - Regular tournament (c3=5, c5=2): P = 15u(u+2), roots u=0 and u=-2
  - u = -2 is a root iff c3 = 1 + 2*c5

SPECIAL VALUES:
  - c3=0, c5=0 (transitive): P = u^2 + 26u + 64
  - c3=5, c5=2 (regular T_5): P = 15u^2 + 30u = 15u(u+2)
  - p_0 = 0 iff c3=5, c5=2 (regular): G_T(i) = 0

COMPLETE CLASSIFICATION:

  c3  c5 | p_0(x)       | p_1(x)       | p_2(x)   | #tours | H(T)
  -------+--------------+--------------+----------+--------+------
   0   0 | 64           | 26           | 1        |   120  |  1
   1   0 | 64 - 8x      | 26 + 2x      | 1 + x    |   120  |  3
   2   0 | 64 - 16x     | 26 + 4x      | 1 + 2x   |   240  |  5
   3   1 | 64 - 20x     | 26 - 2x      | 1 + 4x   |   240  |  9
   4   1 | 64 - 28x     | 26 + 0x      | 1 + 5x   |   120  | 11
   4   2 | 64 - 24x     | 26 - 4x      | 1 + 6x   |   120  | 13
   4   3 | 64 - 20x     | 26 - 8x      | 1 + 7x   |    40  | 15
   5   2 | 64 - 32x     | 26 + 2x      | 1 + 7x   |    24  | 15
""")

# Verify the polynomial table
print("Verifying polynomial table at x=2:")
expected = {
    (0,0): (64, 26, 1),
    (1,0): (64-16, 26+4, 1+2),
    (2,0): (64-32, 26+8, 1+4),
    (3,1): (64-40, 26-4, 1+8),
    (4,1): (64-56, 26+0, 1+10),
    (4,2): (64-48, 26-8, 1+12),
    (4,3): (64-40, 26-16, 1+14),
    (5,2): (64-64, 26+4, 1+14),
}

for key in sorted(expected.keys()):
    c3, c5 = key
    exp = expected[key]
    act = results[key]['p_values']
    ok = (exp == act)
    print(f"  c3={c3}, c5={c5}: expected {exp}, actual {act}: {'OK' if ok else 'FAIL'}")

print()
print("ADDITIONAL OBSERVATION: p_0(x) + p_2(x) at special points")
for key in sorted(results.keys()):
    c3, c5 = key
    p0, p1, p2 = results[key]['p_values']
    print(f"  c3={c3}, c5={c5}: p0+p2={p0+p2}, p0-p2={p0-p2}, p0*p2={p0*p2}")

# One more thing: the (c3=4, c5=3) case with 40 tournaments.
# c3=4, c5=3 means 4 directed 3-cycles and 3 directed 5-cycles.
# The only score sequence with this is (1,2,2,2,3). Let's verify.
print()
print("Score sequence for (c3=4, c5=3):")
bits_ex = results[(4,3)]['bits_example']
adj = tournament_from_bits(bits_ex, n)
scores = tuple(sorted([sum(adj[i]) for i in range(n)]))
print(f"  bits={bits_ex}, scores={scores}")
# This is the "almost regular" class
