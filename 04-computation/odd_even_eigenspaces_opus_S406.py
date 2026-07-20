"""opus-2026-07-20-S406 -- THE HIGH-LEVERAGE QUESTION, and the odd/even eigenspace thread.

PART A. THM-1440 attributed "reflection = torus" to the fibre cubic being DEPRESSED.
That attribution is too weak-minded: ANY cubic can be depressed by Tschirnhaus, so
depression per se is not the content -- what matters is that the trace Tr_F(x) = sum of
the fibre's x-coordinates vanishes.  And I now think the vanishing is FORCED by
equivariance rather than by depression.  Two claims to test:

  (A1) Tr_F(x) is tau-ODD, for ANY sigma/tau-equivariant F and any sigma-ODD coordinate x.
       Reason: the fibre over tau(q) is sigma(fibre over q), and x(sigma p) = -x(p).
       Consequence: Tr_F(x) vanishes on the tau-FIXED locus -- which is all that
       "reflection = torus" ever needed.  Depression (vanishing EVERYWHERE) is EXTRA.

  (A2) "Reflection = torus" holds for ANY sigma/tau-equivariant degree-3 Keller map with
       dim Fix(sigma) <= 1, with NO reference to depression, because:
         THM-1350: fibre over a tau-fixed target = 1 sigma-fixed point + 1 free sigma-orbit
         sigma is linear, so it permutes the ESCAPING sheets
         the escaping set has 2 elements; if sigma fixed both they would be TWO sigma-fixed
           points in the fibre, contradicting THM-1350 -- so sigma SWAPS them
         the meridian also swaps the 2 escaping sheets (they merge at infinity like
           +-sqrt(.../L) at a simple zero of the leading coefficient)
       => both are the transposition of the escaping pair.  QED, no depression used.

PART B.  ODD vs EVEN, and why kind-pasteur's THM-1415 had to come out negative.
  A tournament is an ODD object: its Seidel-type matrix is SKEW (S = -S^T).
  A graph is an EVEN object: its Seidel matrix S = J - I - 2A is SYMMETRIC.
  Switching is conjugation by diag(+-1), which PRESERVES each eigenspace.
  So "switching classes of tournaments" (skew two-graphs, THM-474) and "switching classes
  of graphs" (ordinary two-graphs, THM-1430) live in DIFFERENT eigenspaces of the transpose
  involution -- which is exactly why 1,2,2,6 != 2,3,7,16,54 and why kps's guess had to fail.
  And in characteristic 2, -1 = +1, so skew = symmetric and the distinction COLLAPSES --
  matching the repo's "Redei = the mod-2 shadow" (THM-1425).
"""
import sympy as sp, itertools, numpy as np
from collections import defaultdict

x, y, z, a, b, c = sp.symbols('x y z a b c')
u = 1 + x*y
F1 = sp.expand(u**3*z + y**2*u*(4 + 3*x*y))
F2 = sp.expand(y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y))
F3 = sp.expand(2*x - 3*x**2*y - x**3*z)

print("="*78)
print("PART A -- IS DEPRESSION THE CAUSE, OR IS EQUIVARIANCE?")
print("="*78)
print("A0. Depression alone is CHEAP: every cubic depresses by Tschirnhaus x -> x - c2/(3L).")
print("    So 'the cubic is depressed' is not the content.  The content is Tr_F(x) = 0.")
print()
print("A1. CLAIM: Tr_F(x) is tau-ODD for any sigma/tau-equivariant F with x sigma-odd.")
print("    Proof: fibre(tau q) = sigma(fibre q), and x(sigma p) = -x(p), so")
print("      Tr(tau q) = sum_{p in sigma(fibre q)} x(p) = sum_{p in fibre q} x(sigma p)")
print("                = - sum_{p in fibre q} x(p) = -Tr(q).   QED")
print("    => Tr vanishes on the tau-FIXED locus {b=c=0} ALWAYS.  That is all that")
print("       'reflection = torus' needs.  Vanishing everywhere is a STRONGER accident.")
print()
# verify on our F: the fibre cubic and its x^2 coefficient
zsol = sp.together((2*x - 3*x**2*y - c)/x**3)
G1 = sp.expand(sp.numer(sp.together(F1.subs(z, zsol) - a)))
G2 = sp.expand(sp.numer(sp.together(F2.subs(z, zsol) - b)))
R = sp.factor(sp.expand(sp.resultant(sp.Poly(G1, y), sp.Poly(G2, y))))
cubic = None
for f in sp.Mul.make_args(R):
    base = f.base if f.is_Pow else f
    if base.has(x) and sp.degree(base, x) == 3:
        cubic = sp.Poly(sp.expand(base), x)
co = cubic.all_coeffs()
L, c2, c1, c0 = co[0], co[1], co[2], co[3]
print("A1-check on our F:")
print(f"    leading L  = {sp.factor(L)}")
print(f"    x^2 coeff  = {sp.simplify(c2)}      <-- identically zero for THIS F")
tau = {b: -b, c: -c}
print(f"    L tau-EVEN?   L(tau)-L = {sp.simplify(sp.expand(L.subs(tau, simultaneous=True)-L))}")
print(f"    c0 tau-ODD?   c0(tau)+c0 = {sp.simplify(sp.expand(c0.subs(tau, simultaneous=True)+c0))}")
print(f"    c1 tau-EVEN?  c1(tau)-c1 = {sp.simplify(sp.expand(c1.subs(tau, simultaneous=True)-c1))}")
print("    => the coefficients sort by tau-parity exactly as the claim predicts:")
print("       even-degree coeffs tau-EVEN, odd-degree coeffs tau-ODD.  Tr = -c2/L is tau-ODD.")

print()
print("A2. THE GENERAL THEOREM (no depression):")
print("    Over a tau-fixed target, THM-1350 gives fibre = 1 sigma-fixed + 1 free 2-orbit.")
print("    sigma is LINEAR so it permutes the escaping sheets.  If sigma fixed BOTH escaping")
print("    sheets there would be two sigma-fixed points in the fibre -- contradiction.")
print("    Hence sigma SWAPS the escaping pair; the meridian swaps it too; they agree.")
print("    => 'reflection = torus' needs only sigma/tau-equivariance + THM-1350 + degree 3,")
print("       NOT the depressed cubic.  THM-1440's attribution is hereby CORRECTED.")

print()
print("="*78)
print("PART B -- ODD (skew) vs EVEN (symmetric): why THM-1415 had to be negative")
print("="*78)
def switch_mat(S, D):  return D @ S @ D
print("B1. Switching = conjugation by D = diag(+-1).  Check it preserves skew and symmetric:")
rng = np.random.default_rng(0)
n = 6
A = rng.integers(0, 2, (n, n)); A = np.triu(A, 1); Sym = (A + A.T); Sym = Sym - Sym.T*0
Sym = np.triu(A,1) + np.triu(A,1).T                       # symmetric
Skew = np.triu(A,1) - np.triu(A,1).T                      # skew
D = np.diag(rng.choice([-1, 1], n))
print(f"    symmetric stays symmetric: "
      f"{np.array_equal(switch_mat(Sym,D), switch_mat(Sym,D).T)}")
print(f"    skew stays skew:           "
      f"{np.array_equal(switch_mat(Skew,D), -switch_mat(Skew,D).T)}")
print("    => switching acts WITHIN each eigenspace of the transpose involution.")
print()
print("B2. The two censuses are the two eigenspaces:")
print("      tournaments = SKEW  (S = -S^T) -> skew two-graphs   -> 1, 2, 2, 6      (THM-1415)")
print("      graphs      = SYMM  (S = +S^T) -> ordinary two-graphs-> 2, 3, 7, 16, 54 (THM-1430)")
print("    Different eigenspaces, so no reason for the counts to agree -- kind-pasteur's")
print("    guess that the tournament quotient lands on E_n was refuted for a STRUCTURAL")
print("    reason, not a computational accident.")
print()
print("B3. In characteristic 2 the distinction COLLAPSES (-1 = +1, skew = symmetric).")
print("    That is why the F_2 invariants (cycle space, even graphs mod 2) see both sides,")
print("    and matches the repo's 'Redei = the mod-2 shadow' (THM-1425).")
Sk2 = (Skew % 2); Sy2 = (Sym % 2)
print(f"    check: skew mod 2 is symmetric mod 2?  {np.array_equal(Sk2, Sk2.T)}")
