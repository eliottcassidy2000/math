#!/usr/bin/env python3
"""collatz_procgen_20260924_fixedpt_kawasaki.py

Lane: fixed points / the owner's inequality / arXiv:2502.20642v1
(session collatz-procgen-20260922, 2026-09-24).

Audit of T. Kawasaki, "A proof of the Collatz conjecture", arXiv:2502.20642v1 [math.GM].
The paper's objects are transcribed below from the arXiv HTML (read in full by this lane).

Notation.  X metric space, T: X -> X, coefficient functions (al, be, ga, de, ep, ze) on X x X.
  weighted generalized pseudocontraction (WGP):
     al d(Tx,Ty)^2 + be d(x,Ty)^2 + ga d(Tx,y)^2 + de d(x,y)^2 + ep d(x,Tx)^2 + ze d(y,Ty)^2 <= 0.
  alt1(x,y): al+ze+2min(be,0) > 0,  -(de+ep+2min(be,0)) <= A (al+ze+2min(be,0))   [+ al+be+ze >= B in Thm 2.3]
  alt2(x,y): al+ep+2min(ga,0) > 0,  -(de+ze+2min(ga,0)) <= A (al+ep+2min(ga,0))   [+ al+ga+ep >= B in Thm 2.3]
  condition (5) (pairwise reading, the one Section 3 uses): for every pair, alt1 or alt2; |coeffs| <= M.
  Section 3: T(1) = 1, T(x) = x/2 (x even), (3x+1)/2 (x odd >= 3) on N = {1,2,...}, d = |x-y|,
  lambda = 0, A = 1/2, B = 2, M = 2, with the piecewise table of Theorem 3.1.

Sections (printed):
  K1  Theorem 3.1: the WGP inequality for the paper's T, all 1 <= x,y < 400 (exact), and the paper's
      nine case polynomials re-derived symbolically (one harmless algebra slip found)
  K2  condition (5) of Theorem 2.3 at every pair < 400; census of alternatives; neither global reading holds
  K3  the swap gap: the proof of Thm 2.1(5) needs alt1(p,Tp) OR alt2(Tp,p); census p < 5000
  K4  the paper's intermediate claim d(T^n x,T^(n+1) x)^2 <= A^n d(x,Tx)^2 and its orbit-radius bound
  K5  stretch of consecutive distances at up-steps; the Mersenne family
  K6  counterexamples to Theorems 2.1(5), 2.2(5), 2.3(5) (3-cycles; x -> x+1 on N; the fixed-point step)
  K7  minimality: no 2-point counterexample; any 2-cycle violates (5); the paper's T(1) := 1 is forced
  K8  universality: (5) + WGP are satisfiable whenever T has no near-2-cycles
  K9  SHEET control: the paper's table, verbatim, makes 3n-1 satisfy every hypothesis (false conclusion)
  K10 the corrected Theorem 2.1 = one-step decay; it fails for Collatz for EVERY coefficient choice
  K11 Caristi: the only surviving metric fixed-point principle is equivalent to Collatz
  K12 contraction metrics: |x-y| fails; some metric <=> no cycles (Bessaga); a discrete one <=> Collatz

Every check raises on failure.  Runtime about 1-2 min; peak memory well under 300 MB; one process.
"""
import sys
import time
import itertools
import resource
from fractions import Fraction as Fr
from math import ceil

import numpy as np
import sympy as sp

T0 = time.time()
A_PAPER = Fr(1, 2)
B_PAPER = 2
M_PAPER = 2


def check(cond, msg):
    if not cond:
        raise AssertionError("CHECK FAILED: " + msg)


def hdr(s):
    print()
    print("=" * 100)
    print(s)
    print("=" * 100)


# ----------------------------------------------------------------------------------------------
# the paper's objects
# ----------------------------------------------------------------------------------------------
def T_paper(x):
    """Section 3 map: T(1)=1, x/2 on evens, (3x+1)/2 on odd x >= 3."""
    if x == 1:
        return 1
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


def T_true(x):
    """the shortcut map with its genuine 2-cycle {1,2}"""
    return x // 2 if x % 2 == 0 else (3 * x + 1) // 2


def T_minus(x):
    """the 3n-1 shortcut map; T_-(1) = 1 automatically"""
    return x // 2 if x % 2 == 0 else (3 * x - 1) // 2


def kind(z):
    if z == 1:
        return "one"
    return "even" if z % 2 == 0 else "odd"


TABLE = {  # (al, be, ga, de, ep, ze), Theorem 3.1 of the paper
    ("one", "one"): (1, 0, 0, 0, -1, 1),
    ("one", "even"): (1, 0, 0, -1, 0, 1),
    ("one", "odd"): (0, 0, 0, -2, 1, 2),
    ("even", "one"): (1, 0, 0, -1, 1, 0),
    ("even", "even"): (1, 0, -1, 0, -1, 1),
    ("even", "odd"): (0, 0, -2, 1, -2, 2),
    ("odd", "one"): (0, 0, 0, -2, 2, 1),
    ("odd", "even"): (0, -2, 0, 1, 2, -2),
}


def sub_odd_odd(k, l):
    """beta_0, delta_0, epsilon_0, zeta_0 for x = 2k+1, y = 2l+1 (k, l >= 1)"""
    dkl = k - l
    b0 = -2 if dkl <= -2 else (2 if dkl >= 2 else dkl)
    c1 = dkl <= -2 and 11 * k - 10 * l + 1 <= 0
    c2 = dkl >= 2 and -10 * k + 11 * l + 1 <= 0
    d0 = -2 if (c1 or c2) else -1
    e0 = 2 if c1 else 0
    z0 = 2 if c2 else 0
    return b0, d0, e0, z0


def coeffs_paper(x, y):
    kx, ky = kind(x), kind(y)
    if (kx, ky) in TABLE:
        return TABLE[(kx, ky)]
    k, l = (x - 1) // 2, (y - 1) // 2
    b0, d0, e0, z0 = sub_odd_odd(k, l)
    return (2, b0, -b0, d0, e0, z0)


def dist_vec(x, y, T, d):
    Tx, Ty = T(x), T(y)
    return (d(Tx, Ty) ** 2, d(x, Ty) ** 2, d(Tx, y) ** 2, d(x, y) ** 2, d(x, Tx) ** 2, d(y, Ty) ** 2)


def dabs(a, b):
    return abs(a - b)


def wgp_lhs(c, dv):
    return sum(ci * di for ci, di in zip(c, dv))


def alt1(c, A, B=None):
    al, be, ga, de, ep, ze = c
    m = min(be, 0)
    den = al + ze + 2 * m
    if den <= 0:
        return False
    if -(de + ep + 2 * m) > A * den:
        return False
    if B is not None and al + be + ze < B:
        return False
    return True


def alt2(c, A, B=None):
    al, be, ga, de, ep, ze = c
    m = min(ga, 0)
    den = al + ep + 2 * m
    if den <= 0:
        return False
    if -(de + ze + 2 * m) > A * den:
        return False
    if B is not None and al + ga + ep < B:
        return False
    return True


def cond5(c, A, B, M):
    return (alt1(c, A, B) or alt2(c, A, B)) and max(abs(ci) for ci in c) <= M


def verify_system(X, T, d, coeff, A, B, M, label):
    """all hypotheses of Theorem 2.3(5), lambda = 0, on the finite list X (T values may leave X)."""
    n_alt1 = n_alt2 = 0
    for x in X:
        for y in X:
            c = coeff(x, y)
            dv = dist_vec(x, y, T, d)
            check(wgp_lhs(c, dv) <= 0, f"{label}: WGP inequality at {(x, y)}")
            check(cond5(c, A, B, M), f"{label}: condition (5) at {(x, y)}")
            n_alt1 += alt1(c, A, B)
            n_alt2 += alt2(c, A, B)
    return n_alt1, n_alt2


# ----------------------------------------------------------------------------------------------
hdr("K1  Theorem 3.1 (the paper's WGP inequality for its T), exact, all 1 <= x, y < 400")
# ----------------------------------------------------------------------------------------------
N1 = 400
worst = None
nzero = 0
for x in range(1, N1):
    for y in range(1, N1):
        v = wgp_lhs(coeffs_paper(x, y), dist_vec(x, y, T_paper, dabs))
        check(v <= 0, f"Thm 3.1 inequality at {(x, y)}")
        nzero += (v == 0)
        if worst is None or v > worst[0]:
            worst = (v, x, y)
print(f"  all {(N1 - 1) ** 2} pairs satisfy the inequality; max LHS = {worst[0]} (at {worst[1:]}); "
      f"LHS = 0 at {nzero} pairs")
# numpy extension to 3000 (same table, vectorized, int64 exact)


def np_T_paper(x):
    return np.where(x == 1, 1, np.where(x % 2 == 0, x // 2, (3 * x + 1) // 2))


def np_T_minus(x):
    return np.where(x % 2 == 0, x // 2, (3 * x - 1) // 2)


def np_coeffs(x, y):
    """vectorized table; x array, y scalar"""
    n = x.shape[0]
    out = np.zeros((6, n), dtype=np.int64)
    kx = np.where(x == 1, 0, np.where(x % 2 == 0, 1, 2))
    ky = 0 if y == 1 else (1 if y % 2 == 0 else 2)
    names = ["one", "even", "odd"]
    for kxv in range(3):
        m = kx == kxv
        key = (names[kxv], names[ky])
        if key in TABLE:
            out[:, m] = np.array(TABLE[key], dtype=np.int64)[:, None]
    m = (kx == 2) & (ky == 2)
    if ky == 2 and np.any(m):
        k = (x[m] - 1) // 2
        l = (y - 1) // 2
        dkl = k - l
        b0 = np.where(dkl <= -2, -2, np.where(dkl >= 2, 2, dkl))
        c1 = (dkl <= -2) & (11 * k - 10 * l + 1 <= 0)
        c2 = (dkl >= 2) & (-10 * k + 11 * l + 1 <= 0)
        d0 = np.where(c1 | c2, -2, -1)
        e0 = np.where(c1, 2, 0)
        z0 = np.where(c2, 2, 0)
        out[0, m] = 2
        out[1, m] = b0
        out[2, m] = -b0
        out[3, m] = d0
        out[4, m] = e0
        out[5, m] = z0
    return out


def np_check_table(Tfun, N, label):
    x = np.arange(1, N, dtype=np.int64)
    Tx = Tfun(x)
    mx = -10 ** 18
    for y0 in range(1, N):
        Ty = int(Tfun(np.array([y0]))[0])
        c = np_coeffs(x, y0)
        lhs = (c[0] * (Tx - Ty) ** 2 + c[1] * (x - Ty) ** 2 + c[2] * (Tx - y0) ** 2 + c[3] * (x - y0) ** 2
               + c[4] * (x - Tx) ** 2 + c[5] * (y0 - Ty) ** 2)
        check(int(lhs.max()) <= 0, f"{label}: vectorized WGP inequality at y={y0}")
        mx = max(mx, int(lhs.max()))
    return mx


# cross-check vectorized table against the scalar one on a sample
_x = np.arange(1, 300, dtype=np.int64)
for _y in (1, 2, 3, 7, 50, 51, 299):
    _c = np_coeffs(_x, _y)
    for i, xv in enumerate(_x):
        check(tuple(int(t) for t in _c[:, i]) == coeffs_paper(int(xv), _y), "vectorized table = scalar table")
mx3000 = np_check_table(np_T_paper, 3000, "paper T")
print(f"  vectorized extension: all pairs 1 <= x, y < 3000 satisfy it too (max LHS {mx3000})")

# symbolic re-derivation of the paper's nine case polynomials
k, l = sp.symbols("k l", integer=True, positive=True)


def sym_lhs(c, x, y, Tx, Ty):
    al, be, ga, de, ep, ze = c
    return sp.expand(al * (Tx - Ty) ** 2 + be * (x - Ty) ** 2 + ga * (Tx - y) ** 2 + de * (x - y) ** 2
                     + ep * (x - Tx) ** 2 + ze * (y - Ty) ** 2)


def odd_image(z_half, b):
    # T(2z+1) = 3z + 1 + (b+1)/2 : 3z+2 for b=+1, 3z+1 for b=-1
    return 3 * z_half + 1 + sp.Rational(b + 1, 2)


paper_claims = {
    ("one", "even"): -2 * l ** 2 + 2 * l,
    ("one", "odd"): -6 * l ** 2 + 4 * l + 2,
    ("even", "one"): -2 * k ** 2 + 2 * k,
    ("even", "even"): -k ** 2 + 2 * k * l - 2 * l ** 2,
    ("even", "odd"): -2 * l ** 2 - 4 * k * l - 2 * k - 1,   # as printed in the paper
    ("odd", "one"): -6 * k ** 2 + 4 * k + 2,
    ("odd", "even"): -2 * k ** 2 + 1,
}


def case_poly(kx, ky, b):
    xs = {"one": sp.Integer(1), "even": 2 * k, "odd": 2 * k + 1}[kx]
    ys = {"one": sp.Integer(1), "even": 2 * l, "odd": 2 * l + 1}[ky]
    Txs = {"one": sp.Integer(1), "even": k, "odd": odd_image(k, b)}[kx]
    Tys = {"one": sp.Integer(1), "even": l, "odd": odd_image(l, b)}[ky]
    return sym_lhs(TABLE[(kx, ky)], xs, ys, Txs, Tys)


print("  the paper's case polynomials, re-derived (b = +1):")
slip = None
for key, claim in paper_claims.items():
    got = case_poly(key[0], key[1], +1)
    ok = sp.expand(got - claim) == 0
    print(f"    {key}: derived {got};  paper prints {sp.expand(claim)};  {'agree' if ok else 'DISAGREE'}")
    if not ok:
        slip = (key, got, claim)
check(slip is not None and slip[0] == ("even", "odd") and sp.expand(slip[1] - (-2 * l ** 2 + 1)) == 0,
      "exactly one slip: (even, odd) is -2l^2+1, not -2l^2-4kl-2k-1")
# the slip is harmless: -2 l^2 + 1 <= -1 for l >= 1
print("    -> one algebra slip: for (x even, y odd >= 3) the value is -2l^2 + 1 (<= -1), not the printed")
print("       -2l^2 - 4kl - 2k - 1 (<= -9).  Harmless: the inequality still holds.  Checked numerically:")
check(wgp_lhs(TABLE[("even", "odd")], dist_vec(2, 3, T_paper, dabs)) == -1, "(2,3) gives -1, not -9")
print("       at (x,y) = (2,3) the LHS is -1, not the printed -9.")
# odd/odd generic form and the five sub-cases
b0, d0, e0, z0 = sp.symbols("beta0 delta0 epsilon0 zeta0")
x_, y_ = 2 * k + 1, 2 * l + 1
gen = sym_lhs((2, b0, -b0, d0, e0, z0), x_, y_, 3 * k + 2, 3 * l + 2)
claim_gen = (18 + 4 * d0) * (k - l) ** 2 - 5 * b0 * (k - l) * (k + l + 2) + e0 * (k + 1) ** 2 + z0 * (l + 1) ** 2
check(sp.expand(gen - claim_gen) == 0, "odd/odd generic polynomial")
subcases_plus = [
    ((-2, -2, 2, 0), 2 * (k + 1) * (11 * k - 10 * l + 1)),
    ((-2, -1, 0, 0), 4 * (k - l) * (6 * k - l + 5)),
    ((2, -2, 0, 2), 2 * (l + 1) * (-10 * k + 11 * l + 1)),
    ((2, -1, 0, 0), 4 * (k - l) * (k - 6 * l - 5)),
]
for vals, fac in subcases_plus:
    e = sp.expand(claim_gen.subs({b0: vals[0], d0: vals[1], e0: vals[2], z0: vals[3]}) - fac)
    check(e == 0, f"odd/odd sub-case factorization {vals}")
e = sp.expand(claim_gen.subs({b0: k - l, d0: -1, e0: 0, z0: 0}) - (k - l) ** 2 * (4 - 5 * (k + l)))
check(e == 0, "odd/odd 'other cases' factorization")
print("  odd/odd: generic form (18+4d0)(k-l)^2 - 5b0(k-l)(k+l+2) + e0(k+1)^2 + z0(l+1)^2 and all five")
print("  sub-case factorizations of the paper re-derived symbolically: agree.")
print("  VERDICT K1: Theorem 3.1 is TRUE (the Collatz-specific algebra holds), with one harmless slip.")

# ----------------------------------------------------------------------------------------------
hdr("K2  Condition (5) of Theorem 2.3 (lambda = 0, A = 1/2, B = 2, M = 2), every pair 1 <= x, y < 400")
# ----------------------------------------------------------------------------------------------
cnt = {"alt1 only": 0, "alt2 only": 0, "both": 0}
tuples_gt, tuples_le = set(), set()
where_alt1_fails, where_alt2_fails = None, None
for x in range(1, N1):
    for y in range(1, N1):
        c = coeffs_paper(x, y)
        a1, a2 = alt1(c, A_PAPER, B_PAPER), alt2(c, A_PAPER, B_PAPER)
        check(cond5(c, A_PAPER, B_PAPER, M_PAPER), f"(5) at {(x, y)}")
        cnt["both" if a1 and a2 else ("alt1 only" if a1 else "alt2 only")] += 1
        if kind(x) == "odd" and kind(y) == "odd":
            (tuples_gt if x > y else tuples_le).add(c)
        if not a1 and where_alt1_fails is None:
            where_alt1_fails = (x, y)
        if not a2 and where_alt2_fails is None:
            where_alt2_fails = (x, y)
print(f"  (5) holds at all {(N1 - 1) ** 2} pairs; census: {cnt}")
check(tuples_gt == {(2, 2, -2, -2, 0, 2), (2, 2, -2, -1, 0, 0), (2, 1, -1, -1, 0, 0)},
      "odd/odd x>y combinations as listed by the paper")
check(tuples_le == {(2, -2, 2, -2, 2, 0), (2, -2, 2, -1, 0, 0), (2, -1, 1, -1, 0, 0), (2, 0, 0, -1, 0, 0)},
      "odd/odd x<=y combinations as listed by the paper")
print(f"  odd/odd coefficient tuples: x > y {sorted(tuples_gt)}; x <= y {sorted(tuples_le)} (= the paper's lists)")
check(not alt1(coeffs_paper(3, 2), A_PAPER) and not alt2(coeffs_paper(2, 3), A_PAPER), "global readings fail")
print("  alt1 fails at (3,2) = (odd, even) and alt2 fails at (2,3) = (even, odd), even without the B-clause")
print(f"  (first failures with the B-clause: alt1 at {where_alt1_fails}, alt2 at {where_alt2_fails}):")
print("  so the GLOBAL reading of (5) ('alt1 at every pair' or 'alt2 at every pair') is NOT satisfied;")
print("  Section 3 uses the PAIRWISE reading (its own verification switches alternatives by pair type).")
print("  VERDICT K2: condition (5), pairwise reading, holds for the Collatz table (at least below 400).")

# ----------------------------------------------------------------------------------------------
hdr("K3  The swap gap: the proof of Theorem 2.1(5) needs alt1(p,Tp) OR alt2(Tp,p) (A only, no B)")
# ----------------------------------------------------------------------------------------------
print("  Proof step A (x,y) = (T^(n-1)x, T^n x) in the beta-eliminated inequality uses alt1 at (p, Tp);")
print("  proof step B (x,y) = (T^n x, T^(n-1)x) in the gamma-eliminated inequality uses alt2 at (Tp, p).")
print("  Hypothesis (5) gives [alt1 or alt2] at (p,Tp) and, separately, at (Tp,p).  Uncovered case:")
print("  alt2 only at (p,Tp) and alt1 only at (Tp,p).")
uncovered = []
for p in range(1, 5000):
    q = T_paper(p)
    c_f, c_b = coeffs_paper(p, q), coeffs_paper(q, p)
    check(cond5(c_f, A_PAPER, B_PAPER, M_PAPER) and cond5(c_b, A_PAPER, B_PAPER, M_PAPER),
          "(5) holds at both (p,Tp) and (Tp,p)")
    if not (alt1(c_f, A_PAPER) or alt2(c_b, A_PAPER)):
        uncovered.append(p)
        check(alt2(c_f, A_PAPER, B_PAPER) and not alt1(c_f, A_PAPER), "uncovered: alt2 only at (p,Tp)")
        check(alt1(c_b, A_PAPER, B_PAPER) and not alt2(c_b, A_PAPER), "uncovered: alt1 only at (Tp,p)")
odd_ge3 = [p for p in range(3, 5000, 2)]
check(uncovered == odd_ge3, "uncovered steps = odd p >= 3")
print(f"  steps p -> Tp with 1 <= p < 5000: {len(uncovered)} uncovered = exactly the odd p >= 3 "
      f"({len(odd_ge3)} of them); every even step and p = 1 are covered.")
print("  (5) nevertheless holds at both pairs of every step.  VERDICT K3: the proof of Theorem 2.1(5)")
print("  has a quantifier gap, and for the Collatz table it is hit at EVERY up-step.")

# ----------------------------------------------------------------------------------------------
hdr("K4  The paper's intermediate claim and orbit-radius bound")
# ----------------------------------------------------------------------------------------------


def orbit_paper(x, cap=10 ** 6):
    o = [x]
    while o[-1] != 1:
        o.append(T_paper(o[-1]))
        check(len(o) < cap, "orbit too long")
    return o


o5 = orbit_paper(5)
lhs5 = (o5[1] - o5[2]) ** 2
rhs5 = A_PAPER * (o5[0] - o5[1]) ** 2
check(o5[:5] == [5, 8, 4, 2, 1] and lhs5 == 16 and rhs5 == Fr(9, 2) and lhs5 > rhs5, "x=5, n=1")
print(f"  claim: d(T^n x, T^(n+1) x)^2 <= A^n d(x,Tx)^2 (A = 1/2).  x = 5: orbit {o5};  n = 1: "
      f"d(8,4)^2 = {lhs5} > A d(5,8)^2 = {rhs5}.  FALSE.")
holds = []
for x in range(1, 10001):
    o = orbit_paper(x)
    d0sq = (o[0] - o[1]) ** 2 if len(o) > 1 else 0
    ok = all((o[n] - o[n + 1]) ** 2 <= A_PAPER ** n * d0sq for n in range(1, len(o) - 1))
    if ok:
        holds.append(x)
    if x % 2 == 1 and x >= 3:
        # fails already at n = 1: d(Tx,T^2x)^2 >= (16/9) d(x,Tx)^2 > A d(x,Tx)^2
        check((o[1] - o[2]) ** 2 > A_PAPER * d0sq, "odd x >= 3 fails at n = 1")
check(all(x % 2 == 0 or x == 1 for x in holds) and 96 in holds and 5 not in holds, "holds-set shape")
print(f"  census x <= 10^4: the claim fails at n = 1 for EVERY odd x >= 3; it holds for all n at only "
      f"{len(holds)} starts,")
print("  all 1 or even (a long initial halving phase can pay for the later stretch, e.g. x = 96 = 3*2^5).")
o27 = orbit_paper(27)
bound27 = 27 + (o27[1] - o27[0]) / (1 - 2 ** -0.5)
check(max(o27) == 4616 and 74 < bound27 < 75, "x=27 orbit maximum and bound")
print(f"  orbit radius: the proof gives d(x, T^m x) <= d(x,Tx)/(1 - sqrt(A)), i.e. T^m(27) <= 27 + "
      f"14/(1-1/sqrt2) = {bound27:.2f};  actual max T^m(27) = {max(o27)} (shortcut orbit, {len(o27) - 1} steps).")
viol_radius = 0
for x in range(2, 10001):
    o = orbit_paper(x)
    if max(abs(v - x) for v in o) > abs(o[1] - x) / (1 - 2 ** -0.5) + 1e-9:
        viol_radius += 1
print(f"  the radius bound fails for {viol_radius} of the 9999 starts 2 <= x <= 10^4.")

# ----------------------------------------------------------------------------------------------
hdr("K5  Up-steps stretch consecutive distances; the Mersenne family")
# ----------------------------------------------------------------------------------------------
minr, maxr = None, Fr(0)
for p in range(3, 200001, 2):
    q = T_true(p)
    r = Fr(abs(T_true(q) - q), abs(q - p))
    want = Fr(3, 2) if p % 4 == 3 else Fr(3 * p + 1, 2 * (p + 1))
    check(r == want, f"stretch formula at p={p}")
    if minr is None or r < minr[0]:
        minr = (r, p)
    maxr = max(maxr, r)
for p in range(2, 200001, 2):
    q = T_paper(p)
    if q != p:
        check(Fr(abs(T_paper(q) - q), abs(q - p)) <= 1, "even steps never stretch")
check(minr == (Fr(4, 3), 5) and maxr == Fr(3, 2), "min 4/3 at p=5, sup 3/2")
print("  for odd p >= 3:  d(Tp,T^2 p)/d(p,Tp) = 3/2 (p = 3 mod 4) and (3p+1)/(2(p+1)) (p = 1 mod 4)")
print(f"  (exact, p < 2*10^5); minimum {minr[0]} at p = {minr[1]}, supremum 3/2; even steps never stretch.")
for m in range(2, 81):
    x = 2 ** m - 1
    o = [x]
    for _ in range(m):
        o.append(T_true(o[-1]))
    check(all(o[j] == 3 ** j * 2 ** (m - j) - 1 for j in range(m + 1)), "Mersenne orbit formula")
    check(Fr(o[m] - o[m - 1], o[1] - o[0]) == Fr(3, 2) ** (m - 1), "Mersenne stretch (3/2)^(m-1)")
    check(Fr(o[m] - o[0], o[1] - o[0]) == Fr(3 ** m - 2 ** m, 2 ** (m - 1)), "Mersenne radius ratio")
print("  x = 2^m - 1 (m = 2..80): T^j x = 3^j 2^(m-j) - 1 (j <= m), so")
print("     d(T^(m-1)x, T^m x) / d(x,Tx) = (3/2)^(m-1)  and  d(x, T^m x)/d(x,Tx) = (3^m - 2^m)/2^(m-1).")
print("  Hence no bound d(T^n x,T^(n+1)x) <= g(n) d(x,Tx) uniform in x holds with g(n) < (3/2)^n, and no")
print("  uniform orbit-radius bound d(x,T^m x) <= K d(x,Tx) holds for any K.")

# ----------------------------------------------------------------------------------------------
hdr("K6  Counterexamples to Theorems 2.1(5), 2.2(5), 2.3(5) (pairwise reading)")
# ----------------------------------------------------------------------------------------------
# (a) the coordinator's 3-cycle on {0,1,2}, M = 4
cyc = {0: 1, 1: 2, 2: 0}
Tc = cyc.__getitem__
MC = 4


def coeff_coord(x, y):
    if x == y:
        return (MC, 0, -MC, 0, 0, 0)
    if Tc(x) == y:
        return (1, -MC, MC, 0, 0, 0)
    return (1, MC, -MC, 0, 0, 0)


n1, n2 = verify_system([0, 1, 2], Tc, dabs, coeff_coord, A_PAPER, B_PAPER, MC, "3-cycle coordinator")
check(all(Tc(x) != x for x in cyc), "no fixed point")
# M = 4 is forced for this particular scheme at the pair (1,2): 4 <= M*1
check(wgp_lhs((1, -3, 3, 0, 0, 0), dist_vec(1, 2, Tc, dabs)) > 0, "scheme needs M >= 4")
print(f"  (a) X = {{0,1,2}} (|x-y|), T = 3-cycle 0->1->2->0, lambda=0, A=1/2, B=2, M=4, coefficients")
print("      (M,0,-M,0,0,0) on the diagonal, (1,-M,M,0,0,0) at (p,Tp), (1,M,-M,0,0,0) at (Tp,p):")
print(f"      all 9 WGP inequalities and (5) hold (alt1 at {n1} pairs, alt2 at {n2}); no fixed point; the orbit")
print("      0,1,2,0,... is not Cauchy.  [This scheme needs M >= 4 at the pair (1,2).]")

# (b) the paper's own constants M = 2 on 3-point spaces
SPARSE = {"diag": (2, 0, 0, 0, 0, 0)}


def sparse_table(X, T):
    tab = {}
    grid = sorted(itertools.product(range(-2, 3), repeat=6),
                  key=lambda c: (sum(1 for ci in c if ci), sum(abs(ci) for ci in c), c))
    for x in X:
        for y in X:
            dv = dist_vec(x, y, T, dabs)
            for c in grid:
                if wgp_lhs(c, dv) <= 0 and cond5(c, A_PAPER, B_PAPER, M_PAPER):
                    tab[(x, y)] = c
                    break
            check((x, y) in tab, "sparse table exists")
    return tab


for label, X, T in (("{0,1,2}, 3-cycle", [0, 1, 2], {0: 1, 1: 2, 2: 0}),
                    ("{5,7,10}, the 3n-1 cycle 5->7->10->5", [5, 7, 10], {5: 7, 7: 10, 10: 5})):
    Tm = T.__getitem__
    tab = sparse_table(X, Tm)
    n1, n2 = verify_system(X, Tm, dabs, lambda a, b: tab[(a, b)], A_PAPER, B_PAPER, M_PAPER, label)
    check(all(Tm(x) != x for x in X), "no fixed point")
    print(f"  (b) X = {label}, the PAPER'S constants (lambda=0, A=1/2, B=2, M=2); table:")
    for (x, y), c in tab.items():
        role = "diag" if x == y else ("(p,Tp)" if Tm(x) == y else "(Tp,p)")
        print(f"        {(x, y)} {role:7s} {c}  {'alt1' if alt1(c, A_PAPER, B_PAPER) else 'alt2'}")
    print("      all hypotheses hold; T has no fixed point.")
for x in (5, 7, 10):
    check(T_minus(x) == {5: 7, 7: 10, 10: 5}[x], "{5,7,10} is a T_- cycle")
# equilateral (discrete) 3-point space with the two-vector scheme
ddisc = lambda a, b: 0 if a == b else 1
V1, V2 = (1, 1, -2, 0, 0, 0), (1, -2, 1, 0, 0, 0)


def coeff_disc(x, y):
    return V2 if Tc(x) == y else V1


verify_system([0, 1, 2], Tc, ddisc, coeff_disc, A_PAPER, B_PAPER, M_PAPER, "discrete 3-cycle")
print("      (also: the equilateral 3-point space with V2 = (1,-2,1,0,0,0) at (p,Tp), V1 = (1,1,-2,0,0,0)")
print("       elsewhere satisfies everything with the paper's constants.)")

# (c) translation on N with the paper's constants


def T_shift(x):
    return x + 1


def coeff_shift(x, y):
    return V1 if x >= y else V2


verify_system(range(1, 1201), T_shift, dabs, coeff_shift, A_PAPER, B_PAPER, M_PAPER, "shift")
e = sp.symbols("e", integer=True)
# x - y = e >= 0 with V1:  e^2 + (e-1)^2 - 2(e+1)^2 ;  e < 0 with V2: e^2 - 2(e-1)^2 + (e+1)^2
check(sp.expand(e ** 2 + (e - 1) ** 2 - 2 * (e + 1) ** 2 - (-6 * e - 1)) == 0, "shift closed form, e>=0")
check(sp.expand(e ** 2 - 2 * (e - 1) ** 2 + (e + 1) ** 2 - (6 * e - 1)) == 0, "shift closed form, e<0")
print("  (c) X = N = {1,2,...} (the paper's own space), T(x) = x + 1, lambda=0, A=1/2, B=2, M=2 (the")
print("      paper's own constants): coefficients V1 = (1,1,-2,0,0,0) [alt1] if x >= y, V2 = (1,-2,1,0,0,0)")
print("      [alt2] if x < y.  With e = x - y the WGP left side is -6e-1 (e >= 0) and 6e-1 (e < 0): <= 0.")
print("      Checked exactly for 1 <= x,y <= 1200.  T^n x = x + n is not Cauchy and T has no fixed point.")
print("      => Theorems 2.1(5), 2.2(5), 2.3(5) are FALSE, on (N,|x-y|) with the paper's constants.")

# (d) the fixed-point step: orbits converge, limit not fixed.  X = {0} u {1/n}, T(1/n)=1/(n+1), T(0)=1


def T_harm(x):
    return Fr(1) if x == 0 else Fr(1, x.denominator + 1) if x.numerator == 1 else None


K_H = 5


def coeff_univ(K):
    def cf(x, y, T=T_harm):
        u = abs(T(x) - y)
        v = abs(x - T(y))
        if u >= v:
            return (1, 1, -K, 0, 0, 0)
        return (1, -K, 1, 0, 0, 0)
    return cf


Xh = [Fr(0)] + [Fr(1, n) for n in range(1, 301)]
n1, n2 = verify_system(Xh, T_harm, dabs, coeff_univ(K_H), A_PAPER, B_PAPER, K_H, "harmonic")
check(T_harm(Fr(0)) == 1 and all(T_harm(x) != x for x in Xh), "no fixed point")
# the swap gap in the fixed-point step: need alt1 at (T^n x, u) or alt2 at (u, T^n x) with u = 0
cf = coeff_univ(K_H)
for n in range(2, 301):
    xn = Fr(1, n)
    check(not alt1(cf(xn, Fr(0)), A_PAPER) and not alt2(cf(Fr(0), xn), A_PAPER), "fixed-point step uncovered")
print("  (d) X = {0} u {1/n : n >= 1} (complete), T(1/n) = 1/(n+1), T(0) = 1.  Coefficients (1,1,-5,0,0,0)")
print("      when d(Tx,y) >= d(x,Ty), else (1,-5,1,0,0,0); lambda=0, A=1/2, B=2, M=5.  All hypotheses hold")
print("      (checked exactly on {0} u {1/n: n <= 300}; proof in the note).  Every orbit converges to 0,")
print("      but T(0) = 1: Theorem 2.3(5) fails even where Theorem 2.2's conclusion holds.  The fixed-point")
print("      step of its proof needs alt1 at (T^n x, u) or alt2 at (u, T^n x); here, for every n >= 2, the")
print("      table has alt2 at (1/n, 0) and alt1 at (0, 1/n): the same swap gap.")

# ----------------------------------------------------------------------------------------------
hdr("K7  Minimality: 1 and 2 points are impossible, 3 points suffice; every 2-cycle violates (5)")
# ----------------------------------------------------------------------------------------------
print("  A 1-point space has a fixed point.  A fixed-point-free map on 2 points is the swap a <-> b; at the")
print("  pair (a,b) the distance vector is D^2 (1,0,0,1,1,1), so the WGP inequality reads al+de+ep+ze <= 0,")
print("  while alt1 forces al+de+ep+ze >= (1-A)(al+ze+2min(be,0)) - 4 min(be,0) > 0 and alt2 forces")
print("  al+de+ep+ze >= (1-A)(al+ep+2min(ga,0)) - 4 min(ga,0) > 0.  This uses only A < 1; lambda-mixing")
print("  (Lemma 2.2) produces another coefficient vector obeying the same inequality at the same pair,")
print("  so no lambda helps.  The same computation applies to any 2-cycle {a,b} in any metric space.")
cnt_bad = 0
grid = list(itertools.product(range(-3, 4), repeat=6))
for Aval in (Fr(1, 2), Fr(99, 100)):
    for c in grid:
        if c[0] + c[3] + c[4] + c[5] <= 0 and (alt1(c, Aval) or alt2(c, Aval)):
            cnt_bad += 1
check(cnt_bad == 0, "no grid vector satisfies swap inequality and alt1/alt2")
print(f"  sanity grid: none of the {len(grid)} vectors in {{-3..3}}^6 passes both, for A = 1/2 and A = 99/100.")
# the genuine Collatz 2-cycle {1,2}
dv12 = dist_vec(1, 2, T_true, dabs)
check(dv12 == (1, 0, 0, 1, 1, 1), "true T at (1,2)")
check(wgp_lhs(coeffs_paper(1, 2), dv12) == 1, "paper's table fails the true T at (1,2)")
print(f"  With the true shortcut map (T(1) = 2), the pair (1,2) has distance vector {dv12}: the 2-cycle lemma")
print("  forbids every coefficient choice there (the paper's own table gives LHS = +1 > 0).  Hence the paper's")
print("  redefinition T(1) := 1 is FORCED by condition (5).  3-point counterexamples exist (K6a,b).")
print("  VERDICT K7: the least counterexample has exactly 3 points (a 3-cycle); by K6b it can be taken inside")
print("  N with the paper's constants: {5,7,10} under 3n-1.")

# ----------------------------------------------------------------------------------------------
hdr("K8  Universality: (5) + WGP carry no dynamics beyond excluding (near-)2-cycles")
# ----------------------------------------------------------------------------------------------
print("  Proposition U. If d(Tx,Ty) <= C max(d(Tx,y), d(x,Ty)) for all x,y, then with K = (B/2)(C^2+1) the")
print("  coefficients (B/2,B/2,-K,0,0,0) [alt1] when d(Tx,y) >= d(x,Ty), else (B/2,-K,B/2,0,0,0) [alt2],")
print("  satisfy WGP and (5) (lambda = 0, any A in (0,1)), with M = max(B/2, K).")


def near2_sup(Tnp, N):
    x = np.arange(1, N, dtype=np.int64)
    Tx = Tnp(x)
    best = (0.0, None)
    for y0 in range(1, N):
        Ty = int(Tnp(np.array([y0]))[0])
        u = np.abs(Tx - y0)
        v = np.abs(x - Ty)
        w = np.abs(Tx - Ty)
        m = np.maximum(u, v)
        check(not np.any((m == 0) & (w > 0)), "no exact 2-cycle")
        r = np.where(m > 0, w / np.maximum(m, 1), 0.0)
        i = int(np.argmax(r))
        if r[i] > best[0]:
            best = (float(r[i]), (int(x[i]), y0))
    return best


for label, Tnp in (("paper T (3x+1, T(1)=1)", np_T_paper), ("3x-1", np_T_minus),
                   ("x+1", lambda z: z + 1)):
    s = near2_sup(Tnp, 2500)
    print(f"    C({label}) over 1 <= x,y < 2500: sup = {s[0]:.4f} at {s[1]}")
print("  For T_+ (paper) and T_- on N the ratio is bounded (<= 11 by the elementary argument in the note;")
print("  numerically 5 and < 4).  So BOTH sheets admit bounded coefficient tables satisfying every")
print("  hypothesis of Theorem 2.3(5): the hypothesis carries no dynamical information at all.")

# ----------------------------------------------------------------------------------------------
hdr("K9  SHEET control: the paper's table, VERBATIM, for 3n-1 on N")
# ----------------------------------------------------------------------------------------------
for x in range(1, N1):
    for y in range(1, N1):
        check(wgp_lhs(coeffs_paper(x, y), dist_vec(x, y, T_minus, dabs)) <= 0, f"3n-1 WGP at {(x, y)}")
mxm = np_check_table(np_T_minus, 3000, "3n-1")
print(f"  WGP inequality for T_-(x) = x/2, (3x-1)/2 with the paper's table: all pairs < 400 exact, all")
print(f"  pairs < 3000 vectorized (max LHS {mxm}).  Symbolic case polynomials (b = -1):")
minus_expect = {
    ("one", "even"): -2 * l ** 2 + 2 * l, ("one", "odd"): -6 * l ** 2,
    ("even", "one"): -2 * k ** 2 + 2 * k, ("even", "even"): -k ** 2 + 2 * k * l - 2 * l ** 2,
    ("even", "odd"): -2 * l ** 2 - 4 * l - 1, ("odd", "one"): -6 * k ** 2, ("odd", "even"): -2 * k ** 2 - 4 * k - 1,
}
for key, val in minus_expect.items():
    got = case_poly(key[0], key[1], -1)
    check(sp.expand(got - val) == 0, f"3n-1 case {key}")
    print(f"    {key}: {got}")
genm = sym_lhs((2, b0, -b0, d0, e0, z0), x_, y_, 3 * k + 1, 3 * l + 1)
claim_m = (18 + 4 * d0) * (k - l) ** 2 - 5 * b0 * (k - l) * (k + l) + e0 * k ** 2 + z0 * l ** 2
check(sp.expand(genm - claim_m) == 0, "3n-1 odd/odd generic")
subm = [((-2, -2, 2, 0), 2 * k * (11 * k - 10 * l)), ((-2, -1, 0, 0), 4 * (k - l) * (6 * k - l)),
        ((2, -2, 0, 2), 2 * l * (11 * l - 10 * k)), ((2, -1, 0, 0), 4 * (k - l) * (k - 6 * l))]
for vals, fac in subm:
    check(sp.expand(claim_m.subs({b0: vals[0], d0: vals[1], e0: vals[2], z0: vals[3]}) - fac) == 0, "3n-1 sub")
check(sp.expand(claim_m.subs({b0: k - l, d0: -1, e0: 0, z0: 0}) - (k - l) ** 2 * (14 - 5 * (k + l))) == 0,
      "3n-1 other")
print("    odd/odd: (18+4d0)(k-l)^2 - 5b0(k-l)(k+l) + e0 k^2 + z0 l^2, sub-cases 2k(11k-10l),")
print("    4(k-l)(6k-l), 2l(11l-10k), 4(k-l)(k-6l), (k-l)^2(14-5(k+l)): each <= 0 on its range.")
print("  Condition (5) depends only on the table, so it holds verbatim (K2).  The only fixed point of T_- in N")
print("  is 1.  So Theorem 3.2's argument, verbatim, 'proves' that every positive integer reaches 1 under 3n-1.")
o = [5]
for _ in range(3):
    o.append(T_minus(o[-1]))
check(o == [5, 7, 10, 5], "3n-1 cycle")
print(f"  FALSE: {o}.  The method is SHEET-blind (the foundry's typing, now a theorem about this paper).")

# ----------------------------------------------------------------------------------------------
hdr("K10 The corrected Theorem 2.1 is one-step decay, and Collatz violates it for EVERY coefficient choice")
# ----------------------------------------------------------------------------------------------
print("  Corrected Theorem 2.1. If for every p in the orbit, alt1_lambda(p,Tp) or alt2_lambda(Tp,p) holds, then")
print("  d(Tp,T^2p)^2 <= A d(p,Tp)^2 along the orbit, hence d(T^n x,T^m x) <= A^(n/2)/(1-A^(1/2)) d(x,Tx).")
print("  Converse at a single step: if d(Tp,T^2p)^2 <= A d(p,Tp)^2, the vector (B,0,0,-AB,0,0) at (p,Tp)")
print("  satisfies WGP there and alt1 (with B).  So the corrected hypothesis IS one-step decay.")
decay_ok, cov = [], []
for p in range(1, 5000):
    q = T_paper(p)
    r = T_paper(q)
    dec = (r - q) ** 2 <= A_PAPER * (q - p) ** 2
    decay_ok.append(dec)
    cov.append(alt1(coeffs_paper(p, q), A_PAPER) or alt2(coeffs_paper(q, p), A_PAPER))
    if dec:
        vec = (B_PAPER, 0, 0, -A_PAPER * B_PAPER, 0, 0)
        check(wgp_lhs(vec, dist_vec(p, q, T_paper, dabs)) <= 0 and alt1(vec, A_PAPER, B_PAPER), "converse vector")
    else:
        # no coefficient choice can work: the ratio exceeds 1 > A
        check(Fr((r - q) ** 2, (q - p) ** 2) >= Fr(16, 9), "ratio >= 16/9 at failing steps")
check(decay_ok == cov, "covered steps = decay steps (Collatz table)")
check([p for p in range(1, 5000) if not decay_ok[p - 1]] == odd_ge3, "decay fails exactly at odd p >= 3")
print("  Collatz, p < 5000: one-step decay (A = 1/2) fails exactly at the odd p >= 3, where the ratio is")
print("  >= 16/9 > 1 > A for ANY A < 1; the Collatz table covers exactly the decay steps.  At p = 5 the")
print("  corrected hypothesis would force 16 <= 9A < 9.  Hence NO coefficient functions and NO lambda make")
print("  the corrected Theorem 2.1 apply to Collatz on (N,|x-y|) -- already at the single point p = 5.")
bad_m = []
for p in range(3, 20001, 2):
    q = T_minus(p)
    r_ = T_minus(q)
    if Fr((r_ - q) ** 2, (q - p) ** 2) < Fr(9, 4):
        bad_m.append(p)
check(bad_m == [], "3n-1 up-steps stretch by >= 3/2")
print("  (3n-1, odd p < 20001: every up-step stretches consecutive distances by at least 3/2, so the")
print("   corrected theorem fails at every 3n-1 up-step too.)")

# ----------------------------------------------------------------------------------------------
hdr("K11 Caristi: the one metric fixed-point principle that survives is equivalent to Collatz")
# ----------------------------------------------------------------------------------------------
print("  Caristi: if phi >= 0 (l.s.c.) with d(x,Tx) <= phi(x) - phi(Tx), then T has a fixed point.  On")
print("  (N,|x-y|) with the paper's T such phi exists iff every orbit reaches 1; then phi = the orbit's total")
print("  variation, sum over n < tau(x) of ceil(T^n x / 2), with equality d(x,Tx) = phi(x) - phi(Tx).")
LIM = 100000
phi = {1: 0}


def get_phi(x):
    stack = []
    while x not in phi:
        stack.append(x)
        x = T_paper(x)
    val = phi[x]
    for z in reversed(stack):
        val = val + abs(T_paper(z) - z)
        phi[z] = val
    return phi[stack[0]] if stack else val


for x in range(2, LIM + 1):
    get_phi(x)
    check(abs(x - T_paper(x)) == (x + 1) // 2, "|T(n)-n| = ceil(n/2) for n >= 2")
    check(phi[x] - phi[T_paper(x)] == abs(x - T_paper(x)), "Caristi equality")
print(f"  verified for 2 <= x <= {LIM}: |T(n)-n| = ceil(n/2) and d(x,Tx) = phi(x) - phi(Tx);"
      f" e.g. phi(27) = {phi[27]}, phi(97) = {phi[97]}.")
print("  A Caristi argument therefore needs exactly the finiteness of phi, i.e. Collatz itself.")

# ----------------------------------------------------------------------------------------------
hdr("K12 Contraction metrics for T: |x-y| fails; some metric <=> no cycles (Bessaga); a discrete one <=> Collatz")
# ----------------------------------------------------------------------------------------------
print("  (i)  d = |x-y|: no contraction-type certificate (K5, K10).")
print("  (ii) Bessaga's converse (1959; statement as in Wikipedia 'Banach fixed-point theorem', Converses): if")
print("       every iterate T^n of a self-map of a set has a unique fixed point, then for every q in (0,1) some")
print("       complete metric makes T a q-contraction.  Conversely a contraction's iterates have unique fixed")
print("       points.  For the paper's T on N (T(1) = 1) this says: a Banach-contraction metric exists iff N has")
print("       no T-cycle other than {1} -- the CYCLE half only; divergent orbits would still converge to 1 in it.")
print("  (iii) A UNIFORMLY DISCRETE complete metric making T a contraction exists iff every orbit reaches 1:")
print("       (=>) rho(T^n x, 1) <= q^n rho(x, 1) < inf rho forces T^n x = 1;  (<=) rho(x,y) = K^max(h(x),h(y))")
print("       for x != y (h = steps to 1, K = 1/q > 1) is an ultrametric, >= 1 off the diagonal, and")
print("       rho(Tx,Ty) = q rho(x,y) or 0.")
hgt = {1: 0}


def height(x):
    stack = []
    while x not in hgt:
        stack.append(x)
        x = T_paper(x)
    hv = hgt[x]
    for z in reversed(stack):
        hv += 1
        hgt[z] = hv
    return hgt[stack[0]] if stack else hv


Kq = 2   # q = 1/2
NR = 1500
for x in range(1, NR + 1):
    height(x)
    height(T_paper(x))


def rho(x, y):
    return 0 if x == y else Fr(Kq) ** max(hgt[x], hgt[y])


# rho(x,y) = 2^E(x,y) with E = max(h(x),h(y)); contraction by 1/2 <=> E(Tx,Ty) <= E(x,y) - 1 when Tx != Ty
hx = [0] + [hgt[x] for x in range(1, NR + 1)]
for x in range(1, NR + 1):
    tx = T_paper(x)
    for y in range(1, NR + 1):
        if x == y:
            continue
        ty = T_paper(y)
        if tx != ty:
            check(max(hgt[tx], hgt[ty]) <= max(hx[x], hx[y]) - 1, "rho-contraction with q = 1/2")
rng = np.random.default_rng(924)
for _ in range(20000):
    x, y, z = (int(v) for v in rng.integers(1, NR + 1, size=3))
    check(rho(x, z) <= max(rho(x, y), rho(y, z)), "ultrametric inequality")
    if x != y:
        check(rho(x, y) >= 1, "uniformly discrete")
check(rho(T_paper(27), T_paper(9)) == rho(27, 9) / 2, "equality case example")
print(f"  checked for 1 <= x, y <= {NR} (q = 1/2, heights from the actual orbits): rho >= 1 off the diagonal,")
print("  rho(Tx,Ty) <= rho(x,y)/2, and the ultrametric inequality on 20000 random triples.")
print("  So metric contraction reformulates Collatz exactly as 'the height h is finite everywhere'; the")
print("  divergence half is precisely the uniform discreteness that the paper's last step takes for granted.")

rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
rss_mb = rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024
print()
print(f"[fixedpt_kawasaki] ALL CHECKS PASSED  ({time.time() - T0:.1f} s, peak RSS {rss_mb:.0f} MB)",
      file=sys.stderr)
print("[fixedpt_kawasaki] ALL CHECKS PASSED")
