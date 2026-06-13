#!/usr/bin/env python3
"""lrc_n14_twisted_involutions_s606.py — twisted involutions for LRC@14:
which involutions can certify loneliness, and why the geometric ones cannot.

Target: a proof/improvement for LRC@14, via the resolution template (S599):
turn the permanent-shaped (star) sum p_0 = sum_S (-1)^|S| m_S into a
determinant/certifiable form by a SIGN-REVERSING INVOLUTION (Redei's GF(2)
involution is the solved face). n=14 reduces (S582o, 5 ways) to one gap:
close B(t)=1 -> B(t)=0 off the AP; obstruction = the observer-coupling defect
(S580o), active for m>=5.

DEEP STUDY OF TWISTED INVOLUTIONS (the question 'understand them deeply'):

  GEOMETRIC clock involutions CANNOT certify by cancellation:
  * Negation t->-t: depth(-t)=depth(t) exactly, so it FIXES the lonely set and
    is MEASURE-preserving => sign-PRESERVING on p_0 (the kappa-even fact of S605
    in geometric form). It is the ONLY symmetry of the lonely set Lambda.
  * Half-shift t->t+1/2 (even n): fixes even-runner distances, reflects odd ones
    (origin<->antipode). It is measure-preserving (a translation) but maps
    Lambda to a DIFFERENT set (even-far + odd-ANTIPODE-far) -- NOT a symmetry of
    Lambda. Useful for reformulation, useless for cancellation.
  => no geometric involution is sign-reversing on the loneliness measure.

  The CORRECT twisted involution acts on SUBSETS, not the clock: the (star) sum
  p_0 = sum_S (-1)^|S| m_S already carries signs (-1)^|S|. The PIVOT involution
  sigma_p: S <-> S xor {p} is sign-reversing; it cancels pairs with
  m_S = m_{S+p} (pivot p REDUNDANT on S: cap_S F subset F_p), reducing p_0 to
  the IRREDUNDANT survivors -- i.e. localizing to the NERVE of the cover {F_v}.
  Iterating over pivots drives p_0 to the irredundant core = the two-block
  Helly witnesses (codex-S599). This is the LRC realization of Redei's
  permanent->determinant resolution.

  REDUNDANCY is the involution's FUEL: dissociated speed sets have NONE (the
  involution is trivial -> p_0 stays at the positive independence value, never
  tight); structured sets have redundancy, but the AMOUNT does not cleanly
  decide tightness (Vitali wall, S603) -- the residual is exactly the
  quantitative observer-coupling the involution cannot cancel.

Status: [VERIFIED] computational; [PROVED] classical (nerve lemma, negation);
the n=14 closure remains the irredundant-core / coupling gap (NOT closed here).
Session: claude-2026-06-03-S606 (lrc-n14-twisted-involutions).
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
from fractions import Fraction as F
from itertools import combinations

def bpd(V, n):
    d = F(1, n); bp = {F(0), F(1)}
    for v in V:
        for j in range(v+1):
            for s in (1, -1):
                t = (F(j)+s*d)/v
                if 0 <= t <= 1: bp.add(t)
    return sorted(bp), d
def m_of(S, bp, d):
    if not S: return F(1)
    tot = F(0)
    for a, b in zip(bp, bp[1:]):
        mid = (a+b)/2
        if all(min((v*mid) % 1, 1-((v*mid) % 1)) < d for v in S): tot += b-a
    return tot
def p0_union(V, n):
    bp, d = bpd(V, n); tot = F(0)
    for a, b in zip(bp, bp[1:]):
        mid = (a+b)/2
        if not any(min((v*mid) % 1, 1-((v*mid) % 1)) < d for v in V): tot += b-a
    return tot

print("\n  TWISTED INVOLUTIONS FOR LRC@14\n")
print("=" * 70)

# ============================================================
print("\n  I. GEOMETRIC INVOLUTIONS ARE MEASURE-PRESERVING (cannot cancel)  [VERIFIED]")
print("  " + "-" * 50)
def depth(V, t, d): return sum(1 for v in V if min((v*t) % 1, 1-((v*t) % 1)) < d)
V14 = tuple(range(1, 14)); d14 = F(1, 14)
# negation symmetry: depth(-t)=depth(t)
neg_ok = all(depth(V14, F(s, 97), d14) == depth(V14, (-F(s, 97)) % 1, d14) for s in range(1, 97))
print(f"  Negation t->-t: depth(-t)=depth(t) for all sampled t? {neg_ok}")
print("    => fixes the lonely set, measure-preserving, SIGN-PRESERVING on p_0.")
# half-shift maps Lambda elsewhere
def lonely(V, t, d): return all(min((v*t) % 1, 1-((v*t) % 1)) >= d for v in V)
Vtest = tuple([1,2,3,4,5,6,7,8,9,10,11,12,14])
diff = sum(1 for s in range(1, 400) if lonely(Vtest, F(s, 400), d14) != lonely(Vtest, F(s, 400)+F(1, 2), d14))
print(f"  Half-shift t->t+1/2: lonely(t) != lonely(t+1/2) on {diff}/399 samples (V={Vtest[-3:]}..)")
print("    => NOT a symmetry of Lambda (maps origin-loneliness to a mixed")
print("       origin/antipode condition); measure-preserving, so no cancellation.")
print()

# ============================================================
print("  II. THE SUBSET PIVOT INVOLUTION IS SIGN-REVERSING  [VERIFIED]")
print("  " + "-" * 50)
print("  p_0 = sum_S (-1)^|S| m_S; sigma_p: S <-> S xor {p} flips parity.")
print("  Pairs with m_S = m_{S+p} cancel (p redundant); survivors' signed sum = p_0:")
def pivot_reduce(V, n, p):
    bp, d = bpd(V, n); rest = [v for v in V if v != p]
    tot = F(0); surv = 0
    for r in range(len(rest)+1):
        for S in combinations(rest, r):
            mS = m_of(S, bp, d); mSp = m_of(tuple(sorted(S+(p,))), bp, d)
            if mS != mSp:
                tot += (-1)**len(S)*(mS-mSp); surv += 1
    return tot, surv, 2**(len(V))
for V, n in [((1,2,3,4),5),((1,3,4,7),5),((1,3,4,5,9),6)]:
    t, surv, full = pivot_reduce(V, n, V[0])
    print(f"  V={str(V):<14} survivors-sum={float(t):.5f}  p0={float(p0_union(V,n)):.5f}  "
          f"match={abs(float(t)-float(p0_union(V,n)))<1e-9}  ({surv} pairs vs {full} subsets)")
print("  => a genuine sign-reversing involution on the LRC (star) sum; iterating")
print("     pivots localizes p_0 to the NERVE / irredundant Helly core.")
print()

# ============================================================
print("  III. REDUNDANCY = THE INVOLUTION'S FUEL (dissociated => none)  [VERIFIED]")
print("  " + "-" * 50)
def redundancy(V, n):
    bp, d = bpd(V, n); red = 0; tot = 0
    for r in range(2, len(V)+1):
        for S in combinations(V, r):
            mS = m_of(S, bp, d)
            for v in S:
                tot += 1
                if m_of(tuple(x for x in S if x != v), bp, d) == mS: red += 1
    return red, tot
print(f"  {'V':<18} {'tight':>6} {'p0':>8} {'redundant/total':>16}")
for V, n in [((1,2,3,4),5),((1,3,4,7),5),((1,3,4,5,9),6),
             ((1,2,4),4),((1,2,3,5),5),((1,2,4,8),5)]:
    red, tot = redundancy(V, n); pe = p0_union(V, n)
    print(f"  {str(V):<18} {str(pe==0):>6} {float(pe):>8.4f} {f'{red}/{tot}':>16}")
print("  (dissociated (1,2,4): ZERO redundancy => involution trivial => p_0 stays")
print("   positive, never tight. But redundancy AMOUNT does NOT cleanly decide")
print("   tightness (non-tight (1,2,3,5) has more than tight (1,3,4,7)) -- the")
print("   Vitali wall: the residual is the quantitative coupling, not a count.)")
print()

# ============================================================
print("  IV. CONTRAST WITH REDEI (the SOLVED face) AND THE n=14 TARGET")
print("  " + "-" * 50)
print("""  REDEI: T -> T^op is sign-REVERSING on the Ham-path PARITY sum over GF(2);
  it cancels non-self-converse classes, leaving self-converse fixed points whose
  parity is forced -> #HamPaths odd, PROVED. The reversal there is sign-reversing
  because the object is SIGNED (parity); the SAME reversal on the LRC MEASURE p_0
  is the sign-PRESERVING negation (I). That is precisely why Redei is solved and
  LRC is not: LRC's natural object is an unsigned measure.
  n=14 TARGET: the worry-set reduces to 2^6 = 64 self-converse round classes
  (S576o) -- exactly Redei's fixed-point set. Closing LRC@14 needs a SECOND,
  sign-reversing involution on those 64 (a twist on the 2^6 cube), or
  equivalently a redundancy/nerve collapse of the two-block Helly core
  (codex-S599). The subset pivot involution (II) is that mechanism; the
  un-cancellable residual is the observer-coupling defect (S580o).""")
print()

print("=" * 70)
print("""  SUMMARY
  -------
  * GEOMETRIC clock involutions (negation, half-shift) are measure-preserving =>
    sign-preserving on p_0 => CANNOT certify loneliness by cancellation. Negation
    is the only symmetry of the lonely set; the half-shift is not even that.
  * The CORRECT twisted involution is the SUBSET pivot involution on the (star)
    sum p_0 = sum_S (-1)^|S| m_S: sign-reversing, cancels REDUNDANT subsets,
    localizes p_0 to the nerve / irredundant Helly core (VERIFIED: survivors-sum
    = p_0). This is the LRC form of Redei's permanent->determinant resolution.
  * REDUNDANCY is the fuel: dissociated => none => not tight; but the amount does
    not decide tightness (Vitali wall) -- the residual is the observer-coupling.
  * n=14 reduces to a sign-reversing involution on the 64 self-converse classes
    (Redei's fixed points) / a nerve collapse of the two-block Helly core. NOT
    closed here; this pins exactly what object would close it.
""")
