"""
mac-mini-2026-07-07-S45 (HYP-4947) -- THE WEIGHTED CHERRY LAW (THM-641 candidate):
exact closed form for the fully-weighted, fully-anchored pair mass

  M(q1,q2; alpha,beta; a,b) = meas{x in [0,1): frac(q1 x) in (a, a+alpha],
                                                frac(q2 x) in (b, b+beta]},
  gcd(q1,q2)=1 (general gcd reduces by substitution).

DERIVATION (verified below against brute-exact BEFORE any assertion):
  Pull back both constraints: x-intervals of lengths A=alpha/q1, B=beta/q2 at starts
  (m1+a)/q1, (m2+b)/q2.  With gcd=1 the relative offsets sweep (j+phi)/N, j mod N,
  N = q1 q2, phi = q1 b - q2 a, EXACTLY once (CRT/Bezout).  So
     M = sum_j ov((j+phi)/N),   ov = the circular overlap trapezoid of lengths A,B.
  Periodizing ov via the second Bernoulli polynomial (B2bar(t) = frac(t)^2 - frac(t) + 1/6,
  (B2bar/2)'' = delta - 1 on the circle) and applying the multiplication theorem
  sum_{j mod N} B2bar((t+j)/N) = B2bar(t)/N gives a FOUR-TERM closed form:

     M = alpha*beta + [ s1*B2bar(phi + q1*beta) + s2*B2bar(phi) +
                        s3*B2bar(phi + q1*beta - q2*alpha) + s4*B2bar(phi - q2*alpha) ] / (2 q1 q2)

  with signs (s1..s4) = (+1,-1,-1,+1) up to the offset-direction convention -- the
  script RECONCILES the convention against brute-exact and prints the verified form.

CHECKS: (i) brute-exact (finite Fraction sum over (m1,m2)) == closed form, sweeping
  q <= 24, widths in {1/14, 1/7, 3/14, 2/7}, anchors in {0, 1/2, 1/3}; (ii) at
  (theta,theta,0,0) reproduce THM-638's same-sign law theta^2 + min(r)(s-max(r))/(s^2 q1q2)
  exactly, all residue classes; (iii) the mixed-sign law via anchor a = -alpha (window
  (-alpha, 0]).
APPLICATIONS (this script, after verification):
  (iv) the {0,1/2}-anchored pair table (the 2-anchor tail's exact pair inputs);
  (v) does WIDTH-SWEPT pair data beat klein-S159's 0.1233 pairwise barrier at their
      barrier shape {1,3,5,6,7,9,25,38}?  (multi-width LP, 3 bins/speed.)
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import numpy as np

def B2bar(t: F) -> F:
    ft = t - (t.numerator // t.denominator)   # frac for Fraction (works for negatives)
    return ft*ft - ft + F(1,6)

def brute_exact(q1,q2,al,be,a,b):
    """finite exact sum over (m1,m2) of interval overlaps (circle handled by shifts)."""
    tot = F(0)
    ints1 = [((F(m1)+a)/q1, (F(m1)+a+al)/q1) for m1 in range(q1)]
    ints2 = [((F(m2)+b)/q2, (F(m2)+b+be)/q2) for m2 in range(q2)]
    def clip(lo,hi):
        segs=[]
        for base in (-1,0,1):
            l,h = lo+base, hi+base
            l2,h2 = max(l,F(0)), min(h,F(1))
            if l2<h2: segs.append((l2,h2))
        return segs
    seg1=[]; seg2=[]
    for lo,hi in ints1: seg1+=clip(lo,hi)
    for lo,hi in ints2: seg2+=clip(lo,hi)
    for l1,h1 in seg1:
        for l2,h2 in seg2:
            l,h = max(l1,l2), min(h1,h2)
            if l<h: tot += h-l
    return tot

def closed_form(q1,q2,al,be,a,b,signs=(1,-1,-1,1),phisign=1):
    N = q1*q2
    phi = phisign*(q1*b - q2*a)
    terms = [B2bar(phi + q1*be), B2bar(phi), B2bar(phi + q1*be - q2*al), B2bar(phi - q2*al)]
    corr = sum(s*t for s,t in zip(signs,terms))
    return al*be + corr/(2*N)

# ---------- reconcile the convention ----------
print("=== convention reconciliation (try sign patterns / phi orientation) ===")
tests = [(3,5,F(1,7),F(1,7),F(0),F(0)), (4,7,F(1,7),F(2,7),F(0),F(1,2)),
         (5,9,F(3,14),F(1,7),F(1,3),F(0)), (2,3,F(1,7),F(1,7),F(1,2),F(1,2)),
         (6,7,F(2,7),F(3,14),F(0),F(1,3)),
         (2,9,F(3,14),F(1,14),F(1,2),F(1,3)),   # phi-sign-forcing (asymmetric anchors)
         (2,13,F(1,7),F(1,7),F(1,3),F(1,2))]
best=None
for phisign in (1,-1):
    for signs in [(1,-1,-1,1),(-1,1,1,-1)]:
        ok = all(brute_exact(*t)==closed_form(*t,signs=signs,phisign=phisign) for t in tests)
        print(f"  signs={signs} phisign={phisign}: {'MATCH on 5 probes' if ok else 'no'}")
        if ok: best=(signs,phisign)
assert best, "no convention matched -- derivation bug"
SIGNS, PHISIGN = best
print(f"verified convention: signs={SIGNS}, phi={'+' if PHISIGN==1 else '-'}(q1*b - q2*a)")

# ---------- (i) broad exact verification ----------
print("\n=== (i) broad verification: q<=24 coprime pairs x widths x anchors ===")
import random
random.seed(45)
fails=0; n=0
widths=[F(1,14),F(1,7),F(3,14),F(2,7)]
anchors=[F(0),F(1,2),F(1,3)]
for q1 in range(2,25):
    for q2 in range(q1+1,25):
        if gcd(q1,q2)!=1: continue
        al,be = random.choice(widths), random.choice(widths)
        a,b = random.choice(anchors), random.choice(anchors)
        n+=1
        if brute_exact(q1,q2,al,be,a,b)!=closed_form(q1,q2,al,be,a,b,SIGNS,PHISIGN):
            fails+=1
            if fails<4: print(f"  FAIL {q1},{q2},{al},{be},{a},{b}")
print(f"  {n} random cases: {fails} failures")

# ---------- (ii) THM-638 reproduction ----------
print("\n=== (ii) reproduce THM-638 same-sign law at (theta,theta,0,0), theta=1/7 ===")
th=F(1,7); s=7
bad=0; cnt=0
for q1 in range(1,41):
    for q2 in range(q1+1,41):
        if gcd(q1,q2)!=1: continue
        cnt+=1
        r1,r2 = q1%s or s, q2%s or s
        klein = th*th + F(min(r1,r2)*(s-max(r1,r2)), s*s*q1*q2)
        mine = closed_form(q1,q2,th,th,F(0),F(0),SIGNS,PHISIGN)
        if mine!=klein: bad+=1
print(f"  {cnt} coprime pairs q<=40: {bad} mismatches vs THM-638")

# ---------- (iii) mixed-sign via anchor a = -alpha ----------
print("\n=== (iii) mixed-sign (window (-theta,0] x (0,theta]) via anchoring ===")
bad=0; cnt=0
for q1 in range(1,31):
    for q2 in range(q1+1,31):
        if gcd(q1,q2)!=1: continue
        cnt+=1
        r1,r2 = q1%s or s, q2%s or s
        klein_mixed = th*th - F(min(r1*r2,(s-r1)*(s-r2)), s*s*q1*q2)
        mine = closed_form(q1,q2,th,th,-th,F(0),SIGNS,PHISIGN)
        truth = brute_exact(q1,q2,th,th,-th,F(0))
        if mine!=truth: bad+=1
        if truth!=klein_mixed: cnt2=globals().setdefault('_kx',0); globals()['_kx']=cnt2+1
print(f"  {cnt} pairs: {bad} closed-vs-brute mismatches (THE theorem test); "
      f"{globals().get('_kx',0)} brute-vs-klein-encoding diffs (convention, flagged not asserted)")

# ---------- (iv) the {0,1/2} anchored table (2-anchor tail inputs) ----------
print("\n=== (iv) {0,1/2}-anchor pair table (theta=1/7 windows), sample ===")
print(f"{'(q1,q2)':>9s} {'(0,0)':>12s} {'(0,1/2)':>12s} {'(1/2,1/2)':>12s}  (all exact; corr x 2N shown)")
for (q1,q2) in [(1,2),(2,3),(3,4),(5,7),(6,7),(9,13),(10,13),(11,13),(12,13)]:
    if gcd(q1,q2)!=1: continue
    row=[]
    for (a,b) in [(F(0),F(0)),(F(0),F(1,2)),(F(1,2),F(1,2))]:
        M = closed_form(q1,q2,th,th,a,b,SIGNS,PHISIGN)
        row.append(2*q1*q2*(M-th*th))
    print(f"({q1:2d},{q2:2d})  {str(row[0]):>12s} {str(row[1]):>12s} {str(row[2]):>12s}")

# ---------- (v) width-swept pair LP vs the 0.1233 barrier ----------
print("\n=== (v) width-swept pair-marginal LP at klein-S159's barrier shape ===")
try:
    from scipy.optimize import linprog
    E=[1,3,5,6,7,9,25,38]; k=len(E)
    # bins per speed phase: (0,1/14], (1/14,1/7], (1/7,1]  -> 3 states; atom = state vector
    # constraints: per-speed bin probs (exact: 1/14,1/14,6/7); per-pair joint bin probs
    # from the weighted law: P(phase_i in (0,w1], phase_j in (0,w2]) for w in {1/14,1/7}
    # -> 2x2 cumulative table per pair -> 3x3 joint bin table per pair (inclusion-exclusion).
    def pair_cum(qi,qj,w1,w2):
        g=gcd(qi,qj)
        return closed_form(qi//g,qj//g,w1,w2,F(0),F(0),SIGNS,PHISIGN)
    W=[F(0),F(1,14),F(1,7)]
    nA=3**k
    A_eq=[]; b_eq=[]
    # normalization
    A_eq.append([1.0]*nA); b_eq.append(1.0)
    def atom_states(idx):
        s=[]
        for _ in range(k):
            s.append(idx%3); idx//=3
        return s
    states=[atom_states(i) for i in range(nA)]
    # per-pair joint bin constraints
    for (i,j) in combinations(range(k),2):
        cum={}
        for wi in (1,2):
            for wj in (1,2):
                cum[(wi,wj)] = float(pair_cum(E[i],E[j],W[wi],W[wj]))
        for bi in (0,1):
            for bj in (0,1):
                # P(state_i == bi, state_j == bj) via inclusion-exclusion on cumulatives
                p = cum[(bi+1,bj+1)]
                if bi>0: p-=cum[(bi,bj+1)]
                if bj>0: p-=cum[(bi+1,bj)]
                if bi>0 and bj>0: p+=cum[(bi,bj)]
                row=[1.0 if (states[t][i]==bi and states[t][j]==bj) else 0.0 for t in range(nA)]
                A_eq.append(row); b_eq.append(p)
    # objective: minimize P(no gap>1/7 certificate...) -- klein's LP minimized the
    # empty-atom mass; here: minimize P(all speeds in state 2) LOWER bound proxy for
    # mu via the W-event: mu >= P(exists empty window) >= ... use klein's convention:
    # minimize mass of atoms where NO speed is in bin0 (their 'empty-atom' = the
    # all-far pattern certifying a gap at the anchor).
    obj=[1.0 if all(st==2 for st in states[t]) else 0.0 for t in range(nA)]
    res=linprog(obj, A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                bounds=[(0,None)]*nA, method='highs')
    if res.success:
        print(f"  2-width (3-bin) pair LP floor at barrier shape: {res.fun:.5f} "
              f"(single-width barrier 0.1233; Hunter 6/49 = {6/49:.5f})")
    else:
        print(f"  LP infeasible/failed: {res.message} (grid pair data may be inconsistent -- report)")
except ImportError:
    print("  scipy not available -- LP deferred")
