#!/usr/bin/env python3
"""
lrc14_roundtour_gap_score_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

THE EXACT ROUND-TOURNAMENT <-> GAP DICTIONARY for the LRC residual.

A tie-free difference-winding tournament T(x) is ROUND: order the k cluster phases
p_(0) < p_(1) < ... < p_(k-1) clockwise on R/Z.  Vertex p_(i) BEATS exactly the points
in the OPEN clockwise semicircle (p_(i), p_(i)+1/2).  So the SCORE of p_(i) equals
    s_i = # points strictly inside (p_(i), p_(i)+1/2).

KEY IDENTITIES we prove/verify here:
  (R1)  score s_i = (number of the other k-1 points lying in the clockwise half-circle
        starting just after p_(i)).  Round tournament <=> "interval/semicircle" out-sets.
  (R2)  THE MAX-GAP / SCORE LINK.  A circular gap g_i = p_(i+1) - p_(i) (cyclically).
        Claim:  maxgap > 1/2  <=>  some vertex has score 0 (a SINK = Condorcet loser),
        equivalently some vertex has score k-1 (a SOURCE = Condorcet winner).
        And: maxgap > 1/2  <=>  T(x) is NOT strongly connected (has a dominant vertex).
        Reason: an empty arc of length > 1/2 means the point just clockwise-before it
        sees ALL others in its semicircle (score k-1) and the point just after sees none
        behind => that structure is exactly a transitive split.  TEST EXACTLY.
  (R3)  THE 1/7 LINK (the LRC-relevant one).  maxgap > 1/7 is NOT a pure score event
        (1/7 < 1/2), but it IS a "k-point semicircle-occupancy" event.  We test the
        SHARP reformulation:  for a round tournament, define the SEVENTH-DEGREE vector:
        d^(7)_i = # points in the clockwise 1/7-arc just after p_(i).  Then
            maxgap > 1/7  <=>  some i has d^(7)_i = 0  (an empty 1/7-arc after some point).
        This is a "local sink at scale 1/7" — the LRC 1/7-gap event = existence of a
        SCALE-1/7 LOCAL DOMINATOR.  Verify exactly and tie to OCF/odd-cycle structure.
  (R4)  ODD-CYCLE / REDEI CHECK.  H(T(x)) is always odd (Redei).  We also compute, for
        the round tournament, the number of directed 3-cycles c3 and verify the round
        formula  c3 = C(k,3) - sum_i C(s_i, 2)  (standard: 3-cycles = total triples minus
        transitive triples; transitive triple has a unique 'dominator' counted by C(s_i,2)).
        Then E_x[c3] and E_x[ #odd-cycles-through-fixed-vertex ] vs sector/mu quantities.

Exact Fraction throughout.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb
sys.stdout.reconfigure(line_buffering=True)

def round_tournament(E, x):
    n = len(E); A = [[0]*n for _ in range(n)]; tf = True
    for i in range(n):
        for j in range(n):
            if i == j: continue
            rel = ((E[i]-E[j]) * x) % 1
            if 0 < rel < F(1,2): A[i][j] = 1
            elif rel > F(1,2):   A[i][j] = 0
            else:                A[i][j] = 1 if E[i] < E[j] else 0; tf = False
    return A, tf

def scores(A): return [sum(r) for r in A]
def c3_count(A):
    n=len(A); c=0
    for a,b,cc in itertools.combinations(range(n),3):
        if (A[a][b] and A[b][cc] and A[cc][a]) or (A[a][cc] and A[cc][b] and A[b][a]): c+=1
    return c

def circle_gaps(E, x):
    ps = sorted(set((e*x)%1 for e in E))
    gaps=[]
    for i in range(len(ps)):
        nxt = ps[(i+1)%len(ps)] + (F(1) if i+1==len(ps) else F(0))
        gaps.append(nxt-ps[i])
    return ps, gaps

def breakpoints(E):
    bps=set([F(0),F(1)]); Es=sorted(set(E))
    for e in Es:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
        for m in range(0,2*e+1): bps.add(F(m,2*e))
    diffs=set()
    for a in range(len(Es)):
        for b in range(a+1,len(Es)):
            diffs.add(Es[b]-Es[a]); diffs.add(Es[a]+Es[b])
    for d in diffs:
        if d==0: continue
        for m in range(0,2*d+1): bps.add(F(m,2*d))
    return sorted(x for x in bps if 0<=x<=1)

# ---------- VERIFY identities R2, R3, R4 exactly on sampled cells ----------
def verify_identities(E, label):
    k=len(E); bps=breakpoints(E)
    n_cells=0
    bad_R2_maxgap_half_vs_extreme_score=0
    bad_R3=0
    bad_c3formula=0
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        A,tf=round_tournament(E,xm)
        if not tf: continue
        s=scores(A)
        ps,gaps=circle_gaps(E,xm)
        if len(ps)<k:  # degenerate phase collision at midpoint; skip
            continue
        mg=max(gaps)
        n_cells+=1
        # R2: maxgap>1/2 <=> exists score 0 (and exists score k-1)
        has_sink=(0 in s); has_src=((k-1) in s)
        lhs=(mg>F(1,2))
        if lhs != has_sink or lhs != has_src:
            bad_R2_maxgap_half_vs_extreme_score+=1
        # R3: maxgap>1/7 <=> exists empty 1/7-arc just after some point
        # d7_i = #points in clockwise 1/7 arc after p_(i)
        empty17=False
        P=sorted(set((e*xm)%1 for e in E))
        for a in range(len(P)):
            lo=P[a]
            cnt=0
            for b in range(len(P)):
                if b==a: continue
                d=(P[b]-lo)%1
                if F(0)<d<=F(1,7): cnt+=1
            if cnt==0: empty17=True; break
        if (mg>F(1,7)) != empty17:
            bad_R3+=1
        # R4: c3 = C(k,3) - sum_i C(s_i,2)
        c3_formula=comb(k,3)-sum(comb(si,2) for si in s)
        if c3_count(A)!=c3_formula:
            bad_c3formula+=1
    print(f"  {label:<22} cells={n_cells:>5}  R2_bad={bad_R2_maxgap_half_vs_extreme_score}"
          f"  R3_bad={bad_R3}  c3formula_bad={bad_c3formula}")
    return n_cells

print("="*92)
print("EXACT VERIFICATION of the round-tournament <-> gap dictionary (tie-free cells)")
print("="*92)
print("\nR2: maxgap>1/2 <=> (exists sink score 0) <=> (exists source score k-1)")
print("R3: maxgap>1/7 <=> (exists empty 1/7-arc after some point)  [scale-1/7 local sink]")
print("R4: c3 = C(k,3) - sum_i C(s_i,2)  (round-tournament 3-cycle formula)")
print("-"*92)
shapes = [
  ("k5 consec",[0,1,2,3,4]), ("k5 dissoc",[0,1,3,7,15]),
  ("k6 consec",[0,1,2,3,4,5]),("k6 dissoc",[0,1,3,7,15,31]),
  ("k7 consec",list(range(7))),("k7 dissoc",[0,1,3,7,15,31,63]),
  ("k8 consec",list(range(8))),("k8 perforated",[0,1,2,3,4,5,6,8]),
  ("k8 dissoc",[0,1,3,7,15,31,63,127]),
]
for label,E in shapes:
    verify_identities(E,label)

print("\n" + "="*92)
print("CONSEQUENCE if R3 holds universally:")
print("  mu_{1/7}(E) = P_x[ maxgap > 1/7 ] = P_x[ T(x) has a SCALE-1/7 LOCAL SINK ]")
print("  i.e. some cluster member's clockwise 1/7-arc is empty of other members.")
print("  This is the 'lonely at scale 1/7' tournament event.  The CRUX mu_{1/7}>=thr_k")
print("  becomes: 'a positive fraction of phases have a scale-1/7 local dominator'.")
print("="*92)
print("DONE.")
