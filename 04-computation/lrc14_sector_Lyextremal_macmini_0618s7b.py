#!/usr/bin/env python3
"""
lrc14_sector_Lyextremal_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B)

THM-534 reduced the dangerous rows k=8,9,10 to ONE scalar extremal statement:
   "consec maximizes L_y(E)",   L_y(E) = sum_r y_r S_r(E),   S_r(E) = sum_{|A|=r} J(A,E).
THM-534's optimal duals (integer-root g(t)):
   k=11,12,13 (R=2): L = 1 - (1/2)S_1 + (1/6)S_2,  g(t)=(t-3)(t-4)/12
   k=8,9,10  (R=?): need the exact y. From canon table only R=2 row shown; reconstruct R higher.

GOAL (ANGLE B contribution): attack "consec maximizes L_y" via the cutting-word picture.
Since L_y is a POSITIVE combination only if all y_r>=0 (NOT the case: y_1<0!), per-r extremality
is NOT enough. But test: is each S_r(E) extremized by consec in the RIGHT direction so that
the signed combo L_y is maximized?  And does the Sturmian partial-sum picture explain it?

ALSO test the cleaner Bonferroni route: meas(S7) <= S_0 - S_1 + S_2 (degree-2 Bonferroni,
the k=11,12,13 certificate). Is consec the maximizer of (1 - S_1 + S_2/... )? Reconstruct
the exact L_y for each k by solving the moment-LP dual, then test consec-extremality on the
Sturmian-reparametrized AP vs all bounded-spread competitors.

THE KEY ANGLE-B IDEA: S_r(E) = E_x[ C(N(x), r) ], N(x) = # missed sectors among {1..6}.
For the AP (consec), N(x) is governed by the Sturmian partial-sum walk: at theta=7x in [j,j+1),
the residues hit are partial sums of a slope-frac(theta) mechanical word. The NUMBER of distinct
residues among {floor(e theta) mod 7} = 7 - N. So N(x) = 7 - #distinct partial sums.
For a general E, the residues are a SUBSET (at indices e in E) of a longer walk, so FEWER
distinct => MORE missed => LARGER N pointwise IF E subset {0..maxE} (pointwise domination again).
=> S_r(E) >= ... hmm sign. Let's just COMPUTE and find the structural truth.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
sys.stdout.reconfigure(line_buffering=True)

def factorial_moments(E, R=6):
    """S_r(E) = sum_{|A|=r, A subset {1..6}} J(A,E), via one breakpoint pass.
       At each elementary x-interval, free = # of sectors among {1..6} NOT hit; contributes
       length * C(free, r) to S_r."""
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps)
    S=[F(0)]*(R+1)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; L=x1-x0
        hit=set(int(7*e*xm)%7 for e in E)
        # sector 0 always hit (e=0). free among {1..6}:
        free = 6 - len([s for s in hit if s!=0])
        for r in range(R+1):
            S[r]+=L*comb(free,r)
    return S  # S[0]=1

def measS7_from_moments(S):
    # meas(S7)=p_0 = sum_t (-1)^t S_t  (inclusion-exclusion: p_0 = sum_t (-1)^t C(t,t) ... )
    # Actually p_0 = sum_{t} (-1)^t S_t? Check: meas(S7)=sum_{A}(-1)^|A| J(A)=sum_r (-1)^r S_r.
    return sum((-1)**r*S[r] for r in range(len(S)))

# THM-534 dual functionals (reconstruct by moment-LP). For degree R, find y minimizing
# the worst-case bound, but here just USE the canon degree-2 and try a few degrees.
def Ly_deg2(S):  # 1 - (1/2)S_1 + (1/6)S_2, g(t)=(t-3)(t-4)/12 >=0 on t=0..6, g(0)=1
    return F(1) - F(1,2)*S[1] + F(1,6)*S[2]

def g_deg2(t): return F((t-3)*(t-4),12)

cap = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7), 13:F(1)}

print("="*92)
print("Verify moment identity meas(S7)=sum_r (-1)^r S_r, and the deg-2 Bonferroni bound.")
print("="*92)
for E in [(0,1,2,3,4,5,6,7),(0,1,3,7,12,20,30,44),(0,2,3,4,5,6,7,8)]:
    S=factorial_moments(E)
    s7=measS7_from_moments(S)
    ly=Ly_deg2(S)
    print(f"  E={E}: S=[{', '.join(f'{float(s):.4f}' for s in S)}]")
    print(f"     meas(S7)={float(s7):.5f}, Ly_deg2={float(ly):.5f}, meas(S7)<=Ly? {s7<=ly}")
# confirm g(t)>=1[t=0]
print("  g_deg2(t) for t=0..6:", [str(g_deg2(t)) for t in range(7)], " (need >=0, g(0)=1)")

print()
print("="*92)
print("ANGLE B TARGET: does consec maximize L_y_deg2(E) over bounded-spread competitors?")
print("(the k=11,12,13 certificate; also test if it's a valid bound for k=8,9,10 -- it is a")
print(" valid UPPER bound for all k, just maybe not <=cap_k for the smaller k).")
print("="*92)
def gen(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest; g=0
        for e in E: g=gcd(g,e)
        if g!=1: continue
        out.append(E)
    return out
for k in [8,9,10,11]:
    maxE=k+5
    AP=tuple(range(k)); shapes=gen(k,maxE)
    apLy=Ly_deg2(factorial_moments(AP))
    mx=None
    for E in shapes:
        v=Ly_deg2(factorial_moments(E))
        if mx is None or v>mx[0]: mx=(v,E)
    ck=cap[k]
    print(f"  k={k}: Ly_deg2(AP)={float(apLy):.5f}  max over box={float(mx[0]):.5f} at {mx[1]}  "
          f"consec_is_max={apLy==mx[0]}  Ly(AP)<=cap_k? {apLy<=ck} (cap={float(ck):.4f})")

print()
print("="*92)
print("PER-r factorial moment extremality: which direction does consec extremize each S_r?")
print("L_y = 1 - (1/2)S_1 + (1/6)S_2. To MAXIMIZE L_y: want S_1 MIN, S_2 MAX. Does consec do that?")
print("="*92)
for k in [8,9]:
    maxE=k+5; AP=tuple(range(k)); shapes=gen(k,maxE)
    apS=factorial_moments(AP)
    for r in [1,2,3]:
        vals=[(factorial_moments(E)[r],E) for E in shapes]
        mx=max(vals); mn=min(vals); apv=apS[r]
        print(f"  k={k} S_{r}: AP={float(apv):.5f}  min={float(mn[0]):.5f}  max={float(mx[0]):.5f}  "
              f"AP-is-min={apv==mn[0]} AP-is-max={apv==mx[0]}")
print("\nDONE.")
