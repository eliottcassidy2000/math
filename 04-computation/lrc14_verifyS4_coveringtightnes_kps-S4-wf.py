"""
Adversarial verification of the claimed S3 sub-case closure:
"WINDOW-COLLAPSE LEMMA" + "ALL-MULT7-LARGE sub-case closure".

kind-pasteur S4-wf. EXACT rationals throughout (fractions.Fraction).

Plan:
 (A) Re-derive every load-bearing inequality of the Window-Collapse Lemma exactly.
 (B) Verify the headline validated witness S=[1,2,4,5,6,8,9,10,11,12,13,31,154].
 (C) HUNT for a covering primitive S3 set IN the ALL-MULT7-LARGE sub-case with
     M < 1/14 OR where the hypothesis holds but the claimed global witness fails.
 (D) Stress-test the abstract gap: does a safe u-point ALWAYS survive in (0, 1/(2V*)]?
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r

def safe_components(A, h=F(1,14)):
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j,u); a = (c - h/u) % 1; b = (c + h/u) % 1
            if a < b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged = []
    for a,b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a,b))
    safe = []; prev = F(0)
    for a,b in merged:
        if a > prev: safe.append((prev,a))
        prev = max(prev, b)
    if prev < 1: safe.append((prev,F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc: return F(0)
    ws = [b-a for a,b in sc]
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1:
        ws.append(sc[0][1] + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2): C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b

def is_cov(S):
    return all(any(v % q == 0 for v in S) for q in range(2, 15))

def is_primitive(S):
    g = 0
    for v in S: g = gcd(g, v)
    return g == 1

print("="*70)
print("(A) Re-derive load-bearing inequalities of Window-Collapse Lemma")
print("="*70)

# Claim (a): for tau = k/7 + s, every non-mult-of-7 runner v satisfies
#   ||v*tau|| >= 1/7 - |v|*|s|.
# Proof check: ||v*(k/7 + s)|| = ||v*k/7 + v*s||.
# For 7 not dividing v, v*k/7 mod 1 in {1/7,2/7,3/7,4/7,5/7,6/7}, so ||v*k/7|| in {1/7,2/7,3/7}.
# Min is 1/7. Then ||a+b|| >= ||a|| - ||b|| >= 1/7 - |v*s| = 1/7 - |v||s|.
# Verify the residue rule exactly over k coprime to 7 and v not mult of 7.
print("\n[A1] 7-adic residue rule: ||v*k/7|| in {1/7,2/7,3/7} for 7 not| v, k coprime 7")
ok = True
for k in range(1, 7):
    if gcd(k, 7) != 1: continue
    for v in range(1, 50):
        val = nrm(F(v*k, 7))
        if v % 7 == 0:
            if val != 0: ok = False; print("  FAIL mult7", k, v, val)
        else:
            if val not in (F(1,7), F(2,7), F(3,7)): ok = False; print("  FAIL", k, v, val)
print("  residue rule holds:", ok, "(min nonzero =", F(1,7), ")")

# Claim sharpness: at |s| = 1/(14 V*), v = V*: 1/7 - V*/(14 V*) = 1/7 - 1/14 = 1/14.
print("\n[A2] sharpness bound: 1/7 - V**(1/(14 V*)) =", F(1,7) - F(1,14), "should = 1/14 =", F(1,14))

# Claim (b): for tau = k/7 + s, mult-of-7 runner 7w: ||7w*tau|| = ||7w*k/7 + 7w*s||
#   = ||w*k + w*(7s)|| = ||w*(7s)|| since w*k integer.  Set u = 7s.
print("\n[A3] mult-of-7 reduction: ||7w*(k/7+s)|| = ||w*(7s)||  (since w*k integer)")
ok = True
for k in range(1, 7):
    if gcd(k,7)!=1: continue
    for w in range(1, 30):
        for s in [F(1,100), F(3,200), F(-7,300), F(11,400)]:
            lhs = nrm(F(7*w)*(F(k,7)+s))
            rhs = nrm(F(w)*(7*s))
            if lhs != rhs: ok = False; print("  FAIL", k, w, s, lhs, rhs)
print("  reduction identity holds:", ok)

print("\nWindow: |s| <= 1/(14 V*)  <=>  |u|=|7s| <= 7/(14 V*) = 1/(2 V*).  Confirmed.")
print("So on |u|<=1/(2V*): non-mult-7 runners ALL safe (>=1/14), and the system")
print("reduces to finding u with ||w_i u||>=1/14 for all i, |u|<=1/(2V*).")

print()
print("="*70)
print("[A4] Triangle-ineq subtlety: ||a+b|| >= ||a|| - ||b|| needs ||b||=|v||s|")
print("="*70)
# ||v*s|| = |v*s| only when |v*s| <= 1/2. On the window |s|<=1/(14 V*),
# for v <= V* we have |v*s| <= V*/(14 V*) = 1/14 <= 1/2. OK.
# But the lemma needs it for ALL non-mult-7 v in S, and V* = MAX non-mult-7 runner,
# so every non-mult-7 v <= V*. Hence |v*s| <= 1/14. Good.
# Verify ||v*s|| = |v||s| holds for v<=V* on the window edge.
ok=True
for Vstar in [14,31,100,399]:
    s=F(1,14*Vstar)
    for v in range(1,Vstar+1):
        if v%7==0: continue
        if nrm(F(v)*s) != F(v)*s:
            ok=False; print("  FAIL", Vstar, v, nrm(F(v)*s), F(v)*s)
print("  ||v s|| = |v||s| for all non-mult7 v<=V* on window edge:", ok)

print()
print("="*70)
print("(B) Verify the headline validated witness")
print("="*70)
S = [1,2,4,5,6,8,9,10,11,12,13,31,154]
print("S =", S)
print("  covering:", is_cov(S), " primitive:", is_primitive(S), " size:", len(S))
# S3 check: k=#{v>13}, Vmin, Vmax, is Vmax >= 13 Vmin?
Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
print("  Vmin",Vmin,"Vmax",Vmax,"k=#{v>13}",k,"  S3? (k>=2 and Vmax>=13 Vmin):", k>=2 and Vmax>=13*Vmin)
# mults of 7 in S:
m7=[v for v in S if v%7==0]; nm7=[v for v in S if v%7!=0]
Vstar=max(nm7)
print("  mults of 7:", m7, " V* (max non-mult7):", Vstar)
print("  ALL-MULT7-LARGE? (every mult of 7 > V*):", all(x>Vstar for x in m7), " (7 in S?:", 7 in S, ")")
tau=F(19281,133672)
print("  claimed witness tau =", tau, "=?", F(19281,133672))
mv=min(nrm(F(x)*tau) for x in S)
print("  min_v ||v*tau|| =", mv, "~", float(mv), "  >= 1/14?", mv>=F(1,14))
print("  exact Mval(S) =", Mval(S), "~", float(Mval(S)))

print()
print("="*70)
print("(C) HUNT: covering primitive S3 sets in ALL-MULT7-LARGE with M<1/14")
print("="*70)
# ALL-MULT7-LARGE: every multiple of 7 in S exceeds V* = max non-mult-of-7 runner.
# Equivalent: 7 not in S, and the mult-of-7 runners are the largest elements.
# We need covering -> contains a multiple of 7 AND of 14. So multiples of 14 present.
# Build candidate sets: choose a "small part" of non-mult-of-7 runners (the bulk
# making it covering for q in 2..13 except via 7/14), plus mult-of-7 runners all > V*.
def classify(S):
    S=sorted(set(S))
    if len(S)!=13: return None
    if not is_primitive(S): return None
    if not is_cov(S): return None
    Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
    if not (k>=2 and Vmax>=13*Vmin): return None  # must be S3
    m7=[v for v in S if v%7==0]; nm7=[v for v in S if v%7!=0]
    if not m7 or not nm7: return None
    Vstar=max(nm7)
    allbig = all(x>Vstar for x in m7)  # ALL-MULT7-LARGE
    return dict(Vmin=Vmin,Vmax=Vmax,k=k,m7=m7,nm7=nm7,Vstar=Vstar,allbig=allbig)

import random
random.seed(20260618)
worst=None; tested=0; allbig_tested=0; viol=[]
# Strategy: random covering S3 sets, filter to ALL-MULT7-LARGE, check M>=1/14.
# To get covering we need multiples of all of 2..14. Multiples of 7 and 14 forced large.
def random_allmult7large():
    # small non-mult-7 part: pick from 1..V*  excluding multiples of 7
    Vstar=random.randint(13,120)
    pool=[v for v in range(1,Vstar+1) if v%7!=0]
    if len(pool)<8: return None
    # ensure V* itself in the set (so it's really the max non-mult7)
    nm=set([Vstar])
    while len(nm)<random.randint(8,11):
        nm.add(random.choice(pool))
    nm=sorted(nm)
    if max(nm)!=Vstar:  # enforce
        nm[-1]=Vstar; nm=sorted(set(nm))
    nslots=13-len(nm)
    if nslots<2: return None
    # mult-of-7 runners, all > V*, must include a multiple of 14 for covering
    bigpool=[7*w for w in range(1, 200) if 7*w>Vstar]
    m7=set()
    # force at least one multiple of 14
    pool14=[x for x in bigpool if x%14==0]
    if not pool14: return None
    m7.add(random.choice(pool14))
    while len(m7)<nslots:
        m7.add(random.choice(bigpool))
    if len(m7)!=nslots: return None
    return sorted(nm)+sorted(m7)

for _ in range(400000):
    S=random_allmult7large()
    if S is None: continue
    info=classify(S)
    if info is None: continue
    tested+=1
    if not info['allbig']: continue
    allbig_tested+=1
    M=Mval(S)
    if M<F(1,14):
        viol.append((S,M))
        print("  *** VIOLATION ***", S, "M=",M,"~",float(M))
        if len(viol)>=5: break
    if worst is None or M<worst[1]:
        worst=(S,M)

print(f"  tested covering-primitive-S3: {tested}, of which ALL-MULT7-LARGE: {allbig_tested}")
print("  violations (M<1/14):", len(viol))
if worst:
    print("  worst M in ALL-MULT7-LARGE:", worst[1], "~", float(worst[1]), "  set:", worst[0])

print()
print("="*70)
print("(D) ABSTRACT GAP: does a safe u-point survive in (0, 1/(2V*)] for the")
print("    mult-of-7 subsystem {w_i} with all w_i > V*/7 ?")
print("="*70)
# Reduced subproblem: given V* and a multiset W={w_1,..,w_r} of positive ints
# with each w_i > V*/7 (i.e. 7 w_i > V*), is there u in (0, 1/(2V*)] with
# ||w_i u|| >= 1/14 for all i?  This is the EXACT sub-case-closure claim.
# We can compute the safe set of W on [0, 1/(2V*)] exactly via teeth intervals.

def safe_point_in_window(W, Vstar):
    """Return an exact u in (0,1/(2V*)] with ||w u||>=1/14 all w, or None."""
    R = F(1, 2*Vstar)  # window right endpoint
    h = F(1,14)
    # forbidden intervals on [0,R]: for each w, around j/w radius h/w
    forb=[]
    for w in W:
        j=0
        while True:
            c=F(j,w)
            if c-h/F(w) > R: break
            a=max(F(0), c-h/F(w)); b=min(R, c+h/F(w))
            if a<=b: forb.append((a,b))
            j+=1
    forb.sort()
    # merge
    merged=[]
    for a,b in forb:
        if merged and a<=merged[-1][1]:
            merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else:
            merged.append((a,b))
    # find a gap in (0,R]; we want u>0 strictly (u=0 fails). The first forbidden
    # interval starts at 0 (tooth of every w at j=0). Find first gap after it.
    prev=F(0)
    gaps=[]
    for a,b in merged:
        if a>prev: gaps.append((prev,a))
        prev=max(prev,b)
    if prev<R: gaps.append((prev,R))
    # we need a gap (lo,hi) with hi>0 and a point strictly >0 inside, plus the
    # endpoints: a safe point can be at hi (closed) if ||w hi||>=1/14. Use midpoint
    # if open gap; but exact: pick a rational strictly inside (lo,hi) with lo>=0.
    for lo,hi in gaps:
        if hi<=0: continue
        # pick point: if lo>0 use (lo+hi)/2; if lo==0 we need u>0, still midpoint>0
        u=(lo+hi)/2
        if u<=0: u=hi/2
        if u<=0: continue
        # verify exactly
        if all(nrm(F(w)*u)>=h for w in W):
            return u
    # also boundary points of gaps may be exactly safe (closed teeth boundaries)
    for lo,hi in gaps:
        for u in (lo,hi):
            if u>0 and u<=R and all(nrm(F(w)*u)>=h for w in W):
                return u
    return None

# Exhaustive sweep: V* from 14..399, subsystem sizes r=1..4, all w_i with 7w_i>V*.
# This mirrors the claim's "3.5M adversarial configs, 0 counterexamples". We
# independently reproduce and look for ANY (V*,W) where NO safe point survives.
print("\n[D1] Exhaustive: V*=14..399, |W|<=4, all 7 w_i > V*. Count failures.")
fail_configs=[]
total=0
import itertools
for Vstar in range(14, 400):
    wmin = Vstar//7 + 1  # smallest w with 7w>V*
    # cap w: teeth with period 1/w; only w up to some bound matter on (0,1/(2V*)].
    # If 7w very large, first tooth tiny, easy to dodge. Worst is small w. Cap to
    # keep finite but cover the dangerous small-w regime generously.
    wmax = wmin + 60
    Wpool = list(range(wmin, wmax+1))
    for r in range(1,5):
        for combo in itertools.combinations(Wpool, r):
            total+=1
            if safe_point_in_window(combo, Vstar) is None:
                fail_configs.append((Vstar,combo))
                if len(fail_configs)<=10:
                    print("  *** NO SAFE POINT ***  V*=",Vstar," W=",combo)
    if Vstar%50==0:
        print(f"    ...through V*={Vstar}, configs={total}, failures={len(fail_configs)}")
print(f"  TOTAL configs tested: {total}, configs with NO safe point: {len(fail_configs)}")

print()
print("="*70)
print("[D2] MEASURE ARGUMENT for the ALL-MULT7-LARGE closure (abstract)")
print("="*70)
# Window (0,R], R=1/(2V*). Runner w_i blocks u where ||w_i u||<1/14.
# On (0,R], the number of teeth-centers j/w_i with j>=1 is floor(w_i R). Each tooth
# (interior) has width 2/(14 w_i)=1/(7 w_i). The j=0 tooth at u=0 contributes only
# its right half [0,1/(14 w_i)] inside (0,R].
# Total blocked measure by w_i on (0,R]:
#   <= 1/(14 w_i)  +  floor(w_i R) * 1/(7 w_i)
#   <= 1/(14 w_i)  +  (w_i R) * 1/(7 w_i)  =  1/(14 w_i) + R/7.
# Since w_i > V*/7 => 1/(14 w_i) < 7/(14 V*) = 1/(2 V*) = R. That bound (R) is too
# weak. Use better: 1/(14 w_i) <= R/7 ?  Need w_i >= 2/(R) /14*7... let's just sum.
# SUM over r runners: blocked <= sum_i [1/(14 w_i) + R/7].
# Want < R (total window length) so a safe point survives.
# Each 1/(14 w_i): w_i > V*/7 => 1/(14 w_i) < 1/(14 * V*/7) = 7/(14 V*) = 1/(2V*) = R.
# That's weak. Tighter: smallest possible w_i is wmin = floor(V*/7)+1.
# Let's compute, for each V*, the WORST-CASE total blocked measure UPPER BOUND
# over r runners with distinct w_i >= wmin, and see for which r it stays < R.
def blocked_ub(w, Vstar):
    R=F(1,2*Vstar); h=F(1,14)
    ntooth = (F(w)*R).numerator // (F(w)*R).denominator  # floor(w R), j>=1 teeth
    # right half of j=0 tooth inside window: min(R, h/w)
    half0 = min(R, h/F(w))
    return half0 + ntooth * (F(1,7*w) if F(1,7*w) <= R else R)  # interior tooth width 2h/w=1/(7w)

# For closure we want: for the r SMALLEST allowed w (which maximize blocked measure
# since smaller w -> larger 1/(14w) and 1/(7w)), total blocked < R.
worst_r_needed={}
maxr_safe={}
for Vstar in range(14,400):
    R=F(1,2*Vstar); wmin=Vstar//7+1
    # take r smallest distinct w
    for r in range(1,7):
        ws=list(range(wmin,wmin+r))
        tot=sum(blocked_ub(w,Vstar) for w in ws)
        if tot>=R:
            maxr_safe[Vstar]=r-1
            break
    else:
        maxr_safe[Vstar]=6
mn=min(maxr_safe.values())
print("  min over V* of (max r s.t. measure bound guarantees a safe point):", mn)
# how many mult-of-7 runners can an S3 13-set have? They're all >V* and there are
# at most 13 - (#non-mult7). We need V* present plus enough non-mult7 for covering.
# Show the actual max r encountered: r = #mult-of-7 in S. In ALL-MULT7-LARGE all
# mult7 are the LARGEST. With Vmin small (covering needs small numbers), most
# elements are non-mult7. Empirically r is small. But prove the bound covers r.
print("  V* with measure-guarantee only up to r =", mn, ":",
      [v for v in maxr_safe if maxr_safe[v]==mn][:10], "...")
print("  => measure argument guarantees safe point for all r <=", mn, "at EVERY V*.")
