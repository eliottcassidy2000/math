# opus-2026-07-17-S383 -- HYP-7680: THE CRITICAL-POINT STRUCTURE OF THE GAP FUNCTION.
#
# NEW ANGLE.  LRC(14) is a max-min:  g(V) = sup_t min_v ||v t||,  and the claim
# is g(V) >= 1/14.  Each t -> ||v t|| is a TRIANGULAR WAVE of period 1/v, rising
# and falling with slope +-v, peaking at (2m+1)/(2v).  The lower envelope
# min_v ||v t|| is therefore piecewise linear, and its LOCAL MAXIMA can only
# occur where either
#   (a) two waves CROSS:  v_i t - a = +-(v_j t - b)  =>  t = k/(v_i -+ v_j), or
#   (b) the minimising wave is at its own PEAK:  t = (2m+1)/(2v).
# So the optimum sits at a BEAT FREQUENCY v_i +- v_j (or a half-period 2v), and
# the candidate denominators number at most 13 + 2*C(13,2) = 169 = 13^2.
# If true this is a finite per-family certificate with a fixed, tiny candidate
# set -- and it would explain THM-1105's finding that gaps sit at arithmetically
# special denominators.  TEST IT EXACTLY.
from fractions import Fraction as F
from itertools import combinations
import random
def frac_norm(x):
    """||x|| = distance to nearest integer, exact."""
    r = x - int(x)
    if r < 0: r += 1
    return min(r, 1-r)
def gap_at(V, t):
    return min(frac_norm(t*v) for v in V)
def candidate_denominators(V):
    D=set()
    for v in V: D.add(2*v)
    for a,b in combinations(V,2):
        if a+b>0: D.add(a+b)
        if abs(a-b)>0: D.add(abs(a-b))
    return sorted(D)
def g_candidate(V):
    """max of min_v ||v t|| over the critical-point candidates."""
    best=F(0); arg=None
    for q in candidate_denominators(V):
        for p in range(1, q):
            t=F(p,q)
            g=gap_at(V,t)
            if g>best: best, arg = g, t
    return best, arg
# ---- exact truth via the uncovered set at a given radius ----
def teeth(x, lam):
    w=lam/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered_measure(V, lam):
    allv=[]
    for x in V: allv.extend(teeth(x,lam))
    return 1-sum(b-a for a,b in union(allv))

print("(1) IS THE OPTIMUM ALWAYS AT A CRITICAL POINT?  (candidates: 2v and v_i +- v_j)")
print("    A candidate value g is the true sup iff uncovered(V, lam) is EMPTY for")
print("    every lam > g.  Test at lam = g + tiny.")
random.seed(383)
bad=0; n=0
for trial in range(14):
    V=sorted(random.sample(range(1,60),13))
    g,arg = g_candidate(V)
    eps=F(1,10**7)
    m_at   = uncovered_measure(V, g)          # at g: should be 0 (boundary) or >0
    m_above= uncovered_measure(V, g+eps)      # above g: MUST be 0 if g is the sup
    ok = (m_above==0)
    n+=1
    if not ok: bad+=1
    print(f"    V[:4]={V[:4]}...  g={g} ({float(g):.6f})  at t={arg}"
          f"   uncovered(g+eps)={m_above}  {'OK' if ok else '<-- g is NOT the sup'}")
print(f"    families where a critical point achieved the sup: {n-bad}/{n}")
print()
print(f"(2) CANDIDATE-SET SIZE: at most 13 + 2*C(13,2) = {13+2*len(list(combinations(range(13),2)))} = 13^2 denominators")
V=sorted(random.sample(range(1,60),13))
print(f"    example family {V[:5]}...: {len(candidate_denominators(V))} distinct candidate denominators")
