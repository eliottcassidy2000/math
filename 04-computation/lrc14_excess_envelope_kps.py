import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def N_at(E,x):
    hit=set()
    for e in E:
        v=e*x; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
    return sum(1 for j in range(1,7) if j not in hit)
def dist_p(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); p=[Fraction(0)]*7
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        p[N_at(E,(lo+hi)/2)]+=hi-lo
    return p
def g_poly(k,t):
    if k==8: return Fraction((t-1)*(t-2)*(t-4)*(t-5),40)
    if k in(9,10): return Fraction(-(t-2)*(t-3)*(t-6),36)
    return Fraction((t-3)*(t-4),12)
def L_y(E,k): p=dist_p(E); return sum(p[t]*g_poly(k,t) for t in range(7))
def relrank(E):
    # rank of integer relations sum m_e e=0 with sum m_e=0 (affine), over the k generators.
    # = k - 1 - (rank of {e_i - e_0} over Q) ... but over Q all integers are rank1, so use
    # the ADDITIVE-energy proxy: count independent height-2 relations e_a+e_b=e_c+e_d.
    # Practical proxy: Freiman 'dimension' via |E+E|. We report excess instead.
    E=set(E); return len({a+b for a in E for b in E})-(2*len(E)-1)
caps={8:0.38153,9:0.49426,10:0.6044}
# TEST: is L_y - L_y^inf monotone in relation density? Use excess as inverse-density proxy.
# Compute correlation of L_y with excess across many sets; and the MAX L_y at each excess (envelope).
for k in [8,9]:
    Lyinf={8:0.04927,9:0.15129}[k]; C=list(range(k)); Lc=L_y(C,k)
    print(f"k={k}: L_y(consec)={float(Lc):.4f} L_y^inf={Lyinf} cap={caps[k]}")
    env={}  # excess -> (maxL_y, set)
    box={8:14,9:13}[k]; cnt=0
    for tail in itertools.combinations(range(1,box+1),k-1):
        E=(0,)+tail
        if reduce(gcd,E)!=1: continue
        cnt+=1; ex=relrank(E); L=L_y(E,k)
        if ex not in env or L>env[ex][0]: env[ex]=(L,E)
    print(f"   excess -> max L_y envelope (is it monotone decreasing? the 'dimension penalty'):")
    prev=None; mono=True
    for ex in sorted(env):
        L=env[ex][0]
        tag=""
        if prev is not None and L>prev+Fraction(1,1000): tag=" <-- BUMP UP"; mono=False
        print(f"     exc={ex:2d}: maxL_y={float(L):.4f}")
        prev=L
    print(f"   monotone(strict-ish)? {mono}  scanned {cnt}\n")
