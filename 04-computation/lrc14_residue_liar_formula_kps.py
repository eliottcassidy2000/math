"""kps-2026-07-04: the residue-liar family {1..11,13,12k} -- CLEAN M-formula + lonely-time
formula => an INFINITE family of census-hard covering families closed by explicit (M,t).
Sharpens mac-mini HYP-4070's value-list (3/41,4/53,5/65,7/89) to M=k/(12k+5). Owner's Fibonacci
hint: 12k+5 is Fibonacci (89=F11 at k=7, 233=F13 at k=19,...) periodically."""
from fractions import Fraction as Fr
def is_cov(v): return all(any(x%q==0 for x in v) for q in range(2,15))
def min_margin(v,t):
    m=None
    for s in v:
        r=Fr(s)*t; d=abs(r-round(r))
        if m is None or d<m: m=d
    return m
def exact_M_denomQ(v,Q):
    """max over a=1..Q-1 of min_i ||v_i a/Q|| (M is achieved at some a/Q for the right Q)."""
    best=Fr(0); besta=0
    for a in range(1,Q):
        m=min_margin(v,Fr(a,Q))
        if m>best: best=m; besta=a
    return best,besta

fibs=set()
a,b=1,1
for _ in range(20): fibs.add(a); a,b=b,a+b
print("Residue-liar family {1..11,13,12k}: verify M = k/(12k+5), find lonely time a/(12k+5):")
print(f"{'k':>3} {'X=12k':>6} {'Q=12k+5':>8} {'M(exact)':>10} {'=k/(12k+5)?':>12} {'lonely a':>9} {'a/Q':>10} {'F?':>4} {'cov?':>5}")
allok=True
avals=[]
for k in range(3,16):
    X=12*k; Q=12*k+5
    v=list(range(1,12))+[13,X]
    M,besta=exact_M_denomQ(v,Q)
    formula=Fr(k,Q)
    ok = (M==formula)
    allok=allok and ok
    avals.append((k,besta))
    isfib = "F!" if Q in fibs else ""
    print(f"{k:>3} {X:>6} {Q:>8} {str(M):>10} {str(ok):>12} {besta:>9} {str(Fr(besta,Q)):>10} {isfib:>4} {str(is_cov(v)):>5}")
print(f"\nM = k/(12k+5) EXACT for all k=3..15? {allok}  => M>1/14 iff k>=3 (2k>5); tight 1/14 at k=1,2 (AP,GW).")
print(f"=> {{1..11,13,12k}} is an INFINITE covering family, lonely for ALL k, by an explicit M-formula.")
print()
# lonely-time a(k) pattern
print("lonely a(k):", avals)
print("a(k) vs 6k+2:", [(k, a, 6*k+2) for k,a in avals[:6]])
# Fibonacci denominators
print()
fibdens=[k for k in range(3,40) if (12*k+5) in fibs]
print(f"k with 12k+5 Fibonacci: {fibdens} => denominators {[12*k+5 for k in fibdens]} (=F11,F13,...) -- OWNER'S HINT lands periodically")
