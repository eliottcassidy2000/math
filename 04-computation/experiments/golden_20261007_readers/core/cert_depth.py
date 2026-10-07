"""K(n) = least depth s at which n's class mod 2^s carries a class-decided certificate (descent, or backward-tree branch
with b <= a_s and ratio < 1) that applies to n itself (explicit m < n, m >= 1, T^i(m) = T^s(n)), vs the Terras stopping
time sigma(n) = least t with T^t(n) < n.  Exact integer arithmetic."""
import math, sys
def T(x): return x//2 if x%2==0 else (3*x+1)//2
def lt1(e2,e3):
    if e2<=0 and e3<=0: return not (e2==0 and e3==0)
    if e2>=0 and e3>=0: return False
    if e2<0: return 3**e3 < 2**(-e2)
    return 2**e2 < 3**(-e3)
def bsearch(v, n, e2, e3, p, budget):
    # v: actual integer node value; p: remaining odd budget (b <= a_s); ratio relative to n = 2^e2 3^e3
    if lt1(e2,e3) and 1 <= v < n: return v
    if not lt1(e2+p, e3-p) or budget[0] <= 0: return None
    budget[0] -= 1
    if p > 0 and v % 3 == 2:
        r = bsearch((2*v-1)//3, n, e2+1, e3-1, p-1, budget)
        if r is not None: return r
    if v % 3 == 0 and p > 0:  # no odd preimages ever: only doubling, ratio grows
        return None
    return bsearch(2*v, n, e2+1, e3, p, budget)
def K_of(n, smax=400):
    x = n; a = 0
    for s in range(1, smax+1):
        a += x & 1; x = T(x)
        if 3**a < 2**s and x < n: return s, 'descent'
        budget=[200000]
        m = bsearch(x, n, -s, a, a, budget)
        if m is not None and m != x: return s, f'branch m={m}'
    return None, ''
def sigma(n):
    x=n; t=0
    while True:
        x=T(x); t+=1
        if x<n: return t
records=[27,255,447,639,703,1819,4255,4591,9663,20895,26623,31911,60975,77671,113383,138367,159487,270271,665215,704511,1042431,1212415,1441407,1875711,1988859,2643183,2684647,3041127,3873535,4637979,5656191,6416623,6631675,19638399,38595583,80049391,120080895,210964383,319804831,1410123943]
for n in records:
    s=sigma(n); K,how=K_of(n, s)
    print(f"n={n:>11d} log2n={math.log2(n):5.1f} sigma={s:4d} K(n)={K} {how}  K/sigma={K/s:.2f}  K/log2n={K/math.log2(n):.2f}", flush=True)
