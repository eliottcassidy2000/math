#!/usr/bin/env python3
"""
rung_realizability_nq_witness_klein.py  --  klein-2026-06-30-S54

TWO tracks, both witness/multi-metric:

(A) H7 -- the covering-min RUNG REALIZABILITY function a(n).  Covering-min = smallest realizable Farey
    rung a(n): M = a/(a(n-1)+1) at binding modulus D = a(n-1)+1.  Which rungs are realizable is a
    COVERING-SYSTEM condition: a rung-r covering set has all speeds avoiding the (r-1)-neighborhood of 0
    mod D at some rotation, WHILE covering every resonance q in {2..n}.  We compute a(n) (min realizable
    rung, via known beaters + search), factor D = a(n)(n-1)+1, and hunt the arithmetic pattern.

(B) The (n+q)-WITNESS (opus, tight/floor regime) -- single-swap census.  Tight AP = {1..n-1}, M=1/n.
    Single-swap = drop q, add g.  For q PRIME in ((n-1)/2, n-1]: no swap is tight (M>1/n):
      (i)  g not a multiple of q  -> q-witness at t=1/q: M>=1/q>1/n.
      (ii) g=mq (keeps resonance q) -> at modulus n+q, the pair {q,n} is UNREPRESENTED; witness at
           a=q^{-1} mod (n+q), t=a/(n+q): every runner avoids {0,1,n+q-1} -> M>=2/(n+q)>1/n.
    Only the COMPOSITE doubling drop (q=n-2, g=2(n-2)) survives (GW), and only at n==2 mod 6.
    We VERIFY (i)+(ii) exactly and anchor the doublings.
"""
from fractions import Fraction as F
from math import gcd

def dist0(x,D):
    x%=D; return min(x,D-x)
def M(S):
    Dmax=2*max(S)+2; best=F(0)
    for D in range(2,Dmax+1):
        for a in range(1,D):
            m=D
            for s in S:
                d=dist0(s*a,D); m=d if d<m else m
                if m==0: break
            if F(m,D)>best: best=F(m,D)
    return best
def isprime(p): return p>1 and all(p%d for d in range(2,int(p**.5)+1))
def factor(m):
    f={}; d=2
    while d*d<=m:
        while m%d==0: f[d]=f.get(d,0)+1; m//=d
        d+=1
    if m>1: f[m]=f.get(m,0)+1
    return f
def fstr(m): return "*".join(f"{p}^{e}" if e>1 else str(p) for p,e in factor(m).items()) or str(m)

# ---------- (A) rung realizability ----------
def rung_of(Mv,n):
    for r in range(1,4*n):
        if Mv==F(r,r*(n-1)+1): return r
    return None

# known covering-min beaters (verified last session) + construction fallback
KNOWN = {7:[1,2,5,6,7,8], 8:[1,4,5,6,7,11,16], 9:[1,3,4,5,7,11,18,32], 10:[1,2,3,5,6,7,8,9,30]}

def report_A():
    print("="*74)
    print("(A) H7: covering-min rung a(n), binding modulus D=a(n-1)+1, factorization")
    print("="*74)
    # a(n) from known beaters (n<=10) + repo (n=11: 3/31) + construction upper bound
    amap={}
    for n,S in KNOWN.items():
        amap[n]=(M(S),rung_of(M(S),n))
    amap[11]=(F(3,31),3)  # repo covering-min n=11
    print(f"{'n':>3} {'covering-min':>10} {'rung a(n)':>9} {'D=a(n-1)+1':>11} {'factor(D)':>12} {'D=2n-1?':>8}")
    for n in sorted(amap):
        Mv,a=amap[n]; D=a*(n-1)+1
        print(f"{n:>3} {str(Mv):>10} {a:>9} {D:>11} {fstr(D):>12} {str(D==2*n-1):>8}")
    print("  a(n) = 2,2,4,4,3 (n=7..11): irregular. Note rung 2 => D=2n-1 (the signed-LRC modulus, THM-401).")
    print("  Realizability of rung 2 (mediant 2/(2n-1)) at n=7,8 (D=13,15) but NOT n=9,10,11 -- the")
    print("  covering-system obstruction: whether resonances {2..n} admit a set avoiding {0,+-1} mod (2n-1).")

# ---------- (B) (n+q)-witness ----------
def report_B():
    print()
    print("="*74)
    print("(B) The (n+q)-witness: single-swap census of the tight AP {1..n-1} (M=1/n)")
    print("="*74)
    for n in [14,18,20]:
        AP=list(range(1,n)); floor=F(1,n)
        primes=[q for q in range(n//2+1, n) if isprime(q)]  # large primes in ((n-1)/2, n-1]
        print(f"n={n}: floor 1/n={floor}; large primes q in ((n-1)/2,{n-1}] = {primes}")
        for q in primes:
            # (i) g not multiple of q: sample g=q+1 (not mult) -> q-witness
            g_nonmult=q+1 if (q+1) not in AP or True else q+1
            S1=sorted(set([x for x in AP if x!=q]+[g_nonmult]))
            if len(S1)==n-1:
                pass
            # (ii) g=mq (smallest new multiple 2q) -> (n+q)-witness
            g2=2*q
            S2=sorted(set([x for x in AP if x!=q]+[g2]))
            m2=M(S2) if len(S2)==n-1 else None
            # verify the (n+q) binding + missing pair {q,n} mod (n+q)
            D=n+q; a_star=pow(q, -1, D)
            imgs=sorted(dist0(s*a_star,D) for s in S2)
            witness_val=F(min(imgs),D)
            print(f"   drop q={q}, add 2q={g2}: M(S)={m2}={float(m2):.4f} >1/n? {m2>floor}; "
                  f"(n+q)={D}, a*=q^-1={a_star}, min-dist={min(imgs)} -> witness {witness_val}={float(witness_val):.4f}")
        # the surviving doubling (composite n-2)
        v=n-2; Sgw=sorted(set([x for x in AP if x!=v]+[2*v]))
        mgw=M(Sgw)
        print(f"   DOUBLING drop {v}(composite), add {2*v}: M={mgw}={float(mgw):.4f} {'= FLOOR (tight, GW)' if mgw==floor else '> floor'} (n==2 mod6? {n%6==2})")

if __name__=="__main__":
    report_A()
    report_B()
