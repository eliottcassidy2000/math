#!/usr/bin/env python3
"""
DEFINITIVE exhaustive test of the M3 Walsh-weighted-lead forbidden-class claim,
restricted to PRIMITIVE-REALIZABLE residue tuples (the only legal LRC inputs).

FACTS:
 (1) M3 sign-votes depend ONLY on speeds mod 14 (verified). The only actual-speed dependence is
     reference=slowest's residue and the tie-break by speed order. So every achievable M3 tournament
     for a primitive distinct-speed LRC set arises from an ordered residue-tuple
     (position 0 = slowest = reference; ascending speed positions for tie-break).
 (2) A residue tuple (r_0..r_{n-1}) is realizable by a PRIMITIVE set of distinct positive integers
     (speed_k ≡ r_k mod 14) iff NOT all r_k even AND NOT all r_k ≡ 0 mod 7.
     Reason: 14=2*7; a prime p|14 forces p|speed_k for all t iff p|r_k for all k; primes p∤14 can always
     be dodged. So only p=2 (all even) and p=7 (all ≡0 mod7) can force non-primitivity.
     [Distinctness is free: speeds r_k+14 t_k can always be made distinct by raising some t_k.]
 (3) Position-0 is the slowest. We must also ensure that among the chosen residues, an ASCENDING distinct
     primitive realization exists with position 0 holding residue r_0. This is always possible when (2) holds
     (raise other t_k's to put them above speed_0), EXCEPT we must double-check the ascending/tie-break order
     is honored. We enumerate ALL orderings (full product over positions) so tie-break order is covered.

We enumerate ALL ordered residue-tuples that are primitive-realizable, build the M3 tournament, and check
whether any "forbidden" iso class (by (H,#3cyc,score)) appears.
"""
from fractions import Fraction as F
from math import gcd
from itertools import product, permutations
from functools import reduce
import sys

UNITS=[a for a in range(1,14) if gcd(a,14)==1]
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
FRAC={r: tuple(nrm(F(r*a,14)) for a in UNITS) for r in range(14)}
WPAT={r: tuple(1 if (r*a)%14 in (1,13) else -1 for a in UNITS) for r in range(14)}

def prim_realizable(res):
    if all(r%2==0 for r in res): return False
    if all(r%7==0 for r in res): return False
    return True

def build_adj(res):
    n=len(res); v0=res[0]; w=WPAT[v0]; fr=[FRAC[r] for r in res]
    adj=[[0]*n for _ in range(n)]
    for i in range(n):
        fri=fr[i]
        for j in range(i+1,n):
            frj=fr[j]; D=0
            for k in range(6):
                if fri[k]>frj[k]: D+=w[k]
                elif fri[k]<frj[k]: D-=w[k]
            if D>0: adj[i][j]=1
            elif D<0: adj[j][i]=1
            else: adj[j][i]=1  # tie-break: speed_i<speed_j so x>y false -> y->x
    return adj

def scores(adj,n): return [sum(adj[i][j] for j in range(n)) for i in range(n)]
def c3(adj,n,sc): return n*(n-1)*(n-2)//6 - sum(s*(s-1)//2 for s in sc)
def hpaths(adj,n):
    cnt=0
    for p in permutations(range(n)):
        if all(adj[p[k]][p[k+1]] for k in range(n-1)): cnt+=1
    return cnt
def desc(adj,n):
    sc=scores(adj,n); return (hpaths(adj,n), c3(adj,n,sc), tuple(sorted(sc)))

FORB={
 4:{(5,2,(1,1,2,2))},
 5:{(9,3,(1,1,2,3,3)),(11,4,(1,2,2,2,3)),(13,4,(1,2,2,2,3)),(15,4,(1,2,2,2,3)),(15,5,(2,2,2,2,2))},
}

def realize_primitive_speeds(res_tuple):
    """Find actual ascending distinct primitive speeds whose residues IN ASCENDING SPEED ORDER == res_tuple.
       res_tuple is given in position order = ascending speed order. Returns speeds or None."""
    n=len(res_tuple)
    # greedy ascending; if non-primitive, raise the LAST position by +14 repeatedly (keeps it largest)
    def attempt(bump_last):
        speeds=[]; cur=0
        for idx,r in enumerate(res_tuple):
            c=(r if r>0 else 14)
            while c<=cur: c+=14
            speeds.append(c); cur=c
        speeds[-1]+=14*bump_last
        return speeds
    for bump in range(0,200):
        sp=attempt(bump)
        if len(set(sp))==n and all(s>0 for s in sp) and sp==sorted(sp) \
           and [s%14 for s in sp]==list(res_tuple) and reduce(gcd,sp)==1:
            return sp
    # try raising an interior odd-residue position instead
    for pos in range(n):
        for bump in range(1,200):
            speeds=[]; cur=0
            for idx,r in enumerate(res_tuple):
                c=(r if r>0 else 14)
                while c<=cur: c+=14
                if idx==pos: c+=14*bump
                speeds.append(c); cur=c
            if len(set(speeds))==n and speeds==sorted(speeds) \
               and [s%14 for s in speeds]==list(res_tuple) and reduce(gcd,speeds)==1:
                return speeds
    return None

def run(n):
    forb=FORB[n]; classes={}; found=[]; total=0; skipped_nonprim=0
    for res in product(range(14),repeat=n):
        if not prim_realizable(res):
            skipped_nonprim+=1; continue
        total+=1
        adj=build_adj(list(res)); d=desc(adj,n)
        classes[d]=classes.get(d,0)+1
        if d in forb:
            found.append((res,d))
    print(f"=== n={n}: PRIMITIVE-REALIZABLE ordered residue-tuples = {total} (skipped non-primitive {skipped_nonprim}) ===")
    print(f"distinct (H,c3,score) classes realized: {len(classes)}")
    for d in sorted(classes): print(f"   {d}: {classes[d]}")
    print(f"FORBIDDEN targets: {sorted(forb)}")
    print(f"FORBIDDEN realized over PRIMITIVE inputs? {len(found)>0}")
    for res,d in found[:10]:
        sp=realize_primitive_speeds(res)
        print(f"   WITNESS residues={res} -> {d}; concrete primitive speeds={sp}")
        if sp:
            # double-confirm via independent recompute
            chk=desc(build_adj([s%14 for s in sp]),n)
            print(f"      recheck on residues of concrete speeds -> {chk}  (gcd={reduce(gcd,sp)})")
    print(f"MAX H={max(d[0] for d in classes)}, MAX c3={max(d[1] for d in classes)}")
    sys.stdout.flush()
    return classes,found

if __name__=="__main__":
    for n in [4,5]:
        run(n); print()
