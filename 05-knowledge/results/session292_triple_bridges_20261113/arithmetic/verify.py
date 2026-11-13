"""Exact factorized wheel and its least period; standard library only."""
from math import gcd,lcm,isqrt,prod
from fractions import Fraction
from pathlib import Path
import json,argparse
HERE=Path(__file__).resolve().parent

def check(c,msg):
    if not c: raise RuntimeError(msg)

def prime(n): return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def order(a,p):
    x=a%p; d=1
    while x!=1: x=x*a%p; d+=1
    return d

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--full-period',action='store_true'); args=ap.parse_args()
    P=[p for p in range(2,100) if prime(p)]
    active=[]
    for p in P[1:]:
        d=order(2,p)
        R=[r for r in range(d) if (99*pow(2,r,p)+23)%p==0]
        if R: check(len(R)==1,'single log'); active.append((p,R[0],d))
    A=[t for t in active if t[0] not in (29,37)]
    core=[t for t in A if t[0] in (5,13,17,19,61)]
    petals=[t for t in A if t[0] not in (5,13,17,19,61)]
    M=lcm(*(d for p,r,d in A)); C=lcm(*(d for p,r,d in core))
    check((M,C)==(61209720,360),'period/core')
    good=lambda k: all(k%d!=r for p,r,d in A)
    counts=[]; mod8=[0]*8
    for b in range(C):
        n=0
        if all(b%d!=r for p,r,d in core):
            n=prod(d//gcd(d,C)-int((b-r)%gcd(d,C)==0) for p,r,d in petals)
        counts.append(n); mod8[b%8]+=n
    total=sum(counts)
    check(total==32353278,'count')
    witnesses=[]
    for q,k in zip((2,3,5,11,13,29,41),(5,21,0,56,29,93,69)):
        j=k+M//q
        check(good(k)!=good(j),'least-period witness')
        witnesses.append(dict(q=q,k=k,other=j,good=[good(k),good(j)],rejected_by=[[p for p,r,d in A if v%d==r] for v in (k,j)]))
    for k in range(1,10001):
        check(good(k)==all((99*pow(2,k,p)+23)%p for p in P),'direct modular check')
    four=[149,347,2111,15583]
    check(all(prime(p) for p in four) and prod(four)==99*2**34+23,'four-factor control')
    four_orders=[order(2,p) for p in four]
    four_mod=lcm(M,*four_orders)
    check(four_mod==30670763342539320 and good(34),'four-factor cylinder')
    # Every forbidden 29 or 37 class is contained in the forbidden 5 class.
    for p,r,d in active:
        if p in (29,37): check(d%4==0 and r%4==3,'redundancy')
    # Reproduce counts by an independent generalized-CRT inclusion/exclusion.
    density=Fraction(1)
    from itertools import combinations
    def merge(a,m,b,n):
        g=gcd(m,n)
        if (b-a)%g: return None
        return ((a+m*((b-a)//g*pow(m//g,-1,n//g)%(n//g)))%lcm(m,n),lcm(m,n))
    for size in range(1,len(A)+1):
        for sub in combinations(A,size):
            state=(0,1)
            for p,r,d in sub:
                state=merge(*state,r,d)
                if state is None: break
            if state is not None: density+=Fraction((-1)**size,state[1])
    check(density==Fraction(total,M),'independent IE')
    if args.full_period:
        mask=bytearray(b'\1')*M
        for p,r,d in A: mask[r::d]=b'\0'*len(mask[r::d])
        check(sum(mask)==total,'full periodic census')
        check([sum(mask[r::8]) for r in range(8)]==mod8,'full mod8 census')
    # A genuinely pairwise-compatible but jointly empty toy guard system.
    guards=[{0,1},{1,2},{0,2}]
    check(all(guards[i]&guards[j] for i in range(3) for j in range(i)) and not set.intersection(*guards),'triple hostile')
    out=dict(active=active,redundant=[29,37],minimal_period=M,core_period=C,petals=petals,core_fiber_counts=counts,
      survivors=total,density=str(density),survivors_by_mod8=mod8,minimality_witnesses=witnesses,
      direct_checks=10000,full_period_checked=args.full_period,
      factor419_cylinder=dict(residue=2,period=lcm(M,418)),
      four_factor_cylinder=dict(residue=34,period=four_mod,factors=four,orders=four_orders))
    path=HERE/('receipt.full.json' if args.full_period else 'receipt.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='core_fiber_counts'},indent=2))
if __name__=='__main__': main()
