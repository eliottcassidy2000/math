"""Exact bounded-weight attachment recurrence via integer polynomial squaring."""
from collections import Counter
from itertools import product
from math import gcd
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
def mu(n):
    sign=1;p=2
    while p*p<=n:
        if n%p==0:
            n//=p;sign=-sign
            if n%p==0:return 0
            while n%p==0:n//=p
        p+=1
    return -sign if n>1 else sign
def divisors(n):return [d for d in range(1,n+1) if n%d==0]

def determinant_histogram(L):
    """No coefficient carries: a packed polynomial squared as an integer."""
    f=Counter(a*b for a in range(-L,L+1) for b in range(-L,L+1))
    width=max(1,(((2*L+1)**4).bit_length()+7)//8)
    packed=b''.join(f.get(k,0).to_bytes(width,'little') for k in range(-L*L,L*L+1))
    z=int.from_bytes(packed,'little')**2
    raw=z.to_bytes(width*(4*L*L+1),'little')
    return {k-2*L*L:int.from_bytes(raw[k*width:(k+1)*width],'little') for k in range(4*L*L+1)}

def group_count(d,e,L,s,n,cache):
    if n<1 or d%s or e%s:return 0
    h0=gcd(d//s,e//s);answer=0
    for sign in (-1,1):
        t=(d//s)*(e//s)-sign*n
        for h in divisors(h0):
            if t%(h*h):continue
            ell=L//(s*h)
            if ell not in cache:cache[ell]=determinant_histogram(ell)
            answer+=mu(h)*cache[ell].get(t//(h*h),0)
    return answer

def main():
    cache={};gates=0;cases=0
    for L in range(4):
        D=determinant_histogram(L);cache[L]=D
        brute=Counter(a*f-b*c for a,b,c,f in product(range(-L,L+1),repeat=4))
        if {k:v for k,v in D.items() if v}!=dict(brute):raise RuntimeError('packed convolution')
        gates+=1
        for d,e in product(range(1,7),repeat=2):
            groups=Counter()
            for a,b,c,f in product(range(-L,L+1),repeat=4):
                t=d*e-a*f+b*c
                if not t:continue
                s=gcd(d,e,a,b,c,f);groups[s,abs(t)//(s*s)]+=1;cases+=1
            bound=d*e+2*L*L
            for s in divisors(gcd(d,e)):
                for n in range(1,bound//(s*s)+1):
                    if group_count(d,e,L,s,n,cache)!=groups[s,n]:raise RuntimeError((d,e,L,s,n))
                    gates+=1
    special={n:[group_count(3,3,L,1,n,cache) for L in range(8)] for n in (25,27,57)}
    minima={n:next(L for L,c in enumerate(row) if c) for n,row in special.items()}
    if minima!={25:4,27:4,57:6}:raise RuntimeError('minimum native integer weight')
    high=determinant_histogram(50)
    if sum(high.values())!=101**4:raise RuntimeError('large convolution mass')
    out={'direct_nonsingular_attachments':cases,'gates':gates,'special_counts_by_L_0_through_7':special,
         'minimum_L_for_Hm_from_H3_H3':minima,'L50_total_matrices':sum(high.values()),
         'L50_determinants_in_support':sum(v>0 for v in high.values()),
         'method':'Kronecker polynomial squaring, followed by content Mobius inversion; labelled weighted cuts, no clock constraint.'}
    (ROOT/'weighted_counts.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
