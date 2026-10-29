"""Independent flats-by-span census for the all-m affine observer theorem."""
from itertools import combinations
from collections import Counter
from math import factorial
import json

def need(b,m):
    if not b:raise ArithmeticError(m)

def rank(vs):
    span={0};r=0
    for v in vs:
        if v not in span:span|={x^v for x in span};r+=1
    return r

def choose2(k):return 0 if k<2 else ((2**k-1)*(2**(k-1)-1))//3

def main():
    records=[]
    for m in [3,4,5]:
        N=2**m;d=m-2;spaces=set()
        for basis in combinations(range(1,N),d):
            W={0}
            for b in basis:W|={x^b for x in W}
            if len(W)==2**d:spaces.add(frozenset(W))
        flats={frozenset(x^a for x in W) for W in spaces for a in range(N)}
        need(len(flats)==4*choose2(m),'flat count')
        resp={k:Counter(S for B in flats for S in combinations(sorted(B),k)) for k in range(1,5)}
        for k in range(1,5):
            for S in combinations(range(N),k):
                dimension=rank([s^S[0] for s in S[1:]])
                need(resp[k][S]==choose2(m-dimension),'all subset response formula')
        if m>=4:
            hi=choose2(m-2);lo=choose2(m-3);need(hi>lo,'sharp fourth response')
            for o in range(N):
                for a,b in combinations([x for x in range(N) if x!=o],2):
                    candidates=[c for c in range(N) if c not in [o,a,b] and resp[4][tuple(sorted([o,a,b,c]))]==hi]
                    need(candidates==[o^a^b],'unique recovered sum')
            # Independently enumerate Boolean polynomial evaluations.
            generators=[(1<<N)-1]
            monomials=[(i,) for i in range(m)]+list(combinations(range(m),2))
            for monomial in monomials:generators.append(sum(1<<x for x in range(N) if all(x>>i&1 for i in monomial)))
            words={0}
            for g in generators:words|={w^g for w in words}
            minimum=min(w.bit_count() for w in words if w)
            support={frozenset(x for x in range(N) if w>>x&1) for w in words if w.bit_count()==minimum}
            need(minimum==2**d and support==flats,'all minimum RM supports')
        else:need(flats=={frozenset(S) for S in combinations(range(N),2)},'m3 all pairs hostile')
        gl=1
        for i in range(m):gl*=N-2**i
        geom=factorial(N)//(N*gl);need(geom*N*gl==factorial(N),'orbit-stabilizer integer')
        records.append({'m':m,'flats':len(flats),'labelled_geometries':geom,'fixed_length_bits':(geom-1).bit_length()})
    need([r['fixed_length_bits'] for r in records]==[5,26,90],'bit counts')
    print(json.dumps(records,sort_keys=True))
if __name__=='__main__':main()
