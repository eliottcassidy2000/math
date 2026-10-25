#!/usr/bin/env python3
"""Literal finite-field root fibres, including quadratic extension fields."""
import hashlib
import itertools as it
import json
from collections import Counter
from pathlib import Path


class Field:
    def __init__(self, p, nonsquare=None):
        self.p, self.d = p, nonsquare
        self.q = p if nonsquare is None else p*p

    def add(self, x, y):
        p=self.p
        return ((x%p+y%p)%p)+p*((x//p+y//p)%p)

    def neg(self, x):
        p=self.p
        return (-x%p)+p*(-(x//p)%p)

    def mul(self, x, y):
        p=self.p
        a,b=x%p,x//p;c,d=y%p,y//p
        return (a*c+(self.d or 0)*b*d)%p+p*((a*d+b*c)%p)

    def power(self, x, n):
        ans=1
        while n:
            if n%2:ans=self.mul(ans,x)
            x=self.mul(x,x);n//=2
        return ans

    def prod(self, values):
        ans=1
        for x in values:ans=self.mul(ans,x)
        return ans


def main():
    checks=0;receipts=[]
    matrix=json.loads((Path(__file__).parent/'normal.json').read_text())['C1']
    for p,ns in [(3,None),(5,None),(7,None),(11,None),(13,None),(17,None),(3,2),(5,2)]:
        f=Field(p,ns);q=f.q
        for x in range(1,q):
            if f.power(x,q-1)!=1:
                raise RuntimeError('field multiplication control')
            checks+=1
        roots=Counter(f.mul(x,x) for x in range(q))
        nonsquare=next(x for x in range(1,q) if roots[x]==0)
        counts=[0]*64;Ucounts=[0]*64;Dcounts=[0]*64;V=0
        for x,z in it.product(range(q),repeat=2):
            rad=[f.add(1,x),f.add(1,f.neg(x)),f.add(1,z),f.add(1,f.neg(z)),
                 f.add(f.add(1,x),z),f.add(f.add(1,f.neg(x)),f.neg(z))]
            if 0 in rad:continue
            V+=1
            for mask in range(64):
                b=[f.mul(v,nonsquare) if mask>>i&1 else v for i,v in enumerate(rad)]
                full=1
                for v in b:full*=roots[v]
                wc=[f.prod(v for v,exponent in zip(b,row) if exponent) for row in matrix]
                left=roots[wc[0]]*roots[wc[1]]*roots[wc[2]]
                right=roots[wc[3]]*roots[wc[4]]*roots[wc[5]]
                if full!=left*right:
                    raise RuntimeError('literal fibre product failed')
                checks+=1;counts[mask]+=full;Ucounts[mask]+=left;Dcounts[mask]+=right
        formula=q*q-6*q+(10 if p==3 else 12)
        if V!=formula or sum(counts)!=64*V:
            raise RuntimeError('base count or twist average failed')
        if sum(Ucounts)!=64*V or sum(Dcounts)!=64*V:
            raise RuntimeError('quotient twist marginal average failed')
        checks+=4
        receipts.append({'q':q,'characteristic':p,'quadratic_relation':ns,'base_points':V,
                         'sum_full_twist_counts':sum(counts),
                         'twist_count_histogram':dict(sorted(Counter(counts).items())),
                         'product_of_total_counts_hostiles':sum(Ucounts[i]*Dcounts[i]!=counts[i] for i in range(64))})
    print(json.dumps({'status':'FINITE-EXACT','checks':checks,'fields':receipts,
                      'universe':'Every affine branch-free point and all 64 constant square-class twists over F3,F5,F7,F11,F13,F17,F9,F25. Quotient roots counted directly from field squaring.',
                      'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                     indent=2,sort_keys=True))


if __name__=='__main__':
    main()
