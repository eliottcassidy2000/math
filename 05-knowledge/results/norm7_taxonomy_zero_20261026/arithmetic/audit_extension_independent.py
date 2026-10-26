#!/usr/bin/env python3
"""Independent literal edge-product audit, using different extension moduli."""
from itertools import combinations, product
from pathlib import Path
import json

CHECKS=0
def demand(test,label):
    global CHECKS
    CHECKS+=1
    if not test:
        raise RuntimeError(label)


class FF:
    def __init__(self,p,modulus):
        self.p=p; self.f=modulus; self.n=len(modulus)-1
        self.elems=list(product(range(p),repeat=self.n))
        self.zero=(0,)*self.n; self.one=(1,)+(0,)*(self.n-1)
        sq={self.mul(x,x) for x in self.elems if x!=self.zero}
        self.ch={x:(0 if x==self.zero else 1 if x in sq else -1) for x in self.elems}
        self.inv={x:next(y for y in self.elems if self.mul(x,y)==self.one)
                  for x in self.elems if x!=self.zero}
        self.diff={(x,y):self.ch[self.sub(x,y)] for x in self.elems for y in self.elems}
        self.epsilon=self.diff[self.zero,self.one]
        self.L={l:sum(self.ch[x]*self.diff[x,self.one]*self.diff[x,l] for x in self.elems)
                for l in self.elems}
        demand(len(sq)*2==len(self.elems)-1,'independent field square count')

    def sub(self,x,y):
        return tuple((a-b)%self.p for a,b in zip(x,y))

    def mul(self,x,y):
        c=[0]*(2*self.n-1)
        for i,a in enumerate(x):
            for j,b in enumerate(y):
                c[i+j]=(c[i+j]+a*b)%self.p
        for deg in range(len(c)-1,self.n-1,-1):
            z=c[deg]
            for j,a in enumerate(self.f[:-1]):
                c[deg-self.n+j]=(c[deg-self.n+j]-z*a)%self.p
        return tuple(c[:self.n])


def literal(F,edges):
    adj={i:{b if a==i else a for a,b in edges if i in (a,b)} for i in range(7)}
    cones=next(C for C in combinations(range(7),3)
               if all(len(adj[i])==3 for i in C)
               and all(j not in adj[i] for i,j in combinations(C,2)))
    bases=[i for i in range(7) if i not in cones]
    baseedges=[(a,b) for a,b in edges if a in bases and b in bases]
    domain=[x for x in F.elems if x not in (F.zero,F.one)]
    total=0
    for s,t in product(domain,repeat=2):
        if s==t:
            continue
        X=dict(zip(bases,(F.zero,F.one,s,t)))
        value=1
        for a,b in baseedges:
            value*=F.diff[X[b],X[a]]
        for v in cones:
            term=0
            for x in F.elems:
                edgevalue=1
                for w in adj[v]:
                    edgevalue*=F.diff[X[w],x] if v<w else F.diff[x,X[w]]
                term+=edgevalue
            value*=term
        total+=value
    return total


def main():
    root=Path(__file__).resolve().parents[4]
    graphs=json.loads((root/'04-computation/experiments/pythagorean_hourglass_20261025/tt7/motives/genuine37.json').read_text())
    ids=[2,4,5,6,9,18]
    specs=[(5,[0,1]),(7,[0,1]),(3,[2,1,1]),(3,[1,2,0,1]),(7,[3,1,1])]
    rows=[]
    for p,modulus in specs:
        F=FF(p,modulus)
        q=len(F.elems)
        vals=[literal(F,graphs[i]['edges']) for i in ids]
        domain=[x for x in F.elems if x not in (F.zero,F.one)]
        A={u:sum(F.ch[s]*F.diff[u,s]*F.L[s]**2 for s in domain) for u in domain}
        B={u:sum(F.ch[s]*F.diff[s,F.one]*F.diff[u,s]*F.L[s]**2 for s in domain) for u in domain}
        for u in domain:
            T=F.mul(u,F.inv[F.sub(u,F.one)])
            factor=F.epsilon*F.diff[u,F.one]
            demand(F.L[T]==factor*F.L[u],'Legendre involution sign')
            demand(A[T]==factor*A[u],'A involution sign')
            demand(B[T]==factor*B[u],'B involution sign')
            demand(F.L[F.inv[u]]==F.ch[u]*F.L[u],'Legendre inverse sign')
        S=lambda weights:sum(weights(u)*F.ch[u]*F.L[u] for u in domain)
        transformed=[F.epsilon*S(lambda u:A[u]),S(lambda u:F.diff[u,F.one]*A[u]),
                     S(lambda u:B[u]),S(lambda u:F.diff[u,F.one]*B[u])]
        demand(vals[:4]==transformed,'literal graph versus transformed pair formulas')
        demand(vals[1]==F.epsilon*vals[0] and vals[2]==vals[3],'paired identities')
        rows.append({'p':p,'q':q,'modulus':modulus,'epsilon':F.epsilon,'cone_ids':ids,'literal_sums':vals})
    print(json.dumps({'status':'PASS','checks':CHECKS,'method':'literal cubic products; distinct extension moduli; independent A/B transforms','rows':rows},indent=2,sort_keys=True))


if __name__=='__main__':
    main()
