"""Reproduce the exact row-intersection and residue-variance accounting."""
from collections import Counter, defaultdict
from fractions import Fraction
from math import log

def need(ok, label):
    if not ok:
        raise RuntimeError(label)

def checksum(n):
    out=0
    while n:
        x,b=divmod(n,2)
        if not b:
            out+=x+(x&1)
        n=x
    return out

def run(T):
    D=1
    while D**5<T:
        D+=1                         # delta=1/5, exact integer root
    top=T//(100*D)
    ps=list(map(checksum,range(top+1)))
    counts=Counter()
    rows=[]
    for largest in range(D,top+1):
        for orient in [False,True]:
            row=set()
            previous=-1
            for shorter in range(D,min(largest-int(orient),T//(100*largest))+1):
                x,y=(shorter,largest) if orient else (largest,shorter)
                n=100*x*y+89*x+11*(ps[x]+ps[y])+11*(x&1)+11
                need(n>previous,(T,x,y,n,previous))
                need(n%2==1,(T,x,y,n))
                previous=n
                if T//2<n<=T:
                    row.add(n)
            if row:
                counts.update(row)
                rows.append(row)
    collisions=sum(r*(r-1) for r in counts.values())
    # Independent exact row-intersection algorithm, streamed by row.
    earlier=Counter()
    cross=0
    for row in rows:
        cross+=2*sum(earlier[n] for n in row)
        earlier.update(row)
    need(cross==collisions,'row intersection')
    means=defaultdict(Counter)
    for n in range(T//2+1,T+1,2):
        means[n%22]['total']+=counts[n]
        means[n%22]['size']+=1
    binmeans={b:Fraction(c['total'],c['size']) for b,c in means.items()}
    variance=sum((counts[n]-binmeans[n%22])**2 for n in range(T//2+1,T+1,2))
    parity_baseline=sum(Fraction(c['total']**2,c['size']) for c in means.values())
    need(variance==collisions+sum(counts.values())-parity_baseline,'variance')
    l=0.6*log(T)/100
    print(f'T={T}: rows={len(rows)}, mass={sum(counts.values())}, support={len(counts)}, '
          f'off-diagonal={collisions}, off-diagonal/(T ell^2)={collisions/(T*l*l):.4f}')

if __name__=='__main__':
    c=[Fraction(21,11)]+[Fraction(10,11)]*10
    need(sum(c)/11==1,'mean')
    need(sum(x*x for x in c)/11==Fraction(131,121),'constant')
    for k in (13, 15, 17, 19):
        run(2**k)
