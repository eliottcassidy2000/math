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

def verify_carry(r):
    V=1<<r
    def rev(b):
        return int(f'{b:0{r}b}'[::-1],2)
    vals=set()
    for u in range(V*V,2*V*V):
        for s in range(V//2,V):
            recorded={}
            for b in range(V//8):
                c=s-b
                carry=K=0
                for j in range(r):
                    t=((b>>j)&1)+((c>>j)&1)+carry
                    need(t%2==(s>>j)&1,('carry',r,b,c,j))
                    carry=t//2
                    if j<r-1:
                        K=2*K+carry
                need(carry==0 and rev(b)+rev(c)==rev(s)+3*K,('reverse',r,b,c))
                x,y=V*u+b,V*u+c
                n=100*x*y+89*x+11*(checksum(x)+checksum(y))+11+11*(x%2)
                need(n%2==1 and n<=500*V**6 and n not in vals,('image',r,x,y,n))
                if K in recorded:
                    need(n>recorded[K],('within carry',r,x,y,n))
                recorded[K]=n
                vals.add(n)
    need(len(vals)==V**4//16,('count',r))
    print(f'carry-separated V={V}: {len(vals)} distinct counts, maximum {max(vals)}')

if __name__=='__main__':
    for r in (3,4,5):
        verify_carry(r)

def verify_averaging(r):
    V=1<<r
    def rev(b):
        return int(f'{b:0{r}b}'[::-1],2)
    def K(b,s):
        return (rev(b)+rev(s-b)-rev(s))//3
    def row(u,s,b):
        x,y=V*u+b,V*u+s-b
        return 100*x*y+89*x+11*(checksum(x)+checksum(y))+11+11*(x%2)
    counts=Counter()
    exact_ordered_collisions=0
    for s in range(V//2,V):
        for u in range(V,2*V):
            for b in range(V//8):
                n=row(u,s,b)
                need(n%2==1 and n<=500*V**4, ('bounds averaging',r,n))
                counts[n]+=1
        for b in range(V//8):
            for bb in range(b):
                slope=-33*(K(b,s)-K(bb,s))
                for parity in (0,1):
                    base=V+parity
                    value=row(base,s,b)-row(base,s,bb)
                    need(row(base+2,s,b)-row(base+2,s,bb)-value==2*slope,
                         ('incidence slope',r,s,b,bb,parity))
                    if not slope:
                        need(value!=0, ('same carry',r,s,b,bb,parity))
                        continue
                    quotient,remainder=divmod(-value,2*slope)
                    if remainder==0 and V<=base+2*quotient<2*V:
                        u=base+2*quotient
                        need(row(u,s,b)==row(u,s,bb), ('incidence inverse',r,s,b,bb,u))
                        exact_ordered_collisions+=2
    P=V**3//16
    C=sum(k*(k-1) for k in counts.values())
    need(sum(counts.values())==P and C==exact_ordered_collisions,
         ('all incidence collisions',r,C,exact_ordered_collisions))
    need(C<=V**3//64 and 128*len(counts)>=7*V**3, ('support floor',r,C))
    print(f'averaged V={V}: {P} presentations, {len(counts)} distinct, {C} ordered collisions')

if __name__=='__main__':
    for r in (4,5,6,7):
        verify_averaging(r)
