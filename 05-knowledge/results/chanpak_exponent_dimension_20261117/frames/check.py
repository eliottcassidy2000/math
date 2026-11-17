from fractions import Fraction as Q
from math import comb, factorial
import itertools

def mean(n):
    h=Q(factorial(n),2**(n-1))
    path=closed=Q(1)
    for k in range(2,n,2):
        run=sum(2**c*comb(k//2-1,c-1)*comb(n-k,c) for c in range(1,k//2+1))
        path+=Q(run*factorial(n-k),factorial(n))
        closed+=Q(n*run*factorial(n-k-1),factorial(n))
    return tuple(h*a for a in (path, closed, 2*path-closed))

def literal(n):
    pairs=list((x,y) for x in range(n) for y in range(x+2,n))
    pseq=[]
    index={p:i for i,p in enumerate(pairs)}
    for w in itertools.permutations(range(n)):
        pos=neg=0
        for a,b in zip(w,w[1:]):
            if b==a+1:continue
            if b==a-1:break
            j=index[tuple(sorted((a,b)))]
            if a<b:pos|=1<<j
            else:neg|=1<<j
        else:pseq.append((pos,neg))
    i=index[(0,n-1)];sc=so=numf=ties=nump=0
    for bn in range(1<<len(pairs)):
        if bn>>i&1:continue
        twin=((1<<len(pairs))-1)^bn
        hc=sum(bn&p==p and bn&q==0 for p,q in pseq)
        ho=sum(twin&p==p and twin&q==0 for p,q in pseq)
        sc+=hc; so+=ho
        if hc>ho:nump+=1
        elif hc==ho:ties+=1
        else:numf+=1
    den=1<<(len(pairs)-1)
    return (Q(sc+so,den*2),Q(sc,den),Q(so,den)),(nump,ties,numf)

for n in range(3,8):
    a,signs=literal(n)
    if a!=mean(n):raise ArithmeticError((n,a,mean(n)))
    print(n,tuple(map(str,a)),signs)
