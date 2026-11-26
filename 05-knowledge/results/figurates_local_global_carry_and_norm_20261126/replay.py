from itertools import product, combinations_with_replacement
from math import isqrt

def F(x): return (x**3-x)//6
def P(m,x): return ((m-2)*x*x-(m-4)*x)//2
def tetra(xs): return sum(map(F,xs))
for m,x in product(range(3,24),range(40)):
    assert sum(P(m,i) for i in range(1,x+1))==(m-2)*F(x+1)+(3-m)*x*(x+1)//2
print('polygonal layer identity checked on 840 cases')

for t,a,b,d in product(range(1,9),range(1,9),range(1,9),range(1,9)):
    even=(t,t+a+b,t+a+d,t+b+d)
    odd=(t+a,t+b,t+d,t+a+b+d)
    assert sum(even)==sum(odd)
    assert sum(x*x for x in even)==sum(x*x for x in odd)
    assert tetra(odd)-tetra(even)==a*b*d
print('4096 mixed carries preserve the two lower moments')

residues={}
for xs in combinations_with_replacement(range(4),3):
    key=(sum(xs)%4,sum(x*x for x in xs)%4,sum(xs[i]*xs[j] for i,j in ((0,1),(0,2),(1,2)))%4)
    r=xs[0]*xs[1]*xs[2]%4
    assert residues.setdefault(key,r)==r
print('translation-stable odd-gap obstruction checked exhaustively mod 4')

for k in range(1,9):
    c=2**k
    assert len({F(2*u)%c for u in range(1,c+1)})==c
    qs={a*a+b*b for a in range(c) for b in range(c) if a*a+b*b<=(c-1)**2}
    totals={4*F(c)+F(2*u)+c*q for u in range(1,c+1) for q in qs}
    assert len(totals)==c*len(qs)
    print('dyadic block',c,'norms',len(qs),'targets',len(totals))

for k in range(4):
    m=2*5**k+1; q=3*25**k; n=2*F(m)+2*F(23*m)+m*q
    assert q <= (m-1)**2
    for b in range(isqrt(q//23)+1):
        a2=q-23*b*b
        assert isqrt(a2)**2!=a2
    print('locally soluble unequal-pair obstruction',k,m,n)
N=2*F(3)+2*F(69)+9
atoms={F(x):x for x in range(1,102)}
for ix in range(1,75):
    if (rem:=N-F(ix)) in atoms:
        print('first obstruction target, unrestricted pair witness',ix,atoms[rem]);break
else:
    for a,b,c in combinations_with_replacement(range(1,102),3):
        rem=N-F(a)-F(b)-F(c)
        if rem in atoms:
            print('first obstruction target, unrestricted four-term witness',a,b,c,atoms[rem]);break
