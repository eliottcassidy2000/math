from itertools import product, permutations
from fractions import Fraction
from math import gcd, factorial

def check(b):
    if not b: raise AssertionError('scaffold Burnside check')

def part(n, lo=1):
    if not n:
        yield ()
    else:
        for v in range(lo,n+1,2):
            for p in part(n-v,v):yield (v,)+p

def z(p):
    out=1
    for v in set(p):
        m=p.count(v);out*=v**m*factorial(m)
    return out

def phi(b):
    return sum(gcd(b,x)==1 for x in range(1,b+1))

def exp(p):
    return sum((v-1)//2 for v in p)+sum(gcd(x,y) for i,x in enumerate(p) for y in p[:i])

def by_formula(r, kinds):
    out=Fraction(0)
    for p0 in part(r):
        for aa in product(*(tuple(part(m)) for l,c,m in kinds)):
            den=z(p0)
            groups=[]
            for (l,c,_), p in zip(kinds,aa):
                den*=z(p)
                groups.extend((a,l) for a in p)
            for bb in product(*(tuple(b for b in range(1,l+1) if l%b==0) for a,l in groups)):
                wt=Fraction(1,den)
                pp=[]
                for (a,l),b in zip(groups,bb):
                    wt*=Fraction(phi(b),l);pp.extend([a*b]*(l//b))
                out+=wt*2**(exp(p0+tuple(pp))-len(pp))
    check(out.denominator==1)
    return out.numerator

def literal(r,kinds):
    x=r; scaffold=[]
    for l,c,m in kinds:
        for _ in range(m):
            scaffold.extend((x+i,x+(i+1)%l,c) for i in range(l))
            x+=l
    n=x
    marks=set(scaffold)
    group=[p for p in permutations(range(n))
          if {(p[i],p[j],c) for i,j,c in scaffold}==marks]
    lock={(min(i,j),max(i,j)) for i,j,c in marks}
    pairs=[(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in lock]
    pairid={p:i for i,p in enumerate(pairs)}
    actions=[[(pairid[min(p[i],p[j]),max(p[i],p[j])],p[i]>p[j]) for i,j in pairs] for p in group]
    seen=set(); count=0
    for m in range(1<<len(pairs)):
        if m in seen:continue
        count+=1
        for op in actions:
            y=sum((((m>>i)&1)^k)<<j for i,(j,k) in enumerate(op))
            seen.add(y)
    check(len(seen)==1<<len(pairs))
    return count

fixtures=[(n,()) for n in range(6)]
fixtures.extend((r,((3,c,1),)) for r in range(4) for c in (0,1))
fixtures.extend((r,((5,0,1),)) for r in range(2))
fixtures.extend([(0,((3,0,2),)),(0,((3,0,1),(3,1,1))),
                 (0,((3,1,2),))])
for r,kinds in fixtures:
    check(by_formula(r,kinds)==literal(r,kinds))
check([by_formula(n,()) for n in range(6)]==[1,1,1,2,4,12])
check(sum(by_formula(*x) for x in [(3,()),(0,((3,0,1),)),(0,((3,1,1),))])==4)
# At order three: two path-marked classes; four packing-marked
# classes (one transitive; three for the directed triangle).
print(len(fixtures),'literal scaffold-orbit cases, including indistinguishable and distinct colors.')
