"""Exact finite controls for Pollock's five-contact reductions."""
from itertools import combinations, permutations, combinations_with_replacement, product
from math import isqrt, prod

def F(x):
    return (x**3-x)//6

def primes(n):
    return [p for p in range(2,n+1) if all(p%j for j in range(2,isqrt(p)+1))]

def roots_prime(p):
    # Keep one regular preimage for every finite-field value.
    roots={}
    for x in range(p):
        if (3*x*x-1)%p:
            roots[(x**3-x)%p]=x
    return roots

def pair_pp(r,p,e):
    if p in (2,3):
        x,y,q,n=0,1,p,e+1
    else:
        roots=roots_prime(p)
        v=next(v for v in roots if (6*r-v)%p in roots)
        x,y,q,n=roots[v],roots[(6*r-v)%p],p,e
    target=6*r-(y**3-y)
    for _ in range(1,n):
        assert (target-x**3+x)%q==0
        digit=((target-x**3+x)//q)*pow(3*x*x-1,-1,p)%p
        x+=q*digit
        q*=p
    assert (F(x)+F(y)-r)%(p**e)==0
    return x,y

for p in primes(300):
    if p>=5:
        regular=roots_prime(p)
        assert set(regular)=={(x**3-x)%p for x in range(p)}
        chi=1 if pow((-3)%p,(p-1)//2,p)==1 else -1
        assert 3*len(regular)==2*p+chi

count=0
for p in primes(61):
    for e in range(1,5):
        q=p**e
        for r in set(range(min(q,200)))|{q-1,q//2,q//3}:
            pair_pp(r,p,e)
            count+=1
print('finite-field regular image and prime-power lift checks:',count)

def crt(pairs,mods):
    M=prod(mods)
    return tuple(sum((M//q)*pow(M//q,-1,q)*pair[i] for pair,q in zip(pairs,mods))%M for i in range(2))

factors=[(2,5),(3,4),(5,3),(11,2)]
mods=[p**(e+int(p in (2,3))) for p,e in factors]
M=prod(p**e for p,e in factors)
L=prod(mods)
for k,r in product(range(1,80),range(0,123,7)):
    g=19
    w=[a+j for j in range(5) for a in [(-(k+1))%L+L]]
    x,y=crt([pair_pp(r+F(k+1)-F(g),p,e) for p,e in factors],mods)
    for j,z in enumerate(w):
        assert (F(x+L)+F(y+L)+F(g)+F(z)-r-F(k+1)+F(k+1-j))%M==0
print('five-anchor joint residue phases:',79*len(range(0,123,7)), 'modulus:',M)

# The two real eigenspaces for each inverse pair of pentagonal orderings
# are distinct, as their orthogonal projectors recover their 5-cycle.
cycles=set()
for p in permutations(range(1,5)):
    order=(0,)+p
    cycles.add(frozenset(tuple(sorted((order[i],order[(i+1)%5]))) for i in range(5)))
assert len(cycles)==12
assert all(frozenset(set(combinations(range(5),2))-c) in cycles for c in cycles)
pairs=list(combinations(range(5),2))
edges=[(a,b) for a,b in combinations(pairs,2) if not set(a)&set(b)]
assert len(edges)==15
assert all(sum(a in e for e in edges)==3 for a in pairs)
print('Petersen rational lines:',len(edges),'golden projector cycles:',len(cycles))

def norm(t):
    y=((t-1)*(t+2),(t-2)*(t+1),2,-t*t,2-t*t)
    q=2*t**4-5*t*t+8
    assert sum(y)==sum(v**3 for v in y)==0 and sum(v*v for v in y)==2*q
    return y,q

for t in range(4,2000,121):
    y,q=norm(t)
    assert q%121==77
    c=1-min(y)
    assert sum(F(c+v) for v in y)==5*F(c)+c*q
print('first nonlinear norm witness:',norm(4))

balanced_primes={}
for c in range(1,7):
    for xs in combinations_with_replacement(range(1,5*c-3),5):
        if sum(xs)!=5*c or sum((z-c)**3 for z in xs):
            continue
        N=sum(map(F,xs))
        if N>=2 and all(N%p for p in range(2,isqrt(N)+1)):
            balanced_primes[N]=xs
assert max(balanced_primes)==421
print('largest cubic-balanced prime:',max(balanced_primes),balanced_primes[421])

for xs in combinations_with_replacement(range(1,19),5):
    S=sum(xs)
    N=sum(map(F,xs))
    D=sum((5*x-S)**3 for x in xs)//5
    assert 150*N-D == S*(15*sum(x*x for x in xs)-2*S*S-25)
    assert D%5==2*S**3%5
print('integer defect factorization checked on',len(list(combinations_with_replacement(range(1,19),5))),'tuples')
