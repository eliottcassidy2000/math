#!/usr/bin/env python3
"""Exact controls for Carmichael inheritance and the affine carry boundary.

Universe: every binary word of length <=16; all moduli 5..79 coprime to 6;
all integers 2..1200 for the Fermat-kernel index, every unit base included.
No floating-point claim is certified by this program.
"""
from collections import Counter, deque
from fractions import Fraction
from itertools import product
from math import gcd, prod, lcm
import json
from pathlib import Path

CHECKS = 0

def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)

def factors(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

def phi(n):
    return prod((p-1)*p**(e-1) for p,e in factors(n).items())

def mu(n):
    fs = factors(n)
    return 0 if any(e > 1 for e in fs.values()) else (-1)**len(fs)

def carry(word):
    B = 0
    for i, s in enumerate(word):
        B = 3**s * B + s * 2**i
    return B

def rotation_key(word):
    return min(word[i:] + word[:i] for i in range(len(word)))

def compose(x, y, m):
    a,b=x; c,d=y
    return (a*c % m, (a*d+b) % m)

def inverse(x,m):
    a,b=x; ai=pow(a,-1,m)
    return (ai, -ai*b % m)

def closure(gens,m):
    seen={(1,0)}; todo=deque(seen)
    while todo:
        x=todo.popleft()
        for g in gens:
            y=compose(g,x,m)
            if y not in seen:
                seen.add(y); todo.append(y)
    return seen

def kernel_control():
    examples=[]
    for n in range(2,1201):
        fs=factors(n)
        count=sum(pow(a,n-1,n)==1 for a in range(1,n) if gcd(a,n)==1)
        predicted=prod(gcd(n-1,p-1) for p in fs)
        check(count==predicted, ('kernel',n,count,predicted))
        index=prod(p**(e-1)*(p-1)//gcd(n-1,p-1) for p,e in fs.items())
        check(count*index==phi(n), ('index',n))
        korselt=all(e==1 and (n-1)%(p-1)==0 for p,e in fs.items())
        check((index==1)==korselt, ('Korselt',n))
        if n in (30,91,341,561,1105):
            examples.append({'N':n,'factorization':fs,'unit_count':phi(n),
                             'Fermat_liars':count,'kernel_index':index,
                             'base2_passes':gcd(2,n)==1 and pow(2,n-1,n)==1})
    return examples

def korselt(n):
    fs=factors(n)
    return len(fs)>1 and all(e==1 and (n-1)%(p-1)==0 for p,e in fs.items())

def successor_control():
    ps=[p for p in range(2,500) if factors(p)=={p:1}]
    examples=[]; cases=0
    for N in range(2,500):
        fs=factors(N)
        if any(e>1 for e in fs.values()): continue
        lam=lcm(*(p-1 for p in fs))
        children=[]
        for q in ps:
            if q>N or q in fs: continue
            criterion=(N-1)%(q-1)==0 and (N*q-1)%lam==0
            check(criterion==korselt(N*q),('arbitrary-core successor',N,q))
            if criterion: children.append(q)
            cases+=1
        if N in (21,33,91):
            examples.append({'core':N,'lambda':lam,'gcd_core_lambda':gcd(N,lam),
                             'all_one_prime_children':children})
    return {'comparisons':cases,'examples':examples}

def affine_control():
    rows=[]
    for m in range(5,80):
        if gcd(m,6)!=1: continue
        half=pow(2,-1,m)
        f0=(half,0); f1=(3*half % m,half)
        comm=compose(compose(compose(f0,f1,m),inverse(f0,m),m),inverse(f1,m),m)
        check(comm==(1,-pow(4,-1,m)%m), ('commutator',m,comm))
        group=closure((f0,f1),m)
        multipliers={a for a,b in group}
        expected={(a,b) for a in multipliers for b in range(m)}
        check(group==expected, ('group',m,len(group)))
        check(gcd(comm[1],m)==1,('derived primitive',m))
        rows.append({'modulus':m,'multiplier_order':len(multipliers),'group_order':len(group)})
    return rows

def word_control():
    zero={}; necklaces={}; primitive={}; examples={}
    word_total=0
    for q in range(1,17):
        hits=Counter(); orbits={}; prim={}
        for word in product((0,1),repeat=q):
            p=sum(word)
            if p==0 or p==q: continue
            word_total+=1
            D=2**q-3**p; m=abs(D); B=carry(word)
            pos=[i for i,s in enumerate(word) if s]
            direct=sum(2**i*3**(p-j-1) for j,i in enumerate(pos))
            check(B==direct, ('carry formulas',word))
            rotated=word[1:]+word[:1]
            check(2*carry(rotated)==3**word[0]*B+word[0]*D,
                  ('rotation',word))
            if m>1:
                half=pow(2,-1,m)
                aff=(1,0)
                for s in word:
                    aff=compose((3**s*half % m,s*half % m),aff,m)
                check(aff==(1,B*pow(2**q,-1,m)%m), ('clock carrier',word))
            if B%D==0:
                hits[p]+=1
                key=rotation_key(word)
                orbits.setdefault(p,set()).add(key)
                isprim=all(word!=word[:d]*(q//d) for d in range(1,q) if q%d==0)
                if isprim: prim.setdefault(p,set()).add(key)
                # Integer fixed point has exactly the prescribed parity word.
                n=B//D; x=n
                for s in word:
                    check(x%2==s, ('native parity',word,n,x))
                    x=(3*x+1)//2 if s else x//2
                check(x==n, ('native return',word,n,x))
            if (q,p) in ((3,1),(4,2),(5,3)):
                examples.setdefault(f'{q},{p}',[]).append({'word':''.join(map(str,word)),
                    'carry':B,'modulus':D,'zero':B%D==0})
        for p in range(1,q):
            zero[q,p]=hits[p]
            necklaces[q,p]=len(orbits.get(p,set()))
            primitive[q,p]=len(prim.get(p,set()))
            ds=[d for d in range(1,gcd(q,p)+1) if q%d==0 and p%d==0]
            burnside=sum(phi(d)*zero[q//d,p//d] for d in ds)
            mobius=sum(mu(d)*zero[q//d,p//d] for d in ds)
            check(burnside==q*necklaces[q,p],('Burnside',q,p))
            check(mobius==q*primitive[q,p],('Mobius',q,p))
            if p==1 and q>=3:
                check(zero[q,p]==0,('one-one hostile',q))
    rows=[{'q':q,'p':p,'clock':2**q-3**p,'zero_words':z,
           'necklaces':necklaces[q,p],'primitive_necklaces':primitive[q,p]}
          for (q,p),z in zero.items() if z]
    return {'words_checked':word_total,'nonempty_shapes':rows,'hostiles':examples}

def distribution_dp(q,p,m):
    banks=[Counter() for _ in range(p+1)]; banks[0][0]=1
    for j in range(q):
        nxt=[Counter() for _ in range(p+1)]
        for k in range(p+1):
            for b,n in banks[k].items():
                nxt[k][b]+=n
                if k<p: nxt[k+1][(3*b+2**j)%m]+=n
        banks=nxt
    return banks[p]

def distribution_control():
    tested=0
    for q in range(2,10):
        for p in range(1,q):
            for m in (5,7,11,13):
                direct=Counter(carry(w)%m for w in product((0,1),repeat=q) if sum(w)==p)
                check(distribution_dp(q,p,m)==direct,('DP',q,p,m))
                W=sum(direct.values()); r=W%m
                energy=m*sum(n*n for n in direct.values())-W*W
                check(energy>=r*(m-r),('quantized Fourier energy',q,p,m))
                if W<m:
                    tv=sum(abs(Fraction(direct.get(b,0),W)-Fraction(1,m)) for b in range(m))/2
                    check(tv>=1-Fraction(W,m),('sparse TV barrier',q,p,m))
                tested+=1
    return tested

def occupancy_control():
    # Independent universe, allowing arbitrary occupancies rather than carries.
    cases=0
    for m in range(2,7):
        for cs in product(range(4),repeat=m):
            W=sum(cs)
            if W==0: continue
            a,r=divmod(W,m)
            en=m*sum(c*c for c in cs)-W*W
            check(en>=r*(m-r),('arbitrary occupancy',cs))
            check((en==r*(m-r))==all(c in (a,a+1) for c in cs),('equality boundary',cs))
            cases+=1
    return cases

def main():
    report={'scope':'FINITE-EXACT controls; proofs are in INHERITANCE_AND_BRIDGES.md',
            'kernel_examples':kernel_control(),'affine_groups':affine_control(),
            'carry_necklaces':word_control(),'distribution_DP_instances':distribution_control(),
            'arbitrary_occupancy_instances':occupancy_control(),
            'arbitrary_core_successors':successor_control()}
    report['exact_checks']=CHECKS
    out=Path(__file__).with_suffix('.json')
    out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'exact_checks':CHECKS,'words_checked':report['carry_necklaces']['words_checked'],
                      'DP_instances':report['distribution_DP_instances'],'output':str(out)},sort_keys=True))

if __name__=='__main__': main()
