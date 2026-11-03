#!/usr/bin/env python3
"""Exact controls for the six-state twist and tribonacci marked CA. Stdlib only."""
from collections import Counter, deque
from itertools import product
import json
import random


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


STATES = ((0,0,0,0),(0,0,0,1),(0,0,1,0),(0,1,1,1),(1,0,1,1),(1,1,1,1))
LETTERS = ((0,0),(0,1),(1,0),(1,1))


def transition(q, a, b):
    w, v, A, B = q
    return ((v if b else w) ^ (a & b) ^ (a & B),
            (w if a else v) ^ (a & b) ^ (b & A), A | a, B | b)


def twist(level, a, b):
    q = STATES[0]
    for j in range(level):
        q = transition(q, (a >> j) & 1, (b >> j) & 1)
    return q[0]


def recursive_product(level, a, b):
    """Independent high-digit implementation of the pair multiplication rule."""
    if level == 0:
        return 1, 0
    half = 1 << (level - 1)
    A, B, ah, bh = a % half, b % half, a // half, b // half
    if not ah and not bh:
        return recursive_product(level - 1, A, B)
    if not ah:
        sign, index = recursive_product(level - 1, B, A)
        return sign, half + index
    if not bh:
        sign, index = recursive_product(level - 1, A, B)
        return sign * (-1 if B else 1), half + index
    sign, index = recursive_product(level - 1, B, A)
    return -sign * (-1 if B else 1), index


def perm(level, a, sign=0):
    size = 1 << level
    return tuple((((u // size) ^ sign ^ twist(level, a, u % size)) * size
                  + (a ^ (u % size))) for u in range(2 * size))


def compose(p, q):
    return tuple(p[q[x]] for x in range(len(p)))


def power(p, n):
    ans = tuple(range(len(p)))
    while n:
        if n & 1:
            ans = compose(p, ans)
        p = compose(p, p)
        n //= 2
    return ans


def cycles(p):
    unseen = set(range(len(p)))
    lengths = Counter()
    while unseen:
        x = min(unseen)
        y = p[x]
        length = 1
        unseen.remove(x)
        while y != x:
            unseen.remove(y)
            y = p[y]
            length += 1
        lengths[length] += 1
    return dict(sorted(lengths.items()))


def cd_checks():
    reachable = {STATES[0]}
    while True:
        enlarged = reachable | {transition(q,a,b) for q in reachable for a,b in LETTERS}
        if enlarged == reachable:
            break
        reachable = enlarged
    need(reachable == set(STATES), 'six reachable twist states')
    distinguishing = []
    for i in range(6):
        for j in range(i+1,6):
            found = None
            for length in range(3):
                for word in product(LETTERS, repeat=length):
                    p,q = STATES[i],STATES[j]
                    for a,b in word:
                        p,q = transition(p,a,b),transition(q,a,b)
                    if p[0] != q[0]:
                        found = word
                        break
                if found is not None:
                    break
            need(found is not None, 'pairwise distinguishability')
            distinguishing.append({'states':[i,j], 'suffix':found})
    pair_counts = {}
    for n in range(9):
        for a,b in product(range(1<<n),repeat=2):
            s,c = recursive_product(n,a,b)
            need(((-1)**twist(n,a,b),a^b)==(s,c), 'independent CD multiplication')
        pair_counts[n] = 4**n
    rng = random.Random(168)
    transport_controls=[]
    for n in range(1,9):
        N=1<<n
        identity=tuple(range(2*N))
        generators=[perm(n,a) for a in range(N)]
        for _ in range(200):
            p=identity
            for _ in range(rng.randrange(1,30)):
                p=compose(generators[rng.randrange(N)],p)
            a=p[0]%N
            f=[p[x]//N for x in range(N)]
            need(all(p[s*N+x]==((s^f[x])*N+(x^a)) for s in (0,1) for x in range(N)), 'signed-XOR form')
            need(power(p,4)==identity, 'four-lap transport ceiling')
            if a==0:
                predicted={1:2*f.count(0),2:f.count(1)}
            else:
                h=[f[x]^f[x^a] for x in range(N)]
                predicted={2:h.count(0),4:h.count(1)//2}
            predicted={k:v for k,v in predicted.items() if v}
            need(cycles(p)==predicted,'cycle multiplicity formula')
        need(len({power(generators[1],r)[0] for r in range(4)})==4,'sharp order-four basis orbit')
        transport_controls.append({'level':n,'random_words':200,'degree':2*N})
    L=[perm(3,a) for a in range(8)]
    defect=compose(power(L[3],3),compose(L[1],L[2]))
    need(cycles(defect)=={1:8,2:4},'octonion hostile control')
    return {'transition_table':[[STATES.index(transition(q,a,b)) for a,b in LETTERS] for q in STATES],
            'distinguishing_suffixes':distinguishing,'product_pairs_by_level':pair_counts,
            'transport_controls':transport_controls,'octonion_defect_cycles':cycles(defect)}


# Elements of F_27 are a+3b+9c, representing a+b alpha+c alpha^2.
def coeff(a):
    return (a%3,(a//3)%3,a//9)


def enc(v):
    return sum((x%3)*3**j for j,x in enumerate(v))


def add(a,b):
    return enc([x+y for x,y in zip(coeff(a),coeff(b))])


def mul(a,b):
    out=[0]*5
    for i,x in enumerate(coeff(a)):
        for j,y in enumerate(coeff(b)):
            out[i+j]+=x*y
    for d in (4,3):
        x=out[d]
        for e in (d-1,d-2,d-3):
            out[e]+=x
    return enc(out[:3])


def fpow(a,n):
    ans=1
    while n:
        if n&1:ans=mul(ans,a)
        a=mul(a,a)
        n//=2
    return ans


ALPHA=3
BETA=add(ALPHA,1)
ALPHAMUL=[mul(ALPHA,x) for x in range(27)]


def clock_step(x):
    k=len(x)
    return tuple(add(ALPHAMUL[x[i]],x[(i+1)%k]) for i in range(k))


def guard_to_field(v):
    s,u,d=v
    return enc((d,u-s,s))


def field_to_guard(v):
    d,b,s=coeff(v)
    return s,(b+s)%3,d


def guard_step(v):
    s,u,d=v
    return u,(s+u+d)%3,s


def guard_ca_step(x):
    k=len(x)
    return tuple(tuple((a+b)%3 for a,b in zip(guard_step(x[i]),x[(i+1)%k])) for i in range(k))


def valuation(x):
    """Order at z=1 for X(z)=sum x_i z^i, independent finite differences."""
    from math import comb
    for d in range(len(x)):
        c=0
        for i in range(d,len(x)):
            c=add(c,mul(comb(i,d)%3,x[i]))
        if c:return d
    return len(x)


def scalar_multiply(x,a):
    return tuple(mul(a,v) for v in x)


def clock_checks():
    need(all((x**3-x*x-x-1)%3 for x in range(3)),'irreducible cubic')
    need(fpow(ALPHA,13)==1 and ALPHA!=1,'alpha order thirteen')
    need(fpow(BETA,13)==2 and fpow(BETA,2)!=1,'beta order twenty-six')
    need(all(mul(x,y)!=0 for x in range(1,27) for y in range(1,27)),'field no zero divisors')
    need(all(guard_to_field(guard_step(v))==mul(ALPHA,guard_to_field(v)) for v in product(range(3),repeat=3)), 'literal guard intertwining')
    controls=[]
    for k in (1,3,9,27,81):
        m=26*k
        seed=(1,)+(0,)*(k-1)
        seen=set()
        x=seed
        histogram=Counter()
        for t in range(m):
            need(x not in seen,'delta orbit no premature return')
            seen.add(x)
            histogram[x[0]]+=1
            x=clock_step(x)
        need(x==seed,'delta orbit return')
        need(histogram==Counter({v:k for v in range(1,27)}),'marked-site exact uniform occupation')
        need(histogram[1]%2==1,'negation lift has exactly two laps')
        # Separate matrix-coordinate implementation on the first hundred steps.
        x=seed
        y=tuple(field_to_guard(v) for v in x)
        for t in range(min(m,100)):
            x=clock_step(x)
            y=guard_ca_step(y)
            need(x==tuple(guard_to_field(v) for v in y),'independent guard-register CA')
        # Lift of order four: count events, using explicit level-five basis permutation.
        L1=perm(5,1)
        x,u=seed,0
        for t in range(4*m):
            if t:
                need((x,u)!=(seed,0),'lift no premature return')
            if x[0]==1:u=L1[u]
            x=clock_step(x)
        need((x,u)==(seed,0),'four-lap return')
        controls.append({'sites':k,'clock_period':m,'negation_lift_period':2*m,'basis_lift_period':4*m,'visits_per_nonzero_local_symbol':k})
    # Full dynamical classification at k=3: all 27^3 configurations, no guessed periods.
    k=3
    unseen=set(product(range(27),repeat=k))
    orbit_counts=Counter()
    state_counts=Counter()
    while unseen:
        x=min(unseen)
        root=x
        orbit=[]
        while True:
            need(x in unseen,'disjoint permutation cycles')
            unseen.remove(x)
            orbit.append(x)
            x=clock_step(x)
            if x==root:break
        period=len(orbit)
        orbit_counts[period]+=1
        state_counts[period]+=period
        for z in orbit:
            v=valuation(z)
            if v==k:
                expected=1
            else:
                q=1
                while q<k-v:q*=3
                expected=26*q
            need(period==expected,'all-state valuation-period formula')
    need(orbit_counts==Counter({1:1,26:1,78:252}),'three-site full cycle inventory')
    # Observation ℓ=(1,1,0): determinant two, and reconstruction from time snapshots.
    rng=random.Random(169)
    for k in (1,3,9):
        for _ in range(30):
            x=tuple(tuple(rng.randrange(3) for _ in range(3)) for _ in range(k))
            x1=guard_ca_step(x);x2=guard_ca_step(x1)
            Y0=tuple((s+u)%3 for s,u,d in x)
            Y1=tuple((s+u)%3 for s,u,d in x1)
            Y2=tuple((s+u)%3 for s,u,d in x2)
            for i in range(k):
                u0=Y0[i]
                u1=(Y1[i]-Y0[(i+1)%k])%3
                u2=(Y2[i]-2*Y1[(i+1)%k]+Y0[(i+2)%k])%3
                s,u,d=x[i]
                need((u0,u1,u2)==((s+u)%3,(s+2*u+d)%3,(3*s+3*u+2*d)%3),'three-snapshot observability')
    return {'alpha_powers':{'13':fpow(ALPHA,13)},'beta_powers':{'2':fpow(BETA,2),'13':fpow(BETA,13)},
            'controls':controls,'full_k3_orbit_counts':dict(sorted(orbit_counts.items())),
            'full_k3_state_counts':dict(sorted(state_counts.items()))}


def main():
    result={'cayley_dickson':cd_checks(),'tribonacci_ca':clock_checks(),
            'scope':'All-level statements have proofs in NOTE.md. Finite controls do not replace proofs. The tournament encoding inherits the audited S167 compiler.'}
    print(json.dumps(result,sort_keys=True,indent=2))


if __name__=='__main__':main()
