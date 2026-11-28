"""Replay branch rings on any rooted Moore chart with the repo adjacency layout."""
from collections import Counter
from itertools import product
from math import comb, prod
from fractions import Fraction
from pathlib import Path
import random
import sys

def compose(p, q):
    return tuple(q[y] for y in p)

def fixed(p):
    return sum(x==y for x,y in enumerate(p))

def transfer(q, k, f):
    return (q-1)**k + (-1)**k*(f-1)

def ring_polynomial(m, q, f):
    return [comb(m,k)*transfer(q,k,f) for k in range(m+1)]

def controls():
    # Direct enumeration of departure sheets and parent traversals, including
    # identity, derangement, and mixed fixed-point monodromy.
    for q in range(2,5):
        e=tuple(range(q))
        for m in range(3,7):
            for first in [e,e[1:]+e[:1],(1,0)+e[2:]]:
                p=[first]+[e]*(m-1)
                counts=Counter()
                for word in product(e,repeat=m):
                    k=sum(p[i-1][word[i-1]]!=word[i] for i in range(m))
                    counts[k]+=1
                assert [counts[k] for k in range(m+1)]==ring_polynomial(m,q,fixed(first))
                assert sum(counts.values())==q**m
    for m in range(3,15):
        sets=Counter(mask.bit_count() for mask in range(1<<m)
                     if not (mask&(mask<<1)) and not (mask&1 and mask>>(m-1)))
        assert 1+sets[2]==comb(m-1,2)
        assert all(sets[k]==(1 if not k else m*comb(m-k-1,k-1)//k)
                   for k in range(m//2+1))
    print('Direct transfer and schedule controls pass')

def load(path):
    rows = [list(map(int,l.split())) for l in Path(path).read_text().splitlines()]
    n, = rows.pop(0)
    A = [set(row[1:]) for row in rows]
    assert len(A)==n and all(len(a)==r[0] for a,r in zip(A,rows))
    d=len(A[0]); q=d-1; off=d+1
    assert n==d*d+1
    P={}
    for a in range(d):
        for b in range(d):
            if a==b: continue
            opts=[tuple(y-off-b*q for y in A[off+a*q+x] if off+b*q<=y<off+(b+1)*q)
                  for x in range(q)]
            assert max(map(len,opts))<=1
            if all(opts):
                p=tuple(y for (y,) in opts)
                assert sorted(p)==list(range(q))
                P[a,b]=p
    return d,q,P

def survey(path):
    d,q,P=load(path)
    adj=[sorted(b for b in range(d) if (a,b) in P) for a in range(d)]
    bank_stats={20:Counter(),21:Counter()}
    maxbank=None
    rng=random.Random(20261128)
    for m in (20,21):
        orders=set()
        while len(orders)<10000:
            order=[rng.randrange(d)]
            while len(order)<m:
                opts=[b for b in adj[order[-1]] if b not in order
                      and (len(order)<m-1 or (b,order[0]) in P)]
                if not opts: break
                order.append(rng.choice(opts))
            if len(order)!=m: continue
            j=order.index(min(order)); o=tuple(order[j:]+order[:j])
            o=min(o,(o[0],)+o[:0:-1])
            if o in orders: continue
            orders.add(o)
            a,c=o[-1],o[0]
            R=tuple(range(q))
            for x,y in zip(o,o[1:]):
                R=compose(R,P[x,y])
            f=fixed(compose(R,P[a,c]))
            bank_stats[m][f]+=1
            if m!=20: continue
            claimed={x for x in range(q) if P[a,c][R[x]]==x}
            bank=[]
            for b in sorted(set(adj[a])&set(adj[c])-set(o)):
                hit={x for x in range(q) if P[b,c][P[a,b][R[x]]]==x}
                assert not claimed&hit
                claimed|=hit
                bank.append((b,len(hit)))
            assert sum(v for b,v in bank)+f==len(claimed)<=q
            score=(len(claimed),len(bank),len(o))
            if maxbank is None or score>maxbank[0]:
                maxbank=(score,o,f,bank)
    print('Complete branch pairs',len(P)//2)
    for m,h in bank_stats.items():
        print(m,'branches',sorted(h.items()))
    print('Best sampled 20-branch insertion bank',maxbank)
    _,o,_,_=maxbank
    a,c=o[-1],o[0]
    led=Counter()
    if len(adj[a])==len(adj[c])==d-1:
        for x in range(q):
            depart={c:x}
            for v,w in zip(o,o[1:]):
                depart[w]=P[v,w][depart[v]]
            h=depart[a]
            hits=[b for b in range(d) if b not in (a,c) and P[b,c][P[a,b][h]]==x]
            assert len(hits)+(P[a,c][h]==x)==1
            if not hits: led['shortcut']+=1
            elif hits[0] not in depart: led['fresh branch']+=1
            elif P[a,hits[0]][h]!=depart[hits[0]]: led['repeated branch, simple']+=1
            else: led['repeated vertex']+=1
        print('Complete 56-sheet insertion ledger', sorted(led.items()))
    C=55**21+1
    print('21/63 complement',C)
    witnessed=10000*C-sum(f*num for f,num in bank_stats[21].items())
    print('Distinct 63-cycles from sampled orders alone',witnessed)
    bnd=Fraction(3250,3)*(prod(range(37,58))//42)*(C-Fraction(56,37))
    print('Moore conditional global 63-cycle lower bound',-(-bnd.numerator//bnd.denominator))
    tri_fib={1,3,21,55}
    print('Triangular-Fibonacci fixed counts witnessed at order 21',sorted(tri_fib&set(bank_stats[21])))

if __name__=='__main__':
    controls()
    if len(sys.argv)>1:
        survey(sys.argv[1])
