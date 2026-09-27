"""Finite exact Ackermann carrier controls for reset_20260926_logic.md.
No claim that finite tests prove any axiom scheme or Collatz convergence.
"""
from functools import lru_cache
from itertools import combinations
from collections import Counter
import json
counts=Counter()
def need(ok,label):
    counts[label]+=1
    if not ok:raise RuntimeError(label)
def bits(n):
    out=[]
    while n:
        b=n&-n;out.append(b.bit_length()-1);n-=b
    return out
@lru_cache(None)
def decode(n):return frozenset(decode(i) for i in bits(n))
@lru_cache(None)
def encode(x):return sum(1<<encode(y) for y in x)
def submasks(n):
    m=n
    while True:
        yield m
        if not m:break
        m=(m-1)&n
def powercode(n):return sum(1<<s for s in submasks(n))
def closure_nodes(n):
    result=set();stack=bits(n)
    while stack:
        x=stack.pop()
        if x in result:continue
        result.add(x);stack.extend(bits(x))
    return result
def literal_closure(x):
    answer=set(x);todo=list(x)
    while todo:
        y=todo.pop()
        for z in y:
            if z not in answer:answer.add(z);todo.append(z)
    return frozenset(answer)
@lru_cache(None)
def rank(n):return 0 if not n else 1+max(rank(i) for i in bits(n))
for n in range(4096):
    need(encode(decode(n))==n,'Ackermann roundtrip')
    for i in bits(n):need(i<n,'membership numerical descent')
    tc=sum(1<<i for i in closure_nodes(n))
    need(decode(tc)==literal_closure(decode(n)),'transitive closure independent path')
    need(all(set(bits(i))<=set(bits(tc)) for i in bits(tc)),'closure transitive')
for n in range(64):
  for m in range(64):
    need(decode(n|(1<<m))==decode(n)|{decode(m)},'adjunction graph transport')
for n in range(32):
    xs=list(decode(n))
    literal=frozenset(frozenset(c) for k in range(len(xs)+1) for c in combinations(xs,k))
    need(decode(powercode(n))==literal,'powerset independent path')
    need(powercode(n)>n,'powerset code grows')
for n in range(128):
    B=(1<<(n+1))-1
    need((B>>n)&1,'container contains source')
    for i in range(n+1):
        need(all((B>>j)&1 for j in bits(i)),'container transitive')
        need(all((B>>j)&1 for j in submasks(i)),'container subset closed')
    # Full power-set object exceeds every member code of the finite B.
    need(powercode(n)>n,'container not powerset-object closed')
need(rank(3)==2 and rank(5)==3,'Collatz raises HF rank')
need(not ((3>>5)&1) and not ((5>>3)&1),'Collatz edge is not membership')
need(sum(1<<i for i in closure_nodes(27))==31,'TC may raise ordinary integer')
chain=[];n=0
for k in range(7):
    need(rank(n)==k,'singleton tower rank')
    nodes=closure_nodes(n)|{n}
    need(len(nodes)==k+1,'singleton tower unbounded depth')
    chain.append(dict(rank=k,integer_bitlength=n.bit_length(),graph_nodes=len(nodes)))
    if k<6:n=1<<n
print(json.dumps(dict(status='FINITE-EXACT coding controls only',universes=dict(roundtrip=4096,adjunction_pairs=4096,powersets=32,containers=128),singleton_tower=chain,example_27=dict(top_members=bits(27),closure_members=sorted(closure_nodes(27)),closure_code=31),checks=dict(sorted(counts.items())),total=sum(counts.values()),result='PASS'),indent=2))
