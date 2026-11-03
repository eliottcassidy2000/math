#!/usr/bin/env python3
"""Exact standalone finite verification; no assertions and no third-party dependency."""
from collections import Counter, deque
from itertools import permutations
from pathlib import Path
import json, math

def require(ok, label):
    if not ok: raise RuntimeError(label)
def compose(p,q): return tuple(p[q[i]] for i in range(len(p)))
def cycles(p):
    unseen=set(range(len(p))); out=[]
    while unseen:
        x=min(unseen); c=[]
        while x in unseen:
            c.append(x);unseen.remove(x);x=p[x]
        require(x==c[0], 'not permutation')
        out.append(c)
    return out
A=tuple((5*x+8)%9 for x in range(9))
B=list(range(9)); word=[4,2,5,3,7]
for x,y in zip(word,word[1:]+word[:1]): B[x]=y
B=tuple(B)
L=tuple((527*x+44)%72 for x in range(72))
require(all(L[x]%8==(-x+4)%8 and L[x]%9==A[x%9] for x in range(72)), 'CRT')
require(Counter(map(len,cycles(L)))=={1:2,2:11,6:8},'72 cycle type')
require(all(L[x]==(34+23*(x-34))%72 for x in range(72)),'conjugacy')
require(all(((7*x+4)%8-3)//2%4==(-(x-3)//2-1)%4 for x in [1,3,5,7]),'balanced coordinates')
require(B[4]==2 and B[7]==4 and 4%3==7%3 and 2%3!=4%3,'polynomial obstruction')
def power(p,n):
    r=tuple(range(len(p)))
    for _ in range(n): r=compose(p,r)
    return r
C=power(compose(power(A,3),B),4)
require(C==(5,1,0,3,4,2,6,7,8), 'three-cycle word')
require(power(compose(A,power(B,3)),3)==tuple(range(9)), 'odd-sheet relation')
# Complete independent Cayley-graph enumeration. Composition conventions explicit.
identity=tuple(range(9)); group={identity}; queue=deque([identity])
while queue:
    p=queue.popleft()
    for q in (A,B):
        r=compose(q,p)
        if r not in group: group.add(r);queue.append(r)
require(len(group)==math.factorial(9)//2,'A9 order')
require(all(sum(p[i]>p[j] for i in range(9) for j in range(i+1,9))%2==0 for p in group),'all even')
# Independently enumerate the two-sheet lift, retaining parity of A applications.
lift={(0,identity)}; queue=deque(lift)
while queue:
    e,p=queue.popleft()
    for bit,q in ((1,A),(0,B)):
        r=(e^bit,compose(q,p))
        if r not in lift:lift.add(r);queue.append(r)
require(len(lift)==math.factorial(9), 'C2 x A9 order')
require(all((e,p) in lift for e in (0,1) for p in group),'product structure')
# Common-clock orbit law, independently all small length pairs.
for a in range(1,21):
    for b in range(1,21):
        p=tuple(((x//b+1)%a)*b+(x%b+1)%b for x in range(a*b))
        c=cycles(p)
        require(len(c)==math.gcd(a,b) and all(len(v)==math.lcm(a,b) for v in c),'clock law')
# CRT lift B and joint 72-state reachable sectors.
B72=tuple(next(y for y in range(72) if y%8==x%8 and y%9==B[x%9]) for x in range(72))
unseen=set(range(72)); sectors=[]
while unseen:
    start=min(unseen); reached={start}; q=deque([start])
    while q:
        x=q.popleft()
        for y in (L[x],B72[x]):
            if y not in reached: reached.add(y);q.append(y)
    unseen-=reached; sectors.append(sorted(reached))
require(sorted(map(len,sectors))==[9,9,18,18,18],'72 sectors')
# Literal six-state interpretation cannot descend: 4 and13 coincide but their successors do not.
literal=[4,2,5,3,7,13]; successor=dict(zip(literal,literal[1:]+literal[:1]))
require(4%9==13%9 and successor[4]%9!=successor[13]%9,'six-state nonfactor')
receipt={'status':'PASS','A9_generators':{'A':A,'B':B},'A9_order':len(group),'lift_order':len(lift),'mod72_cycles':cycles(L),'joint_mod72_sector_sizes':sorted(map(len,sectors)),'clock_pairs_tested':400,'sources':'explicit definitions; complete independent permutation enumeration'}
print(json.dumps(receipt,indent=2))
