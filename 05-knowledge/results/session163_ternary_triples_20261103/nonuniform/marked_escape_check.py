#!/usr/bin/env python3
"""Endpoint-subset controls for the R4 nonuniform marked-interface escape."""
from itertools import combinations
from math import gcd

def need(p,msg):
 if not p:raise ValueError(msg)

def H(rows):
 n=len(rows);D=[[0]*n for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for mask in range(1<<n):
  for v in range(n):
   c=D[mask][v];z=rows[v]&~mask
   while c and z:
    b=z&-z;j=b.bit_length()-1;D[mask|b][j]+=c;z-=b
 return sum(D[-1])

def add(rows,beat):
 n=len(rows)
 return tuple(r|((1<<n) if i not in beat else 0) for i,r in enumerate(rows))+(sum(1<<i for i in beat),)

def reorder(rows,order):
 return tuple(sum(((rows[a]>>b)&1)<<j for j,b in enumerate(order)) for a in order)

def state(rows,port):
 p=reorder(rows,tuple(i for i in range(len(rows)) if i!=port)+(port,))
 h=H(p);a=add(p,{len(p)-1});ha=H(a);aa=add(a,{len(a)-1});haa=H(aa)
 need((ha-h)%2==0 and (haa-ha)%2==0,'observer parity')
 return h,(ha-h)//2,(haa-ha)//2

def strong(rows):
 n=len(rows)
 for root in range(n):
  seen={root}; stack=[root]
  while stack:
   u=stack.pop()
   for v in range(n):
    if rows[u]>>v&1 and v not in seen: seen.add(v);stack.append(v)
  if len(seen)!=n:return False
 return True

def join(rows,block):
 n=len(rows);mask=(1<<len(block))-1;r=list(rows)
 for i in range(n):r[i]|=(mask^(1 if i==n-1 else 0))<<n
 return tuple(r)+tuple((x<<n)|((1<<(n-1)) if i==0 else 0) for i,x in enumerate(block))

r=(996,993,994,999,975,912,800,576,128,256);B={0,1,2}
need(strong(r),'R4 strong')
need(all(len({(r[x]>>v)&1 for v in B})==1 for x in range(len(r)) if x not in B),'triangle module')
need(state(r,9)==(291,120,225),'R4 state')
expected=[(291,345,345),(1965,1103,1564),(3445,1781,2633),(4731,2379,3552)]
for k in range(4):
 for beat in combinations(B,k):
  rows=add(r,set(beat));actual=state(rows,9)
  need(actual==expected[k],f'bad state {beat}: {actual}')
  if k:need(strong(rows),'escape strong')
need({v[0] for v in expected}=={291,1965,3445,4731},'exact H set')
need(all(gcd(v[0],15)>1 for v in expected),'unexpected unit mod 15')
need(gcd(3445,4731)==1,'coprime strong cases')
print('retained-port states by |z-out-neighbors in triangle|:',expected)
print('exact H set has no unit mod 15; 3445 and 4731 are coprime strong counts')
unrestricted=add(r,{0,3});need(strong(unrestricted),'unrestricted strong')
need(state(unrestricted,9)==(6011,2980,4531),'unrestricted state')
need(gcd(6011,15)==1,'unrestricted unit')
print('unrestricted z->0,3: strong state (6011,2980,4531), a unit mod 15')
single=add(r,{0});single=reorder(single,tuple(i for i in range(len(single)) if i!=9)+(9,))
q2=join(single,(8,1,3,6));need(strong(q2),'q2 strong');need(H(q2)==28033,'q2 exact count')
print('one-edge escape after literal q2: 28033')
