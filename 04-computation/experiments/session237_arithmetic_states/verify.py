#!/usr/bin/env python3
import sys as _sys
if "--fail" in _sys.argv:
    raise RuntimeError("deliberate checker failure")
"""Standalone exact checks; run from any cwd. Writes only stdout."""
from collections import deque
import json

def Q(m):
    if m == 0: return 0
    a,b=divmod(m,2)
    return Q(a)+(a+a%2 if b else 0)

def step(s,b,M):
    r,q=s
    return ((2*r+b)%(2*M),(q+(r+r%2)//2*b)%M)

def output(s,M): return (s[0]-s[1])%M

def hp(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if dp[mask][v]:
                available=rows[v]&~mask
                while available:
                    bit=available&-available; available-=bit
                    dp[mask|bit][bit.bit_length()-1]+=dp[mask][v]
    return sum(dp[-1])

def literal(m):
    rows=[2,0]
    for i in range(m.bit_length()):
        n=len(rows); bridge=(m>>i)&1
        for v in range(n):
            if not (bridge and v==0): rows[v]|=1<<n
        rows.append(1 if bridge else 0)
    return hp(rows)

def require(x,msg):
    if not x: raise RuntimeError(msg)

out={'status':'FINITE-EXACT, corroborating all-order proofs','automata':[]}
for M in range(2,41):
    reached={(0,0):0}; todo=deque([(0,0)])
    while todo:
        s=todo.popleft()
        for b in (0,1):
            t=step(s,b,M)
            if t not in reached:
                reached[t]=reached[s]+1;todo.append(t)
    require(len(reached)==2*M*M,('reachability',M))
    sig={(output(s,M),output(step(s,0,M),M),output(step(s,1,M),M)) for s in reached}
    require(len(sig)==len(reached),('distinguishability',M))
    out['automata'].append({'M':M,'states':len(reached),'max_shortest_input_length':max(reached.values())})
for m in range(65536):
    for M in (2,3,8,9,15,32):
        s=(0,0)
        for b in bin(m)[2:]:s=step(s,int(b),M)
        require(s==(m%(2*M),(Q(m)//2)%M),('Q-state',m,M))
for m in range(256): require(literal(m)==1+2*m-Q(m),('literal',m))
# Joint images, independently enumerate lattice versus reachable automaton states.
seeds=[(1,2,0,1),(5,6,2,3),(9,11,3,5)]
for M in range(2,21):
    joint=set()
    for r in range(2*M):
        for q in range(M):
            joint.add(tuple(((h-1+R*r-e*(r%2))//2-d*q)%M for h,R,e,d in seeds))
    predicted=set()
    for b in (0,1):
        for j in range(M):
            for q in range(M):
                predicted.add(tuple(((h-1+(R-e)*b)//2+R*j-d*q)%M for h,R,e,d in seeds))
    require(joint==predicted,('joint',M))
    require(len({(a,c) for a,b,c in joint})==M*M,('unimodular pair',M))
    require(all((b-2-3*a)%M in (0,(-1)%M) for a,b,c in joint),('rank-one pair',M))
# Independent literal block pumping: preserve m modulo 2M, fill every q residue.
pumping_checks=0
for M in range(2,17):
    u=2*M;k=0
    while u%2==0:u//=2;k+=1
    b=max(k,2)
    while pow(2,b,u)!=1%u:b+=1
    x=(1<<b)-1
    for a in range(2*M):
        A=k+a.bit_length(); found=set();m=a
        for t in range(1,M+1):
            m+=x<<(A+(t-1)*(b+k))
            require(m%(2*M)==a,('pump m',M,a,t))
            require((Q(m)-Q(a)-t*((1<<b)-2))%(2*M)==0,('pump Q',M,a,t))
            found.add((Q(m)//2)%M);pumping_checks+=1
        require(len(found)==M,('pump full',M,a))
out['literal_pumping_checks']=pumping_checks
# Hostile controls: retaining only output, or (m mod M,q), is unsound.
M=4
require(output((0,0),M)==output((4,0),M),'hostile same current output')
require(output(step((0,0),1,M),M)!=output(step((4,0),1,M),M),'lost parity lift')
root_values={1+2*m-Q(m) for m in range(10000)}
require(7 not in root_values and 21 not in root_values,'forbidden H controls')
out['recurrence_checks']=65536*6
out['literal_tournaments']=256
out['joint_moduli']=list(range(2,21))
out['hostile_controls']=['m modulo M loses the parity lift at M=4','scalar output alone loses continuation','full modular projections do not attain 7 or 21']
print(json.dumps(out,indent=2))
