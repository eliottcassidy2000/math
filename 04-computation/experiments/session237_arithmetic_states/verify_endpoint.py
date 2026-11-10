#!/usr/bin/env python3
import sys as _sys
if "--fail" in _sys.argv:
    raise RuntimeError("deliberate checker failure")
"""Exact endpoint observer; independent physical endpoint DP, no producer imports."""
from collections import deque
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def Q(n):
    ans=0
    while n:
        a,b=divmod(n,2)
        if b:ans+=a+a%2
        n=a
    return ans

def transition(state,b,M):
    r,q=state;L=M if M%2==0 else 2*M
    return ((2*r+b)%L,(q+b*(r+r%2))%M)
def output(state,M):return (state[0]-state[1])%M

def physical(m):
    bits=bin(m)[2:] if m else ''
    l=len(bits);a=l;b=l+1;n=l+2;rows=[0]*n
    for i in range(l+1):
        for j in range(i+1,l+1):rows[i]|=1<<j
    rows[a]|=1<<b
    for i,bit in enumerate(bits):
        if bit=='1':rows[b]|=1<<i
        else:rows[i]|=1<<b
    # Backward endpoint recurrence: dp[S][v] sums dp[S-v][u], u->v.
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        if mask&(mask-1)==0:continue
        for v in range(n):
            if mask>>v&1:
                old=mask^(1<<v)
                dp[mask][v]=sum(dp[old][u] for u in range(n) if old>>u&1 and rows[u]>>v&1)
    return dp[-1][a],sum(dp[-1])

records=[]
for M in range(1,51):
    if M==1:
        records.append({'M':1,'states':1,'max_input_length':0});continue
    states={(0,0):0};todo=deque(states)
    while todo:
        s=todo.popleft()
        for b in (0,1):
            t=transition(s,b,M)
            if t not in states:states[t]=states[s]+1;todo.append(t)
    expect=M*M//2 if M%2==0 else 2*M*M
    need(len(states)==expect,('reachable',M))
    depth=1 if M%2==0 else 2
    signatures={(output(s,M),output(transition(s,0,M),M)) if depth==1 else (output(s,M),output(transition(s,0,M),M),output(transition(s,1,M),M)) for s in states}
    need(len(signatures)==expect,('minimal',M))
    records.append({'M':M,'states':expect,'max_input_length':max(states.values())})
for m in range(65536):
    e=m-Q(m)
    need(2*m-Q(2*m)==m+e,('even recurrence',m))
    need(2*m+1-Q(2*m+1)==e+1-m%2,('odd recurrence',m))
    for M in (2,3,4,6,7,8,15,16):
        s=(0,0)
        for bit in bin(m)[2:]:s=transition(s,int(bit),M)
        need(output(s,M)==e%M,('observer',m,M))
for m in range(256):
    e,h=physical(m);need((e,h)==(m-Q(m),1+2*m-Q(m)),('literal endpoint',m))
# Odd modulus needs the parity lift; same reduced state, different digit-1 outputs.
for M in (3,5,7,9):
    s=(0,0);t=(M,0)
    need(output(s,M)==output(t,M),'hostile initial')
    need(output(transition(s,1,M),M)!=output(transition(t,1,M),M),'hostile odd lift')
# M=2 is genuinely only two states, rather than the normalized-root count eight.
need(records[1]['states']==2,'small even boundary')
print(json.dumps({'status':'PASS','automata':records,'recurrence_inputs':65536,'observer_checks':65536*8,'literal_endpoint_tournaments':256,'hostile_controls':['odd modulus needs parity lift','M=2 has exactly two states','M=1 is constant']},indent=2))
