#!/usr/bin/env python3
"""Exact SCC-profile reachability, with explicit source/clock reconstruction."""
from collections import deque
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent

def act(state,profile):
    F,q,g,h=state
    H,E,p,Q,s,t,f=profile
    return H*F+p*q+(H-s)*g+(f-H)*h,E*q+Q*h,(H-t)*h,H*h

def search(bank,Q=32,Hcap=31,Fcap=64):
    root=(1,0,0,1)
    previous={root:None};queue=deque([root]);transitions=0;eligible={}
    candidates=[(i,b['profile']) for i,b in enumerate(bank)
                if b['profile'][0]<=Hcap and b['profile'][3]<=Q and b['profile'][6]<=Fcap]
    while queue:
        state=queue.popleft()
        _,q,_,h=state
        if (q,h) not in eligible:
            eligible[q,h]=[(i,b) for i,b in candidates if b[0]*h<=Hcap and b[1]*q+b[3]*h<=Q]
        for i,b in eligible[q,h]:
            new=act(state,b);F,q,g,h=new
            if F>Fcap or q>Q or h>Hcap:continue
            if not (0<=g<=h<=F and (q==0 or h<=q)):
                raise AssertionError(('state bounds',new))
            transitions+=1
            if new not in previous:
                previous[new]=(state,i);queue.append(new)
    return previous,transitions

def word(previous,state):
    result=[]
    while previous[state] is not None:
        state,i=previous[state];result.append(i)
    return result[::-1]

def graph(bank,indices):
    # Empty parent plus its marked vertex x=0.
    rows=[0]
    for i in indices:
        b=bank[i];n=len(rows);m=len(b['rows']);rows.extend([0]*m)
        for old in range(1,n):rows[old]|=((1<<m)-1)<<n
        for v,row in enumerate(b['rows']):
            rows[n+v]|=row<<n
            if b['mask']>>v&1:rows[0]|=1<<(n+v)
            else:rows[n+v]|=1
    return rows

def main():
    bank=json.loads((HERE/'strong_profiles.json').read_text())['profiles']
    previous,transitions=search(bank)
    hits={}
    for state in sorted(previous):
        F,q,g,h=state
        if q==0 or q&(q-1):continue
        for r in range(1,7):
            N=q*(1<<r)+F+g-q
            if N<=65 and N not in hits:
                hits[N]={'state':state,'r':r,'word':word(previous,state)}
    targets=[]
    for q,total in [(32,15),(16,31)]:
        targets.append({'q':q,'F_plus_g':total,'states':[s for s in sorted(previous) if s[1]==q and s[0]+s[2]==total]})
    print(json.dumps({'states':len(previous),'transitions':transitions,
       'bounds':{'q':32,'h':31,'F':64},'profiles':len(bank),
       'dyadic_sources':[{'state':s,'word':word(previous,s)} for s in sorted(previous) if s[1] and not s[1]&(s[1]-1)],
       'zero_tail_hits_through65':hits,
       'missing_odd_3_through63':[n for n in range(3,64,2) if n not in hits],
       'target47':targets},indent=2,sort_keys=True))

if __name__=='__main__':main()
