#!/usr/bin/env python3
"""Independent complete binary-continuation reachability, widths 1..8."""
from itertools import product
from collections import deque,defaultdict
import json
counts=defaultdict(int); widths=[]
def check(p,kind,detail=None):
    counts[kind]+=1
    if not p: raise RuntimeError((kind,detail))
def step(x,u):return tuple(abs(x[i]-x[i+1]) for i in range(len(x)-1))+(abs(x[-1]-u),)
def edge(x):
    ans=[]
    while x:
        ans.append(x[0]);x=tuple(abs(x[i]-x[i+1])for i in range(len(x)-1))
    return tuple(ans)
def good(x):return max(edge(x),default=0)<=1
for F in range(1,9):
    states=list(product(range(4),repeat=F));ids={x:i for i,x in enumerate(states)}
    preds=[[]for _ in states]
    for i,x in enumerate(states):
        for u in (0,1):preds[ids[step(x,u)]].append(i)
    bad={i for i,x in enumerate(states)if x[0]>=2};q=deque(bad)
    while q:
        i=q.popleft()
        for j in preds[i]:
            if j not in bad:bad.add(j);q.append(j)
    shield=0;ngood=0
    for i,x in enumerate(states):
        counts['prefixes']+=1
        if not good(x):continue
        ngood+=1
        check(good(x+(1,)),'append_one',x)
        check(good(x+(2,))==any(x),'append_two_iff',x)
        last=next((v for v in reversed(x)if v),0)
        if last<=2:check(good(x+(0,)),'append_zero',x)
        if x[-1] in (1,2):
            check(i not in bad,'all_continuations_safe',x);shield+=1
    widths.append({'width':F,'states':len(states),'good_prefixes':ngood,'shield_prefixes':shield,'can_fail_states':len(bad)})
hostile=(1,2,3,4,1)
check(edge(hostile)==(1,1,0,0,2),'alphabet_four_hostile')
print(json.dumps({'status':'PASS','checks':dict(counts),'widths':widths,'hostile':{'row':hostile,'left_edge':edge(hostile)}},indent=2))
