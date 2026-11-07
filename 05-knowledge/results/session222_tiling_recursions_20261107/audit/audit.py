#!/usr/bin/env python3
"""Independent exact Session222 audit; no producer code imported."""
import itertools,json
from pathlib import Path

def H(n, flips=()):
    F=set(flips)
    incoming=[sum(1<<u for u in range(n) if u!=v and ((u>v)^((max(u,v),min(u,v)) in F))) for v in range(n)]
    dp={1<<v:{v:1} for v in range(n)}
    for S in range(1,1<<n):
        if S not in dp:continue
        for v in range(n):
            if S>>v&1:continue
            count=sum(c for u,c in dp[S].items() if incoming[v]>>u&1)
            if count:dp.setdefault(S|1<<v,{})[v]=count
    return sum(dp[(1<<n)-1].values())

def need(x):
    if not x:raise RuntimeError('audit failure')

out={"scope":"All fixed descending-path tilings n=3..6; all unfiltered labelled tournaments n=1..5; geometric shells m=0..97", "shells":[],"padding":[],"parity":[]}
for m in range(98):
    Q=lambda t:{(a,b) for a in range(t+1) for b in range(t+1-a)}
    old=Q(m); big=Q(m+2); inner={(a+1,b+1) for a,b in old}
    L=big-inner; D=big-old
    # Independently use ordered segment matches, not producer shell_map.
    mapping={p:(p[0]-1,p[1]-1) for p in inner}
    mapping.update({(0,b):(b,m+2-b) for b in range(m+3)})
    mapping.update({(a,0):(a-1,m+2-a) for a in range(1,m+3)})
    need(set(mapping)==big and set(mapping.values())==big)
    need({mapping[p] for p in L}==D and len(L)==2*m+5)
out['shells']={'orders':98,'minimum_n':3,'maximum_n':100}
for n in range(3,7):
    cells=[(x,y) for x in range(n) for y in range(x-1)]
    mismatch=0;first=None
    for mask in range(1<<len(cells)):
        F=[p for i,p in enumerate(cells) if mask>>i&1]
        h=H(n,F); leg=H(n+2,[(x+1,y+1) for x,y in F]); diag=H(n+2,[(x+2,y) for x,y in F])
        need(h==leg and h%2==1)
        if h!=diag:
            mismatch+=1
            if first is None:first={'flips_zero_based':F,'H':[h,leg,diag]}
    out['padding'].append({'n':n,'total':1<<len(cells),'diagonal_H_changes':mismatch,'first_hostile':first})
for n in range(1,6):
    cells=list(itertools.combinations(range(n),2));counts=[]
    for mask in range(1<<len(cells)):
        counts.append(H(n,[(b,a) for i,(a,b) in enumerate(cells) if mask>>i&1]))
    need(all(h>0 and h%2 for h in counts))
    out['parity'].append({'n':n,'labelled_tournaments':len(counts),'min':min(counts),'max':max(counts)})
out['status']='PASS'
Path(__file__).with_name('RESULT.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
