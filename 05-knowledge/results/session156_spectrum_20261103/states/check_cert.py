# Independent subset-endpoint witness verifier. No transfer-automaton import.
import json
from pathlib import Path
w=json.loads(Path(__file__).with_name('witnesses.json').read_text())
w={int(h):v for h,v in w.items()}
expected=set(range(3,302,2))-{7,21}
if set(w)!=expected: raise SystemExit('wrong target set')
for h,(types,ab) in w.items():
    m=len(types); n=m+2
    a=[[0]*n for _ in range(n)]
    for i in range(m):
        for j in range(i+1,m): a[i][j]=1
        for k in range(2):
            if types[i]>>k & 1: a[i][m+k]=1
            else: a[m+k][i]=1
    a[m+(1-ab)][m+ab]=1
    reach=[row[:] for row in a]
    for i in range(n): reach[i][i]=1
    for z in range(n):
      for x in range(n):
       for y in range(n): reach[x][y] |= reach[x][z] and reach[z][y]
    if not all(all(row) for row in reach): raise SystemExit(f'{h}: not strong')
    d=[[0]*n for _ in range(1<<n)]
    for x in range(n):d[1<<x][x]=1
    for mask in range(1,1<<n):
      for last in range(n):
       for nxt in range(n):
        if not(mask>>nxt&1) and a[last][nxt]:d[mask|1<<nxt][nxt]+=d[mask][last]
    if sum(d[-1])!=h:raise SystemExit(f'{h}: count mismatch')
print('checked all 148 strong witnesses; each has FVS size at most 2; max base length',max(len(v[0]) for v in w.values()))
for m in range(1,281):
    if m%2 and m>140: continue
    eligible=set(range(m)) if m%2 else set(range(1,m,2))
    if {h%m for h in w}!=eligible: raise SystemExit(f'{m}: residue coverage failure')
print('checked all odd moduli <=140 and even moduli <=280')
