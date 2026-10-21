import gzip, os
def _open_census(n):
    p=f'census_n{n}.txt'
    return open(p) if os.path.exists(p) else gzip.open(p+'.gz','rt')
# analyze_census.py -- spectra of H, C, tc over all tournament classes n <= 9 (mac-mini-2026-10-21-fanopsp)
import sys, math
from collections import Counter, defaultdict
def load(n):
    rows=[]
    for line in _open_census(n):
        H,C,s,A,AA,c3=map(int,line.split()); rows.append((H,C,s,A,AA,c3))
    return rows
for n in range(3,10):
    R=load(n); fact=math.factorial(n)
    Hc=Counter(r[0] for r in R)
    Hlab=Counter(); 
    for r in R: Hlab[r[0]]+=fact//r[3]
    assert sum(Hlab.values())==2**(n*(n-1)//2)
    Cc=Counter(r[1] for r in R)
    strongH=Counter(r[0] for r in R if r[2])
    strongC=Counter(r[1] for r in R if r[2])
    tc=Counter(r[0]//r[3] for r in R)
    assert all(r[0]%r[3]==0 for r in R)
    maxH=max(Hc); odd_missing=[h for h in range(1,maxH+1,2) if h not in Hc]
    print(f'== n={n}: classes {len(R)}, H values {len(Hc)}, max H {maxH}, max C {max(Cc)}')
    print('   odd H missing up to max:', odd_missing[:40], '...' if len(odd_missing)>40 else '')
    smin=min(strongH) if strongH else None
    print('   strong H range', smin, max(strongH) if strongH else None, ' strong odd missing in range:',
          [h for h in range(smin,max(strongH)+1,2) if h not in strongH][:30] if strongH else None)
    cmissing=[c for c in range(0,max(Cc)+1) if c not in Cc]
    print('   C values missing up to max:', cmissing[:40], '... (#%d)'%len(cmissing))
    print('   tc missing odd up to 60:', [t for t in range(1,61,2) if t not in tc])
    if n<=6:
        print('   H class mult:', sorted(Hc.items()))
        print('   H labeled mult:', sorted(Hlab.items()))
        print('   C class mult:', sorted(Cc.items()))
