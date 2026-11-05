# per uniquely Hamiltonian parent U of order m: number of DISTINCT two-cycle child classes (one-vertex extensions with hc = 2),
# cross-tabulated with U's (sole exits, sole entries); also the total number of (U, child) incidences.
#!/usr/bin/env python3
"""perparent.py (klein-2026-11-04): for every uniquely Hamiltonian parent U of order m (data/cls_{m}.txt, hc = 1 lines) the
number of DISTINCT two-cycle child classes among its one-vertex extensions (tcgen2 with parent index, labelg), cross-tabulated
with U's (sole exits, sole entries); and the number of (parent, child) incidences = number of (Y, v) pairs with Y - v uniquely
Hamiltonian.  Needs scripts/tcgen2 compiled and nauty labelg.  Usage: perparent.py M1 M2."""
import subprocess, sys, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
LABELG='/opt/homebrew/bin/labelg'
def rows_from_d6(s):
    n = ord(s[1]) - 63; bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1): bits.append((v >> k) & 1)
    rows = [0] * n
    for i in range(n):
        for j in range(n):
            if bits[i * n + j]: rows[i] |= 1 << j
    return rows
for m in range(int(sys.argv[1]), int(sys.argv[2])+1):
    U=[l.split()[1] for l in open(os.path.join(DATA, 'cls_%d.txt' % m)) if l.split()[0]=='1']
    raw=subprocess.run([os.path.join(HERE, 'tcgen2'), str(m), '2'],input='\n'.join(U)+'\n',capture_output=True,text=True).stdout.split('\n')
    recs=[l.split() for l in raw if l.strip()]
    can=subprocess.run([LABELG,'-q'],input='\n'.join(r[1] for r in recs)+'\n',capture_output=True,text=True).stdout.split()
    kids={}  # parent index -> set of canonical hc-2 children
    for r,c in zip(recs,can):
        if r[0]=='2': kids.setdefault(int(r[2]),set()).add(c)
    feat={}
    for i,d6 in enumerate(U,1):
        rows=rows_from_d6(d6); od=[bin(r).count('1') for r in rows]; ind=[sum((rows[a]>>b)&1 for a in range(m)) for b in range(m)]
        feat[i]=(od.count(1),ind.count(1))
    tab=Counter(); inc=0; allkids=set()
    for i in range(1,len(U)+1):
        k=len(kids.get(i,())); tab[(feat[i],k)]+=1; inc+=k; allkids|=kids.get(i,set())
    print("parent order %d (children order %d): %d UH parents; distinct two-cycle children overall %d; (U,child) incidences %d"%(m,m+1,len(U),len(allkids),inc))
    for (f,k),v in sorted(tab.items()): print("    parent (sole exits, sole entries)=%s with %2d distinct two-cycle children: %d parents"%(f,k,v))
