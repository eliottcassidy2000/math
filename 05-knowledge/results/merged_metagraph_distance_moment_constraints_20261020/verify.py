"""Producer all-size-small audit. Independent reviewer uses a separate implementation."""
from probe import build
from signed_redei_check import verify
from math import comb
import json, subprocess, sys
from pathlib import Path
root=Path(__file__).resolve().parent
rows=[]
for n in range(3,7):
    R,W,Q,mass,fi=build(n);m=R['m'];f=(n-1)//2;q=(m-f)//2
    fixed=[sum(comb(f,a)*comb(q,(j-a)//2) for a in range(f+1) if 0<=j-a<=2*q and (j-a)%2==0) for j in range(m+1)]
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    perm=[free.index((n-1-j,n-1-i)) for i,j in free]
    for x in range(1<<m):
        sx=sum(((x>>k)&1)<<p for k,p in enumerate(perm))
        assert fi[sx]==fi[x]
    for j in range(m+1):
        for u,s in enumerate(R['mass']):
            assert R['Q'][j][u][u]%8==(comb(m,j) if s%2 else 4*fixed[j])%8
    z=[sum((-1)**x.bit_count() for x in fib) for fib in R['fibres']]
    assert sum(z)==0 and all(z)
    for u,s in enumerate(R['mass']):
        assert z[u]%(2 if s%2 else 4)==(1 if s%2 else 2)
        assert R['Q'][m][u][u]%(8 if s%2 else 32)==(1 if s%2 else 4)
        for v in range(len(z)):
            assert R['Q'][m][u][v]==z[u]*z[v]
            assert R['Q'][m][u][v]==R['Q'][m][0][u]*R['Q'][m][0][v]
    rows.append(R)
    print('PASS',n,'nodes',len(mass),'ranks',R['rank_Q'],'leak_rank',R['leak_rank'],flush=True)
json.dump(rows[:3],open(root/'audit.json','w'),indent=2)
json.dump(rows[-1],open(root/'audit_n6.json','w'),indent=2)
signed=verify(rows);json.dump(signed,open(root/'signed_redei_check.json','w'),indent=2)
print('PASS signed Redei bridge:',len(signed),'merged nodes n=3..6')
subprocess.run([sys.executable,str(root/'hostile.py')],check=True,stdout=open(root/'hostile.out','w'))
print('PASS hostile controls and exact negative principal minor')
