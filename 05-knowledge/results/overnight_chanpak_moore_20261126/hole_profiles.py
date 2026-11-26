"""Literal completion after fixing all fifty already assigned contacts."""
from collections import Counter
from pathlib import Path
from subprocess import check_output
from sys import executable, path
base=Path(__file__).resolve().parent.parent/'moore57_quadratic_tensors_20261125'
path.insert(0,str(base))
from three_stars import SLOTS
it=iter(map(int,check_output([executable,str(base/'nested_core_check.py'),'--board'],text=True).split()))
B=[h for h in SLOTS if h!=1]
assert next(it)==len(B) and B[0]==2
row,ref={},{}
for h in B:
    u,n=next(it),next(it)
    ref[h]=u
    row[h]={next(it):tuple(next(it) for _ in range(3)) for _ in range(n)}
marker=row.pop(2)[7]; ref.pop(2)
for h in row:
    row[h]={y:s for y,s in row[h].items() if all(a!=b for a,b in zip(s,marker))}
holes=[h for h in ref if ref[h]<0]
assert holes==[13,33]
fixed=[Counter(row[h][y][c] for h,y in ref.items() if y>=0) for c in range(3)]
assert all(len(x)==50 for x in fixed)
C=Counter()
w=[]
for x,s in row[13].items():
    for y,t in row[33].items():
        d=tuple(fixed[i][s[i]]+fixed[i][t[i]]+(s[i]==t[i]) for i in range(3))
        C[d]+=1; w.append((sum(d),x,y,d))
best=min(t[0] for t in w)
print('vacancy domain sizes',tuple(len(row[h]) for h in holes), 'profiles',sorted(C.items()))
print('minimum',best,'multiplicity',sum(p[0]==best for p in w),'minimizers',[p[1:] for p in w if p[0]==best])
