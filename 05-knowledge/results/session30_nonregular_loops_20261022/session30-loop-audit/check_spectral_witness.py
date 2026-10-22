import json,itertools
from pathlib import Path
import sympy as sp
HERE=Path(__file__).resolve().parent
rows=json.loads((HERE.parent/'session30-almost-regular/isospectral_witnesses.json').read_text())
out=[]
for row in rows[:2]:
 n=row['n'];pair=[]
 for which in ('first','second'):
  w=row[which];S=sp.zeros(n);bit=0
  for i in range(n):
   for j in range(i+1,n):
    S[i,j]=1 if w['bits'][bit]=='1'else-1;S[j,i]=-S[i,j];bit+=1
  d=S*sp.ones(n,1);assert sorted(d)==[-1]*(n//2)+[1]*(n//2)
  H=sum(all(S[p[j],p[j+1]]==1 for j in range(n-1))for p in itertools.permutations(range(n)))
  cp=list(S.charpoly().all_coeffs());assert cp==row['charpoly'];assert H==int(w['H'])
  M=S/n;delta=d/n;t=(n*delta.T*(sp.eye(n)+M*M).inv()*delta)[0]
  if'tree'in w:assert str(t)==w['tree']
  pair.append({'side':which,'H':H,'charpoly':list(map(int,cp)),'tree':str(t)})
 assert pair[0]['H']!=pair[1]['H']
 out.append({'n':n,'pair':pair,'status':'PASS'})
(HERE/'spectral-witness-results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
