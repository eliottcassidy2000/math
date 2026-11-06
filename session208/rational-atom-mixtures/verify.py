import json,itertools
from fractions import Fraction as F
from pathlib import Path
root=Path(__file__).resolve().parent
source=json.loads((root/'census_session207.json').read_text())
D={}
def matrix(n,code):
 a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  a[i][j]=bool(code>>k&1);a[j][i]=not a[i][j]
 return a
def check(n,code):
 a=matrix(n,code); reach=[sum(1<<j for j in range(n) if a[i][j])|1<<i for i in range(n)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 assert all(x==(1<<n)-1 for x in reach)
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for j in range(n):
   if mask>>j&1:
    for k in range(n):
     if not mask>>k&1 and a[j][k]:dp[mask|1<<k][k]+=dp[mask][j]
 return sum(dp[-1])
for row in source['orders']:
 n=row['n']
 for h in row['strong_spectrum']:
  if h<=1:continue
  t=h;a=b=0
  while t%3==0:t//=3;a+=1
  while t%5==0:t//=5;b+=1
  if t==1 and h not in D:
   code=row['witness_codes'][str(h)]
   assert check(n,code)==h
   D[h]=(a,b,n,code)
assert len(D)==15
slacks={h:11*n-16*a-24*b for h,(a,b,n,_) in D.items()}
assert min(slacks.values())==0
assert {h for h in slacks if slacks[h]==0}=={375,405}
for h in (375,405):
 a,b,n,code=D[h];adj=matrix(n,code)
 count=sum(all(adj[p[i]][p[i+1]] for i in range(n-1)) for p in itertools.permutations(range(n)))
 assert count==h
assert 375**3*405**2==15**11
assert F(3,11)*D[375][2]+F(2,11)*D[405][2]==F(40,11)
f={(0,0):0}
for a in range(12):
 for b in range(12):
  if a+b:f[a,b]=min(f[a-x,b-y]+n for x,y,n,_ in D.values() if x<=a and y<=b)
residues=[f[k,k] for k in range(11)]
assert residues==[0,5,8,13,15,20,23,27,30,35,38]
assert f[11,11]==40
result={'dictionary':D,'dual_price':['16/11','24/11'],'slacks':slacks,'lambda':'40/11','block_exponent':11,'block_order':40,'residue_costs':residues,'scope':'Dictionary exact; global upper bound only.','checks':'All passed'}
(root/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
