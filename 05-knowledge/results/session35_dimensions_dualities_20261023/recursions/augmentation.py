import itertools, math,json
from functools import lru_cache
import sympy as s
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
@lru_cache(None)
def pf(a):
 n=int(len(a)**.5)
 if n==0:return 1
 return sum((-1)**(j+1)*a[j]*pf(tuple(a[r*n+c] for r in range(1,n) if r!=j for c in range(1,n) if c!=j)) for j in range(1,n))
def calc(M):
 n=len(M); K=[[M[i][j]-M[i][-1]+M[j][-1] for j in range(n-1)] for i in range(n-1)]
 p=pf(tuple(sum(K,[])));co=sum((-1)**i*pf(tuple(M[r][c] for r in range(n) if r!=i for c in range(n) if c!=i)) for i in range(n))
 if p!=co:raise RuntimeError('Pf sign')
 d1=math.gcd(*sum(K,[])); ds=[d1]
 if n==5:ds.append(abs(p)//d1)
 if n==7:
  g4=math.gcd(*(pf(tuple(K[i][j] for i in sub for j in sub)) for sub in itertools.combinations(range(6),4)))
  ds.extend([g4//d1,abs(p)//g4])
 out=[sum(1<<j for j in range(n) if M[i][j]>0) for i in range(n)]
 dp=[[0]*n for _ in range(1<<n)];sp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=sp[1<<i][i]=1
 for mask in range(1,1<<n):
  for j in range(n):
   if not mask>>j&1:continue
   for k in range(n):
    if not mask>>k&1 and out[j]>>k&1:
     dp[mask|1<<k][k]+=dp[mask][j]
     sp[mask|1<<k][k]+=(-1)**((mask>>(k+1)).bit_count())*sp[mask][j]
 return {'pf':p,'smith_alt':ds,'H':sum(dp[-1]),'signed_R':sum(sp[-1]),'out':out,'K':K}
for n in [3,5]:
 free=[(i,j) for i in range(n) for j in range(i+2,n)]; seen={};bad=None;profiles={}
 for mask in range(1<<len(free)):
  M=[[0]*n for _ in range(n)]
  for i in range(n-1):M[i][i+1]=1;M[i+1][i]=-1
  for b,(i,j) in enumerate(free):M[i][j]=1 if mask>>b&1 else -1;M[j][i]=-M[i][j]
  c=calc(M);key=tuple(c['smith_alt']);profiles.setdefault(key,set()).add(c['H'])
  if key in seen and seen[key]['H']!=c['H'] and bad is None:bad=[seen[key],c]
  seen.setdefault(key,c)
 print('n',n,'profiles',{str(k):sorted(v) for k,v in profiles.items()},'first_collision',bad)
 if bad:
  for c in bad:
   K=s.Matrix(c['K']);diag=smith_normal_form(K,domain=ZZ)
   print('independent_SNF',list(diag.diagonal()),'detK',K.det())
for n,kind in [(5,'cyclic'),(7,'paley')]:
 residues={1,2} if n==5 else {1,2,4}
 M=[[0 if i==j else (1 if (j-i)%n in residues else -1) for j in range(n)] for i in range(n)]
 c=calc(M);A=s.Matrix([[(int(M[i][j])+1)//2 if i!=j else 0 for j in range(n)] for i in range(n)])
 print(kind,c,'detI2A',(s.eye(n)+2*A).det())
