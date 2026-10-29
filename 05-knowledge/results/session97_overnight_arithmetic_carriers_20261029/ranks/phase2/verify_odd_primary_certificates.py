"""Independent verifier in short Weierstrass coordinates, dense F3 elimination."""
import json,math,copy
from pathlib import Path
from fractions import Fraction as Q
P=Path(__file__).resolve().parent;checks=0

def must(v):
 global checks
 checks+=1
 if not v:raise ValueError('odd-primary certificate rejected')
def pr(p):return p>1 and all(p%d for d in range(2,math.isqrt(p)+1))
def coord(x,p):
 x=Q(x);must(x.denominator%p!=0);return x.numerator*pow(x.denominator,-1,p)%p

def short(a,p):
 a1,a2,a3,a4,a6=a;b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6;c4=b2*b2-24*b4;c6=-b2**3+36*b2*b4-216*b6
 return -c4*pow(48,-1,p)%p,-c6*pow(864,-1,p)%p,b2*pow(12,-1,p)%p

def plus(U,V,A,p):
 if U is None:return V
 if V is None:return U
 x,y=U;z,w=V
 if x==z:
  if (y+w)%p==0:return None
  lam=(3*x*x+A)*pow(2*y%p,-1,p)%p
 else:lam=(w-y)*pow((z-x)%p,-1,p)%p
 xx=(lam*lam-x-z)%p;return(xx,(lam*(x-xx)-y)%p)
def times(k,U,A,p):
 # Left-to-right, unlike producer's right-to-left multiplication.
 V=None
 for bit in bin(k)[2:]:
  V=plus(V,V,A,p)
  if bit=='1':V=plus(V,U,A,p)
 return V

def count(A,B,p):
 freq=[0]*p
 for y in range(p):freq[y*y%p]+=1
 return 1+sum(freq[(x**3+A*x+B)%p]for x in range(p))
def verify(c,z):
 ell=z['ell'];must(pr(ell));a=list(map(int,c['ainvs']));matrix=[]
 for r in z['rows']:
  p,N=r['p'],r['order'];must(pr(p)and p>3 and p!=ell and int(c['discriminant'])%p!=0);A,B,shift=short(a,p);must(count(A,B,p)==N);must(N%ell==0 and N%(ell*ell)!=0)
  def trans(xy):
   x,y=map(lambda x:coord(x,p),xy);return((x+shift)%p,(2*y+a[0]*x+a[2])*pow(2,-1,p)%p)
  G=trans(r['generator']);must(G[1]**2%p==(G[0]**3+A*G[0]+B)%p);must(times(ell,G,A,p)is None)
  for xy,k in zip(c['points'],r['row']):
   U=trans(xy);must(U[1]**2%p==(U[0]**3+A*U[0]+B)%p);must(times(N//ell,U,A,p)==times(k,G,A,p))
  must(len(r['row'])==len(c['points']));matrix.append(r['row'][:])
 rank=0
 for col in range(len(c['points'])):
  pivot=next((j for j in range(rank,len(matrix))if matrix[j][col]%ell),None)
  if pivot is None:continue
  matrix[rank],matrix[pivot]=matrix[pivot],matrix[rank];scale=pow(matrix[rank][col]%ell,-1,ell);matrix[rank]=[v*scale%ell for v in matrix[rank]]
  for j in range(len(matrix)):
   if j!=rank:
    k=matrix[j][col];matrix[j]=[(v-k*w)%ell for v,w in zip(matrix[j],matrix[rank])]
  rank+=1
 must(rank==z['rank']==len(c['points']))
 q=z['torsion_exclusion'];p=q['p'];must(pr(p)and p>3 and p!=ell and int(c['discriminant'])%p!=0);A,B,_=short(a,p);must(count(A,B,p)==q['order'] and q['order']%ell!=0)
cs={c['id']:c for c in json.loads((P/'highrank_inputs.json').read_text())['curves']};zs=json.loads((P/'odd_primary_certificates.json').read_text())
for z in zs:
 c=cs[z['curve']];verify(c,z)
 for hostile in ['duplicate','wrong_order','wrong_scalar']:
  cc=copy.deepcopy(c);zz=copy.deepcopy(z)
  if hostile=='duplicate':cc['points'][0]=cc['points'][1]
  elif hostile=='wrong_order':zz['rows'][0]['order']+=3
  else:zz['rows'][0]['row'][0]=(zz['rows'][0]['row'][0]+1)%3
  try:verify(cc,zz);accepted=True
  except ValueError:accepted=False
  must(not accepted)
print(json.dumps({'checks':checks,'full_ranks':{z['curve']:z['rank']for z in zs},'hostiles_rejected':6},sort_keys=True))
