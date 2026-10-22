"""Independent linear algebra controls, not an asymptotic accuracy certificate."""
import itertools,json,math
import numpy as np

def mat(n,mask):
 s=np.zeros((n,n),dtype=np.int64)
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  s[i,j]=1 if mask>>k&1 else -1;s[j,i]=-s[i,j]
 return s

def rho(m):
 x=np.linalg.eigvalsh(1j*m)
 return math.exp(float(.5*np.sum(np.log1p(x*x)-np.log1p(-x*x))))

def test(s):
 n=len(s);one=np.ones(n);d=s@one;m=s/n
 p=np.eye(n)-np.ones((n,n))/n;c=p@m@p
 a=d/(n*math.sqrt(n));i=np.eye(n)
 assert np.max(np.abs(c@one))<1e-12
 assert np.max(np.abs(m-c-(np.outer(d,one)-np.outer(one,d))/n**2))<1e-12
 plus=float(a@np.linalg.solve(i-c@c,a))
 minus=float(a@np.linalg.solve(i+c@c,a))
 assert math.isclose(rho(m)/rho(c),(1+plus)/(1-minus),rel_tol=2e-12,abs_tol=2e-12)
 bound=(n-1)/(2*n)
 assert np.linalg.norm(m,2)**2<=bound+1e-12
 if np.all(np.abs(d)==1):
  sp=np.block([[s,-d[:,None]],[d[None,:],np.zeros((1,1))]])
  assert np.array_equal(sp@np.ones(n+1),np.zeros(n+1))
  t=float((d@np.linalg.solve(i+m@m,d))/n)
  assert 1-1e-12<=t<=2*n/(n+1)+1e-12
  # rank-two centering alters rho only at second order, uniformly.
  assert 1<=rho(m)/rho(c)+1e-12
  upper=(1+1/n**2)/(1-2*n/(n+1)/n**2)
  assert rho(m)/rho(c)<=upper+1e-12
  return 1
 return 0

out={'scope':'all labelled tournaments n2..5, all almost-regular labelled n6, selected deleted Paley','counts':{},'paley':[]}
for n in range(2,6):
 ar=0
 for mask in range(1<<(n*(n-1)//2)):ar+=test(mat(n,mask))
 out['counts'][str(n)]={'all':1<<(n*(n-1)//2),'almost_regular':ar}
n=6;ar=0
for mask in range(1<<15):
 s=mat(n,mask)
 if np.all(np.abs(s.sum(axis=1))==1):ar+=test(s)
out['counts']['6']={'almost_regular':ar}
for q in [3,7,11,19,23,31,43,47,59,67,71,79,83]:
 residues={x*x%q for x in range(1,q)}
 sp=np.array([[0 if x==y else 1 if (y-x)%q in residues else -1 for y in range(q)] for x in range(q)],dtype=np.int64)
 s=sp[:-1,:-1];n=q-1;d=s.sum(axis=1);j=np.ones((n,n),dtype=np.int64)
 assert np.array_equal(s@s,j-q*np.eye(n,dtype=np.int64)+np.outer(d,d))
 assert np.array_equal(s@d,-np.ones(n,dtype=np.int64))
 test(s)
 predicted=((n*n+1)/(n*n-1))*((n*n+n+1)/(n*n-n-1))**((n-2)/2)
 assert math.isclose(rho(s/n),predicted,rel_tol=2e-12)
 t=float(d@np.linalg.solve(np.eye(n)+(s/n)@(s/n),d)/n)
 assert math.isclose(t,n*n/(n*n-1),rel_tol=2e-12)
 out['paley'].append({'q':q,'rho':predicted,'tree_penalty':t})
out['status']='PASS'
print(json.dumps(out,indent=2))
