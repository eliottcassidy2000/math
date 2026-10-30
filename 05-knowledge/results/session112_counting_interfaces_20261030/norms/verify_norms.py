"""Exact polynomial coefficient tests, not random norm sampling. Python stdlib."""
from itertools import product
from collections import Counter
import json

def require(x,msg):
 if not x: raise RuntimeError(msg)
def lines(q):
 return sorted({tuple(sorted((a,b,a^b))) for a in range(1,q) for b in range(1,q) if a!=b})
def table(q,z):
 f=[[0]*q for _ in range(q)]
 for a in range(1,q):f[a][a]=1
 for (u,v,w),bit in zip(lines(q),z):
  for a,b in ((u,v),(v,w),(w,u)):f[a][b]=bit;f[b][a]=bit^1
 return f

def rectangle(f,t,a,b):
 return f[a][b]^f[a^t][b^t]^f[a][b^t]^f[a^t][b]
def valid(f):
 q=len(f)
 return all(rectangle(f,t,a,b)==1 for t in range(1,q) for a in range(q) for b in range(q))
def polynomial_norm(f):
 q=len(f);coeff=Counter()
 for a,b,c,d in product(range(q),repeat=4):
  if (a^b)==(c^d):
   coeff[tuple(sorted((a,c)))+tuple(sorted((b,d)))]+=(-1)**(f[a][b]+f[c][d])
 for a,b in product(range(q),repeat=2):coeff[(a,a,b,b)]-=1
 return all(c==0 for c in coeff.values())

def line_equation(q,t,a,b):
 ls=lines(q); ix={x:i for i,x in enumerate(ls)};support=set();rhs=1
 for u,v in ((a,b),(a^t,b^t),(a,b^t),(a^t,b)):
  if not u or not v:continue
  if u==v:rhs^=1;continue
  l=tuple(sorted((u,v,u^v)));x,y,z=l
  rhs^=int((u,v) not in ((x,y),(y,z),(z,x)))
  if l in support:support.remove(l)
  else:support.add(l)
 return support,rhs

def main():
 counts={};checked=0
 for q in (2,4,8):
  good=[]
  for z in product(range(2),repeat=len(lines(q))):
   f=table(q,z);v=valid(f);poly=polynomial_norm(f)
   require(v==poly,'rectangle/polynomial disagreement');checked+=1
   if v:good.append(z)
  counts[str(q)]=len(good)
  if q==8:
   formula={(A,B,C,D,1^B^C^D,1^A^C^D,1^A^B^D) for A,B,C,D in product(range(2),repeat=4)}
   require(set(good)==formula,'four-probe classification')
   base=table(8,(0,0,0,0,1,1,1));gauges=set()
   for bits in product(range(2),repeat=4):
    s=[0]*8
    for j,b in zip((3,5,6,7),bits):s[j]=b
    f=[[base[a][b]^s[a]^s[b]^s[a^b] for b in range(8)] for a in range(8)]
    require(valid(f),'gauge invalid')
    gauges.add(tuple(f[u][v] for u,v,w in lines(8)))
   require(gauges==formula,'gauge code incomplete')
 # Exact prime-selection transport to the old sigma Hamming bank.
 primes=(19,2,59,29,5,89,7)
 from math import prod,isqrt
 for z in product(range(2),repeat=7):
  sigma=prod(p+1 for p,b in zip(primes,z) if b)
  square30=sigma%30==0 and isqrt(sigma//30)**2==sigma//30
  require(square30==valid(table(8,z)),'sigma/norm affine-code transport')
 # Five necessary rectangle identities in dimension 16 have XOR 0=1.
 cert=[(2,5,4),(8,1,7),(11,4,1),(12,6,5),(13,4,3)]
 parity=Counter();rhs=0;detail=[]
 for t,a,b in cert:
  support,r=line_equation(16,t,a,b);rhs^=r;parity.update(support)
  detail.append({'t':t,'a':a,'b':b,'line_variables':[list(x) for x in sorted(support)],'rhs':r})
 require(rhs==1 and all(n%2==0 for n in parity.values()),'invalid contradiction')
 # Direct expansion independently checks the reduced line formula for every assignment
 # of the ten line variables appearing in the certificate, not all 35.
 used=sorted(parity); ls=lines(16)
 for bits in product(range(2),repeat=len(used)):
  z=[0]*len(ls)
  for l,b in zip(used,bits):z[ls.index(l)]=b
  f=table(16,z)
  require(not all(rectangle(f,*c)==1 for c in cert),'certificate has solution')
 require(counts=={'2':1,'4':2,'8':16},'wrong classification')
 # Norm defect gives a two-term sparse witness whenever a rectangle fails.
 f=table(8,(0,)*7);bad=next((t,a,b) for t in range(1,8) for a in range(8) for b in range(8) if not rectangle(f,t,a,b))
 t,a,b=bad;x=[0]*8;y=[0]*8;x[a]=x[a^t]=1;y[b]=y[b^t]=1
 xy=[sum((-1)**f[u][v]*x[u]*y[v] for u in range(8) for v in range(8) if u^v==j) for j in range(8)]
 require(sum(v*v for v in xy)!=4,'no sparse defect')
 print(json.dumps({'status':'PASS','sigma_selection_checks':128,'norm_tables':counts,'line_tables_polynomial_checked':checked,'four_probes':['(1,2)','(1,4)','(1,6)','(2,4)'],'sixteen_square_contradiction':detail,'certificate_variables':len(used),'certificate_assignments':2**len(used),'sparse_defect':{'rectangle':bad,'x':x,'y':y,'xy':xy,'norm_product':4,'norm_xy':sum(v*v for v in xy)}},indent=2,sort_keys=True))
if __name__=='__main__':main()
