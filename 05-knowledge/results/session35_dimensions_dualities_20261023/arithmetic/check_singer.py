#!/usr/bin/env python3
import itertools,json,math
from pathlib import Path
from collections import Counter

def require(b,msg):
 if not b:raise RuntimeError(msg)
def pf(n):
 out=[];d=2
 while d*d<=n:
  if n%d==0:
   out.append(d)
   while n%d==0:n//=d
  d+=1
 if n>1:out.append(n)
 return out

def plane(q):
 # Search monic irreducible cubic: a cubic has no proper factor iff no root.
 c=next(c for c in itertools.product(range(q),repeat=3) if c[0] and all((a**3+c[2]*a*a+c[1]*a+c[0])%q for a in range(q)))
 def mul(a,b):
  t=[0]*5
  for i in range(3):
   for j in range(3):t[i+j]=(t[i+j]+a[i]*b[j])%q
  for k in (4,3):
   for j in range(3):t[k-3+j]=(t[k-3+j]-t[k]*c[j])%q
  return tuple(t[:3])
 one=(1,0,0);zero=(0,0,0)
 def power(a,n):
  r=one
  while n:
   if n&1:r=mul(r,a)
   a=mul(a,a);n//=2
  return r
 order=q**3-1
 alpha=next(a for a in itertools.product(range(q),repeat=3) if a!=zero and all(power(a,order//r)!=one for r in pf(order)))
 v=q*q+q+1;D=[]
 for i in range(v):
  a=power(alpha,i);b=power(a,q);cc=power(b,q)
  tr=tuple((a[k]+b[k]+cc[k])%q for k in range(3))
  require(tr[1:]==(0,0),'trace lands in base field')
  if tr==zero:D.append(i)
 require(len(D)==q+1,'line size')
 dif=Counter((a-b)%v for a in D for b in D if a!=b)
 require(dif==Counter({i:1 for i in range(1,v)}),'Singer difference uniqueness')
 lines=[sorted((x+i)%v for x in D) for i in range(v)]
 owners={}
 for i,L in enumerate(lines):
  for u,w in itertools.combinations(L,2):
   require((u,w) not in owners,'unique line edge');owners[u,w]=i
 require(len(owners)==math.comb(v,2),'partition K_v')
 # B rows vertices, columns lines; exact binary rank.
 B=[sum(1<<i for i,L in enumerate(lines) if u in L) for u in range(v)]
 piv={}
 for row in B:
  r=row
  while r:
   k=r.bit_length()-1
   if k in piv:r^=piv[k]
   else:piv[k]=r;break
 rank=len(piv)
 if q%2:
  require(rank==v-1,'odd incidence rank')
  require(all((r.bit_count()%2)==0 for r in B),'all ones kernel')
  for u in range(v):
   for w in range(v):require(((B[u]&B[w]).bit_count()%2)==int(u!=w),'BBt=I+J')
 # Circulant regular baseline, with Fano Paley gauge for q2.
 half={1,2,4} if q==2 else set(range(1,(v+1)//2))
 S=[[0 if u==w else (1 if (w-u)%v in half else -1) for w in range(v)] for u in range(v)]
 h=[[sum(S[u][w] for w in L if w!=u) if u in L else 0 for L in lines] for u in range(v)]
 intersections={(i,j):next(iter(set(lines[i])&set(lines[j]))) for i in range(v) for j in range(i+1,v)}
 # Four times each quadratic Walsh coefficient, integral.
 coeff4={(i,j):-h[w][i]*h[w][j] for (i,j),w in intersections.items()}
 if q%2:require(all(c%2 for c in coeff4.values()),'all quadratic Walsh pairs nonzero')
 if q==2:require(all(x==0 for row in h for x in row),'Fano local cycles regular')
 # Eight times constant term C(n,3)-sum E binom(d_v,2).
 d0=(v-1)//2
 constant8=8*(math.comb(v,3)-v*d0*(d0-1)//2)-sum(x*x for row in h for x in row)
 masks=range(1<<v) if q in (2,3) else [0,(1<<v)-1]+[1<<i for i in range(v)]+[sum(1<<i for i in range(v) if (17*i+11)%23<10)]
 regular=0;triangle_values=[]
 for m in masks:
  z=[-1 if m>>i&1 else 1 for i in range(v)]
  deg=[d0+sum(h[u][i]*z[i] for i in range(v))//2 for u in range(v)]
  # Sum is even since degree integer.
  require(all(sum(h[u][i]*z[i] for i in range(v))%2==0 for u in range(v)),'integer degrees')
  regular+=int(all(d==d0 for d in deg))
  tri=math.comb(v,3)-sum(d*(d-1)//2 for d in deg)
  expansion8=constant8+2*sum(c*z[i]*z[j] for (i,j),c in coeff4.items())
  require(8*tri==expansion8,'triangle Fourier identity')
  if q%2:require(all(((deg[u]-d0)%2)==((B[u]&m).bit_count()%2) for u in range(v)),'score parity syndrome')
  triangle_values.append(tri)
 if q==3:
  require(regular==2,'only converse regular')
  # Independent Walsh transform on all8192 observed values.
  a=triangle_values[:];step=1
  while step<len(a):
   for start in range(0,len(a),2*step):
    for j in range(start,start+step):a[j],a[j+step]=a[j]+a[j+step],a[j]-a[j+step]
   step*=2
  support=[m for m,c0 in enumerate(a) if c0]
  require(len(support)==1+math.comb(v,2),'full exact Walsh support')
  require(all(m==0 or m.bit_count()==2 for m in support),'Walsh degree2')
  # Verify exactly0,all-ones translation stabilizers from independently computed support.
  stabilizers=[m for m in range(1<<v) if all((m&s).bit_count()%2==0 for s in support)]
  require(stabilizers==[0,(1<<v)-1],'continuation stabilizer')
 out={'q':q,'v':v,'polynomial_constant_first':list(c)+[1],'primitive_alpha':alpha,'difference_set':D,'binary_incidence_rank':rank,'edge_blocks':v,'edges_per_block':math.comb(q+1,2),'masks_tested':len(triangle_values),'regular_masks_tested':regular,'triangle_range_tested':[min(triangle_values),max(triangle_values)],'triangle_mean_times8':constant8,'nonzero_quadratic_coefficients':sum(c!=0 for c in coeff4.values())}
 if q==7:
  # Explicit CRT additive transport, and structural collision with ring multiplication.
  require(len({(i%3,i%19) for i in range(v)})==v,'CRT bijective')
  out['difference_set_CRT']=[[i%3,i%19] for i in D]
 return out
out={'status':'FINITE-EXACT','planes':[plane(q) for q in (2,3,7)]}
Path(__file__).with_name('singer_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
