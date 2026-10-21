import sys,json,itertools
from collections import defaultdict
sys.path.insert(0,'session24-pathcovers')
from verify_pathcovers import subset_c,brute_fixed

def mul(a,b,caps):
 d=defaultdict(int)
 for x,v in a.items():
  for y,w in b.items():
   z=tuple(i+j for i,j in zip(x,y))
   if all(i<=j for i,j in zip(z,caps)):d[z]+=v*w
 return d

def walks(A,typ,caps,endpoint=None):
 d=defaultdict(int); zero=(0,)*len(caps)
 if endpoint is not None:d[zero]=1
 def rec(last,x):
  if endpoint is None or A[last][endpoint]:d[x]+=1
  for j in range(len(A)):
   if A[last][j] and x[typ[j]]<caps[typ[j]]:
    z=list(x);z[typ[j]]+=1;rec(j,tuple(z))
 for j in range(len(A)):
  if caps[typ[j]]:
   z=list(zero);z[typ[j]]=1;rec(j,tuple(z))
 return d

def blow(A,B):
 vertices=[(i,j) for i,b in enumerate(B) for j in range(len(b))]
 return [[A[i][k] if i!=k else B[i][j][l] for k,l in vertices] for i,j in vertices]
T=lambda n:[[int(i<j) for j in range(n)] for i in range(n)]
C=[[int((j-i)%3==1) for j in range(3)] for i in range(3)]
rows=[]
for A,B in [(C,[T(1),T(2),C]),(T(3),[C,C,C]),(T(2),[C,T(2)])]:
 caps=[len(b) for b in B];cv=[subset_c(b) for b in B];W=walks(A,list(range(len(A))),caps);P={(0,)*len(A):1};out=[0]
 for k in range(1,sum(caps)+1):
  P=mul(P,W,caps);out.append(sum(v*__import__('math').prod(cv[i][r] for i,r in enumerate(x)) for x,v in P.items()))
 direct=subset_c(blow(A,B))
 if out!=direct:raise RuntimeError(('ordinary mismatch',out,direct))
 rows.append(dict(kind='ordinary unequal or off-family',sizes=caps,cover=out))
A=T(3);B=C;caps=[1,3];typ=[1,0,1];W=walks(A,typ,caps);V=walks(A,typ,caps,1);cv=subset_c(B);dv=brute_fixed(B,[0,2,1]);out=[0]*10
for l in range(5):
 out[2*l+1]=sum(v*dv[2*s+1]*cv[r] for (s,r),v in V.items());V=mul(V,W,caps)
AA=blow(A,[B]*3);tau=[(2-i)*3+(-j)%3 for i in range(3) for j in range(3)];direct=brute_fixed(AA,tau)
if out!=direct:raise RuntimeError(('fixed mismatch',out,direct))
rows.append(dict(kind='full fixed vector off-family transitive outer',d=out))
print(json.dumps(rows,indent=2))
