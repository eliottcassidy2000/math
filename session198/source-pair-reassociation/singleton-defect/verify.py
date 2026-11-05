from pathlib import Path
from collections import deque
from itertools import combinations,permutations
import json,hashlib
P=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def constants(L):
 if L==1:return 1,1,1
 a=(L+1)//2;g,q,b=constants(a)
 return (g+a-a%2,q+b,2*b)if L%2 else(g,q,2*b)

def step(v,L):
 D,E,U,d=v;g,q,b=constants(L)
 return g*D-(L%2)*E+L*U,D-(L==1)*E,q*D+b*U,D

def sink(v):D,E,U,d=v;return D,0,2*U,d

def N(v,w):
 D,E,U,d=v;F,G,V,f=w
 return (D-E+U)*F-(D-E+d)*G+(D+3*U-d)*V-(U+E)*f

def M(v,w):
 D,E,U,d=v;F,G,V,f=w
 return (D+2*U)*F-(U+d)*G+(2*D-E+6*U-d)*V-(U+E)*f

def p(v):return v[0]+v[3]
def rho(v):return v[0]-v[1]+v[2]-v[3]

basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
skews={}
for name,F in [('N',N),('M',M)]:
 A=[[F(v,step(w,1))-F(step(v,1),w)for w in basis]for v in basis]
 for i,v in enumerate(basis):
  for j,w in enumerate(basis):
   expected=rho(v)*p(w)-p(v)*rho(w)+(2*(v[2]*w[0]-v[0]*w[2])if name=='M'else 0)
   need(A[i][j]==expected,'basis defect identity')
 skews[name]=A
expected_N=[[0,1,-1,2],[-1,0,0,-1],[1,0,0,1],[-2,1,-1,0]]
need(skews['N']==expected_N,'N skew matrix')
def determinant(A):
 total=0
 for perm in permutations(range(4)):
  inversions=sum(perm[i]>perm[j]for i in range(4)for j in range(i+1,4))
  term=(-1)**inversions
  for i in range(4):term*=A[i][perm[i]]
  total+=term
 return total
need(determinant(skews['M'])==4,'M rank four')
need(determinant(skews['N'])==0 and skews['N'][0][1]==1,'N nonzero alternating wedge')

cap=501;bank={(1,1,1,1):()};queue=list(bank)
for v in queue:
 for label,w in[(0,sink(v))]+[(L,step(v,L))for L in range(1,cap//v[2]+1)]:
  if w[2]<=cap and w not in bank:bank[w]=bank[v]+(label,);queue.append(w)
t=(1,0,2,1);fiber=[]
for v in bank:
 if v==(1,1,1,1):continue
 D,E,U,d=v
 actual=U-E==2*d
 predicted=v==t or(v==(3*d,d,3*d,d)and any(z[0]==z[2]==d for z in bank))
 need(actual==predicted,('fiber',v,predicted))
 need(N(t,step(v,1))==N(step(t,1),v)if actual else N(t,step(v,1))!=N(step(t,1),v),'fiber equality')
 if actual:fiber.append(list(v))
small=[v for v in bank if v[2]<=41 and v!=(1,1,1,1)];pairs=0
for v in small:
 for w in small:
  delta=rho(v)*p(w)-p(v)*rho(w)
  need(N(v,step(w,1))-N(step(v,1),w)==delta,('N identity',v,w))
  need(M(v,step(w,1))-M(step(v,1),w)==delta+2*(v[2]*w[0]-v[0]*w[2]),('M identity',v,w))
  pairs+=1

def stats(T,root=0):
 n=len(T);F=(1<<n)-1;st=[[0]*n for _ in range(F+1)];en=[[0]*n for _ in range(F+1)];H=[1]+[0]*F
 for S in range(1,F+1):
  for x in range(n):
   if not S>>x&1:continue
   R=S^(1<<x)
   if not R:st[S][x]=en[S][x]=1
   else:
    st[S][x]=sum(st[R][y]for y in range(n)if R>>y&1 and T[x]>>y&1)
    en[S][x]=sum(en[R][y]for y in range(n)if R>>y&1 and T[y]>>x&1)
  H[S]=sum(en[S])
 prof=H[F],en[F][root],sum(H[S]*st[F^S][root]for S in range(F+1)),H[F^(1<<root)]
 cyc=sum(all(T[a]>>b&1 for a,b in zip((0,)+q, q+(0,)))for q in permutations(range(1,n)))if n>=3 else 0
 return prof,st[F][root],cyc

def transitive(n):return[sum(1<<j for j in range(i+1,n))for i in range(n)],0

def J(T,r):
 n=len(T);A=[row|(1<<n)for row in T]+[1<<r];A[r]^=1<<n
 return A,n

def construct(kind,V,c,W,r):
 T,b=J(V,c)
 if kind=='N':entries={b,c}
 else:T,u=J(T,b);entries={u,c}
 n=len(T);G=[((1<<len(W))-1)^(1<<i)^row for i,row in enumerate(W)]
 A=T+[row<<n for row in G]
 for x in range(n):
  for y in range(len(G)):
   if y==r and x in entries:A[n+y]|=1<<x
   else:A[x]|=1<<(n+y)
 return A

rooted={}
for n in(2,3):
 S=set();edges=list(combinations(range(n),2))
 for code in range(1<<len(edges)):
  T=[0]*n
  for k,(i,j)in enumerate(edges):T[i if code>>k&1 else j]|=1<<(j if code>>k&1 else i)
  for r in range(n):
   v,s,c=stats(T,r)
   if s==1:S.add(v)
 rooted[n]=sorted(S)
need(rooted=={2:[t],3:[(1,0,4,1),(3,1,3,1)]},('minimal source universe',rooted))
V,c=transitive(2);W,r=transitive(3);JV,jc=J(V,c);JW,jr=J(W,r);hostile=[]
for kind,expected in[('N',(41,45)),('M',(85,93))]:
 for side,args,h in[('right',(V,c,JW,jr),expected[0]),('left',(JV,jc,W,r),expected[1])]:
  T=construct(kind,*args);v,s,cyc=stats(T)
  need(v[0]==h and cyc==2,('physical hostile',kind,side,v,cyc))
  hostile.append(dict(kind=kind,birth_on=side,H=h,C=cyc,order=len(T),rows=T))
receipt=dict(status='PASS',defect_matrices=skews,source_cap=cap,source_states=len(bank),fiber=fiber,algebra_source_cap=41,ordered_algebra_pairs=pairs,rooted_through_three=rooted,hostile_graphs=hostile,dependencies={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((P/'dependencies').glob('*.md'))})
s=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(P/'VALIDATION.json').write_text(s);print(s,end='')
