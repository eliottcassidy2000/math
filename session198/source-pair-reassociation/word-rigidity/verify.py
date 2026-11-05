from pathlib import Path
from itertools import permutations,product
import hashlib,json
P=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def J(s):D,E,U,d=s;return D-E+U,D-E,D+U,D
def M2(s):D,E,U,d=s;return D+2*U,D,D+2*U,D
def N(s,t):
 D,E,U,d=s;F,G,V,f=t
 return(D-E+U)*F-(D-E+d)*G+(D+3*U-d)*V-(U+E)*f
def delta(v,w,k):
 a=v;b=w
 for _ in range(k):a=J(a);b=J(b)
 return N(v,b)-N(a,w)
def constants(L):
 if L==1:return 1,1,1
 q=(L+1)//2;g,Q,B=constants(q)
 return(g+q-q%2,Q+B,2*B)if L%2 else(g,Q,2*B)
def step(s,L):
 D,E,U,d=s;g,Q,B=constants(L)
 return g*D-(L%2)*E+L*U,D-(L==1)*E,Q*D+B*U,D
basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
for v,w in product(basis,repeat=2):
 D,E,U,d=v;F,G,V,f=w;p=D+d;r=D-E+U-d
 a,b,c=[delta(v,w,k)for k in(1,2,3)]
 need(a==r*(F+f)-p*(F-G+V-f),'first defect')
 need(b-a==2*(U*(F+f)-p*V),'second difference')
 need(c-2*b+a==2*(D*f-d*F),'third difference')
for v in basis:
 js=[v]
 for k in range(4):js.append(J(js[-1]))
 need(js[4]==tuple(js[1][i]+js[2][i]+js[3][i]for i in range(4)),'matrix recurrence')
def det3(A):
 return sum((-1)**sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))*A[0][p[0]]*A[1][p[1]]*A[2][p[2]]for p in permutations(range(3)))
minor_checks=0
for v in product(range(-2,3),repeat=4):
 D,E,U,d=v;rows=[[delta(v,w,k)for w in basis]for k in(1,2,3)]
 signed=[(-1)**j*det3([[r[i]for i in range(4)if i!=j]for r in rows])for j in range(4)]
 need(signed==[-4*(D+d)**2*x for x in v],('signed minors',v))
 minor_checks+=1
cap=101;bank={(1,1,1,1)};queue=list(bank)
for v in queue:
 D,E,U,d=v
 for w in[(D,0,2*U,d)]+[step(v,L)for L in range(1,cap//U+1)]:
  if w[2]<=cap and w not in bank:bank.add(w);queue.append(w)
ss=sorted(bank-{(1,1,1,1)});pairs=0;commuting=0
for i,v in enumerate(ss):
 for w in ss[i:]:
  good=all(delta(v,w,k)==0 for k in(1,2,3))
  proportional=all(v[a]*w[b]==v[b]*w[a]for a in range(4)for b in range(4))
  need(good==proportional,('three test',v,w))
  if good:
   commuting+=1
   need(all(delta(v,w,k)==0 for k in range(4,11)),'longer returning test')
  pairs+=1
t=(1,0,2,1);z=(1,1,1,1);ray=[];trib=[1,1,2]
for k in range(3,20):trib.append(sum(trib[-3:]))
for j in range(8):
 w=M2(z);d=z[0]
 need(w==(3*d,d,3*d,d),'physical balanced ray')
 values=[delta(t,w,k)for k in range(22)]
 need(values[:2]==[0,0]and values[2:]==[4*d*x for x in trib],'tribonacci defect')
 ray.append(dict(d=d,defects=values));z=M2(z)

def birth(T,r):
 n=len(T);A=[x|(1<<n)for x in T]+[1<<r];A[r]^=1<<n
 return A,n

def graft(V,c,W,r):
 T,b=birth(V,c);n=len(T);G=[((1<<len(W))-1)^(1<<i)^x for i,x in enumerate(W)];A=T+[x<<n for x in G]
 for x in range(n):
  for y in range(len(G)):
   if y==r and x in(c,b):A[n+y]|=1<<x
   else:A[x]|=1<<(n+y)
 return A

def HC(T):
 n=len(T);h=sum(all(T[p[i]]>>p[i+1]&1 for i in range(n-1))for p in permutations(range(n)))
 c=sum(all(T[p[i]]>>p[(i+1)%n]&1 for i in range(n))for q in permutations(range(1,n))for p in[(0,)+q])
 return h,c
V=[2,0];c=0;W,r=birth(V,c);JV,jc=birth(*birth(V,c));JW,jr=birth(*birth(W,r))
hostiles=[]
for args,expect in[((V,c,JW,jr),77),((JV,jc,W,r),73)]:
 A=graft(*args);h,cyc=HC(A);need((h,cyc)==(expect,2),'literal word hostile')
 hostiles.append(dict(H=h,C=cyc,order=len(A),rows=A))
receipt=dict(status='PASS',basis_pair_checks=16,signed_minor_checks=minor_checks,source_cap=cap,nontrivial_source_profiles=len(ss),unordered_pairs=pairs,commuting_pairs=commuting,physical_balanced_rays=ray,literal_hostiles=hostiles,dependency_sha256=hashlib.sha256((P/'dependencies/BILINEAR_C2.md').read_bytes()).hexdigest())
s=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(P/'VALIDATION.json').write_text(s);print(s,end='')
