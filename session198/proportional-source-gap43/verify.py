from pathlib import Path
from itertools import permutations,product
import hashlib,json
P=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def J(s):D,E,U,d=s;return D-E+U,D-E,D+U,D
def N(s,t):
 D,E,U,d=s;F,G,V,f=t
 return(D-E+U)*F-(D-E+d)*G+(D+3*U-d)*V-(U+E)*f
def M(s,t):
 D,E,U,d=s;F,G,V,f=t
 return(D+2*U)*F-(U+d)*G+(2*D-E+6*U-d)*V-(U+E)*f
def delta(v,w,k,F=M):
 a=v;b=w
 for _ in range(k):a=J(a);b=J(b)
 return F(v,b)-F(a,w)
def Q(v):
 D,E,U,d=v
 return 2*D*D-D*E-D*U+3*D*d+E*U-E*d+U*U-U*d+d*d
def predicted(v):
 D,E,U,d=v
 return[[-E+3*U-2*d,D+d,-3*D-d,2*D-E+U],[-E+7*U-2*d,D-2*U+d,-7*D+2*E-3*d,2*D-E+3*U],[-E+11*U-4*d,D-2*U+d,-11*D+2*E-5*d,4*D-E+5*U]]
basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
for v in basis:need(predicted(v)==[[delta(v,w,k)for w in basis]for k in(1,2,3)],'M defect coefficient matrix')
def det3(A):
 return sum((-1)**sum(p[i]>p[j]for i in range(3)for j in range(i+1,3))*A[0][p[0]]*A[1][p[1]]*A[2][p[2]]for p in permutations(range(3)))
for v in product(range(4),repeat=4):
 R=predicted(v);minors=[(-1)**j*det3([[r[i]for i in range(4)if i!=j]for r in R])for j in range(4)]
 need(minors==[-4*Q(v)*x for x in v],'cubic maximal minors')
# Degree at most3 in each coordinate:4^4 evaluations certify each identity.
for a,f,e,p in product(range(3),repeat=4):
 v=(a+2*e+p,e,a+f+2*e+2*p,e+p)
 q=2*a*a+10*a*e+a*f+7*a*p+12*e*e+2*e*f+17*e*p+f*f+2*f*p+6*p*p
 need(Q(v)==q and q>=6*(e+p)**2,'quadratic positive expansion')
# Degree at most2 in each coordinate:3^4 evaluations certify the expansion.
def constants(L):
 if L==1:return 1,1,1
 ell=(L+1)//2;g,q,b=constants(ell)
 return(g+ell-ell%2,q+b,2*b)if L%2 else(g,q,2*b)
def step(s,L):
 D,E,U,d=s;g,q,b=constants(L)
 return g*D-L%2*E+L*U,D-(L==1)*E,q*D+b*U,D
cap=101;bank={(1,1,1,1)};queue=list(bank)
for s in queue:
 D,E,U,d=s
 for t in[(D,0,2*U,d)]+[step(s,L)for L in range(1,cap//U+1)]:
  if t[2]<=cap and t not in bank:bank.add(t);queue.append(t)
ss=sorted(bank-{(1,1,1,1)});pairs=0;eq=0
for i,v in enumerate(ss):
 for w in ss[i:]:
  proportional=all(v[a]*w[b]==v[b]*w[a]for a in range(4)for b in range(4))
  tests=all(delta(v,w,k)==0 for k in(1,2,3))
  need(tests==proportional,'M three-test rigidity');pairs+=1;eq+=tests
  if tests:need(all(delta(v,w,k)==0 for k in range(4,9)),'M longer word return')
c1={s[0]for s in bank if s[1]>=1 and 3<=s[0]<=43 and s[2]<=86}
missing=sorted(set(range(3,44,2))-c1)
need(missing==[7,13,21,41,43],'complete C1 small image')
t=(1,0,2,1);q=(3,1,3,1);w=(3,0,6,1)
need([s for s in ss if s[2]<4]==[t,q],'small physical source inventory')
need((N(t,t),N(q,q),M(t,t))==(13,41,29),'small diagonal values')
need(N(t,w)==43,'physical profile43')
defects=[delta(t,w,k,N)for k in(1,2,3)]
need(defects==[-8,-16,-28],('43defects',defects))
def birth(T,r):
 n=len(T);A=[x|(1<<n)for x in T]+[1<<r];A[r]^=1<<n
 return A,n
V=[2,0];c=0;W,r=birth(V,c);n=len(W);W=[x|(1<<n)for x in W]+[0]
T,b=birth(V,c);n=len(T);G=[((1<<len(W))-1)^(1<<i)^x for i,x in enumerate(W)];A=T+[x<<n for x in G]
for x in range(n):
 for y in range(len(G)):
  if y==r and x in(c,b):A[n+y]|=1<<x
  else:A[x]|=1<<(n+y)
size=len(A)
h=sum(all(A[p[i]]>>p[i+1]&1 for i in range(size-1))for p in permutations(range(size)))
cyc=sum(all(A[p[i]]>>p[(i+1)%size]&1 for i in range(size))for q in permutations(range(1,size))for p in[(0,)+q])
need((h,cyc)==(43,2),'literal43')
receipt=dict(status='PASS',coefficient_pair_checks=16,cubic_unisolvent_grid=256,quadratic_unisolvent_grid=81,source_cap=cap,nontrivial_sources=len(ss),unordered_pairs=pairs,commuting_pairs=eq,c1_control=dict(target_cap=43,complete_source_cap=86,missing_strong_C1=missing),diagonal_small_counts=dict(N=[13,41],M=[29]),literal43=dict(rows=A,H=h,C=cyc,order=size,first_three_N_shift_defects=defects),dependencies={f.name:hashlib.sha256(f.read_bytes()).hexdigest()for f in sorted((P/'dependencies').glob('*.md'))})
s=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(P/'VALIDATION.json').write_text(s);print(s,end='')
