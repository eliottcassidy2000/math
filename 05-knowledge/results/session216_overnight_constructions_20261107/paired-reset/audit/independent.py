from pathlib import Path
import itertools,json,random
P=Path(__file__).parent

def graphs(n):
 for bits in itertools.product((0,1),repeat=n*(n-1)//2):
  a=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if bits[k]:a[i]|=1<<j
    else:a[j]|=1<<i
    k+=1
  yield a

def table(a):
 n=len(a);N=1<<n;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for j in range(n):E[1<<j][j]=1
 if n:S[1][0]=1
 for A in range(1,N):
  H[A]=sum(E[A])
  for i in range(n):
   for j in range(n):
    if not A>>j&1 and a[i]>>j&1:E[A|1<<j][j]+=E[A][i];S[A|1<<j][j]+=S[A][i]
 C=[0]*N
 for A in range(N):
  B=A
  while True:
   C[A]+=H[B]*H[A^B]
   if not B:break
   B=(B-1)&A
 if not n:return H,C,None
 F=N-1;st=[sum(s) for s in S]
 p=(H[F],st[F],E[F][0],sum(E[A][0]*H[F^A] for A in range(N)),sum(st[A]*H[F^A] for A in range(N)),H[F^1],C[F],C[F^1],sum(st[A]*C[F^A] for A in range(N)))
 return H,C,p

def seven(p):return tuple(p[i] for i in (0,1,2,3,4,5,7))
def move(p,c):
 h,s,e,L,R,d,b=p
 if c=='O':return h,s,0,L,2*R,d,2*b
 if c=='P':return h,0,e,2*L,R,d,2*b
 if c=='Q':return h+L-s,d,e,2*L-d,R+b,d,2*b
 return h+R-e,s,d,L+b,2*R-d,d,2*b

def birth(a,c):
 n=len(a);a=list(a)+[0]
 for i in range(n):
  out=(c=='O' or c=='Q' and i==0 or c=='Z' and i!=0)
  if out:a[i]|=1<<n
  else:a[n]|=1<<i
 return a

def reduced(p):
 h,s,e,L,R,d,b=p
 if s!=d:raise RuntimeError('reset')
 return h+L-d+2*(d//3),L-d,b,R+b-e,d-e

def reduced_move(p,c):
 x,a,b,w,t=p
 return (x+2*a,2*a,2*b,w+2*b,t) if c=='Q' else (x+w,a+b,2*b,2*w-2*t,0)

def kappa2(p):
 h,s,e,L,R,d,b=p
 return h-L-R+b-d
transfers=reds=stable=0
for n in range(2,6):
 for a in graphs(n):
  p=seven(table(a)[2])
  for c in 'OZPQ':
   if seven(table(birth(a,c))[2])!=move(p,c):raise RuntimeError('physical move')
   transfers+=1
  q=move(p,'Q');z=move(q,'Z')
  for c in 'QZ':
   if reduced(move(q,c))!=reduced_move(reduced(q),c):raise RuntimeError('reduced')
   if kappa2(move(z,c))!=kappa2(z):raise RuntimeError('kappa invariance')
   reds+=1
  h,s,e,L,R,d,b=p
  if s==e==d:
   if kappa2(p)<0 or kappa2(p)%2:raise RuntimeError('stable nonnegative')
   stable+=1
single=seven(table([0])[2])
if kappa2(single)!=-1:raise RuntimeError('singleton hostile')
# Universal integer odd-cycle polynomial, independently evaluated against f^k.
cyclechecks=0;evenchecks=0
samples=[a for n in range(1,6) for a in graphs(n)]
rng=random.Random(216)
for z in range(32):
 a=[0]*6
 for i in range(6):
  for j in range(i+1,6):
   if rng.randrange(2):a[i]|=1<<j
   else:a[j]|=1<<i
 samples.append(a)
for a in samples:
 n=len(a);N=1<<n;F=N-1;H,C,p=table(a);cyclecoef=[0]*(n+1)
 for perm in itertools.permutations(range(n)):
  used=set();cycles=nontriv=0;ok=True
  for v in range(n):
   if v in used:continue
   cycle=[];x=v
   while x not in used:used.add(x);cycle.append(x);x=perm[x]
   cycles+=1
   if len(cycle)>1:
    nontriv+=1
    if len(cycle)%2==0 or not all(a[x]>>y&1 for x,y in zip(cycle,cycle[1:]+cycle[:1])):ok=False;break
  if ok:cyclecoef[cycles]+=2**nontriv
 power=[1]+[0]*(N-1)
 for k in range(5):
  predicted=sum(c*k**j for j,c in enumerate(cyclecoef))
  if power[F]!=predicted:raise RuntimeError('odd cycle polynomial')
  nxt=[0]*N
  for A in range(N):
   B=A
   while True:
    nxt[A]+=H[B]*power[A^B]
    if not B:break
    B=(B-1)&A
  power=nxt;cyclechecks+=1
 if n%2==0 and n>2:
  if C[F]%8:raise RuntimeError('2-adic')
  evenchecks+=1
# Independent physical ray and compiler construction; final graph counts in C++.
A=[108,1,2,6,15,30,62];B=[34,4,9,19,7,30];raw=[];ray=[]
for side,a in enumerate((A,B)):
 initial=table(a)[2]
 expect=[(33,1,1,37,55,1,660,64,589),(31,5,3,41,43,5,328,64,219)][side]
 if initial!=expect:raise RuntimeError('ray profile')
 for k in range(5):
  p=seven(table(a)[2]);x=reduced(p)[0]
  if x!=72*2**k-3:raise RuntimeError('ray value')
  n=len(a);g=list(a)+[0]*3;v=[0,n,n+1,n+2];core=[2,8,3,5]
  for i in range(4):
   for j in range(4):
    if core[i]>>j&1:g[v[i]]|=1<<v[j]
  for i in v[1:]:
   for j in range(1,n):g[i]|=1<<j
  want=216*2**k-7+2*side;raw.append((want,g));ray.append({'side':side,'k':k,'X':x,'target':want})
  a=birth(a,'Q')
(P/'RAY_GRAPHS.txt').write_text(''.join(str(h)+' '+str(len(a))+' '+' '.join(map(str,a))+'\n' for h,a in raw))
print(json.dumps({'physical_seven_transfers':transfers,'reduced_and_stable_invariance_checks':reds,'physical_stable_sources':stable,'singleton_twice_kappa':-1,'odd_cycle_polynomial_evaluations':cyclechecks,'even_div8_controls':evenchecks,'literal_compiler_graphs_prepared':len(raw),'ray':ray},indent=2))
