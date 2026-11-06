import itertools,json,collections,random

def need(q):
 if not q: raise RuntimeError('independent check failed')
def graphs(n):
 pairs=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  a=[0]*n
  for k,(i,j) in enumerate(pairs):
   u,v=(i,j) if mask>>k&1 else (j,i);a[u]|=1<<v
  yield tuple(a)
def profile(a):
 n=len(a);full=(1<<n)-1
 # DP ending at v on each subset, independently using incoming edges.
 dp=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    dp[s][v]=1 if not t else sum(dp[t][u] for u in range(n) if a[u]>>v&1)
  H[s]=sum(dp[s])
 # Start-root path counts use reversed graph's end-root DP.
 rev=tuple(sum(1<<u for u in range(n) if a[u]>>v&1) for v in range(n))
 D=[[0]*n for _ in range(1<<n)]
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v);D[s][v]=1 if not t else sum(D[t][u] for u in range(n) if rev[u]>>v&1)
 R=sum(D[s][0]*H[full^s] for s in range(1,full+1,2))
 L=sum(dp[s][0]*H[full^s] for s in range(1,full+1,2))
 d=H[full^1];beta=sum(H[s]*H[(full^1)^s] for s in range(0,full+1,2))
 return (H[full],R-d,d-dp[full][0],d),(d,beta,L-d,R-d,d-D[full][0],d-dp[full][0],H[full]-d)
def join(a,b):
 n=len(a);r=list(a)+[0]*(len(b)-1)
 for i in range(1,n):r[i]|=((1<<(len(b)-1))-1)<<n
 def idx(i):return 0 if i==0 else n+i-1
 for i in range(len(b)):
  for j in range(len(b)):
   if b[i]>>j&1:r[idx(i)]|=1<<idx(j)
 return tuple(r)
def transfer(s,k):
 h,x,y,d=s;H,B,L,A,S,E,K=k
 return H*h+L*x+S*y+K*d,B*x+A*d,E*d,H*d
small=[a for n in range(1,4) for a in graphs(n)]
cores=[a for n in range(2,5) for a in graphs(n)]
checks=0
for a in small:
 s,_=profile(a)
 for b in cores:
  _,k=profile(b);need(min(k)>=0)
  need(profile(join(a,b))[0]==transfer(s,k));checks+=1
# Physical C3 core constants and all j incidence patterns.
ks=[]
for j in range(4):
 b=[0]*4
 for v in range(3):
  b[v+1]|=1<<((v+1)%3+1)
  if v<j:b[0]|=1<<(v+1)
  else:b[v+1]|=1
 s,k=profile(b);need(k==(3,12,9-3*j,3*j,3-j,j,(3,5,5,3)[j]-3));ks.append(k)
# General integer divisibility controls, including composite c and variable coefficients.
rng=random.Random(7402);grades=0
for c in range(2,13):
 for trial in range(100):
  s=(1,0,0,1);m=0
  for step in range(40):
   h,x,y,d=s;op=rng.randrange(3)
   if op==0:s=(h,2*x+d,d,d)
   elif op==1:s=(h+x+y,2*x,0,d)
   else:
    k=(c*rng.randrange(1,5),c*rng.randrange(1,5),*(rng.randrange(5) for _ in range(5)))
    s=transfer(s,k);m+=1
   h,x,y,d=s;need(d%c**m==0)
   if m:need(x%c**(m-1)==0 and y%c**(m-1)==0)
   if m>=2:need(h%c**(m-2)==0)
   grades+=1
# Independent finite quotient: saturate x at N+1, retain all h<=N states.
# Positive use of saturated x necessarily exceeds N; zero use is exact.
N=67;seed=(1,0,0,1,0);q=collections.deque([seed]);seen={seed};values={1}
while q:
 h,x,y,d,m=q.popleft()
 children=[(h,2*x+d,d,d,m),(h+x+y,2*x,0,d,m)]
 if m<2:children += [(*transfer((h,x,y,d),k),m+1) for k in ks]
 for a,b,c,d,e in children:
  if a>N:continue
  values.add(a);s=(a,min(b,N+1),c,d,e)
  if s not in seen:seen.add(s);q.append(s)
gaps=[i for i in range(1,N+1,2) if i%3 and i not in values];need(gaps==[7,55,67])
print(json.dumps({'status':'PASS','physical_full_state_pair_checks':checks,'core_profiles':ks,'general_grade_steps_c_2_through_12':grades,'independent_saturated_states':len(seen),'coprime_three_odd_gaps':gaps,'coprime_three_odd_values':[i for i in sorted(values) if i%2 and i%3]},indent=2))
