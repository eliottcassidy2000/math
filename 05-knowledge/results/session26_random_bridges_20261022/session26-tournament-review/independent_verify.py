"""Independent controls; no imports from producer. Run python3 independent_verify.py."""
import collections,itertools,json,math,hashlib
from pathlib import Path

def need(test,msg):
 if not test: raise RuntimeError(msg)
def advance(s,b,m=None):
 h,a,p=s
 if b==0: h,a,p=h+a+p,a*2,0
 else: a,p=a*2+1,1
 return (h,a,p) if m is None else (h%m,a%m,p)
def evalword(w):
 s=(1,0,0)
 for b in w:s=advance(s,b)
 return s

def graphcount(w):
 # Build explicit adjacency, then forward subset DP on a dense array.
 n=len(w)+1
 adj=[[False]*n for _ in range(n)]
 for u in range(n-1):
  for v in range(u+1,n-1):adj[u][v]=True
  adj[n-1][u]=bool(w[u]);adj[u][n-1]=not bool(w[u])
 counts=[[0]*n for _ in range(1<<n)]
 for u in range(n):counts[1<<u][u]=1
 for mask in range(1,1<<n):
  for last,c in enumerate(counts[mask]):
   if c:
    for nxt in range(n):
     if not mask>>nxt&1 and adj[last][nxt]:counts[mask|1<<nxt][nxt]+=c
 return sum(counts[-1])

def det(a):
 if not a:return 1
 return sum((-1)**j*x*det([r[:j]+r[j+1:] for r in a[1:]]) for j,x in enumerate(a[0]))
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)]for row in a]
def transpose(a):return list(map(list,zip(*a)))
def minors_gcd(a,k):
 g=0
 for rs in itertools.combinations(range(len(a)),k):
  for cs in itertools.combinations(range(len(a[0])),k):g=math.gcd(g,det([[a[r][c]for c in cs]for r in rs]))
 return g
nwords=0
for n in range(10):
 for w in itertools.product((0,1),repeat=n):
  need(evalword(w)[0]==graphcount(w),('graph',w));nwords+=1
stats=[]
for m in range(1,101):
 start=(1%m,0,0);seen={start};todo=collections.deque([start])
 while todo:
  s=todo.popleft()
  for b in (0,1):
   t=advance(s,b,m)
   if t not in seen:seen.add(t);todo.append(t)
 expected={(h,a,p) for h in range(m) for a in range(m) for p in (0,1) if m%2 or (h%2==1 and a%2==p)}
 need(seen==expected,('reachable set',m))
 # Iterative equivalence partitions retaining complete observable traces.
 states=sorted(seen);index={s:i for i,s in enumerate(states)}
 dest=[[index[advance(s,b,m)]for b in (0,1)]for s in states]
 labels=[s[0]for s in states];rounds=0
 while True:
  sigs=[(s[0],labels[d[0]],labels[d[1]])for s,d in zip(states,dest)]
  ids={};new=[ids.setdefault(t,len(ids))for t in sigs];rounds+=1
  if len(set(new))==len(set(labels)):break
  labels=new
 strong_obs={((h+a+p)%m,(h+2*a+2)%m,(h+4*a+4)%m) for h,a,p in seen}
 want=1 if m<=2 else len(expected)
 need(len(strong_obs)==want,('strong probes',m))
 need(len(ids)==want,('minimum',m));stats.append([m,len(seen),len(ids),rounds])
 # Check constructive reset/suffix certificate on every target, via modular states.
 if m<=30:
  k=0;odd=m
  while odd%2==0:odd//=2;k+=1
  L=max(k,1)
  while pow(2,L,odd)!=1%odd:L+=1
  block=[1]*L+[0]*L;reset=[];s=start
  for _ in range(m//math.gcd(2,m)):
   for b in block:s=advance(s,b,m)
   need(s[1:]==(0,0),('reset',m));reset.append(s)
  built=set()
  for a in range(m):
   for p in (0,1):
    if not m%2 and p!=a%2:continue
    b=a if a%2==p else a+m
    word=list(map(int,bin(b)[2:]))
    for s in reset:
     for digit in word:s=advance(s,digit,m)
     built.add(s)
  need(built==expected,('constructive',m))
crt_pairs=0
for m in range(3,13):
 for n in range(m+1,13):
  if math.gcd(m,n)!=1:continue
  def legal(q):return {(h,a,p)for h in range(q)for a in range(q)for p in (0,1)if q%2 or(h%2 and a%2==p)}
  reduction={((h%m,a%m,p),(h%n,a%n,p))for h,a,p in legal(m*n)}
  fiber={(s,t)for s in legal(m)for t in legal(n)if s[2]==t[2]}
  need(reduction==fiber,('CRT',m,n));crt_pairs+=1
U=[(),(1,),(1,0),(1,1)];V=[(),(0,),(1,0),(0,0)]
B=[[evalword(u+v)[0]for v in V]for u in U]
need([minors_gcd(B,k)for k in range(1,5)]==[1,2,4,16],'minor gcds')
R=[list(evalword(u)) + [1]for u in U]
C=[[1,1,1,1],[0,1,2,3],[0,1,0,1],[0,0,2,0]]
Ts=[[[1,1,1,0],[0,2,0,0],[0,0,0,0],[0,0,0,1]],[[1,0,0,0],[0,2,0,1],[0,0,0,1],[0,0,0,1]]]
row_actions=[[[1,0,0,0],[0,0,1,0],[-1,0,2,0],[-1,-1,2,1]],[[0,1,0,0],[0,0,0,1],[-1,0,1,1],[0,-2,0,3]]]
col_actions=[[[0,0,-1,0],[1,0,0,-2],[0,0,1,0],[0,1,1,3]],[[1,0,-1,-1],[0,0,0,-1],[0,1,2,2],[0,0,0,1]]]
for T,ra,ca in zip(Ts,row_actions,col_actions):
 need(mm(R,transpose(T))==mm(ra,R),'row lattice closure')
 need(mm(transpose(T),C)==mm(C,ca),'column lattice closure')
need(mm(R,C)==B,'Hankel factorization')
print(json.dumps({'status':'PASS','graph_words_lengths_0_through_9':nwords,'modulus_reach_min_rounds':stats,'constructive_moduli':'1..30','strong_probe_moduli':'1..100','CRT_pairs':crt_pairs,'minor':B,'minor_determinant':det(B),'determinantal_divisors':[minors_gcd(B,k)for k in range(1,5)],'integral_lattice_closure':'PASS','producer_proof_sha256':hashlib.sha256((Path(__file__).resolve().parent.parent/'session26-tournament'/'PROOF.md').read_bytes()).hexdigest()},indent=2))
