"""Independent signed-action and tournament checks; no producer imports."""
from itertools import product
from math import gcd
from collections import Counter
import json

def req(ok,msg):
 if not ok: raise RuntimeError(msg)
def factors(n):
 out={};p=2
 while p*p<=n:
  while n%p==0:out[p]=out.get(p,0)+1;n//=p
  p+=1
 if n>1:out[n]=out.get(n,0)+1
 return out
def leg(x,p):
 t=pow(x%p,(p-1)//2,p)
 return 1 if t==1 else -1

def chars(N):
 U=[u for u in range(N) if gcd(u,N)==1]
 squares={u*u%N for u in U};unseen=set(U);cosets=[]
 while unseen:
  a=min(unseen);C={a*s%N for s in squares};cosets.append(C);unseen-=C
 idx={x:i for i,C in enumerate(cosets) for x in C}
 req(idx[1]==0,'identity coset')
 result=[]
 for signs in product((-1,1),repeat=len(cosets)-1):
  v=(1,)+signs
  if all(v[idx[a*b%N]]==v[idx[a]]*v[idx[b]] for a in U for b in U):
   result.append({x:v[idx[x]] for x in U})
 return U,result

def components(N,U,chi):
 # Independent graph: add every unit-action constraint and a separate oddness edge.
 unused=set(range(1,N));free=[];forced=[]
 while unused:
  start=min(unused);labels={start:1};todo=[start];bad=False
  while todo:
   x=todo.pop()
   for y,s in [(u*x%N,chi[u]) for u in U]+[((-x)%N,-1)]:
    value=labels[x]*s
    if y in labels:bad|=labels[y]!=value
    else:labels[y]=value;todo.append(y)
  unused-=set(labels)
  (forced if bad else free).append(labels)
 return free,forced

def poly(sizes):
 c=Counter({0:1})
 for size in sizes:
  nxt=Counter(c)
  for k,v in c.items():nxt[k+size]+=2*v
  c=nxt
 return dict(c)

char_count=literal_words=action_components=0
for N in range(3,106,2):
 U,characters=chars(N);fac=factors(N)
 req(len(characters)==2**len(fac),'character completeness')
 for chi in characters:
  char_count+=1;free,forced=components(N,U,chi);action_components+=len(free)+len(forced)
  E=[]
  # Identify conductor by direct testing of every candidate divisor, not Legendre subsets.
  divisors=[d for d in range(1,N+1) if N%d==0]
  through=[d for d in divisors if all(chi[u]==1 for u in U if u%d==1%d)]
  f=min(through)
  req(all(a==1 for a in factors(f).values()),'squarefree conductor')
  if chi[N-1]==1: req(not free,'even character only zero')
  else:
   eligible=[d for d in divisors if d>1 and d%f==0]
   req(sorted(len(C) for C in free)==sorted(sum(gcd(x,d)==1 for x in range(d)) for d in eligible),'support strata sizes')
   predicted=1
   for p,a in fac.items():predicted*=a if f%p==0 else a+1
   req(len(free)==predicted,'D formula')
   for C in free:
    orders={N//gcd(x,N) for x in C};req(len(orders)==1 and next(iter(orders))%f==0,'literal component stratum')
   if not forced:req(len(fac)==1 and next(iter(fac))%4==3,'full support necessity')
  if N<=11:
   counts=Counter()
   for vals in product((-1,0,1),repeat=(N-1)//2):
    literal_words+=1;sigma={0:0}
    for x,z in enumerate(vals,1):sigma[x]=z;sigma[N-x]=-z
    if all(sigma[u*x%N]==chi[u]*sigma[x] for u in U for x in range(N)):
     counts[sum(sigma[x]!=0 for x in range(1,N))]+=1
   req(dict(counts)==poly([len(C) for C in free]),'literal support polynomial')
# Minimal mixed-prime hostile.
req(11*3%15==3 and leg(11,3)==-1,'N15 hostile')

def symbol(x,p,a,eps):
 N=p**a;x%=N
 if not x:return 0
 v=0
 while x%p==0:x//=p;v+=1
 return eps[v]*leg(x,p)
def mapping(x,p,a,ts):
 out=0;q=1
 for t in ts:out+=(t*(x%p)%p)*q;x//=p;q*=p
 return out
def hp_count(adj):
 n=len(adj);dp={1<<i:{i:1} for i in range(n)}
 for mask in range(1,1<<n):
  if mask not in dp:continue
  for v,count in list(dp[mask].items()):
   for w in range(n):
    if not mask>>w&1 and adj[v][w]:
     nxt=mask|1<<w;d=dp.setdefault(nxt,{})
     d[w]=d.get(w,0)+count
 return sum(dp[(1<<n)-1].values())
frame_pairs=arc_tests=0;ham=[]
for p,a in [(3,1),(7,1),(3,2),(11,1),(3,3),(7,2)]:
 N=p**a;frames=list(product((-1,1),repeat=a))
 for eps,other in product(frames,repeat=2):
  frame_pairs+=1
  ts=[next(t for t in range(1,p) if leg(t,p)==eps[i]*other[i]) for i in range(a)]
  image=[mapping(x,p,a,ts) for x in range(N)]
  req(len(set(image))==N,'digit bijection')
  for x in range(N):
   for y in range(N):
    req(symbol(y-x,p,a,eps)==symbol(image[y]-image[x],p,a,other),'frame isomorphism')
    arc_tests+=1
 if N<=11:
  adj=[[symbol(y-x,p,a,(1,)*a)==1 for y in range(N)] for x in range(N)]
  H=hp_count(adj);order=(p*(p-1)//2)**((N-1)//(p-1))
  req(H%order==0,'Hamiltonian divisibility')
  ham.append({'p':p,'a':a,'vertices':N,'paths':H,'subgroup_order':order})
# Exhaust the actual 3^4-node-choice permutations for the ternary square.
perms=set()
for root,*shifts in product(range(3),repeat=4):
 image=tuple((x%3+root)%3+3*((x//3+shifts[x%3])%3) for x in range(9))
 req(len(set(image))==9,'tree action bijective')
 for x in range(9):
  for y in range(9):req(symbol(y-x,3,2,(1,1))==symbol(image[y]-image[x],3,2,(1,1)),'tree automorphism')
 perms.add(image)
req(len(perms)==81,'faithful subgroup size')
# Eighth-character and Boolean-boundary independent direct checks.
U8=[1,3,5,7];rows=[[1,1,1,1],[1,-1,1,-1],[1,-1,-1,1],[1,1,-1,-1]]
for i,j in product(range(4),repeat=2):req(sum(a*b for a,b in zip(rows[i],rows[j]))==(4 if i==j else 0),'eighth table')
boolean=0
for n in range(2,5):
 N=1<<n
 for word in range(1<<N):
  if word.bit_count()%2!=1:continue
  req(word.bit_count()!=N//2,'no complement permutation')
  for h in range(1,N):
   vals={((word>>x)^(word>>(x^h)))&1 for x in range(N)}
   req(len(vals)==2,'nonconstant derivative');boolean+=1
# Separate odd characters on each additive-order stratum: literal N=21 audit.
N=21;U=[u for u in range(N) if gcd(u,N)==1];divs=[3,7,21]
options=[]
for m in divs:
 Um,cm=chars(m)
 options.append([(chi,e) for chi in cm if chi[m-1]==-1 for e in (-1,1)])
family=[];trace_values=Counter()
for choice in product(*options):
 sigma=[0]*N
 for m,(chi,e) in zip(divs,choice):
  for t in chi:sigma[(N//m)*t%N]=e*chi[t]
 req(all(sigma[x] in (-1,1) and sigma[-x]==-sigma[x] for x in range(1,N)),'stratum odd full support')
 for m,(chi,e) in zip(divs,choice):
  for t in chi:
   x=(N//m)*t%N
   for u in U:req(sigma[u*x%N]==chi[u%m]*sigma[x],'stratum covariance')
 stabilizer={u for u in U if all(sigma[u*x%N]==sigma[x] for x in range(N))}
 reversers={u for u in U if all(sigma[u*x%N]==-sigma[x] for x in range(N))}
 req(stabilizer=={u for u in U if leg(u,3)==leg(u,7)==1},'square stabilizer')
 req(reversers=={u for u in U if leg(u,3)==leg(u,7)==-1},'reverser coset')
 conv=[sum(sigma[t]*sigma[(d-t)%N] for t in range(N)) for d in range(N)]
 tr4=N*sum(z*z for z in conv);trace_values[tr4]+=1
 family.append(tuple(sigma))
req(len(set(family))==16,'separate character count')
remaining=set(family);unit_orbits=[];orbit_traces=[]
while remaining:
 sigma=min(remaining)
 orbit={tuple(sigma[u*x%N] for x in range(N)) for u in U}
 req(orbit<=set(family) and len(orbit)==4,'four-element unit orbit')
 traces=set()
 for sig in orbit:
  conv=[sum(sig[t]*sig[(d-t)%N] for t in range(N)) for d in range(N)]
  traces.add(N*sum(z*z for z in conv))
 req(len(traces)==1,'trace constant on unit orbit')
 orbit_traces.append(next(iter(traces)));unit_orbits.append(orbit);remaining-=orbit
req(len(unit_orbits)==len(set(orbit_traces))==4,'complete isomorphism separation')
print(json.dumps({'status':'PASS','odd_moduli':52,'characters':char_count,'signed_action_components':action_components,'literal_partial_words':literal_words,'frame_pairs':frame_pairs,'frame_arc_tests':arc_tests,'explicit_ternary_subgroup_elements':len(perms),'hamiltonian_controls':ham,'boolean_derivative_checks':boolean,'producer_imports':False,'separate_character_N21':{'symbols':len(family),'trace_four_distribution':dict(sorted(trace_values.items())),'unit_orbits':len(unit_orbits),'orbit_trace_four_values':sorted(orbit_traces),'isomorphism_classes_within_family':4}},indent=2,sort_keys=True))
