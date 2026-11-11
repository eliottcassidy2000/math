"""Exact port interaction and intrinsically decoded tournament CA family."""
import json, math
from collections import Counter, defaultdict
from itertools import permutations, product
from pathlib import Path
from probe import block,H,join

WEIGHTS={(0,0):215,(0,1):179,(1,0):245,(1,1):201}

def check(condition, message):
 """Checks remain active under python -O."""
 if not condition:
  raise RuntimeError(message)

def subset_paths(a):
 n=len(a); dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    if t:dp[s][v]=sum(dp[t][u] for u in range(n) if t>>u&1 and a[u][v])
 return dp

def ports(a,v):
 n=len(a);end=subset_paths(a);start=subset_paths(list(map(list,zip(*a))));full=(1<<n)-1
 h=[sum(x) for x in end];h[0]=1
 L=sum(start[s][v]*h[full^s] for s in range(1,full+1) if s>>v&1)
 R=sum(end[s][v]*h[full^s] for s in range(1,full+1) if s>>v&1)
 return dict(H=h[full],end=end[full][v],start=start[full][v],L=L,R=R)

def strong(a):
 for v in range(len(a)):
  seen={v};todo=[v]
  while todo:
   u=todo.pop()
   for w in range(len(a)):
    if a[u][w] and w not in seen:seen.add(w);todo.append(w)
  if len(seen)!=len(a):return False
 return True

def automorphisms(a):
 groups=defaultdict(list)
 for i,row in enumerate(a):groups[sum(row)].append(i)
 buckets=list(groups.values());domain=[v for b in buckets for v in b];count=0
 for choices in product(*(permutations(b) for b in buckets)):
  p=dict(zip(domain,[v for b in choices for v in b]))
  if all(a[i][j]==a[p[i]][p[j]] for i in range(len(a)) for j in range(len(a))):count+=1
 return count,math.prod(math.factorial(len(b)) for b in buckets)

def covers(a):
 dp=subset_paths(a);h=[sum(x) for x in dp];n=len(a)
 f=[[0]*(n+1) for _ in range(1<<n)];f[0][0]=1
 for mask in range(1,1<<n):
  bit=mask&-mask;sub=mask
  while sub:
   if sub&bit:
    for j in range(1,n+1):f[mask][j]+=h[sub]*f[mask^sub][j-1]
   sub=(sub-1)&mask
 return f[-1][1:]

def pair_counts(x):return Counter((x[i],x[(i+1)%len(x)]) for i in range(len(x)))
def evaluate(x):return math.prod(WEIGHTS[p]**m for p,m in pair_counts(x).items())
def step(x,rule=110):
 n=len(x)
 return tuple((rule>>(4*x[(i-1)%n]+2*x[i]+x[(i+1)%n]))&1 for i in range(n))

def inventory(n):
 d={215**n:1,201**n:1}
 for n0 in range(1,n):
  n1=n-n0
  for k in range(1,min(n0,n1)+1):
   h=215**(n0-k)*201**(n1-k)*(179*245)**k
   numer=n*math.comb(n0-1,k-1)*math.comb(n1-1,k-1)
   check(numer%k==0, 'Failed check: numer%k==0')
   check(h not in d, 'Failed check: h not in d')
   d[h]=numer//k
 return d

if __name__=='__main__':
 seeds=[block(8),block(10)];p=[ports(a,0) for a in seeds]
 check([covers(a) for a in seeds]==[[9,33,31,10,1]]*2, 'Failed check: [covers(a) for a in seeds]==[[9,33,31,10,1]]*2')
 controls=[]
 for x in range(2):
  for y in range(2):
   a=join((x,y),[(0,0,0)]);actual=H(a)
   predicted=p[x]['H']*p[y]['H']-p[x]['end']*p[y]['start']+p[x]['L']*p[y]['R']
   check(predicted==actual==WEIGHTS[x,y] and strong(a), 'Failed check: predicted==actual==WEIGHTS[x,y] and strong(a)')
   aut,candidates=automorphisms(a)
   check(aut==1, 'Failed check: aut==1')
   controls.append(dict(bits=[x,y],H=actual,seam_formula=predicted,strong=True,automorphism_order=aut,degree_preserving_candidates=candidates))
 census=[]
 for n in range(1,13):
  brute=Counter(evaluate(tuple((j>>i)&1 for i in range(n))) for j in range(1<<n))
  exact=inventory(n)
  check(brute==exact and len(exact)==n*n//4+2 and sum(exact.values())==1<<n, 'Failed check: brute==exact and len(exact)==n*n//4+2 and sum(exact.values())==1<<n')
  census.append(dict(cells=n,vertices=10*n,classes=1<<n,distinct_H=len(exact)))
 # The census is not a universality test. A finite-rule trace is only a literal example.
 x=tuple(map(int,'000000000000000000000001'));trace=[]
 for t in range(25):
  trace.append(dict(t=t,state=''.join(map(str,x)),H=str(evaluate(x)),pair_counts={''.join(map(str,k)):v for k,v in sorted(pair_counts(x).items())}))
  x=step(x)
 # Exact single-bit response and finite-cone controls for all rules on width7.
 checks=0
 for rule in [0,30,90,110,255]:
  for integer in range(128):
   x=tuple((integer>>i)&1 for i in range(7));y=list(x);y[3]^=1;y=tuple(y)
   for t in range(4):
    affected={i for i in range(7) if x[i]!=y[i]}
    check(all(min((i-3)%7,(3-i)%7)<=t for i in affected), 'Failed check: all(min((i-3)%7,(3-i)%7)<=t for i in affected)')
    checks+=1;x=step(x,rule);y=step(y,rule)
 # Polynomial aggregate rather than enumerating 2^100 configurations.
 inv100=inventory(100)
 check(len(inv100)==2502 and sum(inv100.values())==2**100, 'Failed check: len(inv100)==2502 and sum(inv100.values())==2**100')
 witnesses=[tuple(map(int,w)) for w in ['101000','100100']]
 witness=[dict(word=''.join(map(str,x)),H=evaluate(x),next_H=evaluate(step(x))) for x in witnesses]
 check(witness[0]['H']==witness[1]['H'] and witness[0]['next_H']!=witness[1]['next_H'], "Failed check: witness[0]['H']==witness[1]['H'] and witness[0]['next_H']!=witness[1]['next_H']")
 # Transfer trace equals the weighted sum of all class Hamilton counts.
 m=[[215,179],[245,201]];power=[[1,0],[0,1]]
 for n in range(1,13):
  power=[[sum(power[i][k]*m[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
  check(power[0][0]+power[1][1]==sum(h*count for h,count in inventory(n).items()), 'Failed check: power[0][0]+power[1][1]==sum(h*count for h,count in inventory(n).items())')
 out=dict(status='FINITE-EXACT',scalar_nonclosure=witness,seed_profiles=[covers(a) for a in seeds],port0=p,seam_controls=controls,
          mixed=-8,multiplicative_coupling=215*201-179*245,census=census,
          width100=dict(vertices=1000,classes=str(2**100),distinct_H=len(inv100),multiplicity_sum=str(sum(inv100.values()))),
          rule110_trace=trace,finite_cone_controls=checks)
 Path(__file__).with_name('ca_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:v for k,v in out.items() if k!='rule110_trace'},indent=2))
