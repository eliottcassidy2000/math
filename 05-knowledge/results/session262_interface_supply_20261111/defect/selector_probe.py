#!/usr/bin/env python3
import itertools,json,pathlib,math
HERE=pathlib.Path(__file__).resolve().parent

def supports(word):return sum(1<<i for i,v in enumerate(word) if v)
def certify(f,k):
 if not f:return None
 J=0
 for s in f:J|=s
 expected={s for s in range(1<<k) if not s&~J}
 return J if f==expected else None
predicates=0;accepted=0;images=0
for k in range(5):
 words=list(itertools.product(range(3),repeat=k))
 for J in range(1<<k):
  image={w for w in words if not supports(w)&~J}
  assert len(image)==3**J.bit_count()
  # Every arbitrary per-block induced deletion has exactly its intact-block mask.
  for sizes in itertools.product(range(4),repeat=k):
   intact=sum(1<<i for i,n in enumerate(sizes) if n==3)
   if intact!=J:continue
   actual=set(itertools.product(*[(0,1,2) if n==3 else (0,) for n in sizes]))
   assert actual==image;images+=1
 if k<=3:
  for code in range(1<<(1<<k)):
   f={s for s in range(1<<k) if code>>s&1};predicates+=1
   direct=[J for J in range(1<<k) if f=={s for s in range(1<<k) if not s&~J}]
   # Separately apply the downward/union axioms.
   down=all(t in f for s in f for t in range(1<<k) if not t&~s)
   union=all(a|b in f for a in f for b in f)
   axioms=bool(f) and down and union
   assert axioms==bool(direct)
   assert (certify(f,k) is not None)==axioms
   if axioms:accepted+=1
assert certify({0,1,2},2) is None

def literal_q(rows,r):
 n=len(rows);total=0
 for assn in itertools.product(range(r),repeat=n):
  prod=1
  for slot in range(r):
   vs=[v for v in range(n) if assn[v]==slot]
   prod*=sum(all(rows[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in itertools.permutations(vs))
  total+=prod
 return total
seed_rows=[]
for r in range(1,6):
 e=literal_q([2,0],r);c=literal_q([2,4,1],r)
 assert e==r*r and c==r**3+2*r
 seed_rows.append({'r':r,'edge':e,'triangle':c})

def join(x,y):
 n=len(x);m=len(y)
 return [v|(((1<<m)-1)<<n) for v in x]+[v<<n for v in y]
def induced(rows,keep):return [sum(1<<j for j,w in enumerate(keep) if rows[v]>>w&1) for v in keep]
def subset_q(rows,r):
 n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   if not mask>>v&1:continue
   prev=mask^(1<<v)
   if prev:dp[mask][v]=sum(dp[prev][u] for u in range(n) if rows[u]>>v&1)
  h[mask]=sum(dp[mask])
 q=[0]*N;q[0]=1
 for _ in range(r):
  nxt=[0]*N
  for mask in range(N):
   s=mask
   while True:
    nxt[mask]+=q[mask^s]*h[s]
    if not s:break
    s=(s-1)&mask
  q=nxt
 return q[-1]
base=[2,4,1];full=join(join(base,base),base);profiles=[]
for J in range(4):
 keep=[v for v in range(9) if v not in [3+3*i for i in range(2) if not J>>i&1]]
 rows=induced(full,keep);vals=[]
 for r in range(1,5):
  val=subset_q(rows,r);expected=(r**3+2*r)**(1+J.bit_count())*(r*r)**(2-J.bit_count())
  assert val==expected;vals.append(val)
 profiles.append({'J':J,'Q1_to_Q4':vals})
assert profiles[1]['Q1_to_Q4']==profiles[2]['Q1_to_Q4']
# The newly constructed correlated pair gate and its leakage boundary.
diamond=[2,12,9,1]
def odd_cycles(rows):
 out=[];n=len(rows)
 for size in range(3,n+1,2):
  for vs in itertools.combinations(range(n),size):
   for tail in itertools.permutations(vs[1:]):
    cyc=(vs[0],)+tail
    if all(rows[cyc[i]]>>cyc[(i+1)%size]&1 for i in range(size)):out.append(cyc)
 return out
assert odd_cycles(diamond)==[(0,1,2),(0,1,3)]
gate_profiles=[]
for r in range(1,6):
 val=literal_q(diamond,r)
 assert val==r**4+4*r*r
 assert subset_q(join(diamond,diamond),r)==val*val
 gate_profiles.append({'r':r,'Q_r':val})
fan_receipt=[]
for m in range(1,6):
 rows=[2, sum(1<<i for i in range(2,m+2))]
 rows += [1+sum(1<<j for j in range(i+1,m+2)) for i in range(2,m+2)]
 cycles=odd_cycles(rows)
 expected=[(0,1)+vs for size in range(1,m+1,2) for vs in itertools.combinations(range(2,m+2),size)]
 assert set(cycles)==set(expected)
 val=subset_q(rows,1)
 assert val==1+2**m
 fan_receipt.append({'m':m,'H':val,'odd_cycles':len(cycles),'cycle_lengths':sorted(map(len,cycles))})
# Gate word decoder has exactly the desired correlated image.
gate_image={(0,0),(1,0),(2,0),(0,1),(0,2)}
assert gate_image=={w for w in itertools.product(range(3),repeat=2) if not(w[0] and w[1])}
receipt={'status':'FINITE-EXACT','shell_residual_size_patterns':images,'support_predicates_tested':predicates,'principal_predicates':accepted,'literal_initial_covers':seed_rows,'two_shell_profiles_on_triangle_core':profiles,'diamond_odd_cycles':odd_cycles(diamond),'diamond_cover_profiles':gate_profiles,'shared_edge_fans':fan_receipt,'first_induced_deletion_hostile':{'supports':[0,1,2],'states':5,'certificate':'supports 1 and 2 accepted, their union 3 rejected'}}
(HERE/'selector_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
