#!/usr/bin/env python3
import functools,itertools,json,math
def require(condition):
 if not condition:raise RuntimeError('verification condition failed')

stats={'tournaments':0,'strong_tournaments':0,'rooted_spines':0,'steps':0,'exceptional_steps':0,'tail_bearing_steps':0,'max_exceptional_stages':0}
witness={}
for n in range(1,7):
 edges=list(itertools.combinations(range(n),2))
 for code in range(1<<len(edges)):
  rows=[0]*n
  for i,(a,b) in enumerate(edges):
   if code>>i&1:rows[a]|=1<<b
   else:rows[b]|=1<<a
  stats['tournaments']+=1
  @functools.lru_cache(None)
  def start(mask,v):
   rest=mask^(1<<v)
   if not rest:return 1
   return sum(start(rest,w) for w in range(n) if rest>>w&1 and rows[v]>>w&1)
  def h(mask):return sum(start(mask,v) for v in range(n) if mask>>v&1) if mask else 1
  def components(mask):
   out=[]
   while mask:
    reach={v:(rows[v]&mask)|(1<<v) for v in range(n) if mask>>v&1}
    for z in range(n):
     if mask>>z&1:
      for v in reach:
       if reach[v]>>z&1:reach[v]|=reach[z]
    C=sum(1<<v for v in reach if reach[v]==mask)
    require(C)
    out.append(C);mask^=C
   return out
  full=(1<<n)-1
  if len(components(full))!=1:continue
  stats['strong_tournaments']+=1
  for root in range(n):
   stats['rooted_spines']+=1
   mask=full;k0=start(mask,root);charge=1;E=Et=Em=B=0
   while mask.bit_count()>1:
    k=start(mask,root);parts=components(mask^(1<<root));C=parts[0]
    q=math.prod(h(D) for D in parts[1:])
    A=[a for a in range(n) if C>>a&1 and rows[root]>>a&1];m=len(A)
    require(A and k==q*sum(start(C,a) for a in A))
    a=min(A,key=lambda a:(start(C,a),a));kp=start(C,a)
    require(q*m*kp<=k)
    charge*=q*m;require(charge*kp<=k0)
    stats['steps']+=1
    if q*m>1:
     E+=1;stats['exceptional_steps']+=1
     if q>1:Et+=1;stats['tail_bearing_steps']+=1
     else:Em+=1
    else:require(kp==k)
    B+=sum(D.bit_count()>1 for D in parts[1:])
    require(2**E<=k0 and 3**Et*2**Em<=k0 and 3**B<=k0)
    require(m*m<=4*(k//q)-3)
    mask,root=C,a
   require(start(mask,root)==1 and charge<=k0)
   stats['max_exceptional_stages']=max(stats['max_exceptional_stages'],E)
   if E and 'exceptional' not in witness:witness['exceptional']={'rows':rows,'start_budget':k0,'stages':E}
   if B and 'tail' not in witness:witness['tail']={'rows':rows,'start_budget':k0,'tail_blocks':B}
print(json.dumps({'scope':'All labelled tournaments of orders 1 through 6, all roots of strong tournaments, full selected minimum-start spines','stats':stats,'controls':witness,'status':'PASS'},sort_keys=True,indent=2))
