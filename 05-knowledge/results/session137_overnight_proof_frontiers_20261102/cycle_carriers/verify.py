from itertools import permutations,product
from collections import Counter,defaultdict
import json, math
from pathlib import Path

def check(c,m):
 if not c: raise RuntimeError(m)
def adj(n,code):
 a=[0]*n
 for k,(i,j) in enumerate((i,j) for i in range(n) for j in range(i+1,n)):
  if code>>k&1:a[i]|=1<<j
  else:a[j]|=1<<i
 return tuple(a)
def canon(a):
 n=len(a); groups=defaultdict(list)
 for i,v in enumerate(a):groups[v.bit_count()].append(i)
 gs=[groups[d] for d in sorted(groups)]
 best=None
 for blocks in product(*(permutations(g) for g in gs)):
  p=sum((tuple(b) for b in blocks),()); code=0;k=0
  for i in range(n):
   for j in range(i+1,n):
    code|=((a[p[i]]>>p[j])&1)<<k;k+=1
  if best is None or code<best:best=code
 return best

def strong(a):
 n=len(a); reach=list(a)
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 return all((reach[i]|1<<i)==(1<<n)-1 for i in range(n))
def cycles(a):
 n=len(a); out=[]
 def rec(path,used):
  last=path[-1]; first=path[0]
  if len(path)>=3 and len(path)%2 and a[last]>>first&1:out.append((tuple(path),used))
  opts=a[last]&~used&~((1<<(first+1))-1)
  while opts:
   bit=opts&-opts;opts-=bit
   rec(path+[bit.bit_length()-1],used|bit)
 for i in range(n):rec([i],1<<i)
 return out

def omega(cs):
 return tuple(sum(1<<j for j,d in enumerate(cs) if i!=j and c[1]&d[1]) for i,c in enumerate(cs))
def graphcanon(a):
 # Use complement, remove universal vertices, and degree-preserving exhaustive permutations.
 n=len(a); full=(1<<n)-1
 comp=[full^((1<<i)|a[i]) for i in range(n)]
 active=[i for i,x in enumerate(comp) if x]
 b=tuple(sum(1<<j for j,v in enumerate(active) if comp[u]>>v&1) for u in active)
 return (n-len(active),len(active),canon(b))
def hp(a):
 n=len(a); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1<<n):
  for i,x in enumerate(dp[s]):
   opts=a[i]&~s
   while x and opts:
    bit=opts&-opts;opts-=bit;dp[s|bit][bit.bit_length()-1]+=x
 return sum(dp[-1])
def independent_poly(a):
 # Include/exclude using memoized masks.
 from functools import lru_cache
 @lru_cache(None)
 def f(s):
  if not s:return (1,)
  v=(s&-s).bit_length()-1;t=s&~(1<<v)
  p=list(f(t));q=f(t&~a[v]);p+=[0]*max(0,len(q)+1-len(p))
  for i,x in enumerate(q):p[i+1]+=x
  return tuple(p)
 return f((1<<len(a))-1)

def main():
 reps={0}; bank={}; allreps={}; clcounts=[]; uq=[]
 for n in range(1,9):
  cand=set()
  for code in reps:
   old=adj(n-1,code)
   for mask in range(1<<(n-1)):
    a=list(old)+[0]
    for i in range(n-1):
     if mask>>i&1:a[i]|=1<<(n-1)
     else:a[-1]|=1<<i
    cand.add(canon(a))
  reps=cand; allreps[n]=sorted(reps); clcounts.append(len(reps))
  stats=Counter(); ustats=Counter()
  for code in sorted(reps):
   a=adj(n,code)
   if not strong(a) or n<3:continue
   cs=cycles(a); om=omega(cs);m=len(cs)
   u=sum(x.bit_count()==m-1 for x in om)
   check(m>=((n-1)**2)//4,'Moon floor');check(u>= (1 if n==3 else 2 if n%2==0 else 4),'universal floor')
   tris=[mask for cyc,mask in cs if len(cyc)==3]
   for i,(cyc,mask) in enumerate(cs):
    transitive=not any(not(mask&t) for t in tris)
    check((om[i].bit_count()==m-1)==transitive,'universal iff transitive complement')
    check(m-1-om[i].bit_count() not in (3,10,14),'height two forbidden disjointness degree')
   stats[m]+=1;ustats[u]+=1
   check(hp(a)==sum(v*(2**i) for i,v in enumerate(independent_poly(om))),'OCF')
   if m<=15:
    key=graphcanon(om)
    value={'key':key,'m':m,'source_n':n,'source_code':code,'source_adj':a,'alpha':independent_poly(om),'H':hp(a),'cycles':[c[0] for c in cs]}
    bank.setdefault(key,value)
  uq.append({'n':n,'strong_cycle_counts':sorted(stats.items()),'universal_counts':sorted(ustats.items())})
 check(clcounts==[1,1,2,4,12,56,456,6880],'class counts')
 out={'class_counts':clcounts,'strong_statistics':uq,'connected_carriers_m_le_15':sorted(bank.values(),key=lambda x:(x['m'],x['key']))}
 out.update(extra_checks(bank,allreps))
 print(json.dumps(out,sort_keys=True,indent=2))
def moon(n):
 a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   if j==i+1:a[i]|=1<<j
   else:a[j]|=1<<i
 return tuple(a)
def addpoly(a,b,shift=0):
 a=list(a)+[0]*max(0,len(b)+shift-len(a))
 for i,x in enumerate(b):a[i+shift]+=x
 return tuple(a)
def context(a):
 n=len(a);b=list(a)+[0,0]
 for i in range(n):b[i]|=1<<n;b[n+1]|=1<<i
 b[n]|=1<<(n+1)
 return tuple(b)
def extra_checks(bank,allreps):
 expected={
 (1,0,0),(2,0,0),(4,0,0),(5,0,0),(4,2,1),(6,0,0),(7,0,0),(8,0,0),
 (5,4,44),(7,2,1),(9,0,0),(6,4,12),(7,3,6),(9,2,1),(11,0,0),
 (4,8,189861952),(8,4,12),(10,2,1),(12,0,0),(10,3,6),(11,2,1),(13,0,0),
 (6,8,297024),(8,6,15632),(10,4,12),(10,4,44),(12,2,1),(9,6,22800),(15,0,0)}
 check(set(bank)==expected,'complete 29-type table')
 check((3,0,0) not in bank and (10,0,0) not in bank and (14,0,0) not in bank and (5,2,1) not in bank,'unrealizable graph controls')
 # Exact bounded-budget semiring closure; each expression compiles by dominance sum.
 bycost=[{} for _ in range(16)];bycost[0][1]=[]
 rows=sorted(bank.values(),key=lambda x:(x['m'],x['key']))
 for cost in range(1,16):
  for r in rows:
   m=r['m']
   if m>cost:continue
   for h,expr in bycost[cost-m].items():
    bycost[cost].setdefault(h*r['H'],expr+[r['key']])
 best={}
 for cost,vals in enumerate(bycost):
  for h,expr in vals.items():best.setdefault(h,{'cost':cost,'components':expr})
 checks={1:0,3:1,5:2,9:2,11:5,13:6,15:3,17:6,19:9,23:9,25:4,27:3,29:10,31:9,35:13,37:14,41:14,53:14}
 for h,cost in checks.items():check(best[h]['cost']==cost,'minimum cycle budget')
 check(7 not in best and 21 not in best,'inherited forbidden outputs')
 # Literal source compilation and direct HP check for every small selected budget.
 constructed=[]
 for h in sorted(checks):
  a=[]
  for key in best[h]['components']:
   old=len(a); b=bank[tuple(key)]['source_adj'];newmask=((1<<len(b))-1)<<old
   a=[row|newmask for row in a]+[row<<old for row in b]
  if not a:a=[0]
  check(hp(a)==h,'compiled source HP')
  check(len(cycles(a))==best[h]['cost'],'compiled source cycle budget')
  constructed.append({'H':h,'cost':best[h]['cost'],'source_adj':a})
 # Complete same-order31 comparison.
 m31={}
 for n in range(1,8):
  vals=[len(cycles(adj(n,c))) for c in allreps[n] if hp(adj(n,c))==31]
  if vals:m31[n]=min(vals)
 check(m31.get(6)==13 and m31.get(7)==9 and all(n>=6 for n in m31),'31 resource tradeoff')
 # Independent interval recurrence and all-order formula controls.
 polys=[(1,)]
 family=[]
 for n in range(1,16):
  p=polys[n-1]
  for length in range(3,n+1,2):p=addpoly(p,polys[n-length],1)
  polys.append(p)
  a=moon(n);cs=cycles(a);om=omega(cs)
  check(all(tuple(range(c[0],c[-1]+1))==c for c,_ in cs),'literal interval cycles')
  check(len(cs)==((n-1)**2)//4,'sharp Moon cycle floor')
  check(independent_poly(om)==p,'interval polynomial recurrence')
  if n<=12:check(hp(a)==sum(x*2**i for i,x in enumerate(p)),'interval OCF direct HP')
  if n%3==0:check(len(p)-1==n//3 and p[-1]==1,'sharp packing budget')
  family.append({'n':n,'alpha':p,'H':sum(x*2**i for i,x in enumerate(p))})
 # Inherited same-order same-Omega context hostile, checked afresh.
 host=[]
 for code in (40,76):
  a=adj(5,code);om=omega(cycles(a));b=context(a)
  check(graphcanon(om)==(7,0,0),'same carrier K7')
  host.append({'code':code,'H':hp(a),'context_H':hp(b),'context_alpha':independent_poly(omega(cycles(b)))})
 check([x['context_H'] for x in host]==[123,135],'inherited context hostile')
 return {'budget_value_counts':[len(x) for x in bycost],
 'exact_budget_values':[sorted(x) for x in bycost],
 'all_minimum_budget_values':{str(h):v for h,v in sorted(best.items())},
 'compiled_selected_witnesses':constructed,'minimum_cycles_for_31_by_order':m31,
 'moon_family':family,'same_carrier_context_hostile':host}
if __name__=='__main__':main()

