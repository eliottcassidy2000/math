from pathlib import Path
import json
P=Path(__file__).parent

def require(condition, message):
 if not condition:
  raise RuntimeError(message)

def adj(code,n):
 a=[0]*n;b=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>b&1:a[i]|=1<<j
   else:a[j]|=1<<i
   b+=1
 return a

def induced(a,vs):
 return [sum(((a[v]>>u)&1)<<j for j,u in enumerate(vs)) for v in vs]
def modules(a):
 n=len(a); full=(1<<n)-1; out=[]
 for s in range(1,full):
  if s.bit_count()<2:continue
  if all((a[v]&s) in (0,s) for v in range(n) if not(s>>v&1)):out.append(s)
 return out
def prime(a):return len(a)>=3 and not modules(a)
def hp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for s in range(1,1<<n):
  for v,x in enumerate(d[s]):
   z=a[v]&~s
   while x and z:
    b=z&-z;z-=b;d[s|b][b.bit_length()-1]+=x
 return sum(d[-1])
def w(k):
 n=2*k+1;a=[0]*n
 for i in range(n-1):
  for j in range(i+1,n-1):a[i]|=1<<j
  if i%2:a[i]|=1<<(n-1)
  else:a[n-1]|=1<<i
 return a
old={1};stats=[];witnesses=[]
for n in range(1,9):
 rows=[]
 for line in (P/f'classes_{n}.txt').read_text().splitlines():
  if line.startswith('#'):continue
  t=line.split();rows.append((int(t[0]),int(t[1])))
 vals={h for c,h in rows};pr=non=crit=0;orph=[]
 for c,h in rows:
  a=adj(c,n)
  require(hp(a)==h, f"Hamilton count mismatch: n={n}, code={c}, expected={h}")
  if prime(a):
   pr+=1
   if h in old:
    non+=1
    pdel=[v for v in range(n) if prime(induced(a,[u for u in range(n) if u!=v]))]
    if not pdel:crit+=1
    if n<=7:witnesses.append(dict(n=n,code=c,H=h,prime_deletions=pdel))
  if n==8 and h==63:orph.append(dict(code=c,modules=modules(a),deletion_H=[hp(induced(a,[u for u in range(n) if u!=v])) for v in range(n)]))
 stats.append(dict(n=n,classes=len(rows),primes=pr,noncreator_primes=non,critical_noncreator_primes=crit))
 if orph:print('ORPHANS',json.dumps(orph))
 old|=vals
print('STATS',json.dumps(stats,indent=2));print('WITNESSES',json.dumps(witnesses,indent=2))
for k in range(1,7):
 a=w(k);h=hp(a);require(h==(4**(k+1)+6*k+5)//9, f"W count formula mismatch: k={k}")
 print('W',len(a),h,'prime',prime(a),'prime_deletions',sum(prime(induced(a,[u for u in range(len(a)) if u!=v])) for v in range(len(a))))
from itertools import permutations
for c in [135528708,36175876]:
 a=adj(c,8);ms=modules(a)
 blocks=[[v for v in range(8) if s>>v&1] for s in ms]
 blocks += [[v] for v in range(8) if not any(s>>v&1 for s in ms)]
 q=induced(a,[b[0] for b in blocks]);base=w(2)
 maps=[perm for perm in permutations(range(5)) if induced(q,perm)==base]
 require(bool(maps) and prime(q) and hp(q)==9, f"Orphan quotient mismatch: code={c}")
 perm=maps[0]
 print('ORPHAN_FACTOR',c,'W5-blocks',[blocks[i] for i in perm],'size_vector',[len(blocks[i]) for i in perm])
a=adj(6241,6)
require(hp(a)==31, "Six-vertex H31 witness mismatch")
print('H31_SMALLER_WITNESS',6241,a)
def inflate_chain_children(q,sizes):
 starts=[];t=0
 for s in sizes:starts.append(t);t+=s
 a=[0]*t
 for i,s in enumerate(sizes):
  for x in range(starts[i],starts[i]+s):
   for y in range(x+1,starts[i]+s):a[x]|=1<<y
   for j,u in enumerate(sizes):
    if q[i]>>j&1:
     for y in range(starts[j],starts[j]+u):a[x]|=1<<y
 return a
for exponent in range(4,13):
 r=exponent//2
 sizes=[r-1,2,r-1,1,1] if exponent%2==0 else [r-1,1,1,1,r-1,1,1]
 q=w(2 if exponent%2==0 else 3)
 a=inflate_chain_children(q,sizes)
 # Direct F_exponent with transitive chain first, hub last.
 n=exponent+2;fan=[0]*n
 for i in range(n-1):
  for j in range(i+1,n-1):fan[i]|=1<<j
  if i+1 in [r,exponent+1-r,exponent+1]:fan[i]|=1<<(n-1)
  else:fan[n-1]|=1<<i
 require(a==fan, f"Fan adjacency mismatch: exponent={exponent}")
 require(hp(a)==2**exponent-1, f"Fan Hamilton count mismatch: exponent={exponent}")
 print('FAN',exponent,'kernel',len(q),'sizes',sizes,'H',hp(a))
def cyclic(r):
 n=2*r+1
 return [sum(1<<j for j in range(n) if 1<=(j-i)%n<=r) for i in range(n)]
def staircase(r):
 n=2*r+1;a=[0]*n
 for i in range(r):
  for j in range(i):a[i]|=1<<j
 for i in range(r+1):
  for j in range(i+1,r+1):a[r+i]|=1<<(r+j)
 for i in range(r):
  for j in range(r+1):
   if j<=i:a[i]|=1<<(r+j)
   else:a[r+j]|=1<<i
 return a
require(sorted([hp(w(2)),hp(cyclic(2)),hp(staircase(2))])==[9,13,15], 'Critical order-five H base mismatch')
for r in range(2,9):
 t=cyclic(r);u=staircase(r)
 require(induced(t,[v for v in range(2*r+1) if v not in [r,2*r]])==cyclic(r-1), f'Cyclic paired deletion mismatch: r={r}')
 require(induced(u,[v for v in range(2*r+1) if v not in [r-1,2*r]])==staircase(r-1), f'Staircase paired deletion mismatch: r={r}')
 require(induced(w(r),list(range(2*r-2))+[2*r])==w(r-1), f'Alternating-hub paired deletion mismatch: r={r}')
 print('CRITICAL_PAIR_DELETIONS',2*r+1,'T,U,W verified')
