"""Independent coefficient/state audit; no producer imports; checks survive -O."""
import itertools as it,json,sys,hashlib
from functools import lru_cache
from pathlib import Path

def check(ok,detail):
 if not ok: raise RuntimeError('AUDIT FAILURE: '+str(detail))
if '--failure-sentinel' in sys.argv: check(False,'deliberate failure sentinel')

def comps(n):
 if n==0: yield (); return
 for c in comps(n-1):
  if c: yield c[:-1]+(c[-1]+1,)
  yield c+(1,)

def plus(a,b):
 c=[0]*max(len(a),len(b))
 for i,v in enumerate(a):c[i]+=v
 for i,v in enumerate(b):c[i]+=v
 while len(c)>1 and c[-1]==0:c.pop()
 return tuple(c)
def times(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,v in enumerate(a):
  for j,w in enumerate(b):c[i+j]+=v*w
 return plus(c,(0,))

def local(m):
 # Count actual local vertex assignments, via parity DP, not character sums.
 M=[[(0,) for j in range(3)] for i in range(3)];M[0][0]=(1,)
 if m==1:return [[(1,),(0,1),(0,)],[(0,),(0,),(1,)],[(1,),(1,),(0,)]]
 subset=[1,0]
 for _ in range(m-1):subset=[sum(subset),sum(subset)]
 subset=subset[::-1] # mandatory endpoint
 for p in range(2):M[p+1][0]=(subset[p],);M[0][p+1]=(0,subset[1-p])
 assignments={(1,1):1} # mandatory first/last in distinct subsets
 parity=[1,0] # through path has two mandatory vertices
 for _ in range(m-2):
  nxt={}
  for (a,b),v in assignments.items():
   for key in ((a,b),(a^1,b),(a,b^1)):nxt[key]=nxt.get(key,0)+v
  assignments=nxt
  parity=[parity[0]+2*parity[1],parity[1]+2*parity[0]]
 for p in range(2):
  for r in range(2):M[r+1][p+1]=(parity[p^r],assignments.get((1-p,r),0))
 return M

def transfer(c):
 v=[(1,),(0,),(0,)]
 for m in c:
  M=local(m);v=[sum_poly(times(a,b) for a,b in zip(row,v)) for row in M]
 return v

def sum_poly(xs):
 r=(0,)
 for x in xs:r=plus(r,x)
 return r

def graph(c):
 n=sum(c);boundaries=set(it.accumulate(c));owner=[]
 for i,m in enumerate(c):owner.extend([i]*m)
 return [[i!=j and ((i<j and (owner[i]==owner[j] or (j==i+1 and j in boundaries))) or (i>j and owner[i]!=owner[j] and not(i==j+1 and i in boundaries))) for j in range(n)] for i in range(n)]

def brute(c):
 A=graph(c);n=len(A);cycles={};paths=[{},{}]
 # Enumerate each path exactly once by its ordered sequence, and each cycle rooted at its least vertex.
 def walk(seq,mask):
  a=seq[-1]
  if a==n-1:paths[len(seq)%2][mask]=paths[len(seq)%2].get(mask,0)+1
  if len(seq)>=3 and len(seq)%2 and seq[0]==min(seq) and A[a][seq[0]]:cycles[mask]=cycles.get(mask,0)+1
  for b in range(n):
   if not(mask>>b&1) and A[a][b]:walk(seq+(b,),mask|(1<<b))
 for i in range(n):walk((i,),1<<i)
 @lru_cache(None)
 def pack(mask):
  if not mask:return (1,)
  bit=mask&-mask;r=pack(mask^bit)
  for S,mult in cycles.items():
   if S&bit and S&mask==S:r=plus(r,times((0,mult),pack(mask^S)))
  return r
 allmask=(1<<n)-1
 return [pack(allmask)]+[sum_poly(times((mult,),pack(allmask^S)) for S,mult in p.items()) for p in paths]

checked=0
for n in range(1,9):
 for c in comps(n):
  check(transfer(c)==brute(c),(c,transfer(c),brute(c)));checked+=1
print('coefficient and full-state checks',checked,flush=True)
# Exhaust all compositions through order18, without pruning or state deduplication.
matrices={m:[[sum(a*2**i for i,a in enumerate(cell)) for cell in row] for row in local(m)] for m in range(1,19)}
values={1}; levels=[]; witnesses={1:(1,)}
def census(c,v,n):
 if n:
  byorder[n].add(v[0])
  if v[0]<=10000:values.add(v[0]);witnesses.setdefault(v[0],c)
 for m in range(1,19-n):
  M=matrices[m];w=tuple(sum(a*b for a,b in zip(row,v)) for row in M)
  census(c+(m,),w,n+m)
byorder=[set() for _ in range(19)]
census((),(1,0,0),0)
check(len(values)==1374,('census size',len(values)))
check(set.union(*byorder[1:6])=={1,3,5,9,13},'small values')
check(19 not in set.union(*byorder[1:7]),'19 obstruction')
check(transfer((1,3,1))[0]==(1,6),'head extension witness')
# Check producer receipt without importing any producer code.
source=Path(__file__).resolve().parent.parent/'transducer'
receipt=json.loads((source/'receipt.json').read_text())
for row in receipt['levels']:
 n=row['n']
 if n<=18:check(row['attained_up_to_10000']==len({v for vals in byorder[1:n+1] for v in vals if v<=10000}),('receipt level',n))
for m,M in enumerate(receipt['matrices_1_to_5'],1):check(M==matrices[m],('matrix receipt',m))
result={'verdict':'PASS','all_three_polynomial_states_checked_through_order':8,'compositions_checked':checked,'unpruned_compositions_through_18':2**18-1,'values_at_most_10000':len(values),'first_gaps':[i for i in range(1,10001,2) if i not in values][:30],'head_extension_example_polynomial':transfer((1,3,1))[0],'producer_receipt_levels_1_to_18':'MATCH','values_sha256':hashlib.sha256(','.join(map(str,sorted(values))).encode()).hexdigest()}
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
Path(__file__).with_name('values.json').write_text(json.dumps(sorted(values))+'\n')
print(json.dumps(result,indent=2))
# Direct Hamilton-path DP, independent of packing identity.
def hp(A):
 n=len(A);D={ (1<<i,i):1 for i in range(n)}
 for S in range(1,1<<n):
  for v in range(n):
   z=D.get((S,v),0)
   if z:
    for w in range(n):
     if not(S>>w&1) and A[v][w]:D[S|1<<w,w]=D.get((S|1<<w,w),0)+z
 return sum(D.get(((1<<n)-1,v),0) for v in range(n))
hp_checks=0
for n in range(1,10):
 for c in comps(n):
  poly=transfer(c)[0];check(hp(graph(c))==sum(a*2**i for i,a in enumerate(poly)),('HP',c));hp_checks+=1
# All eligible original run words through length10: verify complement equality.
run_checks=0
for L in range(1,11):
 for w in it.product((0,1),repeat=L):
  runs=[];s=0
  for j in range(1,L+1):
   if j==L or w[j]!=w[s]:runs.append((s,j-1,w[s]));s=j
  if any(b-a+1<2 for a,b,_ in runs):continue
  windows=[(a,b+2) for a,b,v in runs if v];n=L+2
  c=[];v=0
  while v<n:
   end=next((b for a,b in windows if a==v),v);c.append(end-v+1);v=end+1
  A=graph(c)
  for i in range(n):
   for j in range(i+1,n):
    expected=(j==i+1 or any(a<=i<j<=b for a,b in windows))
    check(A[i][j]==expected,('run-window complement',w,c,i,j))
  run_checks+=1
result.update({'direct_hamilton_path_checks_through_order9':hp_checks,'original_run_window_adjacency_checks_through_word_length10':run_checks})
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print('additional checks',hp_checks,run_checks)
