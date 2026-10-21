"""Independent recurrence and arithmetic controls for the fixed-H multiplicity packet."""
import json,pathlib,collections,math,itertools,csv,hashlib
P=pathlib.Path(__file__).parent
checks=0
def require(c,msg):
 global checks
 checks+=1
 if not c:raise RuntimeError(msg)
d=json.load(open(P/'extended_census.json'))
a={int(h):{int(n):k for n,k in byn.items()} for h,byn in d['strong_atom_inventory_through_order8'].items()}
c={int(h):{int(n):k for n,k in byn.items()} for h,byn in d['self_converse_strong_inventory_through_order8'].items()}
# Entire ordinary census, generated from first SCC and independent convolution.
f=[collections.Counter() for _ in range(9)];f[0][1]=1
for n in range(1,9):
 f[n]+=f[n-1]
 for h,byn in a.items():
  for m,k in byn.items():
   if m>n:continue
   for x,num in f[n-m].items():f[n][h*x]+=k*num
j=[collections.Counter() for _ in range(9)]
center={(0,1):1,(1,1):1}
for h,byn in c.items():
 for m,k in byn.items():center[m,h]=k
for n in range(9):
 for (m,h),k in center.items():
  if m>n or (n-m)%2:continue
  for x,num in f[(n-m)//2].items():j[n][h*x*x]+=k*num
for n in range(3,9):
 actual={int(h):v for h,v in d['counts_by_order'][str(n)].items()}
 require(set(f[n])==set(actual),f'ordinary support at {n}')
 for h,k in f[n].items():
  require(k==actual[h]['ordinary'],f'ordinary recurrence {n},{h}')
  require(j[n][h]==actual[h]['self_converse'],f'center recurrence {n},{h}')
  require((k+j[n][h])//2==actual[h]['merged'],f'converse quotient {n},{h}')
# Independent brute permutation HP counts, all labelled tournaments through n=5.
for n in range(1,6):
 pairs=list(itertools.combinations(range(n),2));hist=collections.Counter()
 perms=list(itertools.permutations(range(n)))
 for mask in range(1<<len(pairs)):
  rows=[0]*n
  for b,(i,k) in enumerate(pairs):
   if mask>>b&1:rows[k]|=1<<i
   else:rows[i]|=1<<k
  hp=sum(all(rows[p[i]]>>(p[i+1])&1 for i in range(n-1)) for p in perms)
  hist[hp]+=1
 if n>=3:
  expected={int(h):v['labelled'] for h,v in d['counts_by_order'][str(n)].items()}
  require(hist==expected,f'all labelled permutation paths n={n}')
# Witnesses for the arithmetic growth rays, with independent Hamilton DP.
rows8=[{k:int(v) for k,v in r.items()} for r in csv.DictReader(open(P/'n8_classes.csv'))]
witnesses={}
for h in (49,63,343):
 r=next(x for x in rows8 if x['H']==h and x['strong']);out=[r['code']>>(8*i)&255 for i in range(8)]
 dp={(1<<v,v):1 for v in range(8)}
 for mask in range(1,256):
  for end in range(8):
   val=dp.get((mask,end),0)
   for v in range(8):
    if out[end]>>v&1 and not mask>>v&1:dp[mask|(1<<v),v]=dp.get((mask|(1<<v),v),0)+val
 require(sum(dp.get((255,v),0) for v in range(8))==h,f'witness H={h}')
 witnesses[h]=r
# Busch's lower bound makes n<=8 complete for odd H<=73: compare exact cubes.
require(5**8>73**3,'Busch cutoff for H<=73')
# Top-length gap ray arithmetic: permissible divisors have no 7 or21.
def maximal_factor_count(a,b):
 N=3**a*7**b
 ds=[3**i*7**k for i in range(a+1) for k in range(b+1) if i+k and 3**i*7**k not in (7,21)]
 memo={1:0}
 def go(n):
  if n not in memo:memo[n]=max((1+go(n//d) for d in ds if n%d==0),default=-1000)
  return memo[n]
 return go(N)
for a0 in range(2,13):require(maximal_factor_count(a0,1)==a0-1,f'gap ray a={a0}')
for b0 in range(2,13):require(maximal_factor_count(0,b0)==b0//2,f'7-power ray b={b0}')
# The last arithmetic check permits every odd prime-power divisor except7,21;
# its upper bound is attained by actual3,49,63,343 witnesses. It does not assume
# that every other divisor occurs as a Hamiltonian-path count.
# Effective Busch onsets, checked directly against truncated binomial coefficients.
formulas=json.load(open(P/'low_H_formulas.json'))
def generalized_binomial(v,r):return math.prod(v-i for i in range(r))//math.factorial(r)
for H in range(1,75,2):
 words=formulas[str(H)]['core_word_terms'];B=0
 while 5**(B+1)<=H**3:B+=1
 for N in range(B,B+10):
  literal=sum(k*math.comb(N-m+r,r) for r,m,k in words if N>=m)
  polynomial=sum(k*generalized_binomial(N-m+r,r) for r,m,k in words)
  require(literal==polynomial,f'ordinary Busch onset H={H},n={N}')
 for N in range(B+1,B+11):
  literal=polynomial=0
  for x in range(1,math.isqrt(H)+1):
   if H%(x*x):continue
   a0=H//(x*x)
   centers0={0:1,1:1} if a0==1 else c.get(a0,{})
   for size,count in centers0.items():
    if (N-size)%2:continue
    half=(N-size)//2
    for r,m,k in formulas[str(x)]['core_word_terms']:
     if half>=m:literal+=count*k*math.comb(half-m+r,r)
     polynomial+=count*k*generalized_binomial(half-m+r,r)
  require(literal==polynomial,f'center Busch onset H={H},n={N}')
require(a[63]=={8:2},'exact global H63 inventory via Busch')
require(c.get(63,{})=={},'no H63 self-converse atom')
require(all(r['aut']==1 and not r['self_converse'] for r in rows8 if r['H']==63),'H63 atoms automorphism-free and non-SC')
require(sum(a[49].values())==8,'exact global H49 inventory via Busch')
receipt={'status':'PASS','checks':checks,'complete_ordinary_classes_by_n':{n:sum(f[n].values()) for n in range(9)},'complete_self_converse_classes_by_n':{n:sum(j[n].values()) for n in range(9)},'witnesses':witnesses,'scope':'Actual H, ordinary and converse-merged classes kept separate. Full n<=8 census; complete strong atom inventory for odd H<=73, using classical Busch lower bound. No completeness theorem for all odd H except7,21 is assumed.'}
print(json.dumps(receipt,indent=2,sort_keys=True))
