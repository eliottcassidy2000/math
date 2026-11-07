"""Exact standalone numeric14/15/27 controls. Python standard library only."""
import collections,itertools,json,math,pathlib
P=pathlib.Path(__file__).resolve().parent
checks=0
def require(test,msg):
 global checks
 checks+=1
 if not test:raise RuntimeError(msg)
def decode(n,code):
 rows=[0]*n
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  a,b=(i,j) if code>>k&1 else (j,i);rows[a]|=1<<b
 return tuple(rows)
def induce(rows,vertices):
 return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(vertices)) for v in vertices)
def canonical(rows):
 n=len(rows);groups=collections.defaultdict(list)
 for v,r in enumerate(rows):groups[r.bit_count()].append(v)
 best=None;ways=0
 for choices in itertools.product(*(itertools.permutations(groups[d]) for d in sorted(groups))):
  perm=sum(choices,());code=0
  for k,(i,j) in enumerate(itertools.combinations(range(n),2)):code|=((rows[perm[i]]>>perm[j])&1)<<k
  if best is None or code<best:best,ways=code,1
  elif code==best:ways+=1
 return best,ways
def hp(rows):
 n=len(rows);size=1<<n;dp=[[0]*n for _ in range(size)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,size):
  for v in range(n):
   value=dp[mask][v]
   if not value:continue
   available=rows[v]&~mask
   while available:
    b=available&-available;w=b.bit_length()-1;available-=b
    dp[mask|b][w]+=value
 return sum(dp[-1])
def strong(rows):
 n=len(rows)
 for start in range(n):
  seen=1<<start;todo=seen
  while todo:
   bit=todo&-todo;todo-=bit;v=bit.bit_length()-1;new=rows[v]&~seen;seen|=new;todo|=new
  if seen!=(1<<n)-1:return False
 return True
raw=[tuple(map(int,line.split())) for line in (P/'n7_representatives.txt').read_text().splitlines() if line.strip()]
require(len(raw)==456,'inherited 456-class source length')
classes={}
for n in (3,4,5,6,7):
 bank={}
 for rows in raw:
  for vv in itertools.combinations(range(7),n):
   rr=induce(rows,vv);key,aut=canonical(rr)
   bank[key]=aut
 require(len(bank)=={3:2,4:4,5:12,6:56,7:456}[n],f'complete inherited-subset class total n{n}')
 classes[n]=bank
atoms=[]
for n,bank in classes.items():
 for code,aut in bank.items():
  rows=decode(n,code);H=hp(rows)
  if H not in (3,5,9,15,27) or not strong(rows):continue
  converse=tuple(((1<<n)-1)^r^(1<<v) for v,r in enumerate(rows))
  sc=canonical(converse)[0]==code
  atoms.append({'n':n,'H':H,'aut':aut,'self_converse':sc,'code':str(code)})
expected={3:{3:1},5:{4:1},9:{5:2},15:{5:2,6:2},27:{6:2,7:3}}
for h,by_n in expected.items():require(dict(collections.Counter(a['n'] for a in atoms if a['H']==h))==by_n,f'global atom inventory H{h}')
require(5**7>27**3,'Busch excludes all strong n>=8 at H<=27')
# Every atom inventory is complete using the inherited n7 atlas and cited Busch bound.
def core_words(target):
 if target==1:yield ();return
 for a in atoms:
  if target%a['H']==0:
   for rest in core_words(target//a['H']):yield (a,)+rest
def count(n,h):
 return sum(math.comb(n-sum(a['n'] for a in w)+len(w),len(w)) for w in core_words(h) if n>=sum(a['n'] for a in w))
def self_count(n,h):
 total=0
 for x in range(1,math.isqrt(h)+1):
  if h%(x*x):continue
  c=h//(x*x)
  centers=[(0,1),(1,1)] if c==1 else [(a['n'],1) for a in atoms if a['H']==c and a['self_converse']]
  for size,mult in centers:
   if n>=size and (n-size)%2==0:total+=mult*count((n-size)//2,x)
 return total
for n in range(5,81):require(count(n,15)==(n-4)*(n-3),f'N15 n{n}')
for n in range(6,81):require(count(n,27)==n*(n-5)*(n-4)//6,f'N27 n{n}')
for n in range(7,81):
 require(self_count(n,15)==2*(n%2),f'S15 n{n}')
 require(self_count(n,27)==((n-5)//2 if n%2 else 0),f'S27 n{n}')
strata={}
for h in (15,27):
 n=14;ordinary_by_degree=collections.Counter();terms=collections.Counter()
 for w in core_words(h):
  size=sum(a['n'] for a in w);r=len(w)
  if size>n:continue
  ways=math.comb(n-size+r,r);aut=math.prod(a['aut'] for a in w)
  ordinary_by_degree[h//aut]+=ways;terms[(r,size)]+=1
 require(self_count(n,h)==0,f'no self-converse H{h} n14')
 require(all(v%2==0 for v in ordinary_by_degree.values()),'converse pair degree divisibility')
 merged_by_degree={str(k):v//2 for k,v in sorted(ordinary_by_degree.items())}
 strata[str(h)]={'ordinary':count(n,h),'self_converse':0,'merged':count(n,h)//2,'paired_black_degree_histogram':merged_by_degree,'ordered_core_terms':[{'components':r,'order':s,'multiplicity':c} for (r,s),c in sorted(terms.items())]}
require(strata['15']['ordinary']==110 and strata['27']['ordinary']==210,'n14 counts')
require(strata['15']['paired_black_degree_histogram'].get('1',0)==0,'H15 no leaves')
require(strata['27']['paired_black_degree_histogram']['1']==28,'H27 exact28 leaves')
require(math.comb(8,3)==56,'H27 three-triangle stack count n14')
require(2*math.comb(9,2)==72,'H15 prime stack count n14')
# Existing strong order14 witness, independently recounted without oracle code.
code=292187853813044375413825338
rows=decode(14,code);value=hp(rows)
require(strong(rows),'order14 witness strong')
require(value==15**6==11390625,'order14 witness count')
# Direct small minimum witnesses and the inherited maximum through five.
require(hp(decode(5,330))==15,'order5 H15 witness')
require(hp(decode(6,10514))==27,'order6 H27 witness')
require(max(hp(decode(4,c)) for c in classes[4])<15,'m15 lower bound')
require(max(hp(decode(5,c)) for c in classes[5])<27,'m27 lower bound')
result={'status':'PASS','active_checks':checks,'scope':'Complete classes inherited via n7 source; exact subset and isomorphism recalculation through n7, cited Busch cutoff, SCC count proofs and current explicit order14 witness. Not an n14 class census or a spectrum-completeness proof.','atom_inventory':atoms,'n14':strata,'minimum_orders':{'15':5,'27':6},'order14_atom':{'code':str(code),'H':value,'power_of_15':6,'strong':True,'rate_upper_bound':'7/3'},'H27_n14_leaf_fraction':'4/15','H15_n14_prime_stack_merged_nodes':36}
print(json.dumps(result,indent=2,sort_keys=True))
