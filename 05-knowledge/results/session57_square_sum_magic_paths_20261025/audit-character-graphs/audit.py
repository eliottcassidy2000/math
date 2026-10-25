import hashlib,json,itertools,pathlib,collections,sys
SRC=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).resolve().parent.parent/'character-graphs'
OUT=pathlib.Path(__file__).parent

def require(v,msg):
 if not v: raise RuntimeError(msg)

def rem(a,b,p):
 a=list(a)
 while len(a)>=len(b):
  k=a[-1]*pow(b[-1],-1,p)%p; d=len(a)-len(b)
  for j in range(len(b)):a[d+j]=(a[d+j]-k*b[j])%p
  a.pop()
 return a
class Field:
 def __init__(self,p,f):
  self.p=p;self.f=f;self.n=len(f)-1;self.q=p**self.n
  require(f[-1]==1,'monic')
  for d in range(1,self.n//2+1):
   for coeff in itertools.product(range(p),repeat=d):
    require(any(rem(f,list(coeff)+[1],p)),f'reducible {p} {f}')
  self.d=[tuple((x//p**i)%p for i in range(self.n)) for x in range(self.q)]
  self.sq={self.mul(x,x) for x in range(self.q)}
  require(len(self.sq)==(self.q+1)//2,'square cardinality')
 def enc(self,cs):return sum((c%self.p)*self.p**i for i,c in enumerate(cs))
 def add(self,*xs):return self.enc(map(sum,zip(*(self.d[x] for x in xs))))
 def neg(self,x):return self.enc([-c for c in self.d[x]])
 def mul(self,x,y):
  cs=[0]*(2*self.n-1)
  for i,a in enumerate(self.d[x]):
   for j,b in enumerate(self.d[y]):cs[i+j]+=a*b
  return self.enc(rem(cs,self.f,self.p))
L=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def matrix(F,c,u,v):
 a=F.add;U=F.neg(u);V=F.neg(v)
 return [a(c,u),a(c,U,V),a(c,v),a(c,U,v),c,a(c,u,V),a(c,V),a(c,u,v),a(c,U)]
def valid(F,m):return len(set(m))==9 and set(m)<=F.sq and len({F.add(*(m[i] for i in l)) for l in L})==1
def adj(F,V,zero):return [[i!=j and F.add(x,y) in F.sq and (zero or F.add(x,y)!=0) for j,y in enumerate(V)] for i,x in enumerate(V)]
def hp(A):
 dp=collections.defaultdict(int)
 for i in range(9):dp[(1<<i,i)]=1
 for mask in range(1,512):
  for i in range(9):
   v=dp[(mask,i)]
   if v:
    for j in range(9):
     if not (mask>>j)&1 and A[i][j]:dp[(mask|(1<<j),j)]+=v
 return sum(dp[(511,i)] for i in range(9))//2
pins={}
expected_pins=json.loads((OUT/'SOURCE_PINS.json').read_text())
for line in (SRC/'SHA256SUMS').read_text().splitlines():
 sha,name=line.split(); actual=hashlib.sha256((SRC/name).read_bytes()).hexdigest();require(sha==actual,'source pin '+name);pins[name]=sha
require(pins==expected_pins,'frozen source pins')
D=json.loads((SRC/'results.json').read_text());report={'pinned_source_files':len(pins),'fields':[], 'affine_planes':[]}
# All normalized grids checked with independently coded field operations.
for r in D['fields']:
 F=Field(r['p'],r['modulus']);count=0
 for c in (0,1):
  for u in ([1] if c==0 else range(F.q)):
   for v in range(F.q):
    m=matrix(F,c,u,v)
    if valid(F,m):
     count+=1
     require(all(F.add(m[i],m[8-i])==F.add(c,c) for i in range(4)), 'opposite identity')
     require(all(F.add(*(m[i] for i in l))==F.add(c,c,c) for l in L),'total identity')
 require(count==r['normalized_magic_count'],f'count {F.q}')
 for g in r['spectral_hp_groups']+([r['first_nine_square_graph']] if 'first_nine_square_graph' in r else []):
  if 'matrix' in g:require(valid(F,g['matrix']),'witness validity')
  for zero in (True,False):
   A=adj(F,g['vertices'],zero);s=g['graphs'][str(zero)]
   require(sum(map(sum,A))//2==s['edges'],'edges')
   require(sorted(map(sum,A))==s['degrees'],'degrees')
   require(hp(A)==s['hp_unoriented'],'HP DP')
 report['fields'].append({'q':F.q,'normalized_magic_count':count})
 print('PASS',F.q,flush=True)
# Exhaustive paths and nonexistence F17.
r=next(r for r in D['fields'] if r['q']==17);F=Field(17,r['modulus']);V=sorted(F.sq)
require(r['normalized_magic_count']==0,'F17 Parker')
report['F17']={}
for zero,want in ((True,1152),(False,112)):
 A=adj(F,V,zero);n=sum(all(A[P[i]][P[i+1]] for i in range(8)) for P in itertools.permutations(range(9)))//2
 require(n==want,'F17 brute HP');report['F17'][str(zero)]=n
# Printed F59 witness, not only representative in census.
F=Field(59,[0,1]);m=[15,35,12,57,1,4,49,26,46]
require(valid(F,m),'printed F59 witness')
for zero in (True,False):
 A=adj(F,m,zero);require(sum(A[3])==0 and hp(A)==0,'F59 isolated57')
report['F59']={'isolated_vertex':57,'hp':0,'magic_total':F.add(*m[:3])}
for q in (81,729,2187):
 r=next(r for r in D['fields'] if r['q']==q) if q==81 else json.loads((SRC/f'witness-{q}.json').read_text())
 g=r['spectral_hp_groups'][0] if q==81 else r;F=Field(r['p'],r['modulus']);V=g['vertices']
 require(valid(F,g['matrix']),'extension magic')
 require((F.neg(1) in F.sq)==(F.n%2==0),'minusone character')
 A=adj(F,V,True);require(sum(map(sum,A))==(72 if F.n%2==0 else 0),'collapse')
 if F.n%2:require(0 not in V,'odd nozero')
 W={F.add(x,F.neg(V[0])) for x in V};require(len(W)==9 and all(F.add(x,y) in W for x in W for y in W),'affineplane closure')
 H={t for t in itertools.combinations(range(9),3) if F.add(*(V[i] for i in t))==0};require(len(H)==12,'H0 lines')
 n=sum(all(tuple(sorted(P[i] for i in l)) in H for l in L) for P in itertools.permutations(range(9)))
 require(n==432,'magic arrangements')
 for zero in (True,False):require(hp(adj(F,V,zero))==g['graphs'][str(zero)]['hp_unoriented'],'extension HP')
 report['affine_planes'].append({'q':q,'H0_edges':12,'magic_arrangements':n,'all_arrangements':362880,'square_sum_edges':sum(map(sum,A))//2})
 print('PASS plane',q,flush=True)
report['status']='PASS';report['total_normalized']=sum(x['normalized_magic_count'] for x in report['fields'])
(OUT/('optimized.json' if not __debug__ else 'normal.json')).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True))
