import pathlib,hashlib,json,itertools as it,math
P=pathlib.Path(__file__).resolve().parent.parent
S=next((p for p in [P/'session43-native-branch-width',P/'native-branch-width'] if (p/'MANIFEST.json').is_file()),None)
if S is None:raise RuntimeError('source packet absent')
checks=0
def need(ok,msg):
 global checks
 if not ok:raise RuntimeError(msg)
 checks+=1
for f in json.loads((S/'MANIFEST.json').read_text())['files']:
 b=(S/f['path']).read_bytes();need(len(b)==f['bytes'] and hashlib.sha256(b).hexdigest()==f['sha256'],'source hash')
def arc(p,t,E,x,y):
 i,u=divmod(x,p);j,v=divmod(y,p)
 if i==j:return 0<(v-u)%p<=(p-1)//2
 if 0<(j-i)%t<=(t-1)//2:return (v-u-1)%p in E
 return (u-v-1)%p not in E
host_records=[]
for p,t,E in [(3,3,{0,1}),(3,5,{0,1}),(3,9,{0,1}),(5,3,{0,2,4}),(5,5,{1}),(7,3,{0,1,3})]:
 n=p*t
 for x in range(n):
  need(not arc(p,t,E,x,x),'no loop')
  need(sum(arc(p,t,E,x,y) for y in range(n))==(n-1)//2,'regular')
  for y in range(n):
   if x!=y:need(arc(p,t,E,x,y)!=arc(p,t,E,y,x),'orientation complement')
 counts=[[0]*p for _ in range(p)];accepted=set()
 for a in it.product(range(p),repeat=t):
  path=[i*p+(a[i]+k)%p for i in range(t) for k in range(p)]
  need(len(set(path))==n,'literal spanning permutation')
  actual=all(arc(p,t,E,x,y) for x,y in zip(path,path[1:]))
  expected=all((v-u)%p in E for u,v in zip(a,a[1:]))
  need(actual==expected,'actual path iff differences')
  if actual:
   accepted.add(tuple(path));counts[a[0]][a[-1]]+=1
   need(all(not arc(p,t,E,y,x) for x,y in zip(path,path[1:])),'reversed path in converse')
 need(len(accepted)==p*len(E)**(t-1),'distinct literal count')
 # Independent endpoint count by the integer coefficient distribution of a product.
 coeff={0:1}
 for _ in range(t-1):
  new={}
  for k,c in coeff.items():
   for e in E:new[k+e]=new.get(k+e,0)+c
  coeff=new
 endpoint=[sum(c for k,c in coeff.items() if k%p==d) for d in range(p)]
 for a in range(p):need(counts[a]==[endpoint[(b-a)%p] for b in range(p)],'endpoint normalization')
 host_records.append({'p':p,'t':t,'E':sorted(E),'paths':len(accepted),'endpoint':endpoint})
# All affine subspaces of F_p^3: points, line parametrizations, plane equations, full space.
coset_records=[]
for p in [3,5]:
 V=list(it.product(range(p),repeat=3));zero=(0,0,0)
 directions=[]
 for v in V:
  if v==zero:continue
  lead=next(x for x in v if x)
  if lead==1:directions.append(v)
 cosets={frozenset([v]) for v in V};cosets.add(frozenset(V))
 for d in directions:
  for v in V:cosets.add(frozenset(tuple((v[j]+k*d[j])%p for j in range(3)) for k in range(p)))
  for b in range(p):cosets.add(frozenset(v for v in V if sum(v[j]*d[j] for j in range(3))%p==b))
 for mask in range(1,2**p-1):
  E={e for e in range(p) if mask>>e&1}
  good={v for v in V if (v[1]-v[0])%p in E and (v[2]-v[1])%p in E}
  contained=[C for C in cosets if C<=good]
  need(max(map(len,contained))==p,'maximum arbitrary affine coset')
  need(len(good)//p==len(E)**2,'cover cardinality lower bound')
  fibres={}
  for v in good:fibres.setdefault(((v[1]-v[0])%p,(v[2]-v[1])%p),set()).add(v)
  need(all(frozenset(F) in cosets for F in fibres.values()) and len(fibres)==len(E)**2,'attaining cover')
 coset_records.append({'p':p,'all_affine_cosets':len(cosets),'proper_nonempty_E_tested':2**p-2})
# Closed form against binomial coefficients, including full (not only even) clock.
for r in range(101):
 for d in range(3):
  actual=sum(math.comb(r,k) for k in range(r+1) if k%3==d)
  num=2**r+[2,1,-1,-2,-1,1][(r-2*d)%6]
  need(num%3==0 and num//3==actual,'period six endpoint clock')
print(json.dumps({'status':'PASS','checks':checks,'hosts':host_records,'affine_tests':coset_records,'source_manifest_sha256':hashlib.sha256((S/'MANIFEST.json').read_bytes()).hexdigest()},sort_keys=True,indent=2))
