import itertools as I,json,hashlib,pathlib,subprocess,sys
P=pathlib.Path(__file__).resolve().parent; R=P.parent
def ck(x):
 if not x: raise RuntimeError('audit failure')
def edges(n,cycle):
 f={tuple(sorted((i,i+1))) for i in range(n-1)}
 if cycle:f.add((0,n-1))
 return [e for e in I.combinations(range(n),2) if e not in f]
result={}
for n in (4,6):
 for cyc in (False,True):
  es=edges(n,cyc);seen=set()
  for bits in I.product((0,1),repeat=len(es)):
   d=[0]*n
   for e,b in zip(es,bits):d[e[b]]+=1
   seen.add(tuple(d))
  good=[]
  for H in I.combinations(range(n),n//2+(not cyc)):
   d=tuple(n//2-2+(v in H) for v in range(n));ok=d in seen
   expect=not(cyc and n==4) or all(len(set(H)&set(e))==1 for e in es)
   ck(ok==expect)
   if ok:good.append(H)
  result[f'orientations_{n}_{cyc}']={'universe':2**len(es),'high_sets':len(good)}
for n in range(3,11):
 for cyc in (False,True):
  es=edges(n,cyc)
  if n%2: ds=[tuple((n-3)//2+(v==w and not cyc) for v in range(n)) for w in range(n)]
  else:ds=[tuple(n//2-2+(v in H) for v in range(n)) for H in I.combinations(range(n),n//2+(not cyc))]
  good=0
  for d in ds:
   ck(sum(d)==len(es))
   feasible=True
   for mask in range(1<<n):
    S={v for v in range(n) if mask>>v&1}
    if sum(d[v] for v in S)<sum(a in S and b in S for a,b in es): feasible=False;break
   if feasible:good+=1
   if n>=5:ck(feasible)
  result[f'cut_{n}_{cyc}']={'patterns':len(ds),'feasible':good}
# Independent enumeration of every modular grid at bad primes, using only lines.
lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
for m in (2,3):
 count=0
 for v in I.product(range(m),repeat=9):
  ss=[sum(v[i] for i in l)%m for l in lines]
  if len(set(ss))==1:
   count+=1;c,a,b=v[4],v[0],v[2]
   ck(v==tuple(x%m for x in (a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a)))
 ck(count==m**3);result[f'all_grids_mod_{m}']=count
pins={}
for folder,names in [('carry-lattice',['PROOF.md','verify.py','certificate.json']),('residual-scores',['THEOREM.md','check.py','README.md'])]:
 for name in names:
  rel=folder+'/'+name;pins[rel]=hashlib.sha256((R/rel).read_bytes()).hexdigest()
result['source_sha256']=pins
print(json.dumps(result,indent=2,sort_keys=True))
