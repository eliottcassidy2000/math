import itertools as I,json,math,hashlib,collections as C,pathlib
D=pathlib.Path(__file__).resolve().parent
F=D/'frozen'
def ck(x,s):
 if not x: raise RuntimeError(s)
V=set(range(7));Q={1,2,4};E=sorted((x,y) for x in V for y in V if (y-x)%7 in Q);idx={e:i for i,e in enumerate(E)}
def mask(es):return sum(1<<idx[e] for e in es)
paths=[]
def dfs(p):
 if len(p)==7:paths.append(tuple(zip(p,p[1:])));return
 for y in V-set(p):
  if not p or (y-p[-1])%7 in Q:dfs(p+[y])
dfs([]);ck(len(paths)==189,'paths')
cover=[sum(1<<j for j,p in enumerate(paths) if e in p) for e in E];full=(1<<len(paths))-1
blocks=[];counts={}
for k in (5,6):
 n=0
 for inds in I.combinations(range(21),k):
  n+=1;c=0
  for i in inds:c|=cover[i]
  if c==full:
   ck(k==6,'minimum');blocks.append(sum(1<<i for i in inds))
 counts[k]=n
ck(len(blocks)==63,'census');bs=set(blocks)
def arcs(b):return {e for i,e in enumerate(E) if b>>i&1}
def action(b,a,c,r=0):
 return mask((((a*(-y if r else x)+c)%7,(a*(-x if r else y)+c)%7) for x,y in arcs(b)))
def classify(b):
 R=set(E)-arcs(b);src=V-{y for x,y in R};snk=V-{x for x,y in R}
 options=[]
 if len(src)==2:options.append(('two-source',src))
 if len(snk)==2:options.append(('two-sink',snk))
 if src&snk: options.append(('star',src&snk))
 for c in V:
  if sum({y for x,y in R if x==z}=={c} for z in V)==3:options.append(('common-out',{c}))
  if sum({x for x,y in R if y==z}=={c} for z in V)==3:options.append(('common-in',{c}))
 ck(len(options)==1,'intrinsic unique classification');return options[0]
def maps(b):
 t,S=classify(b)
 if t.startswith('two'):
  u,v=next((u,v) for u,v in E if {u,v}==S);k={1:0,2:1,4:2}[(v-u)%7]
  if t=='two-source':return frozenset(S),(9*u+16**k)%63
  return frozenset(V-{u,v,(2*v-u)%7}),(9*v-16**k)%63
 c=next(iter(S))
 if t=='star':return frozenset(V-{c}),9*c
 if t=='common-out':return frozenset({c}|{y for x,y in E if x==c}),(9*c+21)%63
 return frozenset({c}|{x for x,y in E if y==c}),(9*c+42)%63
M={b:maps(b) for b in bs};even={frozenset(s) for k in (2,4,6) for s in I.combinations(V,k)}
ck({s for s,e in M.values()}==even,'subset bijection');ck({e for s,e in M.values()}==set(range(63)),'exponent bijection')
ops=list(I.product((1,2,4),range(7),(0,1)));perm=[];eq=0
for a,c,r in ops:
 k={1:0,2:1,4:2}[a];p=[]
 for b in sorted(bs):
  z=action(b,a,c,r);ck(z in bs,'closed');p.append(z)
  s,e=M[b];ss,ee=M[z];ck(ee==(9*c+16**k*(-e if r else e))%63,'G42 equation');eq+=1
  if not r:ck(ss==frozenset((a*x+c)%7 for x in s),'G21 subset equation')
 perm.append(tuple(p))
ck(len(set(perm))==42,'faithful group42')
orbits=[];left=bs.copy()
while left:
 b=min(left);o={action(b,*g) for g in ops};left-=o
 st=[g for g in ops if action(b,*g)==b];ph={9//math.gcd(M[z][1],9) for z in o}
 ck(len(o)*len(st)==42,'stabilizer');orbits.append({'size':len(o),'stabilizer_order':len(st),'phase_orders':sorted(ph)})
ck(sorted(o['size'] for o in orbits)==[7,14,42],'orbit types')
fixed=[b for b in bs if action(b,1,0,1)==b];ck(len(fixed)==1 and classify(fixed[0])==('star',{0}),'fixed blocker')
ck(sum(frozenset(-x%7 for x in s)==s for s in even)==7,'ordinary involution')
R=json.loads((F/'result.json').read_text());P=json.loads((F/'converse.json').read_text())
ck({x['blocker'] for x in R['bijection']}==bs,'producer census')
for x in R['bijection']:ck((frozenset(x['even_subset']),x['type'])==(M[x['blocker']][0],classify(x['blocker'])[0]),'producer subset rows')
for x in P['map']:ck((x['exponent'],x['type'])==(M[x['blocker']][1],classify(x['blocker'])[0]),'producer exponent rows')
# Hostile incorrect converse lift, natural incidence, XOR, and metric.
wrong=sum(M[action(b,1,0,1)][0]!=frozenset(-x%7 for x in M[b][0]) for b in bs)
ck(wrong>0,'reject naive subset converse')
bdry={}
for b in bs:
 s=set()
 for x,y in arcs(b):s.symmetric_difference_update({x,y})
 bdry[b]=frozenset(s)
ck(C.Counter(C.Counter(bdry.values()).values())=={3:21},'boundary loss')
xor=next((a,b,a^b) for a,b in I.combinations(sorted(bs),2) if a^b not in bs);metric=next((a,b,(a^b).bit_count()) for a,b in I.combinations(sorted(bs),2) if (a^b).bit_count()>6)
# Fourier implementation independent polynomial long multiplication.
def mul(a,b):
 z=0
 for k in range(3):
  if b>>k&1:z^=a<<k
 for k in range(5,2,-1):
  if z>>k&1:z^=11<<(k-3)
 return z
powers=[1]
for k in range(6):powers.append(mul(powers[-1],2))
def four(s):
 a=b=0
 for x in s:a^=powers[x];b^=powers[-x%7]
 return a,b
ck(len({four(s) for s in even|{frozenset()}})==64,'Fourier bijection')
for x in R['bijection']:ck(list(four(M[x['blocker']][0]))==x['fourier'],'Fourier rows')
for b in bs:
 a,c=four(M[b][0]);t=classify(b)[0]
 good={'star':mul(a,c)==1,'common-out':a!=0 and c==0,'common-in':a==0 and c!=0,'two-sink':mul(a,c) in {powers[k] for k in (1,2,4)},'two-source':mul(a,c) in {powers[k] for k in (3,5,6)}}[t]
 ck(good,'Fourier type strata')
print(json.dumps({'status':'PASS','enumeration_method':'independent DFS Hamiltonian paths and arc-cover bitsets; no producer imports','paths':len(paths),'candidate_counts':counts,'blockers':len(bs),'G42_equations':eq,'orbits':orbits,'fixed_blockers':fixed,'hostile_naive_converse_mismatches':wrong,'xor_witness':xor,'metric_witness':metric,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(F.iterdir()) if p.is_file()}},sort_keys=True,indent=2))
