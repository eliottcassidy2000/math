import itertools,json,pathlib,hashlib
BASE=pathlib.Path('/workspace/shared/session187/collision')
R=json.loads((BASE/'receipt.json').read_text())
def need(x,s):
 if not x:raise ValueError(s)
def count(rows):
 n=len(rows);size=1<<n;inc=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]; d=[[0]*n for _ in range(size)]
 for v in range(n):d[1<<v][v]=1
 for s in range(1,size):
  b=s
  while b:
   bit=b&-b;b-=bit;v=bit.bit_length()-1;old=s^bit
   if old:
    prev=old&inc[v];z=0
    while prev:
     w=prev&-prev;prev-=w;z+=d[old][w.bit_length()-1]
    d[s][v]=z
 return sum(d[-1])
def port(rows,a,c):
 n=len(rows);out=list(rows)+[1<<(n+1),0]
 for v in range(n):
  for p,m in ((n,a),(n+1,c)):
   if m&(1<<v):out[v]|=1<<p
   else:out[p]|=1<<v
 return out
def strong(rows):
 n=len(rows)
 for v in range(n):
  seen=1<<v
  while True:
   nxt=seen
   for u in range(n):
    if seen>>u&1:nxt|=rows[u]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True
T=R['T'];U=R['U'];n=len(T);O=[((1<<n)-1)^(1<<v)^T[v] for v in range(n)]
need(strong(T) and strong(U),'strong controls')
def homogeneous(rows):
 n=len(rows);return {tuple((u,v)) for u,v in itertools.combinations(range(n),2) if all(((rows[u]>>w)&1)==((rows[v]>>w)&1) for w in range(n) if w not in (u,v))}
need(not homogeneous(T)&homogeneous(U),'no common pair')
tables=[]
for rows in (T,U,O):
 table=[count(port(rows,s&63,s>>6)) for s in range(4096)]; tables.append(table)
contrasts=[]
for k,key in ((1,'arbitrary_difference'),(2,'global_converse_difference')):
 a=[x-y for x,y in zip(tables[0],tables[k])]
 for bit in range(12):
  for s in range(4096):
   if s>>bit&1:a[s]-=a[s^(1<<bit)]
 exact={str(s):v for s,v in enumerate(a) if v}
 need(exact==R[key],'all coefficients '+key)
 contrasts.append({'key':key,'degree':max(s.bit_count() for s,v in enumerate(a) if v),'nonzero':len(exact)})
# Independent literal color assignments, Hamilton words on every subset.
colorchecks=0
for r in range(1,6):
 for i in range(r):
  rows=port([sum(1<<v for v in range(u+1,r)) for u in range(r)],0,1<<i);n=len(rows)
  h={0:1}
  for mask in range(1,1<<n):
   ids=[v for v in range(n) if mask>>v&1]
   h[mask]=sum(all(rows[u]>>v&1 for u,v in zip(p,p[1:])) for p in itertools.permutations(ids))
  for t in range(1,4):
   z=0
   for colors in itertools.product(range(t),repeat=n):
    masks=[0]*t
    for v,col in enumerate(colors):masks[col]|=1<<v
    weight=1
    for mask in masks:weight*=h[mask]
    z+=weight
   need(z==t**(r+2)+t**(r-i+1)*((t+1)**i-(t-1)**i),'literal color formula');colorchecks+=1
result={'status':'PASS','independent_full_truth_queries':12288,'literal_coloring_controls':colorchecks,'strong_T_U':True,'no_common_homogeneous_pair':True,'contrasts':contrasts,'source_hashes':{f:hashlib.sha256((BASE/f).read_bytes()).hexdigest() for f in ('PROOF.md','verify.py','receipt.json')}}
print(json.dumps(result,sort_keys=True,indent=2))
