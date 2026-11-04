import sys,json,itertools,math
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'collision'))
from verify import polynomial,hp,extend,tournament

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def opposite(a):return [((1<<len(a))-1)^(1<<u)^r for u,r in enumerate(a)]
def tensor(a):
 n=len(a);paths=[]
 for k in range(1,n-2):
  for p in itertools.permutations(range(n),k):
   if all(a[u]>>v&1 for u,v in zip(p,p[1:])):paths.append((sum(1<<u for u in p),p[0],p[-1],k))
 out={}
 for A,u,_,ka in paths:
  for B,v,w,kb in paths:
   if A&B or (ka+kb)%2==0 or ka+kb>n-3:continue
   key=((1<<u)|(1<<v),1<<w);out[key]=1-out.get(key,0)
 return out
controls=[]
for n,codes in [(2,range(2)),(4,range(64)),(6,[0,1,2,3,17,31,97,193,777,1234,4095,7171,16383,27171,32767])]:
 for code in codes:
  a=tournament(n,code);p=polynomial(a);q=polynomial(opposite(a));K=tensor(a);tested=0
  for x in range(1<<n):
   if x.bit_count()!=2:continue
   for w in range(n):
    y=1<<w
    for m in [x|(y<<n),y|(x<<n)]:
     c=p.get(m,0)-q.get(m,0);need(c%2==0,'integer evenness');need((c//2)%2==K.get((x,y),0),'independent literal-word cutoff');tested+=1
  controls.append({'order':n,'code':code,'mixed_coordinates':tested})
d=json.loads((BASE/'prime12/cubic.json').read_text());a=d['rows'];q=opposite(a);cache={}
def delta(x,y):
 if (x,y) not in cache:cache[x,y]=hp(extend(a,x,y))-hp(extend(q,x,y))
 return cache[x,y]
coeff={(x,y):v for x,y,v in d['coefficients']};checked=[]
for x,y in [(1,3),(3,4),(3,3)]:
 total=0
 for i in range(1<<12):
  if i&x!=i:continue
  for j in range(1<<12):
   if j&y!=j:continue
   total+=(-1)**((x^i).bit_count()+(y^j).bit_count())*delta(i,j)
 need(total==coeff.get((x,y),0),'independent P literal coefficient');checked.append([x,y,total])
print(json.dumps({'status':'PASS','small_controls':controls,'mixed_coordinates_checked':sum(c['mixed_coordinates'] for c in controls),'P_coefficients':checked,'P_literal_queries':len(cache)},indent=2,sort_keys=True))
