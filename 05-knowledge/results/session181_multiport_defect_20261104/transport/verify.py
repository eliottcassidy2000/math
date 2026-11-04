import itertools,json,sys
from pathlib import Path
from common import response,strong
OUT=Path(__file__).parent

def check(x,msg):
 if not x:raise RuntimeError(msg)
def add(*ps):
 q={}
 for p in ps:
  for m,v in p.items():q[m]=q.get(m,0)+v
 return q
def mul(p,q):
 r={}
 for m,v in p.items():
  for n,w in q.items():
   if not m&n:r[m|n]=r.get(m|n,0)+v*w
 return r
one={0:1}
def matrix(rows):
 n=len(rows);g=[[{} for j in rows]for i in rows]
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(rows[a]>>b&1 for a,b in zip(p,p[1:])):
    m=sum(1<<x for x in p);d=g[p[0]][p[-1]];d[m]=d.get(m,0)+1
 return g
def birth(r,a):
 n=len(r);return [v|((1<<n) if a>>i&1 else 0) for i,v in enumerate(r)]+[((1<<n)-1)^a]
def reduced(rows,w):
 vs=[i for i in range(len(rows))if i!=w]
 return [sum(1<<j for j,b in enumerate(vs)if rows[a]>>b&1)for a in vs],vs
counts={'defect':0,'birth':0,'dominant_prefix':0}
for n in range(1,5):
 pairs=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  r=[0]*n
  for k,(i,j)in enumerate(pairs):
   a,b=(i,j)if bits>>k&1 else(j,i);r[a]|=1<<b
  h,fs=response(r)
  for w in range(n):
   u,vs=reduced(r,w);m=len(u);full=(1<<m)-1;g=matrix(u);I=[i for i,v in enumerate(vs)if r[v]>>w&1]
   e=add(*(g[i][j]for i in range(m)for j in I));s=add(*(g[i][j]for i in I for j in range(m)))
   defect=mul(add(one,e),s).get(full,0)
   check(defect==fs[w]-h,('defect',r,w));counts['defect']+=1
   for cut in range(1<<n):
    aa=[i for i,v in enumerate(vs)if cut>>v&1];bb=[i for i in range(m)if i not in aa];eps=0 if cut>>w&1 else 1
    X=add(*(g[i][j]for i in range(m)for j in aa));Y=add(*(g[i][j]for i in bb for j in range(m)))
    C=add(*(g[i][j]for i in bb for j in I));D=add(*(g[i][j]for i in I for j in aa))
    prediction=add(mul(mul(add(one,X),add(C,{0:eps})),s),mul(mul(add(one,e),add(D,{0:eps})),add(one,Y))).get(full,0)
    nr=birth(r,cut);nh,nf=response(nr);check(prediction==nf[w]-nh,('birth',r,w,cut));counts['birth']+=1
   A=response(u)[0] if u else 1;SI=s.get(full,0);f=fs[w]
   nr=r
   for t in range(1,4):
    nr=birth(nr,1<<w);nh,nf=response(nr)
    check(nh==h+SI+((1<<t)-1)*(f-A),('top_h',r,w,t,nh))
    check(nf[w]==A+(1<<t)*(f-A),('top_f',r,w,t))
    check(nf[w]-nh==defect-SI,('top_b',r,w,t))
    if strong(r) and n>1:check(strong(nr),('strong',r,w,t))
    counts['dominant_prefix']+=1
print(json.dumps(counts,indent=2));(OUT/'verification.json').write_text(json.dumps({'status':'PASS','checks':counts,'scope':'Every labeled tournament orders 1 through 4, every mark and birth cut; dominant-prefix lengths1,2,3. Direct path-polynomial enumeration versus subset DP.'},indent=2)+'\n')
