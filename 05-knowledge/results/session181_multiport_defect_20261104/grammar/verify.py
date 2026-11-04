"""Independent exact checks, standard-library Python; assertions never used."""
import itertools,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parent

def need(x,s):
 if not x:raise RuntimeError(s)
def add(*polys):
 d={}
 for p in polys:
  for k,v in p.items():d[k]=d.get(k,0)+v
 return {k:v for k,v in d.items() if v}
def mul(p,q):
 d={}
 for k,v in p.items():
  for l,u in q.items():
   if not k&l:d[k|l]=d.get(k|l,0)+v*u
 return d
def matrix(a):
 n=len(a);g=[[{} for j in a] for i in a]
 # Independent enumeration of all ordered distinct tuples.
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(a[i]>>j&1 for i,j in zip(p,p[1:])):
    s=sum(1<<i for i in p);q=g[p[0]][p[-1]];q[s]=q.get(s,0)+1
 return g
def direct(a,w):
 n=len(a);h=f=bad=0
 for p in itertools.permutations(range(n)):
  if all(a[i]>>j&1 for i,j in zip(p,p[1:])):h+=1
  t=p.index(w)+1
  if all(a[p[k]]>>p[k+1]&1 for k in range(n-1) if k!=t-1):
   f+=1
   if t<n and a[p[t]]>>w&1:bad+=1
 return h,f,bad
def birth(a,c):
 n=len(a);return tuple(r|(((c>>i)&1)<<n) for i,r in enumerate(a))+(((1<<n)-1)^c,)
def check(a,w,c):
 n=len(a);full=(1<<n)-1;g=matrix(a);one={0:1}
 e=[add(*(g[i][j] for i in range(n))) for j in range(n)]
 s=[add(*row) for row in g]
 h,f,b=direct(a,w);q=add(*(s[v] for v in range(n) if a[v]>>w&1))
 need(mul(e[w],q).get(full,0)==f-h==b,'bad joins identity')
 L=add(*(e[i] for i in range(n) if c>>i&1));R=add(*(s[i] for i in range(n) if not c>>i&1))
 rw=add(*(g[j][w] for j in range(n) if not c>>j&1))
 K=add(*(g[v][i] for v in range(n) if a[v]>>w&1 for i in range(n) if c>>i&1))
 if not c>>w&1:K=add(K,one)
 predh=mul(add(one,L),add(one,R)).get(full,0)
 predb=add(mul(mul(add(one,L),rw),q),mul(mul(e[w],add(one,R)),K)).get(full,0)
 H,F,B=direct(birth(a,c),w)
 need((predh,predb)==(H,F-H)==(H,B),'birth defect recurrence')
 return b

def main():
 count=0;defects=set()
 for n in range(1,5):
  es=list(itertools.combinations(range(n),2))
  for code in range(1<<len(es)):
   a=[0]*n
   for k,(i,j) in enumerate(es):
    if code>>k&1:a[i]|=1<<j
    else:a[j]|=1<<i
   for w in range(n):
    for c in range(1<<n):defects.add(check(tuple(a),w,c));count+=1
 hits=[]
 for line in (ROOT/'search.out').read_text().splitlines():
  if line.startswith('SOURCE'):source=tuple(map(int,line.split()[1:]))
  if line.startswith('HIT'):
   p=line.split();d,x,y=map(int,p[1:4]);a=tuple(map(int,p[5:]));h,f,b=direct(a,0)
   need((h,f)==(3*source[0]+d,3*source[1]),'search hit permutation replay')
   hits.append({'source':source,'digit':d,'cuts':[x,y],'rows':a,'h':h,'F':f,'defect':b})
 markhits=[]
 for line in (ROOT/'search_allmarks.out').read_text().splitlines():
  if line.startswith('SOURCE'):source=tuple(map(int,line.split()[1:]))
  if line.startswith('HIT'):
   p=line.split();d,x,y=map(int,p[1:4]);w=int(p[5]);a=tuple(map(int,p[7:]));h,f,b=direct(a,w)
   need((h,f)==(3*source[0]+d,3*source[1]),'all-mark hit permutation replay')
   markhits.append({'source':source,'digit':d,'cuts':[x,y],'mark':w,'h':h,'F':f,'defect':b})
 out={'all_marks_survivors':markhits,'exhaustive_named_birth_cases':count,'orders':'every labelled tournament of order 1..4; every old mark; every birth cut','defect_range_observed':sorted(defects),'search_hits_permutation_verified':hits}
 (ROOT/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
