"""Independent permutation path counts for two literal carriers; standard library."""
import itertools,json
from pathlib import Path

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sig(rows,mask):
 n=len(rows);N=(1<<n)-1
 h=[0]*(N+1);e=[0]*(N+1);s=[0]*(N+1);h[0]=1
 for k in range(1,n+1):
  for path in itertools.permutations(range(n),k):
   if all(rows[a]>>b&1 for a,b in zip(path,path[1:])):
    S=sum(1<<v for v in path);h[S]+=1
    e[S]+=(mask>>path[-1])&1;s[S]+=(mask>>path[0])&1
 return (h[N],sum(h[A]*h[N^A] for A in range(N+1)),sum(e[A]*h[N^A] for A in range(N+1)),sum(h[A]*s[N^A] for A in range(N+1)),sum(e[A]*s[N^A] for A in range(N+1)),e[N],s[N])
def mul(x,y):
 if x is None:return y
 if y is None:return x
 h,p,F,G,J,e,s=x;H,P,f,g,j,E,S=y
 return(h*H,p*P,F*H+p*f,G*P+h*g,J*H+h*j+e*S+G*f,h*E,s*H)
unit=None
rows7=(2,20,9,3,45,79,31)
rows8=(2,4,33,23,71,155,47,95)
cases=[('24',rows7,42,(35,676,252,267,70,13,14)),('44',rows8,52,(89,2144,504,627,118,17,24))]
contexts=[unit]
for n in range(1,4):
 pairs=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  rows=[0]*n
  for j,(u,v) in enumerate(pairs):
   if bits>>j&1:rows[u]|=1<<v
   else:rows[v]|=1<<u
  for m in range(1<<n):contexts.append(sig(rows,m))
checks=0
for a,r,m,z in cases:
 need(sig(r,m)==z,'literal signature '+a)
 c=z[0]-z[4];need(c<=-3 and min(z[2:4])>-2*c-2,'lemma hypotheses')
 for x,y in itertools.product(contexts,repeat=2):
  zz=mul(mul(x,z),y)
  need(zz[0]-zz[4] not in (3*c+2,3*c+4),'finite context exclusion')
  checks+=1
modules=[]
for k in range(2,7):
 for S in itertools.combinations(range(7),k):
  if all(len({(rows7[v]>>x)&1 for x in S})==1 for v in set(range(7))-set(S)):modules.append(S)
need(modules==[(0,3)],'module list')
fa=min(sum(rows7[p[j]]>>p[i]&1 for i in range(7) for j in range(i+1,7)) for p in itertools.permutations(range(7)))
need(fa==3,'minimum feedback arcs')
result={'status':'PASS','literal_signatures':{a:list(z) for a,r,m,z in cases},'ordered_context_controls':checks,'contexts_each_side':len(contexts),'B7_proper_nontrivial_modules':modules,'B7_feedback_arc_number':fa,'scope':'General all-order context obstruction proved in PROOF.md; finite controls are not its proof.'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
