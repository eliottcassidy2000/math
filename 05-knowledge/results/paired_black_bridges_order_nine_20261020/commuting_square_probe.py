import json,sys,itertools
from pathlib import Path
from collections import Counter
sys.path.insert(0,str(Path('session14-metagraph-census').resolve()))
from census import schema,tournament,reflection

def encode(adj,p,tiles,n):return sum((not(adj[p[n-x]]>>p[n-y]&1))<<k for k,(x,y) in enumerate(tiles))
def path_reverse(adj,p):
 a=adj[:]
 for x,y in zip(p,p[1:]):a[x]^=1<<y;a[y]^=1<<x
 return a

def run(n):
 d=json.load(open(f'session14-metagraph-census/census_n{n}.json'));tiles,sigma=schema(n);full=(1<<len(tiles))-1
 def line(t):return min(t,full^t,reflection(t,sigma),full^reflection(t,sigma))
 def isblack(t):return t!=reflection(t,sigma)
 p=tuple(range(n));c=Counter();fails=[];certs=[]
 for t in range((full+1)//2):
  if not isblack(t) or line(t)!=t:continue
  a,b=d['pi'][t],d['pi'][full^t];key=tuple(sorted((a,b)))
  if a==b:c['loop']+=1;continue
  if d['nodes'][a]['reflection_pairs']==1 or d['nodes'][b]['reflection_pairs']==1:c['leaf']+=1;continue
  adj=tournament(t,n,tiles);adjp=path_reverse(adj,p);q=[];found=False;hasdisjoint=False
  def visit(used):
   nonlocal found,hasdisjoint
   if found:return
   if len(q)==n:
    hasdisjoint=True
    tq=encode(adj,q,tiles,n);tpq=encode(adjp,q,tiles,n)
    adjq=path_reverse(adj,q);tqp=encode(adjq,p,tiles,n)
    ls=[line(x) for x in [t,tq,tpq,tqp]]
    if all(isblack(x) for x in [tq,tpq,tqp]) and ls.count(t)==1:
     found=True;certs.append({'t':t,'q':q[:],'lines':ls});return
   poss=((adj[q[-1]] if q else (1<<n)-1)&~used)
   if q:poss&=~(1<<(q[-1]+1))
   while poss:
    z=poss&-poss;poss^=z;q.append(z.bit_length()-1);visit(used|z);q.pop()
  visit(0)
  if found:c['square']+=1
  else:c['no_square']+=1;c['has_disjoint_no_square' if hasdisjoint else 'no_disjoint']+=1;fails.append({'t':t,'ends':[a,b],'mass':[d['nodes'][x]['mass'] for x in [a,b]],'blue':[d['nodes'][x]['blue'] for x in [a,b]]})
 out={'n':n,'counts':dict(c),'failures':fails,'certificates':certs};Path(f'session15-bridge-classification/square_probe_n{n}.json').write_text(json.dumps(out)+'\n');print(n,c);print(fails[:10])
for n in [4,5,6,7]:run(n)
