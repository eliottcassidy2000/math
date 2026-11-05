"""Independent graph search for the unique-cycle/strong-bridge argument."""
import itertools,json,random,hashlib

def require(x,msg):
 if not x:raise RuntimeError(msg)

def reach(rows,active=None):
 if active is None:active=(1<<len(rows))-1
 rr=[((row|1<<i)&active) if active>>i&1 else 0 for i,row in enumerate(rows)]
 for k in range(len(rows)):
  if active>>k&1:
   for i in range(len(rows)):
    if rr[i]>>k&1:rr[i]|=rr[k]
 return rr

def strong(rows,active=None):
 if active is None:active=(1<<len(rows))-1
 if not active:return False
 rr=reach(rows,active)
 return all(rr[i]==active for i in range(len(rows)) if active>>i&1)

def cycles(rows,active=None):
 if active is None:active=(1<<len(rows))-1
 if active.bit_count()<3:return []
 root=(active&-active).bit_length()-1;ans=[]
 def visit(path,left):
  if len(ans)>=2:return
  if not left:
   if rows[path[-1]]>>root&1:ans.append(tuple(path))
   return
  todo=rows[path[-1]]&left
  while todo and len(ans)<2:
   b=todo&-todo;todo-=b;visit(path+[b.bit_length()-1],left^b)
 visit([root],active^(1<<root));return ans

def canonical(c):
 p=c.index(min(c));return tuple(c[p:]+c[:p])

def components(rows,active):
 rr=reach(rows,active);left=active;ans=[]
 while left:
  b=left&-left;i=b.bit_length()-1
  comp=sum(1<<j for j in range(len(rows)) if left>>j&1 and rr[i]>>j&1 and rr[j]>>i&1)
  ans.append(comp);left^=comp
 return sorted(ans,key=lambda c:-(rows[(c&-c).bit_length()-1]&active&~c).bit_count())

def cycle_ok(rows,c):
 return len(c)==len(rows) and set(c)==set(range(len(rows))) and all(rows[c[i]]>>c[(i+1)%len(c)]&1 for i in range(len(c)))

def graphs(n):
 pairs=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  rows=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1:rows[j]|=1<<i
   else:rows[i]|=1<<j
  yield rows

by_n=[];digest=hashlib.sha256();inheritance=0;cutchecks=0;strongtotal=0;twostrongtotal=0;uniquetotal=0
for n in range(1,7):
 counts=dict(order=n,tournaments=0,strong=0,two_strong=0,unique_cycle=0)
 for rows in graphs(n):
  counts['tournaments']+=1
  if not strong(rows):continue
  counts['strong']+=1;strongtotal+=1
  if n<3:continue
  full=(1<<n)-1;cc=cycles(rows);require(cc,'Camion finite control')
  deletions=[strong(rows,full^(1<<r)) for r in range(n)]
  if all(deletions):
   counts['two_strong']+=1;twostrongtotal+=1;require(len(cc)==2,('two-strong unique cycle',rows))
  if len(cc)!=1:continue
  counts['unique_cycle']+=1;uniquetotal+=1;C=cc[0];require(not all(deletions),'missing cut vertex')
  require(any(row.bit_count()==1 for row in rows),'no outdegree-one vertex')
  witnesses=[]
  for r,remstrong in enumerate(deletions):
   active=full^(1<<r)
   if remstrong:
    child=cycles(rows,active);expected=canonical([x for x in C if x!=r])
    require(len(child)==1 and child[0]==expected,('cycle inheritance',rows,r,C,child));inheritance+=1
    # Explicitly check every insertion gap into the unique deletion cycle.
    insertion=[];D=child[0]
    for i in range(len(D)):
     if rows[D[i]]>>r&1 and rows[r]>>D[(i+1)%len(D)]&1:
      E=list(D[:i+1])+[r]+list(D[i+1:]);insertion.append(canonical(E))
    require(insertion and all(E==C for E in insertion),('insertion',rows,r));continue
   cs=components(rows,active);require(len(cs)>=2,'bad cut components')
   aa=[a for a in range(n) if cs[0]>>a&1 and rows[r]>>a&1]
   bb=[b for b in range(n) if cs[-1]>>b&1 and rows[b]>>r&1]
   require(len(aa)==len(bb)==1,('endpoint choice failure',rows,r,cs,aa,bb))
   a=aa[0];b=bb[0];cut=cs[0]
   backwards=[(x,y) for x in range(n) if not(cut>>x&1) for y in range(n) if cut>>y&1 and rows[x]>>y&1]
   require(backwards==[(r,a)],('not a single reverse arc',rows,r,backwards))
   missing=list(rows);missing[r]&=~(1<<a);require(not strong(missing),('not a bridge',rows,r,a))
   other=list(rows);other[b]&=~(1<<r);require(not strong(other),('terminal bridge failed',rows,r,b))
   witnesses.append((r,a,b));cutchecks+=1
  digest.update(str((rows,C,witnesses)).encode())
 by_n.append(counts)
# Hostile free orientations: the surgeries use only prescribed step1/step2
# arcs and the single even-case chord. Test both choices for that chord.
rng=random.Random(197);surgery=0;odd=0;even_cases=[0,0]
for n in range(5,101):
 for trial in range(20):
  rows=[0]*n
  forced={}
  for i in range(n):
   for step in (1,2):
    j=(i+step)%n;key=tuple(sorted((i,j)));value=(i<j)
    require(key not in forced or forced[key]==value,'inconsistent forced arcs');forced[key]=value
  for i,j in itertools.combinations(range(n),2):
   fw=forced.get((i,j),bool(rng.randrange(2)))
   if n%2==0 and (i,j)==(1,4):fw=bool(trial%2)
   rows[i if fw else j]|=1<<(j if fw else i)
  if n%2:
   C=[2*i%n for i in range(n)];odd+=1
  elif rows[4]>>1&1:
   C=[0,2,4,1,3,5]+list(range(6,n));even_cases[0]+=1
  else:
   C=list(range(4,n,2))+[0,2]+list(range(3,n,2))+[1];even_cases[1]+=1
  require(cycle_ok(rows,C),'surgery failure')
  require(canonical(C)!=tuple(range(n)),'surgery did not change cycle');surgery+=1
print(json.dumps(dict(by_order=by_n,strong_tournaments=strongtotal,two_strong_tournaments=twostrongtotal,unique_cycle_tournaments=uniquetotal,strong_deletion_inheritance_checks=inheritance,cut_vertex_bridge_checks=cutchecks,unique_data_sha256=digest.hexdigest(),hostile_surgery_checks=surgery,odd_surgery_checks=odd,even_surgery_cases=even_cases),sort_keys=True,indent=2))
