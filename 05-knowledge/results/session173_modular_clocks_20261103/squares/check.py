"""Exact checker, no repository imports and no disabled assertions."""
from itertools import product, combinations
from math import isqrt
from pathlib import Path
import json,random
HERE=Path(__file__).resolve().parent

def need(b,msg):
 if not b:raise RuntimeError(msg)
def walk(word,x):
 p=[x]
 for s in word:p.append(s*s-p[-1])
 return p
def offsets(word):return walk(word,0)
def summary(word):
 c=offsets(word)
 return {'sign':(-1)**len(word),'offset':c[-1],'even':c[::2],'odd':c[1::2]}
def allowed(word,L,U):
 c=offsets(word);e=c[::2];o=c[1::2]
 if len(set(e))<len(e) or len(set(o))<len(o):return []
 lo=max([L-v for v in e]+[v-U for v in o]);hi=min([U-v for v in e]+[v-L for v in o])
 bad={(b-a)//2 for a in e for b in o if (b-a)%2==0}
 return [x for x in range(lo,hi+1) if x not in bad]
def valid(p,L,U):return len(p)==len(set(p)) and all(L<=x<=U for x in p)
def square(x):return x>=0 and isqrt(x)**2==x
def checkpath(p,n=None):
 need(len(set(p))==len(p),'collision')
 need(all(square(a+b) for a,b in zip(p,p[1:])),'nonsquare')
 if n:need(sorted(p)==list(range(1,n+1)),'not Hamilton')
 return [isqrt(a+b) for a,b in zip(p,p[1:])]

word_checks=0
for k in range(6):
 for w in product(range(2,8),repeat=k):
  a=allowed(w,1,44);b=[x for x in range(1,45) if valid(walk(w,x),1,44)]
  need(a==b,('placement',w));word_checks+=1

# Independent exhaustive set-partition test of interval-shift rigidity.
shift_checks=0
for n in range(2,15,2):
 whole=set(range(1,n+1))
 for aa in combinations(range(1,n+1),n//2):
  A=set(aa);B=whole-A
  for t in range(1,n):
   actual=({x+t for x in A}|{x-t for x in B}==whole and {x+t for x in A}.isdisjoint({x-t for x in B}))
   predicted=n%(2*t)==0 and A=={x for x in whole if (x-1)//t%2==0}
   need(actual==predicted,('shift rigidity',n,t,A));shift_checks+=1

# Exact twin witness.
twin=json.loads((HERE/'twins.json').read_text());p=twin['path'];q=twin['second'];w=checkpath(p,36)
need(checkpath(q,36)==w,'root twin mismatch');need(allowed(w,1,36)==[3,9],'twin starts')

# Reproducible rotations. The source-supplied seed is independently verified.
r=random.Random(173)
p=[1,15,10,26,23,2,14,22,27,9,16,20,29,7,18,31,5,11,25,24,12,13,3,6,30,19,17,32,4,21,28,8]
checkpath(p,32);need(square(p[0]+p[-1]),'seed not cycle')
paths={32:p[:]};counts={};all_word_starts={}
for n in range(33,528):
 for tries in range(200000):
  if square(p[-1]+n):p.append(n);break
  if square(p[0]+n):p.reverse();p.append(n);break
  options=[i for i in range(len(p)-2) if square(p[-1]+p[i])]
  if options:
   i=r.choice(options);p=p[:i+1]+p[:i:-1]
  if r.randrange(2):p.reverse()
 else:raise RuntimeError(('rotation search stalled',n))
 w=checkpath(p,n);counts[n]=tries
 if n in (44,72,527):
  paths[n]=p[:];all_word_starts[n]=allowed(w,1,n)
  need(len(all_word_starts[n])<=(1 if n%2 else 2),'moment bound')

# Nonzero arithmetic translation whose two modular clocks close.
segment=[230,395,446];need(any(p[i:i+3]==segment for i in range(525)),'segment absent')
need(checkpath(segment)==[25,29],'segment roots');need(segment[-1]-segment[0]==216,'voltage')
for m in (8,9,72):need(segment[0]%m==segment[-1]%m,'clock closure')

# Minimum-connector certificates for the user's named anchors.
base={4,2,5,3,7,13,44,527}
def components(V):
 todo=set(V);out=[]
 while todo:
  v=min(todo);todo.remove(v);C={v};stack=[v]
  while stack:
   x=stack.pop()
   for y in list(todo):
    if square(x+y):todo.remove(y);C.add(y);stack.append(y)
  out.append(sorted(C))
 return out
p0=[527,2,7,317,44,5,4,12,13,3];p1=p0+[1,8]
anchor=[]
for V,P in ((base,p0),(base|{8},p1)):
 roots=checkpath(P);C=components(V);extra=sorted(set(P)-V)
 need(V<=set(P),'missing anchor');need(len(extra)==len(C)-1,'not cardinality optimum')
 anchor.append({'anchors':sorted(V),'components':C,'path':P,'roots':roots,'extra':extra,'minimum_extra':len(C)-1})

# Exhaustive difference-of-squares seams, with explicit occupation exclusions.
def connector_candidates(a,b,occupied=()):
 d=abs(a-b);out=[]
 for f in range(1,isqrt(d)+1):
  if d%f:continue
  g=d//f
  if (f+g)%2:continue
  z=((g-f)//2)**2-min(a,b)
  if z>0 and z not in set(occupied)|{a,b}:out.append(z)
 return sorted(out)
need(connector_candidates(2,7)==[], 'endpoint-repeat hostile')
height_rows=[]
for a in (7,527):
 for b in (3,4,13,44,8):
  C=connector_candidates(a,b,base|{8})
  need(all(square(a+z) and square(b+z) for z in C),'factor candidate')
  height_rows.append({'a':a,'b':b,'connectors':C})
need(min(z for row in height_rows for z in row['connectors'])==317,'height optimum')

print(json.dumps({'minimum_count_height_rows':height_rows,'minimum_count_height':317,'word_interval_checks':word_checks,'shift_partition_checks':shift_checks,'twins':twin,'twin_roots':checkpath(twin['path']), 'twin_start_interval':[3,9],'twin_distinct_counts':[len(set(walk(checkpath(twin['path']),x))) for x in range(3,10)],'Hamilton_paths':paths,'admissible_word_starts':all_word_starts,'rotation_total':sum(counts.values()),'rotation_max':max(counts.values()),'rotation_counts':counts,'closed_clock_nonclosed_integer_path':segment,'clock_roots':[25,29],'integer_translation':216,'mod72_winding':3,'optimal_anchor_augmentations':anchor},sort_keys=True,indent=2))
