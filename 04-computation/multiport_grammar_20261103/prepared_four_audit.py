from itertools import permutations
from functools import cache
from pathlib import Path
import json, hashlib
G=0
def check(v,msg):
 global G
 G+=1
 if not v: raise RuntimeError(msg)
def profile(rows):
 n=len(rows); size=1<<n
 P=[0]*size; P[0]=1
 starts=[[0]*size for _ in rows]; ends=[[0]*size for _ in rows]
 joint=[[ [0]*size for _ in rows] for _ in rows]
 for length in range(1,n+1):
  for order in permutations(range(n),length):
   if all(rows[a]>>b&1 for a,b in zip(order,order[1:])):
    mask=sum(1<<v for v in order); P[mask]+=1
    starts[order[0]][mask]+=1; ends[order[-1]][mask]+=1
    joint[order[0]][order[-1]][mask]+=1
 return P,starts,ends,joint
def mul(A,B):
 size=len(A); out=[0]*size
 for i,a in enumerate(A):
  if not a: continue
  left=(size-1)^i; j=left
  while True:
   out[i|j]+=a*B[j]
   if not j: break
   j=(j-1)&left
 return out
cases=0
for n in range(2,6):
 pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
 for bits in range(1<<len(pairs)):
  rows=[0]*n
  for k,(i,j) in enumerate(pairs):
   if bits>>k&1: rows[i]|=1<<j
   else: rows[j]|=1<<i
  terminals=[u for u in range(n) if rows[u].bit_count()==1]
  if not terminals: continue
  P,A,W,Q=profile(rows); full=(1<<n)-1
  for u in terminals:
   r=rows[u].bit_length()-1; cases+=1
   S=mul(P,A[u])[full]; F=mul(P,W[r])[full]
   L=mul(mul(P,A[u]),W[r])[full]
   K=mul(mul(P,P),Q[u][r])[full]
   check(S-P[full]==W[r][full^(1<<u)],'general terminal start deficit')
   check(Q[u][r][(1<<u)|(1<<r)]==1 and sum(Q[u][r])==1,'only u-r joint path')
   check(F==L+K,'full literal ring joint equality')
# Independent recursive endpoint deletion for the actual graphs and first C2 outputs.
def observer(rows,r):
 n=len(rows); full=(1<<n)-1
 @cache
 def E(mask,v):
  if not mask>>v&1:return 0
  old=mask^(1<<v)
  if not old:return 1
  return sum(E(old,w) for w in range(n) if old>>w&1 and rows[w]>>v&1)
 @cache
 def P(mask):
  return sum(E(mask,v) for v in range(n) if mask>>v&1) if mask else 1
 H=P(full); F=sum(E(mask,r)*P(full^mask) for mask in range(1,full+1))
 return H,F
receipt=json.loads(Path('04-computation/multiport_grammar_20261103/prepared_four.json').read_text())
for rec in receipt['records']:
 r=rec['root']; rows=rec['rows']
 check(observer(rows,r)==(rec['H'],rec['F']),'actual input endpoint-deletion count')
 out=rec['prepared_C2_and_zero_channel'][0]
 check(observer(out['rows'],r)==(out['H'],out['F']),'actual C2 endpoint-deletion count')
 # Match every actual arc against D -> triangle with the same old root.
 k=rec['original_clock']-2; u=k+6; x=len(rows); y=x+1
 output=out['rows']
 check(all(output[k+d]>>v&1 for d in range(6) for v in (u,x,y)),'D dominates triangle')
 check(all(output[a]>>b&1 for a,b in ((u,y),(y,x),(x,u))),'actual triangle orientation')
 check(output[r]>>x&1 and output[r]>>y&1,'root beats both newborns')
# Literal nine-vertex suffix signature via independent endpoint recursion and convolutions.
D=(2,4,57,3,11,27)
rows=tuple(v|448 for v in D)+(128,256,64) # triangle u->x->y->u, isomorphic marked u
ports=(1<<4)|(1<<6); n=len(rows); full=(1<<n)-1
@cache
def end(mask,v):
 if not mask>>v&1:return 0
 rest=mask^(1<<v)
 if not rest:return 1
 return sum(end(rest,w) for w in range(n) if rest>>w&1 and rows[w]>>v&1)
@cache
def path(mask):
 return sum(end(mask,v) for v in range(n) if mask>>v&1) if mask else 1
@cache
def start(mask,v):
 if not mask>>v&1:return 0
 rest=mask^(1<<v)
 if not rest:return 1
 return sum(start(rest,w) for w in range(n) if rest>>w&1 and rows[v]>>w&1)
H=path(full); P2=sum(path(m)*path(full^m) for m in range(full+1))
f=sum(sum(end(m,v) for v in range(n) if ports>>v&1)*path(full^m) for m in range(full+1))
j=sum(sum(end(m,v) for v in range(n) if ports>>v&1)*sum(start(full^m,v) for v in range(n) if ports>>v&1) for m in range(full+1))
e=sum(end(full,v) for v in range(n) if ports>>v&1)
s=sum(start(full,v) for v in range(n) if ports>>v&1)
# g is the ordinary start-port two-slot response.
g=sum(path(m)*sum(start(full^m,v) for v in range(n) if ports>>v&1) for m in range(full+1))
check((H,P2,f,g,j,e,s)==(51,2688,756,579,134,17,12),'independent D-triangle full seven-signature')
print(json.dumps({'gates':G,'terminal_cases_all_labelled_n2_to5':cases,'first_outputs_recounted':3,'signature':[H,P2,f,g,j,e,s]},sort_keys=True))
