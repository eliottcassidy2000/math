from pathlib import Path
import sys,json,collections,hashlib
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'source'))
from twoport_core import transitions,brute

def need(b,m):
 if not b:raise RuntimeError(m)
saved=json.loads((P/'source/CONTROL.json').read_text());V=[tuple(v)for v in saved['vectors']];idx={v:i for i,v in enumerate(V)};need(len(idx)==410,'distinct raw vectors');rej=410
# Independent explicit fixed-port cycle readout and sparse matrix multiplication.
read=[(),(1,2),(3,5),(2,5)]
raw=[]
for v in V:
 row=[]
 for c in range(4):
  C=sum(v[i]for i in read[c])
  if C>=3:row.append((rej,0));continue
  w=[0]*27
  for i,x in enumerate(v):
   for j,m in transitions(2,c,2)[i]:w[j]+=m*x
  w=tuple(min(3,x)for x in w);need(w in idx,'raw transition leaves certificate');row.append((idx[w],int(C==2)))
 raw.append(row)
raw.append([(rej,0)]*4)
need(raw==[[tuple(x)for x in row]for row in saved['raw_edges']],'all raw transitions')
need(V[0]==brute([2,0],[0,1]),'physical seed')
seen={0};todo=[0]
for q in todo:
 for qq,a in raw[q]:
  if qq not in seen:seen.add(qq);todo.append(qq)
need(seen==set(range(411)),'complete raw reachable universe')
# Direct glob NFAs: stars are loops at the current state, no Thompson compiler,
# no epsilon closure, and no import of the source regular-expression functions.
def one(c):return frozenset([c])
def star(*cs):return ('star',frozenset(cs))
A=[one(1),star(1,3)];B=[one(2),star(2,3)];R=[star(0,1),one(1)];S=[star(0,2),one(2)]
patterns=[A+[one(2)]+R,A+[one(0),one(3)],A+[one(0),one(3)]+R]
for pre in ([],B):
 for suf in ([],S):patterns.append(pre+[one(0),one(2),one(3)]+suf)
for suf in ([],S):patterns.append(B+[one(1),one(3)]+suf)
patterns += [A+[one(0),star(0),one(3)],B+[one(1),one(0),star(0),one(3)],[one(0),one(1),one(3)]]
phi=[3,1,2,0]
def reverse(p):
 out=[]
 for t in p[::-1]:
  if isinstance(t,tuple):out.append(('star',frozenset(phi[x]for x in t[1])))
  else:out.append(frozenset(phi[x]for x in t))
 return out
patterns=patterns+[reverse(p)for p in patterns];nfa=[];initial=[];final=set()
def state():nfa.append([set()for _ in range(4)]);return len(nfa)-1
for p in patterns:
 q=state();initial.append(q)
 for t in p:
  if isinstance(t,tuple):
   for c in t[1]:nfa[q][c].add(q)
  else:
   qq=state()
   for c in t:nfa[q][c].add(qq)
   q=qq
 final.add(q)
start=(0,0,frozenset(initial));todo=[start];seen={start};parent={start:None};labels={};edgecount=0
for node in todo:
 q,a,X=node
 need(bool(a)==bool(final.intersection(X)),('all-length mismatch',node))
 for c in range(4):
  qq,aa=raw[q][c];Y=frozenset(z for x in X for z in nfa[x][c]);nxt=(qq,aa,Y);edgecount+=1
  if nxt not in seen:seen.add(nxt);todo.append(nxt);parent[nxt]=node;labels[nxt]=c
# Check source minimization is an exact edge-output quotient, independently.
cls=saved['class_map'];quot=saved['quotient']
for q in range(411):
 for c,(qq,a)in enumerate(raw[q]):need(quot[cls[q]][c]==[cls[qq],a],'quotient congruence')
print(json.dumps({'status':'PASS','raw_states_with_reject':411,'raw_edges':1644,'glob_patterns':len(patterns),'independent_nfa_states':len(nfa),'complete_raw_product_states':len(todo),'complete_raw_product_edges':edgecount,'source_controller_classes':len(quot),'source_control_sha256':hashlib.sha256((P/'source/CONTROL.json').read_bytes()).hexdigest()},indent=2,sort_keys=True))
