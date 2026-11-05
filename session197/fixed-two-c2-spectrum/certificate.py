import sys,json
from collections import deque
from control_core import brute,cycle_readout,capped_step

def need(x,s):
 if not x:raise RuntimeError(s)
V=[brute([2,0],[0,1])];seen={V[0]:0};E=[]
for v in V:
 row=[]
 for S in range(4):
  C=cycle_readout(v,S)
  if C>2:row.append((-1,0));continue
  w=capped_step(v,S,2,3)
  if w not in seen:seen[w]=len(V);V.append(w)
  row.append((seen[w],int(C==2)))
 E.append(row)
reject=len(V);E=[[(q if q>=0 else reject,a)for q,a in row]for row in E]+[[(reject,0)]*4]
cls=[0]*len(E)
while True:
 D={};new=[]
 for row in E:
  key=tuple((a,cls[q])for q,a in row)
  if key not in D:D[key]=len(D)
  new.append(D[key])
 if new==cls:break
 cls=new
Q=[None]*len(set(cls))
for p,row in enumerate(E):
 x=tuple((cls[q],a)for q,a in row)
 if Q[cls[p]] is not None:need(Q[cls[p]]==x,'unstable quotient')
 Q[cls[p]]=x
# Simple regular-expression syntax: lit, eps, concat, union, star.
def lit(w):return ('lit',tuple(w))
def eps():return('eps',)
def cat(*r):return('cat',r)
def alt(*r):return('alt',r)
def star(r):return('star',r)
A=cat(lit([1]),star(alt(lit([1]),lit([3]))))
B=cat(lit([2]),star(alt(lit([2]),lit([3]))))
R=cat(star(alt(lit([0]),lit([1]))),lit([1]))
S=cat(star(alt(lit([0]),lit([2]))),lit([2]))
N=cat(A,lit([2]),R)
M=cat(A,lit([0,3]),alt(eps(),R))
P=cat(alt(eps(),B),lit([0,2,3]),alt(eps(),S))
QQ=cat(B,lit([1,3]),alt(eps(),S))
F=cat(A,lit([0]),star(lit([0])),lit([3]))
G=cat(B,lit([1,0]),star(lit([0])),lit([3]))
T=alt(N,M,P,QQ,F,G,lit([0,1,3]))
phi=(3,1,2,0)
def dual(r):
 k=r[0]
 if k=='lit':return lit([phi[c]for c in r[1][::-1]])
 if k=='eps':return r
 if k=='cat':return cat(*(dual(x)for x in r[1][::-1]))
 if k=='alt':return alt(*(dual(x)for x in r[1]))
 return star(dual(r[1]))
T=alt(T,dual(T));nfa=[]
def state():nfa.append([set()for _ in range(5)]);return len(nfa)-1
def build(r):
 a,b=state(),state();k=r[0]
 if k=='eps':nfa[a][4].add(b)
 elif k=='lit':
  p=a
  for c in r[1]:q=state();nfa[p][c].add(q);p=q
  nfa[p][4].add(b)
 elif k=='cat':
  p=a
  for x in r[1]:s,t=build(x);nfa[p][4].add(s);p=t
  nfa[p][4].add(b)
 elif k=='alt':
  for x in r[1]:s,t=build(x);nfa[a][4].add(s);nfa[t][4].add(b)
 else:
  s,t=build(r[1]);nfa[a][4].update((b,s));nfa[t][4].update((b,s))
 return a,b
start,end=build(T)
def closure(X):
 X=set(X);todo=list(X)
 for p in todo:
  for q in nfa[p][4]:
   if q not in X:X.add(q);todo.append(q)
 return frozenset(X)
initial=(cls[0],0,closure([start]));queue=[initial];visited={initial};parent={initial:None};letters={}
for node in queue:
 q,a,X=node
 if bool(a)!=(end in X):
  w=[];p=node
  while parent[p]is not None:w.append(letters[p]);p=parent[p]
  raise RuntimeError(('language mismatch',w[::-1],q,a,end in X))
 for c in range(4):
  qq,aa=Q[q][c];Y=closure(set().union(*(nfa[p][c]for p in X)))
  nxt=(qq,aa,Y)
  if nxt not in visited:visited.add(nxt);queue.append(nxt);parent[nxt]=node;letters[nxt]=c
need((len(V),len(Q),len(nfa),len(queue))==(410,15,450,47),'frozen dimensions')
product_index={x:i for i,x in enumerate(queue)}
product_edges=[]
for q,a,X in queue:
 row=[]
 for c in range(4):
  qq,aa=Q[q][c];Y=closure(set().union(*(nfa[p][c]for p in X)))
  row.append(product_index[(qq,aa,Y)])
 product_edges.append(row)
result=dict(initial_class=cls[0],reject_class=cls[reject],regex_nfa=[[sorted(z)for z in row]for row in nfa],regex_start=start,regex_end=end,product_states=[[q,a,sorted(X)]for q,a,X in queue],product_edges=product_edges,status='PASS',raw_vectors=len(V),control_classes=len(Q),regex_nfa_states=len(nfa),complete_product_states=len(queue),vectors=[list(v)for v in V],raw_edges=E,class_map=cls,quotient=Q)
if __name__=='__main__':print(json.dumps(result,sort_keys=True,separators=(',',':')))
