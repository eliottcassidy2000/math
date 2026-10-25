#!/usr/bin/env python3
import ast,json
from pathlib import Path
from itertools import combinations,permutations

def require(x,msg):
 if not x:raise RuntimeError(msg)
p=Path(__file__).resolve().parent.parent/'cycle-moves/check.py'
tree=ast.parse(p.read_text());tree.body=[node for node in tree.body if isinstance(node,(ast.Import,ast.ImportFrom,ast.FunctionDef))]
ns={};exec(compile(tree,str(p),'exec'),ns)

def replay(n,E,m,moves):
 idx={e:k for k,e in enumerate(E)}
 for c in moves:
  require(len(set(c))==len(c),'repeated vertex')
  require(sum(u in c and v in c for u,v in E)==len(c),'not induced')
  flip=0
  for u,v in zip(c,c[1:]+c[:1]):
   require((min(u,v),max(u,v)) in idx,'missing edge')
   bit=1<<idx[min(u,v),max(u,v)]
   require(bool(m&bit)==(u<v),'illegal cycle')
   flip|=bit
  m^=flip
 return m
routes=0
for n in range(1,8):
 E=ns['edges'](n);fibres={}
 for m in range(1<<len(E)):
  s=ns['score'](n,E,m)
  if s not in fibres:fibres[s]=m
  else:
   moves=ns['route'](n,E,m,fibres[s]);require(replay(n,E,m,moves)==fibres[s],'bad route');routes+=1
 print('routes',n,'PASS',flush=True)
checks=0
for n in range(1,6):
 P=list(combinations(range(n),2))
 for graph in range(1<<len(P)):
  E=[e for j,e in enumerate(P) if graph>>j&1];idx={e:j for j,e in enumerate(E)}
  for C in ns['cycles'](n,E):
   fixed={idx[min(u,v),max(u,v)]:u<v for u,v in zip(C,C[1:]+C[:1])};free=[j for j in range(len(E)) if j not in fixed]
   cm=sum(1<<j for j in fixed)
   for x in range(1<<len(free)):
    m=sum(1<<j for j,yes in fixed.items() if yes)|sum(1<<j for i,j in enumerate(free) if x>>i&1)
    moves=[];end=ns['reverse_cycle'](E,idx,m,tuple(C),moves)
    require(end==m^cm,'wrong net cycle');require(len(moves)<=len(C)-2,'move bound');require(replay(n,E,m,moves)==end,'bad chord simulation');checks+=1
 print('generic cycle replay',n,'PASS',flush=True)
result={'status':'PASS','producer_routes_replayed':routes,'generic_cycle_reversals_replayed':checks,'explicit_exception_checks':True}
Path(__file__).with_name('replay-result'+('-O' if not __debug__ else '')+'.json').write_text(json.dumps(result,indent=2)+'\n')
print(result)
