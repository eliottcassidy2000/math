"""Promised-input physical two-port compiler for complete C1 and C2 classes.

Inputs retain their original vertex labels. At most two current ports are exposed;
no forgotten old vertex can ever be re-exposed. No scalar-coverage claim is made.
"""
from itertools import permutations


def need(x,s):
 if not x:raise RuntimeError(s)

def induced(rows,vs):
 return tuple(sum(1<<j for j,w in enumerate(vs) if rows[v]>>w&1) for v in vs)

def scc(rows,vs):
 vs=tuple(sorted(vs));reach={v:{v}|{u for u in vs if rows[v]>>u&1} for v in vs}
 for k in vs:
  for v in vs:
   if k in reach[v]:reach[v]|=reach[k]
 left=set(vs);parts=[]
 while left:
  v=min(left);C={u for u in left if u in reach[v] and v in reach[u]};parts.append(tuple(sorted(C)));left-=C
 snapshot=tuple(parts)
 parts.sort(key=lambda C:sum(bool(rows[C[0]]>>D[0]&1) for D in snapshot if C!=D),reverse=True)
 return parts


class Compiler:
 def __init__(self,rows):
  self.rows=tuple(rows);self.ports=[];self.born=[];self.steps=[]
 def add(self,v,keep):
  need(v not in self.born,'vertex born twice')
  beaten={u for u in self.born if self.rows[v]>>u&1}
  need(beaten<=set(self.ports),('unexposed old outneighbor',v,beaten,self.ports))
  mask=sum(1<<i for i,u in enumerate(self.ports) if u in beaten)
  nextports=self.ports+[v];need(set(keep)<=set(nextports),'requested old vertex was forgotten')
  drop=None
  if len(nextports)>2:
   drop=next(i for i,u in enumerate(nextports) if u not in keep);del nextports[drop]
  self.born.append(v);self.ports=nextports
  self.steps.append(dict(vertex=v,mask=mask,drop=drop,ports=nextports[:]))
 def one(self,V,requested):
  """C1 graph (or singleton), exposing any given one/two start-one vertices."""
  V=set(V);requested=set(requested)
  need(1<=len(requested)<=2 and requested<=V,'bad requested start-one set')
  if len(V)==1:
   need(requested==V,'singleton marks');self.add(next(iter(V)),V);return
  b=min(requested);parts=scc(self.rows,V-{b});U=set(parts[0]);tail=[v for C in parts[1:] for v in C]
  need(all(len(C)==1 for C in parts[1:]),'start-one tail is not transitive')
  entry={v for v in U if self.rows[b]>>v&1};need(len(entry)==1,'start-one entry size')
  a=next(iter(entry));other=next(iter(requested-{b}),None)
  old_requested={a}
  if not tail and other is not None and other!=a:old_requested.add(other)
  if tail and other in U:need(other==a,'nonempty-tail old marked vertex')
  self.one(U,old_requested)
  keep={b}|({other} if other in self.born else set());self.add(b,keep)
  for w in tail:
   keep={b}|({other} if other in self.born or other==w else set());self.add(w,keep)
  need(requested<=set(self.ports),'final start-one pair lost')
 def two(self,V,root):
  """Strong graph with two Hamilton paths starting at root; expose root."""
  V=set(V);parts=scc(self.rows,V-{root});U=set(parts[0]);tail=[v for C in parts[1:] for v in C]
  need(len(U)>=3 and all(len(C)==1 for C in parts[1:]),'start-two component form')
  entry={v for v in U if self.rows[root]>>v&1};need(len(entry) in (1,2),'start-two entry size')
  if len(entry)==1:self.two(U,next(iter(entry)))
  else:self.one(U,entry)
  self.add(root,{root})
  for w in tail:self.add(w,{root})
  need(root in self.ports,'start-two root lost')
 def regular_five(self):
  R=(6,12,24,17,3)
  order=next(p for p in permutations(range(5)) if induced(self.rows,p)==R)
  self.add(order[0],{order[0]});self.add(order[1],set(order[:2]))
  for v in order[2:]:self.add(v,set(order[:2]))
 def normalize_seed(self):
  # The standard two-port language starts from a forward ordered edge.
  # Swap the first two births if necessary, then re-index every later mask/drop.
  # Their exposed SET is unchanged before every later birth.
  if len(self.born)<2:return
  old=self.steps[:];a,b=self.born[:2]
  if self.rows[b]>>a&1:a,b=b,a
  born=[a,b];ports=[a,b]
  steps=[dict(vertex=a,mask=0,drop=None,ports=[a]),dict(vertex=b,mask=0,drop=None,ports=[a,b])]
  for step in old[2:]:
   v=step['vertex'];keep=set(step['ports']);mask=sum(1<<i for i,u in enumerate(ports) if self.rows[v]>>u&1)
   ports.append(v);drop=next(i for i,u in enumerate(ports) if u not in keep);del ports[drop]
   born.append(v);steps.append(dict(vertex=v,mask=mask,drop=drop,ports=ports[:]))
  self.born,self.ports,self.steps=born,ports,steps
 def verify(self):
  self.normalize_seed()
  need(sorted(self.born)==list(range(len(self.rows))),'graph vertices missing')
  ports=[];born=[];out=[0]*len(self.rows)
  for s in self.steps:
   v,mask,drop=s['vertex'],s['mask'],s['drop']
   need(0<=mask<1<<len(ports),'invalid mask')
   chosen={u for i,u in enumerate(ports) if mask>>i&1}
   for u in born:out[v if u in chosen else u]|=1<<(u if u in chosen else v)
   born.append(v);ports.append(v)
   if len(ports)==3:
    need(drop in (0,1,2),'missing drop');del ports[drop]
   else:need(drop is None,'unexpected drop')
   need(ports==s['ports'] and len(ports)<=2,'port receipt')
  need(tuple(out)==self.rows,'physical adjacency mismatch')
  return dict(rows=list(self.rows),birth_order=self.born[:],ports=self.ports[:],steps=self.steps[:])


def compile_one(rows,requested):
 c=Compiler(rows);c.one(set(range(len(rows))),set(requested));return c.verify()
def compile_two(rows,root):
 c=Compiler(rows);c.two(set(range(len(rows))),root);return c.verify()
def compile_regular_five(rows=(6,12,24,17,3)):
 c=Compiler(rows);c.regular_five();return c.verify()
