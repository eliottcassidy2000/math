def partition(edges,m):
 p=list(range(m+1))
 def root(a):
  while p[a]!=a:a=p[a]
  return a
 for a,b in edges:
  a,b=root(a),root(b)
  if a!=b:p[b]=a
 groups={}
 for a in range(m+1):groups.setdefault(root(a),[]).append(a)
 return tuple(sorted(tuple(v) for v in groups.values()))
def independent(edges,m):return len(partition(edges,m))==m+1-len(edges)
def enumerate_fiber(r,m,F):
 out=[]
 def visit(d,t):
  if d>r:out.append(t.copy());return
  for a in range(m-d+1):
   e=(a,a+d)
   if independent(list(t.values())+[e],m):t[d]=e;visit(d+1,t);del t[d]
 visit(1,F.copy());return out

