from itertools import combinations,permutations
def native(n):
 edges=list(combinations(range(n),2)); M=len(edges); ids={e:i for i,e in enumerate(edges)}
 cells=[set() for _ in range(3)]
 for w in permutations(range(n)):
  forced=ones=0
  for a,b in zip(w,w[1:]):
   bit=1<<ids[min(a,b),max(a,b)]; forced|=bit
   if a<b: ones|=bit
  opts=[1<<i for i in range(M) if forced>>i&1==0]
  for k in range(1,3):
   for frees in combinations(opts,k):
    F=sum(frees); O=((1<<M)-1)^forced^F
    d=O
    while True:
     cells[k].add((ones|d,F))
     if not d:break
     d=(d-1)&O
 parent=list(range(1<<M))
 def root(a):
  while parent[a]!=a: a=parent[a]
  return a
 ev={}
 for u,F in sorted(cells[1]):
  a,b=root(u),root(u|F)
  if a==b:ev[u,F]=1<<len(ev)
  else:parent[a]=b
 Z={}
 for u,F in sorted(cells[2]):
  a=F&-F;b=F^a
  w=ev.get((u,a),0)^ev.get((u|b,a),0)^ev.get((u,b),0)^ev.get((u|a,b),0)
  while w:
   lead=w.bit_length()
   if lead in Z:w^=Z[lead]
   else:Z[lead]=w;break
 print('n,V,E,faces,components,cycle rank,rank(square),H1F2:',n,1<<M,len(cells[1]),len(cells[2]),len({root(i) for i in range(1<<M)}),len(ev),len(Z),len(ev)-len(Z))

 H=[]
 for mask in range(1<<M):
  arc=lambda a,b:a!=b and (bool(mask>>ids[min(a,b),max(a,b)]&1)==(a<b))
  reach=[[a==b or arc(a,b) for b in range(n)] for a in range(n)]
  for k in range(n):
   for a in range(n):
    for b in range(n):reach[a][b]|=reach[a][k] and reach[k][b]
  single=[sum(reach[a][b] and reach[b][a] for b in range(n))==1 for a in range(n)]
  for i,(a,b) in enumerate(edges):
   if not arc(a,b):a,b=b,a
   forbidden=single[a] and single[b] and not any(arc(a,c) and arc(c,b) for c in range(n))
   F=1<<i
   if ((mask&~F,F) in cells[1]) == forbidden:raise RuntimeError("native arc classification")
  tr=[[0]*n for _ in range(1<<n)]
  for a in range(n):tr[1<<a][a]=1
  for S in range(1,1<<n):
   for a in range(n):
    if not tr[S][a]:continue
    for b in range(n):
     if S>>b&1:continue
     bit=mask>>ids[min(a,b),max(a,b)]&1
     if bit==(a<b):tr[S|(1<<b)][b]+=tr[S][a]
  H.append(sum(tr[-1]))
 start,finish=H.index(min(H)),H.index(max(H))
 parent=list(range(1<<M)); sizes=[1]*(1<<M)
 for w,u,v in sorted((abs(H[u]-H[u|F])//2,u,u|F) for u,F in cells[1]):
  a,b=root(u),root(v)
  if a!=b:
   if sizes[a]<sizes[b]:a,b=b,a
   sizes[a]+=sizes[b];parent[b]=a
  if root(start)==root(finish):break
 else:
  w=None
 print("  spectrum:",sorted(set(H)),"anchor bottleneck:",w)
if __name__=='__main__':
 for n in range(3,6):native(n)
