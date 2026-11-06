"""Independent exact overlap-injection and forward-sector controls."""
import functools,itertools,json

def need(b,m):
 if not b:raise RuntimeError(m)
def rows_from_mask(n,mask):
 out=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if (mask>>k)&1:out[i]|=1<<j
   else:out[j]|=1<<i
   k+=1
 return tuple(out)
def paths_function(rows):
 @functools.cache
 def paths(S):
  if not S:return ((),)
  out=[];todo=S
  while todo:
   bit=todo&-todo;todo-=bit;first=bit.bit_length()-1
   for tail in paths(S^bit):
    if not tail or rows[first]>>tail[0]&1:out.append((first,)+tail)
  return tuple(out)
 return paths

def merged(rows,P,Q,shared):
 common=[v for v in P if shared>>v&1]
 need(common==[v for v in Q if shared>>v&1],'incompatible')
 answer=[];p0=q0=0
 for anchor in common+[None]:
  p1=P.index(anchor) if anchor is not None else len(P)
  q1=Q.index(anchor) if anchor is not None else len(Q)
  a=list(P[p0:p1]);b=list(Q[q0:q1])
  while a and b:
   if rows[a[0]]>>b[0]&1:answer.append(a.pop(0))
   else:answer.append(b.pop(0))
  answer+=a+b
  if anchor is not None:answer.append(anchor)
  p0=p1+1;q0=q1+1
 return tuple(answer)
def scc_sizes(rows):
 n=len(rows);reach=[row|(1<<v) for v,row in enumerate(rows)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 unseen=set(range(n));parts=[]
 while unseen:
  i=min(unseen);C={j for j in unseen if reach[i]>>j&1 and reach[j]>>i&1};parts.append(C);unseen-=C
 return parts

tournaments=cases=pairchecks=authorcases=authorpairs=forward1=forward2=0
for n in range(6):
 for code in range(1<<(n*(n-1)//2)):
  rows=rows_from_mask(n,code);getpaths=paths_function(rows);allmask=(1<<n)-1;full=set(getpaths(allmask));tournaments+=1
  parts=scc_sizes(rows);nontrivial=[C for C in parts if len(C)>1]
  for u in range(n):
   for v in range(n):
    if not(rows[u]>>v&1):continue
    f=sum(p.index(u)<p.index(v) for p in full)
    transitive=not nontrivial
    triangle_padding=len(nontrivial)==1 and len(nontrivial[0])==3 and {u,v}<=nontrivial[0]
    need((f==1)==transitive,'iff forward1 structure')
    need((f==2)==triangle_padding,'iff forward2 structure')
    if f==1:need(len(full)==1,'forward1 reverse0');forward1+=1
    if f==2:need(len(full)==3,'forward2 reverse1');forward2+=1
  for labels in itertools.product(range(3),repeat=n):
   X=sum(1<<v for v in range(n) if labels[v]!=1);Y=sum(1<<v for v in range(n) if labels[v]!=0);S=X&Y
   image=set();expected=0;private=0 in labels and 1 in labels
   for P in getpaths(X):
    order=tuple(v for v in P if S>>v&1)
    for Q in getpaths(Y):
     if order!=tuple(v for v in Q if S>>v&1):continue
     R=merged(rows,P,Q,S)
     need(R in full,'merged path legality')
     need(tuple(v for v in R if X>>v&1)==P,'X recovery')
     need(tuple(v for v in R if Y>>v&1)==Q,'Y recovery')
     need(R not in image,'input collision')
     image.add(R);expected+=1;pairchecks+=1
     if private:authorpairs+=1
   need(len(image)==expected<=len(full),'overlap count')
   if S.bit_count()<=1:need(expected==len(getpaths(X))*len(getpaths(Y)),'small overlap product')
   cases+=1
   if private:authorcases+=1
need((tournaments-1,authorcases,authorpairs)==(1099,187620,515524),'author universe agreement')
sharp=(10,28,25,16,1);sp=paths_function(sharp)
need(len(sp(31))==9 and len(sp(7))==len(sp(25))==3,'strong sharp product')
need(len(scc_sizes(sharp))==1,'sharp graph strong')
hostile=rows_from_mask(6,91);pp=paths_function(hostile)(63)
f=sum(p.index(0)<p.index(1) for p in pp);r=len(pp)-f
need(hostile==(54,8,2,5,14,30) and hostile[0]&2,'hostile graph')
need((f,r)==(7,8),'majority hostile')
# Independently count every path in the entire rigid strong-block family through m=10.
rigid=[]
for m in range(1,11):
 n=m+2;rows=[0]*n;rows[0]=2;rows[1]=sum(1<<i for i in range(2,n))
 for i in range(2,n):rows[i]=1|sum(1<<j for j in range(i+1,n))
 hp=paths_function(tuple(rows))((1<<n)-1);f=sum(p.index(0)<p.index(1) for p in hp)
 need((f,len(hp)-f)==(2**m,1),'rigid blowup sectors')
 rigid.append({'m':m,'forward':f,'reverse':len(hp)-f})
print(json.dumps({'status':'PASS','scope':'finite controls support separately audited all-order proof','tournaments_including_empty':tournaments,'all_ordered_covering_pairs':cases,'compatible_path_pairs':pairchecks,'author_private_nonempty_cases':authorcases,'author_compatible_path_pairs':authorpairs,'arc_forward1_structure_checks':forward1,'arc_forward2_structure_checks':forward2,'sharp_H':9,'sharp_strong':True,'hostile_forward_reverse':[f for f in (7,8)],'rigid_family':rigid},indent=2,sort_keys=True))
