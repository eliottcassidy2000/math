"""Two actual start-one sources; no formal inverse closure is assumed."""
from core import constants,glue,gadget,birth,require,inventories,adjacency_map

def source_bank(cap):
 if cap<1:return {}
 base=(1,1,1,1);out={base:()};queue=[base]
 for h,e,R,d in queue:
  moves=[('O',None,(h,0,2*R,d))]
  for L in range(1,cap//R+1):
   g,Q,B=constants(L)
   moves.append(('M',L,(g*h-(L%2)*e+L*R,h-int(L==1)*e,Q*h+B*R,h)))
  for kind,L,s in moves:
   if s[2]<=cap and s not in out:
    out[s]=out[h,e,R,d]+((kind,L),);queue.append(s)
 return out

def build_source(recipe):
 T=[0];r=0
 for kind,L in recipe:
  if kind=='M':T,r=glue(T,{r},gadget(L))
  else:T,_=birth(T,[r],0,1)
 return T,r

def converse(T):
 n=len(T);return [((1<<n)-1)^(1<<v)^row for v,row in enumerate(T)]

def left_response(kind,s):
 D,E,U,d=s
 if kind=='N':return D-E+U,D-E+d,D+3*U-d,U+E
 require(kind=='M','pair kind')
 return D+2*U,U+d,2*D-E+6*U-d,U+E

def value(kind,left,right):
 h,E,R,K=left_response(kind,left);g,s,L,d=right
 return h*g-E*s+R*L-K*d

MATRICES={
 'N':((6,4,5,5),(4,3,4,3),(5,4,5,3),(5,3,3,4)),
 'M':((11,8,10,10),(8,6,7,7),(10,7,9,8),(10,7,8,9)),
}

def positive_coordinates(s):
 D,E,U,d=s;return D-E-d,U-D,E,d

def matrix_value(kind,v,w):
 x=positive_coordinates(v);y=positive_coordinates(w);A=MATRICES[kind]
 return sum(x[i]*A[i][j]*y[j]for i in range(4)for j in range(4))

def construct(kind,V,c,W,a):
 require(kind in('N','M'),'pair kind')
 require(kind!='N'or len(V)>=2,'N excludes singleton left')
 require(len(W)>=2,'right source must be nontrivial')
 T,b=glue(V,{c},[0])
 if kind=='N':entries={b,c}
 else:
  T,u=glue(T,{b},[0]);entries={u,c}
 G=converse(W);return graft(T,entries,G,a)

def graft(T,entries,G,a):
 n=len(T);A=list(T)+[row<<n for row in G]
 for x in range(n):
  for y in range(len(G)):
   if y==a and x in entries:A[n+y]|=1<<x
   else:A[x]|=1<<(n+y)
 return A,n+a

def induced(T,vertices):
 return [sum(((T[x]>>y)&1)<<j for j,y in enumerate(vertices))for x in vertices]

def components(T,vertices):
 S=set(vertices);reach={v:({v}|{u for u in S if T[v]>>u&1})for v in S}
 for k in S:
  for v in S:
   if k in reach[v]:reach[v]|=reach[k]
 todo=set(S);out=[]
 while todo:
  v=min(todo);C=sorted(u for u in todo if u in reach[v]and v in reach[u]);out.append(C);todo.difference_update(C)
 return sorted(out,key=lambda C:-sum(bool(T[C[0]]>>D[0]&1)for D in out if D!=C))

def factor_graph(T):
 """Promised strong C2 input: exact two-source graph factorization, or R5."""
 n=len(T);ones=[v for v in range(n)if T[v].bit_count()==1]
 if not ones:
  require(n==5 and all(row.bit_count()==2 for row in T),'non-R5 C2 graph lacks degree-one vertex')
  return dict(kind='R5')
 v=ones[0];root=T[v].bit_length()-1
 U=components(T,[x for x in range(n)if x!=v])[0]
 require(root in U,'closing root is not in initial SCC')
 while True:
  C=components(T,[x for x in U if x!=root])[0]
  entries=[x for x in C if T[root]>>x&1]
  require(len(entries)in(1,2),'canonical start-two entry count')
  if len(entries)==2:break
  U=C;root=entries[0]
 left=list(C);first_entry=root;right=[x for x in range(n)if x not in left]
 for x in left:
  for y in right:require(bool(T[y]>>x&1)==(y==first_entry and x in entries),'bundled cross incidences')
 one=lambda x:[y for y in left if T[x]>>y&1]
 adjacent=next(((a,b)for a in entries for b in entries if a!=b and one(a)==[b]),None)
 if adjacent:
  a,c=adjacent;keep=[x for x in left if x!=a];order=keep+[a]+right;kind='N'
 else:
  u=next(a for a in entries if len(one(a))==1);b=one(u)[0];c=next(x for x in entries if x!=u)
  require(b not in entries and one(b)==[c],'three-vertex unique-start block')
  keep=[x for x in left if x not in(u,b)];order=keep+[b,u]+right;kind='M'
 V=induced(T,keep);W=converse(induced(T,right));c=keep.index(c);a=right.index(first_entry)
 A,newroot=construct(kind,V,c,W,a);require(adjacency_map(A,T,order),'exact reconstructed target graph')
 return dict(kind=kind,left_rows=V,left_root=c,right_rows=W,right_root=a,map_to_target=order,bundled_entry=first_entry)

def support(bound):
 bank=source_bank(max(0,(bound-1)//6));base=(1,1,1,1)
 sources=sorted((s for s in bank if s!=base),key=lambda s:(s[2],s))
 out={15:dict(kind='R5')}if bound>=15 else{}
 pairs=0
 for kind in('N','M'):
  for i,v in enumerate(sources):
   rv=bank[v]
   h,E,R,K=left_response(kind,v)
   if R>bound:continue
   for w in sources[i:]:
    lower=(3*v[2]*w[2]+1)if kind=='N'else(6*v[2]*w[2]+v[2]+w[2]+1)
    if lower>bound:break
    pairs+=1;N=value(kind,v,w)
    require(N>=R and N>=h*w[2]and N>=lower and N==matrix_value(kind,v,w)==value(kind,w,v),'positive bilinear bounds')
    if N<=bound:out.setdefault(N,dict(kind=kind,left=list(v),right=list(w),left_recipe=list(rv),right_recipe=list(bank[w])))
 return out,bank,pairs

def realize(record):
 if record['kind']=='R5':return [6,12,24,17,3],0
 V,c=build_source(record['left_recipe']);W,a=build_source(record['right_recipe'])
 return construct(record['kind'],V,c,W,a)

def inverse(target):
 out,_,_=support(target)
 return out.get(target)
