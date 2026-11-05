"""Physical complete C3 source union. No claim of presentation uniqueness."""
from core import require,glue,inventories,cycles
from source_oracle import source_bank,build_source,converse,graft,induced,components

ONE=(1,1,1,1)
TYPE_I=('A','B','C','T')
MATRICES=(
 ((0,0,0,3),(0,0,0,-3),(0,0,12,-3),(3,0,-12,3)),
 ((0,0,0,3),(0,0,0,-2),(0,0,12,-6),(3,-1,-9,5)),
 ((0,0,0,3),(0,0,0,-1),(0,0,12,-9),(3,-2,-6,5)),
 ((0,0,0,3),(0,0,0,0),(0,0,12,-12),(3,-3,-3,3)))
POSITIVE=(
 ((12,12,12,12),(12,12,12,9),(12,12,12,9),(3,0,3,6)),
 ((12,12,12,9),(12,12,12,6),(12,12,12,7),(6,3,5,8)),
 ((12,12,12,6),(12,12,12,3),(12,12,12,5),(9,6,7,8)),
 ((12,12,12,3),(12,12,12,0),(12,12,12,3),(12,9,9,6)))

def response(kind,s):
 D,E,U,d=s
 if kind=='A':return D-E+3*U-2*d,D-E+2*U-2*d,3*D-2*E+7*U-2*d,2*D-2*E+2*U
 if kind=='B':return D-E+3*U-2*d,U,2*D-E+8*U-4*d,U+E
 if kind=='C':return D+4*U,U+d,4*D-E+14*U-d,U+E+2*d
 require(kind=='T','unknown Type I')
 return D+2*U,D+U+d,3*D-E+8*U-d,2*D+4*U

def left_graph(kind,V,c):
 if kind=='T':
  T,b=glue(V,{c},[0]);T,u=glue(T,{b},[0]);return T,{u,b,c}
 if kind=='C':
  x=c;T,c=glue(V,{c},[0]);T,a=glue(T,{c},[0]);T,r=glue(T,{a,c},[0]);return T,{r,x}
 require(kind in('A','B'),'unknown Type I')
 T,a=glue(V,{c},[0]);T,r=glue(T,{a,c},[0]);return T,{r,a if kind=='A' else c}

def actual_response(T,S):
 p=inventories(T);F=(1<<len(T))-1
 K=sum(p['end'][M][a]*p['start'][F^M][b]for M in range(F+1)for a in S for b in S if a!=b)
 return p['H'],sum(p['e'][a]for a in S),sum(p['R'][a]for a in S),K

def profile(T,c):
 p=inventories(T);return p['H'],p['e'][c],p['R'][c],p['allH'][((1<<len(T))-1)^(1<<c)]

def admissible(kind,v,w):
 if kind in TYPE_I:return w!=ONE and (kind not in('A','B')or v!=ONE)
 k=int(kind[1:]);return (v!=ONE and w!=ONE)or(k==0 and v!=ONE and w==ONE)or(k==3 and v==ONE and w!=ONE)

def value(kind,v,w):
 if kind in TYPE_I:
  h,E,R,K=response(kind,v);D,e,U,d=w
  return h*D-E*e+R*U-K*d
 k=int(kind[1:]);return sum(v[i]*MATRICES[k][i][j]*w[j]for i in range(4)for j in range(4))

def amalgamate(V,c,W,a,k):
 # vertex order: shared root, V minus root, cyclic triangle, converse(W) minus root.
 vl=[i for i in range(len(V))if i!=c];wl=[i for i in range(len(W))if i!=a]
 n=len(V)+len(W)+2;out=[0]*n;tri=list(range(1+len(vl),4+len(vl)))
 vmap={c:0,**{x:i+1 for i,x in enumerate(vl)}}
 wmap={a:0,**{x:i+4+len(vl) for i,x in enumerate(wl)}}
 for T,mp in((V,vmap),(converse(W),wmap)):
  for x in mp:
   for y in mp:
    if T[x]>>y&1:out[mp[x]]|=1<<mp[y]
 blocks=([vmap[x]for x in vl],tri,[wmap[x]for x in wl])
 for i in range(3):
  for j in range(i+1,3):
   for x in blocks[i]:
    for y in blocks[j]:out[x]|=1<<y
 for i,x in enumerate(tri):
  out[x]|=1<<tri[(i+1)%3]
  if i<k:out[0]|=1<<x
  else:out[x]|=1
 return out,0

def construct(kind,V,c,W,a):
 if kind in TYPE_I:
  T,S=left_graph(kind,V,c);return graft(T,S,converse(W),a)
 return amalgamate(V,c,W,a,int(kind[1:]))

def realize(record):
 V,c=build_source(record['left_recipe']);W,a=build_source(record['right_recipe'])
 return construct(record['kind'],V,c,W,a)

def support(bound):
 bank=source_bank(bound//3);sources=sorted(bank,key=lambda s:(s[2],s));out={};tested=0
 for kind in TYPE_I+tuple('J'+str(k)for k in range(4)):
  for v in sources:
   for w in sources:
    if not admissible(kind,v,w):continue
    if kind in TYPE_I:
     h,E,R,K=response(kind,v)
     if h*w[2]>bound:break
    elif v!=ONE and w!=ONE and 3*v[2]*w[2]>bound:break
    tested+=1;h=value(kind,v,w)
    if h<=bound:
     out.setdefault(h,dict(kind=kind,left=list(v),right=list(w),left_recipe=list(bank[v]),right_recipe=list(bank[w])))
 return out,bank,tested

def inverse(target):
 # The support engine retains constructive predecessor recipes, not arbitrary formal states.
 return support(target)[0].get(target)
