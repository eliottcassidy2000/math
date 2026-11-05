import importlib.util,hashlib,itertools,json,math
from pathlib import Path
p=Path(__file__).resolve().parent.parent/'ribbon-terminals'/'verify.py'
if hashlib.sha256(p.read_bytes()).hexdigest()!='9ada1cc1cfaf33030a89122dcf4ed65103aa23f5a2ad01e5dd5c288cde198a3e':raise RuntimeError('dependency changed')
spec=importlib.util.spec_from_file_location('terminals',p);t=importlib.util.module_from_spec(spec);spec.loader.exec_module(t)
A=[[5,-1,7],[1,0,1],[7,-1,14]];D=[[1,0,4],[0,0,1],[4,-1,6]]
cells=[([8,5,1,6],0,3,'2213',A,[3,1,2,0]),([2,0,3],1,2,'013',D,[2,0,1])];seed=[4,1,2];initial=(3,1,3)
need=t.require

def valid(h,e,q):return h>=3 and h%2==1 and q in(-2,2) and 1<=e<=h and e<=h+q

def decode(h,e,q):
 out=[]
 while (h,e,q)!=(3,1,-2):
  if not valid(h,e,q):return None
  if 7*e>h and 5*e<h:
   a=0;nn=(-h+8*e-q,-3*h+19*e+2*q);dd=(5,5);qq=q
  elif 5*e>=h and 3*e<h:
   a=1;nn=(h-4*e,2*h-9*e-q);dd=(1,1);qq=-q
  else:return None
  if any(n%d for n,d in zip(nn,dd)):return None
  hh,ee=(n//d for n,d in zip(nn,dd))
  if not valid(hh,ee,qq) or hh>=h:return None
  out.append(a);h,e,q=hh,ee,qq
 return tuple(reversed(out))

for i,(G,b,c,w,M,path) in enumerate(cells):
 need(t.matrix(G,b,c)==M,('matrix',i));need(t.span(G)<=2 and b<2 and c==len(G)-1,('ordered cell',i))
 need(sorted(path)==list(range(len(G))) and path[-1]==b and all(G[x]>>y&1 for x,y in zip(path,path[1:])),('entry-ending path',i))
need(t.rm([-2,1,1],A)==[-2,1,1] and t.rm([-2,1,1],D)==[2,-1,-1],'charge rows')
need(t.mm(D,D)==[[17,-4,28],[4,-1,6],[28,-6,51]],'squared flip')
physical=[]
for n in range(4):
 for word in itertools.product(range(2),repeat=n):
  rows=seed;v=initial;mask='022'
  for a in word:
   G,b,c,w,M,path=cells[a];rows=t.glue(rows,len(rows)-1,G,b);v=t.mv(M,v);mask+=w
  need(rows==t.wordgraph(mask,2),('literal word',word));need(t.span(rows)<=2,('literal span',word))
  need(t.profile(rows,len(rows)-1)==v,('literal profile',word))
  physical.append(dict(cells=''.join('AD'[x] for x in word),order=len(rows),profile=v,word=mask))
frontier=[((),initial)];states={};digest=hashlib.sha256()
for depth in range(11):
 nxt=[]
 for word,v in frontier:
  h,e,R=v;q=-2*h+e+R
  need(valid(h,e,q) and h%4==3,('physical cone',word,v))
  need(decode(h,e,q)==word,('decode',word,v));need(v not in states,('collision',word));states[v]=word
  ha=t.mv(A,v)[0];haa=t.mv(A,t.mv(A,v))[0];hdda=t.mv(A,t.mv(D,t.mv(D,v)))[0]
  need(haa-12*ha-13*h==48*R and hdda-61*ha+28*h==64*R,'content identities')
  need(math.gcd(h,ha,haa,hdda)==1,'observable content')
  digest.update(str((word,v)).encode())
  if depth<10:
   for a,cell in enumerate(cells):
    nv=t.mv(cell[4],v);need(nv[0]>=[11,5][a]*h,'growth');nxt.append((word+(a,),nv))
 frontier=nxt
print(json.dumps(dict(physical_graph_checks=physical,decoded_words=len(states),state_data_sha256=digest.hexdigest(),first_D_state=t.mv(D,initial),second_D_state=t.mv(D,t.mv(D,initial))),sort_keys=True,indent=2))
