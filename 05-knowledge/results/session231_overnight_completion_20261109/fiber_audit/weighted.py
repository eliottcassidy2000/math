import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit import literal,check,decode
from physical import ordinal,close,ham,strong
from collections import deque
from pathlib import Path
import json,math
BASE=Path(__file__).resolve().parent.parent/'fibers';OUT=Path(__file__).parent

def mul(x,y):
 h,q,a,g,t,s,e=x;j,k,b,f,u,v,w=y
 return (h*j,q*k,a*k+h*b,f*q+j*g,a*f+t*j+h*u+e*v,s*j,h*w)
def fold(parts):
 v=parts[0]
 for z in parts[1:]:v=mul(v,z)
 return v
one=[(1,2,b,b,0,b,b) for b in (0,1)]
out={}
C5=[8,1,3,6];C9=[16,1,3,7,14];small={
 'moving_H5':[[([0],1),(C5,4),([0],1),([0],1),([0],0),([0],0)],[([0],1),([0],0),([0],1),([0],1),(C5,2),([0],0)]],
 'moving_H9':[[([0],1),(C9,13),([0],0),([0],1),([0],0)],[([0],1),([0],1),(C9,4),([0],1),([0],0)],[([0],1),([0],1),([0],0),(C9,22),([0],0)]]}
for name,graphs in small.items():
 rr=[]
 for parts in graphs:
  a,m=ordinal(parts);s=literal(a)[m];b=close(a,m);count=ham(b);check(strong(b),'not strong');check(count==2*s[2]+s[0]-s[4],'count mismatch');rr.append(dict(state=s,count=count,adjacency=b))
 check(len({tuple(r['state'][:4]+r['state'][5:]) for r in rr})==1,'not fiber');out[name]=rr

receipt=json.loads((BASE/'PHYSICAL_CYCLE.json').read_text());adj=receipt['adj'];mask=receipt['mask'];s=literal(adj)[mask];check(list(s)==receipt['core'],'physical core mismatch');h,q,A,G,T,S,E=s;check(math.gcd(h,q)==1,'saturated')
path=receipt['path'];a,g,p=receipt['start'];start=a,g,p;gain=0;edges=[]
for x,y in path:
 check(x in (0,1) and y in (0,1),'bits');aa=2*a+h*y-q*x;gg2=g+q*y-h*x;check(gg2%2==0,'nonintegral');gg=gg2//2
 check(-h<aa<q and -h<gg<q,'out of bounds');w=a*y+h*p*y-x*gg-h*x*(gg%2);gain+=w;edges.append((a,g,p,x,y,aa,gg,w));a,g,p=aa,gg,y
check((a,g,p)==start,'not closed');check(gain==190,'gain');check(len(path)==312,'length');X=[one[x] for x,y in path];Y=[one[y] for x,y in path];sx=fold(X);sy=fold(Y);u=path[-1][0];v=path[0][1]
left=fold([one[u],s]+Y+[one[v]]);right=fold([one[u]]+X+[s,one[v]])
check(left[:4]+left[5:]==right[:4]+right[5:],'not full interface');delta=left[4]-right[4];check(delta==gain,'cocycle mismatch')
# Evaluate the core-independent cancellation identity, without using local weights.
D=A*sy[3]-sx[2]*G+h*(sy[4]-sx[4])-h*u*sx[5]+h*v*sy[6]
check(D==delta,'boundary identity mismatch')
for k in (1,2,3,6,10):
 values=[];interfaces=set()
 for j in range(k+1):
  value=fold([one[1],one[u]]+X*j+[s]+Y*(k-j)+[one[v]])
  values.append(value[4]);interfaces.add(value[:4]+value[5:])
 check(len(interfaces)==1,'iterated interface');check(all(values[j]-values[j+1]==190 for j in range(k)),'iterated drift')
out['cycle']=dict(core=s,length=len(path),gain=delta,first_x=path[0][0],last_x=u,first_y=v,last_y=path[-1][1],X_state=sx,Y_state=sy,left=left,right=right,independent_global_identity=D,edges=edges,verified_k=[1,2,3,6,10])
# Independently compute the SCC of the cycle's physical node within the whole bounded graph.
# Forward search from this physical root, then reverse-reachable subset.
root=start;graph={};reverse={};todo=deque([root]);seen={root}
while todo:
 a,g,p=todo.popleft();node=(a,g,p);graph[node]=[]
 for x in (0,1):
  for y in (0,1):
   aa=2*a+h*y-q*x;gg2=g+q*y-h*x
   if gg2%2:continue
   gg=gg2//2
   if not(-h<aa<q and -h<gg<q):continue
   dest=(aa,gg,y);graph[node].append(dest);reverse.setdefault(dest,[]).append(node)
   if dest not in seen:seen.add(dest);todo.append(dest)
co={root};todo=deque([root])
while todo:
 for node in reverse.get(todo.popleft(),[]):
  if node not in co:co.add(node);todo.append(node)
check(len(co)==2232,'SCC size');supplied={tuple(v) for c in json.loads((BASE/'NONZERO_SCC.json').read_text()) for v in c};check(co==supplied,'SCC mismatch')
out['cycle']['SCC_size']=len(co)
# All h95 classes through order9: forward endpoint DP computes q and every A/G.
def interfaces(a):
 n=len(a);N=1<<n;F=N-1;en=[[0]*n for _ in range(N)];st=[[0]*n for _ in range(N)];hp=[0]*N;hp[0]=1
 for v in range(n):en[1<<v][v]=st[1<<v][v]=1
 for m in range(1,N):
  hp[m]=sum(en[m])
  for v in range(n):
   if not m>>v&1:continue
   for u in range(n):
    if m>>u&1:continue
    if a[v]>>u&1:en[m|1<<u][u]+=en[m][v]
    if a[u]>>v&1:st[m|1<<u][u]+=st[m][v]
 q=sum(hp[m]*hp[F^m] for m in range(N));aw=[sum(st[m][v]*hp[F^m] for m in range(N)) for v in range(n)];gw=[sum(en[m][v]*hp[F^m] for m in range(N)) for v in range(n)]
 return hp[F],q,[(sum(aw[v] for v in range(n) if m>>v&1),sum(gw[v] for v in range(n) if m>>v&1)) for m in range(N)]
bridge=[]
for n in range(1,10):
 selected=[];hits=[];h95=0
 for line in (BASE/f'classes9_{n}.txt').read_text().splitlines():
  row=line.split()
  if int(row[1])!=95:continue
  code=int(row[0]);h95+=1;H,Q,ag=interfaces(decode(code,n));check(H==95,'H95 mismatch')
  if Q!=1276:continue
  selected.append(code)
  for m,(a,g) in enumerate(ag):
   if (a,g,0) in co or (a,g,1) in co:hits.append((code,m,a,g))
 bridge.append(dict(n=n,H95_classes=h95,H95_Q1276_classes=selected,hits=hits))
out['bridge']=bridge
(OUT/'WEIGHTED.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(length=len(path),gain=delta,SCC=len(co),bridge=bridge),indent=2))
