"""Independent audit: no imports from producer or its verifier."""
import json,hashlib,itertools
from pathlib import Path
import sympy as sp
BASE=Path(__file__).resolve().parents[1]
frontier_path=BASE/'math/05-knowledge/results/session18_native_annulus/frontier.json'
V=set(json.loads(frontier_path.read_text())['frontier'])
checks=0

def require(p,msg):
 global checks
 checks+=1
 if not p:raise RuntimeError(msg)

def Q(w):return w[0]**2-w[0]*w[1]+w[1]**2

def add(*ws):return tuple(sum(w[j] for w in ws) for j in range(2))
def scale(k,w):return tuple(k*x for x in w)
def beta(a,b):return Q(add(a,b))-Q(a)-Q(b)
def increment(w,h):return Q(add(w,h))-Q(w)
L0=(1,2);L2=(1,-1);U1=(-1,-1)
x,y,a,b=sp.symbols('x y a b',real=True)
w=add(scale(x+a,L0),scale(y+b,L2));h=add(scale(x,U1),scale(-y,L2))
F=-2*x*x+3*x*y-3*y*y-3*a*x+3*a*y-6*b*y
require(sp.expand(increment(w,h)-F)==0,'boundary identity')
M=sp.Matrix([[2,-sp.Rational(3,2)],[-sp.Rational(3,2),3]])
l=sp.Matrix([-3*a,3*a-6*b]);bound=(l.T*M.inv()*l)[0]/4
require(sp.expand(bound-sp.Rational(6,5)*(a*a-a*b+4*b*b))==0,'boundary Schur')
require(M.det()==sp.Rational(15,4) and M[0,0]>0,'positive M')
u,v=sp.symbols('u v',real=True)
require(sp.expand(4*Q((u,v))-9*Q((0,(u+v)/3))-3*(u-v)**2)==0,'sharp contraction')
ox,oy,hx,hy=sp.symbols('ox oy hx hy',real=True)
o=(ox,oy);ex=(hx,hy);B=add(scale(a,L0),scale(b,L2));cw=add(o,B)
l1=beta(ex,U1)+beta(cw,U1)+beta(L0,ex)
l2=beta(ex,scale(-1,L2))+beta(cw,scale(-1,L2))+beta(L2,ex)
full=increment(add(cw,scale(x,L0),scale(y,L2)),add(ex,scale(x,U1),scale(-y,L2)))
require(sp.expand(full-(increment(o,ex)+beta(B,ex)-2*x*x+3*x*y-3*y*y+l1*x+l2*y))==0,'exterior exact')
lv=sp.Matrix([l1,l2]);schur=(lv.T*M.inv()*lv)[0]/4
require(sp.expand(schur-(3*l1*l1+3*l1*l2+2*l2*l2)/15)==0,'exterior Schur')

def primitive(s):return s>0 and s%2==1 and s%8!=5

def P(s):
 require(primitive(s),'source primitive')
 clock=2 if s%8==1 else 1
 t=(3*s+1)//(1<<clock)
 while t%8==5:t=(t-1)//4;clock+=2
 require(primitive(t),'target primitive')
 return t,clock

def charge(N,s):
 vectors=[(1,2),(-2,-1),(1,-1)]
 uppers=[(0,1),(-1,-1),(1,0)]
 acc=[0,0]
 while s<4*N:
  vec=(vectors if s<2*N else uppers)[s%3]
  acc[0]+=vec[0];acc[1]+=vec[1];s=4*s+1
 return tuple(acc)

def family(N):
 cut=(3*N+1)//2
 for tag,lo,hi,res in [('A',cut,2*N,{9,27,63}),('B',N,cut,{17,35,71}),('C',cut,2*N,{3,39,57}),('D',cut,2*N,{23,41,59})]:
  for r in res:
   for parent in range(lo+(r-lo)%72,hi,72):
    if tag=='A':numerators=[(16*parent+3,4)]
    elif tag=='B':numerators=[(64*parent+19,6)]
    elif tag=='C':numerators=[(32*parent+9,5),(16*parent+3,4)]
    else:numerators=[(128*parent+41,7)]
    require(all(num%3==0 for num,_ in numerators),'integrality')
    yield tag,parent,[(num//3,clock) for num,clock in numerators]

def arrival(s):
 clock=0
 for depth in range(501):
  if s in V:return (s,clock,depth)
  if depth==500:return None
  s,c=P(s);clock+=c

digest=hashlib.sha256();rows=children=arrivals=0;summary=[]
require(len(V)==627 and max(V)==63665,'frontier')
for N in list(range(1,513))+[144*10,4096,7958,7959,8192,63666,65536,100000]:
 seen=set();counts={x:0 for x in 'ABCD'}
 for tag,p,cs in family(N):
  counts[tag]+=1;rows+=1
  ids=[p]+[s for s,c in cs]
  require(len(set(ids))==len(ids) and not seen.intersection(ids),'original source unique');seen.update(ids)
  require(primitive(p) and N<=p<2*N,'parent carrier')
  require(charge(N,p)==charge(64*N,p),'old source unchanged')
  h=(0,0)
  for s,c in cs:
   children+=1
   require(P(s)==(p,c),'literal native edge and full clock')
   require(4*N<=s<256*N and s>8*N,'new source carrier')
   if N>=7959:require(s not in V,'stopping guard')
   h=add(h,scale((-1)**c,charge(64*N,s)))
   digest.update(f'{N}:{tag}:{p}:{s}:{c}:{p%15}:{p%35}:{s%15}:{s%35};'.encode())
  target=(L0,U1) if tag in ('A','C') else (L2,scale(-1,L2))
  require((charge(N,p),h)==target,'literal vector block')
  if N==63666:
   rp=arrival(p);require(rp is not None,'finite control parent')
   for s,c in cs:
    require(arrival(s)==(rp[0],rp[1]+c,rp[2]+1),'actual synchronized arrival');arrivals+=1
 for tag,lo,hi,res in [('A',(3*N+1)//2,2*N,{9,27,63}),('B',N,(3*N+1)//2,{17,35,71}),('C',(3*N+1)//2,2*N,{3,39,57}),('D',(3*N+1)//2,2*N,{23,41,59})]:
  expected=sum((hi-1-r)//72-(lo-1-r)//72 for r in res)
  require(counts[tag]==expected,'exact residue count')
 if N%144==0:require(all(c==N//48 for c in counts.values()),'aligned count')
 if N>=4096:summary.append({'N':N,'clusters':counts,'distinct_sources':len(seen)})
# All-scale symbolic linear substitutions B^j(u), rather than only arithmetic samples.
r=sp.symbols('r',integer=True)
for multiplier,offset,initial,j in [(16,3,2,1),(64,19,2,2),(32,9,1,2),(128,41,1,3)]:
 child=(multiplier*r+offset)/3
 post=(3*child+1)/(2**initial)
 Bj=(4**j)*r+sp.Rational(4**j-1,3)
 require(sp.expand(post-Bj)==0,'exact inverse chain identity')
require(increment(add((-10,0),L2),scale(-1,L2))==27,'hostile exterior full increment')
require(increment(L2,scale(-1,L2))==-3,'hostile isolated pair increment')
require(Q((0,2))*9==Q((3,3))*4,'sharp attainable bank lattice equality')
out={'status':'PASS','symbolic_identities':11,'checks':checks,'literal_clusters':rows,'literal_children':children,'actual_first_arrival_child_checks':arrivals,'source_clock_register_digest':digest.hexdigest(),'frontier_sha256':hashlib.sha256(frontier_path.read_bytes()).hexdigest(),'selected_scales':summary,'uniform_153N':'NOT PROVED','producer_imports':False}
print(json.dumps(out,indent=2,sort_keys=True))
