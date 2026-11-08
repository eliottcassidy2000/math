"""Independent audit. Python 3 standard library only. Assertions are not used."""
import itertools as it, json, math, hashlib
from pathlib import Path

def check(x,msg):
 if not x: raise RuntimeError(msg)
def mm(a,b):return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]
def mv(a,v):return [sum(x*y for x,y in zip(r,v)) for r in a]
def sub(a,b,h):return [[x-h*y for x,y in zip(r,s)] for r,s in zip(a,b)]
def zero(a):return all(not x for r in a for x in r)
def local(adj,out):
 n=len(adj); m=[[0]*6 for _ in range(6)]
 for perm in it.permutations(range(n)):
  for cut in range(n+1):
   p,r=perm[:cut],perm[cut:]
   if any(not(adj[a]>>b&1) for seq in (p,r) for a,b in zip(seq,seq[1:])):continue
   for j in range(6):
    pp,rr=divmod(j,2)
    if not rr and r and r[0] not in out:continue
    np=(2 if p[-1] in out else 1) if p else pp
    nr=int(bool(rr or r));m[2*np+nr][j]+=1
 return m

def stats(adj):
 n=len(adj);full=(1<<n)-1;end=[[0]*n for _ in range(full+1)];h=[1]+[0]*full
 for mask in range(1,full+1):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    end[mask][v]=sum(end[rest][u] for u in range(n) if adj[u]>>v&1) if rest else 1
  h[mask]=sum(end[mask])
 # Compute first endpoint counts independently by reversing arcs.
 rev=[sum(1<<j for j in range(n) if adj[j]>>i&1) for i in range(n)]
 first=[[0]*n for _ in range(full+1)]
 for mask in range(1,full+1):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v);first[mask][v]=sum(first[rest][u] for u in range(n) if adj[v]>>u&1) if rest else 1
 alpha=[sum(first[s][v]*h[full^s] for s in range(1,full+1)) for v in range(n)]
 return h, end, first, alpha

def module_stats(adj):
 h,e,f,a=stats(adj);full=len(h)-1
 return h[-1],sum(h[s]*h[full^s] for s in range(full+1)),a
S=[local([0],set() if b==0 else {0}) for b in range(2)]
tri=[2,4,1];adj=[18,8,3,5,14]
H,Q,alpha=module_stats(adj)
check((H,Q,alpha)==(13,104,[18,14,18,15,26]),'module')
F=[sub(local(tri,set(range(k))),mm(S[k%2],S[k//2]),3) for k in range(4)]
G=[sub(local(adj,out),mm(S[b],mm(S[1-b],S[b])),13) for b,out in [(0,{4}),(1,{0,1,2,3})]]
W=[[0,0],[1,0],[0,0],[-1,1],[0,0],[0,-1]]
def in_w(a):
 return all(a[i][j]==0 for i in (0,2,4) for j in range(len(a[0]))) and all(a[1][j]+a[3][j]+a[5][j]==0 for j in range(len(a[0])))
# Universal category-difference basis, in order x,y,a,b.
U=[]
for t in range(4):
 x,y,a,b=[int(t==i) for i in range(4)]
 weights={(0,1):x,(0,2):-x,(1,0):y,(2,0):-y,(1,1):a,(1,2):b,(2,1):-x-a,(2,2):x-b}
 C=[[0]*6 for _ in range(6)]
 for j in range(6):
  p,r=divmod(j,2)
  for (pp,rr),value in weights.items():
   if not r and rr==2:continue
   C[2*(pp or p)+int(bool(r or rr))][j]+=value
 expected_C=[[0]*6,[x,0,0,0,0,0],[y,0,y,0,y,0],
             [a,a+b+y,a+x,a+b+y,a,a+b+y],[-y,0,-y,0,-y,0],
             [-a-x,-a-b-y,-a-x,-a-b-y,-a,-a-b-y]]
 check(C==expected_C,'UNIVERSAL.md displayed family')
 U.append(C)
for D in U:
 check(zero(mm(D,W)),'universal annihilation')
 for B in S:check(in_w(mm(B,D)),'universal entry')
for B in S:check(in_w(mm(B,W)),'invariance')
for D in F+G:
 check(zero(mm(D,W)),'annihilation')
 for B in S:check(in_w(mm(B,D)),'entry to W')
A=mm(S[1],S[0]);powers=[[[int(i==j) for j in range(6)] for i in range(6)]]
for i in range(4):powers.append(mm(A,powers[-1]))
check(zero([[powers[4][i][j]-6*powers[3][i][j]+9*powers[2][i][j]-4*powers[1][i][j] for j in range(6)] for i in range(6)]),'recurrence')
v=[1,0,0,0,0,0];ell=[1,1,1,1,0,0]
for r in range(1,5):check(sum(a*b for a,b in zip(ell,mv(powers[r],v)))==(4**r+6*r-1)//9,'baseline grid')
# G has rank exactly three: all rows span rows 1,2,3, with a nonzero minor.
g=G[0]
check(all(g[4][j]==-g[2][j] and g[5][j]==-g[3][j]-g[1][j] for j in range(6)) and not any(g[0]),'rank upper bound')
a,b,c=[g[i][:3] for i in (1,2,3)]
check(a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])!=0,'rank lower bound')
def scalar(x):return sum(a*b for a,b in zip(ell,x))
grid=[]
for j in range(4):
 row=[]
 for k in range(4):
  value=scalar(mv(powers[k],mv(S[1],mv(G[0],mv(powers[j],v)))))
  check(value==4**(j+1)-10+4*(j==0),'defect grid');row.append(value)
 grid.append(row)
# literal tournament construction, independent of local transfer

def graph(r,p,trim=False):
 bits='01'*r;blocks=[];outs=[];i=0
 while i<len(bits):
  if i==p:blocks.append(adj);outs.append({4});i+=3
  else:blocks.append([0]);outs.append({0} if bits[i]=='1' else set());i+=1
 n=sum(map(len,blocks));a=[0]*(n+1);off=0
 for block,out in zip(blocks,outs):
  for u,mask in enumerate(block):
   a[off+u]=(mask<<off)|sum(1<<j for j in range(off+len(block),n))
   if u in out:a[n]|=1<<(off+u)
   else:a[off+u]|=1<<n
  off+=len(block)
 if trim:
  keep=list(range(1,n-1))+[n]
  a=[sum(1<<j for j,w in enumerate(keep) if a[u]>>w&1) for u in keep]
 return a

def state(a):
 n=len(a);h,e,f,alpha=stats(a);root=n-1;small=a[:-1];h2,_,_,al=stats(small)
 return [h[-1],e[-1][root],h2[-1]+sum(al[i] for i in range(root) if a[root]>>i&1),h2[-1]]
def strong(a):
 for start in range(len(a)):
  reach=1<<start
  while True:
   new=reach
   for i in range(len(a)):
    if reach>>i&1:new|=a[i]
   if new==reach:break
   reach=new
  if reach!=(1<<len(a))-1:return False
 return True
witnesses=[]
for r,js,trim in [(4,range(3),False),(7,range(1,5),True)]:
 for j in js:
  a=graph(r,2*j,trim);x=state(a)
  d=13*(4**r+6*r-1)//9+4**(j+1)-10+4*(j==0)
  expected=[d,13 if trim else 0,13*(2**(2*r-1)+1)//3 if trim else 13*(4**r+2)//3,13]
  check(x==expected,'literal state');check(not trim or strong(a),'strong')
  witnesses.append({'r':r,'j':j,'trim':trim,'adj':a,'state':x})
check(math.gcd(*[w['state'][0] for w in witnesses[:3]])==1,'gcd small')
check(math.gcd(*[w['state'][0] for w in witnesses[3:]])==1,'gcd strong')
census={}
for n in range(1,7):
 counts={};signatures=set();edges=list(it.combinations(range(n),2))
 for code in range(1<<len(edges)):
  a=[0]*n
  for i,(u,vv) in enumerate(edges):
   if code>>i&1:a[u]|=1<<vv
   else:a[vv]|=1<<u
  h,q,al=module_stats(a)
  if q%h==0 and q//h & (q//h-1)==0:
   key=f'{h},{q//h}';counts[key]=counts.get(key,0)+1
   sol=tuple((mask,sum(al[v] for v in range(n) if mask>>v&1)//h) for mask in range(1,1<<n) if sum(al[v] for v in range(n) if mask>>v&1)%h==0 and 0<sum(al[v] for v in range(n) if mask>>v&1)//h<q//h-1)
   signatures.add((h,q//h,tuple(sorted(al)),tuple(sorted(set(v for _,v in sol)))))
 bank=json.loads((Path(__file__).parent/'inputs'/f'bank{n}.json').read_text())
 bank_signatures=set()
 for rec in bank:
  hh,qq,aa=module_stats(rec['adj'])
  check((hh,qq,aa)==(rec['H'],rec['Q'],rec['alpha']),'bank statistics')
  ss=[[mask,sum(aa[v] for v in range(n) if mask>>v&1)//hh] for mask in range(1,1<<n) if sum(aa[v] for v in range(n) if mask>>v&1)%hh==0 and 0<sum(aa[v] for v in range(n) if mask>>v&1)//hh<qq//hh-1]
  check(ss==rec['sol'],'bank subsets')
  bank_signatures.add((hh,qq//hh,tuple(sorted(aa)),tuple(sorted(set(v for _,v in ss)))))
 check(len(bank)==len(bank_signatures) and signatures==bank_signatures,'bank completeness')
 census[str(n)]=counts
expected={'1':{'1,2':1},'2':{'1,4':2},'3':{'1,8':6,'3,4':2},'4':{'1,16':24,'3,8':16},'5':{'1,32':120,'3,16':120,'13,8':120,'15,8':24},'6':{'1,64':720,'3,32':960,'9,16':80,'13,16':1440,'15,16':288}}
check(census==expected,'census')
print(json.dumps({'passed':True,'module':{'H':H,'Q':Q,'alpha':alpha},'G':G,'universal_basis':U,'input_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((Path(__file__).parent/'inputs').iterdir())},'invariant_basis':W,'defect_grid':grid,'literal_witnesses':witnesses,'census':census,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
