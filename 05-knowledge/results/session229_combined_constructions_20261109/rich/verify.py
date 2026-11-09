"""Portable finite-exact checks; no repository inputs, assertions, or dependencies."""
import itertools as it,math,json
B=[26,44,57,16,34,9];MS=[9,10,17,18]
def require(c,s):
 if not c:raise RuntimeError(s)
def categories(a,m):
 w={(p,r):0 for p in range(3) for r in range(3) if (p,r)!=(0,0)}
 for seq in it.permutations(range(len(a))):
  for k in range(len(a)+1):
   p,r=seq[:k],seq[k:]
   if any(not(a[u]>>v&1) for z in (p,r) for u,v in zip(z,z[1:])):continue
   pp=(2 if m>>p[-1]&1 else 1) if p else 0;rr=(1 if m>>r[0]&1 else 2) if r else 0
   w[pp,rr]+=1
 return w
def coords(w):
 H=w[0,1]+w[0,2];Q=sum(w.values());A=w[0,1]+w[1,1]+w[2,1];G=w[2,0]+w[2,1]+w[2,2]
 return H,Q,A,G,w[0,1],w[2,0],w[2,1]
def graph(m):
 # ordered p,B,q and final root
 a=[0]*9
 for i in range(8):a[i]|=sum(1<<j for j in range(i+1,8))
 for i,w in enumerate(B):a[i+1]=(w<<1)|(1<<7)
 for i in range(8):
  if i==0 or 1<=i<=6 and m>>(i-1)&1:a[8]|=1<<i
  else:a[i]|=1<<8
 return a
def ham(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for u,c in enumerate(dp[s]):
   if c:
    for v in range(n):
     if a[u]>>v&1 and not s>>v&1:dp[s|1<<v][v]+=c
 return sum(dp[-1])
cs=[coords(categories(B,m)) for m in MS]
expected=[(29,320,81,129,7,15,22),(29,320,81,117,7,13,16),(29,320,81,93,7,9,12),(29,320,81,81,7,7,14)]
require(cs==expected,'full category enumeration');require(math.gcd(29,320)==1,'unsaturated')
aw=[coords(categories(B,1<<i))[2] for i in range(6)];require(aw==[50,50,86,31,31,43],'alpha weights')
fiber=[m for m in range(64) if sum(a for i,a in enumerate(aw) if m>>i&1)==81];require(fiber==MS,'complete alpha81 fiber')
ds=[ham(graph(m)) for m in MS];require(ds==[673,691,719,729],'independent literal nine-vertex DP')
for c,d in zip(cs,ds):
 H,Q,A,G,S,E,T=c;require(d==2*Q+H+2*A-G-S-T,'strong flank scalar formula')
# Exhaust all 4^k words through k=7 for G injectivity.
levels=[];gs={0}
for k in range(1,8):
 gs={29*g+320**(k-1)*c[3] for g in gs for c in cs}
 require(len(gs)==4**k,'G word injectivity');levels.append({'length':k,'distinct_G':len(gs)})
# Aut(B) trivial, independently by all 720 vertex permutations.
autos=[]
for p in it.permutations(range(6)):
 if all(bool(B[i]>>j&1)==bool(B[p[i]]>>p[j]&1) for i in range(6) for j in range(6)):autos.append(p)
require(len(autos)==1,'rigidity')
rich_result=dict(block=B,alpha_weights=aw,coordinates=cs,strong_scalars=ds,amplification=levels,automorphism_order=len(autos))
# Unsaturated interface-preserving gap-two switches, independently enumerated.
SW=[([64,121,35,37,13,81,28],[12,20],[(95,1276,312,270,24,18,38),(95,1276,312,270,24,18,36)]),([20,85,72,195,204,159,33,71],[47,181],[(389,5600,3369,3456,250,271,1858),(389,5600,3369,3456,250,271,1860)])]
def star(a,m):
 n=len(a);z=[w|(1<<n) for w in a]+[1<<(n+1),0]
 for i in range(n):
  if m>>i&1:z[n+1]|=1<<i
  else:z[i]|=1<<(n+1)
 return z
def compose(c,d):
 h,q,A,G,S,E,T=c;k,r,B,J,V,F,U=d
 return h*k,q*r,A*r+h*B,J*q+k*G,S*k,h*F,A*J+T*k+h*U+E*V
def transfer(w):
 M=[[0]*6 for _ in range(6)]
 for j in range(6):
  p,r=divmod(j,2)
  for (pp,rr),v in w.items():
   if not r and rr==2:continue
   M[2*(pp or p)+int(bool(r or rr))][j]+=v
 return M
def strong(a):
 for start in range(len(a)):
  seen=1<<start
  while True:
   nxt=seen
   for i in range(len(a)):
    if seen>>i&1:nxt|=a[i]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<len(a))-1:return False
 return True
switch_results=[]
for a,ms,expected in SW:
 ws=[categories(a,m) for m in ms];cs=[coords(w) for w in ws];require(cs==expected,'independent unsaturated interface switch')
 require(strong(a),'strong switch core')
 mats=[transfer(w) for w in ws];delta=cs[1][-1]-cs[0][-1];w=[0,0,0,1,0,-1];rho=[1,0,1,0,1,0]
 require(all(mats[1][i][j]-mats[0][i][j]==-delta*w[i]*rho[j] for i in range(6) for j in range(6)),'central defect formula')
 require(math.gcd(cs[0][0],cs[0][1])==1,'switch coprimality')
 require(cs[0][:-1]==cs[1][:-1] and abs(cs[0][-1]-cs[1][-1])==2,'full interface equality')
 ds=[ham(star(a,m)) for m in ms];require(ds==[2*c[2]+c[0]-c[-1] for c in cs],'literal star closure DP')
 # All 16 length-four singleton mask contexts, split before/after at every cut.
 for bits in it.product((0,1),repeat=4):
  for cut in range(5):
   out=[]
   for c in cs:
    for b in bits[:cut][::-1]:c=compose((1,2,b,b,b,b,0),c)
    for b in bits[cut:]:c=compose(c,(1,2,b,b,b,b,0))
    out.append(c)
   require(out[0][:-1]==out[1][:-1] and out[0][-1]-out[1][-1]==cs[0][-1]-cs[1][-1],'all short singleton contexts')
 switch_results.append(dict(block=a,masks=ms,coordinates=cs,star_scalars=ds))
print(json.dumps(dict(status='PASS',four_choice_module=rich_result,unsaturated_switches=switch_results),indent=2))
