"""Standalone stdlib-only audit; run next to UNIT_WITNESS.json."""
import itertools,json,math,hashlib
from pathlib import Path
ROOT=Path(__file__).parent

def require(b,message):
 if not b:raise RuntimeError(message)

def literal_core(adj,mask):
 n=len(adj);require(all(not(adj[v]>>v&1) and adj[v]<1<<n for v in range(n)),'invalid row');require(all(((adj[u]>>v)^(adj[v]>>u))&1 for u in range(n) for v in range(u)),'invalid tournament')
 h=q=A=G=T=S=E=0
 for p in itertools.permutations(range(n)):
  bad=[i+1 for i in range(n-1) if not(adj[p[i]]>>p[i+1]&1)]
  if not bad:
   h+=1;S+=(mask>>p[0])&1;E+=(mask>>p[-1])&1;splits=range(n+1)
  elif len(bad)==1:splits=bad
  else:continue
  for k in splits:
   q+=1;a=k<n and bool(mask>>p[k]&1);g=k>0 and bool(mask>>p[k-1]&1);A+=a;G+=g;T+=a and g
 return h,q,A,G,T,S,E

def multiply(z,w):
 h,q,a,g,t,s,e=z;j,r,b,f,v,u,l=w
 return h*j,q*r,a*r+h*b,f*q+j*g,a*f+t*j+h*v+e*u,s*j,h*l

def fold(parts):
 z=parts[0]
 for w in parts[1:]:z=multiply(z,w)
 return z

def word(bits):
 require(bool(bits),'empty word has no endpoint state')
 # Direct singleton-word formulas, independently of multiply().
 # T sums selected ordered pairs with gap-weight 2^(j-i-1),
 # with one additional unit for each selected adjacent pair.
 a=int(bits,2);g=int(bits[::-1],2);prefix=0;prev=0;t=0
 for c in bits:
  b=int(c);t+=b*(prefix+prev);prefix=2*prefix+b;prev=b
 return 1,1<<len(bits),a,g,t,int(bits[0]),int(bits[-1])

w=json.loads((ROOT/'UNIT_WITNESS.json').read_text());adj=w['adj'];mask=w['mask'];core=literal_core(adj,mask);require(list(core)==w['core'],'literal core');h,q,A,G,T,S,E=core
# Verify the class encoding independently.
code=0;bit=0
for i in range(len(adj)):
 for j in range(i+1,len(adj)):
  code|=((adj[i]>>j)&1)<<bit;bit+=1
require(code==1147012 and mask==30,'class code/cut');require(math.gcd(h,q)==37,'saturation')
X,Y=w['X'],w['Y'];require(len(X)==len(Y)==1188 and set(X+Y)<={'0','1'},'binary words')
start=(A,G,int(Y[-1]));a,g,p=start;gain=0;segments=[];trace=[]
for i,(cx,cy) in enumerate(zip(X,Y)):
 x,y=int(cx),int(cy);aa=2*a+h*y-q*x;two_g=g+q*y-h*x;require(two_g%2==0,'nonintegral forward edge');gg=two_g//2;require(-h<aa<q and -h<gg<q,'out-of-bounds forward edge')
 wt=a*y+h*p*y-x*gg-h*x*(gg%2);gain+=wt;trace.append((a,g,p,x,y,aa,gg,wt));a,g,p=aa,gg,y
 if (i+1)%108==0:
  require((a,g,p)==start,'108-edge segment not rooted closed');segments.append(sum(z[-1] for z in trace[-108:]))
require(segments==[34]*8+[-90]*3,'directed loop concatenation');require((a,g,p)==start and gain==2,'closed total gain')
require(X[:108]*8+X[864:972]*3==X and Y[:108]*8+Y[864:972]*3==Y,'not literal eight-plus-three repetition')

one=[(1,2,b,b,0,b,b) for b in (0,1)];sx,sy=word(X),word(Y);u,v=int(X[-1]),int(Y[0]);require((u,v)==(0,0),'boundary bits')
left=fold([one[u],core,sy,one[v]]);right=fold([one[u],sx,core,one[v]])
require(left[:4]+left[5:]==right[:4]+right[5:],'flanked full interface');require(left[4]-right[4]==2,'actual flanked unit gap')
# Independent global boundary identity. It does not use any local edge weights.
global_delta=A*sy[3]-sx[2]*G+h*(sy[4]-sx[4])-h*u*sx[5]+h*v*sy[6]
require(global_delta==2,'global word cancellation')
for pos,wanted in ((0,34),(864,-90)):
 xx,yy=word(X[pos:pos+108]),word(Y[pos:pos+108]);require(A*yy[3]-xx[2]*G+h*(yy[4]-xx[4])==wanted,'global constituent loop gain')
checks=[]
for k in range(9):
 vals=[];ifaces=set()
 for j in range(k+1):
  parts=[one[1],one[0]]
  if j:parts.append(word(X*j))
  parts.append(core)
  if k-j:parts.append(word(Y*(k-j)))
  parts.append(one[0]);z=fold(parts);ifaces.add(z[:4]+z[5:]);vals.append(2*z[2]+z[0]-z[4])
 require(len(ifaces)==1,'iterated full interface');require(all(n==vals[0]+2*j for j,n in enumerate(vals)),'not consecutive');require(vals[0]>0 and vals[0]%2==1,'not positive odd')
 yy=word(Y*k) if k else None
 # Closed baseline formula uses only one repeated singleton word.
 baseline=6996*(1<<(1188*k))-939
 if yy:baseline+=148*yy[2]-1749*yy[3]-37*yy[4]
 require(baseline==vals[0],'baseline formula');qtotal=740*(1<<(1188*k+3));require(math.gcd(37,qtotal)==37,'iterated saturation')
 checks.append(dict(k=k,star_order=1188*k+12,size=k+1,baseline_bitlength=baseline.bit_length(),baseline_mod10000=baseline%10000,last_mod10000=vals[-1]%10000))
# Literal physical strongness control at k=0; all k follows from selected first singleton.
def ordinal(parts):
 n=sum(len(a) for a,m in parts);out=[];cut=0;offset=0
 for a,m in parts:
  later=((1<<n)-1)^((1<<(offset+len(a)))-1)
  out.extend((r<<offset)|later for r in a);cut|=m<<offset;offset+=len(a)
 return out,cut

a,cut=ordinal([([0],1),([0],0),(adj,mask),([0],0)]);n=len(a);star=[row|(1<<n)|((not(cut>>i&1))<<(n+1)) for i,row in enumerate(a)]+[1<<(n+1),cut]
for v0 in range(len(star)):
 reached=1<<v0
 while True:
  nex=reached
  for i,row in enumerate(star):
   if reached>>i&1:nex|=row
  if nex==reached:break
  reached=nex
 require(reached==(1<<len(star))-1,'literal nonstrong base')
N=1<<len(star);dp=[{} for _ in range(N)]
for v in range(len(star)):dp[1<<v][v]=1
for m in range(1,N):
 for v,count in dp[m].items():
  avail=star[v]&~m
  while avail:
   b=avail&-avail;avail-=b;u=b.bit_length()-1;dest=dp[m|b];dest[u]=dest.get(u,0)+count
require(sum(dp[-1].values())==6057,'literal base Hamilton count')
result=dict(status='PASS',core=core,class_code=code,mask=mask,length=len(X),root=start,rooted_segment_lengths=[108]*11,rooted_segment_gains=segments,total_gain=gain,global_flanked_gap=global_delta,saturation_gcd=37,first_x=int(X[0]),last_x=int(X[-1]),first_y=int(Y[0]),last_y=int(Y[-1]),checks=checks,base_star_adjacency=star,base_direct_H=6057,trace_sha256=hashlib.sha256(repr(trace).encode()).hexdigest())
(ROOT/'UNIT_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
