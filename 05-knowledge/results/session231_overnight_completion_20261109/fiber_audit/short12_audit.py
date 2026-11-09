"""Standalone stdlib-only audit; run next to SHORT12_WITNESS.json."""
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

w=json.loads((ROOT/'SHORT12_WITNESS.json').read_text());core=literal_core(w['adj'],w['mask']);require(list(core)==w['core'],'literal core');h,q,A,G,T,S,E=core;X,Y=w['X'],w['Y'];require((X,Y)==('101000111000','111011100110'),'words');require(math.gcd(h,q)==45,'saturation');sx,sy=word(X),word(Y);require(sx==(1,4096,2616,453,288,1,0),'X state');require(sy==(1,4096,3814,1655,1540,1,0),'Y state')
a,g,p=A,G,int(Y[-1]);root=a,g,p;trace=[];gain=0
for cx,cy in zip(X,Y):
 x,y=int(cx),int(cy);aa=2*a+h*y-q*x;gg2=g+q*y-h*x;require(gg2%2==0,'nonintegral');gg=gg2//2;require(-h<aa<q and -h<gg<q,'out of bounds');wt=a*y+h*p*y-x*gg-h*x*(gg%2);trace.append([a,g,p,x,y,aa,gg,wt]);gain+=wt;a,g,p=aa,gg,y
require((a,g,p)==root and gain==-2,'closed loop gain')
one=[(1,2,b,b,0,b,b) for b in (0,1)];u=int(X[-1]);v=int(Y[0]);left=fold([one[u],core,sy,one[v]]);right=fold([one[u],sx,core,one[v]]);require(left[:4]+left[5:]==right[:4]+right[5:],'flanked interface');require(left[4]-right[4]==-2,'actual unit identity')
global_delta=A*sy[3]-sx[2]*G+h*(sy[4]-sx[4])-h*u*sx[5]+h*v*sy[6];require(global_delta==-2,'global word cancellation')
checks=[];old_baseline=None
for k in list(range(9))+[20,100]:
 vals=[];ifaces=set()
 for j in range(k+1):
  ps=[one[1],one[0]]
  if j:ps.append(word(X*j))
  ps.append(core)
  if k-j:ps.append(word(Y*(k-j)))
  ps.append(one[1]);z=fold(ps);ifaces.add(z[:4]+z[5:]);vals.append(2*z[2]+h-z[4])
 require(len(ifaces)==1,'full interface');require(all(c==vals[0]-2*j for j,c in enumerate(vals)),'interval');require(all(c>0 and c%2 for c in vals),'positive odds');require(next(iter(ifaces))[1]==1440*(1<<(12*k+3)),'q formula')
 yy=word(Y*k) if k else None;baseline=11274*(1<<(12*k))-1365
 if yy:baseline+=135*yy[2]-3758*yy[3]-45*yy[4]
 require(baseline==vals[0],'baseline formula')
 closed_num=245051264*(4096**k)+1602510*k+1118023
 require(closed_num%24843==0 and closed_num//24843==baseline,'closed baseline')
 if old_baseline is not None and k<=8:require(baseline==4096*old_baseline-264150*(k-1)-184225,'baseline recurrence')
 old_baseline=baseline
 checks.append(dict(k=k,vertices=12*k+13,size=k+1,baseline=baseline,minimum=vals[-1],gcd=math.gcd(h,1440*(1<<(12*k+3)))))
# Direct actual 13-vertex closure for k=0: full forward DP and reachability.
adj=w['adj'];mask=w['mask'];parts=[([0],1),([0],0),(adj,mask),([0],1)];n=sum(len(a) for a,m in parts);a=[];cut=0;offset=0
for rows,m in parts:
 later=((1<<n)-1)^((1<<(offset+len(rows)))-1);a.extend((r<<offset)|later for r in rows);cut|=m<<offset;offset+=len(rows)
star=[r|(1<<n)|((not(cut>>i&1))<<(n+1)) for i,r in enumerate(a)]+[1<<(n+1),cut]
for rootv in range(len(star)):
 reached=1<<rootv
 while True:
  nex=reached
  for i,row in enumerate(star):
   if reached>>i&1:nex|=row
  if nex==reached:break
  reached=nex
 require(reached==(1<<len(star))-1,'strongness')
N=1<<len(star);dp=[{} for _ in range(N)]
for i in range(len(star)):dp[1<<i][i]=1
for m in range(1,N):
 for i,count in dp[m].items():
  avail=star[i]&~m
  while avail:
   bit=avail&-avail;avail-=bit;j=bit.bit_length()-1;row=dp[m|bit];row[j]=row.get(j,0)+count
require(sum(dp[-1].values())==9909,'base actual count')
result=dict(status='PASS',core=core,length=12,root=root,forward_trace=trace,gain=gain,X_state=sx,Y_state=sy,flanked_left=left,flanked_right=right,global_delta=global_delta,checks=checks,base_star_adjacency=star,base_direct_H=9909)
(ROOT/'SHORT12_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
