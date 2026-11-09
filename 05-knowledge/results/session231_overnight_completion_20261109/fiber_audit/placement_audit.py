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

w=json.loads((ROOT/'SHORT12_WITNESS.json').read_text());core=literal_core(w['adj'],w['mask']);one=[(1,2,b,b,0,b,b) for b in (0,1)];X,Y=w['X'],w['Y'];contexts=values=0
for a in range(5):
 for P in (range(1<<(a-1),1<<a) if a else [0]):
  prefix=[] if not a else [word(format(P,f'0{a}b'))]
  for b in range(5):
   for R in range(1<<b):
    suffix=[] if not b else [word(format(R,f'0{b}b'))];rev=0 if not b else int(format(R,f'0{b}b')[::-1],2);baseline=[]
    for k in range(7):
     counts=[]
     for j in range(k+1):
      parts=prefix+[one[1],one[0]]
      if j:parts.append(word(X*j))
      parts.append(core)
      if k-j:parts.append(word(Y*(k-j)))
      parts+=[one[1]]+suffix;s=fold(parts);N=2*s[2]+s[0]-s[4];counts.append(N);values+=1
      require(N>=5760*(4096**k)*(1<<(a+b)),'path-cover bound')
     require(counts==[counts[0]-2*j for j in range(k+1)],'context interval');require((counts[0]-((1<<(b+1))-rev-13*P+7+5*k))%45==0,'mod45');baseline.append(counts[0])
    require(all(baseline[k+2]==4097*baseline[k+1]-4096*baseline[k]-264150 for k in range(5)),'fixed-context recurrence');contexts+=1
# Independent control of h95 five-coordinate period and scalar residues.
r=json.loads((ROOT/'WEIGHTED.json').read_text())['cycle'];c=tuple(r['core']);bits=''.join(str(e[4]) for e in r['edges']);yy=word(bits);require(yy[5]==0,'h95 Y start');z=fold([one[1],one[0],c]);initial=tuple(v%190 for v in z[:5]);residues=[];period=[]
for k in range(7):
 residues.append((4*z[2]+z[0]-z[4])%190);period.append(tuple(v%190 for v in z[:5]));z=multiply(z,yy)
require(residues==[149,169,119,149,169,119,149],'h95 residues');require(period[3]==initial and period[6]==initial,'h95 five-coordinate period');require(all(math.gcd(r,190)==1 for r in residues),'h95 gcd')
result=dict(status='PASS',contexts=contexts,scalar_values=values,k_max=6,h95_residues=residues,h95_five_coordinate_period3=True)
(ROOT/'PLACEMENT_AUDIT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
