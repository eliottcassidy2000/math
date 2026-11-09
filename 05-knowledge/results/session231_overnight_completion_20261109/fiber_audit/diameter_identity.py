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

checks=0
for fn in ('UNIT_WITNESS.json','SHORT12_WITNESS.json'):
 w=json.loads((ROOT/fn).read_text());h,q,A,G,*_=w['core'];a,g,p=A,G,0;gain=0
 for i,(cx,cy) in enumerate(zip(w['X'],w['Y'])):
  x,y=int(cx),int(cy);aa=2*a+h*y-q*x;gg=(g+q*y-h*x)//2;weight=a*y+h*p*y-x*gg-h*x*(gg%2);gain+=weight;a,g,p=aa,gg,y
  sx=word(w['X'][:i+1]);sy=word(w['Y'][:i+1]);identity=A*sy[3]-sx[2]*g+h*(sy[4]-sx[4])-h*sx[6]*(g%2)
  require(identity==gain,'open-path telescoping');require(abs(weight)<=q+2*h-2,'edge bound');checks+=1
result=dict(status='PASS',open_path_prefix_identities=checks,edge_bound='q+2h-2')
(ROOT/'DIAMETER_IDENTITY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
