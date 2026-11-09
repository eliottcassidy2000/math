import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit import literal,validate,check,canonical,decode
from pathlib import Path
import json,itertools as it
OUT=Path(__file__).parent

def ordinal(parts):
 a=[];mask=0;total=sum(len(g) for g,m in parts);offset=0
 for g,m in parts:
  for row in g:a.append((row<<offset)|(((1<<total)-1)^((1<<(offset+len(g)))-1)))
  mask|=m<<offset;offset+=len(g)
 validate(a);return a,mask

def close(a,m):
 n=len(a);b=[v|(1<<n)|((not(m>>i&1))<<(n+1)) for i,v in enumerate(a)]
 b += [1<<(n+1),m];validate(b);return b

def ham(a):
 # Independently append vertices to paths, rather than deletion-endpoint recurrence.
 n=len(a);full=(1<<n)-1;dp=[{} for _ in range(full+1)]
 for v in range(n):dp[1<<v][v]=1
 for m in range(1,full+1):
  for v,c in dp[m].items():
   left=a[v]&~m
   while left:
    bit=left&-left;left-=bit;u=bit.bit_length()-1;d=dp[m|bit];d[u]=d.get(u,0)+c
 return sum(dp[full].values())
def strong(a):
 for root in range(len(a)):
  reached=1<<root
  while True:
   nex=reached
   for v in range(len(a)):
    if reached>>v&1:nex|=a[v]
   if nex==reached:break
   reached=nex
  if reached!=(1<<len(a))-1:return False
 return True

if __name__ == "__main__":
 one=([0],1);zero=([0],0)
 groups={
  'H29': [([28,41,58,48,34,1],19),([20,41,58,49,34,1],19)],
  'H45':[([48,49,3,7,12,28],9),([48,49,3,7,12,28],3),([40,37,17,6,11,28],3),([48,49,3,7,12,28],6),([40,37,17,6,11,28],10)],
  'H5':[(decode(2050,7),50),(decode(1114112,7),70)],
  'H9':[(decode(32772,7),49),(decode(1056768,7),42),(decode(262400,7),54)]}
 out={}
 for name,cores in groups.items():
  rows=[]
  for core in cores:
   a,m=ordinal([one,core,zero]);state=literal(a)[m];b=close(a,m);H=ham(b)
   check(strong(b),'not strong');check(H==2*state[2]+state[0]-state[4],'star formula fails')
   rows.append(dict(core=core,state=state,star_adjacency=b,star_count=H,order=len(b)))
  check(len({tuple(r['state'][:4]+r['state'][5:]) for r in rows})==1,'fiber mismatch')
  out[name]=rows;print(name,[r['star_count'] for r in rows],flush=True)
 check(len({canonical(a)[0] for a in ([28,41,58,48,34,1],[20,41,58,49,34,1],[26,44,57,16,34,9])})==3,'H29 isomorphism')
 # Verify the published 29-vertex H95 graphs explicitly exist. Their count
 # evaluation uses the composition law after independent literal small-graph controls.
 def compose(x,y):
  h,q,A,G,T,S,E=x;j,k,B,J,U,V,W=y
  return h*j,q*k,A*k+h*B,J*q+j*G,A*J+T*j+h*U+E*V,S*j,h*W
 ss=literal([64,121,35,37,13,81,28]);ones=[(1,2,b,b,0,b,b) for b in (0,1)]
 # Test both orderings for representative extreme and ambiguous cuts.
 for m in (0,14,22,46,54,127):
  for b in (0,1):
   for parts,want in [([([64,121,35,37,13,81,28],m),([0],b)],compose(ss[m],ones[b])),([([0],b),([64,121,35,37,13,81,28],m)],compose(ones[b],ss[m]))]:
    a,k=ordinal(parts);check(literal(a)[k]==want,'literal composition control fails')
 path=[(1,1),(0,1),(0,1),(1,0),(0,0),(0,0),(0,1),(0,1),(0,0),(1,1),(0,0),(1,1),(0,1),(0,1),(1,1),(1,1),(1,0),(0,1),(0,1),(0,0)]
 rows=[]
 for side,m in [(0,14),(0,22),(1,46),(1,54)]:
  parts=[one]+([([64,121,35,37,13,81,28],m)]+[([0],y) for x,y in path] if side==0 else [([0],x) for x,y in path]+[([64,121,35,37,13,81,28],m)])+[one]
  a,k=ordinal(parts);b=close(a,k);check(strong(b),'H95 nonstrong');state=None
  for core,mask in parts:
   term=ss[mask] if len(core)>1 else ones[mask]
   state=term if state is None else compose(state,term)
  rows.append(dict(state=state,star_count=2*state[2]+state[0]-state[4],star_adjacency=b,order=len(b)))
 check({r['state'][:4]+r['state'][5:] for r in rows}=={(95,5351931904,4184607987,3843289869,95,95)},'H95 interface')
 check(sorted(r['star_count'] for r in rows)==[5364196511,5364196519,5364196551,5364196559],'H95 counts')
 out['H95']=rows
 (OUT/'PHYSICAL.json').write_text(json.dumps(out,indent=2))
