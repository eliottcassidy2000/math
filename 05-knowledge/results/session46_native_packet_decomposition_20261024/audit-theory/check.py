"""Independent elementary group-algebra and exact-encoding proof controls."""
import argparse,hashlib,itertools,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--packet',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();checks=0

def need(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
def rank(vs):
 piv={}
 for v in vs:
  while v:
   k=v.bit_length()-1
   if k in piv:v^=piv[k]
   else:piv[k]=v;break
 return len(piv),list(piv.values())
reports=[]
for mods in ((2,),(4,),(2,4),(2,2,4)):
 G=list(itertools.product(*(range(m) for m in mods)));I={g:i for i,g in enumerate(G)};d=1+sum(m-1 for m in mods)
 def shift(v,g):
  out=0
  while v:
   bit=v&-v;i=bit.bit_length()-1;v^=bit
   out^=1<<I[tuple((x+y)%m for x,y,m in zip(G[i],g,mods))]
  return out
 basis=[1<<i for i in range(len(G))];dims=[len(basis)]
 for k in range(1,d+1):
  r,basis=rank([v^shift(v,g) for v in basis for g in G]);dims.append(r)
  expected=sum(sum(e)>=k for e in G);need(r==expected,'augmentation power dimension from actual multiplication')
 need(dims[-1]==0 and dims[-2]==1,'nilpotency degree exact')
 reports.append({'moduli':mods,'nilpotency_index':d,'ideal_dimensions':dims})
# Actual Q minimal-zero sequence: four order2 basis vectors,127 copies of order128 basis,then their sum.
zero=[]
for bits in itertools.product(range(2),repeat=4):
 for copies in range(128):
  for last in range(2):
   if all((b+last)%2==0 for b in bits) and (copies+last)%128==0:zero.append(sum(bits)+copies+last)
need(zero==[0,132],'explicit length132 sequence has only empty and full zero subsequences')
# Distinguish fixed-pair minimality from balancedcross-subset minimality.
P=(0,1);M=(1,0);diff=[(x-y)%3 for x,y in zip(P,M)]
need(diff==[2,1] and sum(diff)%3==0 and all(diff),'chosen pairing minimallyzero')
need(P[0]==M[1] and P[1]==M[0],'cross-pair balancedsubtrades neverthelessexist')
# A lattice combination can be impossible in a binary box.
states=[(1,0),(0,1)];moves=[(2,-2),(-2,2),(3,-3),(-3,3)]
need(tuple(moves[2][i]-moves[0][i] for i in range(2))==(1,-1),'integerlatticegeneration')
need(not any(tuple(s[i]+m[i] for i in range(2)) in states for s in states for m in moves),'no listed binaryfeasibleevent')
# Multi-targetmargin equations, trivial group: only synchronizedchoices.
for n in range(3,10):
 legal=[b for b in itertools.product(range(2),repeat=n) if all(1-b[j]+b[j-1]==1 for j in range(n))]
 need(legal==[(0,)*n,(1,)*n],'targetincidencealone forcesunboundedatom')
m=json.loads((a.packet/'MANIFEST.json').read_text())
for e in m['files']:need(hashlib.sha256((a.packet/e['path']).read_bytes()).hexdigest()==e['sha256'],'reviewedtheorymanifest '+e['path'])
result={'status':'PASS','checks':checks,'group_algebra_controls':reports,'Q_minimal_zero_lengths':zero,'manifest_sha256':hashlib.sha256((a.packet/'MANIFEST.json').read_bytes()).hexdigest(),'scope':'Independent finite proof controls plus source-hash verification. General theorems are reviewed by proofs; producer and inherited native replays are separately retained.'}
a.out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2))
