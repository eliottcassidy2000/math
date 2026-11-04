"""Staircase prime quotient: exact physical sixteen-query controls."""
from pathlib import Path
import sys,json,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'marked_pair_quotient'))
from verify import strong,prime,battery,need

def quotient(r):
 n=2*r+1;T=[0]*n
 for u in range(n):
  for v in range(u+1,n):
   edge=(u>=r) or (u<r<=v and v-r<=u)
   T[u if edge else v]|=1<<(v if edge else u)
 return T

def parts(r,z):
 keep=[u for u in range(2*r+1) if u!=z]
 if z<r:return [r+z],[r+z+1]
 if z==r:return [u for u in keep if u!=0],[0]
 if z==2*r:return [r-1],[u for u in keep if u!=r-1]
 j=z-r;return [j],[j-1]

def main():
 structures=[];deletions=0
 for r in range(1,8):
  Q=quotient(r);n=len(Q)
  cyclic=[sum(1<<((u+d)%n) for d in range(1,r+1)) for u in range(n)]
  for u in range(r):
   for v in range(u+1,r):cyclic[u]^=1<<v;cyclic[v]^=1<<u
  need(Q==cyclic,'cyclic first-block reversal convention');need(strong(Q) and prime(Q),'strong prime staircase')
  for z in range(n):
   I,J=parts(r,z);M=sum(1<<u for u in I+J);need(I and J and len(I+J)==len(set(I+J)) and z not in I+J,'disjoint nonempty exterior sides')
   need(all(Q[u]>>v&1 for u in I for v in J),'ordered side direction')
   need(all(Q[u]&M in(0,M) for u in range(n) if u!=z and not M>>u&1),'deleted module homogeneity');deletions+=1
  structures.append({'r':r,'rows':Q})
 positive=[];boundary=[]
 for r in (1,2):
  Q=quotient(r)
  for z in range(len(Q)):
   I,J=parts(r,z);rec=battery(Q,[[2,0] for _ in Q],z,I,J);need(rec['quartic_coefficient']>0,'every root positive with nontrivial children');positive.append(rec)
 for r in (3,4):
  Q=quotient(r)
  for z in (0,r,r+1,2*r):
   blocks=[[0] for _ in Q];blocks[z]=[2,0];I,J=parts(r,z)
   if z not in(r,2*r):blocks[I[0]]=[2,0]
   rec=battery(Q,blocks,z,I,J);need(rec['quartic_coefficient']>0,'minimal sufficient child sizes');positive.append(rec)
 # Nontransitive interior and reversed pair away from the designated pair child.
 for z in (0,3):
  r=2;Q=quotient(r);I,J=parts(r,z);blocks=[[0] for _ in Q];blocks[z]=[2,0];blocks[I[0]]=[2,4,1];blocks[J[0]]=[0,1]
  rec=battery(Q,blocks,z,I,J);need(rec['quartic_coefficient']>0,'arbitrary child interiors');positive.append(rec)
 for r in (1,2):
  Q=quotient(r)
  for z in range(len(Q)):
   if r>1 and z in(r,2*r):continue
   blocks=[[0] for _ in Q];blocks[z]=[2,0];I,J=parts(r,z);rec=battery(Q,blocks,z,I,J);need(rec['quartic_coefficient']==0,'two-singleton exposed module boundary');boundary.append(rec)
 root=Path(__file__).resolve().parent;paths=[root/'PROOF.md',root/'verify.py',root.parent/'marked_pair_quotient'/'PROOF.md',root.parent/'marked_pair_quotient'/'verify.py',root.parent/'ordered_module_skew'/'PROOF.md']
 print(json.dumps({'status':'PASS','quotient_structure':structures,'deletion_module_checks':deletions,'positive_literal_batteries':positive,'zero_quartic_boundary_batteries':boundary,'source_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2,sort_keys=True))
if __name__=='__main__':main()
