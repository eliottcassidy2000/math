"""Exact alternating-hub quotient and marked-pair controls."""
from pathlib import Path
import hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'marked_pair_quotient'))
from verify import strong,prime,induced,battery,need


def quotient(r):
 n=2*r+1;T=[0]*n
 for u in range(2*r):
  for v in range(u+1,2*r):T[u]|=1<<v
  if u%2:T[u]|=1<<(2*r)
  else:T[2*r]|=1<<u
 return T

def parts(r,z):
 h=2*r;remaining=[u for u in range(h+1) if u!=z]
 if z==h:return [0],remaining[1:]
 if z==0:return [1],[u for u in remaining if u!=1]
 if z==h-1:return [u for u in remaining if u!=h-2],[h-2]
 return [z-1],[z+1]

def main():
 structure=[];deletions=0
 for r in range(1,8):
  Q=quotient(r);need(strong(Q) and prime(Q),'strong prime alternating quotient')
  for z in range(len(Q)):
   I,J=parts(r,z);keep=[u for u in range(len(Q)) if u!=z];M=set(I+J)
   need(I and J and not set(I)&set(J),'nonempty disjoint sides')
   need(all(Q[i]>>j&1 for i in I for j in J),'ordered module direction')
   need(all(all(bool(Q[v]>>w&1)==bool(Q[v]>>next(iter(M))&1) for w in M) for v in keep if v not in M),'deleted module homogeneity')
   deletions+=1
  structure.append({'r':r,'rows':Q})
 controls=[]
 for r in (1,2):
  Q=quotient(r)
  for z in range(len(Q)):
   I,J=parts(r,z);w=battery(Q,[[2,0]]*len(Q),z,I,J);need(w['quartic_coefficient']>0,'all-nontrivial child separator');controls.append(w)
 for r in (3,4):
  Q=quotient(r)
  for z in (0,2*r,2):
   blocks=[[0] for _ in Q];blocks[z]=[2,0]
   if z not in (0,2*r,2*r-1):blocks[z-1]=[2,0]
   I,J=parts(r,z);w=battery(Q,blocks,z,I,J);need(w['quartic_coefficient']>0,'large quotient minimal-size condition');controls.append(w)
 boundary=[]
 for z in range(3):
  Q=quotient(1);blocks=[[0] for _ in Q];blocks[z]=[2,0];I,J=parts(1,z);w=battery(Q,blocks,z,I,J);need(w['quartic_coefficient']==0,'two-singleton exterior zero boundary');boundary.append(w)
 root=Path(__file__).resolve().parent
 print(json.dumps({'status':'PASS','quotient_structure':structure,'deletion_module_checks':deletions,'positive_literal_batteries':controls,'zero_quartic_boundary_batteries':boundary,'source_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'PROOF.md',root/'verify.py',root.parent/'marked_pair_quotient'/'verify.py',root.parent/'ordered_module_skew'/'PROOF.md')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
