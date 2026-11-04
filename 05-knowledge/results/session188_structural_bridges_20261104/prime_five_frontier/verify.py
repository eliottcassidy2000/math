from pathlib import Path
from itertools import combinations
import hashlib,json
HERE=Path(__file__).resolve().parent
DEPS=['prime_five_all_ambient/PROOF.md','prime_five_all_profiles/PROOF.md','prime_five_even_ambient/PROOF.md','cyclic_ambient_synthesis/PROOF.md','module_decomposition/PROOF.md','highest_quotient_reduction/PROOF.md','monotone_dual_queries/PROOF.md']
def need(c,m):
 if not c:raise RuntimeError(m)
def op(T):
 f=(1<<len(T))-1
 return [f^(1<<i)^r for i,r in enumerate(T)]
def join(A,B):
 n=len(A);m=len(B)
 return [r|(((1<<m)-1)<<n) for r in A]+[r<<n for r in B]
def H(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for M in range(1,1<<n):
  for v in range(n):
   if M>>v&1:
    S=M^(1<<v);d[M][v]=1 if S==0 else sum(d[S][u] for u in range(n) if S>>u&1 and T[u]>>v&1)
 return sum(d[-1])
def aug(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for v in range(n):
  for mask,w in [(a,n),(c,n+1)]:
   if mask>>v&1:R[v]|=1<<w
   else:R[w]|=1<<v
 return R
def mods(T):
 n=len(T);out=[]
 for S in range(1,1<<n):
  if all((T[v]&S) in (0,S) for v in range(n) if not S>>v&1):out.append(S)
 return out
def prime(T):return not any(1<S.bit_count()<len(T) for S in mods(T))
def tournament(n,c):
 R=[0]*n
 for b,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if c>>b&1 else(j,i);R[u]|=1<<v
 return R
def main():
 p4=[c for c in range(64) if prime(tournament(4,c))];need(not p4,'no prime quotient order4')
 reps=[[20,1,2,7,14],[18,8,3,5,14],[24,17,3,6,12]];C=[2,4,1];controls=[]
 masks=[x for x in range(32) if x.bit_count()<=2]
 for Q in reps:
  need(prime(Q),'prime5 source');T=join(Q,C);V=join(op(Q),C);U=join(op(Q),op(C));need(31 in mods(T) and 224 in mods(T) and 31 in mods(U) and 224 in mods(U),'disjoint modules on both sides');need(H(T)==H(V)==H(U),'scalar equality for converse children')
  witness=None
  for a in masks:
   if witness is not None:break
   for c in masks:
    vals=[H(aug(W,a,c)) for W in (T,V,U)]
    if vals[0]!=vals[1]:
     need(vals[1]==vals[2],'exact other-node restoration under physical uniform ports');witness={'a':a,'c':c,'values':vals};break
  need(witness is not None,'prime5 ambient physical witness');controls.append({'quotient':Q,'source':T,'intermediate':V,'target':U,'original_family_masks':[31,224],'hamilton':H(T),'query':witness})
 pins={p:hashlib.sha256((HERE.parent/p).read_bytes()).hexdigest() for p in DEPS}
 pins['PROOF.md']=hashlib.sha256((HERE/'PROOF.md').read_bytes()).hexdigest();pins['verify.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 print(json.dumps({'status':'PASS','dependency_pins':pins,'prime_four_codes':p4,'multiple_highest_node_controls':controls,'scope':'Original-module XOR only; residual highest prime quotient order at least six.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
