"""Independent original-tree and literal restoration audit of the prime-five integration."""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
from math import factorial
import hashlib,json,subprocess,sys
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'prime_five_frontier'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def opposite(T):return [r^((1<<len(T))-1)^(1<<i) for i,r in enumerate(T)]
def modules(T):
 n=len(T);out=[]
 for S in range(1,1<<n):
  vs=[i for i in range(n) if S>>i&1];a=vs[0]
  if all(all(((T[x]>>v)&1)==((T[x]>>a)&1) for v in vs) for x in range(n) if not S>>x&1):out.append(S)
 return out
def strong(T):
 ms=modules(T)
 return [A for A in ms if not any(A&B and A&~B and B&~A for B in ms)]
def tree(T,U):
 ss=strong(T);nodes=[];changed=[]
 for M in ss:
  if M.bit_count()<2:continue
  subs=[A for A in ss if A!=M and A&~M==0];children=[A for A in subs if not any(A!=B and A&~B==0 for B in subs)];children.sort(key=lambda x:(x&-x).bit_length());need(sum(children)==M,'strong children partition original names');reps=[(a&-a).bit_length()-1 for a in children];qs=[]
  for V in (T,U):
   need(all(A in modules(V) for A in children),'persistent child modules')
   qs.append([sum(1<<j for j,v in enumerate(reps) if V[u]>>v&1) for u in reps])
  if qs[0]!=qs[1]:changed.append(M)
  nodes.append({'mask':M,'children':children,'quotients':qs,'changed':qs[0]!=qs[1]})
 high=[A for A in changed if not any(A!=B and A&~B==0 for B in changed)];need(all(not(A&B) for A,B in combinations(high,2)),'highest nodes disjoint')
 for i,j in combinations(range(len(T)),2):
  if (T[i]^U[i])>>j&1:need(any(A>>i&1 and A>>j&1 for A in high),'every changed edge contained in a highest node')
 return {'all_original_strong_masks':ss,'nodes':nodes,'highest_changed':sorted(high)}
def count(T):
 n=len(T)
 @lru_cache(None)
 def finish(S,u):
  if not S:return 1
  return sum(finish(S^(1<<v),v) for v in range(n) if S>>v&1 and T[u]>>v&1)
 return sum(finish(((1<<n)-1)^(1<<u),u) for u in range(n))
def profiles(T):
 n=len(T)
 @lru_cache(None)
 def finish(S,u):
  if not S:return 1
  return sum(finish(S^(1<<v),v) for v in range(n) if S>>v&1 and T[u]>>v&1)
 H=[0]*(1<<n);H[0]=1
 for S in range(1,1<<n):H[S]=sum(finish(S^(1<<v),v) for v in range(n) if S>>v&1)
 @lru_cache(None)
 def partitions(S,k):
  if not S:return int(k==0)
  if k==0:return 0
  pivot=S&-S;free=S^pivot;B=free;z=0
  while True:
   A=B|pivot;z+=H[A]*partitions(S^A,k-1)
   if not B:break
   B=(B-1)&free
  return z
 return [[partitions(S,k)*factorial(k) for k in range(n+1)] for S in range(1<<n)]
def aug(T,a,c):
 n=len(T);R=T[:]+[1<<(n+1),0]
 for v in range(n):
  if a>>v&1:R[v]|=1<<n
  else:R[n]|=1<<v
  if c>>v&1:R[v]|=1<<(n+1)
  else:R[n+1]|=1<<v
 return R
def blow(Q,blocks):
 offs=[];N=0
 for B in blocks:offs.append(N);N+=len(B)
 R=[0]*N
 for i,B in enumerate(blocks):
  for u,row in enumerate(B):R[offs[i]+u]=(row<<offs[i])|sum(((1<<len(blocks[j]))-1)<<offs[j] for j in range(len(Q)) if Q[i]>>j&1)
 return R
def xor(T,masks):
 V=T[:]
 for A in masks:
  for i,j in combinations(range(len(T)),2):
   if A>>i&1 and A>>j&1:V[i]^=1<<j;V[j]^=1<<i
 return V
def main():
 src=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired source receipts')
 for f,h in src['dependency_pins'].items():need(digest((SRC if f in ['PROOF.md','verify.py'] else ROOT.parent)/f)==h,'source and dependency pin '+f)
 rec=[];hcounts=0;subsets=0
 for r in src['multiple_highest_node_controls']:
  T,V,U=r['source'],r['intermediate'],r['target'];need(xor(T,r['original_family_masks'])==U,'literal original-family XOR');tr=tree(T,U);need(tr['highest_changed']==[31,224],'genuine two highest changed nodes');P,Q=profiles(T),profiles(U);need(P==Q,'complete hereditary profile equality');subsets+=len(P);a,c=r['query']['a'],r['query']['c'];values=[count(aug(W,a,c)) for W in (T,V,U)];need(values==r['query']['values'] and values[0]!=values[1]==values[2],'literal source/baseline/target restoration');hcounts+=3;need(count(T)==count(V)==count(U)==r['hamilton'],'scalar equality');hcounts+=3;rec.append({'source_control':r,'original_tree':tr,'literal_values':values,'hereditary_profile_rows':len(P)})
 # Prime parent, non-singleton prime5 child, changed descendant pair and another highest node.
 Q=[18,8,3,5,14];C=[2,4,1];children=[[2,0],[0],[0],[0],[0]];M=blow(Q,children);B=blow(opposite(Q),children);T=blow(C,[M,C,[0]]);V=blow(C,[B,C,[0]]);W=blow(C,[opposite(M),C,[0]]);U=blow(C,[opposite(M),opposite(C),[0]]);need(xor(T,[63,448])==U,'nested original-family XOR');tr=tree(T,U);need(tr['highest_changed']==[63,448],'prime-parent highest-node antichain');need(next(r for r in tr['nodes'] if r['mask']==3)['changed'],'descendant pair changed but not highest');P,Qp=profiles(T),profiles(U);need(P==Qp,'nested complete hereditary equality');subsets+=len(P);a,c=515,512;vv=[count(aug(Z,a,c)) for Z in (T,V,W,U)];hcounts+=4;need(vv[0]!=vv[1]==vv[2]==vv[3],'nested child restoration and disjoint-node restoration');need((vv[0]-vv[1])%4==2,'physical endpoint residue2');nested={'source':T,'baseline':V,'child_restored':W,'target':U,'original_family_masks':[63,448],'original_tree':tr,'query':[a,c],'literal_values':vv,'hereditary_profile_rows':len(P)}
 p4=[]
 for code in range(64):
  T=[0]*4
  for bit,(i,j) in enumerate(combinations(range(4),2)):
   v,w=(i,j) if code>>bit&1 else(j,i);T[v]|=1<<w
  if not any(1<A.bit_count()<4 for A in modules(T)):p4.append(code)
 need(not p4,'independent no-prime-four census')
 replay=[]
 for flags in [[],['-O']]:
  out=subprocess.check_output([sys.executable,*flags,str(SRC/'verify.py')]);need(out==(SRC/'normal.json').read_bytes(),'fresh source replay');replay.append(hashlib.sha256(out).hexdigest())
 return {'status':'PASS','source_pins':{f:digest(SRC/f) for f in ['PROOF.md','README.md','verify.py','normal.json','optimized.json']},'dependency_pins':src['dependency_pins'],'independent_source_controls':rec,'nested_prime_parent_control':nested,'literal_H_counts':hcounts,'complete_hereditary_subset_rows':subsets,'prime_four_codes':p4,'producer_replay_hashes':replay,'scope':'Highest changed prime quotients have order at least6 in any possible nontrivial full-response collision under simultaneous XOR of original modules. This is source-relative, not arbitrary-output injectivity or current-module inversion sequences.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
