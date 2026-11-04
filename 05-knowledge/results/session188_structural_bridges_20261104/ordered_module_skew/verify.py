"""Independent literal cover injection and prime-five grouped quartic controls."""
from itertools import combinations,permutations,product
from pathlib import Path
import hashlib,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'marked_pair_quotient'))
from verify import rows,opposite,hp,covers,q2,induced,blow,strong,prime,battery,need


def module(T,M):
 mask=sum(1<<u for u in M)
 return all(T[u]&mask in (0,mask) for u in range(len(T)) if u not in M)
def paths_cover(T):
 n=len(T);out=[]
 for word in permutations(range(n)):
  for k in range(1,n):
   P,Q=word[:k],word[k:]
   if all(T[u]>>v&1 for Z in (P,Q) for u,v in zip(Z,Z[1:])):out.append((P,Q))
 return out

def injection(T,A,B):
 M=A|B;allcovers=paths_cover(T)
 left={(P,Q) for P,Q in allcovers if P[-1] in A and Q[0] in B}
 right={(P,Q) for P,Q in allcovers if P[-1] in B and Q[0] in A}
 image=set();reentries=0
 for P,Q in left:
  i=len(P)-1
  while i>0 and P[i-1] in M:i-=1
  j=1
  while j<len(Q) and Q[j] in M:j+=1
  need(set(P[i:])<=A and set(Q[:j])<=B,'pure marked runs')
  R,S=P[:i]+Q[:j],P[i:]+Q[j:]
  need((R,S) in right,'invalid swapped cover');need((R,S) not in image,'not injective');image.add((R,S))
  if any(x in M for x in P[:i]) or any(x in M for x in Q[j:]):reentries+=1
  ii=len(R)-1
  while ii>0 and R[ii-1] in M:ii-=1
  jj=1
  while jj<len(S) and S[jj] in M:jj+=1
  need(R[:ii]+S[:jj]==P and R[ii:]+S[jj:]==Q,'inverse swap')
 # Contract M to its first name; its exterior incidences are uniform.
 keep=[min(M)]+[u for u in range(len(T)) if u not in M];K=induced(T,keep)
 AA=induced(T,sorted(A));BB=induced(T,sorted(B));lower=hp(K)*(q2(AA)*hp(BB)+hp(AA)*q2(BB)+q2(AA)*q2(BB))
 need(len(right)-len(left)>=lower,'quantified positive skew')
 need((len(right)>len(left))==(len(M)>=3),'sharp singleton sides')
 return len(left),reentries


def ordered_parts(T):
 n=len(T);N=(1<<n)-1
 for am in range(1,N):
  A={u for u in range(n) if am>>u&1};B=set(range(n))-A
  if all(T[u]>>v&1 for u in A for v in B):return sorted(A),sorted(B)
 for u,v in combinations(range(n),2):
  if module(T,{u,v}):return ([u],[v]) if T[u]>>v&1 else ([v],[u])
 raise RuntimeError('four-vertex ordered module missing')


def main():
 cases=pairs=reentries=0
 corpus=[rows(n,c) for n in range(2,5) for c in range(1<<(n*(n-1)//2))]
 # M={0,1,2} is an ordered module; exterior cycle allows genuine reentry.
 T,_=blow([2,4,1],[[6,4,0],[2,0],[0]]);corpus.append(T)
 for T in corpus:
  n=len(T)
  for labels in product((0,1,2),repeat=n):
   A={u for u in range(n) if labels[u]==1};B={u for u in range(n) if labels[u]==2};M=A|B
   if not A or not B or not module(T,M) or not all(T[u]>>v&1 for u in A for v in B):continue
   a,b=injection(T,A,B);cases+=1;pairs+=a;reentries+=b
 need(reentries>0,'missing reentry control')
 selections=0;representatives={}
 for code in range(1024):
  Q=rows(5,code)
  for z in range(5):
   keep=[u for u in range(5) if u!=z];R=induced(Q,keep);I,J=ordered_parts(R)
   M=set(I)|set(J);need(module(R,M) and all(R[i]>>j&1 for i in I for j in J),'classified ordered module');selections+=1
  if prime(Q) and strong(Q):
   canon=min(tuple(induced(Q,p)) for p in permutations(range(5)))
   representatives.setdefault(canon,Q)
 controls=[]
 for Q in representatives.values():
  for z in range(5):
   keep=[u for u in range(5) if u!=z];I,J=ordered_parts(induced(Q,keep));I=[keep[i] for i in I];J=[keep[j] for j in J]
   w=battery(Q,[[2,0]]*5,z,I,J);need(w['quartic_coefficient']>0,'prime-five full query separator');controls.append(w)
 cyclic_structure=[];cyclic_controls=[]
 for r in range(2,8):
  n=2*r+1;Q=[sum(1<<((u+d)%n) for d in range(1,r+1)) for u in range(n)]
  need(strong(Q) and prime(Q),'cyclic strong prime quotient')
  for z in range(n):
   keep=[u for u in range(n) if u!=z];R=induced(Q,keep);index={u:i for i,u in enumerate(keep)}
   i,j=(z+r)%n,(z+r+1)%n
   need(module(R,{index[i],index[j]}) and Q[i]>>j&1,'cyclic deleted homogeneous pair')
  cyclic_structure.append({'order':n,'quotient':Q})
  if r in (3,4):
   blocks=[[0] for _ in range(n)];blocks[0]=[2,0];blocks[r]=[2,0]
   w=battery(Q,blocks,0,[r],[r+1]);need(w['quartic_coefficient']>0,'larger cyclic prime literal separator');cyclic_controls.append(w)
 root=Path(__file__).resolve().parent
 print(json.dumps({'status':'PASS','literal_ordered_module_cases':cases,'injected_ordered_covers':pairs,'covers_with_other_module_visits':reentries,'five_vertex_root_classifications':selections,'strong_prime_five_isomorphism_types':len(representatives),'literal_prime_five_batteries':controls,'cyclic_prime_structure':cyclic_structure,'larger_cyclic_prime_batteries':cyclic_controls,'source_sha256':{str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'PROOF.md',root/'verify.py',root.parent/'marked_pair_quotient'/'PROOF.md',root.parent/'marked_pair_quotient'/'verify.py')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
