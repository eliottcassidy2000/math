#!/usr/bin/env python3
"""Explicit checks, valid under Python -O; no assertions are relied upon."""
from collections import Counter
from fractions import Fraction
from itertools import product
from functools import reduce
from pathlib import Path
import hashlib,json,sys
sys.path.insert(0,str(Path(__file__).parent/'sources'))
import session55_verify as old

def check(ok,msg):
 if not ok: raise RuntimeError(msg)
def row(d):return tuple(sorted(d.items()))
def fuse(a,b):
 d=dict(a)
 for i,s in b:
  if i in d and d[i]!=s:return None
  d[i]=s
 return row(d)
def join(A,B):return {c for a in A for b in B if (c:=fuse(a,b)) is not None}
def project(A,scope):return {tuple((i,s) for i,s in a if i in scope) for a in A}
def semi(A,B):return {a for a in A if any(fuse(a,b) is not None for b in B)}
def relation(scope,flag):
 return {tuple(zip(scope,v)) for n,v in enumerate(product((0,1),repeat=len(scope))) if flag>>n&1}
def weight(a,W):return reduce(lambda x,is_:x*W[is_],a,Fraction(1))
def local(q):
 Q=old.roots(q)
 return sorted({old.magic(q,c,a,b) for c,a,b in product(Q,repeat=3) if all(x in Q for x in old.magic(q,c,a,b))})
def hist(A,states,W=None):
 H=Counter()
 for a in A:
  m=reduce(lambda x,is_:x&old.mask(states[is_[0]][is_[1]]),a,old.ALL)
  H[m]+=1 if W is None else weight(a,W)
 return H
def tree_count(R,W):
 # Three-node running-intersection chain ab--bc--cd.
 # Own a,b at root ab; c at bc; d at cd. Weight each variable once.
 msg2=Counter()
 for r in R[2]:
  d=dict(r);msg2[d[2]]+=W[(3,d[3])]
 msg1=Counter()
 for r in R[1]:
  d=dict(r);msg1[d[1]]+=W[(2,d[2])]*msg2[d[2]]
 return sum(W[(0,dict(r)[0])]*W[(1,dict(r)[1])]*msg1[dict(r)[1]] for r in R[0])
def main():
 checks={}
 partial=[row({i:s for i,s in enumerate(v) if s!=-1}) for v in product((-1,0,1),repeat=3)]
 for a,b,c in product(partial,repeat=3):
  ab=fuse(a,b);bc=fuse(b,c)
  check((None if ab is None else fuse(ab,c))==(None if bc is None else fuse(a,bc)),'partial assignment associativity')
 checks['partial_assignment_associativity_triples']=len(partial)**3
 W={(i,s):Fraction((i+1)*s,3) for i in range(4) for s in (0,1)}
 # Includes many zero weights; no division is used.
 for flags in product(range(16),repeat=3):
  R=[relation(scope,f) for scope,f in zip(((0,1),(1,2),(2,3)),flags)]
  J=join(join(R[0],R[1]),R[2])
  check(J==join(R[0],join(R[1],R[2])),'relation associativity')
  check(tree_count(R,W)==sum(weight(a,W) for a in J),'tree weighted count')
  # Bottom-up then top-down full semijoin reduction.
  S=list(R);S[1]=semi(S[1],S[2]);S[0]=semi(S[0],S[1]);S[1]=semi(S[1],S[0]);S[2]=semi(S[2],S[1])
  check(all(S[i]==project(J,set((i,i+1))) for i in range(3)),'tree global extension')
 checks['binary_chain_relation_triples']=16**3
 mods=(2,5,3);states={i:local(q) for i,q in enumerate(mods)}
 leaves=[{((i,s),) for s in range(len(states[i]))} for i in range(3)]
 left=join(leaves[0],leaves[1]);right=join(leaves[1],leaves[2]);good=join(left,right)
 true=join(join(leaves[0],leaves[1]),leaves[2])
 check(good==true,'overlap exact assignment set')
 check(join(left,left)==left,'relation self-copy')
 H=hist(good,states);direct=Counter(old.mask(a) for a in local(30))
 check(H==direct and sum(H.values())==112 and 0 not in H,'mod30 exact histogram')
 crt_arrays={tuple(old.crt(v,mods) for v in zip(*(states[i][s] for i,s in a))) for a in good}
 check(crt_arrays==set(local(30)),'CRT array bijection')
 check(all(old.valid(30,a) for a in crt_arrays),'CRT line equations')
 check(set(H)==set(old.independent_masks(30)),'independent line-constraint support')
 bad=Counter()
 for a,b in product(left,right):
  mask=reduce(lambda m,is_:m&old.mask(states[is_[0]][is_[1]]),a+b,old.ALL);bad[mask]+=1
 check(sum(bad.values())==784 and bad[0]>0,'hostile independent duplication')
 WF={(i,s):Fraction((s%3),s+2) for i in states for s in range(len(states[i]))}
 check(hist(good,states,WF)==hist(true,states,WF),'fractional once-only weights')
 # Every normalized nonzero integral vector on four masks with entries 0,1,2.
 count=0
 for v in product(range(3),repeat=4):
  if not sum(v):continue
  f={s:Fraction(w,sum(v)) for s,w in enumerate(v) if w}
  diagonal={(s,s):w for s,w in f.items()};push=Counter()
  for (s,t),w in diagonal.items():push[s&t]+=w
  check(push==f,'meet diagonal identity');count+=1
 checks['rational_diagonal_laws']=count
 # Evaluation is not a multiplicative map on choice projectors.
 p=Fraction(1,3);check(p!=p*p,'nonhomomorphism hostile')
 # Equal aggregate masks still need distinct original state IDs for exact counts.
 buckets={}
 for s,a in enumerate(states[1]):buckets.setdefault(old.mask(a),[]).append(s)
 collision=next((ss for ss in buckets.values() if len(ss)>1),None)
 check(collision is not None,'same-mask distinct choices exist')
 coarse_pairs=sum(len(ss)**2 for ss in buckets.values())
 check(coarse_pairs>len(states[1]),'matching only masks overcounts')
 # Complete marginal overlap interfaces do not imply lossless decomposition.
 parity={row(dict(enumerate(v))) for v in product((0,1),repeat=3) if sum(v)%2==0}
 marginal_join=join(project(parity,{0,1}),project(parity,{1,2}))
 check(len(parity)==4 and len(marginal_join)==8,'lossy marginal decomposition hostile')
 eq01=relation((0,1),9);eq12=relation((1,2),9);ne02=relation((0,2),6)
 cyc=[eq01,eq12,ne02]
 check(all(semi(A,B)==A for A,B in product(cyc,repeat=2)),'cycle pairwise consistent')
 check(not join(join(eq01,eq12),ne02),'cycle globally inconsistent')
 result={'schema':'correlated-choice-gluing-v1','checks':checks,'hostile':{'moduli':mods,'local_state_counts':[len(states[i]) for i in states],'left_rows':len(left),'right_rows':len(right),'repaired_rows':len(good),'independent_rows':sum(bad.values()),'repaired_distinct_count':H[0],'fabricated_distinct_count':bad[0],'histogram':{hex(k):v for k,v in sorted(H.items())},'same_mask_F5_pair_count':coarse_pairs,'F5_state_count':len(states[1])},'boundary_controls':{'parity_original_rows':len(parity),'marginal_join_rows':len(marginal_join),'cycle_pairwise_consistent':True,'cycle_global_rows':0,'zero_weight_safe':True,'scalar_evaluation_not_multiplicative':True}}
 p=Path(__file__).with_name('certificate.json')
 if '--verify' in sys.argv:check(json.loads(p.read_text())==json.loads(json.dumps(result)),'certificate changed')
 else:p.write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2));print('PASS')
if __name__=='__main__':main()
