"""Countercontrols for inference routes, not tournament counterexamples."""
from explore import tournament,c3_blowup,H,strong,extend,endpoint_Q_dp
from itertools import combinations
from functools import reduce
from math import gcd
import json,hashlib

def check(x,msg):
 if not x:raise RuntimeError(msg)
def induced(a,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
def field(a,C):
 n=len(a);A=[sum(C[i][j] for j in range(n) if a[i]>>j&1) for i in range(n)];f=[H(a)]*(1<<n)
 for S in range(1,1<<n):
  bit=S&-S;i=bit.bit_length()-1;P=S^bit;f[S]=f[P]+A[i]-sum(C[i][j] for j in range(n) if P>>j&1)
 return f,reduce(gcd,A+[C[i][j] for i,j in combinations(range(n),2)])
T=c3_blowup(tournament(4,4));B=induced(T,list(range(1,12)));s,e,Q=endpoint_Q_dp(B);C=[[Q[i][j]+Q[j][i] for j in range(11)] for i in range(11)];D=[row[:] for row in C]
for j in [2,3,4]:D[0][j]-=498;D[j][0]-=498
for j in [5,6,7]:D[0][j]+=498;D[j][0]+=498
f,g=field(B,C);ff,gg=field(B,D);old=1793;other=480;target=H(T)
check(H(B)==24651 and target==241137,'actual base and child values')
check(f[other]==239643 and ff[other]==target,'formal repair changes targeted level')
check(f[old]==ff[old]==target,'original count preserved')
check(g==gg==6,'exact ear lattice preserved')
check(sum(sum(row) for row in C)==sum(sum(row) for row in D),'total symmetric capacity preserved')
for m in range(12):check(sum(f[S] for S in range(2048) if S.bit_count()==m)==sum(ff[S] for S in range(2048) if S.bit_count()==m),'every cardinality-layer mean preserved')
maxbound=0
for i,j in combinations(range(11),2):
 lower=2*H(induced(B,[v for v in range(11) if v not in [i,j]]));maxbound=max(maxbound,lower);check(D[i][j]>=lower,'every actual pair-deletion lower bound preserved')
check(all(v%2==1 and v>=H(B) for v in ff),'parity and positivity')
realzeros=[S for S,v in enumerate(f) if v==target];formalzeros=[S for S,v in enumerate(ff) if v==target]
check(realzeros==[old] and len(formalzeros)>1,'real versus formal zero uniqueness')
# Strongness is needed in the prospective uniqueness conjecture.
NT=c3_blowup(tournament(2,1));NB=induced(NT,list(range(1,6)));nz=[S for S in range(32) if H(extend(NB,S))==H(NT)]
check(not strong(NT) and nz==[14,22,26,29],'nonstrong quotient countercontrol')
row={'test':'formal network zero-closure firewall','actual_parent_masks':B,'parent_H':H(B),'original_H':target,'original_cut':old,'alternative_cut':other,'actual_alternative_H':f[other],'modified_alternative_H':ff[other],'ear_gcd_both':g,'total_unordered_capacity':sum(C[i][j] for i,j in combinations(range(11),2)),'global_maximum_pair_lower_bound':maxbound,'modified_orbits':{'local0_to_local2_3_4':[7182,6684],'local0_to_local5_6_7':[5616,6114]},'actual_zero_cuts':realzeros,'formal_zero_cuts':formalzeros,'nonstrong_countercontrol':{'masks':NT,'H':H(NT),'parent_masks':NB,'zero_cuts':nz}}
print(json.dumps(row,sort_keys=True,separators=(',',':')))
print('PASS semantic_sha256='+hashlib.sha256(json.dumps(row,sort_keys=True,separators=(',',':')).encode()).hexdigest())
