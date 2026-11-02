"""Formal network still has a zero escape after preserving twin-lift constraints.

This is deliberately NOT an actual tournament kernel or counterexample.
"""
import sys,itertools,json,hashlib,math,functools
from reference_paths import tournament,c3_blowup,H,endpoint_Q_dp

def ck(z,msg):
 if not z:raise RuntimeError(msg)
def induce(a,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
def field(a,C):
 n=len(a);A=[sum(C[i][j] for j in range(n) if a[i]>>j&1) for i in range(n)];f=[H(a)]*(1<<n)
 for S in range(1,1<<n):
  bit=S&-S;i=bit.bit_length()-1;P=S^bit;f[S]=f[P]+A[i]-sum(C[i][j] for j in range(n) if P>>j&1)
 return A,f,functools.reduce(math.gcd,A+[C[i][j] for i,j in itertools.combinations(range(n),2)])
T=c3_blowup(tournament(4,4));B=induce(T,list(range(1,12)));s,e,K=endpoint_Q_dp(B);C=[[K[i][j]+K[j][i] for j in range(11)] for i in range(11)];D=[r[:] for r in C]
for i in [2,3,4]:
 for j in [8,9,10]:D[i][j]+=498;D[j][i]+=498
for i in [5,6,7]:
 for j in [8,9,10]:D[i][j]-=498;D[j][i]-=498
A,f,g=field(B,C);AA,ff,gg=field(B,D);old=1793;other=480;target=H(T)
ck(f[old]==ff[old]==target,'canonical target')
ck(f[other]==239643 and ff[other]==target,'formal zero escape')
ck(g==gg==6,'exact gcd')
ck(sum(map(sum,C))==sum(map(sum,D)),'total capacity')
for m in range(12):ck(sum(f[S] for S in range(2048) if S.bit_count()==m)==sum(ff[S] for S in range(2048) if S.bit_count()==m),'layer means')
for i,j in itertools.combinations(range(11),2):ck(D[i][j]>=2*H(induce(B,[v for v in range(11) if v not in [i,j]])),'pair-deletion lower bound')
for i in range(2,11):ck(D[0][i]-D[1][i]==C[0][i]-C[1][i],'twin gradient')
ck(D[0][1]==C[0][1],'twin pair capacity')
for a in [2,5,8]:
 for i,j in itertools.combinations(range(a,a+3),2):ck(C[i][j]==D[i][j]==2*H(induce(B,[v for v in range(11) if v!=i])),'internal triangle identity')
low=2*H(induce(B,list(range(2,11))));mingap=min(ff[S|1]-ff[S|2] for S in range(0,2048,4))
ck(mingap==432 and mingap>=low,'every cyclic-middle lower bound')
row={'parent_masks':B,'target':target,'canonical_cut':old,'other_cut':other,'actual_other':f[other],'formal_other':ff[other],'gcd_both':g,'total_capacity_both':sum(map(sum,C))//2,'modified_cross_orbits':{'blocks1_3':[13770,14268],'blocks2_3':[10002,9504]},'formal_zero_cuts':[S for S in range(2048) if ff[S]==target],'minimum_twin_gap_both':mingap,'twin_lower_bound':low,'preserved':['positivity','C3 symmetries','all pair-deletion lower bounds','exact gcd','total capacity','all cardinality-layer means','parent H','canonical child H','all twin capacity differences','twin pair capacity','all homogeneous-C3 internal capacities','every cyclic-middle lower bound']}
print(json.dumps(row,sort_keys=True))
print('PASS sha256='+hashlib.sha256(json.dumps(row,sort_keys=True).encode()).hexdigest())
