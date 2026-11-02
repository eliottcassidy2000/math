"""Homogeneous-pair strict cyclic-vs-middle comparison; finite controls."""
import sys,itertools,json,hashlib
from reference_paths import tournament,extend,H,endpoint_Q_dp,c3_blowup

def ck(z,msg):
 if not z:raise RuntimeError(msg)
def induce(a,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
def field(a,C):
 n=len(a);A=[sum(C[i][j] for j in range(n) if a[i]>>j&1) for i in range(n)]
 return A
def ext_pair(a,R):
 # Make u,w labelled 0,1 with u->w and identical externally.
 n=len(a);b=[2,0]+[v<<2 for v in a]
 for i in range(n):
  if R>>i&1:b[0]|=1<<(i+2);b[1]|=1<<(i+2)
  else:b[i+2]|=3
 return tuple(b)
rows=0;least=None
for n in range(5):
 for code in range(1<<(n*(n-1)//2)):
  D=tournament(n,code);low=2*H(D)
  for R in range(1<<n):
   B=ext_pair(D,R)
   for S in range(1<<n):
    cyc=H(extend(B,(S<<2)|1));mid=H(extend(B,(S<<2)|2));gap=cyc-mid
    ck(gap>=low,'strict twin comparison')
    rows+=1
    if least is None or gap-low<least[0]:least=(gap-low,n,code,R,S,gap,low)
T=c3_blowup(tournament(4,4));B=induce(T,list(range(1,12)));s,e,K=endpoint_Q_dp(B);C=[[K[i][j]+K[j][i] for j in range(11)] for i in range(11)];D=[r[:] for r in C]
for j in [2,3,4]:D[0][j]-=498;D[j][0]-=498
for j in [5,6,7]:D[0][j]+=498;D[j][0]+=498
A=field(B,C);AA=field(B,D);diff=[C[0][j]-C[1][j] for j in range(2,11)];dd=[D[0][j]-D[1][j] for j in range(2,11)]
S=sum(1<<j for j in range(2,11) if dd[j-2]>0)
def f(S,C,A):return H(B)+sum(A[i] for i in range(11) if S>>i&1)-sum(C[i][j] for i in range(11) for j in range(i) if S>>i&1 and S>>j&1)
row={'literal_rows':rows,'equality_control':least,'formal_bad_external_cut':S,'formal_cyclic_cut':S|1,'formal_middle_cut':S|2,'actual_cyclic':f(S|1,C,A),'actual_middle':f(S|2,C,A),'formal_cyclic':f(S|1,D,AA),'formal_middle':f(S|2,D,AA),'actual_min_twin_gap':A[0]-A[1]-sum(max(x,0) for x in diff),'formal_min_twin_gap':AA[0]-AA[1]-sum(max(x,0) for x in dd),'lower_bound':2*H(induce(B,list(range(2,11))))}
ck(row['formal_min_twin_gap']==-90,'formal negative gap')
ck(row['actual_min_twin_gap']==432,'actual positive gap')
print(json.dumps(row,sort_keys=True))
print('PASS sha256='+hashlib.sha256(json.dumps(row,sort_keys=True).encode()).hexdigest())
