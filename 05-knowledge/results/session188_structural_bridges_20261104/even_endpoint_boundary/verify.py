"""Exact parity hostile to an even-order endpoint-only shortcut."""
from pathlib import Path
import hashlib,json
R=(598,508,344,865,40,69,656,317,113,438)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows_from_start(start):
 n=len(R);dp=[0]*(1<<n);dp[1<<start]=1<<start;all_end=0
 for used in range(1,1<<n):
  all_end^=dp[used];lasts=dp[used]
  while lasts:
   bit=lasts&-lasts;lasts-=bit;u=bit.bit_length()-1;remaining=R[u]&~used
   while remaining:
    b=remaining&-remaining;remaining-=b;dp[used|b]^=b
 return all_end,dp[-1]
def main():
 n=len(R);need(all(not (R[i]>>i&1) for i in range(n)),'no loops')
 need(all(((R[i]>>j)^(R[j]>>i))&1 for i in range(n) for j in range(i+1,n)),'tournament')
 rows=[rows_from_start(i) for i in range(n)]
 matrices=[[r[k] for r in rows] for k in range(2)];vectors=[]
 for P in matrices:
  starts=sum((p.bit_count()%2)<<u for u,p in enumerate(P));ends=0
  for p in P:ends^=p
  need(starts==ends==296,'equal start/end parity');vectors.append({'start':starts,'end':ends})
 P=matrices[0];witness=next((i,j) for i in range(n) for j in range(i+1,n) if ((P[i]>>j)^(P[j]>>i))&1)
 need(sum(p.bit_count() for p in P)%2==1 and sum(P[i]>>i&1 for i in range(n))%2==0,'odd total versus even diagonal')
 modules=[M for M in range(1,1<<n) if M.bit_count()>1 and all(r&M in(0,M) for u,r in enumerate(R) if not M>>u&1)]
 need(modules==[(1<<n)-1],'prime finite quotient')
 out={'status':'PASS','rows':R,'all_simple_path_parity_rows':matrices[0],'hamilton_path_parity_rows':matrices[1],'all_path_endpoint_vectors':vectors[0],'hamilton_endpoint_vectors':vectors[1],'asymmetric_pair':witness,'nontrivial_modules':modules,'scope':'Counterexample to even-order endpoint-parity nonvanishing and to zero symmetric row sums implying a symmetric all-path matrix. No complete-response collision is asserted.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 print(json.dumps(out,sort_keys=True,indent=2))
if __name__=='__main__':main()
