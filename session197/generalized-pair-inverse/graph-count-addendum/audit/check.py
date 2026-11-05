from fractions import Fraction as F
from itertools import product
import json
N=300
B=[0]*(N+1);B[1]=1
for n in range(3,N+1):B[n]=B[n-1]+sum((1<<(n-m-2))*B[m]for m in range(1,n-1))
W=[sum(B[:n+1])for n in range(N+1)]
N_counts=[0]*(N+1);M_counts=[0]*(N+1);P_counts=[0]*(N+1);Q_counts=[0]*(N+1)
for source in range(1,N+1):
 for length in range(N+1-source):
  choices=1 if length==0 else 1<<(length-1)
  if source>=2:
   if length and source+2+length<=N:N_counts[source+2+length]+=W[source]*choices
   if source+3+length<=N:
    M_counts[source+3+length]+=W[source]*choices
    Q_counts[source+3+length]+=W[source]*choices
  if source+4+length<=N:P_counts[source+4+length]+=W[source]*choices
presentations=[sum(x)for x in zip(N_counts,M_counts,P_counts,Q_counts)]
# Recover the claimed rational numerator by multiplication, not its recurrence.
def product_coeff(seq,pol):return[sum(c*seq[n-k]for k,c in enumerate(pol)if k<=n)for n in range(len(seq))]
if product_coeff(presentations,[1,-5,7,-2])!=[4 if n==5 else -8 if n==6 else 4 if n==7 else 0 for n in range(N+1)]:raise RuntimeError('exact presentation rational identity')
# Independent convolution of two Delta reciprocals, then the C2 numerator.
a=[1]
for n in range(1,N+1):a.append(3*a[-1]-(a[-2]if n>=2 else 0))
square=[sum(a[k]*a[n-k]for k in range(n+1))for n in range(N+1)]
C2=[]
for n in range(N+1):
 C2.append((n==5)+sum(c*square[n-k]for k,c in ((5,1),(6,-1),(7,-1),(8,1))if n>=k))
if(presentations[14],C2[14])!=(41736,44394):raise RuntimeError('order14 strict count gap')
if any(presentations[n]>=C2[n]for n in range(14,N+1)):raise RuntimeError('finite continued gap control')
# Direct tiny rooted-source counts are a guard against confusing roots with vertices.
if B[1:8]!=[1,0,1,3,8,21,55]or W[1:8]!=[1,1,2,5,13,34,89]:raise RuntimeError('rooted orbit and weak-tail bases')
result=dict(status='PASS',orders_checked=N+1,rooted_strong_first7=B[1:8],all_rooted_sources_first7=W[1:8],order14=dict(N=N_counts[14],M=M_counts[14],P=P_counts[14],Q=Q_counts[14],presentations=presentations[14],C2=C2[14],missing_at_least=C2[14]-presentations[14]),first_gap=next(n for n in range(N+1)if C2[n]>presentations[n]),finite_gap_controls=N-13,scope='Presentation upper bound only; no scalar omission.')
print(json.dumps(result,indent=2,sort_keys=True))
