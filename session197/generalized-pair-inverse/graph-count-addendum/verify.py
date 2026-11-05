from pathlib import Path
import json
N=100
F=[0,1]
for _ in range(2,2*N):F.append(F[-1]+F[-2])
B=[0]*(N+1);B[1]=1
for n in range(3,N+1):B[n]=F[2*n-4]
W=[sum(B[:n+1])for n in range(N+1)];S=[1]+[2**(m-1)for m in range(1,N+1)]
P=[0]*(N+1)
for n in range(1,N+1):
 for m in range(N+1-n):
  if n>=2 and n+2+m<=N and m>=1:P[n+2+m]+=W[n]*S[m]
  if n>=2 and n+3+m<=N:P[n+3+m]+=2*W[n]*S[m]
  if n+4+m<=N:P[n+4+m]+=W[n]*S[m]
def rational(num,den):
 out=[]
 for n in range(N+1):out.append(num.get(n,0)-sum(den[i]*out[n-i]for i in range(1,min(n,len(den)-1)+1)))
 return out
R=rational({5:4,6:-8,7:4},[1,-5,7,-2])
if P!=R:raise RuntimeError('presentation convolution vs rational function')
A=rational({5:1,6:-1,7:-1,8:1},[1,-6,11,-6,1]);A[5]+=1
if(P[14],A[14])!=(41736,44394):raise RuntimeError('order14 separation')
receipt=dict(status='PASS',coefficient_controls=N+1,first_coefficient_separation=next(n for n in range(N+1)if A[n]>P[n]),order14=dict(presentations=P[14],all_C2=A[14],minimum_missing=A[14]-P[14]),rows=[dict(order=n,presentations=P[n],all_C2=A[n])for n in range(5,25)])
text=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(Path(__file__).parent/'VALIDATION.json').write_text(text);print(text,end='')
