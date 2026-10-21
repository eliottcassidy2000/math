from math import comb
N=16
T=[2**(j*(j-1)//2) for j in range(N+1)]; s=[0]*(N+1)
for k in range(1,N+1): s[k]=T[k]-sum(comb(k,i)*s[i]*T[k-i] for i in range(1,k))
c=lambda m: 2**((m-1)//2) if m%2 else 0
import subprocess
for r in range(1,N+1):
    D=s[r]-c(r)
    f=subprocess.run(['gp','-q'],input=f'print(factor({D}))',capture_output=True,text=True).stdout.strip() if D else '0'
    print(r, s[r], D, f)
