# exact s(n) (strong labeled tournaments) for n <= 80 and per-prime-power test of s(n) = 2^((n-1)/2)
from math import comb
N=200
T=[2**(j*(j-1)//2) for j in range(N+1)]
s=[0]*(N+1)
for k in range(1,N+1):
    s[k]=T[k]-sum(comb(k,i)*s[i]*T[k-i] for i in range(1,k))
print('s(1..10)=',s[1:11])
def factor(n):
    f={};d=2;x=n
    while d*d<=x:
        while x%d==0: f[d]=f.get(d,0)+1; x//=d
        d+=1
    if x>1: f[x]=f.get(x,0)+1
    return f
rows=[]
for n in range(9,N+1,2):
    f=factor(n)
    if len(f)==1 and list(f.values())[0]==1: continue
    circ=pow(2,(n-1)//2)
    res={}
    for p,a in f.items():
        q=p**a; res[f'{p}^{a}']=( (s[n]-circ)%q==0 )
    ok=all(res.values())
    rows.append((n,ok,res))
for n,ok,res in rows:
    print(n, 'PASS' if ok else 'fail', ' '.join(f"{k}:{'y' if v else 'n'}" for k,v in res.items()))
