from itertools import product
from math import isqrt, prod
from pathlib import Path
from runpy import run_path
bank=run_path(str(Path(__file__).with_name("third-slot-bank.py")))
data=[[d,q,u[0],v[0]] for d,q,u,v in bank["data"]]
cert={460648142166497107: [3, {2: 1, 3: 1, 56951: 1, 1348083270901: 1}], 56951: [13, {2: 1, 5: 2, 17: 1, 67: 1}], 1348083270901: [2, {2: 2, 3: 2, 5: 2, 19: 1, 78835279: 1}], 78835279: [3, {2: 1, 3: 1, 13139213: 1}], 13139213: [2, {2: 2, 743: 1, 4421: 1}], 461036408494562707: [2, {2: 1, 3: 1, 13: 1, 19: 1, 109: 1, 2854043064137: 1}], 2854043064137: [3, {2: 3, 11: 1, 32432307547: 1}], 32432307547: [3, {2: 1, 3: 1, 107: 1, 50517613: 1}], 50517613: [5, {2: 2, 3: 2, 163: 1, 8609: 1}], 461046213093554707: [2, {2: 1, 3: 1, 61: 1, 1259689106812991: 1}], 1259689106812991: [7, {2: 1, 5: 1, 7: 1, 79: 1, 89: 1, 467: 1, 733: 1, 7477: 1}], 461025215748156307: [2, {2: 1, 3: 1, 13: 1, 103: 1, 389: 1, 147517400581: 1}], 147517400581: [2, {2: 2, 3: 1, 5: 1, 11: 1, 41: 1, 5451493: 1}], 5451493: [2, {2: 2, 3: 1, 17: 1, 26723: 1}], 26723: [2, {2: 1, 31: 1, 431: 1}], 461413800957885907: [3, {2: 1, 3: 1, 13: 1, 2093323: 1, 2825919149: 1}], 2093323: [5, {2: 1, 3: 1, 7: 1, 11: 1, 23: 1, 197: 1}], 2825919149: [2, {2: 2, 13: 1, 54344599: 1}], 54344599: [6, {2: 1, 3: 1, 7: 1, 11: 1, 19: 1, 41: 1, 151: 1}], 461423613609357907: [3, {2: 1, 3: 1, 13: 1, 17: 1, 65777: 1, 5290323503: 1}], 65777: [3, {2: 4, 4111: 1}], 5290323503: [5, {2: 1, 109: 1, 24267539: 1}], 24267539: [2, {2: 1, 12133769: 1}], 12133769: [3, {2: 3, 73: 1, 79: 1, 263: 1}]}
def prime(n):
    if n<10000:return n>=2 and all(n%p for p in range(2,isqrt(n)+1))
    g,F=cert[n]
    return (prod(q**e for q,e in F.items())==n-1 and all(prime(q) for q in F)
        and pow(g,n-1,n)==1 and all(pow(g,(n-1)//q,n)!=1 for q in F))
B=[data[i] for i in [0,8,11]]
D=prod(x[0] for x in B); Q=prod(x[1] for x in B)
assert (D,Q,-3*D-6*Q)==(7683910325,36380429405056,-218305628161311)
for s in product(range(2), repeat=3):
    R=prod(x[2+j] for x,j in zip(B,s)); h=11*D+42*Q+20*R
    if s==(0,1,1):assert h%109==0
    elif s==(1,1,1):assert h%17==0
    else:assert prime(h)
def pf(n):
    return [p for p in range(2,n+1) if n%p==0 and prime(p)]
def rank(rows,p):
    A=[[v%p for v in row] for row in rows]
    i=0
    for j in range(len(A[0])):
        b=next((k for k in range(i,len(A)) if A[k][j]),None)
        if b is None:continue
        A[i],A[b]=A[b],A[i]
        for k in range(i+1,len(A)):
            x,y=A[i][j],A[k][j]
            A[k]=[(x*u-y*v)%p for u,v in zip(A[k],A[i])]
        i+=1
        if i==len(A):break
    return i
for P,I,G in [
  ([7,13,17,19,29,37,41,47],[3,4,6,8,9,10,12,13],[3,2,3,2,2,2,6,5]),
  ([13,23,31,37,41,47],[1,2,3,6,7,13],[2,5,3,2,6,5]),
]:
    L=[]
    for p,g in zip(P,G):
        f=p*(p-1);m=p*p
        log={pow(g,k,m): k for k in range(f)}
        assert len(log)==f
        L.append([log[data[i][2]%m]-log[data[i][3]%m] for i in I])
    for q in set(q for p in P for q in pf(p*(p-1))):
        rows=[v for p,v in zip(P,L) if p*(p-1)%q==0]
        assert rank(rows,q)==len(rows)
print('six certified primes; both prime-power wheels verified')
