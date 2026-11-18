from math import isqrt
from selector import dp, connected

D=(42,116,121,114,97,64,1)
R=tuple(sum(1<<j for j in range(7) if (j-i)%7 in (1,2,3))|128 for i in range(7))+(0,)

def read(rows,I):
    n=len(rows);f=(1<<n)-1;A=dp(rows)
    B=dp([sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)])
    G=Q=0
    for s in range(1,f+1):
        w=sum(A[s][i] for i in I)
        t=f^s
        G+=w*(sum(B[t]) if t else 1)
        Q+=w*(sum(B[t][i] for i in range(n) if i not in I) if t else 1)
    return sum(A[f]),G,Q

def wrap(D,I,a):
    n=len(D);b=n+a
    return [r|((i in I)<<b) for i,r in enumerate(D)]+[
        sum(1<<j for j in range(i+1,b))|((1<<n)-1) for i in range(n,b)]+[
        sum(1<<i for i in range(b) if i not in I)]
for rows,I,expected in [(D,{0,5},(43,284,266)),(R,{7},(175,1988,1988))]:
    assert read(rows,I)==expected
    print('core h,G,Q:',expected)
for a in range(1,7):
    L,U=wrap(D,{0,5},a),wrap(R,{7},a)
    low,high=sum(dp(L)[-1]),sum(dp(U)[-1])
    assert connected(L) and connected(U)
    assert high==7*low==7*(284*2**a+25)
    isp=all(low%d for d in range(2,isqrt(low)+1))
    print(a,low,high,'prime lower target' if isp else '')
