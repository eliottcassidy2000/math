"""Independent finite calibrations of the color-cover transfer."""
from functools import cache


def dp(T):
    n=len(T);full=(1<<n)-1
    F=[[0]*n for _ in range(full+1)]
    for i in range(n):F[1<<i][i]=1
    for S in range(1,full+1):
        for i,p in enumerate(F[S]):
            if not p:continue
            A=T[i]&(full^S)
            while A:
                j=A&-A;A^=j
                F[S|j][j.bit_length()-1]+=p
    return F



@cache
def profile(T,k):
    H=[sum(a) if i else 1 for i,a in enumerate(dp(T))]
    F=[1]+[0]*(len(H)-1)
    for _ in range(k):
        for S in range(len(H)-1,-1,-1):
            V=S; z=0
            while True:
                z+=F[S^V]*H[V]
                if not V:break
                V=(V-1)&S
            F[S]=z
    return F[-1]
