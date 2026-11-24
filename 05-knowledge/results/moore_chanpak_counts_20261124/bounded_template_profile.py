"""Independent finite calibrations of the color-cover transfer."""
from functools import cache
from math import prod
from counts import tournament, dp, cycle5

one=(0,)
two=(2,0)
tri=(2,4,1)
c5=cycle5()

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

def sub(T,*args):
    L=[(g,v) for g,X in enumerate(args) for v in range(len(X))]
    def f(i,j):
        a,u=L[i]; b,v=L[j]
        return ((T[a]>>b)&1) if a!=b else ((args[a][u]>>v)&1)
    return tournament(len(L),f)

