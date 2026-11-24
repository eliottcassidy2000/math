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

done=0
for T,A in [(tri,(one,one,one)),(tri,(two,one,tri)),
            (tri,(tri,two,two)),(c5,(one,one,one,one,one)),
            (c5,(one,two,one,tri,one)),(c5,(two,one,two,one,two))]:
    R=sub(T,*A)
    for k in (1,2,3,4):
        U=[profile(x,k+1)-2*profile(x,k)+profile(x,k-1) if len(x)>1 else 1 for x in A]
        assert profile(R,k)>=prod(U)
        assert all(u*3*(k+1)**2>=2*profile(a,k+1) for a,u in zip(A,U))
        done+=1
print(done,'literal strong-template color-transfer inequalities checked')
