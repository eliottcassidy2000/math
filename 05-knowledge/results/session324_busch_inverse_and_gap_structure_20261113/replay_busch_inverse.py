"""Check exact seam formulas by direct path recursion and explicit cover sums."""
from functools import lru_cache
from random import Random

def need(cond):
    if not cond:
        raise ValueError("replay failed")

def data(M, root):
    n=len(M); full=(1<<n)-1; active=full^(1<<root)
    @lru_cache(None)
    def start(S,i):
        if S==(1<<i): return 1
        S^=1<<i
        return sum(start(S,j) for j in range(n) if S>>j&1 and M[i][j])
    @lru_cache(None)
    def end(S,i):
        if S==(1<<i): return 1
        S^=1<<i
        return sum(end(S,j) for j in range(n) if S>>j&1 and M[j][i])
    def h(S):
        return sum(start(S,i) for i in range(n) if S>>i&1) if S else 1
    L=0
    for Q in range(full+1):
        if Q&~active: continue
        s=sum(start(Q,i) for i in range(n) if Q>>i&1 and M[root][i]) if Q else 1
        L+=h(active^Q)*s
    return h(full),h(active),end(full,root),L

def NT(j):
    n=j+2
    return [[int(i<k and (i,k)!=(0,n-1) or (i,k)==(n-1,0)) for k in range(n)] for i in range(n)]

def join(A,B):
    a,b=len(A),len(B)
    T=[[0]*(a+b-1) for _ in range(a+b-1)]
    for i in range(a):
        for j in range(a): T[i][j]=A[i][j]
    for i in range(b):
        for j in range(b):
            T[0 if i==0 else a+i-1][0 if j==0 else a+j-1]=B[i][j]
    for i in range(1,a):
        for j in range(a,a+b-1): T[i][j]=1
    return T

def reroot(T,r):
    order=[r]+[i for i in range(len(T)) if i!=r]
    return [[T[i][j] for j in order] for i in order]

def strong(T):
    for row in (T, tuple(zip(*T))):
        reach={0}
        while True:
            nxt=reach|{j for i in reach for j,x in enumerate(row[i]) if x}
            if len(nxt)==len(reach): break
            reach=nxt
        if len(reach)!=len(T): return False
    return True

def main():
    rng=Random(324)
    def draw(n):
        while True:
            M=[[0]*n for _ in range(n)]
            for i in range(n):
                for j in range(i+1,n):
                    M[i][j]=rng.randrange(2); M[j][i]=1-M[i][j]
            if strong(M): return M
    for _ in range(160):
        A,B=draw(rng.randrange(3,7)),draw(rng.randrange(3,7))
        H,h,E,L=data(A,0)
        J,k,S,R=data(tuple(zip(*B)),0)
        T=join(A,B)
        need(data(T,0)[0]==(L-h)*(R-k)+(H-E)*k+h*(J-S)+E*S)
        need(data(T,0)[0]-H*J==(H-E)*(R-J)+(J-S)*(L-H)+(L-H-h+E)*(R-J-k+S))
        degree=sum(A[0]); A0=[r[1:] for r in A[1:]]
        need((L==H)==(degree==1 and not strong(A0)))
        need(L-H>=2 or degree==1)
        j=rng.randrange(1,5); T=join(A,NT(j))
        need(data(T,0)[::3]==(H+(1<<j)*L,(1<<(j+1))*L-h))
        N,_,_,Lnew=data(reroot(T,len(T)-1),0)
        need(N==Lnew)
    for j in range(1,10):
        X=NT(j)
        for i in range(j+2):
            T=[[0]+[int(i==k) for k in range(j+2)]]+[[int(i!=k)]+row for k,row in enumerate(X)]
            H,_,_,L=data(T,0)
            if i==0: want=(3*2**j-1,1)
            elif i==j+1: want=(2**(j+1)+1,1)
            else: want=(2**j+1-2**(i-1)+2**(j+1-i)+2**i*3**(j-i)-int(i==1),2**(i-1))
            need((H,L-H)==want)
    print("Busch seam: 160 direct overlapping pairs, 160 root resets, 72 dyadic menu roots passed")
if __name__=="__main__": main()
