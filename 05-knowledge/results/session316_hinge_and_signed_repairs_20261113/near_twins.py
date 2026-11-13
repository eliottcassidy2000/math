from itertools import product

def rows(P,Q,cyclic):
    m=len(P); n=m+3
    A=[0]*n
    def put(i,j,b=True):
        if not b:i,j=j,i
        A[i]|=1<<j
    for i in range(m):
        for j in range(i+1,m):
            put(i,j)
        put(m,i,P[i]); put(m+1,i,Q[i]); put(m+2,i,Q[i])
    a,p,q=m,m+1,m+2
    put(p,a); put(a,q); put(p,q,not cyclic)
    return A

def path_count(A):
    n=len(A); N=(1<<n)-1
    D=[[0]*n for _ in range(N+1)]
    for i in range(n):
        D[1<<i][i]=1
    for S in range(1,N+1):
        for j,val in enumerate(D[S]):
            if val:
                rest=A[j]&~S
                while rest:
                    b=rest&-rest; rest^=b
                    D[S|b][b.bit_length()-1]+=val
    return sum(D[N])

def mixed(P,Q):
    m=len(P)
    if not m:return 1
    def w(i,j):
        if j==i+1:return 1
        if j<i:return 2**(i-j-1)
        return 0
    return P[0]+(1-Q[-1])+sum(w(i,j) for i in range(m) for j in range(m) if P[j] and not Q[i])

def reduced(P,Q):
    P=list(P);Q=list(Q)
    while P and P[0]==Q[0]==0:
        P.pop(0);Q.pop(0)
    while P and P[-1]==Q[-1]==1:
        P.pop();Q.pop()
    return tuple(P),tuple(Q)

if __name__=='__main__':
    controls={}
    for m in range(6):
        for P in product((0,1),repeat=m):
            for Q in product((0,1),repeat=m):
                a,b=(path_count(rows(P,Q,c)) for c in (False,True))
                d=mixed(P,Q)+mixed(Q,P)
                if b-a!=d:
                    raise ValueError((P,Q,a,b,d))
                if d==2:
                    R=reduced(P,Q)
                    if len(R[0])>2 or (R in controls and controls[R]!=(a,b)):
                        raise ValueError((P,Q,a,b,R))
                    controls[R]=(a,b)
    print(controls)
