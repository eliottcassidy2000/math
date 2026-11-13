def rows(w,k=2):
    n=len(w)+k
    def arc(i,j):
        if i<k and j<k:return i<j
        if i<k:return w[j-k]>>i&1
        if j<k:return not(w[i-k]>>j&1)
        return i<j
    return [sum(1<<j for j in range(n) if i!=j and arc(i,j)) for i in range(n)]
def subtable(A):
    n=len(A);D=[[0]*n for _ in range(1<<n)]
    for j in range(n):D[1<<j][j]=1
    for S in range(1,1<<n):
        for j,c in enumerate(D[S]):
            if c:
                N=A[j]&~S
                while N:
                    b=N&-N;N^=b;D[S|b][b.bit_length()-1]+=c
    return D
def strong(A):
    n=len(A)
    for i in range(n):
        S=1<<i
        for _ in A:
            for j,b in enumerate(A):
                if S>>j&1:S|=b
        if S!=(1<<n)-1:return False
    return True
def count(w):
    return sum(subtable(rows(w))[-1])
def profile(w):
    A=rows(w); W=[r>>1 for r in A[1:]]
    N=(1<<len(W))-1; U=A[0]>>1; D=subtable(W)
    h=[1]+list(map(sum,D[1:]))
    E=sum(h[N^S]*sum(x for j,x in enumerate(D[S]) if not U>>j&1)
          for S in range(1,N+1))
    R=[sum(1<<j for j,r in enumerate(W) if r>>i&1) for i in range(len(W))]
    s=sum(x for j,x in enumerate(subtable(R)[N]) if U>>j&1)
    return h[-1],E,s,count(w)
def K(B):
    r=B//2
    return (K(r)+(B%2)*(r+(r%2))) if B else 0
for word,want in (
    ("021102",(9,245,0,91)),("0102330",(9,245,0,93)),("0102101",(9,245,0,95)),
    ("120",(3,12,3,13)),("312",(3,12,3,15))):
    w=tuple(map(int,word)); got=profile(w)
    if got!=want or not strong(rows((1,)+w)):raise ValueError((word,got,want))
    h,E,s,c=got
    for B in range(16):
        P=tuple(map(int,bin(B)[2:]))
        if count(P+w)!=c+(h+E)*B-h*K(B)-s*(B%2):raise ValueError((word,B))
    print(word,got)
