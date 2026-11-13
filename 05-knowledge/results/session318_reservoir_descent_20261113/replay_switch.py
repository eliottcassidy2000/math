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

rowspecs = {
 '2':(1,2,1,3), '30':(3,3,2,5), '320':(5,16,3,15),
}
for word,want in rowspecs.items():
    w=tuple(map(int,word))
    if profile(w)!=want or not strong(rows((1,)+w)): raise ValueError(word)
def check(word,B,N,direct):
    h,e,s,c=rowspecs[word]
    F=c+(h+e)*B-h*K(B)-s*(B%2)
    if F!=N:raise ValueError((word,B,F,N))
    if direct:
        w=tuple(map(int,bin(B)[2:]+word))
        if count(w)!=N or not strong(rows(w)):raise ValueError(w)
for k in range(4,100):
    v=9*2**(k-2)+9
    check('2',2**k+3,v,k<=5)
    check('30',2**(k-1)+1,v,k<=5)
    v=21*2**k
    check('320',2**k,v+15,k<=5)
    check('30',2**(k+2)+2,v+17,k<=5)
    check('2',2**(k+3)+6,v+19,k<=5)
print('three physical profiles; 480 exact chart checks; 10 literal strong path counts')
