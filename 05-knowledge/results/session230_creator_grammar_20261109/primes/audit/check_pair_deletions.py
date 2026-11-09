# Independent narrow check of section 7 definitions; active under python -O.
def require(ok, message):
    if not ok:
        raise ValueError(message)
def cyclic(r):
    n=2*r+1
    return [[0<(j-i)%n<=r for j in range(n)] for i in range(n)]
def staircase(r):
    V=[('A',i) for i in range(r)]+[('B',j) for j in range(r+1)]
    def arc(x,y):
        a,i=x;b,j=y
        if a==b:return i>j if a=='A' else i<j
        return j<=i if a=='A' else i>j
    return [[arc(x,y) for y in V] for x in V]
def alternating(r):
    n=2*r+1
    def arc(i,j):
        if i==j:return False
        if i==2*r:return j%2==0
        if j==2*r:return i%2==1
        return i<j
    return [[arc(i,j) for j in range(n)] for i in range(n)]
def restrict(A,V):return [[A[i][j] for j in V] for i in V]
for r in range(2,21):
    n=2*r+1
    for name,make,deleted in [('T',cyclic,{r,2*r}),('U',staircase,{r-1,2*r}),('W',alternating,{2*r-2,2*r-1})]:
        require(restrict(make(r),[v for v in range(n) if v not in deleted])==make(r-1),f'{name}, r={r}')
print('T/U/W PAIR DELETIONS PASS for r=2,...,20')
