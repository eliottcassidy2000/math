from itertools import combinations

P = (2,3,5,7,11,13,17,19,23,29)
R = (3,5,11,13,17,19,23,29,31,37)
F = ((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),
     (77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))

def require(ok):
    if not ok:
        raise AssertionError

def out_rows(n, mask):
    A = [0]*n
    for j,(a,b) in enumerate(combinations(range(n),2)):
        a,b = (a,b) if mask>>j&1 else (b,a)
        A[a] |= 1<<b
    return A

def ham(A):
    n=len(A)
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        d[1<<i][i]=1
    for s in range(1,1<<n):
        for v,h in enumerate(d[s]):
            b=A[v]&~s
            while h and b:
                l=b&-b;b-=l
                d[s|l][l.bit_length()-1]+=h
    return sum(d[-1])

def strong(A):
    for v in range(len(A)):
        s=1<<v
        for _ in A:
            s2=s
            for u,a in enumerate(A):
                if s>>u&1: s2 |= a
            s=s2
        if s!=(1<<len(A))-1: return False
    return True

def fac(n):
    V = [0]*len(P)
    for j,p in enumerate(P):
        while n%p==0:
            n//=p
            V[j]+=1
    require(n==1)
    return tuple(V)

def product(vals, exp):
    n=1
    for a,b in zip(vals,exp):
        n*=a**b
    return n

def find():
    W={}
    for n in range(3,7):
        for mask in range(1<<(n*(n-1)//2)):
            a=out_rows(n,mask)
            if strong(a):
                h=ham(a)
                if h in R: W.setdefault(h,(n,mask))
                if len(W)==len(R): return W
    raise AssertionError(W)

def run():
    W=find()
    print('tokens',[(p,R[i],*W[R[i]]) for i,p in enumerate(P)])
    tape=fac(2)
    raw=2
    rules=[(fac(a),fac(b)) for a,b in F]
    clean=[]
    for t in range(1,6000):
        j=next(j for j,(a,b) in enumerate(F) if raw*a%b==0)
        k=next(k for k,(a,b) in enumerate(rules) if all(x>=y for x,y in zip(tape,b)))
        require(j==k)
        a,b=F[j]
        raw=raw*a//b
        r,s=rules[j]
        tape=tuple(e+x-y for e,x,y in zip(tape,r,s))
        require(product(P,tape)==raw)
        if all(x==0 for x in tape[1:]):
            require(product(R,tape)==3**tape[0])
            clean.append((t,tape[0]))
    print('checkpoints through step 5999',clean)
    require(clean==[(19,2),(69,3),(281,5),(710,7),(2375,11),(3893,13)])
if __name__=='__main__':run()
