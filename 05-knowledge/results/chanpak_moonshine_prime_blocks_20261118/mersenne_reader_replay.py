"""Literal replay of cyclic sandwich and ordered-slot reader identities."""
from itertools import combinations
from math import comb

def h(rows,S):
    if not S:
        return 1
    D=[{} for _ in range(S+1)]
    for i in range(len(rows)):
        if (1<<i)&S:
            D[1<<i][i]=1
    for p,row in enumerate(D):
        for i,w in row.items():
            A=rows[i]&S&~p
            while A:
                b=A&-A
                A^=b
                new=D[p|b]
                j=b.bit_length()-1
                new[j]=new.get(j,0)+w
    return sum(D[S].values())

def p2(rows):
    F=(1<<len(rows))-1
    counts=[h(rows,S) for S in range(F+1)]
    return counts[F],sum(counts[S]*counts[F^S] for S in range(F+1))

def gluing(A,r):
    n=len(A)
    b=((1<<r)-1)<<n
    c=1<<(n+r)
    return [a|b for a in A]+[c|(b&~((1<<(i+1))-1)) for i in range(n,n+r)]+[(1<<n)-1]

def gen(n):
    P=list(combinations(range(n),2))
    for k in range(1<<len(P)):
        A=[0]*n
        for t,(i,j) in enumerate(P):
            if not ((k>>t)&1): i,j=j,i
            A[i]|=1<<j
        yield A

def check(b,msg):
    if not b:
        raise AssertionError(msg)

m=0
for n in range(1,6):
    for rows in gen(n):
        a,b=p2(rows)
        check((b-(-1)**n*a)%3==0,("mod3",rows,a,b))
        if n>=4 and n%2==0:
            check(b%8==0,("mod8",rows,b))
        for r in range(1,5):
            C=gluing(rows,r)
            t=h(C,(1<<len(C))-1)
            check(t==2*a+((1<<r)-1)*(b-a),("gluing", rows,r,a,b,t))
            m+=1
print("verified glued orders and transfer",m)
for A in ([0],[2,0],[2,4,1],[8,1,3,6]):
    print("source",A,p2(A),"outputs",[h(C,(1<<len(C))-1) for C in (gluing(A,r) for r in range(1,5))])
