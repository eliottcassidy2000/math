from math import comb, ceil
from collections import defaultdict
from itertools import permutations
from fractions import Fraction

def sides(d):
    N=d*d+1;m=N//2
    return [(a,d-a) for a in range(d+1) if (2*a-d)**2+(2*a-d)==d-1]

def quotient_rows():
    trace=[s for s in range(126) if (125+2*s-57+8*124)%15==0]
    assert trace==list(range(10,126,15))
    return trace

def path_counts_5():
    p=list(permutations(range(5)))
    def arc(x,y):return (y-x)%5 in (1,2)
    c=[0]*6
    for p0 in p:
        bad=sum(not arc(a,b) for a,b in zip(p0,p0[1:]))
        # Partition the permutation into an ordered list of r nonempty paths.
        for r in range(1,6):
            if bad<=r-1:c[r]+=comb(4-bad,r-1-bad)
    assert c[1]==15
    return c

def lex25():
    c=path_counts_5()
    cur={}
    for i in range(5):
        s=[0]*5;s[i]=1
        cur[tuple(s),i]=1
    total=0
    for k in range(1,26):
        nxt=defaultdict(int)
        for (s,i),v in cur.items():
            if all(s):
                w=v
                for r in s:w*=c[r]
                total+=w
            for j in ((i+1)%5,(i+2)%5):
                if s[j]==5:continue
                t=list(s);t[j]+=1
                nxt[tuple(t),j]+=v
        cur=nxt
    assert total%25==0
    return total

def rank(A):
    A=[list(map(Fraction,row)) for row in A]; r=0
    for j in range(len(A[0])):
        i=next((i for i in range(r,len(A)) if A[i][j]),None)
        if i is None: continue
        A[i],A[r]=A[r],A[i]
        for k in range(r+1,len(A)):
            if A[k][j]:
                f=A[k][j]/A[r][j]
                A[k]=[x-f*y for x,y in zip(A[k],A[r])]
        r+=1
    return r

def controls():
    for q,inc in [(1,1),(5,5)]:
        if q==1:
            N=10; p=5
            A=[set() for _ in range(N)]
            for x in range(5):
                for y in ((x+1)%5,(x-1)%5,5+x):A[x].add(y)
                for y in (x,5+(x+2)%5,5+(x-2)%5):A[5+x].add(y)
        else:
            N=50;p=25
            A=[set() for _ in range(N)]
            V=[(z,i,x) for z in range(2) for i in range(5) for x in range(5)]
            for v,(z,i,x) in enumerate(V):
                for w,(z1,j,y) in enumerate(V):
                    if z==z1 and i==j and (y-x)%5 in (1+z,4-z):A[v].add(w)
                    if z==0 and z1==1 and (y-x-i*j)%5==0:
                        A[v].add(w);A[w].add(v)
        C=[[int(j in A[i]) for j in range(p,2*p)] for i in range(p)]
        assert rank(C)==(5 if q==1 else 21)
        d=inc+2
        for i,row in enumerate(A):
            assert len(row)==d and sum((j<p)!=(i<p) for j in row)==inc
            for j in range(i):assert len(row&A[j])==int(j not in row)

if __name__=='__main__':
    controls()
    assert sides(3)==[(2,1)]
    assert sides(7)==[(2,5)]
    assert sides(57)==[(32,25)]
    assert [120-(2*c)//15+1 for c in (30,120,210)]==[117,105,93]
    assert ceil(Fraction(40000*4,225))==712
    print('125-cell quotient: possible numbers of cubic cells',quotient_rows())
    print('cyclic T5 path-cover counts:',path_counts_5())
    print('lexicographic T25 exact H:',lex25())
