from itertools import product
from counter36 import H
from language_certificate import tr, acc, start

def rows(w):
    n=len(w)+2
    M=[2,0]+[sum(1<<j for j in range(i+1,n)) for i in range(2,n)]
    for i,d in enumerate(map(int,w),2):
        for h in (0,1):
            if d>>(1-h)&1: M[h]|=1<<i
            else: M[i]|=1<<h
    return M
def literal(w):
    M=rows(w);n=len(M)
    V=[[0]*n for _ in range(1<<n)]
    for j in range(n):V[1<<j][j]=1
    for mask,v in enumerate(V):
        for i,c in enumerate(v):
            if c:
                for j in range(n):
                    if mask>>j&1==0 and M[i]>>j&1:V[mask|1<<j][j]+=c
    return sum(V[-1])
assert all(literal('1'*a+'2'+'1'*b)==1+3*2**(a+b) for a,b in product((1,2,3),repeat=2))
assert all(literal('1'*a+'20')==9*2**a-3 for a in (1,2,3))
print('independent subset DP: 9 A-board and 3 B-board coefficients')

from itertools import permutations,groupby
def canonical(w):
    M=rows(w)
    perm=sorted(range(len(M)),key=lambda i:M[i].bit_count())
    ans=None
    for p0 in product(*(tuple(permutations(g)) for _,g in groupby(perm,lambda i:M[i].bit_count()))):
        p=sum(p0,())
        u=tuple(sum((M[i]>>j&1)<<k for k,j in enumerate(p)) for i in p)
        if ans is None or u<ans:ans=u
    return ans
def boards(m):
    L={}
    for w in map(''.join,product('0123',repeat=m)):
        s=start
        for d in map(int,w):
            s=tr(s,d)
            if s is None:break
        if acc(s):L[canonical(w)]=H(w)
    return L
from collections import Counter
for n,out in [(5,{13:1,15:2}),(6,{25:2,33:2}),(8,{97:4,141:2}),(10,{385:6,573:2})]:
    assert Counter(boards(n-2).values())==out
    print(n, 'hub-free class multiplicities by count:',out)
