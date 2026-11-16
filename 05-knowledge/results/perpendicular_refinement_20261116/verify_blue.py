from itertools import combinations, permutations
from collections import Counter

def hub(n,i,j):
    return i>j if (j-i)%3==0 else (j-i)%3==1

def graph(n,edits=(),comp=False):
    edits={tuple(sorted(p)) for p in edits}
    return [[i!=j and (hub(n,i,j)^(tuple(sorted((i,j))) in edits)^
                        ((not comp) and abs(i-j)>1))
             for j in range(n)] for i in range(n)]

def H(A):
    n=len(A)
    return sum(all(A[i][j]for i,j in zip(p,p[1:]))for p in permutations(range(n)))

def signature(A,d):
    ds=list(map(sum,A))
    return Counter(tuple(sorted(ds[j]for j in range(len(A))if A[i][j]))
                   for i in range(len(A))if ds[i]==d)

for m in range(1,21):
    n=6*m; M=2*m
    A=(0,3);D=(n-4,n-1);B=(3*(m-1)+1,3*m+1)
    line=[[],[A],[D],[A,B],[D,B]]
    hub0=graph(n,comp=True)
    for E in line:
        h=graph(n,E,True)
        p=list(range(n))
        for i,j in E:p[i],p[j]=p[j],p[i]
        assert all(h[i][j]==hub0[p[i]][p[j]] for i in range(n) for j in range(n))
    u,p,q,r,s=[graph(n,E)for E in line]
    assert all(u[i][j]==u[n-1-j][n-1-i] and
               p[i][j]==q[n-1-j][n-1-i] and r[i][j]==s[n-1-j][n-1-i]
               for i in range(n)for j in range(n))
    ex=[(min(map(sum,t)),max(map(sum,t)))for t in (hub0,u,p,q,r,s)]
    assert len(set(ex))==4 and ex[2]==ex[4] and ex[3]==ex[5]
    if m==1:
        assert list(map(H,[hub0,u,p,q,r,s]))==[45,33,37,37,33,33]
    elif m==2:
        sig=(3,4,5,5,6,6,7)
        assert (signature(p,7)[sig],signature(r,7)[sig])==(2,1)
    else:
        rows=[]
        for t in [p,r]:
            ids=[i for i in range(n)if sum(t[i])==2*M-m]
            assert len(ids)==3
            rows.append(sorted(sum(t[i][j]for j in ids)for i in ids))
        assert rows==[[0,1,2],[1,1,1]]
print('order 6 through 120 verified, and literal order-six Hamilton counts 45; 33; 37,37,33,33')
