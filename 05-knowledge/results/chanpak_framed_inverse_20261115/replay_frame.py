"""Standalone sanity control for the framed graph note."""
from itertools import permutations, combinations
from math import comb


def H(rows):
    return sum(all(rows[u][v] for u, v in zip(p, p[1:])) for p in permutations(range(len(rows))))


def add_outlet(rows, root):
    return [row+[i != root] for i,row in enumerate(rows)] + [[i == root for i in range(len(rows))]+[0]]


def chart(rows):
    n = len(rows)
    return {(i,j) for i,j in combinations(range(n-1),2) if rows[j+1][i]}


def Ea(rows, u):
    n = len(rows)
    answer = 0
    for r in range(2, n+1, 2):
        for p in permutations((i for i in range(n) if i != u),r-1):
            left = (u,)+p
            if all(rows[a][b] for a,b in zip(left,left[1:])):
                right = sorted(set(range(n))-set(left))
                answer += H([[rows[i][j] for j in right] for i in right])
    return answer


for n in range(3,8):
    for b in range(3,n+1):
        rows = [[i<j for j in range(n)] for i in range(n)]
        rows[0][b-1],rows[b-1][0]=False,True
        assert chart(rows) == {(0,b-2)}
        assert chart(add_outlet(rows,0)) == {(0,b-2),(0,n-1)}
        e = 2**(n-2) - (2**(n-b-1) if n>b else 1)
        h = 1+2**(b-2)
        assert Ea(rows,0)==e
        assert H(rows)==h
        assert H(add_outlet(rows,0))==h+2*e


def fence(beta):
    return sum(all((a>b)==bit for (a,b), bit in zip(zip(p,p[1:]),beta))
               for p in permutations(range(len(beta)+1)))


for m in range(1,5):
    pairs=list(combinations(range(m),2))
    for code in range(1<<len(pairs)):
        G={frozenset(pair) for i,pair in enumerate(pairs) if code>>i&1}
        def wf(seq):
            return fence(tuple(frozenset({a,b}) in G for a,b in zip(seq,seq[1:])))
        mh,me=0,0
        for order in permutations(range(m)):
            order_index={vertex:i for i,vertex in enumerate(order)}
            rows=[[i!=j and ((order_index[i]<order_index[j]) != (frozenset({i,j}) in G))
                   for j in range(m)] for i in range(m)]
            mh+=H(rows)
            me+=Ea(rows,0)
        assert mh==sum(wf(seq) for seq in permutations(range(m)))
        expect=0
        for r in range(2,m+1,2):
            for p in permutations(range(1,m),r-1):
                left=(0,)+p
                for right in permutations(sorted(set(range(m))-set(left))):
                    expect+=comb(m,r)*wf(left)*wf(right)
        assert me==expect
print("Passed through staircase order 7 and all fence graphs through order 4.")
