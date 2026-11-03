from functools import lru_cache
from itertools import combinations, product

def need(x, message):
    if not x: raise RuntimeError(message)

def hp(rows):
    n=len(rows)
    @lru_cache(None)
    def walk(v, unused):
        if not unused: return 1
        can=rows[v]&unused
        ans=0
        while can:
            bit=can&-can; can-=bit
            ans += walk(bit.bit_length()-1, unused^bit)
        return ans
    full=(1<<n)-1
    return sum(walk(v,full^(1<<v)) for v in range(n))

def singleton(rows):
    n=len(rows)
    return tuple(row | ((1<<n) if j != n-1 else 0) for j,row in enumerate(rows)) + (1<<(n-1),)

def state(rows):
    h=hp(rows); a=hp(singleton(rows)); aa=hp(singleton(singleton(rows)))
    need((a-h)%2==0 and (aa-a)%2==0,'observer parity')
    return h,(a-h)//2,(aa-a)//2

def append(rows, block):
    n=len(rows); k=len(block)
    return tuple(row | sum(1<<(n+j) for j in range(k) if not (i==n-1 and j==0)) for i,row in enumerate(rows)) + tuple((r<<n)|((1<<(n-1)) if j==0 else 0) for j,r in enumerate(block))

def times(v,P):
    return tuple(sum(v[i]*P[i][j] for i in range(3)) for j in range(3))

P=((5,3,4),(8,7,8),(6,6,7)); B=((1,2,2),(4,3,4),(4,2,3))
blocks=[((8,1,3,6),P),((0,1,3),B)]

count=0
for n in range(1,5):
    edges=list(combinations(range(n),2))
    for code in range(1<<len(edges)):
        rows=[0]*n
        for bit,(i,j) in enumerate(edges):
            if (code>>bit)&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
        rows=tuple(rows); v=state(rows)
        for block,M in blocks:
            observed=state(append(rows,block))
            need(observed==times(v,M),('transfer',n,code,v,observed))
        count+=1
print('PASS: exact independent state DP and both transfer matrices on',count,'labelled sources n<=4')

t=(2,-1,0)
need(times(t,P)==t and times(t,B)==tuple(-x for x in t),'signed gain')
for n in range(7):
    for word in product((0,1),repeat=n):
        z=t
        for bit in word:z=times(z,(P,B)[bit])
        want=t if sum(word)%2==0 else tuple(-x for x in t)
        need(z==want,('suffix gain',word))
print('PASS: symbolic single-letter gains and all mixed words length <=6')

graphs=[(112,97,115,103,74,16,32),(252,249,194,228,204,148,32,64),
        (496,497,499,487,456,400,288,64,128),(1008,1017,1011,997,776,912,816,592,128,256),
        (996,993,994,999,975,912,800,576,128,256)]
want=[(31,14,25),(95,41,75),(97,40,75),(285,123,225),(291,120,225)]
def verify_graph(rows):
    n=len(rows)
    reach=[rows[i]|(1<<i) for i in range(n)]
    for i in range(n):
        need(rows[i]>>n==0 and not (rows[i]>>i&1),('row invalid',i))
        for j in range(i): need(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1,('arc invalid',i,j))
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i] |= reach[k]
    need(all(x==(1<<n)-1 for x in reach),'not strong')
for rows,expected in zip(graphs,want):
    verify_graph(rows)
    need(state(rows)==expected,('literal full source state',len(rows)))
    for block,M in blocks:
        g=append(rows,block)
        verify_graph(g)
        need(state(g)==times(expected,M),('literal full source extension',len(rows)))
need(want[1] == tuple(3*want[0][i]+t[i] for i in range(3)), 'first offset')
need(want[2] == tuple(3*want[0][i]+2*t[i] for i in range(3)), 'second offset')
need(want[3] == tuple(3*x for x in want[1]), 'first exact return')
need(want[4] == tuple(3*x for x in want[2]), 'second exact return')
print('PASS: five literal strong graphs, attached strong graphs, all 3-state and exact return identities')
