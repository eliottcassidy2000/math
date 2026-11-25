"""Three saturated branch pairs, with a nonlinear order-56 Latin board on each."""
from collections import Counter
from itertools import combinations
from math import gcd

N = 56
PIVOT, A, B = 56, 0, 28
SLOTS = [h for h in range(N) if h not in (A,B)]
IDENTITY = tuple(range(N))
TRIPLES = ((9,10,19),(8,12,20),(7,11,18),(6,15,21),(5,22,27),
           (4,13,17),(3,23,26),(2,14,16),(1,24,25))
SHORT = ((1,5,0),(3,10,0),(8,29,0))
TRADES = (
    (40,2,14),(19,29,40),(24,38,33),(18,26,42),(45,11,8),
    (11,1,16),(10,38,44),(10,38,15),(7,53,3),(15,10,1),
    (25,37,25),(22,24,47),(46,32,33),(50,52,52),(12,48,51),
    (35,32,53),(12,17,2),(2,9,29),(53,19,35),(30,18,44),(1,50,0),(25,21,24),(51,46,32),
    (40,14,43),(14,10,22),(54,30,41),(39,24,41),(29,45,25))


def inverse(row):
    result=[None]*N
    for x,y in enumerate(row):
        result[y]=x
    assert None not in result
    return tuple(result)


def compose(g,f):
    return tuple(g[y] for y in f)


def translation(t):
    return tuple((x+t)%N for x in IDENTITY)


def make_chart(trades=TRADES):
    s={}
    for u,v,w in TRIPLES:
        triple=(u,v,-w)
        for x,y in zip(triple,triple[1:]+triple[:1]):
            for sign in (-1,1):
                s[sign*x%N]=-sign*y%N
    assert set(s)==set(SLOTS)==set(s.values())=={(h-s[h])%N for h in s}
    p={h:list(translation(-h)) for h in SLOTS}
    q={h:list(translation(-s[h])) for h in SLOTS}
    for h,k,start in trades:
        ip,iq=inverse(p[k]),inverse(q[k])
        orbit={start}
        queue=[start]
        while queue:
            x=queue.pop()
            for y in (ip[p[h][x]],iq[q[h][x]]):
                if y not in orbit:
                    orbit.add(y)
                    queue.append(y)
        assert len(orbit)<N
        for x in orbit:
            p[h][x],p[k][x]=p[k][x],p[h][x]
            q[h][x],q[k][x]=q[k][x],q[h][x]
    maps={(PIVOT,j):IDENTITY for j in IDENTITY}
    maps[A,B]=translation(B)
    for h in SLOTS:
        maps[h,A]=tuple(p[h])
        maps[h,B]=tuple(q[h])
    maps.update({(j,i):inverse(f) for (i,j),f in maps.items()})
    return maps


def routes(maps,i,j):
    rows=[maps[i,j]] if (i,j) in maps else []
    rows.extend(compose(maps[k,j],maps[i,k]) for k in range(57) if (i,k) in maps and (k,j) in maps)
    return rows


def rank(board,prime):
    # Columns represent values of f on the sheets; constants are automatically in the kernel.
    base={}
    for row in board:
        for x in IDENTITY:
            a=[0]*N
            for j,sign in ((row[x],1),(x,-1),(row[0],-1),(0,1)):
                a[j]=(a[j]+sign)%prime
            for j in IDENTITY:
                coef=a[j]
                if not coef:
                    continue
                if j not in base:
                    inv=pow(coef,-1,prime)
                    base[j]=[(inv*v)%prime for v in a]
                    break
                a=[(v-coef*w)%prime for v,w in zip(a,base[j])]
            if len(base)==55:
                return 55
    return len(base)


def symmetric_group_certificate(board,word):
    assert any(row==translation(row[0]) and gcd(56,row[0])==1 for row in board)
    product=IDENTITY
    for r in word:
        product=compose(board[r],product)
    unseen=set(IDENTITY)
    cycles=[]
    while unseen:
        y=min(unseen)
        cycle=[]
        while y in unseen:
            unseen.remove(y)
            cycle.append(y)
            y=product[y]
        cycles.append(cycle)
    special=[cycle for cycle in cycles if len(cycle)%3==0]
    assert len(special)==1 and len(special[0])==3
    i,j,k=sorted(special[0])
    assert gcd(N,j-i,k-i)==1
    return sorted(map(len,cycles)),(i,j,k)


def check(trades):
    maps=make_chart(trades)
    before=make_chart(())
    results=[]
    for i,j in combinations(range(57),2):
        board=routes(maps,i,j)
        full=i in (PIVOT,A,B) and j in (PIVOT,A,B)
        assert len(board)==(56 if full else 3)
        assert all(len(set(column))==len(board) for column in zip(*board))
        if full:
            rs=[rank(board,p) for p in (2,7)]
            assert rs==[55,55]
            raw=routes(before,i,j)
            edits=sum(x!=y for f,g in zip(raw,board) for x,y in zip(f,g))
            if trades==TRADES:
                words={(0,28):(1,12,16),(0,56):(1,10,12),(28,56):(2,11,32)}
                print((i,j),"symmetric group witness",symmetric_group_certificate(board,words[i,j]))
                defects=[sum(N-max(Counter((r[x]-x)%v for x in IDENTITY).values()) for r in board)
                         for v in (2,7,56)]
                print((i,j),"minimum standard-coordinate defects at 2,7,56",defects)
                assert min(defects)>=320
            results.append(((i,j),rs,edits))
    print(len(trades),"trades; pair, character ranks and changed cells",results)


if __name__=="__main__":
    check(SHORT)
    check(TRADES)
