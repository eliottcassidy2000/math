"""Repunit 57 -> PBD(54;5,7,8) -> six-column OA. Python stdlib."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path

N=54
C=list(combinations(range(6),2))


def projective(v):
    first=next(x for x in v if x)
    return tuple(x*pow(first,-1,7)%7 for x in v)


def blocks():
    P=sorted({projective(x) for x in product(range(7),repeat=3) if any(x)})
    assert len(P)==57
    # The three removed points are on z=0.
    D={(1,0,0),(1,1,0),(0,1,0)}
    V=[p for p in P if p not in D]
    B=[tuple(i for i,p in enumerate(V) if sum(a*b for a,b in zip(p,w))%7==0)
       for w in P]
    assert Counter(map(len,B))=={5:1,7:21,8:35}
    assert Counter(pair for block in B for pair in combinations(block,2))=={
        pair:1 for pair in combinations(range(N),2)}
    return B


def gf8(a,b):
    ret=0
    while b:
        if b&1:ret^=a
        b>>=1
        a<<=1
        if a&8:a^=0b1011  # x^3+x+1
    return ret


def array54():
    B=blocks()
    core,=(b for b in B if len(b)==5)
    result=[(x,)*6 for x in range(N) if x not in core]
    for b in B:
        q=len(b)
        for a,t in product(range(q),repeat=2):
            if q==5:
                # Include the projective direction infinity, and all 25 tuples.
                vals=[(a+i*t)%5 for i in range(5)]+[t]
            elif t:
                # Delete exactly the diagonal tuples; the common global diagonal is above.
                vals=[(a+i*t)%7 if q==7 else a^gf8(i,t) for i in range(6)]
            else:continue
            result.append(tuple(b[x] for x in vals))
    assert len(result)==N*N
    for i,j in C:assert len({(r[i],r[j]) for r in result})==N*N
    return result


def read_tables(path):
    seq=iter(int(w,16) for w in Path(path).read_text().split())
    table={ij:tuple(next(seq) for _ in range(N)) for ij in C}
    assert next(seq,None) is None
    return table


def checked(rows,table):
    for i,j in C:
        used={(r[i],r[j]) for r in rows}
        assert len(used)==len(rows) and all(table[i,j][x]>>y&1 for x,y in used)


if __name__=='__main__':
    import sys, random
    arr=array54()
    print("PBD block inventory", sorted(Counter(map(len,blocks())).items()))
    print("OA rows",len(arr),"verified ordered pair coverage",len(arr)*len(C))
    for case,file in enumerate(sys.argv[1:]):
        p=Path(file)
        table=read_tables(p)
        rr=random.Random(5700+case)
        best=[]
        for t in range(250):
            ps=[rr.sample(range(N),N) if t else list(range(N)) for _ in range(6)]
            rows=[tuple(ps[i][r[i]] for i in range(6)) for r in arr]
            bad=[sum(not(table[i,j][r[i]]>>r[j]&1) for i,j in C) for r in rows]
            good=[r for r,k in zip(rows,bad) if not k]
            if len(good)>len(best):best=good
        checked(best,table)
        print(case,"best 250 OA relabelings",len(best))
        p.with_suffix('.seed').write_text(''.join(' '.join(map(str,r))+'\n' for r in best))
