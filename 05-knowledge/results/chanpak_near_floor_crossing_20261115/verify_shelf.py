"""Literal controls for the crossing-shelf family."""
from itertools import pairwise

OLD = [4030,4024,4090,1008,4064,4032,3971,3840,3584,3072,2056,8]
MARK = 4008

def need(ok, detail):
    if not ok:
        raise RuntimeError(detail)

def family(r):
    order = list(range(7)) + [7]*r + [10,11]
    rows = [sum(1<<j for j,w in enumerate(order)
                if (i<j if u==w else bool(OLD[u]>>w&1)))
            for i,u in enumerate(order)]
    mark = sum(1<<i for i,u in enumerate(order) if MARK>>u&1)
    return rows, mark

def ext(rows, x):
    n = len(rows)
    return [y | (1<<n if not(x>>i&1) else 0)
            for i,y in enumerate(rows)] + [x]

def paths(rows):
    n = len(rows)
    def visit(last, todo, pref):
        if not todo:
            yield pref
        legal = rows[last] & todo
        while legal:
            bit = legal & -legal
            legal ^= bit
            nxt = bit.bit_length()-1
            yield from visit(nxt, todo^bit, pref+(nxt,))
    for i in range(n):
        yield from visit(i, ((1<<n)-1)^(1<<i), (i,))

def strong(rows):
    n=len(rows)
    for i in range(n):
        got=1<<i
        while True:
            nxt=got
            for j in range(n):
                if got>>j&1:
                    nxt|=rows[j]
            if nxt==got:
                break
            got=nxt
        if got!=(1<<n)-1:
            return False
    return True

def check():
    for r in range(1,6):
        rows,X=family(r)
        need(strong(rows),'strong')
        need(strong(ext(rows,X)),'strong extension')
        ps=list(paths(rows))
        double=next((p for p in ps if sum((X>>w&1) and not(X>>u&1)
                                         for w,u in pairwise(p))>=2),None)
        need(double, ('double',r))
        a,b=len(ps),sum(1 for _ in paths(ext(rows,X)))
        need((a,b)==(267*2**r+1,528*2**r+45),(r,a,b))
        print(r,a,b,b-2*a,','.join(map(str,double)))
if __name__=='__main__':
    check()
