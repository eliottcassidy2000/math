"""Exact observer identities in the inherited six-state funnel compiler.

This is a bounded algebra search, not a tournament census.  Coordinates are
(H,E(1),S(2),P(2),E(2),S(3)); every nonempty word ends in 1.
"""
from functools import lru_cache
from itertools import product
from pathlib import Path
import json


@lru_cache(None)
def coeff(word, m):
    t = len(word)
    L = (m-1)**t + sum(m**i*(m-1)**(t-i-1) for i,b in enumerate(word) if b)
    D = (m-1)**t + sum((m-1)**i*m**(t-i-1) for i,b in enumerate(word) if not b)
    states = {(0,0): 1}
    for b in word:
        nxt = {}
        for (last,started), z in states.items():
            for key,v in (((last,started),(m-1)*z),
                          ((1 if b else 2,started),z)):
                nxt[key] = nxt.get(key,0)+v
            if started or not b:
                key = (last,1)
                nxt[key] = nxt.get(key,0)+z
        states = nxt
    return L,D,sum(z for (last,_),z in states.items() if last != 2)


@lru_cache(None)
def matrix(word):
    t = len(word)
    l1,d1,b1 = coeff(word,1)
    l2,d2,b2 = coeff(word,2)
    l3,d3,_ = coeff(word,3)
    return ((b1,-d1,l2,0,0,0),
            (l1,-int(not t),0,0,0,0),
            (d2,0,2**t,0,0,0),
            (0,0,0,2*b2,-2*d2,2*l3),
            (0,0,0,l2,-1,0),
            (0,0,0,d3,0,3**t))


def rowmul(row, M):
    return tuple(sum(row[i]*M[i][j] for i in range(6)) for j in range(6))


def ready(word):
    M = matrix(word)
    return tuple(M[0][i]+M[4][i] for i in range(6))


def words(bound):
    return [()] + [p+(1,) for n in range(1,bound+1) for p in product((0,1), repeat=n-1)]


def run(bound=6):
    W=words(bound)
    bank={}
    for inside in W:
        for outside in W:
            row=rowmul(ready(outside),matrix(inside))
            recipe=(inside,outside)
            if row not in bank or sum(map(len,recipe)) < sum(map(len,bank[row])):
                bank[row]=recipe
    hits=[]
    for row,src in bank.items():
        dst=bank.get(tuple(3*x for x in row))
        if dst is not None:
            hits.append(dict(source=src,target=dst,source_row=row,
                             target_row=tuple(3*x for x in row)))
    hits.sort(key=lambda z: sum(map(len,z['source']))+sum(map(len,z['target'])))
    result=dict(status='FINITE-EXACT',scope=f'Every two-funnel address with each word length <= {bound}; empty or ending1',
                words=len(W),distinct_observers=len(bank),factor3_identities=hits)
    path=Path(__file__).with_name('intertwiner_receipt.json')
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(words=len(W),observers=len(bank),hits=len(hits),first=hits[:8]),indent=2))


if __name__=='__main__':
    run()
