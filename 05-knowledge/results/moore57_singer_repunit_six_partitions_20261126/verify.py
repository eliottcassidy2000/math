"""Checks stored packings against the exact fifteen boundary boards."""
from collections import Counter
from pathlib import Path
from itertools import combinations
from singer_seed import C, N, read_tables, checked
from highholes import high, bits

PATH=Path(__file__).parent
ABC="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"


def residual(p, rows):
    a={ij:list(b) for ij,b in p.items()}
    for z in rows:
        for i,j in C:a[i,j][z[i]]^=1<<z[j]
    return a


def remaining_cliques(a):
    def dfs(z,S):
        i=len(z)
        if i==6:return 1
        return sum(dfs(z+[x],[w&a[i,j][x] for j,w in zip(range(i+1,6),S[1:])])
                   for x in bits(S[0]))
    return dfs([],[(1<<N)-1]*6)


for i, fn in enumerate(sorted(PATH.glob('witness*.b64'))):
    code=fn.read_text().strip()
    assert len(code)%6==0
    rows=[tuple(ABC.index(v) for v in code[i:i+6]) for i in range(0,len(code),6)]
    assert all(0<=x<N for r in rows for x in r)
    t=read_tables(PATH/f'b{i}.txt')
    checked(rows,t)
    a=residual(t,rows)
    left={sum(x.bit_count() for x in b) for b in a.values()}
    assert left=={(2754 if i%2==0 else 2755)-len(rows)}
    # Every stored packing is admissible, not necessarily globally maximum.
    print(i, len(rows), "verified six-way tuples",
          "uncovered pairs per board",left.pop(), "direct addable rows",remaining_cliques(a))
