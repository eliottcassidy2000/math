"""Collapse transitive twin blocks and canonically identify the actual lower hosts."""
from collections import Counter
from itertools import combinations, permutations, product
import math
import sys
from seam_widths import solve

E = 2**60*(2**61-1)
FOUR_SOURCE = "110111111111111111111111101111111111101111101"
FOUR_CASES = (
    ((0,4,6,7), (29,29,30,30)),
    ((0,4,6,9), (14,14,30,61)),
    ((0,4,6,9), (14,14,61,30)),
)


def canonical(rows, sizes):
    rows=[list(row) for row in rows]
    sizes=list(sizes)
    while True:
        twins=next(((i,j) for i,j in combinations(range(len(rows)),2)
                    if all(rows[i][v]==rows[j][v] for v in range(len(rows)) if v not in (i,j))),None)
        if twins is None: break
        i,j=twins
        sizes[i]+=sizes.pop(j)
        rows.pop(j)
        for row in rows: row.pop(j)
    ordered=sorted(set(sizes))
    blocks=[[v for v,w in enumerate(sizes) if w==s] for s in ordered]
    def word(perm):
        return tuple(rows[i][j] for i,j in combinations(perm,2))
    bits=min(word(sum(parts,())) for parts in product(*(permutations(b) for b in blocks)))
    return tuple(sorted(sizes)),bits


def decode(line):
    n=(1+math.isqrt(1+8*len(line)))//2
    rows=[[0]*n for i in range(n)]
    for bit,(i,j) in zip(line,combinations(range(n),2)):
        rows[i][j]=int(bit)
        rows[j][i]=1-int(bit)
    return rows


def main():
    sample=(2,8,3,37,15,151,63,95)
    reference=[[mask>>i&1 for i in range(8)] for mask in sample]
    sizes=(1,29,1,1,29,1,61,1)
    names={canonical(reference,sizes):"prior",canonical(list(zip(*reference)),sizes):"prior converse"}
    board=decode(FOUR_SOURCE)
    for label,(sites,widths) in zip(("half","quarter A","quarter B"),FOUR_CASES):
        sizes=[1]*len(board)
        for site,width in zip(sites,widths): sizes[site]=width
        names[canonical(board,sizes)]=label
        names[canonical(list(zip(*board)),sizes)]=label+" converse"
    assert len(names)==8
    counts=Counter()
    for record in sys.stdin:
        line,*fields=record.split()
        k=next(i for i in range(1,10) if i+2**i==len(fields))
        ports=tuple(map(int,fields[:k]))
        f=list(map(int,fields[k:]))
        rows=decode(line)
        for widths in solve(f,E-1):
            sizes=[1]*len(rows)
            for i,w in zip(ports,widths): sizes[i]=w
            key=canonical(rows,sizes)
            label=names.get(key,str(key))
            counts[len(rows),label]+=1
    print(counts)


if __name__=="__main__":
    main()
