"""Emit one representative of each both-strong missing-arc quotient of a small order."""
import argparse
from itertools import combinations, permutations
from seed_lifts import edit, strong


def bank(n):
    free_edges=list(combinations(range(n),2))[1:]
    orderings=[head+tail for head in ((0,1),(1,0))
               for tail in permutations(range(2,n))]
    seen=set()
    for code in range(1<<len(free_edges)):
        if code in seen:
            continue
        rows=[0]*n
        rows[0]|=2
        for k,(i,j) in enumerate(free_edges):
            if code>>k&1:
                rows[i]|=1<<j
            else:
                rows[j]|=1<<i
        for p in orderings:
            seen.add(sum((rows[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(free_edges)))
        if strong(rows) and strong(edit(rows,0,1)):
            yield rows


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('order',type=int)
    args=ap.parse_args()
    for r in bank(args.order): print(len(r),*r)
