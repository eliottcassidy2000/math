"""Independent check of physical candidates, projected completions, and two-branch hinges."""
from collections import Counter
from itertools import combinations
from pathlib import Path
import sys

from nested_core_check import graph, legal, orbit_chart
from three_stars import SLOTS

BRANCHES = tuple(h for h in SLOTS if h != 1)
CENTERS = (56, 0, 28)
SOURCE = 58 + 56 + 1


def read_board(path):
    items = iter(map(int, Path(path).read_text().split()))
    assert next(items) == len(BRANCHES) == 53
    ref, rows = [], []
    for _ in BRANCHES:
        ref.append(next(items))
        rows.append({})
        for _ in range(next(items)):
            y, *slot = (next(items) for _ in range(4))
            assert y not in rows[-1]
            rows[-1][y] = tuple(slot)
    assert next(items, None) is None
    return ref, rows


def leaf(h, y):
    return 58 + 56*h + y


def certify_board(rows):
    chart, stars = orbit_chart('period_fourteen_site.txt', 4, ((10, 45, 2), (27, 49, 0)))
    A = graph(chart, stars)
    assert sum(map(len, A)) == 2*13231
    for R in A:
        # The original graph is independently verified to be triangle- and square-free.
        assert len(R) <= 57
    star = {v for v,R in enumerate(A) if len(R)==57}
    assert len(star)==240 and len(A[SOURCE])==4
    for v,R in enumerate(A):
        reached = {v} | R
        for w in R:
            D=A[w]-{v}
            assert not (D & reached)
            reached.update(D)
    marked = leaf(BRANCHES[0], 7)
    assert BRANCHES[0] == 2 and legal(A, SOURCE, marked)
    existing_slots = [set(A[SOURCE]) | {SOURCE}]
    existing_slots.append(set().union(*(A[v] for v in existing_slots[0])) - {SOURCE})
    allowed = tuple(set(range(56)) - {y for y in range(56)
                                       if any(leaf(c,y) in R for R in existing_slots)}
                    for c in CENTERS)
    assert tuple(map(len, allowed)) == (53, 53, 53)
    for r,h in enumerate(BRANCHES):
        expected = {y:tuple(chart[h,c][y] for c in CENTERS)
                    for y in range(56)
                    if legal(A, SOURCE, leaf(h,y)) and (r!=0 or y==7)}
        assert rows[r] == expected, (r, h)
    A[SOURCE].add(marked); A[marked].add(SOURCE)
    R = [{y: t for y,t in dom.items() if legal(A, SOURCE, leaf(h,y))}
         if i else dom for i,(h,dom) in enumerate(zip(BRANCHES,rows))]
    m = rows[0][7]
    assert R[0] == {7:m}
    assert all(R[i] == {y:t for y,t in dom.items() if all(a!=b for a,b in zip(t,m))}
               for i,dom in enumerate(rows) if i)
    assert all({t[k] for D in R[1:] for t in D.values()} == allowed[k]-{m[k]}
               for k in range(3))
    V = [(i,y,t,leaf(BRANCHES[i],y)) for i,dom in enumerate(R) if i
         for y,t in dom.items()]
    assert len(V)==1741
    N = [sum(1<<u for u in X) for X in A]
    for (i,y,t,u),(j,z,s,v) in combinations(V,2):
        actual = v in A[u] or bool(N[u]&N[v])
        encoded = i==j or any(a==b for a,b in zip(t,s))
        assert actual == encoded, (i,y,j,z)
    print(f"Rebuilt {len(V)} contacts and checked {len(V)*(len(V)-1)//2} physical pairs.")
    return R


def projected_covers(rows, path):
    C = [(h,y) for h in range(1,53) for y in rows[h]]
    seen = set()
    for line in Path(path).read_text().splitlines():
        g,k,*ys = map(int,line.split())
        assert 0<=g<3 and 0<=k<len(C) and (g,k) not in seen
        seen.add((g,k))
        h,y=C[k]
        assert len(ys)==53 and ys[h]==y and ys[0]==7
        picks=[dom[y] for dom,y in zip(rows,ys)]
        for z in range(3):
            if z!=g:
                assert len({P[z] for P in picks})==53
    assert len(seen)==3*len(C)==5223
    print("Verified 5,223 completions: every physical contact in every two-center projection.")


def critical_rectangles(rows):
    # The common first coordinate is the physical sheet. A 2x2 rectangle
    # sensitive to deleting each of the 3 margins must have two equal-sheet
    # diagonal pairs. The remaining two conflicts then use different margins.
    assert all(t[0]==y for D in rows for y,t in D.items())
    count=Counter()
    example=None
    for h,k in combinations(range(1,53),2):
        for a,b in combinations(sorted(rows[h].keys()&rows[k].keys()),2):
            bits = [sum(1<<g for g,(v,w) in enumerate(zip(rows[h][y],rows[k][z])) if v==w)
                    for y in (a,b) for z in (a,b)]
            if all(bits) and {1,2,4} <= set(bits):
                count[BRANCHES[h],BRANCHES[k]]+=1
                example=example or (BRANCHES[h],BRANCHES[k],a,b,bits)
    assert example==(3,12,11,28,[1,4,2,1])
    print("Critical two-branch rectangles:",sum(count.values()),"on",len(count),"branch pairs;",
          "first",example,"maximum per branch pair",max(count.values()))
    return count


if __name__=='__main__':
    assert len(sys.argv)>1
    _, R=read_board(sys.argv[1])
    R=certify_board(R)
    critical_rectangles(R)
    if len(sys.argv)>2:
        projected_covers(R,sys.argv[2])
