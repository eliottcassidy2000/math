"""Necessary intersection test for the seven-branch partitions, and actual small charts."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from itertools import combinations_with_replacement
from math import ceil
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/"moore57_thirteen_20261125"))
from moore_chart import hoffman_edge


def partitions(n, lower=6):
    if n == 0:
        yield ()
    for first in range(lower, n+1):
        for tail in partitions(n-first, first):
            yield (first,)+tail


def compatible_cuts(sizes):
    n = len(sizes)
    all_mask = (1<<n)-1
    widths = [0]*(all_mask+1)
    for m in range(1,all_mask+1):
        lone = m&-m
        widths[m] = widths[m^lone] + sizes[lone.bit_length()-1]
    cuts = [m for m in range(1,all_mask) if widths[m] in (14,21,28)
            and (widths[m] < 28 or m < (m^all_mask))]
    rows = [0]*len(cuts)
    for i,j in combinations(range(len(cuts)),2):
        good = True
        for a in (cuts[i],cuts[i]^all_mask):
            for b in (cuts[j],cuts[j]^all_mask):
                x, y = widths[a]//7, widths[b]//7
                left = x*y-widths[a&b]
                if not 0 <= left <= min(x,y):
                    good = False
        if good:
            rows[i] |= 1<<j
            rows[j] |= 1<<i
    return cuts, rows


def max_cliques(rows, all_at_target=None):
    best = 0
    saved = []

    def extend(mask, selected):
        nonlocal best
        if all_at_target is not None and len(selected) == all_at_target:
            saved.append(selected)
            return
        best = max(best,len(selected))
        required = all_at_target if all_at_target is not None else best+1
        while mask and len(selected)+mask.bit_count() >= required:
            v = mask&-mask
            mask ^= v
            i = v.bit_length()-1
            extend(mask & rows[i], selected+(i,))
            required = all_at_target if all_at_target is not None else best+1
    extend((1<<len(rows))-1, ())
    return best,saved


def extremal_cuts():
    hist = Counter()
    for widths in partitions(56):
        masks, rows = compatible_cuts(widths)
        cap = max_cliques(rows)[0]
        hist[cap] += 1
        if cap == 7:
            assert widths == (7,)*8
            geometries = Counter()
            for indices in max_cliques(rows,7)[1]:
                masks_used = [masks[j] for j in indices]
                k = masks_used[0].bit_count()
                assert all(m.bit_count()==k for m in masks_used)
                if k == 3:
                    used = 0
                    for mask in masks_used:
                        used |= mask
                    assert used.bit_count()==7
                    assert all((a&b).bit_count()==1 for a,b in combinations(masks_used,2))
                else:
                    assert k == 4
                    # The seven balanced bipartitions, all oriented away from the
                    # eighth bit, form the nonzero binary simplex code of dimension 3.
                    assert all(a^b in masks_used for a,b in combinations(masks_used,2))
                geometries[k] += 1
            assert geometries == {3:240,4:30}
    assert hist == {0:407,1:675,2:218,3:41,4:8,5:2,6:1,7:1}
    return dict(sorted(hist.items()))


def bound(m, s):
    frac = Fraction(m*s,56)%1
    initial = Fraction((m-8)*(m+7)*m*s*(56-s),56)
    noise = (57-m)*56*frac*(1-frac)
    return ceil((initial+noise)/(4*m-6))


def simultaneous_contact_profiles():
    # Four off-line Fano atoms, allowing any nonnegative counts with total <= 7.
    fano = [counts for counts in combinations_with_replacement(range(8),4)
            if sum(counts) <= 7 and all(a+b in (2,3) for a,b in combinations(counts,2))]
    assert fano == [(1,1,1,1),(1,1,1,2)]
    assert 35*max(map(sum,fano)) < 4*7*7

    # Seven core neighbors counted on eight point labels. Branch 1 is used
    # for the exceptional control by symmetry of nonzero characters.
    accepted = Counter()
    for labels in combinations_with_replacement(range(8),7):
        coeff = [sum((-1)**((x&j).bit_count()%2) for x in labels) for j in range(1,8)]
        ordinary = all(abs(v)==1 for v in coeff)
        exceptional = abs(coeff[0]) == 7 and all(abs(v)==1 for v in coeff[1:])
        if ordinary:
            assert len(set(labels)) == 7
            accepted["ordinary"] += 1
        if exceptional:
            assert sorted(Counter(labels).values()) == [1,2,2,2]
            accepted["fixed exceptional branch"] += 1
    assert accepted == {"ordinary":8,"fixed exceptional branch":8}
    print("Surviving affine neighbor multisets per type",dict(accepted))


def check_known():
    p = list(combinations(range(5),2))
    A = {v:{w for w in p if set(v).isdisjoint(w)} for v in p}
    H = {v:{w for w in range(50) if hoffman_edge(v,w)} for v in range(50)}
    for graph in (A,H):
        root = next(iter(graph))
        leaves = [graph[v]-{root} for v in graph[root]]
        ans = {}
        bipartite = Counter()
        for m in range(2,len(leaves)+1):
            shapes = set()
            for fibers in combinations(leaves,m):
                selected = set.union(*fibers)
                counts = []
                two_colorable = True
                while selected:
                    chunk = set()
                    fringe = {next(iter(selected))}
                    while fringe:
                        if any(graph[x] & fringe for x in fringe):
                            two_colorable = False
                        chunk.update(fringe)
                        fringe = set.union(*(graph[x] for x in fringe)) & (selected-chunk)
                    selected -= chunk
                    assert not len(chunk)%m
                    counts.append(len(chunk)//m)
                shapes.add(tuple(sorted(counts)))
                if two_colorable:
                    bipartite[m] += 1
            ans[m] = shapes
        print("Degree",len(leaves),"component sizes per branch",ans)
        print("Bipartite subsets by branch count",dict(bipartite))
        if len(graph)==50:
            assert ans[3] == {(2,2,2)} and ans[4] == {(6,)}
            assert bipartite == {2:21, 3:35, 4:35}
        else:
            assert ans[3] == {(2,)}


if __name__=="__main__":
    print("Seven-branch partitions by maximum allowed exceptional extensions",extremal_cuts())
    simultaneous_contact_profiles()
    assert [bound(8,r) for r in range(7)] == [0,13,22,26,26,22,13]
    assert (bound(8,13),bound(9,13)) == (13,56)
    print("Eight/nine branches, thirteen tracks:", bound(8,13), bound(9,13))
    assert 52*49 + Fraction(49**2,18) < 2688
    assert 52*50 + 10**2 - 2*10 < 2688
    check_known()
