# When the switching shadow stops identifying distinct frame moves

In the switching observer of README, let P range over undirected
Hamilton paths on n labelled vertices. Two frames have the same move
exactly when their symmetric difference is a cut. We obtain:

| n | distinct moves | number of frames per move |
|---:|---:|---:|
| 3 | 1 | 3 |
| 4 | 6 | 2 |
| 5 | 15 | 4 |
| 6 | 180 | 2 |
| 7 | 1260 | 2 |
| ≥8 | n!/2 | 1 |

The last row is an all-order recovery theorem, not a census. P△Q has
maximum degree four and at most 2(n−1) edges. Every nonzero cut on
n≥9 has some degree at least five (take the smaller shore). At n=8
the only nonzero cut with degree at most four is K4,4, which has
sixteen edges instead of at most fourteen. So the projection from
frames to moves is injective from n=8 onwards.

At n=7 a nonzero cut can only be K3,4 and it has exactly twelve
edges. Thus P,Q must be disjoint and exhaust this cut. The only
possible cut is the path's alternating bipartition. For the path
0123456 its mate is 2503614, and all other frames are relabellings.
At n=6 a cut of maximum degree four and at most ten edges can
only have two vertices on its smaller shore. Both must be
nonadjacent interior vertices of P. For P=012345 the possible
pairs of positions are {1,3},{1,4},{2,4}; only {1,4} gives a
Hamilton path after taking the symmetric difference (042315).
The small remaining rows are checked by the independent short
enumeration below. In particular the five-label Petersen shadow's
four-to-one phenomenon is an exceptional small-order compression,
not a recursive identification persisting in large size.

On odd n a switching move's Eulerian representative is K_n△P△δ(E),
where E is the pair of frame endpoints. For n≥5 its degrees are
2 at the endpoints, n−3 at their two neighbours along the path,
and n−5 at the remaining vertices. Thus for odd n≥9 one can
locate the frame endpoints from the shadow and recover its full
edge set immediately.

For n=1 mod 4 (n≥5) even the converse quotient is a genuine
double cover of a simple graph, not necessarily complete. Write
z=K_n. If two distinct moves differed by z, their Eulerian
representatives would be complements. They have the same weight
binom(n,2)−3n+9, which equals half the complete graph's weight
only at n=9 in this range. At n=9, complementing changes the
two degree-6 vertices (the two endpoint neighbours in P) into
the new degree-2 endpoints, so if a frame Q existed their combined
endpoint set E△E' would have four vertices. The interior vertex
two steps along P is adjacent in P to one of those four. Its
degree in Q=K_n△P△δ(E△E') would be 6−4+2=4 (or higher if
it touches two), impossible. At n=5 the no-collision result
also follows from the different weights or README's quadratic bit.
No move equals z because P is not itself a cut: a spanning path
on n≥5 has degree at most two and every nonzero cut has degree ≥3.
Hence the projection is locally bijective, with no loops introduced.
Only at n=5 is the converse-quotient base K16; n≥9 has far fewer
moves than its number of other vertices.

Replay for the initial rows:

    from itertools import combinations, permutations
    from collections import Counter
    for n, expected in [(3,(1,3)),(4,(6,2)),(5,(15,4)),(6,(180,2)),(7,(1260,2)),(8,(20160,1))]:
        es=list(combinations(range(n),2))
        ids={e:j for j,e in enumerate(es)}
        cuts=[sum((((x>>i)^(x>>j))&1)<<k for k,(i,j) in enumerate(es))
              for x in range(1<<(n-1))]
        C=Counter()
        for p in permutations(range(n)):
            if p[0]>p[-1]: continue
            a=sum(1<<ids[tuple(sorted(e))] for e in zip(p,p[1:]))
            C[min(a^c for c in cuts)]+=1
        assert (len(C),next(iter(C.values())))==expected
        assert len(set(C.values()))==1
