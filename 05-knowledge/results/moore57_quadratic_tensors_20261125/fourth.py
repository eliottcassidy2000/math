"""Joint voltage constraint and a literal extension of the 28-trade seed."""
from collections import Counter
from itertools import combinations
from three_stars import A, B, IDENTITY, N, PIVOT, SLOTS, TRADES, make_chart, routes

# h in increasing order through SLOTS except the source branch 1.
PICKS = (53,27,14,42,13,2,29,37,17,20,16,26,45,7,12,13,36,32,
         38,35,44,52,31,8,51,11,19,15,54,5,50,9,10,25,48,3,
         22,47,24,18,55,41,49,6,33,39,1,34,40,30,47,28,46)


def selected(source=1, slot=0):
    maps = make_chart()
    branches = [h for h in SLOTS if h != source]
    picks = dict(zip(branches, PICKS))
    centers = (PIVOT, A, B)
    blocked = [{row[slot] for row in routes(maps, source, c)} for c in centers]
    triples = {h: tuple(maps[h, c][y] for c in centers) for h, y in picks.items()}
    assert all(all(value not in bad for value, bad in zip(t, blocked)) for t in triples.values())
    clashes = [(i,j) for i,j in combinations(branches,2)
               if any(a == b for a,b in zip(triples[i],triples[j]))]
    assert clashes == [(6,17), (40,53)]
    for h in (17,53):
        del picks[h]
    for t in range(3):
        assert len({triples[h][t] for h in picks}) == len(picks)
    assert len(picks) == 51
    print("source (1,0), 51 new edges, unfilled branches 17 and 53")


def minimum_changes():
    # Sets of residues are rows of 56-bit integers; rotation implements addition.
    maps = make_chart()
    base = make_chart(())
    mask = (1 << N)-1
    def rot(bits, k):
        return ((bits << k) | (bits >> (N-k))) & mask

    profile = {}
    for source in SLOTS:
        costs = [[0]*N for _ in range(8)]
        costs[0][0] = 1
        for h in SLOTS:
            if h == source:
                continue
            steps = {((y-maps[h,A][y])%N, (y-maps[h,B][y])%N) for y in IDENTITY}
            original = ((-base[h,A][0])%N, (-base[h,B][0])%N)
            nxt = [[0]*N for _ in range(8)]
            for a,b in steps:
                price = int((a,b) != original)
                for cost in range(8-price):
                    for i,bits in enumerate(costs[cost]):
                        if bits:
                            nxt[cost+price][(i+a)%N] |= rot(bits,b)
            costs = nxt
        profile[source] = next(k for k, block in enumerate(costs) if block[28] & (1<<28))
    print("minimum nonoriginal rows",dict(sorted(Counter(profile.values()).items())))
    assert Counter(profile.values()) == {3:1,4:4,5:45,6:4}
    assert [h for h,k in profile.items() if k == 6] == [2,12,35,42]
    assert [h for h,k in profile.items() if k == 3] == [38]


if __name__ == "__main__":
    assert len(TRADES) == 28
    import sys
    if sys.argv[1:] == ["--board"]:
        maps = make_chart()
        centers = (PIVOT, A, B)
        forbidden = [{r[0] for r in routes(maps, 1, c)} for c in centers]
        for h, selected_y in zip((v for v in SLOTS if v != 1), PICKS):
            for y in IDENTITY:
                triple = tuple(maps[h, c][y] for c in centers)
                if all(z not in bad for z, bad in zip(triple, forbidden)):
                    print(h, *triple, int(y == selected_y))
    else:
        selected()
        minimum_changes()
