"""Independent audit of decoder_pair_repair section 5: minimum reversals to reach
H5[TT2] (pair modules + regular 5-vertex quotient), per core class.
Claimed (transitive, source/sink singleton inside, strong+transitive triple,
strong+cyclic triple, source/sink singleton at apex) -> (15, 13, 12, 7, 10)."""
from itertools import combinations, permutations
block = [0,0,0,1,1,1,2,2,2,3]
cyc = {(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(6,7),(7,8),(8,6)}
core_pairs = list(combinations(range(4), 2))
def matchings(rem):
    if not rem:
        yield []; return
    a = rem[0]
    for i in range(1, len(rem)):
        rest = rem[1:i] + rem[i+1:]
        for m in matchings(rest):
            yield [(rem[0], rem[i])] + m
Ms = list(matchings(list(range(10))))
qpairs = list(combinations(range(5), 2))
regs = []
for mask in range(1 << 10):
    out = [0]*5
    for i, (a, b) in enumerate(qpairs):
        if (mask >> i) & 1: out[a] += 1
        else: out[b] += 1
    if all(o == 2 for o in out): regs.append(mask)
assert len(regs) == 24
res = {}
for mask in range(64):
    core = {}
    for i, (a, b) in enumerate(core_pairs):
        core[(a,b)] = 1 if (mask >> i) & 1 else -1
        core[(b,a)] = -core[(a,b)]
    def arc(u, v):
        bu, bv = block[u], block[v]
        if bu == bv: return 1 if (u, v) in cyc else -1
        return core[(bu, bv)]
    # core class descriptors
    outd = [sum(1 for j in range(4) if j != i and core[(i,j)] == 1) for i in range(4)]
    scores = tuple(sorted(outd))
    root_out = outd[3]
    tri = [0,1,2]
    tri_cyc = (core[(0,1)] == core[(1,2)] == core[(2,0)])
    best_unres = best_can = None
    for M in Ms:
        k = {}
        for (p, r) in qpairs:
            k[(p,r)] = sum(1 for u in M[p] for v in M[r] if arc(u, v) == 1)
        un = sum(min(x, 4 - x) for x in k.values())
        if best_unres is None or un < best_unres: best_unres = un
        for rm in regs:
            c = 0
            for i, (p, r) in enumerate(qpairs):
                c += (4 - k[(p,r)]) if (rm >> i) & 1 else k[(p,r)]
            if best_can is None or c < best_can: best_can = c
    if scores == (0,1,2,3): cls = "transitive"
    elif scores == (1,1,2,2): cls = "strong, deleted-singleton triple " + ("cyclic" if tri_cyc else "transitive")
    else: cls = "source/sink over C3, singleton " + ("at apex" if root_out in (0,3) else "inside triangle")
    res.setdefault(cls, set()).add((best_unres, best_can))
for k, v in sorted(res.items()): print(k, "-> (unrestricted, canonical H5[TT2]) minima:", sorted(v))
