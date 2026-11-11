"""a10_verify_witness.py (AUDIT) - standalone check of the counterexample to HYP-9664 (C) at n = 5.
Reads a10_witness_n5.txt (465-bit word; pairs (u < v) of points 1..31 in lexicographic order, bit i = 1 iff u -> v)
and counts, with no linear algebra: cyclic lines, and for each of the 155 planes its even zero-sum 4-sets
(a 4-set is even iff it contains 0 or 2 cyclic triangles).  Also the same counts for the Cayley-Dickson sign
tournament CD_5 (own doubling) for comparison."""
import itertools, os
HERE = os.path.dirname(os.path.abspath(__file__))
PTS = list(range(1, 32))
PAIRS = [(u, v) for u in PTS for v in PTS if u < v]
t = int(open(os.path.join(HERE, "a10_witness_n5.txt")).read())
beats = {}
for i, (u, v) in enumerate(PAIRS):
    b = (t >> i) & 1
    beats[(u, v)] = bool(b)
    beats[(v, u)] = not b
def conj(x):
    if len(x) == 1: return x[:]
    h = len(x) // 2
    return conj(x[:h]) + [-z for z in x[h:]]
def mul(x, y):
    if len(x) == 1: return [x[0] * y[0]]
    h = len(x) // 2
    a, b, c, d = x[:h], x[h:], y[:h], y[h:]
    return [p - q for p, q in zip(mul(a, c), mul(conj(d), b))] + [p + q for p, q in zip(mul(d, a), mul(b, conj(c)))]
cd = {}
for i in PTS:
    for j in PTS:
        if i != j:
            x = [0] * 32; y = [0] * 32; x[i] = 1; y[j] = 1
            cd[(i, j)] = mul(x, y)[i ^ j] == 1
planes = sorted({tuple(sorted({a, b, c, a ^ b, a ^ c, b ^ c, a ^ b ^ c})) for a in PTS for b in PTS for c in PTS if a < b < c and c != a ^ b})
lines = sorted({tuple(sorted((a, b, a ^ b))) for a in PTS for b in PTS if a < b})
def stats(arc):
    def cyc(u, v, w):
        return (arc[(u, v)] and arc[(v, w)] and arc[(w, u)]) or (arc[(v, u)] and arc[(w, v)] and arc[(u, w)])
    cl = sum(1 for (a, b, c) in lines if cyc(a, b, c))
    odd = ev = 0
    for pl in planes:
        bad = 0
        for q in itertools.combinations(pl, 4):
            if q[0] ^ q[1] ^ q[2] ^ q[3]: continue
            k = sum(1 for tr in itertools.combinations(q, 3) if cyc(*tr))
            assert k <= 2
            bad += (k != 1)
        ev += bad; odd += (bad == 0)
    sc = sorted({sum(1 for v in PTS if v != u and arc[(u, v)]) for u in PTS})
    return cl, odd, ev, sc
res = []
for name, arc in (("audit witness (planes skew to the line {1,2,3} odd)", beats), ("Cayley-Dickson sign tournament CD_5", cd)):
    cl, odd, ev, sc = stats(arc)
    res.append("%s: cyclic lines %d/155, odd planes %d/155, even zero-sum 4-sets %d, scores %s" % (name, cl, odd, ev, sc))
print("\n".join(res))
open(os.path.join(HERE, "a10_verify_witness.out"), "w").write("\n".join(res) + "\n")
