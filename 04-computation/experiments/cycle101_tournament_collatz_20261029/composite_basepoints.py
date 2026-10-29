"""HYP-9446 (4): tournament base points per necklace for composite odd q (the two-base-point lemma of THM-4849 (B)
uses that 2t generates Z_q).  For q in 9, 15, 21, 25, 27: number of connection sets, tournament cycles, the
distribution of base points per cycle, and whether every cycle with >= 2 base points is m*{1..k} for a unit m."""
import itertools
from math import gcd
def canon(S, q):
    return min(tuple(sorted((x - t) % q for x in S)) for t in range(q) if t in S or True)
for q in (9, 15, 21, 25, 27):
    k = (q - 1) // 2
    neck = {}
    for choice in itertools.product((0, 1), repeat=k):
        S = frozenset(i if c == 0 else q - i for i, c in zip(range(1, k + 1), choice))
        neck.setdefault(canon(S, q), []).append(S)
    I = frozenset(range(1, k + 1))
    Rcls = {canon(frozenset((m * x) % q for x in I), q) for m in range(1, q) if gcd(m, q) == 1}
    dist = {}
    bad = []
    for c, sets in neck.items():
        dist[len(sets)] = dist.get(len(sets), 0) + 1
        if len(sets) >= 2 and c not in Rcls:
            bad.append((c, len(sets)))
    print('q=%d: %d connection sets, %d tournament cycles, base points per cycle %s, cycles with >= 2 base points not of the form m*I: %d %s'
          % (q, 2 ** k, len(neck), dict(sorted(dist.items())), len(bad), bad[:3]))

# structure check (THM-4849 (B')): if S and S + t are both connection sets, g = gcd(t, q), then every coset of <t>
# other than <t> lies inside S or is disjoint from it, and S cap <t> = -2t * {1, ..., (q/g - 1)/2}:
# T(S) is the lexicographic product C[R_(q/g)] of a circulant tournament C on Z_q/<t> with R_(q/g).
print()
for q in (9, 15, 21, 25, 27, 7, 11, 13):
    k = (q - 1) // 2
    sets = [frozenset(i if c == 0 else q - i for i, c in zip(range(1, k + 1), ch)) for ch in itertools.product((0, 1), repeat=k)]
    setset = set(sets)
    pairs = ok = 0
    for S in sets:
        for t in range(1, q):
            St = frozenset((x + t) % q for x in S)
            if St in setset:
                pairs += 1
                g = gcd(t, q)
                sub = set(range(0, q, g))
                good = all(len(S & {(c + y) % q for y in sub}) in (0, len(sub)) for c in range(1, g))
                half = (q // g - 1) // 2
                good &= (S & sub) == {(-2 * t * j) % q for j in range(1, half + 1)}
                ok += good
    print('q=%d: ordered pairs (S, S+t) of connection sets: %d, with the lexicographic-product structure: %d' % (q, pairs, ok))
