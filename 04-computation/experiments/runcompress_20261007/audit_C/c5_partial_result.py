#!/usr/bin/env python3
"""audit_C item 5: HYP-9242 'Partial result: the upper bound for random t' - checks of the ingredients.

(a) exact uniformity: for K in {2, 3, 5}, L = 16, enumerate ALL odd L-bit t in the residual class; the child h_1's parity
    bits from y_1 = T^(K-2)(h_1) on (the D = 1 chain's driving bits) of length L are: a fixed 3-bit prefix, then all
    2^(L-3) patterns exactly once; and conditionally on (J, r) the bits after the two-run and the forced post-run prefix
    are again exactly uniform.
(b) merge value: whenever the D = 1 chain (state (1, 2) at x) is absorbed at chain time s <= L - K - 2, the common value
    is >= 2^(K+2) >= 2 (so it is a deletion certificate), and the (odd count, Terras difference) criterion holds for D = 1.
(c) the implication 'absorbed by L - K - 2  =>  not an orphan' holds sample-wise, and P(D = 1 chain not absorbed by
    L - K - 2) is printed with sqrt(L) x P (the claimed upper-bound rate L^-1/2 up to logs), next to P(orphan).
(d) the post-run state of the D = 1 chain is (J, 1) (THM-4601 (v)) and J is geometric (P(J >= m) = 4^(1-m)).
"""
import random, math
from collections import Counter, defaultdict
exec(open('c4_general_sources.py').read().split('N = int(sys.argv[1])')[0])   # ot(), v2()


def T(x):
    return (3 * x + 1) >> 1 if x & 1 else x >> 1


def bits(y, s):
    out = []
    for _ in range(s):
        out.append(y & 1)
        y = T(y)
    return tuple(out)


# (a)
L = 16
for K in (2, 3, 5):
    ts = [t for t in range(1 << (L - 1), 1 << L) if t & 1 and v2(3 ** K * t - 1) == 1]
    vecs = Counter()
    byJ = defaultdict(Counter)
    for t in ts:
        y1 = 2 * 3 ** (K - 2) * t - 1 if K >= 2 else None
        b = bits(y1, L)
        vecs[b] += 1
        x = 2 * 3 ** (K - 1) * t - 1
        j = v2(x - 1)
        J, r = (j - 1) // 2, j - 2 * ((j - 1) // 2)
        pre = 2 * J + (2 if r == 1 else 3)
        if pre + 4 <= L:
            byJ[(J, r)][b[pre:]] += 1
    prefixes = Counter(b[:3] for b in vecs)
    assert len(prefixes) == 1 and all(c == 1 for c in vecs.values()) and len(vecs) == 2 ** (L - 3)
    unif = all(len(c) == 2 ** (L - pre_) and len(set(c.values())) == 1
               for (J, r), c in byJ.items()
               for pre_ in [2 * J + (2 if r == 1 else 3)])
    print(f"(a) K={K}: {len(ts)} t; driving bits = fixed prefix {list(prefixes)[0]} + all {2**(L-3)} patterns once; "
          f"conditionally uniform after the two-run and post-run prefix for every (J, r) with room: {unif}")

# (b), (c), (d)
rnd = random.Random(5)
for K in (2, 5, 9):
    for Lb in (100, 200, 400, 800, 1600):
        N = 600 if Lb <= 400 else 300
        fail = orph = 0
        Jc = Counter()
        tot = 0
        while tot < N:
            t = rnd.getrandbits(Lb) | (1 << (Lb - 1)) | 1
            if v2(3 ** K * t - 1) != 1:
                continue
            tot += 1
            n = (t << K) - 1
            x = 2 * 3 ** (K - 1) * t - 1
            y = 2 * 3 ** (K - 2) * t - 1 if K >= 2 else None
            j = v2(x - 1)
            J = (j - 1) // 2
            Jc[J] += 1
            u, v, k = x, y, 1
            Tmax = Lb - K - 2
            absorbed = None
            for s in range(1, Tmax + 1):
                k += (u & 1) - (v & 1)
                u, v = T(u), T(v)
                if s == 2 * J:
                    assert k == J and u - 3 ** J * v == 1, (K, Lb, J)      # post-run state (J, 1)
                if u == v and k == 0:
                    absorbed = (s, u)
                    break
            on, tn = ot(n)
            oh, th = ot(((n + 1) >> 1) - 1)
            d1 = (on == oh and tn - th == 1)
            if absorbed:
                assert absorbed[1] >= 2 ** (K + 2) and d1
            else:
                fail += 1
            # orphan?
            found = d1
            D = 2
            while not found and D < K:
                oh, th = ot(((n + 1) >> D) - 1)
                found = (on == oh and tn - th == D)
                D += 1
            orph += not found
            assert (not absorbed) or found
        pf, po = fail / N, orph / N
        print(f"K={K} L={Lb:5d}: P(D=1 chain not absorbed by L-K-2) = {pf:.3f} (x sqrt(L) = {pf*math.sqrt(Lb):.2f}); "
              f"P(orphan) = {po:.3f}; J >= 3 share {sum(c for J, c in Jc.items() if J >= 3)/N:.3f} (4^-2 = 0.0625)")
print("ALL CHECKS PASSED")
