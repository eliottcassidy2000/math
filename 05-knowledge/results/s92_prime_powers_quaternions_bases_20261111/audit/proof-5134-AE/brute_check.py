# brute_check.py -- auditor proof-5134-AE.
# Slow, pure-Python, definition-level census of F1 for a small X, used to validate census_f1.py.
# Enumerates ORDERED triples (h1, h2, s) directly from the definitions, no symmetry, no pruning tricks
# other than the trivial bound 2^s y1 y2 <= M.
import sys, math, json
from functools import lru_cache


@lru_cache(maxsize=None)
def psi(h):
    if h == 0:
        return 0
    k = h // 2
    if h % 2 == 0:
        return psi(k) + k
    return psi(k) + k // 2 + 1


def main():
    L = int(sys.argv[1])
    mods = [int(t) for t in sys.argv[2].split(",")]
    X = 1 << L
    # q recursion check: q(2h) = q(h), q(2h+1) = q(h) + ceil(h/2)
    for h in range(0, 5000):
        q = h - psi(h)
        assert (2 * h) - psi(2 * h) == q
        assert (2 * h + 1) - psi(2 * h + 1) == q + (h + 1) // 2
        assert 0 <= q <= h / 2
    out = {}
    cnt = {m: [[dict() for _ in range(2)] for _ in range(2)] for m in mods}  # cnt[m][sheet][diag?][(d, r)]
    tot = [0, 0]
    hmax = X + 2
    for s in range(0, L + 2):
        sh = 0 if s == 0 else 1
        Xp = X + 1 if s == 0 else X     # N = M - [s = 0] <= X
        for h1 in range(0, hmax):
            x1, y1 = 2 * psi(h1) + 1, 2 * h1 + 1
            if (1 << s) * y1 > Xp:
                break
            for h2 in range(0, hmax):
                x2, y2 = 2 * psi(h2) + 1, 2 * h2 + 1
                if (1 << s) * y1 * y2 > Xp:
                    break
                M = x1 * x2 + (1 << s) * y1 * y2
                if M <= Xp:
                    tot[sh] += 1
                    hs = min(h1, h2)
                    g = math.gcd(2 * psi(hs) + 1, 2 * hs + 1)
                    for m in mods:
                        key = (math.gcd(g, m), M % m)
                        dd = cnt[m][sh][1 if h1 == h2 else 0]
                        dd[key] = dd.get(key, 0) + 1
    print("X = 2^%d ordered triples: sheet0 %d sheet1 %d" % (L, tot[0], tot[1]))
    res = {"log2X": L, "tot": tot, "data": {}}
    for m in mods:
        res["data"][str(m)] = [[sorted([[k[0], k[1], v] for k, v in cnt[m][sh][dg].items()]) for dg in range(2)] for sh in range(2)]
    json.dump(res, open(sys.argv[3], "w"))


if __name__ == "__main__":
    main()
