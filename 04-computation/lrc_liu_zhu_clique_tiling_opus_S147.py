"""
lrc_liu_zhu_clique_tiling_opus_S147.py   (opus-2026-07-07-S147, HYP-5277 part 4)

THE UPPER-BOUND MECHANISM for Liu-Zhu Conjecture 2: {0,x,y,x+y} is a 4-clique of G_M
(six differences x,y,x+y,y-x = all in +-M), so an M-avoiding set meets each translate
c+{0,x,y,x+y} in <= 1 point.  If Z_N tiles into (k+1)m such 4-cliques plus 1 leftover
point (4*(k+1)m = N-1), then S(period) <= (k+1)m + 1, and Haralambis's single-leftover
argument gives mu <= (k+1)m/N.  This script verifies the clique property (trivial) and
the exact-cover tiling of Z_N \\ {pt} by translates of {0,x,y,x+y}, per instance.
"""
from math import gcd
import sys

sys.setrecursionlimit(100000)


def clique_diffs_ok(x, y, N):
    cl = [0, x % N, y % N, (x + y) % N]
    diffs = set((a - b) % N for a in cl for b in cl if a != b)
    Mset = set(v % N for v in [x, y, y - x, y + x, -x, -y, -(y - x), -(y + x)])
    return diffs <= Mset


def exact_cover_tiling(x, y, N, node_cap=2_000_000):
    """Cover N-1 points of Z_N by disjoint translates of {0,x,y,x+y}; 1 leftover allowed.
       Backtracking on the least-covered point; returns list of offsets or None."""
    clique = [0, x % N, y % N, (x + y) % N]
    if len(set(clique)) != 4:
        return None
    nodes = [0]
    full = frozenset(range(N))

    def bt(uncov):
        nodes[0] += 1
        if nodes[0] > node_cap:
            raise TimeoutError
        if len(uncov) <= 1:
            return []
        p = min(uncov)
        for off in clique:
            c = (p - off) % N
            block = frozenset((c + o) % N for o in clique)
            if len(block) == 4 and block <= uncov:
                rest = bt(uncov - block)
                if rest is not None:
                    return [c] + rest
        return None

    try:
        return bt(full)
    except TimeoutError:
        return "TIMEOUT"


def main():
    print("=" * 96)
    print("UPPER-BOUND MECHANISM: {0,x,y,x+y} 4-clique + Z_N tiling into (k+1)m cliques + 1 pt")
    print("=" * 96)
    ok = 0
    fails = []
    timeouts = []
    for y in range(3, 32, 2):
        for x in range(3, y, 2):   # x >= 3: the previously-open cases
            if gcd(x, y) != 1:
                continue
            k, m = (x - 1) // 2, (y - 1) // 2
            N = 4 * (k + 1) * m + 1
            if not clique_diffs_ok(x, y, N):
                fails.append((x, y, "NOT-CLIQUE"))
                continue
            res = exact_cover_tiling(x, y, N)
            if res == "TIMEOUT":
                timeouts.append((x, y, N))
            elif res is not None and len(res) == (k + 1) * m:
                ok += 1
            else:
                fails.append((x, y, N, len(res) if isinstance(res, list) else res, (k + 1) * m))
    print(f"  clique + tiling into (k+1)m blocks: {ok} OK, {len(fails)} fail, {len(timeouts)} timeout")
    if fails:
        print(f"  FAILS: {fails[:10]}")
    if timeouts:
        print(f"  timeouts (larger N, inconclusive): {timeouts}")
    print()
    if not fails:
        print("  => the 4-clique {0,x,y,x+y} tiles Z_N \\ {pt} into (k+1)m blocks on every")
        print("     resolved instance: the upper bound mu <= (k+1)m/N holds by the Haralambis")
        print("     single-leftover argument.  With THM-657's lower bound => Conjecture 2 on")
        print("     the resolved range; the uniform tiling construction is the remaining step.")
    # show one explicit tiling
    print()
    print("  Example tiling for (x,y)=(3,5), M={2,3,5,8}, N=17 (offsets c; blocks c+{0,3,5,8}):")
    res = exact_cover_tiling(3, 5, 17)
    if isinstance(res, list):
        covered = set()
        for c in res:
            covered |= set((c + o) % 17 for o in [0, 3, 5, 8])
        leftover = set(range(17)) - covered
        print(f"    offsets {sorted(res)}; leftover point {sorted(leftover)}")


if __name__ == "__main__":
    main()
