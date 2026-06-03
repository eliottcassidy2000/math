#!/usr/bin/env python3
"""
lrc_round_count_m89_s574.py    monad-compute-2026-06-02-S574

HYP-1998 open task (A): confirm ROUND tournament iso-class counts via the DIRECT
round generator for m=8,9 (closed form A000016 predicts 16, 30). Prior work
(oracle-S523, lrc_realizable_isoclasses_s523.py) only reached m=7 because that
generator iterated all m^m d-vectors and canonicalized with full m! permutations.

INDEPENDENT METHOD (this script):
  1. Generate ONLY valid round d-vectors via pruned backtracking. A round
     tournament puts vertices on a circle; out-set of vertex i is the clockwise
     arc {i+1,...,i+d_i}. Validity = the arc rule yields a tournament:
         for every pair i<j with k=j-i:  (k<=d_i) XOR ((m-k)<=d_j).
     Backtracking assigns d_0,...,d_{m-1} and checks each new vertex against all
     earlier ones, pruning inconsistent prefixes -> only a few hundred leaves
     even at m=11 (vs m^m brute force).
  2. Count tournament isomorphism classes with an EXACT individualization-
     refinement canonical labeling (no m! blowup; handles vertex-transitive
     cases via individualization). Score(i) = d_i, so most round tournaments
     split immediately under refinement.
  3. Cross-validate: at m<=6 also build the FULL tournament set via brute force
     and check (a) the I-R canon reproduces A000568(m), and (b) round ==
     locally-transitive. This pins the canon's correctness before trusting it
     at m=8,9.
  4. Compare the round count to the A000016 closed form for all m.

No claim is trusted without an independent check: the closed form A000016 and
the direct generator are two separate computations that must agree.
"""
from itertools import combinations, product
from collections import defaultdict
from math import gcd


def _totient(d):
    return sum(1 for k in range(1, d + 1) if gcd(k, d) == 1)


def _divisors(m):
    return [d for d in range(1, m + 1) if m % d == 0]


# ---------- closed form (the value we are independently confirming) ----------
def A000016(m):
    return sum(_totient(d) * 2 ** (m // d)
               for d in _divisors(m) if d % 2 == 1) // (2 * m)


A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}


# ---------- exact individualization-refinement canonical labeling ----------
def _refine(adj, m, colors):
    """1-WL refinement to an equitable coloring; returns compacted int colors."""
    while True:
        sig = []
        for v in range(m):
            outc = sorted(colors[w] for w in range(m) if w != v and adj[v][w])
            inc = sorted(colors[w] for w in range(m) if w != v and adj[w][v])
            sig.append((colors[v], tuple(outc), tuple(inc)))
        ranks = {s: i for i, s in enumerate(sorted(set(sig)))}
        new = [ranks[sig[v]] for v in range(m)]
        if new == colors:
            return colors
        colors = new


def canon(adj, m):
    """Exact canonical flat-adjacency tuple via individualization-refinement."""
    best = [None]

    def rec(colors):
        colors = _refine(adj, m, colors)
        cells = defaultdict(list)
        for v in range(m):
            cells[colors[v]].append(v)
        target = next((c for c in sorted(cells) if len(cells[c]) > 1), None)
        if target is None:                       # discrete coloring -> a labeling
            order = sorted(range(m), key=lambda v: colors[v])
            flat = tuple(adj[order[a]][order[b]]
                         for a in range(m) for b in range(m) if a != b)
            if best[0] is None or flat < best[0]:
                best[0] = flat
            return
        for v in cells[target]:                  # individualize each vertex in cell
            nc = [c * (m + 1) for c in colors]
            nc[v] -= 1                            # rank v just before its cell-mates
            rec(nc)

    # seed with score (= out-degree) coloring
    scores = [sum(adj[v]) for v in range(m)]
    rec(scores)
    return best[0]


def build_adj(d, m):
    adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for k in range(1, d[i] + 1):
            adj[i][(i + k) % m] = 1
    return adj


# ---------- pruned backtracking generator of valid round d-vectors ----------
def valid_dvectors(m):
    """Yield every valid round d-vector on the labeled circle Z_m."""
    d = [0] * m

    def ok(i):
        # check vertex i against all earlier vertices j < i
        di = d[i]
        for j in range(i):
            k = i - j                            # clockwise gap j -> i
            if (k <= d[j]) == ((m - k) <= di):   # both or neither => not a tournament
                return False
        return True

    res = []

    def bt(i):
        if i == m:
            res.append(tuple(d))
            return
        for val in range(m):                     # d_i in {0,...,m-1}
            d[i] = val
            if ok(i):
                bt(i + 1)
    bt(0)
    return res


def round_classes(m):
    """Iso-class count of round tournaments, plus #valid labeled d-vectors."""
    seen = set()
    dv = valid_dvectors(m)
    for d in dv:
        seen.add(canon(build_adj(d, m), m))
    return len(seen), len(dv), seen


# ---------- cross-validation helpers (small m only) ----------
def all_classes(m):
    pairs = list(combinations(range(m), 2))
    seen = {}
    for bits in product((0, 1), repeat=len(pairs)):
        adj = [[0] * m for _ in range(m)]
        for (i, j), b in zip(pairs, bits):
            if b:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        c = canon(adj, m)
        if c not in seen:
            seen[c] = adj
    return seen


def locally_transitive(adj):
    m = len(adj)

    def trans(S):
        S = list(S)
        return sorted(sum(1 for b in S if a != b and adj[a][b]) for a in S) \
            == list(range(len(S)))
    for v in range(m):
        if not trans([w for w in range(m) if w != v and adj[v][w]]):
            return False
        if not trans([w for w in range(m) if w != v and adj[w][v]]):
            return False
    return True


def main():
    print("HYP-1998(A): confirm ROUND == A000016 via the direct generator")
    print("monad-compute-2026-06-02-S574\n")

    # --- correctness pin: I-R canon must reproduce A000568 on full sets (m<=6) ---
    print("CANON CORRECTNESS CHECK (I-R canon vs A000568 on full tournament set):")
    for m in range(3, 7):
        cls = all_classes(m)
        lt = {c for c, a in cls.items() if locally_transitive(a)}
        cnt, ndv, rd = round_classes(m)
        ok_total = (len(cls) == A000568[m])
        ok_lt = (rd == lt)
        print(f"  m={m}: |all-iso|={len(cls):3d} (A000568={A000568[m]:3d}, "
              f"match={ok_total}) | round={cnt} loc-trans={len(lt)} "
              f"round==loctrans:{ok_lt}")
        assert ok_total, f"I-R canon WRONG at m={m}: {len(cls)} != {A000568[m]}"
        assert ok_lt, f"round != locally-transitive at m={m}"
    print("  -> canon verified exact; round == locally-transitive confirmed.\n")

    # --- the main extension: direct round counts vs A000016, m=3..11 ---
    print(f"{'m':>3} {'round(direct)':>13} {'A000016':>8} {'#valid-d':>9} "
          f"{'match':>6}")
    rows = []
    for m in range(3, 12):
        cnt, ndv, _ = round_classes(m)
        a16 = A000016(m)
        rows.append((m, cnt, a16))
        print(f"{m:>3} {cnt:>13d} {a16:>8d} {ndv:>9d} {str(cnt == a16):>6}")

    print("\n ROUND (direct) m=3..11 :", [r[1] for r in rows])
    print(" A000016         m=3..11 :", [r[2] for r in rows])
    all_match = all(r[1] == r[2] for r in rows)
    print(f"\n ALL MATCH: {all_match}")
    if all_match:
        print(" => HYP-1998(A) CONFIRMED: round(8)=16, round(9)=30 reproduced by the")
        print("    direct generator; A000016 closed form holds m=3..11.")


if __name__ == "__main__":
    main()
