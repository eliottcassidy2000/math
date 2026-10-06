#!/usr/bin/env python3
"""
chessboard_weave_20261006_hp_blocking_theory.py

Checks for the proof attempts on HYP-9168 (beta(T) = hall(T)).

 A. Gutin's theorem (SIAM J. Discrete Math. 6 (1993), Corollary 1: a complete multipartite
    digraph has a Hamiltonian path iff it has a spanning path + disjoint cycles) on every
    multipartite tournament obtainable as T - (arcs inside the parts of a partition), all
    classes T with N <= 7 and all set partitions.  Consequences checked: cluster-type
    blocking sets are never cheaper than hall; the cheapest cluster-type blocking set vs beta.
 B. MERGE LEMMA.  P = p_1..p_m a path, C = c_0..c_{k-1} a cycle of a tournament T, disjoint.
    The "simple merges" of P and C are: c_i -> p_1 (front), p_m -> c_i (back),
    p_j -> c_i & c_{i-1} -> p_{j+1} (middle).  For the diagonal class delta (pairs {p_j,c_i}
    with i+j = delta mod k) let x_j = [p_j -> c_{delta-j} in T].  Claim: the minimum number of
    arcs of E(V(P),V(C)) whose deletion blocks every simple merge is
         MB(P,C) = sum_delta (1 + asc_delta),   asc_delta = #{j<m : x_j = 0, x_{j+1} = 1}
    (>= k).  Verified against brute force on random instances.
 C. The single-factor merge argument would prove HYP-9168 if
         (MB inequality)  hall(T) <= sum_i MB(P, C_i)  for every spanning 1-path-cycle factor
    P + C_1 + ... + C_r (r >= 1) of T.  Tested on all classes N <= 7; failures show where the
    argument breaks (explicit witness: a set D_F of size sum MB blocking every simple merge of F
    while T - D_F still has a Hamiltonian path).
 D. Census of the minimum blocking sets D that leave a 1-path-cycle factor (N <= 7) by shape.
 E. k-path generalization: beta_k (fewest deletions leaving no cover by <= k vertex-disjoint
    paths) vs hall_k = min_{|A|+|B| = N+k+1} e(A,B) (fewest deletions leaving no path-cycle
    subgraph with <= k paths), k = 2, 3, all classes N <= 7 (brute force cross-check N <= 6).
 F. The analogue for semicomplete digraphs (2-cycles allowed), all labeled ones with N <= 5.
 G. On the exotic minima (N <= 7): the merge count of a fewest-cycle (hence stuck) factor vs |D|.
"""
import itertools
import os
import random
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "chessboard_weave_20261006_hp_blocking_engine.c")
BIN = os.path.join(tempfile.gettempdir(), "chessboard_weave_20261006_hp_blocking_theory_%d.bin" % os.getpid())
T0 = time.time()


def build():
    subprocess.run(["clang", "-O3", "-march=native", "-o", BIN, SRC], check=True)


def engine(args, lines):
    r = subprocess.run([BIN] + list(args), input="\n".join(lines) + "\n", capture_output=True, text=True, check=True)
    return [ln for ln in r.stdout.splitlines() if ln.strip()]


def field(line, key):
    for tok in line.split():
        if tok.startswith(key + "="):
            return tok.split("=", 1)[1]
    return None


def classes(n):
    return subprocess.run(["gentourng", str(n)], capture_output=True, text=True).stdout.split()


def from_string(s):
    L = len(s)
    N = 1
    while N * (N - 1) // 2 < L:
        N += 1
    adj = [[False] * N for _ in range(N)]
    k = 0
    for i in range(N):
        for j in range(i + 1, N):
            if s[k] == "1":
                adj[i][j] = True
            else:
                adj[j][i] = True
            k += 1
    return N, adj


def hall(N, adj):
    best = 10 ** 9
    for A in range(1 << N):
        a = bin(A).count("1")
        b = N + 2 - a
        if a < 2 or b > N:
            continue
        d = sorted(sum(1 for u in range(N) if (A >> u) & 1 and adj[u][v]) for v in range(N))
        best = min(best, sum(d[:b]))
    return best


def has_hp(N, adj, removed=frozenset()):
    ends = {}
    for v in range(N):
        ends[1 << v] = 1 << v
    full = (1 << N) - 1
    for S in range(1, full + 1):
        e = ends.get(S, 0)
        if not e:
            continue
        for w in range(N):
            if (S >> w) & 1:
                continue
            for u in range(N):
                if (e >> u) & 1 and adj[u][w] and (u, w) not in removed:
                    ends[S | (1 << w)] = ends.get(S | (1 << w), 0) | (1 << w)
                    break
    return ends.get(full, 0) != 0


def find_hp(N, adj, removed=frozenset()):
    for perm in itertools.permutations(range(N)):
        if all(adj[perm[i]][perm[i + 1]] and (perm[i], perm[i + 1]) not in removed for i in range(N - 1)):
            return perm
    return None


# ------------------------------------------------------------------ merge lemma
def mb(adj, P, C):
    m, k = len(P), len(C)
    tot = 0
    for delta in range(k):
        x = [1 if adj[P[j - 1]][C[(delta - j) % k]] else 0 for j in range(1, m + 1)]
        asc = sum(1 for j in range(m - 1) if x[j] == 0 and x[j + 1] == 1)
        tot += 1 + asc
    return tot


def simple_merge_exists(adj, P, C, removed):
    m, k = len(P), len(C)

    def arc(u, v):
        return adj[u][v] and (u, v) not in removed

    if any(arc(c, P[0]) for c in C):
        return True
    if any(arc(P[-1], c) for c in C):
        return True
    for j in range(m - 1):
        for i in range(k):
            if arc(P[j], C[i]) and arc(C[(i - 1) % k], P[j + 1]):
                return True
    return False


def min_block_simple(adj, P, C, limit):
    pairs = [(p, c) if adj[p][c] else (c, p) for p in P for c in C]
    for s in range(0, limit + 1):
        for D in itertools.combinations(pairs, s):
            if not simple_merge_exists(adj, P, C, frozenset(D)):
                return s, D
    return None, None


def all_cycles(N, adj):
    cyc = []
    for L in range(3, N + 1):
        for S in itertools.combinations(range(N), L):
            first = S[0]
            for rest in itertools.permutations(S[1:]):
                c = (first,) + rest
                if all(adj[c[i]][c[(i + 1) % L]] for i in range(L)):
                    cyc.append(c)
    return cyc


def hps_of(adj, verts):
    out = []
    for perm in itertools.permutations(verts):
        if all(adj[perm[i]][perm[i + 1]] for i in range(len(perm) - 1)):
            out.append(perm)
    return out


def min_mb_over_factors(N, adj):
    """min over spanning 1-path-cycle factors with >= 1 cycle of sum_i MB(P, C_i); returns (value, factor)."""
    cyc = all_cycles(N, adj)
    best = (10 ** 9, None)
    hp_cache = {}

    def rec(start, used, chosen):
        nonlocal best
        if chosen:
            rest = tuple(v for v in range(N) if not (used >> v) & 1)
            if rest:
                if rest not in hp_cache:
                    hp_cache[rest] = hps_of(adj, rest)
                for P in hp_cache[rest]:
                    val = sum(mb(adj, P, C) for C in chosen)
                    if val < best[0]:
                        best = (val, (P, tuple(chosen)))
        for idx in range(start, len(cyc)):
            c = cyc[idx]
            mask = sum(1 << v for v in c)
            if mask & used:
                continue
            rec(idx + 1, used | mask, chosen + [c])

    rec(0, 0, [])
    return best


def blocking_set_for_factor(adj, P, C):
    """one deletion per descent on every diagonal (realizes MB)."""
    m, k = len(P), len(C)
    D = []
    for delta in range(k):
        x = [1] + [1 if adj[P[j - 1]][C[(delta - j) % k]] else 0 for j in range(1, m + 1)] + [0]
        for j in range(0, m + 1):
            if x[j] == 1 and x[j + 1] == 0:
                pos = j + 1 if j + 1 <= m else j  # delete the non-virtual position
                if j + 1 > m:
                    pos = j
                p = P[pos - 1]
                c = C[(delta - pos) % k]
                D.append((p, c) if adj[p][c] else (c, p))
    return D


def section(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)
    sys.stdout.flush()


def main(SECTIONS="ABCDEFG"):
    build()
    print(__doc__.strip())
    print(f"sections run: {SECTIONS}")
    rng = random.Random(4242)

    # ------------------------------------------------------------------ A
    if "A" in SECTIONS:
        section("A. Gutin's Corollary 1 on all multipartite tournaments T - (intra-part arcs), N <= 7")
        for n in range(3, 8):
            res = engine(["gutin"], classes(n))
            mism = sum(int(field(ln, "hp_vs_1pcf_mismatch")) for ln in res)
            parts = sum(int(field(ln, "partitions")) for ln in res)
            below = sum(1 for ln in res if int(field(ln, "min_cluster_blocking")) < int(field(ln, "hall")))
            strict = sum(1 for ln in res if int(field(ln, "min_cluster_blocking")) > int(field(ln, "beta")))
            print(f"N={n}: {len(res)} classes x all set partitions = {parts} multipartite tournaments; "
                  f"HP <=> 1-path-cycle factor fails on {mism}; cheapest cluster-type blocking set < hall: {below} classes; "
                  f"> beta (so no minimum blocking set is cluster-type): {strict} classes")
        sys.stdout.flush()

    # ------------------------------------------------------------------ B
    if "B" in SECTIONS:
        section("B. merge lemma: min deletions blocking all simple merges of (P, C) == MB(P,C) = sum_delta (1 + asc_delta)")
        ok = tot = 0
        hist = {}
        for trial in range(1500):
            while True:
                m = rng.randint(1, 4)
                k = rng.randint(3, 5)
                if m * k <= 12:
                    break
            N = m + k
            # random tournament containing the cycle C on vertices m..m+k-1 and path P on 0..m-1
            adj = [[False] * N for _ in range(N)]
            for i in range(N):
                for j in range(i + 1, N):
                    if rng.random() < 0.5:
                        adj[i][j] = True
                    else:
                        adj[j][i] = True
            P = list(range(m))
            C = list(range(m, N))
            rng.shuffle(P)
            rng.shuffle(C)
            for i in range(m - 1):
                a, b = P[i], P[i + 1]
                adj[a][b], adj[b][a] = True, False
            for i in range(k):
                a, b = C[i], C[(i + 1) % k]
                adj[a][b], adj[b][a] = True, False
            val = mb(adj, P, C)
            s, D = min_block_simple(adj, P, C, val)
            tot += 1
            ok += (s == val)
            hist[val - k] = hist.get(val - k, 0) + 1
        print(f"{tot} random (T, P, C) with |P| = 1..4, |C| = 3..5, |P||C| <= 12: brute-force minimum == MB formula on {ok}/{tot}; "
              f"distribution of MB - |C| (= total ascents): {dict(sorted(hist.items()))}")
        sys.stdout.flush()

    # ------------------------------------------------------------------ C
    if "C" in SECTIONS:
        section("C. MB inequality hall(T) <= min_F sum_i MB(P, C_i) over spanning 1-path-cycle factors F (>= 1 cycle)")
        example = None
        for n in range(3, 8):
            cl = classes(n)
            fails = 0
            nocyc = 0
            tight = 0
            worst = 0
            for s in cl:
                N, adj = from_string(s)
                h = hall(N, adj)
                val, F = min_mb_over_factors(N, adj)
                if F is None:
                    nocyc += 1
                    continue
                if val < h:
                    fails += 1
                    worst = max(worst, h - val)
                    if example is None or (N, h - val) > (example[0], example[3] - example[2]):
                        example = (N, s, val, h, F)
                elif val == h:
                    tight += 1
            print(f"N={n}: {len(cl)} classes; no 1-path-cycle factor with a cycle: {nocyc}; "
                  f"min_F sum MB < hall (merge count alone insufficient): {fails} (max deficit {worst}); == hall: {tight}")
            sys.stdout.flush()
        if example:
            N, s, val, h, (P, Cs) = example
            _, adj = from_string(s)
            D = []
            for C in Cs:
                D += blocking_set_for_factor(adj, P, C)
            Dset = frozenset(D)
            stuck = all(not simple_merge_exists(adj, P, C, Dset) for C in Cs)
            hp = find_hp(N, adj, Dset)
            print(f"witness: T = {s} (N={N}, hall={h}); F = path {P} + cycles {list(Cs)}; sum MB = {val}")
            print(f"  D_F = {sorted(D)} (|D_F| = {len(D)}) blocks every simple merge of F: {stuck}; "
                  f"yet T - D_F has the Hamiltonian path {hp}")

    # ------------------------------------------------------------------ D
    if "D" in SECTIONS:
        section("D. minimum blocking sets leaving a 1-path-cycle factor (exotic minima), N <= 8")
        for n in range(3, 9):
            res = engine(["exotic"], classes(n))
            shapes = {}
            regular_like = set()
            for ln in res:
                Nn, adj = from_string(ln.split()[0])
                sc = sorted(sum(adj[u]) for u in range(Nn))
                regular_like.add((ln.split()[0], tuple(sc), field(ln, "hall")))
            for ln in res:
                Dstr = field(ln, "D").split(",")
                arcs = [tuple(map(int, a.split(">"))) for a in Dstr]
                N = int(field(ln, "N"))
                deg = [0] * N
                for (u, v) in arcs:
                    deg[u] += 1
                    deg[v] += 1
                beta = int(field(ln, "beta"))
                if max(deg) == N - 1 and beta == N - 1:
                    shape = "isolated vertex (all N-1 arcs at one vertex)"
                else:
                    # source/sink created?
                    Nn, adj = from_string(ln.split()[0])
                    rem = set(arcs)
                    srcs = [v for v in range(N) if not any(adj[u][v] and (u, v) not in rem for u in range(N))]
                    snks = [v for v in range(N) if not any(adj[v][w] and (v, w) not in rem for w in range(N))]
                    shape = f"{len(srcs)} source(s) + {len(snks)} sink(s) in T-D"
                shapes[shape] = shapes.get(shape, 0) + 1
            ncls = len(set(ln.split()[0] for ln in res))
            print(f"N={n}: {len(res)} exotic minimum blocking sets in {ncls} classes; shapes: {shapes}")
            for (t, sc, h) in sorted(regular_like):
                print(f"     class {t}: score sequence {list(sc)}, hall = beta = {h}")
        sys.stdout.flush()

    # ------------------------------------------------------------------ E
    if "E" in SECTIONS:
        section("E. k-path generalization beta_k = hall_k (k = 2, 3)")
        for n in range(3, 8):
            cl = classes(n)
            for K in (2, 3):
                if K >= n:   # beta_K undefined (every vertex alone is a cover by n <= K paths)
                    continue
                brute = "1" if n <= 6 else "0"
                res = engine(["kpath", "0", str(K), brute], cl)
                eq = sum(1 for ln in res if ln.endswith(" EQ"))
                bm = sum(1 for ln in res if n > 6 or field(ln, "beta_K_brute") == field(ln, "beta_K"))
                print(f"N={n} k={K}: beta_k == hall_k on {eq}/{len(res)} classes"
                      + (f"; B&B == brute force on {bm}/{len(res)}" if n <= 6 else ""))
                for ln in res:
                    if not ln.endswith(" EQ"):
                        print("   NEQ: " + ln)
            sys.stdout.flush()
    # ------------------------------------------------------------------ F
    if "F" in SECTIONS:
        section("F. the statement is special to tournaments: semicomplete digraphs (2-cycles allowed)")
        for n in (3, 4, 5):
            E = n * (n - 1) // 2
            lines = ["".join(t) for t in itertools.product("012", repeat=E) if "2" in t]
            res = engine(["bnb"], lines)
            ce = [ln for ln in res if "COUNTEREXAMPLE" in ln]
            print(f"N={n}: all {len(lines)} labeled semicomplete digraphs with >= 1 two-cycle: beta < hall for {len(ce)}"
                  + (f"; e.g. {ce[0].split()[0]} (pairs i<j: 1 = i->j, 0 = j->i, 2 = both) hall={field(ce[0], 'hall')} "
                     f"beta={field(ce[0], 'beta')} blocking set {field(ce[0], 'witness')}" if ce else ""))
            sys.stdout.flush()
    # ------------------------------------------------------------------ G
    if "G" in SECTIONS:
        section("G. merge count on the exotic minima: fewest-cycle 1-path-cycle factor F of T - D (stuck), |D| vs sum MB(F)")
        for n in range(5, 8):
            res = engine(["exotic"], classes(n))
            stats = {}
            for ln in res:
                t = ln.split()[0]
                N, adj = from_string(t)
                D = set(tuple(map(int, a.split(">"))) for a in field(ln, "D").split(","))
                gadj = [[adj[u][v] and (u, v) not in D for v in range(N)] for u in range(N)]
                cyc = all_cycles(N, gadj)
                best_r = None
                vals = []

                def rec(start, used, chosen):
                    nonlocal best_r, vals
                    rest = tuple(v for v in range(N) if not (used >> v) & 1)
                    if chosen and rest:
                        for P in hps_of(gadj, rest):
                            r = len(chosen)
                            val = sum(mb(adj, P, C) for C in chosen)
                            if best_r is None or r < best_r:
                                best_r, vals = r, [val]
                            elif r == best_r:
                                vals.append(val)
                    for idx in range(start, len(cyc)):
                        c = cyc[idx]
                        m = sum(1 << v for v in c)
                        if not (m & used):
                            rec(idx + 1, used | m, chosen + [c])

                rec(0, 0, [])
                key = (len(D), best_r, min(vals), max(vals))
                stats[key] = stats.get(key, 0) + 1
            print(f"N={n}: (|D|, fewest cycles r, min sum MB, max sum MB over such factors) -> count: {dict(sorted(stats.items()))}")
            sys.stdout.flush()
    print(f"\ntotal time {time.time() - T0:.1f}s")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "ABCDEFG")
