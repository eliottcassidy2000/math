#!/usr/bin/env python3
"""
chessboard_weave_20261006_hp_blocking_run.py

HYP-9168: beta(T) (fewest arcs whose deletion kills every Hamiltonian path) equals
hall(T) (Hall-deficiency bound) for every tournament T.

Parts
  1. Konig identity: hall(T) = min_{|A|+|B|=N+2} e_T(A,B)  (= min deletions killing every
     1-path-cycle factor) against the literal (X,Y) definition, in- AND out-versions,
     Y allowed to meet X.  Pure-Python reference on all labeled N <= 5, C engine on all
     labeled N <= 6, all classes N = 7, random labeled N = 7, 8.
  2. Exact beta: branch and bound (C engine) against brute-force subset enumeration on all
     classes N <= 7 (infeasibility of hall-1 AND iterative deepening), and on 5000
     lightly damaged tournaments (general oriented graphs, so feasible answers occur).
  3. beta = hall at N = 10..14: uniform random labeled tournaments, random regular /
     almost-regular tournaments (score-preserving 3-cycle-reversal walk), and structured
     families (all circulants on Z_11 and Z_13, Paley QR_11, QR_11 - v, transitive,
     one/two reversed arcs of the transitive tournament, compositions with cyclic
     triangles and regular blocks, lexicographic products).
  4. Gutin route data: for every class N <= 7, ALL minimum blocking sets D (|D| = beta):
     how many are of "multipartite type" (deleted pairs = disjoint union of cliques, so
     T - D is a multipartite tournament) and how many leave a 1-path-cycle factor.
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
BIN = os.path.join(tempfile.gettempdir(), "chessboard_weave_20261006_hp_blocking_run_%d.bin" % os.getpid())

T0 = time.time()


def build():
    subprocess.run(["clang", "-O3", "-march=native", "-o", BIN, SRC], check=True)


def engine(mode, lines, extra=()):
    inp = "\n".join(lines) + "\n"
    r = subprocess.run([BIN, mode] + list(extra), input=inp, capture_output=True, text=True, check=True)
    return [ln for ln in r.stdout.splitlines() if ln.strip()]


def field(line, key):
    for tok in line.split():
        if tok.startswith(key + "="):
            return tok.split("=", 1)[1]
    return None


# ---------------------------------------------------------------- tournaments as strings
def to_string(N, out):
    """out: set of arcs (u,v). upper-triangle gentourng string."""
    return "".join("1" if (i, j) in out else "0" for i in range(N) for j in range(i + 1, N))


def from_string(s):
    L = len(s)
    N = 1
    while N * (N - 1) // 2 < L:
        N += 1
    arcs = set()
    k = 0
    for i in range(N):
        for j in range(i + 1, N):
            arcs.add((i, j) if s[k] == "1" else (j, i))
            k += 1
    return N, arcs


def adj_lists(N, arcs):
    outl = {v: sorted(w for (u, w) in arcs if u == v) for v in range(N)}
    return outl


def random_tournament(N, rng):
    return "".join(rng.choice("01") for _ in range(N * (N - 1) // 2))


def circulant(N, S):
    S = set(x % N for x in S)
    assert all(((-s) % N) not in S for s in S) and len(S) == (N - 1) // 2
    arcs = {(i, j) for i in range(N) for j in range(N) if i != j and (j - i) % N in S}
    return arcs


def transitive(N):
    return {(i, j) for i in range(N) for j in range(i + 1, N)}


def random_regularish(N, rng, steps=None):
    """score-preserving random walk (reverse directed 3-cycles) from a circulant / near-regular
    start; for odd N regular, for even N almost regular (scores N/2, N/2-1)."""
    if N % 2 == 1:
        arcs = circulant(N, range(1, (N - 1) // 2 + 1))
    else:
        # regular tournament on N-1 (odd) vertices plus a vertex beating N/2 - 1 of them:
        # scores are N/2 - 1 or N/2 (almost regular)
        base = circulant(N - 1, range(1, (N - 2) // 2 + 1))
        arcs = set(base)
        others = list(range(N - 1))
        rng.shuffle(others)
        for idx, v in enumerate(others):
            if idx < N // 2 - 1:
                arcs.add((N - 1, v))
            else:
                arcs.add((v, N - 1))
    out = {v: set() for v in range(N)}
    for (u, v) in arcs:
        out[u].add(v)
    steps = steps or 50 * N * N
    verts = list(range(N))
    for _ in range(steps):
        a, b, c = rng.sample(verts, 3)
        # orient so that we check a->b->c->a or reverse
        if b in out[a] and c in out[b] and a in out[c]:
            out[a].remove(b); out[b].remove(c); out[c].remove(a)
            out[b].add(a); out[c].add(b); out[a].add(c)
        elif a in out[b] and b in out[c] and c in out[a]:
            out[b].remove(a); out[c].remove(b); out[a].remove(c)
            out[a].add(b); out[b].add(c); out[c].add(a)
    arcs = {(u, v) for u in out for v in out[u]}
    return arcs


def automorphisms(N, arcs, limit=256):
    import networkx as nx
    from networkx.algorithms.isomorphism import DiGraphMatcher
    G = nx.DiGraph()
    G.add_nodes_from(range(N))
    G.add_edges_from(arcs)
    auts = []
    for m in DiGraphMatcher(G, G).isomorphisms_iter():
        auts.append([m[i] for i in range(N)])
        if len(auts) > limit:
            return None
    return auts


def line_with_auts(N, arcs, use_aut=True):
    s = to_string(N, arcs)
    if not use_aut:
        return s
    auts = automorphisms(N, arcs)
    if auts is None or len(auts) <= 1:
        return s
    return s + " " + " ".join("A:" + ",".join(map(str, a)) for a in auts)


# ---------------------------------------------------------------- pure-python reference (N <= 5)
def py_hall_brute(N, arcs):
    inn = {v: {u for (u, w) in arcs if w == v} for v in range(N)}
    outn = {v: {w for (u, w) in arcs if u == v} for v in range(N)}
    V = range(N)
    best_in = best_out = 10 ** 9
    subsets = [frozenset(c) for r in range(N + 1) for c in itertools.combinations(V, r)]
    for X in subsets:
        if len(X) < 2:
            continue
        for Y in subsets:
            if len(Y) > len(X) - 2:
                continue
            ci = sum(len(inn[x] - Y) for x in X)
            co = sum(len(outn[x] - Y) for x in X)
            best_in = min(best_in, ci)
            best_out = min(best_out, co)
    return best_in, best_out


def py_hall_konig(N, arcs):
    V = range(N)
    subsets = [frozenset(c) for r in range(N + 1) for c in itertools.combinations(V, r)]
    best = 10 ** 9
    for A in subsets:
        for B in subsets:
            if len(A) + len(B) != N + 2:
                continue
            best = min(best, sum(1 for (u, v) in arcs if u in A and v in B))
    return best


def py_min_deletions_kill_1pcf(N, arcs):
    """min |D| such that T-D has no 1-path-cycle factor, by brute force over D (N <= 4)."""
    import networkx as nx
    arcl = sorted(arcs)
    for k in range(len(arcl) + 1):
        for D in itertools.combinations(arcl, k):
            rest = set(arcl) - set(D)
            B = nx.Graph()
            B.add_nodes_from([("o", v) for v in range(N)] + [("i", v) for v in range(N)])
            B.add_edges_from((("o", u), ("i", v)) for (u, v) in rest)
            m = nx.bipartite.maximum_matching(B, top_nodes=[("o", v) for v in range(N)])
            if len(m) // 2 <= N - 2:
                return k
    return None


def section(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)
    sys.stdout.flush()


def main():
    build()
    rng = random.Random(20261006)
    print(__doc__.strip())
    allok = True

    # ------------------------------------------------------------------ PART 1
    section("PART 1. Konig identity hall = min_{|A|+|B|=N+2} e(A,B); in/out versions agree")
    # 1a pure python, all labeled N <= 5 (and the 1-pcf brute force for N <= 4)
    cnt = 0
    bad = 0
    for N in range(2, 6):
        E = N * (N - 1) // 2
        for bits in range(2 ** E):
            s = format(bits, "0%db" % E) if E else ""
            if N == 2:
                s = s or "1"
            N_, arcs = from_string(s) if E else (N, set())
            hi, ho = py_hall_brute(N, arcs)
            hk = py_hall_konig(N, arcs)
            ok = (hi == ho == hk)
            if N <= 4:
                ok &= (py_min_deletions_kill_1pcf(N, arcs) == hk)
            cnt += 1
            bad += (not ok)
    print(f"1a pure Python, all labeled tournaments N=2..5 ({cnt}): hall_in == hall_out == Konig formula"
          f" (and == brute-force min deletions killing all 1-path-cycle factors for N <= 4): "
          f"{'ALL AGREE' if bad == 0 else str(bad) + ' DISAGREE'}")
    allok &= bad == 0
    # 1b C engine, all labeled N <= 6, classes N = 7, random labeled N = 7, 8
    for N in range(3, 7):
        E = N * (N - 1) // 2
        lines = [format(b, "0%db" % E) for b in range(2 ** E)]
        res = engine("hall", lines)
        ag = sum(1 for ln in res if ln.endswith("AGREE") and not ln.endswith("DISAGREE") and field(ln, "obstruction_ok") == "1")
        print(f"1b C engine, all {len(lines)} labeled tournaments N={N}: formula == in-version == out-version"
              f" and optimal (A,B) kills every HP and every 1-path-cycle factor: {ag}/{len(res)}")
        allok &= ag == len(lines)
    g7 = subprocess.run(["gentourng", "7"], capture_output=True, text=True).stdout.split()
    res = engine("hall", g7)
    ag = sum(1 for ln in res if ln.endswith(" AGREE") and field(ln, "obstruction_ok") == "1")
    print(f"1b C engine, all {len(g7)} classes N=7: agree {ag}/{len(res)}")
    allok &= ag == len(g7)
    for N, M in ((7, 3000), (8, 1000)):
        lines = [random_tournament(N, rng) for _ in range(M)]
        res = engine("hall", lines)
        ag = sum(1 for ln in res if ln.endswith(" AGREE") and field(ln, "obstruction_ok") == "1")
        print(f"1b C engine, {M} random labeled N={N}: agree {ag}/{len(res)}")
        allok &= ag == M

    # ------------------------------------------------------------------ PART 2
    section("PART 2. exact beta: branch and bound vs brute-force subset enumeration")
    for N in range(3, 8):
        gl = subprocess.run(["gentourng", str(N)], capture_output=True, text=True).stdout.split()
        r1 = engine("both", gl)
        m1 = sum(1 for ln in r1 if " EQ " in ln and ln.endswith("MATCH") and not ln.endswith("MISMATCH"))
        r2 = engine("deep", gl)
        m2 = sum(1 for ln in r2 if ln.endswith(" MATCH"))
        eq = sum(1 for ln in r1 if field(ln, "beta_brute") == field(ln, "hall"))
        print(f"N={N}: {len(gl)} classes; brute-force beta == hall: {eq}/{len(gl)};"
              f" B&B(hall-1 infeasible) == brute: {m1}/{len(gl)}; B&B deepening == brute: {m2}/{len(gl)}")
        allok &= (m1 == m2 == eq == len(gl))
    # root symmetry breaking (arc orbits of Aut(T)) validated: deepening WITH symmetry vs brute force
    tot = msym = 0
    for N in range(3, 8):
        gl = subprocess.run(["gentourng", str(N)], capture_output=True, text=True).stdout.split()
        lines = []
        for s in gl:
            _, arcs = from_string(s)
            ln = line_with_auts(N, arcs)
            if " A:" in ln:
                lines.append(ln)
        if lines:
            r = engine("deepsym", lines)
            tot += len(lines)
            msym += sum(1 for ln in r if ln.endswith(" MATCH"))
    print(f"classes N<=7 with nontrivial Aut: {tot}; deepening WITH root symmetry breaking == brute: {msym}/{tot}")
    allok &= msym == tot
    # damaged tournaments (general oriented graphs): feasible answers occur
    lines = []
    for N, M in ((6, 2000), (7, 2000), (8, 1000)):
        for _ in range(M):
            s = random_tournament(N, rng)
            _, arcs = from_string(s)
            m = rng.randint(0, N // 2)
            D = rng.sample(sorted(arcs), m)
            lines.append(s + ((" D:" + ",".join("%d>%d" % a for a in D)) if D else ""))
    res = engine("deep", lines)
    m = sum(1 for ln in res if ln.endswith(" MATCH"))
    hist = {}
    neq = []
    for ln, src in zip(res, lines):
        b = int(field(ln, "beta_deep"))
        hist[b] = hist.get(b, 0) + 1
        if field(ln, "hall") != field(ln, "beta_deep"):
            neq.append((src, field(ln, "hall"), field(ln, "beta_deep")))
    print(f"damaged tournaments (T minus 0..N/2 random arcs), N=6,7,8, {len(lines)} digraphs: "
          f"B&B deepening == brute force on {m}/{len(lines)}; beta histogram {dict(sorted(hist.items()))}")
    allok &= m == len(lines)
    one_arc = [x for x in neq if x[0].count(">") == 1]
    print(f"  (side remark) beta < hall on {len(neq)} of these NON-tournaments, e.g. T minus ONE arc: "
          + "; ".join(f"{a} hall={h} beta={b}" for a, h, b in one_arc[:3]))

    # ------------------------------------------------------------------ PART 3
    section("PART 3. beta = hall at N = 10..14")
    summary = {}

    def run_batch(name, lines, extra=()):
        t = time.time()
        res = engine("bnb", lines, extra)
        eq = sum(1 for ln in res if ln.endswith(" EQ") or " EQ " in ln)
        ce = [ln for ln in res if "COUNTEREXAMPLE" in ln or "WITNESS_FAIL" in ln or "ABORTED" in ln]
        ob = sum(1 for ln in res if field(ln, "obstruction_ok") == "1")
        halls = {}
        for ln in res:
            h = int(field(ln, "hall"))
            halls[h] = halls.get(h, 0) + 1
        mx = max(res, key=lambda ln: int(field(ln, "nodes")))
        lt = sum(1 for ln in res if int(field(ln, "hall")) < int(field(ln, "sigma")))
        print(f"{name}: {len(lines)} tournaments; beta == hall: {eq}/{len(lines)}; hall obstruction verified: "
              f"{ob}/{len(lines)}; hall<sigma: {lt}; hall histogram {dict(sorted(halls.items()))}; "
              f"max nodes {field(mx, 'nodes')}; {time.time() - t:.1f}s")
        for ln in ce:
            print("   !!! " + ln)
        sys.stdout.flush()
        summary[name] = (eq, len(lines), ce)
        return res

    # 3a uniform random labeled
    for N, M in ((10, 10000), (11, 10000), (12, 10000), (13, 1000), (14, 200)):
        lines = [random_tournament(N, rng) for _ in range(M)]
        run_batch(f"3a uniform random labeled N={N}", lines)
    # 3b random regular / almost regular
    for N, M in ((10, 300), (11, 300), (12, 200), (13, 20)):
        lines = [to_string(N, random_regularish(N, rng)) for _ in range(M)]
        run_batch(f"3b random {'regular' if N % 2 else 'almost-regular'} N={N}", lines)
    # 3c structured
    lines = []
    for S in itertools.product(*[(s, -s) for s in range(1, 6)]):
        lines.append(line_with_auts(11, circulant(11, S)))
    run_batch("3c all 32 circulant tournaments on Z_11 (regular; prediction beta = 10)", lines)
    lines = []
    for S in itertools.product(*[(s, -s) for s in range(1, 7)]):
        lines.append(line_with_auts(13, circulant(13, S)))
    run_batch("3c all 64 circulant tournaments on Z_13 (regular; prediction beta = 12)", lines)
    qr11 = circulant(11, [1, 3, 4, 5, 9])
    qr11mv = {(u - 1, v - 1) for (u, v) in qr11 if u != 0 and v != 0}
    qr7 = circulant(7, [1, 2, 4])
    lines = [line_with_auts(11, qr11), line_with_auts(10, qr11mv)]
    run_batch("3c Paley QR_11 and QR_11 minus a vertex", lines)
    # transitive and near-transitive
    lines = []
    for N in (10, 11, 12):
        TT = transitive(N)
        lines.append(to_string(N, TT))
        for (i, j) in sorted(TT):
            lines.append(to_string(N, (TT - {(i, j)}) | {(j, i)}))
        for _ in range(200):
            a1, a2 = rng.sample(sorted(TT), 2)
            lines.append(to_string(N, (TT - {a1, a2}) | {(a1[1], a1[0]), (a2[1], a2[0])}))
        for _ in range(200):
            k = rng.randint(3, 8)
            R = rng.sample(sorted(TT), k)
            lines.append(to_string(N, (TT - set(R)) | {(b, a) for (a, b) in R}))
    run_batch("3c transitive TT_N and TT_N with 1 (all), 2, 3..8 (random) reversed arcs, N=10,11,12", lines)

    # compositions
    def relabel(arcs, off):
        return {(u + off, v + off) for (u, v) in arcs}

    C3 = circulant(3, [1])
    R5 = circulant(5, [1, 2])
    R5b = circulant(5, [1, 3])
    R7 = circulant(7, [1, 2, 4])
    R7b = circulant(7, [1, 2, 3])
    blocks = {"C3": (3, C3), "R5": (5, R5), "R7": (7, R7), "R7b": (7, R7b), "TT1": (1, set()),
              "TT2": (2, {(0, 1)}), "TT3": (3, transitive(3))}

    def hub_comp(nameA, nameB, hubs):
        """hubs h: h -> A, A => B, B -> h ; hubs among themselves transitive."""
        nA, aA = blocks[nameA]
        nB, aB = blocks[nameB]
        N = hubs + nA + nB
        arcs = set()
        for i in range(hubs):
            for j in range(i + 1, hubs):
                arcs.add((i, j))
        arcs |= relabel(aA, hubs)
        arcs |= relabel(aB, hubs + nA)
        for h in range(hubs):
            for a in range(hubs, hubs + nA):
                arcs.add((h, a))
            for b in range(hubs + nA, N):
                arcs.add((b, h))
        for a in range(hubs, hubs + nA):
            for b in range(hubs + nA, N):
                arcs.add((a, b))
        return N, arcs

    lines = []
    names = []
    for A in ("C3", "R5", "R7", "R7b", "TT2", "TT3"):
        for B in ("C3", "R5", "R7", "R7b", "TT1", "TT2", "TT3"):
            for hubs in (1, 2, 3):
                N, arcs = hub_comp(A, B, hubs)
                if 7 <= N <= 14:
                    lines.append(line_with_auts(N, arcs))
                    names.append(f"{hubs}hub {A}=>{B}")
    run_batch("3c compositions  h -> A => B -> h  (A,B in C3,R5,R7,TT1..3; 1-3 hubs; N=7..14)", lines)

    # chains of cyclic/regular blocks with one back arc set and random perturbations
    def chain(blist, back=True):
        arcs = set()
        off = 0
        offs = []
        for b in blist:
            n, a = blocks[b]
            arcs |= relabel(a, off)
            offs.append((off, n))
            off += n
        N = off
        for x in range(len(offs)):
            for y in range(x + 1, len(offs)):
                for u in range(offs[x][0], offs[x][0] + offs[x][1]):
                    for v in range(offs[y][0], offs[y][0] + offs[y][1]):
                        arcs.add((u, v))
        if back:  # reverse the arcs between the first vertex of the first and last block
            u, v = 0, N - 1
            arcs.discard((u, v))
            arcs.add((v, u))
        return N, arcs

    lines = []
    for bl in (("C3", "C3", "C3"), ("C3", "C3", "C3", "TT1"), ("C3", "R5", "C3"), ("R5", "R5"),
               ("C3", "C3", "C3", "C3"), ("R5", "C3", "TT3"), ("C3", "R7"), ("R5", "R7"), ("C3", "C3", "R5"),
               ("TT1", "C3", "C3", "C3", "TT1"), ("C3", "TT2", "C3", "TT2", "C3")):
        for back in (False, True):
            N, arcs = chain(bl, back)
            if N <= 14:
                lines.append(line_with_auts(N, arcs))
    # lexicographic products
    def lex(outer, inner):
        no, ao = outer
        ni, ai = inner
        arcs = set()
        for (x, y) in ao:
            for i in range(ni):
                for j in range(ni):
                    arcs.add((x * ni + i, y * ni + j))
        for x in range(no):
            for (i, j) in ai:
                arcs.add((x * ni + i, x * ni + j))
        return no * ni, arcs

    for o, i in ((C3, C3), (C3, transitive(3)), (transitive(3), C3), (C3, R5), (R5, C3), (C3, transitive(4)),
                 (transitive(4), C3), ({(0, 1)}, R5), ({(0, 1)}, R7), (R5, {(0, 1)}), (C3, {(0, 1), }), (R7, {(0, 1)})):
        no = 1 + max(max(a) for a in o) if o else 1
        ni = 1 + max(max(a) for a in i) if i else 1
        N, arcs = lex((no, o), (ni, i))
        if N <= 14:
            lines.append(line_with_auts(N, arcs))
    # random perturbations of the above (reverse 1-3 arcs)
    base = list(lines)
    for s in base:
        tour = s.split()[0]
        N, arcs = from_string(tour)
        for _ in range(15):
            R = rng.sample(sorted(arcs), rng.randint(1, 3))
            lines.append(to_string(N, (arcs - set(R)) | {(b, a) for (a, b) in R}))
    run_batch("3c block chains, lexicographic products (C3[C3], C3[TT3], R5[C3], TT2[R7], ...) + 1-3 random arc reversals", lines)

    # ------------------------------------------------------------------ PART 4
    section("PART 4. Gutin route: ALL minimum blocking sets, classes N <= 7")
    print("cluster = deleted pairs form a disjoint union of cliques (T - D is a multipartite tournament)")
    print("1pcf    = T - D still has a 1-path-cycle factor (bipartite matching >= N-1)")
    for N in range(3, 8):
        gl = subprocess.run(["gentourng", str(N)], capture_output=True, text=True).stdout.split()
        res = engine("minsets", gl)
        tot = len(res)
        every_cl = some_cl = no_cl = 0
        any_pcf = 0
        sum_min = sum_cl = 0
        examples_none = []
        for ln in res:
            nmin = int(field(ln, "nmin"))
            ncl = int(field(ln, "cluster"))
            npcf = int(field(ln, "with_1pcf"))
            sum_min += nmin
            sum_cl += ncl
            if ncl == nmin:
                every_cl += 1
            if ncl > 0:
                some_cl += 1
            else:
                no_cl += 1
                if len(examples_none) < 3:
                    examples_none.append(ln.split()[0] + " beta=" + field(ln, "beta"))
            if npcf > 0:
                any_pcf += 1
        print(f"N={N}: {tot} classes; minimum blocking sets in total {sum_min} (cluster-type {sum_cl});"
              f" classes with EVERY min set cluster-type: {every_cl}; with SOME: {some_cl}; with NONE: {no_cl}"
              f" (e.g. {', '.join(examples_none)});"
              f" classes having a min blocking set D with a 1-path-cycle factor in T-D: {any_pcf}")
        sys.stdout.flush()

    section("VERDICT")
    tot_eq = sum(v[0] for v in summary.values())
    tot_n = sum(v[1] for v in summary.values())
    print(f"Part 3: beta == hall on {tot_eq}/{tot_n} tournaments (N = 7..14); counterexamples: "
          f"{sum(len(v[2]) for v in summary.values())}")
    print("validation parts 1-2 all OK:", allok)
    print(f"total time {time.time() - T0:.1f}s")


if __name__ == "__main__":
    main()
