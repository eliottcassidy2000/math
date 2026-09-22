#!/usr/bin/env python3
"""Lane square_sum_hamiltonicity (session collatz-mod6-20260922, wave 6).

Q_n = square-sum graph on {1..n}: x ~ y iff x != y and x + y is a perfect
square.  Exact Hamiltonian path / cycle decisions for 1 <= n <= 40, exact
counts (A071983 / A090460 / A071984 comparison), forced-edge obstruction
proofs for n = 18..22 and 24, the n = 15 edge census at vertex 4, a
verification of Gerbicz's 25-fold blow-up (the family (71*25^m-1)/2) and
the correction of the pasted "horizons" table.

Deterministic: no timing output, no network access unless SQSUM_LIVE=1.
Output must be identical under python3 and python3 -O (no bare asserts).
"""
import os
import sys
from math import isqrt

sys.setrecursionlimit(10000)

NMAX = 40
COUNT_PATH_MAX = 36      # exact undirected Hamiltonian path counts for n <= this
COUNT_CYCLE_MAX = 40     # exact undirected Hamiltonian cycle counts for n <= this


def fail(msg):
    raise RuntimeError(msg)


def is_square(m):
    r = isqrt(m)
    return r * r == m


def build(n):
    adj = {v: set() for v in range(1, n + 1)}
    for x in range(1, n + 1):
        for y in range(x + 1, n + 1):
            if is_square(x + y):
                adj[x].add(y)
                adj[y].add(x)
    return adj


def components(adj):
    seen = set()
    comps = []
    for v in sorted(adj):
        if v in seen:
            continue
        stack = [v]
        seen.add(v)
        comp = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for w in adj[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        comps.append(sorted(comp))
    return comps


# ---------------------------------------------------------------------------
# Hamiltonian path / cycle search (bitmask backtracking with sound pruning)
# ---------------------------------------------------------------------------

class Ham:
    def __init__(self, adj):
        self.n = len(adj)
        self.adj = adj
        self.nb = {v: sum(1 << (w - 1) for w in adj[v]) for v in adj}
        self.full = (1 << self.n) - 1

    def prune_ok(self, cur, visited):
        """Sound necessary condition for the remaining path cur -> all unvisited.

        Every unvisited vertex needs >= 1 available neighbour (unvisited or cur);
        at most one unvisited vertex may have exactly 1 available neighbour
        (it must be the terminal vertex); the unvisited set plus cur must be
        connected.
        """
        unv = self.full & ~visited
        if unv == 0:
            return True
        ones = 0
        curbit = 1 << (cur - 1)
        m = unv
        while m:
            b = m & -m
            m ^= b
            v = b.bit_length()
            av = self.nb[v] & (unv | curbit)
            c = bin(av).count("1")
            if c == 0:
                return False
            if c == 1:
                ones += 1
                if ones > 1:
                    return False
        # connectivity of unv | cur
        region = unv | curbit
        frontier = curbit
        reach = curbit
        while frontier:
            b = frontier & -frontier
            frontier ^= b
            v = b.bit_length()
            new = self.nb[v] & region & ~reach
            reach |= new
            frontier |= new
        return reach == region

    def count_from(self, start, count_all, target=None, want_path=False):
        """Directed Hamiltonian paths from start.  If target is given, only
        paths ending at target adjacent-closing to start are counted (cycles).
        Returns (count, one witness path or None)."""
        n = self.n
        nb = self.nb
        full = self.full
        result = [0, None]
        path = [start]

        def rec(cur, visited):
            if visited == full:
                if target is not None and not (nb[cur] >> (target - 1)) & 1:
                    return False
                result[0] += 1
                if result[1] is None:
                    result[1] = list(path)
                return not count_all
            if not self.prune_ok(cur, visited):
                return False
            cand = nb[cur] & ~visited
            # Warnsdorff order: fewest onward options first
            lst = []
            while cand:
                b = cand & -cand
                cand ^= b
                v = b.bit_length()
                lst.append((bin(nb[v] & ~visited).count("1"), v))
            lst.sort()
            for _, v in lst:
                path.append(v)
                if rec(v, visited | (1 << (v - 1))):
                    return True
                path.pop()
            return False

        rec(start, 1 << (start - 1))
        return result[0], result[1]

    def path_exists(self):
        # try starts in increasing degree order (leaves first)
        for v in sorted(self.adj, key=lambda u: (len(self.adj[u]), u)):
            c, w = self.count_from(v, count_all=False)
            if c:
                return w
        return None

    def count_paths_undirected(self):
        tot = 0
        for v in self.adj:
            c, _ = self.count_from(v, count_all=True)
            tot += c
        if tot % 2:
            fail("odd directed path count")
        return tot // 2

    def cycle_exists(self):
        if self.n < 3:
            return None
        c, w = self.count_from(1, count_all=False, target=1)
        return w if c else None

    def count_cycles_undirected(self):
        if self.n < 3:
            return 0
        c, _ = self.count_from(1, count_all=True, target=1)
        if c % 2:
            fail("odd directed cycle count")
        return c // 2


def check_path(n, seq):
    if sorted(seq) != list(range(1, n + 1)):
        return False
    return all(is_square(seq[i] + seq[i + 1]) for i in range(len(seq) - 1))


# ---------------------------------------------------------------------------
# Forced-edge reduction (Hamiltonian PATH obstruction certificates)
# ---------------------------------------------------------------------------

def forced_reduction(n, adj0, log):
    """Iterated necessary conditions for a Hamiltonian path in the graph.

    Rules (each is a theorem about any Hamiltonian path P of a graph H):
      R0  an isolated vertex (n >= 2) => no path.
      R1  a degree-1 vertex is an endpoint of P and its edge is in P (forced).
      R2  a degree-2 vertex has both edges in P.
      R3  a vertex with 3 forced edges => contradiction (P has max degree 2).
      R4  three leaves => contradiction (P has exactly two endpoints).
      R5  if v has two forced edges, its other edges are not in P: delete them.
      R6  forced edges form vertex-disjoint paths (fragments); a non-forced edge
          joining the two ends of one fragment would close a cycle on < n
          vertices: delete it.  A forced cycle on < n vertices => contradiction.
      R7  if two leaves exist, they are the endpoints; every other vertex is
          interior, so an interior vertex's fragment-end may not be joined to
          a leaf-fragment-end ... (not needed here; not implemented).
    Returns ('LEAVES>=3', leaves) / ('FORCED_CYCLE', cycle) / ('DEG3', v) /
    ('ISOLATED', v) / ('OPEN', None), with the log of steps.
    """
    adj = {v: set(adj0[v]) for v in adj0}
    forced = set()

    def frag_ends():
        # forced edges -> fragments; returns dict end -> (other_end, vertices)
        fadj = {v: set() for v in adj}
        for (a, b) in forced:
            fadj[a].add(b)
            fadj[b].add(a)
        seen = set()
        frags = []
        for v in adj:
            if v in seen or not fadj[v]:
                continue
            # walk to an end
            comp = []
            stack = [v]
            seen.add(v)
            while stack:
                u = stack.pop()
                comp.append(u)
                for w in fadj[u]:
                    if w not in seen:
                        seen.add(w)
                        stack.append(w)
            ends = [u for u in comp if len(fadj[u]) == 1]
            if ends and all(len(fadj[u]) <= 2 for u in comp):
                # order comp as a path from its smallest end (audit fix: the DFS
                # order printed a non-path 1-8-...-4-24-12 at n=24)
                ordered = [min(ends)]
                prev = None
                while len(ordered) < len(comp):
                    nxt = [w for w in fadj[ordered[-1]] if w != prev]
                    prev = ordered[-1]
                    ordered.append(nxt[0])
                comp = ordered
                ends = [comp[0], comp[-1]]
            frags.append((comp, ends))
        return fadj, frags

    rnd = 0
    while True:
        rnd += 1
        changed = False
        leaves = [v for v in sorted(adj) if len(adj[v]) == 1]
        if len(leaves) >= 3:
            log.append(f"  round {rnd}: leaves {leaves} -> contradiction (R4)")
            return ("LEAVES>=3", leaves)
        for v in sorted(adj):
            d = len(adj[v])
            if d == 0 and n >= 2:
                log.append(f"  round {rnd}: vertex {v} isolated -> contradiction")
                return ("ISOLATED", v)
            if d in (1, 2):
                for w in adj[v]:
                    e = (min(v, w), max(v, w))
                    if e not in forced:
                        forced.add(e)
                        changed = True
                        log.append(f"  round {rnd}: deg({v})={d} forces edge {e[0]}-{e[1]}")
        fadj, frags = frag_ends()
        for v in sorted(adj):
            if len(fadj[v]) >= 3:
                log.append(f"  round {rnd}: vertex {v} has forced edges to "
                           f"{sorted(fadj[v])} -> contradiction (R3)")
                return ("DEG3", v)
        for comp, ends in frags:
            if not ends and len(comp) < n:
                # forced cycle: order it
                cyc = [comp[0]]
                prev = None
                while True:
                    nxt = [w for w in fadj[cyc[-1]] if w != prev]
                    prev = cyc[-1]
                    if nxt[0] == cyc[0]:
                        break
                    cyc.append(nxt[0])
                log.append(f"  round {rnd}: forced edges close the cycle "
                           f"{'-'.join(map(str, cyc))}-{cyc[0]} on {len(cyc)} < {n} vertices"
                           f" -> contradiction (R6)")
                return ("FORCED_CYCLE", cyc)
        # R5: delete non-forced edges at vertices with two forced edges
        for v in sorted(adj):
            if len(fadj[v]) == 2:
                for w in sorted(adj[v]):
                    if w not in fadj[v]:
                        adj[v].discard(w)
                        adj[w].discard(v)
                        changed = True
                        log.append(f"  round {rnd}: {v} saturated by forced edges "
                                   f"{sorted(fadj[v])} -> delete {v}-{w} (R5)")
        # R6: delete non-forced edges between the two ends of a fragment
        for comp, ends in frags:
            if len(ends) == 2 and len(comp) < n:
                a, b = ends
                if b in adj[a] and (min(a, b), max(a, b)) not in forced:
                    adj[a].discard(b)
                    adj[b].discard(a)
                    changed = True
                    log.append(f"  round {rnd}: fragment {'-'.join(map(str, comp))} "
                               f"has ends {a},{b}; delete {a}-{b} (R6)")
        if not changed:
            return ("OPEN", None)


def cut_search(adj, kmax=3):
    """Smallest S (|S| <= kmax) with c(Q - S) > |S| + 1, else None."""
    from itertools import combinations
    verts = sorted(adj)
    for k in range(1, kmax + 1):
        best = None
        for S in combinations(verts, k):
            sub = {v: adj[v] - set(S) for v in verts if v not in S}
            c = len(components(sub))
            if c > k + 1 and (best is None or c > best[1]):
                best = (S, c)
        if best:
            return best
    return None


# ---------------------------------------------------------------------------
# Gerbicz 25-fold blow-up (Mersenneforum post #19, 2018-01-17), verified here
# ---------------------------------------------------------------------------

GLUE = [1, ("T", -1), ("T", 1), ("R", -7), ("T", 6), ("T", -6), ("R", 0), 11,
        ("R", -5), 5, 4, 12, ("T", -12), ("T", 12), ("R", 7), ("T", -8),
        ("R", 2), ("T", -3), 9, 7, ("T", 4), ("T", -4), 10, 6, ("T", 5),
        ("R", -11), 2, ("T", -2), 8, ("T", 3), ("R", -9), ("R", 9), ("T", -10),
        ("T", 10), ("T", 11), ("R", 8), 3]


def blow_up(a):
    n = len(a)
    if n % 2 == 0 or a[0] != 1 or a[-1] != 3:
        fail("blow_up needs odd n, a[0]=1, a[-1]=3")
    out = []
    for item in GLUE:
        if isinstance(item, int):
            out.append(item)
        else:
            kind, c = item
            blk = [25 * a[i] + (c if i % 2 == 0 else -c) for i in range(n)]
            if kind == "R":
                blk = blk[::-1]
            out.extend(blk)
    return out


V35 = [1, 8, 28, 21, 4, 32, 17, 19, 6, 30, 34, 15, 10, 26, 23, 13, 12, 24, 25, 11,
       5, 20, 29, 35, 14, 2, 7, 18, 31, 33, 16, 9, 27, 22, 3]

# OEIS data fetched 2026-09-22 (curl, fmt=json); embedded so the run is offline.
A090460_DATA = [1, 1, 1, 0, 0, 0, 0, 0, 3, 0, 10, 12, 35, 52, 19, 20, 349, 361, 637,
                3678, 15237, 11875, 13306, 10964, 27223, 37054]        # offset 15
A071983_DATA = [1, 1, 1, 0, 0, 0, 0, 0, 3, 0, 10, 12, 35, 52, 19, 20, 349, 392, 669,
                4041, 17175, 12960, 14026, 11889, 29123, 39550]         # offset 15
A071984_DATA = [1, 1, 11, 57, 31, 20, 25, 50, 64]                        # offset 32
A090461_DATA = [15, 16, 17, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
                38, 39, 40]
A078107_DATA = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 20, 21, 22, 24]

OEIS_COMMENTS = {
    "A090461": [
        "Conjecture: sequence includes all integers k > 24. See A090460 for the "
        "number of essentially different solutions.",
        "It is now known that 25..299 are in the sequence, see the Numberphile 2 "
        "link. - _Jud McCranie_, Jan 11 2018",
        "Every 25 <= k <= 2^20 is in the sequence and (71*25^m-1)/2 is also in the "
        "sequence for every m, hence this sequence is infinite, see Mersenneforum "
        "link for the proof; we give Hamiltonian cycle for these k values if "
        "k >= 32. - _Robert Gerbicz_, Jan 17 2017",
        "The conjecture has been proved: every k >= 25 is in the sequence, moreover "
        "for k >= 32 there is a Hamiltonian cycle; see Mersenneforum topic for a "
        "code and deterministic algorithm to find a sequence. - _Robert Gerbicz_, "
        "Jan 21 2018",
    ],
    "A090460": [
        "For n > 31, some solutions are circular; that is, the first and last "
        "numbers also sum to a square. Note that A071983 counts each circular "
        "solution n times. This sequence counts each circular solution only once. "
        "The Mathematica program uses backtracking to find all solutions, which "
        "can be printed by removing the comment symbols.",
        "This sequence counts all essentially different square chains with "
        "reversals not counted as different: the non-circular (linear) chains "
        "counted by A398909, together with the square loops in A071984, each "
        "counted only once. Linear chains first occur for n = 15, while square "
        "loops first occur for n = 32. - _Bernard Schott_, Sep 01 2026",
    ],
    "A071984": [
        "From _Bert Dobbelaere_, Dec 28 2018: (Start)",
        "It is easy to see that no solutions for n <= 30 can exist: for each value "
        "of n <= 30 at least one number exists that can only be paired with at "
        "most one other number to form a square (e.g., 18 for n=30 can only be "
        "paired with 7). No Hamiltonian cycle can exist if the graph contains a "
        "vertex of degree less than 2.",
        "For the case n=31, the nonexistence of a Hamiltonian cycle is less "
        "trivial but can be shown by hand.",
        "(End)",
    ],
    "A078107": [
        "It has been proven that there are no more terms. See A090461 for details. "
        "- _Paolo Xausa_, May 29 2024",
    ],
}
OEIS_LINKS = {
    "A090461": "Mersenneforum, The Square-Sum problem, "
               "http://mersenneforum.org/showthread.php?p=477787 (live URL returned "
               "404 on 2026-09-22; readable via web.archive.org/web/2023/ prefix)",
    "A071983": "R. Gerbicz, Proof that there is a square chain for all n > 31, "
               "https://mersenneforum.org/showthread.php?p=477787",
}


def main():
    P = print
    P("# collatz_mod6_20260922_w6_square_sum_hamiltonicity")
    P(f"# Q_n: x~y iff x!=y and x+y is a square; n = 1..{NMAX}")
    P()

    graphs = {n: build(n) for n in range(1, NMAX + 1)}

    # ---- S1: components, degree <= 1 vertices, self-loop values ----------
    P("== S1  components and degree<=1 vertices (FINITE-EXACT)")
    P("n   #comp  components (if >1)                     deg<=1 vertices")
    comp_count = {}
    leaf_table = {}
    for n in range(1, NMAX + 1):
        adj = graphs[n]
        comps = components(adj)
        comp_count[n] = len(comps)
        low = [v for v in sorted(adj) if len(adj[v]) <= 1]
        leaf_table[n] = low
        cs = " ".join("{" + ",".join(map(str, c)) + "}" for c in comps) if len(comps) > 1 else "-"
        if len(cs) > 40:
            cs = cs[:37] + "..."
        P(f"{n:<3} {len(comps):<6} {cs:<40} {low}")
    for n in range(4, 13):
        if comp_count[n] != 3:
            fail(f"components at n={n}")
    if comp_count[13] != 2:
        fail("components at 13")
    for n in range(14, NMAX + 1):
        if comp_count[n] != 1:
            fail(f"connected at n={n}")
    P("check: 3 components for 4<=n<=12, 2 at n=13, connected for 14<=n<=40 : OK")
    P(f"deg<=1 at n=18: {leaf_table[18]}  n=19: {leaf_table[19]}  "
      f"n=20..30 all == [18]: {all(leaf_table[n] == [18] for n in range(20, 31))}  "
      f"n=31,32: {leaf_table[31]},{leaf_table[32]}")
    P("neighbours of 18 in Q_n for n<=30: " + str(sorted(graphs[30][18])) +
      "  (18+18=36 is excluded by x!=y; 31 arrives at n=31)")
    loops = [x for x in range(1, NMAX + 1) if is_square(2 * x)]
    P(f"values x<={NMAX} with 2x a square (self-loop in a graph without x!=y): {loops}")
    P("  -> the pasted Lean `Adj x y := exists k, (x+1)+(y+1)=k^2` is not irreflexive "
      "(vertex 2: 2+2=4), so SimpleGraph's loopless field cannot be proved; REFUTED as stated")
    lad = [3, 7, 11, 17]
    P(f"pasted 'Delta=4 prime ladder' {lad}: differences "
      f"{[lad[i+1]-lad[i] for i in range(3)]} -> not an arithmetic progression (REFUTED)")
    P()

    # ---- S2: Hamiltonian path / cycle existence 1..40 --------------------
    P("== S2  Hamiltonian path / cycle existence, 1 <= n <= 40 (FINITE-EXACT)")
    P("n   path  cycle  #leaves  witness (path; cycle witness for n>=32 printed separately)")
    path_yes = []
    cycle_yes = []
    witness = {}
    cyc_witness = {}
    for n in range(1, NMAX + 1):
        adj = graphs[n]
        H = Ham(adj)
        if n == 1:
            w = [1]
        else:
            w = H.path_exists()
        cw = H.cycle_exists() if n >= 3 else None
        if w is not None:
            path_yes.append(n)
            witness[n] = w
            if not check_path(n, w):
                fail(f"bad witness n={n}")
        if cw is not None:
            cycle_yes.append(n)
            cyc_witness[n] = cw
            if not (check_path(n, cw) and is_square(cw[0] + cw[-1])):
                fail(f"bad cycle witness n={n}")
        ws = ",".join(map(str, w)) if w else "-"
        if len(ws) > 70:
            ws = ws[:67] + "..."
        P(f"{n:<3} {'yes' if w else 'no ':<5} {'yes' if cw else 'no ':<6} "
          f"{len(leaf_table[n]):<8} {ws}")
    P(f"n with a Hamiltonian path : {path_yes}")
    P(f"n without (2<=n<=40)      : {[n for n in range(2, NMAX + 1) if n not in path_yes]}")
    P(f"n with a Hamiltonian cycle: {cycle_yes}")
    exp_path = [1] + A090461_DATA
    if path_yes != exp_path:
        fail(f"path set {path_yes} != {exp_path}")
    if [n for n in range(2, NMAX + 1) if n not in path_yes] != A078107_DATA[1:]:
        fail("A078107 mismatch")
    if cycle_yes != list(range(32, NMAX + 1)):
        fail("cycle set")
    P("check: path set == {1} u A090461 (A078107 lists 1 by its chain convention; "
      "Q_1 has the trivial one-vertex path), cycle set == 32..40 == support of A071984 : OK")
    P("cycle witnesses:")
    for n in range(32, NMAX + 1):
        P(f"  n={n}: {','.join(map(str, cyc_witness[n]))}")
    # session lead paths
    p15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
    p23 = [18, 7, 9, 16, 20, 5, 11, 14, 2, 23, 13, 12, 4, 21, 15, 10, 6, 19, 17, 8, 1, 3, 22]
    P(f"session lead path n=15 valid: {check_path(15, p15)}; n=23 valid: {check_path(23, p23)}")
    P()

    # ---- S3: exact counts vs OEIS ----------------------------------------
    P(f"== S3  exact counts of Hamiltonian paths (up to reversal) for n <= {COUNT_PATH_MAX}"
      f" and cycles (up to rotation/reversal) for n <= {COUNT_CYCLE_MAX} (FINITE-EXACT)")
    P("n   paths(=A071983)  cycles(=A071984)  A090460 = paths-(n-1)*cycles   OEIS A071983 A090460")
    ok = True
    counts = {}
    for n in range(15, COUNT_CYCLE_MAX + 1):
        H = Ham(graphs[n])
        cyc = H.count_cycles_undirected()
        if n <= COUNT_PATH_MAX:
            pth = H.count_paths_undirected()
            a90 = pth - (n - 1) * cyc
            o83 = A071983_DATA[n - 15]
            o90 = A090460_DATA[n - 15]
            flag = "OK" if (pth == o83 and a90 == o90) else "MISMATCH"
            ok = ok and flag == "OK"
            counts[n] = (pth, cyc, a90)
            P(f"{n:<3} {pth:<16} {cyc:<17} {a90:<30} {o83:<8} {o90:<8} {flag}")
        else:
            o84 = A071984_DATA[n - 32] if n >= 32 else 0
            flag = "OK" if cyc == o84 else "MISMATCH"
            ok = ok and flag == "OK"
            P(f"{n:<3} {'(not counted)':<16} {cyc:<17} {'':<30} {'':<8} A071984={o84:<4} {flag}")
    if not ok:
        fail("OEIS count mismatch")
    P("check: all counted terms agree with A071983 / A090460 / A071984 : OK")
    P("n=32: 392 paths up to reversal = 361 essentially different + 31 cuts of the unique cycle")
    P()

    # ---- S4: obstruction certificates for the failing n --------------------
    P("== S4  forced-edge obstruction certificates for n = 18..22, 24 (PROVED by the printed reduction)")
    P("Degree lists of the failing graphs (vertex: neighbours):")
    for n in [18, 19, 20, 21, 22, 24]:
        adj = graphs[n]
        P(f"  n={n}: " + "; ".join(f"{v}:{sorted(adj[v])}" for v in sorted(adj)))
    verdicts = {}
    for n in [18, 19, 20, 21, 22, 24]:
        log = []
        verdict, data = forced_reduction(n, graphs[n], log)
        verdicts[n] = (verdict, data)
        P(f"-- n={n}: verdict {verdict} {data}")
        for line in log:
            P(line)
    n24_forced = sum(1 for line in log if "round 1" in line and "forces edge" in line)
    P(f"  n=24: forced edges in round 1: {n24_forced}")
    for n in [19, 20, 21, 22, 24]:
        if verdicts[n][0] == "OPEN":
            fail(f"reduction inconclusive at n={n}")
    if verdicts[18][0] != "LEAVES>=3":
        fail("n=18 should be the three-leaf case")
    P("cut-set search: smallest S with c(Q_n - S) > |S| + 1, |S| <= 3")
    for n in [18, 19, 20, 21, 22, 24]:
        r = cut_search(graphs[n], 3)
        P(f"  n={n}: {('S=' + str(list(r[0])) + ', c(Q-S)=' + str(r[1])) if r else 'none with |S|<=3'}")
    P("n=23 control: leaf 18 (neighbour 7) is an endpoint; forced-edge reduction verdict:")
    log = []
    v23 = forced_reduction(23, graphs[23], log)
    P(f"  n=23: verdict {v23[0]} (a path exists: {','.join(map(str, witness[23]))})")
    # which vertices are endpoints in n=23 paths
    H23 = Ham(graphs[23])
    ends = {}
    for s in range(1, 24):
        c, _ = H23.count_from(s, count_all=True)
        if c:
            ends[s] = c
    P(f"  n=23: directed path counts by start vertex: {ends}  (every path has 18 at one end)")
    P()

    # ---- S5: n=15 edge census at vertex 4 ----------------------------------
    P("== S5  n = 15: the vertex 4 and the pasted '4 is the edge to avoid' (FINITE-EXACT)")
    adj15 = graphs[15]
    P(f"neighbours of 4 in Q_15: {sorted(adj15[4])} (4+5=9, 4+12=16; 4+21=25 needs 21>15)")
    H15 = Ham(adj15)
    tot = 0
    thru45 = thru412 = both = 0
    seen = set()
    for s in range(1, 16):
        # enumerate all directed paths from s explicitly
        paths = []

        def rec(cur, visited, path):
            if visited == H15.full:
                paths.append(list(path))
                return
            if not H15.prune_ok(cur, visited):
                return
            for v in sorted(adj15[cur]):
                if not (visited >> (v - 1)) & 1:
                    path.append(v)
                    rec(v, visited | (1 << (v - 1)), path)
                    path.pop()
        rec(s, 1 << (s - 1), [s])
        for p in paths:
            key = tuple(p) if p[0] < p[-1] else tuple(reversed(p))
            if key in seen:
                continue
            seen.add(key)
            tot += 1
            es = {(min(p[i], p[i + 1]), max(p[i], p[i + 1])) for i in range(14)}
            a = (4, 5) in es
            b = (4, 12) in es
            thru45 += a
            thru412 += b
            both += a and b
    P(f"Hamiltonian paths of Q_15 up to reversal: {tot}; through 4-5: {thru45}; "
      f"through 4-12: {thru412}; through both: {both}")
    P("  -> deg(4)=2 so BOTH edges at 4 are forced in every Hamiltonian path; "
      "'4 is the edge to avoid' is REFUTED (4 is an interior vertex of the unique path)")
    deg2_15 = [v for v in sorted(adj15) if len(adj15[v]) == 2]
    P(f"degree-2 vertices of Q_15 (all edges forced): {deg2_15}; leaves: {leaf_table[15]}")
    P()

    # ---- S6: Gerbicz 25-fold blow-up verified ------------------------------
    P("== S6  Gerbicz's 25-fold blow-up (Mersenneforum post #19, 2018-01-17): "
      "finite check + verified instances")
    # finite check of the induction step: glue sums depend only on a[0]=1, a[-1]=3, n odd
    # T(c) block: first term 25*1+c = 25+c, last term (n odd) 25*3+c = 75+c
    # R(T(c)) block: first 75+c, last 25+c
    ends = []
    for item in GLUE:
        if isinstance(item, int):
            ends.append((item, item))
        else:
            kind, c = item
            ends.append((25 + c, 75 + c) if kind == "T" else (75 + c, 25 + c))
    glue_sums = [(ends[i][1], ends[i + 1][0], ends[i][1] + ends[i + 1][0]) for i in range(len(ends) - 1)]
    bad = [g for g in glue_sums if not is_square(g[2])]
    cs = sorted(c for it in GLUE if not isinstance(it, int) for c in [it[1]])
    singles = sorted(it for it in GLUE if isinstance(it, int))
    P(f"  c-values used (each once): {cs}")
    P(f"  singleton glue integers: {singles}")
    P(f"  junction sums (last of block, first of next, sum): {[(a, b, s) for a, b, s in glue_sums]}")
    P(f"  number of junctions: {len(glue_sums)}")
    P(f"  all junction sums square: {not bad}; first element {ends[0][0]}, last element {ends[-1][1]}, "
      f"closing sum {ends[0][0] + ends[-1][1]} square: {is_square(ends[0][0] + ends[-1][1])}")
    if bad or cs != list(range(-12, 13)) or singles != list(range(1, 13)):
        fail("glue check failed")
    P("  => PROVED (finite check above): for every odd n and every Hamiltonian path a of Q_n with "
      "a(1)=1, a(n)=3, the glued word is a Hamiltonian path of Q_{25n+12} with the same "
      "properties, hence (1+3=4) a Hamiltonian cycle; iterating from n=35 gives (71*25^m-1)/2.")
    a = V35
    if not (check_path(35, a) and a[0] == 1 and a[-1] == 3):
        fail("V35")
    P(f"  base n=35 path valid, starts 1 ends 3: True; closing sum {a[0] + a[-1]}")
    for m in range(1, 4):
        a = blow_up(a)
        N = len(a)
        if N != (71 * 25 ** m - 1) // 2:
            fail("length")
        okp = check_path(N, a) and a[0] == 1 and a[-1] == 3 and is_square(a[0] + a[-1])
        if not okp:
            fail(f"blow-up m={m}")
        P(f"  m={m}: N={N}=(71*25^{m}-1)/2, Hamiltonian cycle of Q_N with a(1)=1, a(N)=3 verified: {okp}")
    P()

    # ---- S7: OEIS comments verbatim (fetched 2026-09-22) ------------------
    P("== S7  OEIS comments (CITED; fetched 2026-09-22 with curl fmt=json; reproduced under OEIS CC BY-SA)")
    for key in ["A090461", "A090460", "A071984", "A078107"]:
        P(f"-- {key}")
        for c in OEIS_COMMENTS[key]:
            P("   | " + c)
    for key in ["A090461", "A071983"]:
        P(f"-- {key} link: {OEIS_LINKS[key]}")
    P("-- authors: A090460, A090461: _T. D. Noe_, Dec 01 2003; A071983, A071984: _William Rex Marshall_, "
      "Jun 16 2002; A078107: _R. K. Guy_, Dec 06 2002; A398909: _Bernard Schott_, Aug 14 2026")
    P(f"A090461 data (first 20): {A090461_DATA}")
    P(f"A078107 data (complete): {A078107_DATA}")
    P(f"A071984 offset 32, data: {A071984_DATA}")
    if os.environ.get("SQSUM_LIVE") == "1":
        import json
        import subprocess
        for key in ["A090461"]:
            try:
                raw = subprocess.run(["curl", "-s", "-A", "Mozilla/5.0 (research)",
                                      f"https://oeis.org/search?q=id:{key}&fmt=json"],
                                     capture_output=True, text=True, timeout=60).stdout
                d = json.loads(raw)
                r = d["results"][0] if isinstance(d, dict) else d[0]
                P(f"live {key}: comments match embedded: {r.get('comment') == OEIS_COMMENTS[key]}")
            except Exception as e:  # noqa
                P(f"live {key}: unavailable ({type(e).__name__})")
    P()
    P("== S7b  Mersenneforum thread 'The Square-Sum problem' (web.archive.org copy read 2026-09-22): "
      "method summary in the lane's words (CITED, not verified beyond S6)")
    P("  post #9  (R. Gerbicz, 2018-01-12): solutions for n=15,16,17,23 and all 25<=n<=2^20=1048576, "
      "Hamiltonian cycles for 32<=n<=1048576; extension step F(i,j): from a chain of length n-1, "
      "reverse the segment i..j and insert n at position i (three new adjacent sums); 'S' = fresh "
      "search, last needed at n=6109; heuristic failure chance ~ exp(-sqrt(n)).")
    P("  post #19 (2018-01-17): the 25-fold blow-up verified in S6: blocks T(c)=25*a(k)+(-1)^(k+1)*c, "
      "c=-12..12, glued with 1..12; from n=35 gives (71*25^m-1)/2 = 887, 22187, ...")
    P("  post #21 (2018-01-17): why 25 (odd square, so c and -c give distinct residues; 9 fails to glue); "
      "single-square blow-ups cannot cover all n (density sum(k>1) 1/k^2 < 1); proposal: 'nice pairs' "
      "of chains of lengths n and n+1 with equal position parity.")
    P("  post #22 (2018-01-21): THEOREM: every n>=25 has a chain, every n>=32 a Hamiltonian cycle. "
      "Method: a nice pair (seq0 of length n, seq1 of length n+1, same-parity positions) is blown up "
      "by 49 with c=-24..24 (choice of T(c,0) or T(c,1) forced by parity) and glued with 1..24, giving "
      "nice pairs for 49*n+res, res=24..72 (a complete residue system mod 49); base table n=41..2032; "
      "smaller n by lookup; O(n log n) time, n=10000000 in about 5 seconds; all chains have a(1)=1 "
      "and a(n)=8 (n even) or 3 (n odd).")
    P("  lane status: the post #19 step is PROVED here by S6's finite check; the post #22 induction is "
      "CITED only (its glue tables live in the linked squares.c, not fetched).")
    P()
    P("== S9  inheritance and scope lines")
    P("  inherited by path, not re-derived: THM-2422 (summand closure P minus {1,4,6}; M_t = 27*2^(t-4)+1), "
      "THM-2433, THM-362; summand reflection three modes 1/4/6; collatz_mod6_20260917_synthesis.md sections 1-12; "
      "counterexample_portrait (SHEET-BLIND vs SIGN-SPECIFIC).")
    P("  the three components {1,3,6,8,10}, {2,7,9}, {4,5,11,12} of Q_12 are NOT the summand module {1,4,6}: "
      "1 and 6 share a component, 4 sits in a third; no map found beyond the cardinality 3.")
    P("  untouched paste items (SCOPE, no square-sum content): 196 = 14^2 horizon, 2363 = 17*139 clock, "
      "octonion monodromy, 0.700 bits, 'martingale decay per sheet'.")
    P()

    # ---- S8: corrected horizons table ---------------------------------------
    P("== S8  the pasted 'horizons' table corrected against S1-S5")
    sq15 = sorted({p15[i] + p15[i + 1] for i in range(14)})
    rows = [
        ("N<=13 disconnected", f"3 components for 4<=N<=12, 2 at N=13 (N=1,2,3: {comp_count[1]},{comp_count[2]},{comp_count[3]}); Q_1 is connected, so true for 2<=N<=13 only", "CORRECT for 2<=N<=13 (audit: fails at N=1)"),
        ("N=14 first unification, 3 components join", "2 components at N=13 join at N=14 (edge 14-2, 14-11); the 3rd component vanished at N=13 (13+3=16, 13+12=25)", "REFUTED as stated"),
        ("N=15 first Hamiltonian path 'braiding 9,16,25'", f"first path at N=15, unique up to reversal ({counts[15][0]}); squares used by it: {sq15}; the edge 1-3 (sum 4) exists in Q_15 but is unused", "CORRECT"),
        ("N=16,17 connected", "connected from N=14 on; paths exist (1 each)", "CORRECT (weak)"),
        ("N in [18,22] 'parity desert'", "N=18: three leaves 16,17,18; N=19,21,22: vertex 7 acquires three forced edges 2-7, 7-9, 7-18; N=20: vertex 5 acquires 4-5, 5-11, 5-20 (S4); no parity statement is involved", "REFUTED (mechanism)"),
        ("N=23 'the 23 valve'", f"N=23: leaf 18 is an endpoint with neighbour 7; {counts[23][0]} paths, all start or end at 18", "CORRECT (data), name unearned"),
        ("N=24 stalled", f"N=24: verdict {verdicts[24][0]} {verdicts[24][1]}: after the forced-edge saturations delete 1-3, 1-15, 4-5, 3-6, 2-7, 2-14 and the fragment-closing edge 4-12, the vertices 2, 4, 18 are all leaves", "CORRECT (data); obstruction is a second-round triple of leaves"),
        ("N>=25 conjectured connected", "connected from 14; Hamiltonian PATH for all N>=25 is a THEOREM (Gerbicz 2018, CITED), not a conjecture", "REFUTED (wrong notion, wrong status)"),
        ("'4 is the edge to avoid' at N=15", "4 has degree 2; both edges 4-5, 4-12 are forced and used", "REFUTED"),
    ]
    for a_, b_, c_ in rows:
        P(f"  [{c_}] {a_}\n      -> {b_}")
    P()
    P("== done")


if __name__ == "__main__":
    main()
