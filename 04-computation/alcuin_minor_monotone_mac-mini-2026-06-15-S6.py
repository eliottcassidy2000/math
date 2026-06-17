#!/usr/bin/env python3
"""
alcuin_minor_monotone_mac-mini-2026-06-15-S6.py

Question: Is the Alcuin number minor-monotone?  Is it subgraph-monotone?
Compare against tau = n - alpha (vertex cover number), which IS minor-monotone.

Definitions
-----------
tau(G)        = minimum vertex cover = n - alpha(G)  (alpha = max independent set; brute force).

Alcuin(G)     = the Alcuin number: minimum boat capacity b such that the whole vertex
                set can be ferried across a river, where at every moment, on each bank,
                the set of vertices NOT in the company of the ferryman must be an
                INDEPENDENT set of G ("nothing eats anything when unsupervised").
                Computed via a BFS / reachability search over states.

State model (standard Alcuin / Csorba-Hurkens-Woeginger formulation)
--------------------------------------------------------------------
We want to move all n vertices from the LEFT bank to the RIGHT bank.
A "configuration" is (L, ferryman_side) where L = set of vertices on the left bank
and ferryman_side in {0=left, 1=right}.  The right bank is V \ L.

Safety: when the ferryman is on the LEFT bank, the RIGHT bank (V\L) must be independent.
        when the ferryman is on the RIGHT bank, the LEFT bank (L) must be independent.
(The bank the ferryman is on is "supervised"; the opposite bank must be conflict-free.)

A boat trip of capacity b: the ferryman picks a subset S, |S| <= b, of vertices on his
current bank, moves them and himself to the other bank.  Both resulting banks must be
safe (the bank he just left -- now unsupervised -- must be independent; the bank he
arrives at is supervised so it's automatically fine, but for cleanliness we check the
unsupervised one).

Start: L = V, ferryman on left (side 0).  For the start to be safe the right bank is
empty -> trivially independent.
Goal:  L = {} , ferryman on right (side 1).

Alcuin(G) = min b in {1,...,n} for which the goal is reachable.  (b=n always works:
            carry everyone at once if alpha trivial... but we still need each
            intermediate bank safe; b=n always succeeds because we can move everyone
            in one trip if needed only when one bank empty -- handled by BFS.)
Known closed form (Csorba, Hurkens, Woeginger 2008/2012):
   Alcuin(G) = tau(G)            if G has a "spanning structure" allowing it, more precisely
   tau(G) <= Alcuin(G) <= tau(G) + 1
   and Alcuin(G) in {tau(G), tau(G)+1}.  We compute it directly by BFS, not by the formula,
   so the BFS is the ground truth.

Minors (single step)
--------------------
edge deletion:      remove one edge.
vertex deletion:    remove one vertex (and incident edges).
edge contraction:   pick edge uv, merge u,v into one vertex; drop resulting loops and
                    parallel edges (simple graph).

We enumerate, for each connected/all non-iso graph G up to n=7, every single-step minor H
(via each of the three operations), and test:
   tau(H)    <= tau(G)      (expected: always true)
   Alcuin(H) <= Alcuin(G)   (question)
and separately the SUBGRAPH (deletion-only) version: edge & vertex deletion only.

Graph source: nauty `geng` for ALL non-isomorphic simple graphs on n vertices.
"""

import sys
import subprocess
import itertools
from functools import lru_cache

GENG = "/opt/homebrew/bin/geng"

# ----------------------------------------------------------------------
# Graph representation: n = number of vertices (0..n-1), edges = frozenset of frozenset{u,v}
# We canonicalize via a simple but correct isomorphism canonical form (small n, fine).
# ----------------------------------------------------------------------

def g6_to_graph(line):
    """Decode graph6 string -> (n, set of edges as tuples (i<j))."""
    line = line.strip()
    data = [ord(c) - 63 for c in line]
    # number of vertices
    if data[0] <= 62:
        n = data[0]
        rest = data[1:]
    else:
        # extended (not needed for n<=7) -- handle minimal case
        if data[1] <= 62:
            n = (data[1] << 12) + (data[2] << 6) + data[3]
            rest = data[4:]
        else:
            n = (data[2] << 30) + (data[3] << 24) + (data[4] << 18) + (data[5] << 12) + (data[6] << 6) + data[7]
            rest = data[8:]
    # bit vector
    bits = []
    for x in rest:
        for k in range(5, -1, -1):
            bits.append((x >> k) & 1)
    edges = set()
    idx = 0
    for j in range(1, n):
        for i in range(j):
            if idx < len(bits) and bits[idx]:
                edges.add((i, j))
            idx += 1
    return n, edges


def gen_all_graphs(n):
    """Yield all non-isomorphic simple graphs on n vertices via geng."""
    proc = subprocess.run([GENG, str(n)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("geng failed: " + proc.stderr)
    for line in proc.stdout.splitlines():
        if line.strip():
            yield g6_to_graph(line)


# ----------------------------------------------------------------------
# Canonical form for isomorphism (brute-force over permutations; n<=7 ok)
# ----------------------------------------------------------------------

def canonical(n, edges):
    """Return a canonical (frozenset of edges) under vertex relabeling. n<=7."""
    best = None
    verts = list(range(n))
    # reduce work: only permute, n<=7 -> 5040 perms worst case
    for perm in itertools.permutations(verts):
        mapped = frozenset(
            (min(perm[u], perm[v]), max(perm[u], perm[v])) for (u, v) in edges
        )
        key = tuple(sorted(mapped))
        if best is None or key < best:
            best = key
    return (n, best)


# ----------------------------------------------------------------------
# alpha (max independent set) and tau
# ----------------------------------------------------------------------

def alpha(n, edges):
    """Max independent set size, brute force over vertex subsets (n<=7)."""
    adj = [[False] * n for _ in range(n)]
    for (u, v) in edges:
        adj[u][v] = adj[v][u] = True
    best = 0
    for r in range(n, best, -1):
        found = False
        for combo in itertools.combinations(range(n), r):
            ok = True
            for a, b in itertools.combinations(combo, 2):
                if adj[a][b]:
                    ok = False
                    break
            if ok:
                found = True
                break
        if found:
            best = r
            break
    return best


def tau(n, edges):
    return n - alpha(n, edges)


# ----------------------------------------------------------------------
# Alcuin number via state BFS
# ----------------------------------------------------------------------

def is_independent(vset, adj):
    vs = list(vset)
    for i in range(len(vs)):
        for j in range(i + 1, len(vs)):
            if adj[vs[i]][vs[j]]:
                return False
    return True


def alcuin_reachable(n, edges, b):
    """Can we ferry everyone across with boat capacity b? BFS over states.

    State = (left_bank frozenset, ferryman_side)  side 0=left, 1=right.
    Start = (full set, 0).  Goal = (empty, 1).
    Safety: the UNSUPERVISED bank (the one without ferryman) must be independent.
    """
    adj = [[False] * n for _ in range(n)]
    for (u, v) in edges:
        adj[u][v] = adj[v][u] = True

    full = frozenset(range(n))
    start = (full, 0)
    goal = (frozenset(), 1)

    def safe(left, side):
        right = full - left
        if side == 0:
            # ferryman on left; right unsupervised
            return is_independent(right, adj)
        else:
            # ferryman on right; left unsupervised
            return is_independent(left, adj)

    if not safe(full, 0):
        # start config unsafe -> impossible (right empty so always safe actually)
        pass

    from collections import deque
    seen = {start}
    dq = deque([start])
    while dq:
        left, side = dq.popleft()
        if (left, side) == goal:
            return True
        # ferryman picks subset S of his current bank, |S|<=b, moves to other side.
        if side == 0:
            bank = left  # ferryman on left, moves people from left to right
        else:
            bank = full - left  # ferryman on right, moves people from right to left
        bank_list = list(bank)
        # subsets of size 0..b   (size 0 = ferryman crosses alone)
        for k in range(0, min(b, len(bank_list)) + 1):
            for S in itertools.combinations(bank_list, k):
                Sset = set(S)
                if side == 0:
                    new_left = left - Sset
                    new_side = 1
                else:
                    new_left = left | Sset
                    new_side = 0
                if not safe(new_left, new_side):
                    continue
                ns = (new_left, new_side)
                if ns not in seen:
                    seen.add(ns)
                    dq.append(ns)
    return False


def alcuin(n, edges):
    """Minimum boat capacity b in 1..n such that everyone can cross.
    For n==0 return 0; for graphs with no edges alpha=n, tau=0, alcuin=1 (need to move at least someone)."""
    if n == 0:
        return 0
    for b in range(1, n + 1):
        if alcuin_reachable(n, edges, b):
            return b
    return n  # should always succeed at b=n


# ----------------------------------------------------------------------
# Minor operations (single step)
# ----------------------------------------------------------------------

def delete_edge(n, edges, e):
    return n, set(edges) - {e}


def delete_vertex(n, edges, v):
    # remove v, relabel remaining vertices 0..n-2
    remaining = [u for u in range(n) if u != v]
    relabel = {old: new for new, old in enumerate(remaining)}
    new_edges = set()
    for (a, b) in edges:
        if a != v and b != v:
            na, nb = relabel[a], relabel[b]
            new_edges.add((min(na, nb), max(na, nb)))
    return n - 1, new_edges


def contract_edge(n, edges, e):
    """Contract edge e=(u,v): merge v into u, drop loops & parallels (simple graph)."""
    u, v = e
    if u > v:
        u, v = v, u
    # merge v into u: every edge touching v now touches u
    new_edges = set()
    for (a, b) in edges:
        if (a, b) == (u, v):
            continue  # the contracted edge -> loop, drop
        na = u if a == v else a
        nb = u if b == v else b
        if na == nb:
            continue  # loop, drop
        new_edges.add((min(na, nb), max(na, nb)))
    # now vertex v is gone; relabel to remove the gap
    remaining = [w for w in range(n) if w != v]
    relabel = {old: new for new, old in enumerate(remaining)}
    final = set()
    for (a, b) in new_edges:
        na, nb = relabel[a], relabel[b]
        final.add((min(na, nb), max(na, nb)))
    return n - 1, final


def single_step_minors(n, edges):
    """Yield (op_name, (n', edges')) for each single-step minor."""
    edges = set(edges)
    # edge deletion
    for e in list(edges):
        yield ("delete_edge", delete_edge(n, edges, e))
    # vertex deletion
    if n >= 1:
        for v in range(n):
            yield ("delete_vertex", delete_vertex(n, edges, v))
    # edge contraction
    for e in list(edges):
        yield ("contract_edge", contract_edge(n, edges, e))


# ----------------------------------------------------------------------
# Main driver
# ----------------------------------------------------------------------

def main():
    out = []
    def log(*a):
        s = " ".join(str(x) for x in a)
        out.append(s)
        print(s, flush=True)

    MAXN = 7
    log("=" * 78)
    log("Alcuin number minor-monotonicity study  (mac-mini-2026-06-15-S6)")
    log("tau(G)=n-alpha (vertex cover);  Alcuin(G)=min boat capacity via state-BFS")
    log("=" * 78)

    # caches keyed by canonical form
    tau_cache = {}
    alc_cache = {}

    def get_tau(n, edges):
        key = canonical(n, edges)
        if key not in tau_cache:
            tau_cache[key] = tau(n, edges)
        return tau_cache[key]

    def get_alc(n, edges):
        key = canonical(n, edges)
        if key not in alc_cache:
            alc_cache[key] = alcuin(n, edges)
        return alc_cache[key]

    # global counters
    tau_violations = []
    alc_violations_minor = []          # any single-step minor
    alc_violations_deletion = []       # subgraph (deletion only)
    op_counts = {"delete_edge": 0, "delete_vertex": 0, "contract_edge": 0}
    alc_viol_by_op = {"delete_edge": 0, "delete_vertex": 0, "contract_edge": 0}

    # also track Alcuin-tau relationship
    alc_minus_tau_hist = {}

    total_graphs = 0
    for n in range(1, MAXN + 1):
        graphs = list(gen_all_graphs(n))
        log(f"\n--- n={n}: {len(graphs)} non-isomorphic graphs ---")
        n_tau_v = 0
        n_alc_v_minor = 0
        n_alc_v_del = 0
        for (gn, gedges) in graphs:
            total_graphs += 1
            tG = get_tau(gn, gedges)
            aG = get_alc(gn, gedges)
            diff = aG - tG
            alc_minus_tau_hist[diff] = alc_minus_tau_hist.get(diff, 0) + 1

            for (op, (hn, hedges)) in single_step_minors(gn, gedges):
                op_counts[op] += 1
                tH = get_tau(hn, hedges)
                aH = get_alc(hn, hedges)
                # tau monotonicity
                if tH > tG:
                    n_tau_v += 1
                    if len(tau_violations) < 20:
                        tau_violations.append((n, op, sorted(gedges), tG, tH))
                # alcuin minor monotonicity
                if aH > aG:
                    n_alc_v_minor += 1
                    alc_viol_by_op[op] += 1
                    rec = (n, op, hn, sorted(gedges), sorted(hedges), tG, aG, tH, aH)
                    alc_violations_minor.append(rec)
                    if op in ("delete_edge", "delete_vertex"):
                        n_alc_v_del += 1
                        alc_violations_deletion.append(rec)
        log(f"    tau violations (minor): {n_tau_v}")
        log(f"    Alcuin violations (any single-step minor): {n_alc_v_minor}")
        log(f"    Alcuin violations (deletion-only / subgraph): {n_alc_v_del}")

    log("\n" + "=" * 78)
    log("SUMMARY")
    log("=" * 78)
    log(f"Total graphs tested (n=1..{MAXN}): {total_graphs}")
    log(f"Single-step minor checks by op: {op_counts}")
    log("")
    log("(1) tau minor-monotone?  -> " +
        ("CONFIRMED (no violations)" if not tau_violations else f"VIOLATIONS: {len(tau_violations)}"))
    if tau_violations:
        for v in tau_violations[:10]:
            log("    tau-violation:", v)

    log("")
    log("(2) Alcuin minor-monotone (any single-step minor)?  -> " +
        ("MONOTONE (no violations)" if not alc_violations_minor else
         f"NOT MONOTONE: {len(alc_violations_minor)} violations"))

    log("")
    log("(4) Alcuin SUBGRAPH-monotone (deletion-only: delete_edge + delete_vertex)?  -> " +
        ("MONOTONE (no deletion violations)" if not alc_violations_deletion else
         f"NOT MONOTONE: {len(alc_violations_deletion)} deletion violations"))

    log("")
    log("(3) Violations by operation (which op breaks Alcuin monotonicity):")
    for op in ("delete_edge", "delete_vertex", "contract_edge"):
        log(f"    {op}: {alc_viol_by_op[op]} violations")

    # smallest counterexample for Alcuin minor-monotonicity
    if alc_violations_minor:
        # sort by n (of G), then by number of edges in G
        def keyf(rec):
            n, op, hn, gedges, hedges, tG, aG, tH, aH = rec
            return (n, len(gedges), op)
        alc_violations_minor.sort(key=keyf)
        log("")
        log("SMALLEST Alcuin minor-monotonicity counterexample(s):")
        smallest_n = alc_violations_minor[0][0]
        shown = 0
        for rec in alc_violations_minor:
            n, op, hn, gedges, hedges, tG, aG, tH, aH = rec
            if n != smallest_n:
                break
            log(f"    G: n={n}, edges={gedges}, tau={tG}, Alcuin={aG}")
            log(f"    --[{op}]--> H: n={hn}, edges={hedges}, tau={tH}, Alcuin={aH}")
            log(f"      => Alcuin(H)={aH} > Alcuin(G)={aG}  (tau ok: {tH}<={tG})")
            shown += 1
            if shown >= 8:
                break

        # smallest counterexample specifically via contraction
        contr = [r for r in alc_violations_minor if r[1] == "contract_edge"]
        if contr:
            contr.sort(key=lambda rec: (rec[0], len(rec[3])))
            r = contr[0]
            n, op, hn, gedges, hedges, tG, aG, tH, aH = r
            log("")
            log("SMALLEST CONTRACTION counterexample:")
            log(f"    G: n={n}, edges={gedges}, tau={tG}, Alcuin={aG}")
            log(f"    --[contract_edge]--> H: n={hn}, edges={hedges}, tau={tH}, Alcuin={aH}")
            log(f"      => Alcuin jumps {aG} -> {aH}")

    # smallest deletion counterexample if any
    if alc_violations_deletion:
        alc_violations_deletion.sort(key=lambda rec: (rec[0], len(rec[3])))
        r = alc_violations_deletion[0]
        n, op, hn, gedges, hedges, tG, aG, tH, aH = r
        log("")
        log("SMALLEST DELETION (subgraph) counterexample:")
        log(f"    G: n={n}, edges={gedges}, tau={tG}, Alcuin={aG}")
        log(f"    --[{op}]--> H: n={hn}, edges={hedges}, tau={tH}, Alcuin={aH}")

    log("")
    log("Alcuin - tau distribution over all graphs (should be {0,1}):")
    for d in sorted(alc_minus_tau_hist):
        log(f"    Alcuin-tau = {d}:  {alc_minus_tau_hist[d]} graphs")

    log("")
    log("=" * 78)
    log("DONE")

    # write output file
    with open("/Users/e/Documents/GitHub/math/05-knowledge/results/alcuin_minor_monotone_mac-mini-2026-06-15-S6.out", "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
