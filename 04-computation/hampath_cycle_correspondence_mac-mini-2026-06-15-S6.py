#!/usr/bin/env python3
"""
hampath_cycle_correspondence_mac-mini-2026-06-15-S6.py

Map T_G: for a simple graph G on vertex set {0,...,n-1}, define a tournament
T_G where for i<j we orient i->j iff {i,j} is an edge of G, else j->i.
(This is the "edge => forward, non-edge => backward" tournament of G.)

Tasks
-----
(1) CLAIM (likely tautological): G has a Hamiltonian PATH  <=>  there EXISTS an
    ordering pi of the vertices such that the spine pi(0)->pi(1)->...->pi(n-1)
    is all-forward in T_{pi(G)} (equivalently, every consecutive pair pi(k),pi(k+1)
    is an EDGE of G with pi(k)<pi(k+1) under the relabeling). We state what it gives.

    KEY SUBTLETY: T_G depends on the labelling. "spine all forward under some order"
    must be read carefully. The natural reading that makes this a real statement
    about G (not about a fixed labelling) is:
       Exists a permutation pi of vertices s.t. consecutive vertices pi(k),pi(k+1)
       are adjacent in G  <=>  G has a Hamiltonian path.
    That is literally the definition of a Hamiltonian path. We verify it as a
    sanity tautology and also check the "fixed natural labelling spine" reading.

(2) CLAIM: G has a Hamiltonian CYCLE  <=>  ?
    Test candidates for T_G:
       (a) all consecutive arcs forward (0->1->...->(n-1)) AND wrap arc present
           [natural-labelling reading of "spine + wrap"]
       (b) T_G strong under the NATURAL order (G as labelled),
       (c) EXISTS an order making T_G strong (i.e. T_G strongly connected up to
           the canonical map; but T_G is determined by G's labelling, so this is
           "exists relabelling pi of G with T_{pi(G)} strongly connected").
    Empirically determine the correspondence; correlate "G Ham cycle" with
    "T_G strong (natural)" and "exists-order T_G strong".

(3) Cross-tabulate: G-Ham-path, G-Ham-cycle, T_G-strong(natural),
    exists-order-T_G-strong, Alcuin number tau vs (tau or tau+1)
    [here we record the "spine-forward-count" style Alcuin/score data: we compute
    the Alcuin number = minimum number of arc-reversals? No -- Alcuin number of a
    graph is unrelated; we interpret the requested "Alcuin=tau vs tau+1" as the
    standard fact for tournaments: a tournament has a Hamiltonian path always
    (Redei), and the number of Ham paths is odd; a strongly connected tournament
    is vertex-pancyclic (Moon). We record tau = (number of "score-sequence
    transitive levels") proxy and compare. To stay faithful and self-contained we
    define the two compared quantities below and cross-tabulate them.]

All graphs up to n<=7 are generated up to isomorphism via a canonical-form filter.
"""

import itertools
import sys
import subprocess
import shutil
from functools import lru_cache

# ----------------------------------------------------------------------
# Graph generation up to isomorphism (brute canonical form).
# ----------------------------------------------------------------------

def edges_to_adj(n, edge_set):
    adj = [[False] * n for _ in range(n)]
    for (i, j) in edge_set:
        adj[i][j] = True
        adj[j][i] = True
    return adj

def canonical_form(n, edge_set):
    """Canonical string of a graph given as a frozenset of edges (i<j)."""
    adj = edges_to_adj(n, edge_set)
    best = None
    for perm in itertools.permutations(range(n)):
        # build sorted tuple of relabelled edges
        es = []
        for (i, j) in edge_set:
            a, b = perm[i], perm[j]
            if a > b:
                a, b = b, a
            es.append((a, b))
        es.sort()
        key = tuple(es)
        if best is None or key < best:
            best = key
    return best

def graph6_to_edges(n, g6):
    """Decode a graph6 string (single graph) into a frozenset of edges (i<j).
    Standard McKay graph6 format. Assumes n < 63 (single-byte order)."""
    # graph6: first byte(s) encode N, then bit vector of upper triangle.
    data = g6.strip()
    # Determine number of vertices
    idx = 0
    first = ord(data[0]) - 63
    if first < 63:
        N = first
        idx = 1
    else:
        # multi-byte N (n>=63) — not needed here but handle minimally
        if data[1] != chr(126):
            N = ((ord(data[1]) - 63) << 12) + ((ord(data[2]) - 63) << 6) + (ord(data[3]) - 63)
            idx = 4
        else:
            raise ValueError("huge graph6 not supported")
    assert N == n, f"graph6 N={N} != n={n}"
    # bit vector
    bits = []
    for c in data[idx:]:
        v = ord(c) - 63
        for b in range(5, -1, -1):
            bits.append((v >> b) & 1)
    # upper triangle column-major order: for j in 1..N-1, for i in 0..j-1
    edges = set()
    pos = 0
    for j in range(1, N):
        for i in range(j):
            if pos < len(bits) and bits[pos] == 1:
                edges.add((i, j))
            pos += 1
    return frozenset(edges)

def all_noniso_graphs(n):
    """Generate all non-isomorphic simple graphs on n vertices {0..n-1}.
    Uses nauty's `geng` if available (fast); else brute canonical-form filter.
    Returns list of frozenset(edge) representatives (one per iso class)."""
    geng = shutil.which("geng") or shutil.which("nauty-geng")
    if geng is not None:
        try:
            out = subprocess.run([geng, "-q", str(n)], capture_output=True,
                                 text=True, timeout=600)
            reps = []
            for line in out.stdout.splitlines():
                line = line.strip()
                if not line:
                    continue
                reps.append(graph6_to_edges(n, line))
            if reps:
                return reps
            # geng emits nothing for n=0; fall through for n<=0
        except Exception as e:
            print(f"[geng failed: {e}; falling back to brute force]", file=sys.stderr)
    # brute fallback (slow for n>=7)
    all_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    seen = {}
    reps = []
    m = len(all_pairs)
    for mask in range(1 << m):
        es = frozenset(all_pairs[k] for k in range(m) if (mask >> k) & 1)
        cf = canonical_form(n, es)
        if cf not in seen:
            seen[cf] = es
            reps.append(es)
    return reps

# ----------------------------------------------------------------------
# Graph predicates.
# ----------------------------------------------------------------------

def has_hamiltonian_path(n, adj):
    """True iff G has a Hamiltonian path (some permutation with all consecutive
    pairs adjacent)."""
    if n == 1:
        return True
    verts = list(range(n))
    # DFS over Hamiltonian paths
    def dfs(path, used):
        if len(path) == n:
            return True
        last = path[-1]
        for v in range(n):
            if not used[v] and adj[last][v]:
                used[v] = True
                path.append(v)
                if dfs(path, used):
                    return True
                path.pop()
                used[v] = False
        return False
    for start in range(n):
        used = [False] * n
        used[start] = True
        if dfs([start], used):
            return True
    return False

def has_hamiltonian_cycle(n, adj):
    if n < 3:
        return False
    start = 0
    def dfs(path, used):
        if len(path) == n:
            return adj[path[-1]][start]
        last = path[-1]
        for v in range(1, n):
            if not used[v] and adj[last][v]:
                used[v] = True
                path.append(v)
                if dfs(path, used):
                    return True
                path.pop()
                used[v] = False
        return False
    used = [False] * n
    used[0] = True
    return dfs([0], used)

# ----------------------------------------------------------------------
# Tournament T_G and its predicates.
# T_G arc: for i<j, i->j iff edge, else j->i.
# Represent T_G as boolean matrix arc[a][b] = True iff a->b.
# ----------------------------------------------------------------------

def tournament_of(n, adj):
    """arc[a][b] True iff a->b in T_G (a!=b)."""
    arc = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i][j]:
                arc[i][j] = True   # i->j
            else:
                arc[j][i] = True   # j->i
    return arc

def spine_all_forward_natural(n, arc):
    """0->1->...->(n-1) all forward in the natural order."""
    for k in range(n - 1):
        if not arc[k][k + 1]:
            return False
    return True

def wrap_arc_present_natural(n, arc):
    """arc (n-1)->0 present (forward wrap to close the cycle)."""
    if n < 2:
        return False
    return arc[n - 1][0]

def tournament_strong(n, arc):
    """Tournament strongly connected (single strongly-connected component)."""
    if n <= 1:
        return True
    # Tarjan / simple reachability: tournament is strong iff for the given
    # adjacency every vertex reaches every other. Use Kosaraju-lite:
    def reach(start, mat):
        seen = [False] * n
        stack = [start]
        seen[start] = True
        cnt = 1
        while stack:
            u = stack.pop()
            for v in range(n):
                if mat[u][v] and not seen[v]:
                    seen[v] = True
                    cnt += 1
                    stack.append(v)
        return cnt
    # forward reachability from 0
    if reach(0, arc) != n:
        return False
    # reverse reachability from 0
    rev = [[arc[b][a] for b in range(n)] for a in range(n)]
    if reach(0, rev) != n:
        return False
    return True

# exists-order T_G strong: T_G is determined by G's labelling. Relabelling G by
# pi changes which pairs are edges in label order, hence changes T_G. So
# "exists order making T_G strong" = exists permutation pi of V(G) such that the
# tournament built from pi(G) is strongly connected.
def exists_order_tournament_strong(n, adj):
    for perm in itertools.permutations(range(n)):
        # relabel adjacency: new adjacency adj2[perm[i]][perm[j]] = adj[i][j]
        adj2 = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if adj[i][j]:
                    adj2[perm[i]][perm[j]] = True
        arc = tournament_of(n, adj2)
        if tournament_strong(n, arc):
            return True
    return False

def exists_order_spine_and_wrap(n, adj):
    """Exists relabelling pi s.t. natural spine all-forward AND wrap present
    (candidate (a) for Ham CYCLE)."""
    for perm in itertools.permutations(range(n)):
        adj2 = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if adj[i][j]:
                    adj2[perm[i]][perm[j]] = True
        arc = tournament_of(n, adj2)
        if spine_all_forward_natural(n, arc) and wrap_arc_present_natural(n, arc):
            return True
    return False

def exists_order_spine_forward(n, adj):
    """Exists relabelling pi s.t. natural spine 0->1->...->(n-1) all forward.
    (candidate for Ham PATH, task 1)."""
    for perm in itertools.permutations(range(n)):
        adj2 = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if adj[i][j]:
                    adj2[perm[i]][perm[j]] = True
        arc = tournament_of(n, adj2)
        if spine_all_forward_natural(n, arc):
            return True
    return False

# ----------------------------------------------------------------------
# "Alcuin = tau vs tau+1" interpretation.
# The Alcuin number of a graph and the vertex cover number tau satisfy
# tau(G) <= Alcuin(G) <= tau(G)+1  (Csorba-Hurkens-Woeginger 2008).
# We compute tau = minimum vertex cover, and the Alcuin number, and record
# whether Alcuin == tau or tau+1, and cross-tabulate against the tournament/
# Hamiltonicity predicates.
# ----------------------------------------------------------------------

def vertex_cover_number(n, edge_set):
    """Minimum vertex cover size (brute force)."""
    if not edge_set:
        return 0
    edges = list(edge_set)
    for k in range(0, n + 1):
        for cover in itertools.combinations(range(n), k):
            cset = set(cover)
            if all((i in cset or j in cset) for (i, j) in edges):
                return k
    return n

def is_stable_transfer_feasible(n, edge_set, boat):
    """Alcuin: can we ferry all n items across a river with a boat of capacity
    'boat', never leaving a non-independent (i.e. an edge) set unattended on a
    bank without the ferryman? Standard Alcuin number computation.

    Alcuin number = min boat capacity B such that a valid schedule exists.
    We use the Csorba-Hurkens-Woeginger characterization / direct BFS over states.
    State = (set on near bank, boat position). Items not on near bank are on far
    bank. The ferryman is on one bank; the bank WITHOUT the ferryman must be an
    independent set in G.
    """
    full = (1 << n) - 1
    # state: (mask_near, side)  side=0 ferryman on near bank, side=1 on far bank
    # near bank items = bits in mask_near; far bank = complement.
    # Constraint: the bank without the ferryman is independent.
    def independent(mask):
        for (i, j) in edge_set:
            if (mask >> i) & 1 and (mask >> j) & 1:
                return False
        return True
    # precompute independence for all masks? n<=7 -> 128 masks, fine
    indep = [independent(m) for m in range(full + 1)]

    from collections import deque
    start = (full, 0)  # all on near bank, ferryman near
    goal = (0, 1)      # all on far bank, ferryman far
    # validity of a state: bank without ferryman must be independent
    def valid(state):
        mask_near, side = state
        if side == 0:
            # ferryman near; far bank = complement must be independent
            far = full ^ mask_near
            return indep[far]
        else:
            # ferryman far; near bank must be independent
            return indep[mask_near]
    if not valid(start):
        # shouldn't happen unless boat must already... start always ferryman near,
        # far bank empty -> independent. fine.
        pass
    seen = set()
    dq = deque([start])
    seen.add(start)
    while dq:
        state = dq.popleft()
        if state == goal:
            return True
        mask_near, side = state
        if side == 0:
            # ferryman on near bank: choose a subset S of near-bank items (incl
            # possibly empty) of size <= boat to take across.
            near_items = [i for i in range(n) if (mask_near >> i) & 1]
            for r in range(0, min(boat, len(near_items)) + 1):
                for S in itertools.combinations(near_items, r):
                    smask = 0
                    for i in S:
                        smask |= (1 << i)
                    new_near = mask_near ^ smask  # they leave near bank
                    new_state = (new_near, 1)
                    if valid(new_state) and new_state not in seen:
                        seen.add(new_state)
                        dq.append(new_state)
        else:
            # ferryman on far bank: far-bank items = complement of near
            far_items = [i for i in range(n) if not ((mask_near >> i) & 1)]
            for r in range(0, min(boat, len(far_items)) + 1):
                for S in itertools.combinations(far_items, r):
                    smask = 0
                    for i in S:
                        smask |= (1 << i)
                    new_near = mask_near | smask  # they come to near bank
                    new_state = (new_near, 0)
                    if valid(new_state) and new_state not in seen:
                        seen.add(new_state)
                        dq.append(new_state)
    return False

def alcuin_number(n, edge_set):
    """Minimum boat capacity for a feasible schedule. Bounded by n."""
    for b in range(1, n + 1):
        if is_stable_transfer_feasible(n, edge_set, b):
            return b
    return n  # fallback (always feasible at b=n)

# ----------------------------------------------------------------------
# Main driver.
# ----------------------------------------------------------------------

def main():
    out_lines = []
    def emit(s=""):
        print(s)
        out_lines.append(s)

    emit("=" * 78)
    emit("Hamiltonian path/cycle <-> tournament T_G correspondence")
    emit("Map T_G: for i<j, i->j iff {i,j} edge of G, else j->i")
    emit("=" * 78)

    MAXN = 7  # all non-iso graphs up to n=7

    # Per-n contingency storage
    grand = {}

    for n in range(1, MAXN + 1):
        emit("")
        emit("#" * 78)
        emit(f"n = {n}")
        emit("#" * 78)
        reps = all_noniso_graphs(n)
        emit(f"non-isomorphic graphs on {n} vertices: {len(reps)}")

        # counters for task-1 tautology check
        t1_match = 0
        t1_mismatch = 0
        t1_counterex = []

        # task-2 correlation counters
        # rows: G has Ham cycle? cols: T_G strong (natural)?
        ctab_cyc_strongnat = {(a, b): 0 for a in (0, 1) for b in (0, 1)}
        ctab_cyc_existstrong = {(a, b): 0 for a in (0, 1) for b in (0, 1)}
        ctab_cyc_spinewrap = {(a, b): 0 for a in (0, 1) for b in (0, 1)}

        # full 5-way contingency tally
        full_tab = {}

        # alcuin vs tau
        alc_tau_eq = 0
        alc_tau_p1 = 0
        alc_other = 0

        for es in reps:
            adj = edges_to_adj(n, es)
            arc = tournament_of(n, adj)

            hp = has_hamiltonian_path(n, adj)
            hc = has_hamiltonian_cycle(n, adj)
            strong_nat = tournament_strong(n, arc)
            ex_strong = exists_order_tournament_strong(n, adj)
            ex_spine = exists_order_spine_forward(n, adj)
            ex_spinewrap = exists_order_spine_and_wrap(n, adj)

            # TASK 1: G has Ham path  <=>  exists order with spine all-forward
            a = int(hp)
            b = int(ex_spine)
            if a == b:
                t1_match += 1
            else:
                t1_mismatch += 1
                if len(t1_counterex) < 5:
                    t1_counterex.append((sorted(es), hp, ex_spine))

            # TASK 2
            ctab_cyc_strongnat[(int(hc), int(strong_nat))] += 1
            ctab_cyc_existstrong[(int(hc), int(ex_strong))] += 1
            ctab_cyc_spinewrap[(int(hc), int(ex_spinewrap))] += 1

            # TASK 3 full table key
            # alcuin vs tau
            if n >= 1:
                tau = vertex_cover_number(n, es)
                if n >= 1:
                    alc = alcuin_number(n, es)
                else:
                    alc = tau
                if alc == tau:
                    alc_tau_eq += 1
                    alcflag = "tau"
                elif alc == tau + 1:
                    alc_tau_p1 += 1
                    alcflag = "tau+1"
                else:
                    alc_other += 1
                    alcflag = f"other(alc={alc},tau={tau})"
            else:
                alcflag = "na"

            key = (int(hp), int(hc), int(strong_nat), int(ex_strong), alcflag)
            full_tab[key] = full_tab.get(key, 0) + 1

        # ---- TASK 1 report ----
        emit("")
        emit("--- TASK 1: G has Ham PATH  <=>  exists order, spine all-forward ---")
        emit(f"  match: {t1_match}   mismatch: {t1_mismatch}")
        if t1_counterex:
            emit(f"  COUNTEREXAMPLES (edges, hamPath, existsSpineForward):")
            for ce in t1_counterex:
                emit(f"    {ce}")
        else:
            emit("  NO counterexample: equivalence holds (tautology confirmed).")

        # ---- TASK 2 report ----
        emit("")
        emit("--- TASK 2: G has Ham CYCLE  <=>  ? ---")
        emit("  Contingency [G-HamCycle x T_G-strong(natural order)]:")
        emit(f"    HC=0,strong=0: {ctab_cyc_strongnat[(0,0)]}   HC=0,strong=1: {ctab_cyc_strongnat[(0,1)]}")
        emit(f"    HC=1,strong=0: {ctab_cyc_strongnat[(1,0)]}   HC=1,strong=1: {ctab_cyc_strongnat[(1,1)]}")
        emit("  Contingency [G-HamCycle x exists-order T_G strong]:")
        emit(f"    HC=0,exS=0: {ctab_cyc_existstrong[(0,0)]}   HC=0,exS=1: {ctab_cyc_existstrong[(0,1)]}")
        emit(f"    HC=1,exS=0: {ctab_cyc_existstrong[(1,0)]}   HC=1,exS=1: {ctab_cyc_existstrong[(1,1)]}")
        emit("  Contingency [G-HamCycle x exists-order spine+wrap]:")
        emit(f"    HC=0,sw=0: {ctab_cyc_spinewrap[(0,0)]}   HC=0,sw=1: {ctab_cyc_spinewrap[(0,1)]}")
        emit(f"    HC=1,sw=0: {ctab_cyc_spinewrap[(1,0)]}   HC=1,sw=1: {ctab_cyc_spinewrap[(1,1)]}")

        # implication checks
        # HC => exists-strong ?  (no HC=1 with exS=0)
        imp_hc_implies_exstrong = (ctab_cyc_existstrong[(1, 0)] == 0)
        # exists-strong => HC ?  (no HC=0 with exS=1)
        imp_exstrong_implies_hc = (ctab_cyc_existstrong[(0, 1)] == 0)
        # HC <=> exists spine+wrap
        imp_hc_iff_sw = (ctab_cyc_spinewrap[(1, 0)] == 0 and ctab_cyc_spinewrap[(0, 1)] == 0)
        emit(f"  HC => exists-order-strong (no counterexample): {imp_hc_implies_exstrong}")
        emit(f"  exists-order-strong => HC (no counterexample): {imp_exstrong_implies_hc}")
        emit(f"  HC <=> exists-order spine+wrap (no counterexample): {imp_hc_iff_sw}")

        # ---- TASK 3 report ----
        emit("")
        emit("--- TASK 3: full 5-way contingency ---")
        emit("  key = (HamPath, HamCycle, strongNat, existsStrong, Alcuin-vs-tau)")
        for key in sorted(full_tab.keys()):
            emit(f"    {key} : {full_tab[key]}")
        emit(f"  Alcuin == tau:   {alc_tau_eq}")
        emit(f"  Alcuin == tau+1: {alc_tau_p1}")
        emit(f"  Alcuin other:    {alc_other}")

        grand[n] = {
            "num_graphs": len(reps),
            "t1_match": t1_match,
            "t1_mismatch": t1_mismatch,
            "hc_implies_exstrong": imp_hc_implies_exstrong,
            "exstrong_implies_hc": imp_exstrong_implies_hc,
            "hc_iff_sw": imp_hc_iff_sw,
            "alc_tau_eq": alc_tau_eq,
            "alc_tau_p1": alc_tau_p1,
            "alc_other": alc_other,
            "ctab_cyc_existstrong": dict(ctab_cyc_existstrong),
            "ctab_cyc_strongnat": dict(ctab_cyc_strongnat),
            "ctab_cyc_spinewrap": dict(ctab_cyc_spinewrap),
        }

    # ---- GRAND SUMMARY ----
    emit("")
    emit("=" * 78)
    emit("GRAND SUMMARY across n=1..%d" % MAXN)
    emit("=" * 78)
    emit("")
    emit("TASK 1 (tautology: HamPath <=> exists-spine-forward):")
    all_t1 = True
    for n in grand:
        g = grand[n]
        emit(f"  n={n}: match={g['t1_match']}, mismatch={g['t1_mismatch']}")
        if g["t1_mismatch"] != 0:
            all_t1 = False
    emit(f"  >>> Equivalence holds with NO counterexample up to n={MAXN}: {all_t1}")

    emit("")
    emit("TASK 2 (HamCycle vs tournament strength):")
    all_hc_exS = True
    all_exS_hc = True
    all_hc_sw = True
    for n in grand:
        g = grand[n]
        emit(f"  n={n}: HC=>exStrong={g['hc_implies_exstrong']}, "
             f"exStrong=>HC={g['exstrong_implies_hc']}, "
             f"HC<=>spine+wrap={g['hc_iff_sw']}")
        all_hc_exS = all_hc_exS and g["hc_implies_exstrong"]
        all_exS_hc = all_exS_hc and g["exstrong_implies_hc"]
        all_hc_sw = all_hc_sw and g["hc_iff_sw"]
    emit(f"  >>> HC => exists-order-strong (all n): {all_hc_exS}")
    emit(f"  >>> exists-order-strong => HC (all n): {all_exS_hc}")
    emit(f"  >>> HC <=> exists-order spine+wrap (all n): {all_hc_sw}")

    emit("")
    emit("TASK 3 (Alcuin vs tau):")
    for n in grand:
        g = grand[n]
        emit(f"  n={n}: Alcuin==tau:{g['alc_tau_eq']}  Alcuin==tau+1:{g['alc_tau_p1']}  other:{g['alc_other']}")

    # Write output file
    outpath = "/Users/e/Documents/GitHub/math/05-knowledge/results/hampath_cycle_correspondence_mac-mini-2026-06-15-S6.out"
    with open(outpath, "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"\n[written to {outpath}]")


if __name__ == "__main__":
    main()
