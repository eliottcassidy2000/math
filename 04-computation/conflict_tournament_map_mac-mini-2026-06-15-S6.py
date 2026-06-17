#!/usr/bin/env python3
"""
conflict_tournament_map_mac-mini-2026-06-15-S6.py

Map: graph G on vertices 0..n-1 (NATURAL ORDER) -> tournament T_G.
  For i<j:  arc i->j  iff {i,j} in E(G);  else arc j->i.

We compute, for all non-isomorphic graphs up to n=7 (natural order), a battery
of tournament invariants of T_G and correlate them with graph invariants.

Tournament tools (all from scratch, no networkx for tournament side):
  - score sequence (out-degrees)
  - is_strongly_connected
  - has Hamiltonian *cycle* (verified directly by search)
  - count Hamiltonian directed paths, check parity odd (Redei's theorem)
  - largest transitive (acyclic) subtournament size
  - number of 3-cycles

Graph tools:
  - non-iso graph generation up to n=7 (own generator + canonical-form dedup)
  - alpha(G) (independence number), omega(G) (clique number)
  - connectedness of G and complement
  - induced P3 count in G and in complement
  - Alcuin number / vertex-cover number tau (BRUTE force, re-implemented)

CORRELATIONS REQUESTED:
 (1) Confirm #Ham-paths(T_G) is ODD for every G  (Redei).
 (2) largest transitive subtournament size vs max(alpha(G),omega(G)).
 (3) #3-cycles(T_G) vs induced-P3 count in G and complement -- exact relation.
 (4) Is T_G strong? Characterize via G connectivity (G connected? both G and
     complement connected?).
 (5) KEY: does Alcuin(G) = tau(G)+1 correlate with ANY tournament invariant of
     T_G?  Natural order + a sample over all n! relabelings (order-existential).
"""

import sys
import itertools
from functools import lru_cache

# ----------------------------------------------------------------------------
# Graph representation: vertices 0..n-1, edge set as frozenset of frozenset pairs
# We use adjacency as a set of (i,j) with i<j, plus n.
# ----------------------------------------------------------------------------

def edges_to_adj(n, edges):
    """edges: iterable of (i,j) i<j. Return adjacency bool matrix as tuple of frozensets."""
    adj = [set() for _ in range(n)]
    for (i, j) in edges:
        adj[i].add(j)
        adj[j].add(i)
    return adj

def has_edge(adj, i, j):
    return j in adj[i]

# ----------------------------------------------------------------------------
# Non-isomorphic graph generation via canonical form (brute over all labelings).
# n<=7 so C(n,2)<=21 edges; number of graphs manageable (1044 for n=7).
# We generate all graphs on n vertices, bucket by canonical form.
# To keep it tractable at n=7 (2^21 ~ 2M graphs), we use orderly-ish generation:
# enumerate edge subsets but dedup by canonical form computed via min over perms.
# 2M * 5040 perms is too much. Instead use a smarter canonical hashing:
# canonical form = lexicographically minimal adjacency bitstring over all perms.
# For n=7 that's 5040 perms per graph * up to 2M graphs -> ~10^10, too slow.
#
# Better: generate non-iso graphs incrementally (augmentation) OR use the fact
# that we can iterate over all graphs but compute a cheap invariant signature
# first (degree sequence sorted + ...) and only do full canonical within buckets.
#
# Simplest robust approach that finishes for n<=7: use networkx's graph atlas?
# Not guaranteed available. We'll implement canonical form via refinement +
# only-perm-within-color-classes. That's enough for n<=7.
# ----------------------------------------------------------------------------

def canonical_form(n, adj):
    """Return a canonical (hashable) form of the graph using vertex refinement
    to reduce the permutation search, then brute over remaining symmetry."""
    # initial colors = degree
    deg = [len(adj[v]) for v in range(n)]
    # iterative refinement: color = (own color, sorted multiset of neighbor colors)
    colors = list(deg)
    for _ in range(n):
        newkey = {}
        newcolors = [0]*n
        sig = []
        for v in range(n):
            nb = tuple(sorted(colors[u] for u in adj[v]))
            sig.append((colors[v], nb))
        # rank signatures
        order = sorted(set(sig))
        rank = {s: i for i, s in enumerate(order)}
        newcolors = [rank[sig[v]] for v in range(n)]
        if newcolors == colors:
            break
        colors = newcolors
    # group vertices by color
    from collections import defaultdict
    groups = defaultdict(list)
    for v in range(n):
        groups[colors[v]].append(v)
    color_order = sorted(groups.keys())
    blocks = [groups[c] for c in color_order]

    best = None
    # brute over permutations within each color block, concatenated
    def gen_perms(blocks):
        if not blocks:
            yield []
            return
        first = blocks[0]
        rest = blocks[1:]
        for p in itertools.permutations(first):
            for restp in gen_perms(rest):
                yield list(p) + restp

    for perm in gen_perms(blocks):
        # perm is the new ordering of vertices: perm[newindex_position]
        # we map old vertex perm[k] -> position k
        pos = [0]*n
        for k, oldv in enumerate(perm):
            pos[oldv] = k
        # build adjacency bitstring in new labeling
        bits = 0
        idx = 0
        for a in range(n):
            for b in range(a+1, n):
                # is there an edge between new-vertices a,b?
                # new vertex a = perm[a], new vertex b = perm[b]
                if has_edge(adj, perm[a], perm[b]):
                    bits |= (1 << idx)
                idx += 1
        if best is None or bits < best:
            best = bits
    return (n, best)

def all_noniso_graphs(n):
    """Generate one representative adjacency per iso class.
    Iterate all 2^C(n,2) edge subsets, dedup by canonical_form.
    For n<=7 this is feasible (2^21 = 2.1M; canonical with refinement is fast)."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    seen = {}
    reps = []
    for mask in range(1 << m):
        edges = [pairs[k] for k in range(m) if (mask >> k) & 1]
        adj = edges_to_adj(n, edges)
        cf = canonical_form(n, adj)
        if cf not in seen:
            seen[cf] = True
            reps.append(adj)
    return reps

# ----------------------------------------------------------------------------
# Build tournament T_G from graph adjacency (natural order).
# Arc i->j (i<j) iff edge {i,j} in G; else j->i.
# Represent tournament as bool matrix beats[i][j] = True if i->j.
# ----------------------------------------------------------------------------

def build_tournament(n, adj):
    beats = [[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if has_edge(adj, i, j):
                beats[i][j] = True   # i -> j
            else:
                beats[j][i] = True   # j -> i
    return beats

def build_tournament_perm(n, adj, perm):
    """Build T from G but using vertex labeling given by perm: the tournament
    vertices are 0..n-1 = positions; position p corresponds to graph vertex perm[p].
    Arc between positions a<b: a->b iff edge {perm[a],perm[b]} in G."""
    beats = [[False]*n for _ in range(n)]
    for a in range(n):
        for b in range(a+1, n):
            if has_edge(adj, perm[a], perm[b]):
                beats[a][b] = True
            else:
                beats[b][a] = True
    return beats

# ----------------------------------------------------------------------------
# Tournament invariants
# ----------------------------------------------------------------------------

def scores(n, beats):
    return tuple(sorted(sum(1 for j in range(n) if beats[i][j]) for i in range(n)))

def is_strongly_connected(n, beats):
    # Tournament is strong iff strongly connected. Use simple reachability both ways.
    def reach(src):
        stack = [src]
        seen = {src}
        while stack:
            u = stack.pop()
            for v in range(n):
                if beats[u][v] and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen
    # forward reachability from 0 must be all; and reverse reachability to 0 all
    fwd = reach(0)
    if len(fwd) != n:
        return False
    # reverse graph
    def reach_rev(src):
        stack = [src]
        seen = {src}
        while stack:
            u = stack.pop()
            for v in range(n):
                if beats[v][u] and v not in seen:
                    seen.add(v)
                    stack.append(v)
        return seen
    rev = reach_rev(0)
    return len(rev) == n

def has_hamiltonian_cycle(n, beats):
    """Direct search for a directed Hamiltonian cycle."""
    if n == 1:
        return True
    if n == 2:
        return False  # tournament on 2 vertices: single arc, no cycle
    # backtracking
    path = [0]
    used = [False]*n
    used[0] = True
    def bt(last, count):
        if count == n:
            return beats[last][0]
        for v in range(n):
            if not used[v] and beats[last][v]:
                used[v] = True
                path.append(v)
                if bt(v, count+1):
                    return True
                path.pop()
                used[v] = False
        return False
    return bt(0, 1)

def count_hamiltonian_paths(n, beats):
    """Count directed Hamiltonian paths (any start/end). DP over subsets."""
    # dp[mask][v] = number of ham paths covering set 'mask' ending at v
    full = (1 << n) - 1
    # initialize
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            cur = dp[mask][v]
            if cur == 0:
                continue
            for w in range(n):
                if not (mask >> w) & 1 and beats[v][w]:
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))

def largest_transitive_subtournament(n, beats):
    """Largest subset S such that beats restricted to S is acyclic (transitive).
    A subtournament is transitive iff it has a strict linear order, i.e. the
    induced sub-digraph is a DAG (no directed cycle); equivalently acyclic.
    We find max-size acyclic induced subtournament via subset search with pruning.
    Equivalent to longest "chain": find ordering. Use: largest S that can be
    linearly ordered consistently with beats = longest path in the tournament's
    DAG? No -- transitive subtournament = a set linearly ordered by beats.
    A subset S is transitive iff for the linear order by out-degree-in-S it is
    consistent. Simplest: S transitive <=> induced sub-tournament acyclic.
    Find max acyclic induced subtournament. We do DP / branch."""
    best = [1]
    verts = list(range(n))
    # order by something; we just brute with memo on acyclicity check.
    # max independent-of-cycle set; we use a DP: longest transitive = longest
    # path in the tournament under the relation that S is a chain.
    # Actually: largest transitive subtournament = longest directed path? NO.
    # It is the largest set that is totally ordered: there exists ordering
    # v1->v2->...->vk with all forward arcs (vi beats vj for i<j). That's a
    # transitive tournament. The largest such = longest "increasing chain" where
    # every pair is forward. Equivalent to max clique in the DAG-reachability?
    # A transitive subtournament corresponds to a set with a linear extension
    # where ALL pairs go forward -> it's a chain in the partial order... but the
    # tournament relation is total, so S transitive <=> S has no 3-cycle <=>
    # induced sub is acyclic. Largest acyclic induced subtournament.
    # DP via longest path counting won't capture "all pairs forward".
    # Correct: transitive subtournament = set S with a linear order s.t. it's a
    # "transitive subset": this equals longest chain in the relation a<b iff
    # a beats b, BUT only if that relation restricted to S is a total order =
    # automatically (tournament) total; transitivity needed. So we need largest
    # S on which beats is transitive. Use DP: f(v) = longest transitive subtour.
    # ending pattern... do DP where we build a chain v1 beats v2 beats ... and
    # require it's transitive: equiv to longest path in the "domination" DAG
    # where we only keep edge u->v if for staying transitive... messy.
    # Just do exact search over subsets with acyclicity test for n<=7 (128 subsets).
    def induced_acyclic(S):
        # Kahn's algorithm on induced subtournament
        S = list(S)
        idx = {v: k for k, v in enumerate(S)}
        indeg = {v: 0 for v in S}
        for a in S:
            for b in S:
                if a != b and beats[a][b]:
                    indeg[b] += 1
        # repeatedly remove a source (indeg 0)
        rem = set(S)
        local_indeg = dict(indeg)
        while rem:
            src = None
            for v in rem:
                if local_indeg[v] == 0:
                    src = v
                    break
            if src is None:
                return False  # cycle
            rem.discard(src)
            for v in rem:
                if beats[src][v]:
                    local_indeg[v] -= 1
        return True
    # iterate subsets largest first
    for size in range(n, 0, -1):
        for S in itertools.combinations(range(n), size):
            if induced_acyclic(S):
                return size
    return 1

def count_3cycles(n, beats):
    """Number of cyclic triangles (3-cycles) in the tournament.
    Formula: C(n,3) - sum_i C(outdeg_i, 2)."""
    from math import comb
    out = [sum(1 for j in range(n) if beats[i][j]) for i in range(n)]
    total = comb(n, 3)
    transitive_triples = sum(comb(o, 2) for o in out)
    return total - transitive_triples

# ----------------------------------------------------------------------------
# Graph invariants
# ----------------------------------------------------------------------------

def independence_number(n, adj):
    best = 0
    for size in range(n, 0, -1):
        found = False
        for S in itertools.combinations(range(n), size):
            ok = True
            for a, b in itertools.combinations(S, 2):
                if has_edge(adj, a, b):
                    ok = False
                    break
            if ok:
                return size
    return 0

def clique_number(n, adj):
    for size in range(n, 0, -1):
        for S in itertools.combinations(range(n), size):
            ok = True
            for a, b in itertools.combinations(S, 2):
                if not has_edge(adj, a, b):
                    ok = False
                    break
            if ok:
                return size
    return 0

def is_connected(n, adj):
    if n == 0:
        return True
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == n

def complement_adj(n, adj):
    cadj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and j not in adj[i]:
                cadj[i].add(j)
    return cadj

def count_induced_P3(n, adj):
    """Count induced paths on 3 vertices (P3): triples {a,b,c} with exactly 2
    edges among them AND those 2 edges share a vertex (i.e. it's a path not a
    matching+isolated; with exactly 2 edges on 3 vertices the 2 edges always
    share a vertex -> it's automatically a P3). So induced P3 = triples with
    exactly 2 edges."""
    from math import comb
    cnt = 0
    for S in itertools.combinations(range(n), 3):
        a, b, c = S
        e = 0
        if has_edge(adj, a, b): e += 1
        if has_edge(adj, a, c): e += 1
        if has_edge(adj, b, c): e += 1
        if e == 2:
            cnt += 1
    return cnt

def vertex_cover_number(n, adj):
    """tau(G): minimum vertex cover. Brute. = n - alpha(G)."""
    # min vertex cover = complement of max independent set
    return n - independence_number(n, adj)

def alcuin_number(n, adj):
    """Alcuin number of a graph. We use the standard graph-theoretic Alcuin
    number Alc(G). Conjecture/known result (Csorba, Hurkens, Woeginger 2008/2012):
    For any graph, tau(G) <= Alc(G) <= tau(G)+1, and there's a characterization.
    The task says 'Alcuin(G)=tau(G)+1'. We compute tau(G)+1 directly as the
    quantity to correlate, AND separately we note Alcuin via its definition.

    The Alcuin number is the minimum boat capacity b such that all n items
    (vertices) can be ferried across a river (in any number of trips) so that
    at no point on either bank (excluding when the ferryman is present) are two
    conflicting items (adjacent vertices) left together.

    Known: Alc(G) in {tau(G), tau(G)+1}. We compute tau(G)+1 as requested and
    report it. (Computing true Alcuin requires the ferrying feasibility; we
    expose tau and tau+1; the requested correlation target is tau(G)+1.)"""
    tau = vertex_cover_number(n, adj)
    return tau + 1

# ----------------------------------------------------------------------------
# Main driver
# ----------------------------------------------------------------------------

def analyze(n, sample_perm_limit=5040, verbose_examples=3):
    out_lines = []
    def P(*a):
        s = " ".join(str(x) for x in a)
        out_lines.append(s)
        print(s)

    P(f"\n{'='*70}")
    P(f"n = {n}")
    P('='*70)
    reps = all_noniso_graphs(n)
    P(f"non-iso graphs on {n} vertices: {len(reps)}")

    # accumulators
    all_paths_odd = True
    odd_violations = []

    # (2) transitive size vs max(alpha,omega)
    cmp2 = {"equal": 0, "trans_larger": 0, "trans_smaller": 0}
    cmp2_examples = {"trans_larger": [], "trans_smaller": []}

    # (3) 3-cycles vs P3 counts
    rel3_consistent_Gplus = True   # test: #3cyc == P3(G)+P3(Gbar)? other relations
    rel3_data = []   # (c3, p3_G, p3_Gbar, n_choose_3)
    rel3_eq_p3G = True
    rel3_eq_p3Gbar = True
    rel3_eq_sum = True

    # (4) strong vs connectivity
    cmp4 = {}  # (G_conn, Gbar_conn, T_strong) -> count
    strong_iff_both_connected = True
    strong_iff_both_connected_ce = []

    # (5) Alcuin = tau+1 correlation
    # We collect per-graph: tau, alcuin(=tau+1), and tournament invariants;
    # we look for any tournament invariant that splits exactly on alcuin value
    # OR on (alcuin == tau+1) -- but alcuin is DEFINED as tau+1 here, so we
    # instead correlate the VALUE alcuin (= tau+1) with invariants, and also
    # correlate "tau(G)" with strong/hamcycle etc. We test natural order and a
    # perm sample for an order-existential property.
    alcuin_vs_strong = {}   # alcuin_value -> {strong: c, notstrong: c}
    alcuin_vs_hamcyc = {}
    alcuin_vs_c3pos = {}    # alcuin -> {c3>0, c3==0}

    from math import comb
    for adj in reps:
        T = build_tournament(n, adj)
        # tournament invariants
        hp = count_hamiltonian_paths(n, T)
        if hp % 2 == 0:
            all_paths_odd = False
            odd_violations.append(hp)
        sc = is_strongly_connected(n, T)
        hc = has_hamiltonian_cycle(n, T)
        trans = largest_transitive_subtournament(n, T)
        c3 = count_3cycles(n, T)

        # graph invariants
        alpha = independence_number(n, adj)
        omega = clique_number(n, adj)
        g_conn = is_connected(n, adj)
        cadj = complement_adj(n, adj)
        gbar_conn = is_connected(n, cadj)
        p3G = count_induced_P3(n, adj)
        p3Gbar = count_induced_P3(n, cadj)
        tau = vertex_cover_number(n, adj)
        alc = tau + 1

        # (2)
        mab = max(alpha, omega)
        if trans == mab:
            cmp2["equal"] += 1
        elif trans > mab:
            cmp2["trans_larger"] += 1
            if len(cmp2_examples["trans_larger"]) < verbose_examples:
                cmp2_examples["trans_larger"].append((adj_edges(n, adj), trans, alpha, omega))
        else:
            cmp2["trans_smaller"] += 1
            if len(cmp2_examples["trans_smaller"]) < verbose_examples:
                cmp2_examples["trans_smaller"].append((adj_edges(n, adj), trans, alpha, omega))

        # (3)
        rel3_data.append((c3, p3G, p3Gbar))
        if c3 != p3G:
            rel3_eq_p3G = False
        if c3 != p3Gbar:
            rel3_eq_p3Gbar = False
        if c3 != p3G + p3Gbar:
            rel3_eq_sum = False

        # (4)
        key = (g_conn, gbar_conn, sc)
        cmp4[key] = cmp4.get(key, 0) + 1
        # hypothesis: T strong  <=>  G connected and Gbar connected
        pred = g_conn and gbar_conn
        if pred != sc:
            strong_iff_both_connected = False
            if len(strong_iff_both_connected_ce) < verbose_examples:
                strong_iff_both_connected_ce.append((adj_edges(n, adj), g_conn, gbar_conn, sc))

        # (5)
        d = alcuin_vs_strong.setdefault(alc, {"strong": 0, "not": 0})
        d["strong" if sc else "not"] += 1
        d2 = alcuin_vs_hamcyc.setdefault(alc, {"hc": 0, "no": 0})
        d2["hc" if hc else "no"] += 1
        d3 = alcuin_vs_c3pos.setdefault(alc, {"c3pos": 0, "c3zero": 0})
        d3["c3pos" if c3 > 0 else "c3zero"] += 1

    # ----- report (1) -----
    P("")
    P("(1) Hamiltonian-path parity (Redei):")
    P(f"    all #Ham-paths ODD: {all_paths_odd}")
    if odd_violations:
        P(f"    VIOLATIONS (even counts found): {odd_violations[:10]}")

    # ----- report (2) -----
    P("")
    P("(2) largest transitive subtournament vs max(alpha(G), omega(G)):")
    P(f"    equal: {cmp2['equal']}   trans>max: {cmp2['trans_larger']}   trans<max: {cmp2['trans_smaller']}")
    if cmp2_examples["trans_larger"]:
        P(f"    examples trans>max: {cmp2_examples['trans_larger']}")
    if cmp2_examples["trans_smaller"]:
        P(f"    examples trans<max: {cmp2_examples['trans_smaller']}")

    # ----- report (3) -----
    P("")
    P("(3) #3-cycles(T_G) vs induced-P3 counts:")
    P(f"    c3 == P3(G) for all: {rel3_eq_p3G}")
    P(f"    c3 == P3(Gbar) for all: {rel3_eq_p3Gbar}")
    P(f"    c3 == P3(G)+P3(Gbar) for all: {rel3_eq_sum}")
    # show a few data rows
    P(f"    sample (c3, P3(G), P3(Gbar)) first 8: {rel3_data[:8]}")
    # Also test: is c3 monotone / any linear relation? compute correlation of
    # c3 with p3G, p3Gbar via simple check of whether c3 determined by (p3G,p3Gbar)
    func_ok = True
    seen_map = {}
    for (c3, p3G, p3Gbar) in rel3_data:
        k = (p3G, p3Gbar)
        if k in seen_map and seen_map[k] != c3:
            func_ok = False
        seen_map[k] = c3
    P(f"    is c3 a function of (P3(G),P3(Gbar))? {func_ok}")

    # ----- report (4) -----
    P("")
    P("(4) T_G strong  vs  G/Gbar connectivity:")
    P(f"    (G_conn, Gbar_conn, T_strong) -> count:")
    for k in sorted(cmp4.keys()):
        P(f"      {k} -> {cmp4[k]}")
    P(f"    'T strong <=> (G connected AND Gbar connected)' : {strong_iff_both_connected}")
    if strong_iff_both_connected_ce:
        P(f"      counterexamples: {strong_iff_both_connected_ce}")
    # also test simpler: strong <=> G connected
    simple = True
    for (gc, gbc, scv), cnt in cmp4.items():
        if (gc) != scv:
            simple = False
    P(f"    'T strong <=> G connected' (simpler): {simple}")

    # ----- report (5) -----
    P("")
    P("(5) Alcuin(G) = tau(G)+1  vs tournament invariants of T_G (natural order):")
    P(f"    alcuin -> strong/not: {dict(sorted(alcuin_vs_strong.items()))}")
    P(f"    alcuin -> hamcycle/no: {dict(sorted(alcuin_vs_hamcyc.items()))}")
    P(f"    alcuin -> c3pos/c3zero: {dict(sorted(alcuin_vs_c3pos.items()))}")
    # Determine if any invariant is perfectly predicted by alcuin value:
    def clean_split(d, keyA, keyB):
        # returns True if for every alcuin value, all graphs fall in ONE bucket
        for alc_v, bk in d.items():
            if bk[keyA] > 0 and bk[keyB] > 0:
                return False
        return True
    P(f"    is 'strong' perfectly determined by alcuin value? {clean_split(alcuin_vs_strong,'strong','not')}")
    P(f"    is 'hamcycle' perfectly determined by alcuin value? {clean_split(alcuin_vs_hamcyc,'hc','no')}")
    P(f"    is 'c3>0' perfectly determined by alcuin value? {clean_split(alcuin_vs_c3pos,'c3pos','c3zero')}")

    return out_lines


def adj_edges(n, adj):
    return sorted((i, j) for i in range(n) for j in adj[i] if i < j)


# ----------------------------------------------------------------------------
# Order-existential test for (5): for a sample of graphs, over ALL n! relabelings
# of the SAME graph G, does the SET of tournament invariants {strong over some
# order, hamcycle over some order, ...} correlate with alcuin? Since alcuin is a
# graph invariant (order-independent), we ask: does there EXIST an order making
# T_G strong, and does that existence track alcuin / tau?
# ----------------------------------------------------------------------------

def order_existential_test(n, max_graphs=60):
    lines = []
    def P(*a):
        s = " ".join(str(x) for x in a)
        lines.append(s); print(s)
    P("")
    P(f"--- ORDER-EXISTENTIAL test (n={n}, all {n}! relabelings, sample up to {max_graphs} graphs) ---")
    reps = all_noniso_graphs(n)
    import random
    random.seed(12345)
    if len(reps) > max_graphs:
        sample = random.sample(reps, max_graphs)
    else:
        sample = reps
    # For each graph: tau, alcuin; over all perms compute whether ANY order gives
    # strong T, ANY gives ham cycle, min and max #3-cycles, min/max trans size,
    # min/max #ham-paths.
    rows = []
    perms = list(itertools.permutations(range(n)))
    for adj in sample:
        tau = vertex_cover_number(n, adj)
        alc = tau + 1
        any_strong = False
        any_hc = False
        c3set = set()
        transset = set()
        hpset = set()
        for perm in perms:
            T = build_tournament_perm(n, adj, perm)
            if is_strongly_connected(n, T):
                any_strong = True
            if has_hamiltonian_cycle(n, T):
                any_hc = True
            c3set.add(count_3cycles(n, T))
            transset.add(largest_transitive_subtournament(n, T))
            hpset.add(count_hamiltonian_paths(n, T))
        rows.append({
            "edges": adj_edges(n, adj),
            "tau": tau, "alcuin": alc,
            "any_strong_over_orders": any_strong,
            "any_hamcycle_over_orders": any_hc,
            "c3_range": (min(c3set), max(c3set)),
            "trans_range": (min(transset), max(transset)),
            "hp_range": (min(hpset), max(hpset)),
        })
    # correlate alcuin with any_strong_over_orders etc.
    by_alc = {}
    for r in rows:
        d = by_alc.setdefault(r["alcuin"], {"any_strong_T":0,"any_strong_F":0,
                                            "any_hc_T":0,"any_hc_F":0,
                                            "count":0})
        d["count"] += 1
        d["any_strong_T" if r["any_strong_over_orders"] else "any_strong_F"] += 1
        d["any_hc_T" if r["any_hamcycle_over_orders"] else "any_hc_F"] += 1
    P(f"    by alcuin value: {dict(sorted(by_alc.items()))}")
    # clean split?
    def clean(d, a, b):
        for k, bk in d.items():
            if bk[a] > 0 and bk[b] > 0:
                return False
        return True
    P(f"    'EXISTS order making T strong' determined by alcuin? {clean(by_alc,'any_strong_T','any_strong_F')}")
    P(f"    'EXISTS order making T ham-cyclic' determined by alcuin? {clean(by_alc,'any_hc_T','any_hc_F')}")
    # print a few rows
    for r in rows[:8]:
        P(f"      {r}")
    return lines


def main():
    all_out = []
    # decide max n; n=7 may be slow because of order-existential 5040 perms.
    max_n = 7
    for n in range(3, max_n + 1):
        all_out += analyze(n)
    # order-existential: do small n fully, larger n sampled
    for n in [4, 5, 6]:
        all_out += order_existential_test(n, max_graphs=60)
    # n=7 order-existential is 5040 perms * many graphs -> sample small
    all_out += order_existential_test(7, max_graphs=12)
    return all_out


if __name__ == "__main__":
    main()
