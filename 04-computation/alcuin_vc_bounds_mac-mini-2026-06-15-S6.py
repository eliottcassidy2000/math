#!/usr/bin/env python3
"""
Alcuin number vs vertex cover number bounds.

Background
----------
Alcuin's "wolf, goat, cabbage" river-crossing puzzle generalizes to an arbitrary
conflict graph G. Items are vertices; an edge {u,v} means u and v cannot be left
together on a bank without the ferryman. The ferryman shuttles across a river,
each trip carrying a subset M of items (|M| <= boat capacity b) drawn from items
on his current bank. After each move, the bank WITHOUT the ferryman must be an
INDEPENDENT set in G (no conflicting pair left unsupervised).

Alcuin(G) = minimum boat capacity b for which all items can be moved from the
left bank to the right bank.

tau(G) = vertex cover number = n - alpha(G), alpha = independence number.

Known theorem (Csorba, Hurkens, Woeginger 2008): tau(G) <= Alcuin(G) <= tau(G)+1.

This script:
  (1) Verifies tau <= Alcuin <= tau+1 for all non-iso graphs up to n=7
      (using networkx graph_atlas_g, n<=7). Reports any violations.
  (2) Reports the fraction Alcuin=tau vs Alcuin=tau+1 by n.
  (3) Tests structural characterizations of the Alcuin=tau+1 case.
"""

import sys
import itertools
from collections import deque

import networkx as nx

try:
    from networkx.generators.atlas import graph_atlas_g
    HAVE_ATLAS = True
except Exception:
    HAVE_ATLAS = False


# ---------------------------------------------------------------------------
# Graph invariants
# ---------------------------------------------------------------------------

def independence_number_and_set(adj, n):
    """Brute-force maximum independent set over vertices 0..n-1.
    adj: list of int bitmasks, adj[v] = bitmask of neighbors of v.
    Returns (alpha, one_max_independent_set_as_frozenset)."""
    best = 0
    best_set = frozenset()
    # iterate over all 2^n subsets; n<=7 so at most 128
    full = 1 << n
    for s in range(full):
        # check independence: no vertex in s has a neighbor in s
        indep = True
        ss = s
        while ss:
            v = (ss & -ss).bit_length() - 1
            if adj[v] & s:
                indep = False
                break
            ss &= ss - 1
        if indep:
            c = bin(s).count("1")
            if c > best:
                best = c
                best_set = s
    items = frozenset(i for i in range(n) if (best_set >> i) & 1)
    return best, items


def is_independent(adj, mask):
    """Is the vertex-set given by bitmask `mask` independent?"""
    ss = mask
    while ss:
        v = (ss & -ss).bit_length() - 1
        if adj[v] & mask:
            return False
        ss &= ss - 1
    return True


# ---------------------------------------------------------------------------
# Alcuin number via BFS over states for a given boat capacity
# ---------------------------------------------------------------------------

def feasible_with_capacity(adj, n, b):
    """Is the crossing feasible with boat capacity b?
    State = (left_mask, ferry_side) where ferry_side in {0=left, 1=right}.
    left_mask = bitmask of items currently on the LEFT bank.
    Start: (full_mask, 0). Goal: (0, 1).
    Move: ferryman picks subset M (0<=|M|<=b) of items on HIS current bank,
          crosses to other side carrying M. After the move, the bank WITHOUT
          the ferryman must be independent.
    """
    full = (1 << n) - 1
    start = (full, 0)        # everything left, ferry on left
    goal = (0, 1)            # everything right, ferry on right
    if start == goal:
        return True

    # Precompute, for each bank-mask, the list of subsets of size<=b.
    # We do this lazily via combinations during BFS.
    seen = set([start])
    dq = deque([start])

    while dq:
        left_mask, side = dq.popleft()
        if side == 0:
            here_mask = left_mask          # ferry on left -> items he can carry
        else:
            here_mask = full & ~left_mask  # ferry on right -> items on right
        here_items = [i for i in range(n) if (here_mask >> i) & 1]

        new_side = 1 - side
        # choose subset M of here_items with |M| in 0..b
        maxk = min(b, len(here_items))
        for k in range(0, maxk + 1):
            for combo in itertools.combinations(here_items, k):
                m = 0
                for c in combo:
                    m |= (1 << c)
                if side == 0:
                    new_left = left_mask & ~m   # M leaves the left bank
                else:
                    new_left = left_mask | m    # M arrives on the left bank
                # After move: ferry is on new_side. The bank WITHOUT ferry:
                if new_side == 0:
                    # ferry now on left -> unattended bank is RIGHT
                    unattended = full & ~new_left
                else:
                    # ferry now on right -> unattended bank is LEFT
                    unattended = new_left
                if not is_independent(adj, unattended):
                    continue
                st = (new_left, new_side)
                if st == goal:
                    return True
                if st not in seen:
                    seen.add(st)
                    dq.append(st)
    return False


def alcuin_number(adj, n):
    """Minimum boat capacity b in 1..n for which crossing is feasible.
    Returns b, or None if infeasible for all b<=n (shouldn't happen)."""
    for b in range(1, n + 1):
        if feasible_with_capacity(adj, n, b):
            return b
    return None


# ---------------------------------------------------------------------------
# Build adjacency bitmasks from a networkx graph
# ---------------------------------------------------------------------------

def adj_from_nx(G):
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    adj = [0] * n
    for u, v in G.edges():
        if u == v:
            continue
        iu, iv = idx[u], idx[v]
        adj[iu] |= (1 << iv)
        adj[iv] |= (1 << iu)
    return adj, n


# ---------------------------------------------------------------------------
# Structural tests for the Alcuin = tau + 1 case
# ---------------------------------------------------------------------------

def has_isolated_vertex(adj, n):
    return any(adj[v] == 0 for v in range(n))


def min_degree(adj, n):
    if n == 0:
        return 0
    return min(bin(adj[v]).count("1") for v in range(n))


def has_perfect_matching(G):
    """True if G has a perfect matching (n even and matching covers all)."""
    if G.number_of_nodes() % 2 != 0:
        return False
    M = nx.max_weight_matching(G, maxcardinality=True)
    return 2 * len(M) == G.number_of_nodes()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if HAVE_ATLAS:
        atlas = graph_atlas_g()
        print(f"Using networkx graph_atlas_g: {len(atlas)} graphs (n=0..7)")
    else:
        print("graph_atlas_g unavailable; would need to generate. Aborting.")
        sys.exit(1)

    MAXN = 7  # atlas covers up to 7

    # per-n tallies
    by_n = {}            # n -> dict with counts
    violations = []      # (graph_index, n, tau, alcuin)
    plus1_examples = []  # store (n, tau, alcuin, edges, properties) for Alcuin=tau+1
    eq_examples = []     # a few Alcuin=tau examples for contrast

    for gi, G in enumerate(atlas):
        n = G.number_of_nodes()
        if n < 1 or n > MAXN:
            continue
        adj, nn = adj_from_nx(G)
        # n may differ from G node count only if isolated relabeling; use nn
        n = nn
        if n == 0:
            continue

        alpha, _ = independence_number_and_set(adj, n)
        tau = n - alpha
        alc = alcuin_number(adj, n)

        d = by_n.setdefault(n, {
            "total": 0, "eq": 0, "plus1": 0, "viol": 0,
        })
        d["total"] += 1

        # bound check
        if alc is None or alc < tau or alc > tau + 1:
            d["viol"] += 1
            violations.append((gi, n, tau, alc))
            continue

        if alc == tau:
            d["eq"] += 1
            if len(eq_examples) < 5 and n >= 3:
                eq_examples.append((n, tau, alc, sorted(G.edges())))
        else:  # alc == tau + 1
            d["plus1"] += 1
            props = {
                "n": n,
                "tau": tau,
                "alpha": alpha,
                "isolated": has_isolated_vertex(adj, n),
                "min_deg": min_degree(adj, n),
                "n_eq_2tau": (n == 2 * tau),       # tau = n/2
                "perfect_matching": has_perfect_matching(G),
                "edges": sorted(G.edges()),
                "n_edges": G.number_of_edges(),
            }
            plus1_examples.append(props)

    # -----------------------------------------------------------------
    # Report (1): bounds verification
    # -----------------------------------------------------------------
    print("\n" + "=" * 70)
    print("(1) BOUND VERIFICATION  tau <= Alcuin <= tau+1")
    print("=" * 70)
    total_graphs = sum(d["total"] for d in by_n.values())
    total_viol = sum(d["viol"] for d in by_n.values())
    print(f"Checked {total_graphs} non-isomorphic graphs, n=1..{MAXN}")
    print(f"Violations of tau <= Alcuin <= tau+1: {total_viol}")
    if violations:
        for (gi, n, tau, alc) in violations:
            print(f"   VIOLATION atlas#{gi} n={n} tau={tau} Alcuin={alc}")
    else:
        print("   NONE -- bound holds for every graph checked.")

    # -----------------------------------------------------------------
    # Report (2): fraction by n
    # -----------------------------------------------------------------
    print("\n" + "=" * 70)
    print("(2) FRACTION  Alcuin=tau  vs  Alcuin=tau+1  BY n")
    print("=" * 70)
    print(f"{'n':>3} {'total':>7} {'A=tau':>7} {'A=tau+1':>9} "
          f"{'frac(=tau)':>11} {'frac(+1)':>10}")
    for n in sorted(by_n):
        d = by_n[n]
        tot = d["total"]
        eq = d["eq"]
        p1 = d["plus1"]
        fe = eq / tot if tot else 0.0
        fp = p1 / tot if tot else 0.0
        print(f"{n:>3} {tot:>7} {eq:>7} {p1:>9} {fe:>11.4f} {fp:>10.4f}")

    # -----------------------------------------------------------------
    # Report (3): structural characterization of Alcuin = tau+1
    # -----------------------------------------------------------------
    print("\n" + "=" * 70)
    print("(3) STRUCTURAL PATTERN FOR  Alcuin = tau+1")
    print("=" * 70)
    P = plus1_examples
    nP = len(P)
    print(f"Total Alcuin=tau+1 graphs (n=1..{MAXN}): {nP}")

    if nP == 0:
        print("None found.")
    else:
        # (a) tau = n/2 / perfect matching
        c_n2tau = sum(1 for p in P if p["n_eq_2tau"])
        c_pm = sum(1 for p in P if p["perfect_matching"])
        # (b) min degree / isolated vertices
        c_iso = sum(1 for p in P if p["isolated"])
        c_mindeg0 = sum(1 for p in P if p["min_deg"] == 0)
        c_mindeg_ge1 = sum(1 for p in P if p["min_deg"] >= 1)
        # (d) alpha vs tau:  n == 2*tau  <=>  alpha == tau
        c_alpha_eq_tau = sum(1 for p in P if p["alpha"] == p["tau"])

        print("\n  Property prevalence among Alcuin=tau+1 graphs:")
        print(f"    tau = n/2  (n == 2*tau)           : {c_n2tau}/{nP}")
        print(f"    has perfect matching             : {c_pm}/{nP}")
        print(f"    alpha == tau                     : {c_alpha_eq_tau}/{nP}")
        print(f"    has isolated vertex (deg 0)      : {c_iso}/{nP}")
        print(f"    min degree == 0                  : {c_mindeg0}/{nP}")
        print(f"    min degree >= 1                  : {c_mindeg_ge1}/{nP}")

        # CONVERSE test: among ALL graphs, does the candidate property
        # exactly predict Alcuin=tau+1?  Recompute over all graphs.
        print("\n  CONVERSE / EXACTNESS test of candidate characterizations:")
        # re-scan atlas for confusion matrices
        cands = {
            "tau==n/2 (perfect VC split)": lambda p: p["n_eq_2tau"],
            "has perfect matching": lambda p: p["pm"],
            "min_deg>=1 AND tau==n/2": lambda p: (p["min_deg"] >= 1 and p["n_eq_2tau"]),
        }
        # Build full record list (recompute lightweight props for all graphs)
        recs = []
        for gi, G in enumerate(atlas):
            n = G.number_of_nodes()
            if n < 1 or n > MAXN:
                continue
            adj, nn = adj_from_nx(G)
            n = nn
            if n == 0:
                continue
            alpha, _ = independence_number_and_set(adj, n)
            tau = n - alpha
            alc = alcuin_number(adj, n)
            recs.append({
                "n": n, "tau": tau, "alpha": alpha, "alc": alc,
                "n_eq_2tau": (n == 2 * tau),
                "min_deg": min_degree(adj, n),
                "pm": has_perfect_matching(G),
            })

        for name, pred in cands.items():
            tp = fp = fn = tn = 0
            for r in recs:
                want = (r["alc"] == r["tau"] + 1)
                got = pred(r)
                if want and got:
                    tp += 1
                elif (not want) and got:
                    fp += 1
                elif want and (not got):
                    fn += 1
                else:
                    tn += 1
            exact = (fp == 0 and fn == 0)
            print(f"    [{name}]")
            print(f"       TP={tp} FP={fp} FN={fn} TN={tn}  "
                  f"{'EXACT CHARACTERIZATION' if exact else 'not exact'}")

        # The CHW theorem's refined criterion: examine triangle-free /
        # "spider/matching" structure.  Empirically, test the cleanest:
        # Alcuin=tau+1  <=>  the graph has a UNIQUE minimum vertex cover?
        # We instead test the known sharp criterion: Alcuin=tau+1 iff
        # tau == n/2 AND G has no isolated vertices is NOT exact; the true
        # CHW result: Alcuin=tau+1 iff every maximum independent set is also
        # a "dominating-free" ... -> we test min-vertex-cover uniqueness.
        print("\n  Min-vertex-cover uniqueness test:")
        def num_min_vertex_covers(adj, n, tau):
            # count vertex sets of size tau that are covers
            cnt = 0
            full = (1 << n) - 1
            for s in range(1 << n):
                if bin(s).count("1") != tau:
                    continue
                comp = full & ~s
                if is_independent(adj, comp):  # complement independent => s is a cover
                    cnt += 1
            return cnt
        # tally for plus1 vs eq
        uniq_in_plus1 = 0
        for gi, G in enumerate(atlas):
            n = G.number_of_nodes()
            if n < 1 or n > MAXN:
                continue
            adj, nn = adj_from_nx(G)
            n = nn
            if n == 0:
                continue
            alpha, _ = independence_number_and_set(adj, n)
            tau = n - alpha
            alc = alcuin_number(adj, n)
            if alc == tau + 1:
                nm = num_min_vertex_covers(adj, n, tau)
                if nm == 1:
                    uniq_in_plus1 += 1
        print(f"    Alcuin=tau+1 graphs with a UNIQUE min vertex cover: "
              f"{uniq_in_plus1}/{nP}")

        # smallest +1 examples
        print("\n  SMALLEST Alcuin=tau+1 examples (by n, then edge count):")
        P_sorted = sorted(P, key=lambda p: (p["n"], p["n_edges"]))
        for p in P_sorted[:12]:
            print(f"    n={p['n']} tau={p['tau']} alpha={p['alpha']} "
                  f"min_deg={p['min_deg']} iso={p['isolated']} "
                  f"PM={p['perfect_matching']} edges={p['edges']}")

    print("\nDONE.")


if __name__ == "__main__":
    main()
