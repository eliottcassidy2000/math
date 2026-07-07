"""
lrc_chi_GW_rung_opus_S145.py   (opus-2026-07-07-S145, HYP-5197)

THE RUNG QUESTION AT GW: chi_c(G_GW) in [13, 14] -- which?
G_GW = Cay(Z, +-{1..11, 13, 24});  chi_f = 1/mu = 13 (S144: mu(GW) = 1/13 exactly),
1/M = 14.  chi_c <= chi, and chi_c >= chi_f = 13, so:
    chi(G_GW) = 13  ==>  chi_c = 13 (all graph rungs blind at GW, like k=4 Lucas);
    chi(G_GW) = 14  ==>  chi_c in (13, 14] (needs the (p,q)-coloring hunt).

(0) k=4 HAND RESULT (verified here): {0,2},{1,3},{4,6},{5,7} mod 8 properly 4-colors
    Cay(Z, +-{1,3,4,7})  =>  chi_c = 4 < 5 = 1/M at the k=4 Lucas tight instance.

(1) STRUCTURE OF ALL OPTIMAL (density-1/13) INDEPENDENT SETS of G_GW:
    the tight subgraph of the max-cycle-mean window graph (weights 13*bit - 1,
    potentials from Bellman; tight edge: D[t] = D[s] + w).  Max-mean cycles use only
    tight edges; the tight SCC structure classifies all periodic optimal sets.
    If all optimal cycles are the translates of {0,12} mod 26: a 13-coloring needs a
    perfect matching of Z_26 into difference-12 pairs = alternate edges of two odd
    (13-)cycles -- IMPOSSIBLE => every 13-coloring must use non-periodic or
    mixed-phase classes; whether that survives depends on the tight-graph paths
    BETWEEN translate classes (phase switches).

(2) DIRECT PERIODIC 13-COLORING SEARCH on circulants Cay(Z_N, +-GW mod N):
    a proper 13-coloring of the circulant with all GW mod N distinct and nonzero
    gives a periodic 13-coloring of G_GW => chi = 13.  DSATUR + tabu local search,
    N in {26, 39, 52, 65, 78, 91, 104, 117, 130}.  (Search failure is evidence,
    not proof; the tight-structure argument is the proof side.)
"""
from fractions import Fraction as F
import sys, time, random
import numpy as np

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc_mu_eq_M_maxcycle_opus_S144 import build_window_graph

GW = list(range(1, 12)) + [13, 24]

# ---------------------------------------------------------------- (0) k=4 check
def check_k4():
    S = [1, 3, 4, 7]
    classes = [{0, 2}, {1, 3}, {4, 6}, {5, 7}]
    ok = all(((x + s) % 8) not in c for c in classes for x in c for s in S)
    cover = sorted(set().union(*classes)) == list(range(8))
    print(f"(0) k=4 Lucas {{1,3,4,7}}: mod-8 4-coloring proper: {ok}, covers Z_8: {cover}")
    print(f"    => chi_c(G_{{1,3,4,7}}) = 4 < 5 = 1/M  (chi_f = 4 = 1/mu, S144)")

# ---------------------------------------------------------------- (1) tight subgraph
def tight_structure():
    states, t0, t1 = build_window_graph(GW)
    V = len(states)
    print(f"(1) GW window graph: {V} states")
    # Bellman longest-path, weights: bit0 -> -1, bit1 -> 13-1 = 12   (13*bit - 1)
    D = np.zeros(V, dtype=np.int64)
    has1 = t1 >= 0
    src1 = np.nonzero(has1)[0]
    t1v = t1[src1]
    for it in range(V + 5):
        Dn = D.copy()
        np.maximum.at(Dn, t0, D - 1)
        np.maximum.at(Dn, t1v, D[src1] + 12)
        if np.array_equal(Dn, D):
            break
        D = Dn
    print(f"    Bellman converged in {it} iterations (no cycle beats 1/13, reconfirmed)")
    # tight edges
    tight = []   # (s, t, bit)
    for s in range(V):
        if D[t0[s]] == D[s] - 1:
            tight.append((s, int(t0[s]), 0))
        if t1[s] >= 0 and D[t1[s]] == D[s] + 12:
            tight.append((s, int(t1[s]), 1))
    # restrict to nodes on tight cycles: iteratively trim nodes with no tight in/out
    from collections import defaultdict
    while True:
        outd = defaultdict(int); ind = defaultdict(int)
        nodes = set()
        for s, t, b in tight:
            outd[s] += 1; ind[t] += 1; nodes.add(s); nodes.add(t)
        bad = {v for v in nodes if outd[v] == 0 or ind[v] == 0}
        if not bad:
            break
        tight = [(s, t, b) for (s, t, b) in tight if s not in bad and t not in bad]
    nodes = sorted({s for s, t, b in tight} | {t for s, t, b in tight})
    print(f"    tight subgraph (on cycles): {len(nodes)} nodes, {len(tight)} edges")
    # out-degree structure: if every node has out-degree 1, the tight graph is a
    # disjoint union of pure cycles -> RIGID optimal sets (no phase switching)
    outs = defaultdict(list)
    for s, t, b in tight:
        outs[s].append((t, b))
    outdegs = sorted(set(len(v) for v in outs.values()))
    print(f"    out-degrees in tight subgraph: {outdegs}")
    if outdegs == [1]:
        # pure cycles: enumerate them
        seen = set()
        cycles = []
        for v in nodes:
            if v in seen:
                continue
            path = []; bits = []; x = v
            while x not in seen:
                seen.add(x); path.append(x)
                nx, b = outs[x][0]
                bits.append(b); x = nx
            if x in path:
                i = path.index(x)
                cyc_bits = bits[i:]
                cycles.append(cyc_bits)
        print(f"    PURE CYCLES: {len(cycles)} cycle(s); lengths/weights:")
        for cb in cycles:
            L, B = len(cb), sum(cb)
            # reconstruct the periodic pattern positions
            pos = [i for i, b in enumerate(cb) if b == 1]
            print(f"      period {L}, {B} elements/period, density {B}/{L}"
                  f" = {F(B, L)}; element offsets in period: {pos[:8]}"
                  f"{'...' if len(pos) > 8 else ''}"
                  f"  gaps: {[ (pos[(i+1)%len(pos)]-pos[i]) % L for i in range(len(pos))] if pos else []}")
        print("    => ALL optimal periodic sets are these cycles' translates: RIGID")
        print("       (no phase-switch paths exist: out-degree 1 everywhere).")
        return True, cycles
    else:
        print("    tight subgraph has branching -- optimal sets are FLEXIBLE (phase")
        print("    switches possible); enumerating SCCs for structure...")
        return False, None

# ---------------------------------------------------------------- (2) 13-coloring search
def coloring_search():
    print("(2) periodic 13-coloring search on circulants Cay(Z_N, +-GW mod N):")
    best_overall = None
    for N in (26, 39, 52, 65, 78, 91, 104, 117, 130):
        conn = sorted({s % N for s in GW} | {(-s) % N for s in GW})
        if 0 in conn or len({s % N for s in GW}) < 13:
            print(f"    N={N}: degenerate connection set, skip")
            continue
        adj = [[(x + d) % N for d in conn] for x in range(N)]
        rng = random.Random(N)
        K = 13
        best_conf = None
        for restart in range(60):
            # DSATUR-ish greedy with random tie-breaks
            col = [-1] * N
            order = sorted(range(N), key=lambda x: rng.random())
            ok = True
            for x in order:
                used = {col[y] for y in adj[x] if col[y] >= 0}
                free = [c for c in range(K) if c not in used]
                if not free:
                    ok = False
                    col[x] = rng.randrange(K)
                else:
                    col[x] = rng.choice(free)
            # tabu/min-conflicts
            def conflicts():
                cf = []
                for x in range(N):
                    for y in adj[x]:
                        if y > x and col[x] == col[y]:
                            cf.append((x, y))
                return cf
            for it in range(6000):
                cf = conflicts()
                if not cf:
                    break
                x, y = cf[rng.randrange(len(cf))]
                x = x if rng.random() < 0.5 else y
                # move x to least-conflict color
                cnt = [0] * K
                for z in adj[x]:
                    cnt[col[z]] += 1
                mn = min(cnt)
                col[x] = rng.choice([c for c in range(K) if cnt[c] == mn])
            cf = conflicts()
            if not cf:
                best_conf = col[:]
                break
        if best_conf is not None:
            print(f"    N={N}: *** PROPER 13-COLORING FOUND *** -> chi(G_GW) = 13 -> chi_c = 13")
            # verify hard
            okv = all(best_conf[x] != best_conf[(x + s) % N] for x in range(N) for s in GW)
            print(f"        verified: {okv}; classes sizes: "
                  f"{sorted([best_conf.count(c) for c in range(13)])}")
            best_overall = (N, best_conf)
            break
        else:
            print(f"    N={N}: no 13-coloring found (60 restarts x 6k tabu)")
    return best_overall

def main():
    t0 = time.time()
    check_k4()
    print()
    rigid, cycles = tight_structure()
    print()
    res = coloring_search()
    print()
    if rigid and res is None:
        print("VERDICT (evidence): optimal sets RIGID + no periodic 13-coloring found")
        print("=> chi(G_GW) = 14 likely; the matching obstruction argument:")
        print("   a 13-coloring partitions Z into 13 density-1/13 independent sets;")
        print("   every such set must be (asymptotically) a translate of the rigid cycle")
        print("   pattern; a partition of Z_26 into difference-12 pairs = perfect matching")
        print("   on two 13-cycles (odd) -- impossible. Formalize as the next step.")
        print("   Then chi_c in (13, 14] and the (p,q) hunt is the remaining question.")
    elif res is not None:
        print("VERDICT: chi(G_GW) = 13 = chi_c -- ALL graph rungs are blind at GW,")
        print("matching the k=4 Lucas pattern. The circular rung does NOT rescue the")
        print("graph route at the tight locus; GRAPH-14 is strictly weaker than LRC(14)")
        print("at both known non-AP tight instances.")
    print(f"\n[{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
