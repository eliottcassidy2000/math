# n = 9 only (191,536 classes); same checks as verify_reconstruction.py
# Exhaustive check of the legal-path reconstruction theorem (all tournaments, n <= nmax).
#   (S) structure: in a strong T, an arc u -> v on no Hamilton path has T - u non-strong with v in the TERMINAL
#       component of T - u, and T - v non-strong with u in the INITIAL component of T - v;
#   (W) witness: such an arc is the chord of a Hamilton path from u to v, and no Hamilton path runs from v to u;
#   (R) decoding: from the path set alone (used arcs; cut positions -> ordered blocks, forward between blocks;
#       an unused pair inside a block oriented from the start to the end of a block-restricted path), the
#       labelled tournament is recovered exactly -- checked on every labelled representative of every class;
#   (X) exact characterisation of missed arcs: T - u non-strong (components C_1..C_k), v in C_k, u beats all of
#       C_{k-1}, and C_k admits no split into a path ending in N^-(u) and a path starting at v.
import sys, itertools
sys.path.insert(0, '../joint_port_bisimulation_20261026')
from rooted import tournaments, ham_paths, is_strong
from ham_recovery import components

def sub(adj, keep):
    idx = {v: i for i, v in enumerate(keep)}
    return tuple(sum(1 << idx[u] for u in keep if (adj[v] >> u) & 1) for v in keep)

def comps_minus(adj, n, x):
    keep = [y for y in range(n) if y != x]
    s = sub(adj, keep)
    return [[keep[i] for i in c] for c in components(s)]

def decode(n, HP):
    used = {(a, b) for P in HP for a, b in zip(P, P[1:])}
    cuts = [m for m in range(1, n) if len({frozenset(P[:m]) for P in HP}) == 1]
    P0 = HP[0]
    bounds = [0] + cuts + [n]
    blocks = [frozenset(P0[bounds[i]:bounds[i + 1]]) for i in range(len(bounds) - 1)]
    bi = {v: i for i, b in enumerate(blocks) for v in b}
    ends = {}
    for i, b in enumerate(blocks):
        ends[i] = {(P[bounds[i]], P[bounds[i + 1] - 1]) for P in HP}
    rows = [0] * n
    for x in range(n):
        for y in range(n):
            if x == y: continue
            if (x, y) in used: rows[x] |= 1 << y; continue
            if (y, x) in used: continue
            if bi[x] != bi[y]:
                if bi[x] < bi[y]: rows[x] |= 1 << y
                continue
            fx, fy = (x, y) in ends[bi[x]], (y, x) in ends[bi[x]]
            if fx == fy: return None          # decoding rule ambiguous
            if fx: rows[x] |= 1 << y
    return tuple(rows)

def can_split(adj, Ck, v, targets):
    # partition C_k into a path alpha ending in targets (nonempty) and a path beta starting at v
    Ck = list(Ck); rest = [x for x in Ck if x != v]
    for r in range(0, len(rest)):
        for beta_tail in itertools.permutations(rest, r):
            beta = (v,) + beta_tail
            if any(not ((adj[a] >> b) & 1) for a, b in zip(beta, beta[1:])): continue
            others = [x for x in rest if x not in beta_tail]
            for alpha in itertools.permutations(others):
                if alpha[-1] not in targets: continue
                if all((adj[a] >> b) & 1 for a, b in zip(alpha, alpha[1:])): return True
    return False

nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
labelled_max = int(sys.argv[2]) if len(sys.argv) > 2 else 6
for n in range(9, nmax + 1):
    classes = strong = missed = S_ok = W_ok = X_ok = X_tot = R_ok = 0
    for adj in tournaments(n):
        classes += 1
        HP = ham_paths(adj)
        R_ok += decode(n, HP) == tuple(adj)
        if not is_strong(adj): continue
        strong += 1
        used = {(a, b) for P in HP for a, b in zip(P, P[1:])}
        ends = {(P[0], P[-1]) for P in HP}
        for u in range(n):
            cu = comps_minus(adj, n, u)
            inn = {x for x in range(n) if (adj[x] >> u) & 1}
            for v in range(n):
                if not ((adj[u] >> v) & 1): continue
                is_missed = (u, v) not in used
                # exact characterisation (X)
                pred = False
                if len(cu) >= 2 and v in cu[-1]:
                    beats_prev = all((adj[u] >> x) & 1 for x in cu[-2])
                    pred = beats_prev and not can_split(adj, cu[-1], v, inn)
                X_tot += 1; X_ok += pred == is_missed
                if not is_missed: continue
                missed += 1
                cv = comps_minus(adj, n, v)
                s_ok = len(cu) >= 2 and v in cu[-1] and len(cv) >= 2 and u in cv[0]
                S_ok += s_ok
                W_ok += (u, v) in ends and (v, u) not in ends
    print(f"n = {n}: classes {classes}, strong {strong}; missed arcs in strong classes {missed}; "
          f"(S) {S_ok}/{missed}; (W) {W_ok}/{missed}; (X) exact on {X_ok}/{X_tot} arcs; "
          f"(R) decoded exactly {R_ok}/{classes} class representatives", flush=True)

# (R) on every labelled tournament for small n: decode(HP(T)) == T
for n in range(2, labelled_max + 1):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    ok = tot = 0
    for bits in range(1 << len(pairs)):
        rows = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: rows[i] |= 1 << j
            else: rows[j] |= 1 << i
        adj = tuple(rows)
        tot += 1; ok += decode(n, ham_paths(adj)) == adj
    print(f"labelled n = {n}: decode(HP(T)) == T for {ok}/{tot}", flush=True)
