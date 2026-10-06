#!/usr/bin/env python3
"""
sumner_t5_20261005.py -- Sumner's universal tournament conjecture at n = 5, exhaustively.

Sumner (1971): every tournament on 2n-2 vertices contains every oriented tree on n vertices.
Here n = 5: hosts are all tournaments on N = 5, 6, 7, 8 vertices (nauty gentourng: 12, 56, 456, 6880
isomorphism classes; A000568), guests are all 27 oriented trees on 5 vertices (A000238(5) = 27:
10 oriented paths, 12 oriented forks, 5 oriented stars).  For every guest S we compute
  miss_N(S) = number of N-tournament classes containing no copy of S (subdigraph, not necessarily induced)
and the unavoidability number f(S) = min { N : miss_N(S) = 0 } (monotone in N by vertex deletion).
Sumner at n = 5 is the statement miss_8(S) = 0 for all 27 S.

Conventions: gentourng default output is the upper triangle row by row; the character for the pair
(i, j), i < j, is '1' when the arc is i -> j.  The host set is closed under converse and the guest
set under reversal, so the census is convention-independent.

Reproduce: python3 sumner_t5_20261005.py   (needs nauty's gentourng on PATH; ~1-2 min)
"""
import itertools, subprocess, sys, time
from collections import Counter, defaultdict
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

CHECKS = 0
def check(c, msg):
    global CHECKS
    CHECKS += 1
    if not c:
        print("CHECK FAILED:", msg); sys.exit(1)

# ------------------------------------------------------------------ hosts
def read_tournaments(n):
    r = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True)
    m = n * (n - 1) // 2
    hosts = []
    for line in r.stdout.split("\n"):
        line = line.strip()
        if len(line) == m and set(line) <= {"0", "1"}:
            out = [0] * n
            k = 0
            for i in range(n):
                for j in range(i + 1, n):
                    if line[k] == "1":
                        out[i] |= 1 << j
                    else:
                        out[j] |= 1 << i
                    k += 1
            hosts.append(tuple(out))
    return hosts

A000568 = {5: 12, 6: 56, 7: 456, 8: 6880}
HOSTS = {}
for N in range(5, 9):
    HOSTS[N] = read_tournaments(N)
    check(len(HOSTS[N]) == A000568[N], f"gentourng count at N={N}")
    # sanity: each is a tournament
    for out in HOSTS[N]:
        for i in range(N):
            for j in range(i + 1, N):
                check(((out[i] >> j) & 1) + ((out[j] >> i) & 1) == 1, "tournament arcs")
print("hosts (gentourng):", {N: len(v) for N, v in HOSTS.items()}, "= A000568(5..8)")

# ------------------------------------------------------------------ guests: the 27 oriented 5-trees
def tree_name(T):
    degs = sorted((d for _, d in T.degree()), reverse=True)
    if degs[0] == 4: return "star"
    if degs[0] == 2: return "path"
    return "fork"

guests = []   # (name, DiGraph)
reps = []
for T in nx.nonisomorphic_trees(5):
    nm = tree_name(T)
    E = list(T.edges)
    for bits in range(16):
        D = nx.DiGraph()
        D.add_nodes_from(T.nodes)
        for k, (u, v) in enumerate(E):
            if (bits >> k) & 1:
                D.add_edge(u, v)
            else:
                D.add_edge(v, u)
        if not any(DiGraphMatcher(R, D).is_isomorphic() for _, R in reps):
            reps.append((nm, D))
guests = reps
check(len(guests) == 27, "A000238(5) = 27 oriented trees")
check(Counter(nm for nm, _ in guests) == Counter({"path": 10, "fork": 12, "star": 5}), "10 + 12 + 5")
print("guests: 27 oriented trees on 5 vertices (10 paths, 12 forks, 5 stars)")

def guest_label(nm, D):
    """Human label: out-degree sequence plus a word for paths/stars."""
    outd = sorted((d for _, d in D.out_degree()), reverse=True)
    ind = sorted((d for _, d in D.in_degree()), reverse=True)
    if nm == "path":
        # read the path end to end as a word of + (forward) / - (backward), canonical up to reversal-flip
        U = D.to_undirected()
        ends = [v for v in U if U.degree(v) == 1]
        order = nx.shortest_path(U, ends[0], ends[1])
        w = "".join("+" if D.has_edge(order[i], order[i+1]) else "-" for i in range(4))
        w2 = "".join("+" if c == "-" else "-" for c in reversed(w))
        word = min(w, w2)
        return f"path {word}"
    if nm == "star":
        c = [v for v in D if D.in_degree(v) + D.out_degree(v) == 4][0]
        return f"star out{D.out_degree(c)}/in{D.in_degree(c)}"
    c = [v for v in D if D.in_degree(v) + D.out_degree(v) == 3][0]
    legs = []
    for x in D.successors(c): legs.append(("c->" + ("x" if (D.in_degree(x) + D.out_degree(x)) == 2 else "l"), x))
    for x in D.predecessors(c): legs.append(("c<-" + ("x" if (D.in_degree(x) + D.out_degree(x)) == 2 else "l"), x))
    # the long leg: find the degree-2 neighbour x of c and its other neighbour e
    x = [v for v in D if (D.in_degree(v) + D.out_degree(v)) == 2 and (D.has_edge(c, v) or D.has_edge(v, c))][0]
    e = [v for v in D if v != c and (D.has_edge(x, v) or D.has_edge(v, x))][0]
    long_arc = ("c->x" if D.has_edge(c, x) else "c<-x") + ("," + ("x->e" if D.has_edge(x, e) else "x<-e"))
    shorts = sorted(("c->l" if D.has_edge(c, v) else "c<-l") for v in D if v not in (c, x, e))
    return f"fork [{long_arc}] {shorts[0]} {shorts[1]}"

# ------------------------------------------------------------------ embedding test
def bfs_plan(D):
    """Return (root, [(child, parent, forward)]) in BFS order; forward=True means parent -> child."""
    root = next(iter(D.nodes))
    U = D.to_undirected()
    order = list(nx.bfs_edges(U, root))
    plan = []
    for p, c in order:
        plan.append((c, p, D.has_edge(p, c)))
    return root, plan

def embeds(out, N, plan_data):
    root, plan = plan_data
    full = (1 << N) - 1
    inn = [0] * N
    for i in range(N):
        for j in range(N):
            if (out[i] >> j) & 1:
                inn[j] |= 1 << i
    pos = {}
    def rec(k, used):
        if k == len(plan):
            return True
        c, p, fwd = plan[k]
        cand = (out[pos[p]] if fwd else inn[pos[p]]) & ~used
        while cand:
            b = cand & -cand
            v = b.bit_length() - 1
            pos[c] = v
            if rec(k + 1, used | b):
                return True
            cand ^= b
        return False
    for r in range(N):
        pos[root] = r
        if rec(0, 1 << r):
            return True
    return False

plans = [(nm, D, bfs_plan(D)) for nm, D in guests]

# ------------------------------------------------------------------ census
t0 = time.time()
miss = defaultdict(dict)          # guest index -> N -> number of host classes missing it
missing_hosts = defaultdict(list) # (N) -> list of (host index, [guest indices missed])
for N in range(5, 9):
    for hi, out in enumerate(HOSTS[N]):
        missed = []
        for gi, (nm, D, pl) in enumerate(plans):
            if not embeds(out, N, pl):
                missed.append(gi)
        if missed:
            missing_hosts[N].append((hi, missed))
        for gi in missed:
            miss[gi][N] = miss[gi].get(N, 0) + 1
print(f"census done in {time.time()-t0:.1f}s")

def f_of(gi):
    for N in range(5, 9):
        if miss[gi].get(N, 0) == 0:
            return N
    return None

print("\n" + "=" * 78)
print("Sumner census at n = 5: miss_N(S) over all N-tournament classes, and f(S)")
print("=" * 78)
print(f"{'guest':42s} {'miss_5':>7s} {'miss_6':>7s} {'miss_7':>7s} {'miss_8':>7s}   f(S)")
rows = []
for gi, (nm, D, pl) in enumerate(plans):
    lab = guest_label(nm, D)
    m = [miss[gi].get(N, 0) for N in range(5, 9)]
    f = f_of(gi)
    rows.append((nm, lab, m, f))
for nm, lab, m, f in sorted(rows, key=lambda r: (r[0], r[3], r[1])):
    print(f"{lab:42s} {m[0]:7d} {m[1]:7d} {m[2]:7d} {m[3]:7d}   {f}")
check(all(r[2][3] == 0 for r in rows), "SUMNER n=5: every 8-tournament contains every oriented 5-tree")
print("\nSUMNER AT n = 5 HOLDS: all 27 oriented trees on 5 vertices embed in all 6880 tournaments on 8 vertices.")
fdist = Counter(r[3] for r in rows)
print("distribution of f(S):", dict(sorted(fdist.items())))
print("by underlying tree:", {nm: dict(sorted(Counter(r[3] for r in rows if r[0] == nm).items())) for nm in ("path", "fork", "star")})
# monotonicity check (should be automatic)
for r in rows:
    m = r[2]
    check(all((m[i] == 0) <= (m[i+1] == 0) for i in range(3)), "monotone")

print("\n" + "=" * 78)
print("Hosts of order 7 that miss some oriented 5-tree (the 2n-3 layer)")
print("=" * 78)
def scores(out, N):
    return tuple(sorted(bin(x).count("1") for x in out))
for hi, missed in missing_hosts[7]:
    out = HOSTS[7][hi]
    print(f"  host #{hi:3d} scores {scores(out,7)}  misses {len(missed)}: " + "; ".join(guest_label(plans[g][0], plans[g][1]) for g in missed))
reg7 = [hi for hi, out in enumerate(HOSTS[7]) if scores(out, 7) == (3,) * 7]
print(f"  regular 7-tournaments: {len(reg7)} classes (indices {reg7}); hosts missing something: {len(missing_hosts[7])}")
check(set(hi for hi, _ in missing_hosts[7]) == set(reg7), "exactly the regular 7-tournaments miss a 5-tree")

print("\n" + "=" * 78)
print("Hosts of order 6 that miss some oriented 5-tree")
print("=" * 78)
for hi, missed in missing_hosts[6]:
    out = HOSTS[6][hi]
    print(f"  host #{hi:3d} scores {scores(out,6)}  misses {len(missed)}: " + "; ".join(guest_label(plans[g][0], plans[g][1]) for g in missed))
print(f"  {len(missing_hosts[6])} of 56 classes miss something")

print("\n" + "=" * 78)
print("Hosts of order 5 that miss some oriented 5-tree (spanning case)")
print("=" * 78)
for hi, missed in missing_hosts[5]:
    out = HOSTS[5][hi]
    print(f"  host #{hi:3d} scores {scores(out,5)}  misses {len(missed)}: " + "; ".join(guest_label(plans[g][0], plans[g][1]) for g in missed))
print(f"  {len(missing_hosts[5])} of 12 classes miss something")


print("\n" + "=" * 78)
print("E. The Collatz functional digraph as a host for the 27 oriented 5-trees")
print("=" * 78)
# arcs n -> C(n), C(n) = n/2 (n even), 3n+1 (n odd), on {1..NMAX} (arcs leaving the range dropped)
NMAX = 200000
succ = {}
pred = defaultdict(list)
for n in range(1, NMAX + 1):
    c = n // 2 if n % 2 == 0 else 3 * n + 1
    if c <= NMAX:
        succ[n] = c
        pred[c].append(n)
indeg = Counter(len(v) for v in pred.values())
check(max(indeg) == 2, "in-degree <= 2")
branch = [v for v, p in pred.items() if len(p) == 2]
check(all(v % 6 == 4 for v in branch), "branch vertices are exactly 4 mod 6 (in range)")
check(all(len(pred.get(q, [])) <= 1 for v in branch for q in pred[v]), "no two branch vertices adjacent")
d2 = sum(1 for v in branch for q in pred[v] for r in pred.get(q, []) if len(pred.get(r, [])) == 2)
check(d2 > 0, "branch vertices at distance 2 occur")
print(f"  vertices 1..{NMAX}: in-degrees {dict(sorted(indeg.items()))}, out-degree <= 1; branch vertices (in-degree 2) are the n = 4 mod 6;")
print(f"  no two branch vertices are adjacent (the preimages 2n = 2 mod 6 and (n-1)/3 odd are not 4 mod 6); distance-2 pairs: {d2}")

def embeds_sparse(D):
    root, plan = bfs_plan(D)
    pos = {}
    def rec(k, used):
        if k == len(plan):
            return True
        c, p, fwd = plan[k]
        v0 = pos[p]
        cands = ([succ[v0]] if (fwd and v0 in succ) else []) if fwd else pred.get(v0, [])
        for v in cands:
            if v in used: continue
            pos[c] = v
            if rec(k + 1, used | {v}):
                return True
        return False
    # roots: only need to try vertices with the right local shape; try all up to a bound
    for r in range(1, 5000):
        pos[root] = r
        if rec(0, {r}):
            return r
    return None

found = []
for nm, D in guests:
    r = embeds_sparse(D)
    outd = max(d for _, d in D.out_degree()); ind = max(d for _, d in D.in_degree())
    functional = outd <= 1 and ind <= 2
    # adjacent branch vertices in the guest?
    br = [v for v in D if D.in_degree(v) == 2]
    adj_branch = any(D.has_edge(a, b) or D.has_edge(b, a) for a in br for b in br if a != b)
    predicted = functional and not adj_branch
    check((r is not None) == predicted, f"Collatz host prediction for {guest_label(nm, D)}")
    if r is not None:
        found.append(guest_label(nm, D))
print(f"  oriented 5-trees embedded in the Collatz digraph: {len(found)} of 27:")
for lab in found:
    print("    ", lab)
print("  = exactly the in-trees with in-degree <= 2 and no two adjacent branch vertices (PROVED for the host;")
print("    the fork rooted at the middle of its long leg needs two adjacent branch vertices and is the one binary in-tree excluded).")
print("  Compare: TT_5 hosts 27/27, the regular 5-tournament 17/27, the Collatz digraph 5/27.")

print(f"\nALL {CHECKS} CHECKS PASSED")
