"""
Check whether Cayley tournaments on Z/7 : Z/3 (Frobenius group, order 21)
are self-converse using full backtracking search.
"""
import os, sys
os.environ['PYTHONIOENCODING'] = 'utf-8'

def mul(x, y):
    a1, b1 = x
    a2, b2 = y
    return ((a1 + pow(2, b1, 7) * a2) % 7, (b1 + b2) % 3)

def inv_g(x):
    a, b = x
    return ((-pow(2, (3 - b) % 3, 7) * a) % 7, (3 - b) % 3)

e = (0, 0)
elements = [(a, b) for a in range(7) for b in range(3)]
n = 21
idx = {g: i for i, g in enumerate(elements)}

# Inversion pairs
non_id = [g for g in elements if g != e]
pairs = []
seen_pairs = set()
for g in non_id:
    fg = frozenset([g, inv_g(g)])
    if fg not in seen_pairs:
        pairs.append((g, inv_g(g)))
        seen_pairs.add(fg)

# Build a tournament with S = first element of each pair
S = set(p[0] for p in pairs)

# Build adjacency bitmasks
adj = [0] * n
for g in elements:
    gi = idx[g]
    for h in elements:
        if g != h and mul(inv_g(g), h) in S:
            adj[gi] |= (1 << idx[h])

# T^op bitmasks
adj_op = [0] * n
for i in range(n):
    for j in range(n):
        if adj[i] & (1 << j):
            adj_op[j] |= (1 << i)

# We want sigma such that: adj[i] has bit j set iff adj[sigma(j)] has bit sigma(i) set
# i.e., i->j in T iff sigma(j)->sigma(i) in T (anti-isomorphism = isomorphism T -> T^op)

# Since T is vertex-transitive, we can fix sigma(0) = any vertex.
# By VT of T^op, we can fix sigma(0) = 0.

# Key observation: if i->j in T, then sigma(j)->sigma(i) in T.
# So 0->j in T means sigma(j)->sigma(0)=0 in T, i.e., sigma(j) is in in-nbrs of 0.
# And j->0 in T means 0->sigma(j) in T, i.e., sigma(j) is in out-nbrs of 0.
# So sigma maps out-nbrs(0) to in-nbrs(0) and vice versa.

out_0 = [j for j in range(n) if j > 0 and (adj[0] & (1 << j))]
in_0 = [j for j in range(n) if j > 0 and (adj[j] & (1 << 0))]
print(f"out-degree of 0: {len(out_0)}, in-degree of 0: {len(in_0)}")

# Profile each vertex by its score within out_0 and in_0
def get_profile(v, source_set, adj_bits):
    return sum(1 for w in source_set if w != v and (adj_bits[v] & (1 << w)))

out_in_out = {}  # for u in out_0: (score in out_0, score in in_0)
for u in out_0:
    key = (get_profile(u, out_0, adj), get_profile(u, in_0, adj))
    out_in_out.setdefault(key, []).append(u)

# For T^op: sigma maps out_0 -> in_0.
# In T^op, sigma(u) must satisfy: for sigma(u) in in_0,
# "score of sigma(u) among in_0 in T^op" should match "score of u among out_0 in T"
# T^op score of v among in_0 = number of w in in_0 with adj_op[v][w] = adj[w][v]
in_in_in = {}
for v in in_0:
    s1 = sum(1 for w in in_0 if w != v and (adj[w] & (1 << v)))  # T^op: v beats w in in_0
    s2 = sum(1 for w in out_0 if (adj[w] & (1 << v)))  # T^op: v beats w in out_0
    key = (s1, s2)
    in_in_in.setdefault(key, []).append(v)

print("Profiles of out_0 (in T):", {k: len(v) for k, v in sorted(out_in_out.items())})
print("Profiles of in_0 (in T^op):", {k: len(v) for k, v in sorted(in_in_in.items())})

# Backtracking search
# sigma[0] = 0 (fixed)
# sigma maps out_0 -> in_0 and in_0 -> out_0, preserving profiles

found_iso = [False]

def check_pair(sigma, i, j):
    """Check if the pair (i,j) is consistent with sigma being an anti-isomorphism."""
    if i not in sigma or j not in sigma:
        return True  # can't check yet
    si, sj = sigma[i], sigma[j]
    i_to_j = bool(adj[i] & (1 << j))
    sj_to_si = bool(adj[sj] & (1 << si))
    return i_to_j == sj_to_si

def backtrack(sigma, remaining_src, remaining_tgt, assigned_src_list):
    """Try to extend sigma to map remaining_src -> remaining_tgt."""
    if found_iso[0]:
        return

    if not remaining_src:
        # Check full consistency (should already be consistent from pruning)
        found_iso[0] = True
        print("FOUND ANTI-ISOMORPHISM!")
        for k in sorted(sigma.keys()):
            g = elements[k]
            sg = elements[sigma[k]]
            # print(f"  {g} -> {sg}")
        return

    # Pick the source vertex with fewest compatible targets (MRV heuristic)
    best_src = None
    best_compatible = None
    for src in remaining_src:
        compatible = []
        for tgt in remaining_tgt:
            # Check if sigma[src] = tgt is consistent with all existing assignments
            ok = True
            for assigned in assigned_src_list:
                # Check (src, assigned) pair
                src_to_a = bool(adj[src] & (1 << assigned))
                sa = sigma[assigned]
                sa_to_tgt = bool(adj[sa] & (1 << tgt))
                if src_to_a != sa_to_tgt:
                    ok = False
                    break
                # Check (assigned, src) pair
                a_to_src = bool(adj[assigned] & (1 << src))
                tgt_to_sa = bool(adj[tgt] & (1 << sa))
                if a_to_src != tgt_to_sa:
                    ok = False
                    break
            if ok:
                compatible.append(tgt)
        if best_src is None or len(compatible) < len(best_compatible):
            best_src = src
            best_compatible = compatible
        if len(compatible) == 0:
            return  # Dead end

    # Try each compatible target
    remaining_src_new = remaining_src - {best_src}
    for tgt in best_compatible:
        sigma[best_src] = tgt
        assigned_src_list.append(best_src)
        remaining_tgt_new = remaining_tgt - {tgt}
        backtrack(sigma, remaining_src_new, remaining_tgt_new, assigned_src_list)
        if found_iso[0]:
            return
        assigned_src_list.pop()
        del sigma[best_src]

# Phase 1: map out_0 -> in_0
sigma = {0: 0}
print("\nPhase 1: mapping out_0 -> in_0...")
backtrack(sigma, set(out_0), set(in_0), [0])

if found_iso[0]:
    print("Phase 1 succeeded. Extending to in_0 -> out_0...")
    # Phase 2: map in_0 -> out_0
    found_iso[0] = False
    backtrack(sigma, set(in_0), set(out_0), [0] + out_0)
    if found_iso[0]:
        print("\nSUCCESS: Tournament IS self-converse!")
    else:
        print("\nNo extension found for in_0 -> out_0 from this partial.")
        print("Trying full search...")
else:
    print("Phase 1 found no valid mapping out_0 -> in_0")

# If phase-by-phase failed, do full unified search
if not found_iso[0]:
    print("\nFull unified backtracking search...")
    sigma = {0: 0}
    all_others = [j for j in range(1, n)]
    backtrack(sigma, set(out_0 + in_0), set(out_0 + in_0), [0])
    if found_iso[0]:
        print("SUCCESS: Tournament IS self-converse!")
    else:
        print("FAILURE: Tournament is NOT self-converse!")

print("\nDone.")
