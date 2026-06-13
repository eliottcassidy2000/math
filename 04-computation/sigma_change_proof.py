"""
sigma_change_proof.py -- kind-pasteur-2026-03-13-S61

PROOF that Vitali atoms change sigma for exactly 4*(n-4) pairs.

sigma(u,v) = #{common successors of u,v} + #{common predecessors of u,v}
           = #{w: (u->w AND v->w) OR (w->u AND w->v)}

Under reversal of S = {a,b,c,d} (4-vertex subset with scores (1,1,2,2)):
  For arcs within S: all directions flip.
  For arcs between S and outside: unchanged.

Which sigma(u,v) values change?

Case 1: u,v both outside S (XX pair).
  sigma(u,v) = #{w: u->w,v->w OR w->u,w->v}
  For w outside S: arcs u-w, v-w unchanged. Contribution unchanged.
  For w in S: arcs u-w and v-w are NOT arcs within S. They are arcs
  between S and outside, which are UNCHANGED.
  So sigma(u,v) is UNCHANGED for XX pairs. PROVED.

Case 2: u in S, v outside S (SX pair).
  sigma(u,v) = #{w: u->w,v->w OR w->u,w->v}
  For w outside S (w != v): arcs u-w and v-w are unchanged (u-w is S-outside,
    v-w is outside-outside). Contribution unchanged.
  For w in S (w != u): arc u-w is WITHIN S and FLIPS. Arc v-w is S-outside,
    unchanged. So the contribution of w changes.
    Before: (u->w means A[u][w]=1) and (v->w means A[v][w]=1): both out.
            (w->u means A[w][u]=1) and (w->v means A[w][v]=1): both in.
    After: u->w flips to w->u (or vice versa). v->w unchanged.
    So the "common successor/predecessor" status of w w.r.t. (u,v) changes.

  Does EVERY w in S\{u} cause a sigma change? Not necessarily.
  It changes when the u-w arc direction matters for the common successor status.

  Let's count precisely.

Case 3: u,v both in S (SS pair).
  sigma(u,v) = #{w: u->w,v->w OR w->u,w->v}
  For w outside S: arcs u-w and v-w are S-outside, unchanged. Contribution unchanged.
  For w in S: arcs u-w, v-w, and u-v all flip simultaneously.
    Before: u->w?, v->w?. After: w->u?, w->v? (flipped arcs).
    But: common successor means u->w AND v->w. After flip: w->u AND w->v.
    So "common successor" becomes "common predecessor" and vice versa.
    The sum sigma = #(common succ) + #(common pred) is PRESERVED under swap!
    (because we just swap the two sets)

  Wait: but u-v itself also flips. sigma(u,v) doesn't involve the arc u-v
  directly (sigma counts witnesses w != u,v). So sigma(u,v) for SS pairs:
  - Outside w: unchanged
  - Inside w (in S, != u,v): arcs u-w and v-w both flip. The "both out" becomes
    "both in" and vice versa. But sigma counts EITHER, so sigma is preserved!

  PROOF: For w in S\{u,v}: after flipping all S-arcs:
    "u->w AND v->w" becomes "w->u AND w->v" (common successor -> common predecessor)
    "w->u AND w->v" becomes "u->w AND v->w" (common predecessor -> common successor)
    "u->w AND w->v" stays as "w->u AND v->w" — but wait, v->w also flips to w->v.
    Actually ALL arcs u-w, v-w, w-u, w-v flip.
    So A'[u][w] = A[w][u] and A'[v][w] = A[w][v].
    Common successor: A'[u][w]*A'[v][w] = A[w][u]*A[w][v] = common predecessor before.
    Common predecessor: A'[w][u]*A'[w][v] = A[u][w]*A[v][w] = common successor before.
    sigma' = #(new common succ) + #(new common pred)
           = #(old common pred) + #(old common succ) = sigma.

  So sigma(u,v) is UNCHANGED for SS pairs. PROVED.

Therefore: sigma changes ONLY for SX pairs.
There are 4*(n-4) SX pairs. But does EVERY SX pair change?
"""

import numpy as np
from itertools import combinations
from collections import Counter

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

# Test the proof at multiple n values
for n in [7, 8, 9]:
    total_bits = n * (n-1) // 2
    print(f"\n{'='*60}")
    print(f"SIGMA CHANGE ANALYSIS AT n={n}")
    print(f"{'='*60}")

    np.random.seed(42)
    n_SX_pairs = 4 * (n - 4)
    change_counts = []

    for trial in range(2000):
        if total_bits <= 31:
            bits = np.random.randint(0, 1 << total_bits)
        else:
            bits = int(np.random.randint(0, 2**31)) | (int(np.random.randint(0, 2**(total_bits-31))) << 31)
        A = bits_to_adj(bits, n)
        L = lambda_graph(A, n)

        for subset in combinations(range(n), 4):
            ss = sub_scores(A, n, list(subset))
            if ss != (1, 1, 2, 2):
                continue
            B = reverse_subtournament(A, n, list(subset))
            if not np.array_equal(L, lambda_graph(B, n)):
                continue

            S = set(subset)
            A2_A = A @ A
            A2_B = B @ B

            # Count changes by pair type
            xx_changes = 0
            sx_changes = 0
            ss_changes = 0

            for u in range(n):
                for v in range(u+1, n):
                    sig_A = n - 2 - int(A2_A[u][v]) - int(A2_A[v][u])
                    sig_B = n - 2 - int(A2_B[u][v]) - int(A2_B[v][u])
                    if sig_A != sig_B:
                        u_in = u in S
                        v_in = v in S
                        if u_in and v_in:
                            ss_changes += 1
                        elif u_in or v_in:
                            sx_changes += 1
                        else:
                            xx_changes += 1

            change_counts.append({
                'xx': xx_changes,
                'sx': sx_changes,
                'ss': ss_changes,
                'total': xx_changes + sx_changes + ss_changes,
            })
            break

    if not change_counts:
        print(f"  No Vitali pairs found at n={n}")
        continue

    # Non-trivial (sigma actually changes)
    nontrivial = [c for c in change_counts if c['total'] > 0]
    print(f"  Vitali pairs: {len(change_counts)}, non-trivial: {len(nontrivial)}")
    print(f"  Expected SX pairs: {n_SX_pairs}")

    if nontrivial:
        xx_dist = Counter(c['xx'] for c in nontrivial)
        sx_dist = Counter(c['sx'] for c in nontrivial)
        ss_dist = Counter(c['ss'] for c in nontrivial)
        total_dist = Counter(c['total'] for c in nontrivial)

        print(f"  XX changes: {dict(xx_dist)}")
        print(f"  SX changes: {dict(sx_dist)}")
        print(f"  SS changes: {dict(ss_dist)}")
        print(f"  Total changes: {dict(total_dist)}")
        print(f"  XX=0 always? {all(c['xx'] == 0 for c in nontrivial)}")
        print(f"  SS=0 always? {all(c['ss'] == 0 for c in nontrivial)}")
        print(f"  SX={n_SX_pairs} always? {all(c['sx'] == n_SX_pairs for c in nontrivial)}")

# Now prove WHY every SX pair changes sigma
print(f"\n{'='*60}")
print("WHY EVERY SX PAIR CHANGES SIGMA")
print(f"{'='*60}")

n = 7
total_bits = n * (n-1) // 2

# For SX pair (s, x) where s in S, x outside S:
# sigma(s, x) = #{w: common succ/pred of s,x}
# The witnesses w can be: (1) other S-vertices, (2) other outside vertices.
#
# For w outside S (w != x): arcs s-w (S-outside, unchanged) and x-w (outside-outside, unchanged).
#   So these contributions DON'T change.
#
# For w in S (w != s): arc s-w (WITHIN S, FLIPS) and x-w (S-outside, unchanged).
#   Before flip: s->w or w->s. After flip: w->s or s->w (opposite).
#   The "common successor" status of w for pair (s,x):
#     Before: s->w AND x->w (both out to w) = A[s][w]*A[x][w]
#     After:  A'[s][w]*A[x][w] = A[w][s]*A[x][w]  (since s,w both in S)
#   The "common predecessor" status of w for pair (s,x):
#     Before: w->s AND w->x = A[w][s]*A[w][x]
#     After:  A'[w][s]*A[w][x] = A[s][w]*A[w][x]
#
# sigma_before(s,x) from w: A[s][w]*A[x][w] + A[w][s]*A[w][x]
# sigma_after(s,x) from w:  A[w][s]*A[x][w] + A[s][w]*A[w][x]
#
# delta_sigma from w = (A[w][s]*A[x][w] + A[s][w]*A[w][x]) - (A[s][w]*A[x][w] + A[w][s]*A[w][x])
#                     = A[w][s]*(A[x][w] - A[w][x]) + A[s][w]*(A[w][x] - A[x][w])
#                     = (A[w][s] - A[s][w])*(A[x][w] - A[w][x])
#                     = -sign(s->w) * sign(x->w)
#
# where sign(u->v) = A[u][v] - A[v][u] in {-1, +1}.
#
# So delta_sigma(s,x) = sum_{w in S, w!=s} -sign(s->w)*sign(x->w)
#
# This is nonzero unless ALL w in S\{s} have the same product sign(s->w)*sign(x->w).
# But in a (1,1,2,2) tournament on S, not all sign(s->w) are the same.

print("Algebraic formula for delta_sigma:")
print("  For SX pair (s,x): delta_sigma = -sum_{w in S, w!=s} sign(s->w)*sign(x->w)")
print("  where sign(u->v) = A[u][v] - A[v][u] in {-1,+1}")
print()

# Verify this formula
np.random.seed(42)
formula_verified = 0
formula_tested = 0

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        S = set(subset)
        A2_A = A @ A
        A2_B = B @ B

        for s in sorted(S):
            for x in range(n):
                if x in S:
                    continue
                # Direct computation
                sig_A = n - 2 - int(A2_A[s][x]) - int(A2_A[x][s])
                sig_B = n - 2 - int(A2_B[s][x]) - int(A2_B[x][s])
                dsig_direct = sig_B - sig_A

                # Formula
                dsig_formula = 0
                for w in sorted(S):
                    if w == s:
                        continue
                    sign_sw = int(A[s][w]) - int(A[w][s])  # +1 if s->w, -1 if w->s
                    sign_xw = int(A[x][w]) - int(A[w][x])  # +1 if x->w, -1 if w->x
                    dsig_formula -= sign_sw * sign_xw

                formula_tested += 1
                if dsig_direct == dsig_formula:
                    formula_verified += 1
                else:
                    print(f"  FORMULA FAIL: s={s}, x={x}, direct={dsig_direct}, formula={dsig_formula}")

        break

print(f"Formula verified: {formula_verified}/{formula_tested}")

# Now: when is delta_sigma(s,x) = 0?
# delta_sigma = -sum_{w in S\{s}} sign(s->w)*sign(x->w)
# |S\{s}| = 3. So sum has 3 terms, each +1 or -1. Sum in {-3,-1,1,3}.
# delta_sigma = 0 iff sum is 0, which needs 3 terms summing to 0 — IMPOSSIBLE.
# Therefore delta_sigma(s,x) != 0 for ALL SX pairs!

print(f"\nKey insight: sum of 3 terms in {{-1,+1}} is ALWAYS odd, so NEVER 0.")
print(f"Therefore delta_sigma(s,x) is ALWAYS nonzero for SX pairs.")
print(f"This PROVES that sigma changes for ALL |S|*(n-|S|) = 4*(n-4) SX pairs.")

# Verify: delta_sigma values for SX pairs
print(f"\ndelta_sigma values for SX pairs:")
dsig_sx_vals = []
np.random.seed(42)
for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        S = set(subset)
        A2_A = A @ A
        A2_B = B @ B

        for s in sorted(S):
            for x in range(n):
                if x in S:
                    continue
                sig_A = n - 2 - int(A2_A[s][x]) - int(A2_A[x][s])
                sig_B = n - 2 - int(A2_B[s][x]) - int(A2_B[x][s])
                dsig_sx_vals.append(sig_B - sig_A)

        break

dsig_dist = Counter(dsig_sx_vals)
print(f"  delta_sigma distribution for SX pairs: {dict(sorted(dsig_dist.items()))}")
print(f"  Always +/-1? {all(abs(v) == 1 for v in dsig_sx_vals)}")

print("\nDone.")
