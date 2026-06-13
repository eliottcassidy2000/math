#!/usr/bin/env python3
"""
beta2_delta_ratio_proof.py - Prove delta_|A_3| = (n-3)*delta_|A_2|

THEOREM: For any tournament T on n vertices and any arc (u->v),
let T' be the tournament obtained by flipping (u->v) to (v->u). Then:
  |A_3(T')| - |A_3(T)| = (n-3) * (|A_2(T')| - |A_2(T)|)

where A_p(T) = set of allowed p-paths (sequences of p+1 distinct vertices
following arcs).

PROOF IDEA:
Let out(v) = {w : v->w} \ {u}, in(u) = {w : w->u} \ {v}.
(These are the OUT neighbors of v excluding u, and IN neighbors of u excluding v.)

After flip: u becomes an in-neighbor of v (was out), v becomes out-neighbor of u.
Define: a = |out(u) cap in(v)| \ {u,v} etc. -- the ABCD partition.

LOST 2-PATHS through u->v:
  Position 0: (u, v, w) for w in out(v)\{u} -- |out(v)|-[v->u] paths
  Position 1: (w, u, v) for w in in(u)\{v} -- |in(u)|-[u->v doesn't apply, w->u] paths
  Total lost_2 = |out(v)\{u}| + |in(u)\{v}|

GAINED 2-PATHS through v->u in T':
  Position 0: (v, u, w) for w in out_T'(u)\{v} = (out_T(u)\{v})\{u} = out_T(u)\{v}
  Position 1: (w, v, u) for w in in_T'(v)\{u} = (in_T(v)\{u}) union {u}\{u} wait...

After flip: T' has v->u instead of u->v.
  out_T'(u) = out_T(u) \ {v}  (u lost the arc to v)
  in_T'(u) = in_T(u) union {v}  (u gained incoming from v)
  out_T'(v) = out_T(v) union {u}  (v gained arc to u... wait, v->u means v has arc TO u?)

No: flip (u->v) to (v->u) means:
  In T: u->v. In T': v->u.
  So out_T'(v) = out_T(v) (v's other arcs unchanged)... wait.
  In T, A[u][v]=1 (u->v). A[v][u]=0.
  In T', A'[u][v]=0, A'[v][u]=1 (v->u).
  So out_T'(u) = out_T(u) \ {v} (u lost v from outgoing)
  And out_T'(v) = out_T(v) union {u} (v gained u in outgoing)
  in_T'(u) = in_T(u) union {v}
  in_T'(v) = in_T(v) \ {u}

GAINED 2-PATHS through v->u in T':
  Position 0: (v, u, w) for w in out_T'(u)\{v} = out_T(u)\{v}. Count = |out_T(u)| - 1
    (since v was in out_T(u), we remove it)
    Wait: out_T(u) includes v (since u->v in T). So out_T'(u) = out_T(u)\{v}.
    Paths: (v, u, w) for w in out_T'(u) = out_T(u)\{v}. These w != v, w != u.
    Count = |out_T(u)| - 1  (removing v from out_T(u))

  Position 1: (w, v, u) for w in in_T'(v)\{u} = (in_T(v)\{u})\{u}...
    in_T(v) = {w : w->v in T, w != v}. Does u -> v in T? Yes. So u in in_T(v)?
    No: in_T(v) means w->v, i.e., A[w][v]=1. We have A[u][v]=1, so u IS in in_T(v).
    in_T'(v) = in_T(v) \ {u}  (since in T', u no longer points to v)
    Paths: (w, v, u) for w in in_T'(v) = in_T(v)\{u}. w != u, w != v.
    Count = |in_T(v)| - 1  (removing u from in_T(v))

So:
  lost_2 = |out_T(v)\{u}| + |in_T(u)\{v}|
         = (n-1 - |in_T(v)\{u}| - 1) + (n-1 - |out_T(u)\{v}| - 1)
         Hmm, let me be more careful.

In tournament T on n vertices:
  out_T(v) includes all w != v with v->w. u may or may not be in out_T(v).
  Since u->v in T, v does NOT have arc to u. So u is NOT in out_T(v).

  out_T(v) = {w != v : v->w}. Since u->v, we know v->u is false, so u not in out_T(v).
  |out_T(v)| = out-degree of v in T.
  out_T(v) \ {u} = out_T(v) (since u not in it anyway).

  in_T(u) = {w != u : w->u}. Since u->v, v->u is false, so v not in in_T(u).
  in_T(u) \ {v} = in_T(u).

So:
  lost_2 = |out_T(v)| + |in_T(u)|

  gained_2 = |out_T(u)\{v}| + |in_T(v)\{u}|
           = (|out_T(u)| - 1) + (|in_T(v)| - 1)

  delta_2 = gained_2 - lost_2
          = (|out_T(u)| - 1 + |in_T(v)| - 1) - (|out_T(v)| + |in_T(u)|)
          = |out_T(u)| + |in_T(v)| - |out_T(v)| - |in_T(u)| - 2

In a tournament: |out_T(w)| + |in_T(w)| = n-1 for each w.
  So |in_T(u)| = n-1 - |out_T(u)| and |in_T(v)| = n-1 - |out_T(v)|.

  delta_2 = |out_T(u)| + (n-1-|out_T(v)|) - |out_T(v)| - (n-1-|out_T(u)|) - 2
          = 2*|out_T(u)| - 2*|out_T(v)| - 2
          = 2*(|out_T(u)| - |out_T(v)| - 1)

Let d_u = |out_T(u)|, d_v = |out_T(v)|. Then:
  delta_2 = 2*(d_u - d_v - 1)

Note: d_u >= 1 (since u->v), and d_v can be anything from 0 to n-2 (but v->u is false,
so v's out-degree is among arcs to other n-2 vertices minus u... well, it's just the
out-degree.)

Now for 3-paths:
LOST 3-paths through u->v: arc at position 0, 1, or 2 in a length-3 path.

Position 0: (u, v, w, x) where v->w, w->x, all distinct, {w,x} disjoint from {u,v}
  = # of allowed 2-paths (v, w, x) in V\{u} starting from v
  But w,x must also avoid u.
  = # of allowed 2-paths in the sub-tournament T[V\{u}] starting from v
  Wait: v->w and w->x with w,x not in {u,v} and w!=x.
  The 2-path (v,w,x) must use vertices in V\{u} and the arcs are from T.
  Since T is a tournament, all arcs between V\{u} are present in T.
  So this equals: # of 2-paths starting at v in the tournament T - {u} (on n-1 vertices).

Hmm, this is getting complicated for a general formula. Let me verify computationally
that delta_2 = 2*(d_u - d_v - 1) and then find the formula for delta_3.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, os, random
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import enumerate_allowed_paths
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

# Verify delta_2 = 2*(d_u - d_v - 1)
print("=" * 70)
print("VERIFYING DELTA_2 FORMULA")
print("=" * 70)

for n in [5, 6, 7, 8]:
    violations = 0
    total = 0
    for _ in range(500):
        A = random_tournament(n)
        for _ in range(3):
            u = random.randint(0, n-1)
            v = random.randint(0, n-1)
            while u == v or A[u][v] == 0:
                u = random.randint(0, n-1)
                v = random.randint(0, n-1)

            d_u = sum(A[u])
            d_v = sum(A[v])

            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1

            a2_before = len(enumerate_allowed_paths(A, n, 2))
            a2_after = len(enumerate_allowed_paths(B, n, 2))
            delta_2 = a2_after - a2_before
            formula = 2 * (d_u - d_v - 1)

            total += 1
            if delta_2 != formula:
                violations += 1
                if violations <= 3:
                    print(f"  VIOLATION at n={n}: d_u={d_u}, d_v={d_v}, delta_2={delta_2}, formula={formula}")

    print(f"  n={n}: {violations}/{total} violations")

# Now find formula for delta_3
print(f"\n{'='*70}")
print("FINDING DELTA_3 FORMULA")
print(f"{'='*70}")

# We know delta_3 = (n-3) * delta_2 = (n-3) * 2 * (d_u - d_v - 1)
# But let's also express it directly in terms of the tournament structure.

# Try: delta_3 in terms of d_u, d_v, and the ABCD partition
# A = out(u) cap in(v) (w: u->w, w->v)
# B = in(u) cap out(v)  (w: w->u, v->w)
# C = out(u) cap out(v) (w: u->w, v->w)
# D = in(u) cap in(v)   (w: w->u, w->v)
# Note: u->v in T, so v not counted in out(u) in the A,B,C,D partition
# Actually let's be careful about what's in these sets.

# For the ABCD partition, consider w in V\{u,v}:
# u->w or w->u (since tournament), and w->v or v->w (since tournament)
# Combined with u->v in T:
#   A_set: u->w AND w->v
#   B_set: w->u AND v->w
#   C_set: u->w AND v->w
#   D_set: w->u AND w->v

# |A| + |C| = |out(u) cap V\{u,v}| = d_u - 1  (d_u includes v)
# |B| + |D| = |in(u) cap V\{u,v}| = n-1-d_u  (= in(u), which excludes v since v->u false)
# |A| + |D| = |in(v) cap V\{u,v}| = n-1-d_v - 1 = n-2-d_v
#   Wait: in(v) = {w: w->v, w != v}. u->v so u in in(v). But A,D are subsets of V\{u,v}.
#   |A| + |D| = |{w in V\{u,v} : w->v}| = |in(v)| - 1 = (n-1-d_v) - 1 = n-2-d_v
# |B| + |C| = |{w in V\{u,v} : v->w}| = d_v

# So:
#   |A| = (d_u-1) - |C|
#   |B| = d_v - |C|
#   |D| = (n-2-d_v) - |A| = (n-2-d_v) - (d_u-1) + |C| = n-1-d_v-d_u+|C|
# And |A|+|B|+|C|+|D| = n-2. Check: (d_u-1-|C|) + (d_v-|C|) + |C| + (n-1-d_v-d_u+|C|) = n-2. Yes.

print("\nVerifying delta_3 = 2*(n-3)*(d_u - d_v - 1):")
for n in [5, 6, 7, 8]:
    violations = 0
    total = 0
    for _ in range(500):
        A = random_tournament(n)
        for _ in range(3):
            u = random.randint(0, n-1)
            v = random.randint(0, n-1)
            while u == v or A[u][v] == 0:
                u = random.randint(0, n-1)
                v = random.randint(0, n-1)

            d_u = sum(A[u])
            d_v = sum(A[v])

            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1

            a3_before = len(enumerate_allowed_paths(A, n, 3))
            a3_after = len(enumerate_allowed_paths(B, n, 3))
            delta_3 = a3_after - a3_before
            formula = 2 * (n - 3) * (d_u - d_v - 1)

            total += 1
            if delta_3 != formula:
                violations += 1
                if violations <= 3:
                    print(f"  VIOLATION at n={n}: d_u={d_u}, d_v={d_v}, delta_3={delta_3}, formula={formula}")

    print(f"  n={n}: {violations}/{total} violations")

# And verify that both formulas hold simultaneously
print(f"\n{'='*70}")
print("COMBINED FORMULA VERIFICATION")
print(f"{'='*70}")
print("delta_|A_2| = 2*(d_u - d_v - 1)")
print("delta_|A_3| = 2*(n-3)*(d_u - d_v - 1)")
print()

for n in [4, 5, 6, 7, 8, 9]:
    v2 = 0
    v3 = 0
    total = 0
    for _ in range(300):
        A = random_tournament(n)
        for _ in range(3):
            u = random.randint(0, n-1)
            v = random.randint(0, n-1)
            while u == v or A[u][v] == 0:
                u = random.randint(0, n-1)
                v = random.randint(0, n-1)

            d_u = sum(A[u])
            d_v = sum(A[v])

            B = [row[:] for row in A]
            B[u][v] = 0
            B[v][u] = 1

            a2_b = len(enumerate_allowed_paths(A, n, 2))
            a2_a = len(enumerate_allowed_paths(B, n, 2))
            a3_b = len(enumerate_allowed_paths(A, n, 3))
            a3_a = len(enumerate_allowed_paths(B, n, 3))

            total += 1
            if (a2_a - a2_b) != 2*(d_u - d_v - 1):
                v2 += 1
            if (a3_a - a3_b) != 2*(n-3)*(d_u - d_v - 1):
                v3 += 1

    print(f"  n={n}: delta_2 violations={v2}/{total}, delta_3 violations={v3}/{total}")

print("\nDone.")
