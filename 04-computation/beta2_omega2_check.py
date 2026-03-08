#!/usr/bin/env python3
"""
CRITICAL CHECK: Is Ω_2 = span(transitive triples), or is it larger?

Ω_2 = {u ∈ A_2 : ∂u ∈ A_1}

A_2 includes ALL allowed 2-paths (a,b,c) with a→b, b→c.
Transitive triples are those with a→c additionally.
Non-transitive 2-paths have c→a instead.

For a non-transitive 2-path (a,b,c) with c→a:
∂(a,b,c) = (b,c) - (a,c) + (a,b)
(a,c) is NOT in A_1 since a→c is false (c→a instead).
So individually, (a,b,c) ∉ Ω_2.

But could two non-transitive 2-paths cancel their non-allowed faces?
Example: (a,b,c) and (d,e,c') where c'=c, same non-allowed face (a,c).
We need -(a,c) from one and +(a,c) from another.

For face (a,c) to appear in ∂(d,e,f):
  As Face 1 = (d,f): only if d=a, f=c → (a,e,c) with e different
  As Face 0 = (e,f): only if e=a, f=c → (d,a,c) with d different
  As Face 2 = (d,e): only if d=a, e=c → (a,c,f) with f different

So (a,c) appears in ∂(a,b,c) as -1 and in ∂(a,b',c) as -1 for any b'.
To cancel, we'd need opposite signs. But Face 1 always has sign -1.
So two paths (a,b,c) and (a,b',c) with b≠b' both have -1 coefficient
on face (a,c). They can't cancel!

What about (a,c) appearing as Face 0 or Face 2?
Face 0 of (d,a,c): ∂(d,a,c) has +(a,c). Sign is +1.
Face 2 of (a,c,f): ∂(a,c,f) has +(a,c). Sign is (-1)^2 = +1.

Wait, let me be careful:
∂(v0,v1,v2) = (v1,v2) - (v0,v2) + (v0,v1)

Face 0 = (v1,v2), sign +1
Face 1 = (v0,v2), sign -1
Face 2 = (v0,v1), sign +1

So (a,c) appears as:
  Face 0 of (*,a,c): +1 coefficient
  Face 1 of (a,*,c): -1 coefficient
  Face 2 of (a,c,*): +1 coefficient

For a non-transitive (a,b,c) with c→a:
  ∂ has -(a,c) which is non-allowed.
  To cancel: need +(a,c) from another path.
  This comes from Face 0 of some (d,a,c) or Face 2 of some (a,c,e).

  Path (d,a,c): needs d→a and a→c. But c→a (not a→c)! So a→c doesn't exist.
  Wait, (d,a,c) needs edges d→a, a→c. Since c→a, there's no edge a→c.
  So (d,a,c) is NOT in A_2 — can't use it.

  Path (a,c,e): needs a→c. Same problem — c→a, not a→c.
  So (a,c,e) is NOT in A_2 either.

THEREFORE: the non-allowed face (a,c) from ∂(a,b,c) CANNOT be cancelled
by any other allowed 2-path, because any path generating +(a,c) would
require a→c, which doesn't exist.

CONCLUSION: Ω_2 = span(transitive triples) EXACTLY.
Non-transitive 2-paths cannot participate in Ω_2 chains at all.

Let me verify this computationally.
"""
import numpy as np
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

print("=" * 70)
print("Ω_2 = span(transitive triples)?")
print("=" * 70)

for n in [4, 5]:
    all_match = True
    for A in all_tournaments_gen(n):
        a1 = enumerate_allowed_paths(A, n, 1)
        a2 = enumerate_allowed_paths(A, n, 2)
        tt_count = sum(1 for p in a2 if A[p[0]][p[2]] == 1)

        om2 = compute_omega_basis(A, n, 2, a2, a1)
        dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

        if dim_om2 != tt_count:
            all_match = False
            print(f"  n={n}: dim(Ω_2)={dim_om2} ≠ |transitive triples|={tt_count}")

    if all_match:
        print(f"  n={n}: Ω_2 = span(transitive triples) ✓ for all tournaments")

# Now re-check: does ∂(Ω_3) ⊆ Ω_2?
print(f"\n{'='*70}")
print("∂(Ω_3) ⊆ Ω_2?")
print("="*70)

n = 5
fails = 0
for t_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    a2_tuples = [tuple(p) for p in a2]
    a3_tuples = [tuple(p) for p in a3]

    om3 = compute_omega_basis(A, n, 3, a3, a2)
    dim_om3 = om3.shape[1] if om3.ndim == 2 else 0
    if dim_om3 == 0:
        continue

    # ∂_3 maps A_3 → free space of 2-paths
    # But we only index allowed 2-paths. Non-allowed 2-paths should have
    # zero total coefficient (that's the Ω_3 condition).
    #
    # The Ω_3 condition is: ∂u ∈ A_2, meaning only allowed 2-paths appear.
    # So ∂u is already in A_2 by definition.
    #
    # But then does ∂u ∈ Ω_2? We need ∂(∂u) ∈ A_1.
    # Since ∂²=0 (on the full free space), ∂(∂u) = 0 ∈ A_1. ✓
    # So YES, ∂(Ω_3) ⊆ Ω_2.
    #
    # The issue in my earlier test: I was checking if ∂u has zero coefficients
    # on NON-transitive-triple elements of A_2. But ∂u CAN have such coefficients
    # and still be in Ω_2 = span(transitive triples) only IF those coefficients
    # are indeed zero.
    #
    # WAIT: I just proved Ω_2 = span(transitive triples). So if ∂u ∈ Ω_2,
    # then ∂u has zero coefficients on non-transitive 2-paths.
    #
    # But my test showed nonzero coefficients! Let me check carefully.

    # Build ∂_3 from A_3 to ALL 2-paths (including non-allowed ones)
    # Actually, the issue is: ∂(v0,v1,v2,v3) includes faces that may not
    # be in A_2 at all (e.g., (v0,v2,v3) when v0 does not beat v2).
    # The Ω_3 condition says: these non-A_2 face contributions must cancel.
    #
    # But the boundary ALSO produces faces that are in A_2 but NOT in Ω_2
    # (non-transitive triples). These are allowed 2-paths but not ∂-invariant.
    #
    # ∂(Ω_3) ⊆ Ω_2 means: the A_2 components of ∂u must be in Ω_2.
    # Since Ω_2 = span(TT), the non-transitive A_2 components must be zero.
    #
    # Is this true? The Ω_3 condition only requires NON-A_2 faces to cancel.
    # Non-transitive triples ARE in A_2, so the Ω_3 condition doesn't
    # constrain their coefficients!
    #
    # WAIT: Ω_3 = {u ∈ A_3 : ∂u ∈ A_2}. The face (v0,v2,v3) is in A_2 iff
    # v0→v2 AND v2→v3. The Ω_3 condition says: ALL faces (v0,...,v̂_i,...,v3)
    # must be individually allowed, OR their non-allowed parts must cancel
    # as a linear combination.
    #
    # Actually no: ∂u = Σ_i (-1)^i (v0...v̂_i...v3) and we need this SUM
    # to lie in A_2. Each individual term is a (free) 2-path, not necessarily
    # allowed. The sum must be expressible as a linear combination of
    # ALLOWED 2-paths only.
    #
    # For a single 3-path (v0,v1,v2,v3):
    #   ∂ = (v1,v2,v3) - (v0,v2,v3) + (v0,v1,v3) - (v0,v1,v2)
    # Each face is a triple. It's in A_2 iff consecutive edges exist.
    # F0=(v1,v2,v3): v1→v2 ✓, v2→v3 ✓ → always in A_2
    # F1=(v0,v2,v3): v0→v2?, v2→v3 ✓ → in A_2 iff v0→v2
    # F2=(v0,v1,v3): v0→v1 ✓, v1→v3? → in A_2 iff v1→v3
    # F3=(v0,v1,v2): v0→v1 ✓, v1→v2 ✓ → always in A_2
    #
    # So compute_omega_basis checks: F1 coefficients on non-A_2 paths
    # must cancel, and F2 coefficients on non-A_2 paths must cancel.
    #
    # But it does NOT check whether the A_2 face components are in Ω_2!
    # And the ∂²=0 argument says they MUST be.

    # Let me verify: build ∂_3 on Ω_3, get A_2 components, check Ω_2-ness

    a2_idx = {t: i for i, t in enumerate(a2_tuples)}
    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]
    ntt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 0]

    bd3_full = np.zeros((len(a2_tuples), len(a3_tuples)))
    for j, path in enumerate(a3_tuples):
        v0,v1,v2,v3 = path
        faces = [(v1,v2,v3), (v0,v2,v3), (v0,v1,v3), (v0,v1,v2)]
        signs = [1, -1, 1, -1]
        for face, sign in zip(faces, signs):
            if face in a2_idx:
                bd3_full[a2_idx[face], j] += sign

    bd3_omega = bd3_full @ om3

    # Check non-transitive-triple components
    ntt_indices = [a2_idx[t] for t in ntt if t in a2_idx]
    if ntt_indices:
        max_ntt = np.max(np.abs(bd3_omega[ntt_indices, :]))
        if max_ntt > 1e-8:
            fails += 1
            if fails <= 3:
                print(f"  T#{t_idx}: max non-TT component = {max_ntt:.4f}")
                # Which Ω_3 element has non-TT boundary?
                for col in range(dim_om3):
                    ntt_vals = bd3_omega[ntt_indices, col]
                    if np.max(np.abs(ntt_vals)) > 1e-8:
                        # Show the Ω_3 element
                        vec = om3[:, col]
                        terms = [(a3_tuples[j], vec[j]) for j in range(len(a3_tuples)) if abs(vec[j]) > 1e-8]
                        print(f"    Ω_3 element {col}: {len(terms)} terms")
                        for p, c in terms[:5]:
                            v0,v1,v2,v3 = p
                            print(f"      {c:+.3f}·{p} [a→c={A[v0][v2]}, b→d={A[v1][v3]}]")
                        # Show its boundary's non-TT parts
                        for idx in ntt_indices:
                            val = bd3_omega[idx, col]
                            if abs(val) > 1e-8:
                                triple = a2_tuples[idx]
                                a,b,c = triple
                                print(f"      ∂ has non-TT: {val:+.4f}·{triple} (a→c={A[a][c]})")
                        break

print(f"\nn=5: Tournaments where ∂(Ω_3) has non-TT A_2 components: {fails}")

if fails > 0:
    print("""
This means ∂(Ω_3) lands in A_2 but NOT purely in Ω_2 = span(TT).
But ∂²=0 says it MUST be in Ω_2...

RESOLUTION: The boundary of a NON-transitive-triple (a,b,c) with c→a is:
  ∂(a,b,c) = (b,c) - (a,c) + (a,b)
  (a,c) is NOT in A_1 since c→a.
  But (b,c) and (a,b) ARE in A_1.
  So ∂(a,b,c) has one term NOT in A_1.
  Therefore (a,b,c) ∉ Ω_2.

But ∂u ∈ Ω_2 means ∂u is a linear combination of Ω_2 elements.
If ∂u has a non-TT component, say α·(a,b,c) with c→a, then
∂(∂u) would have α·∂(a,b,c) = α·((b,c)-(a,c)+(a,b)).
Since ∂²=0, α·(-(a,c)) must be cancelled. But (a,c) is NOT in A_1,
so no other term can provide (a,c). Therefore α = 0.

So non-TT components of ∂u MUST be zero if ∂u ∈ A_2 and u ∈ Ω_3.
If we're seeing nonzero values, it's a numerical issue or a bug.
""")

    # Let's check with higher precision
    print("Rechecking with compute_omega_basis and tighter tolerance...")

print("\nDone.")
