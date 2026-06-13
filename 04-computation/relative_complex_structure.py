"""
relative_complex_structure.py — Study the relative chain complex C_*(T)/C_*(T\\v)

The quotient complex C_*(T)/C_*(T\\v) has chain groups = span of through-v paths.
Its homology H_*(T, T\\v) is related to H_*(T) and H_*(T\\v) via LES.

Key insight: the "through-v-only cycles" in ker(d_3^T) that are supported
entirely on through-v paths — these are NOT the same as cycles in the
RELATIVE complex C_3(T)/C_3(T\\v), because the boundary in the quotient
complex maps differently from the boundary in C_*(T).

Wait — actually in the quotient complex, d^rel is induced by d^T:
  d^rel([σ]) = [d^T(σ)]
For a through-v chain σ, d^T(σ) has both through-v and non-through-v components.
In the quotient, the non-through-v components project to 0.
So d^rel([σ]) = [through-v part of d^T(σ)].

A through-v-only cycle z in C_3(T) has d_3^T(z) = 0.
In the quotient, [z] has d^rel([z]) = [through-v part of d_3(z)] = [0] = 0
since d_3(z) is entirely zero.

So through-v-only cycles in C_3(T) map to CYCLES in the relative complex.
But they might map to BOUNDARIES in the relative complex.

From the LES:
  H_4(T) → H_4(T,T\\v) → H_3(T\\v) → H_3(T) → H_3(T,T\\v) → H_2(T\\v)

When b4(T)=0, b2(T\\v)=0:
  0 → H_4(T,T\\v) → H_3(T\\v) → H_3(T) → H_3(T,T\\v) → 0

H_4(T,T\\v) = ker(i_*). When i_* is injective, H_4(T,T\\v) = 0.
H_3(T,T\\v) = coker(i_*). When i_* is surjective, H_3(T,T\\v) = 0.

So at n=7 (i_* always injective for BAD v):
  H_4(T,T\\v) = 0
  H_3(T,T\\v) = 0  (since b3T = b3Tv = 1 and rank(i_*)=1)

H_4(T,T\\v) = 0 means: every relative 4-cycle is a relative 5-boundary.
H_3(T,T\\v) = 0 means: every relative 3-cycle is either from H_3(T) (via j_*)
or a relative 4-boundary.

Actually, H_3(T,T\\v) = 0 means the sequence
  H_3(T) → H_3(T,T\\v) → H_2(T\\v)
has H_3(T,T\\v)=0, so j_*: H_3(T) → H_3(T,T\\v) is zero.
This means every element of H_3(T) maps to zero in relative homology.
Which means i_*: H_3(T\\v) → H_3(T) is surjective.
Combined with injectivity: i_* is an ISOMORPHISM when both work.

Now, through-v-only cycles in C_3(T) represent classes in H_3(T).
If they're not boundaries, they'd give nonzero H_3(T) classes.
But then j_*(class) should be nonzero in H_3(T,T\\v) — unless j_*
maps it to zero.

Wait — j_* is the zero map when H_3(T,T\\v)=0. So j_*(h) = 0 for
all h ∈ H_3(T). But what does j_* look like on chain level?

j_#: C_*(T) → C_*(T)/C_*(T\\v) is the quotient map.
A through-v-only cycle z maps to [z] ∈ C_3(T)/C_3(T\\v) with
d^rel([z]) = 0, so [z] represents a class in H_3(T,T\\v).

But we know H_3(T,T\\v) = 0. So [z] must be a boundary in the
relative complex: [z] = d^rel([w]) for some 4-chain w.

This means z = d_4(w) + c where c ∈ C_3(T\\v) (a non-through-v chain).
If z is through-v-only, then the non-through-v part of d_4(w) must equal -c.
And z equals the through-v part of d_4(w).

So z = π_new(d_4(w)) where π_new is the projection to through-v paths.
And d_4(w) = z + c with c purely non-through-v.

This doesn't immediately say z ∈ im(d_4). It says z differs from
d_4(w) by a non-through-v chain. So z = d_4(w) - c.

Hmm, let me reconsider. If z is a through-v-only cycle in C_3(T),
and if H_3(T,T\\v) = 0, then z is in im(d_4^rel) in the quotient complex.
This means there exists w such that d_4(w) ≡ z (mod C_3(T\\v)).
So d_4(w) = z + c for some c ∈ C_3(T\\v).

But z ∈ ker(d_3) and d_4(w) ∈ ker(d_3), so c = d_4(w) - z ∈ ker(d_3) too.
So c is a non-through-v 3-cycle.

Now: z = d_4(w) - c. If c = d_4(w') for some other 4-chain, then
z = d_4(w - w'), so z IS a boundary.

But c might not be a boundary — c could be a nonzero homology class!
In that case, z = d_4(w) - c where c represents [c] ∈ H_3(T).
But z also represents [z] ∈ H_3(T). So [z] = [d_4(w)] - [c] = -[c].
So [z] = -[c], meaning the through-v-only class equals minus the
non-through-v class.

Wait, c is a non-through-v CYCLE. If c represents a nontrivial H_3 class,
then [z] = -[c], which means both are nonzero. But the through-v-only
class [z] and the non-through-v class [c] are NEGATIVES of each other.
This means the generator has representatives in both "sides".

For z to be a boundary, we need [z] = 0, which means [c] = 0, which
means c is also a boundary. And c = d_4(w) - z, so if z is a boundary
then c is too.

This is circular. The LES argument gives:
  If H_3(T,T\\v) = 0, then j_*([z]) = 0, meaning [z] = i_*([h]) for
  some h ∈ H_3(T\\v).

But [z] is in H_3(T), and it's the image of some H_3(T\\v) class under i_*.
Since i_* is injective (rank=1), [z] = i_*(gen_{T\\v}) or [z] = 0.

If [z] = i_*(gen_{T\\v}), then z = embedded_gen + boundary.
The embedded gen is NOT through-v-only (it's supported on non-v paths).
So z - (embedded gen) is a boundary, meaning z has a non-through-v component
from the embedded gen. But z is through-v-only, contradiction!

Unless z = 0 in homology, i.e., z is a boundary. ∎

THIS IS THE PROOF!

Let me verify this logic carefully and implement a test.

Author: opus-2026-03-09-S58
"""
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

print("=" * 70)
print("PROOF OF THE GHOST CYCLE THEOREM")
print("=" * 70)
print()
print("Claim: If T is a tournament with b3(T)=1, b3(T\\v)=1,")
print("and rank(i_*)=1 (i.e., H_3(T,T\\v)=0 and H_4(T,T\\v)=0),")
print("then every through-v-only cycle in ker(d_3^T) is a boundary.")
print()
print("Proof:")
print("  Let z ∈ ker(d_3^T) be supported entirely on through-v 3-paths.")
print("  Then z represents a class [z] ∈ H_3(T).")
print()
print("  Case 1: [z] = 0 in H_3(T). Then z ∈ im(d_4^T), done. ■")
print()
print("  Case 2: [z] ≠ 0 in H_3(T). Since b3(T)=1, [z] generates H_3(T).")
print("  Since rank(i_*)=1 (i_* is an isomorphism F_p → F_p),")
print("  [z] = i_*(h) for a unique nonzero h ∈ H_3(T\\v).")
print()
print("  i_*(h) is represented by the embedded cycle: a chain in C_3(T)")
print("  supported ENTIRELY on non-through-v paths (since i embeds T\\v")
print("  paths into old paths of T).")
print()
print("  So [z] has a representative on non-through-v paths, but z itself")
print("  is supported on through-v paths. The difference")
print("  z - (embedded representative) is a boundary (both represent [z]).")
print()
print("  But z is through-v-only, and the embedded representative is")
print("  non-through-v-only. Their sum can only be a boundary if BOTH")
print("  individually are... No wait, that's not right.")
print()
print("  The issue: z - emb_rep ∈ im(d_4), but z and emb_rep have disjoint")
print("  support (different coordinate subsets). So d_4(w) = z - emb_rep")
print("  must have nonzero components on both through-v and non-through-v")
print("  coordinates. This is fine — 4-boundaries CAN mix coordinates.")
print()
print("  REVISED ARGUMENT:")
print("  [z] ≠ 0 means z ∉ im(d_4). But we assumed z is through-v-only.")
print("  Also, i_*(gen) = [z] means the embedded gen represents [z].")
print("  The embedded gen is NON-through-v-only.")
print()
print("  Now: [z] can be represented by BOTH a through-v-only chain z")
print("  AND a non-through-v-only chain (the embedded gen).")
print("  This means z - emb_gen ∈ im(d_4).")
print("  So d_4(w) = z - emb_gen for some 4-chain w.")
print()
print("  This is a consistent equation. So Case 2 IS possible in principle.")
print("  The through-v-only [z] = i_*(gen) doesn't lead to contradiction.")
print()
print("  CONCLUSION: The proof DOESN'T work via pure LES argument.")
print("  The Ghost Cycle Theorem (tv-only cycles are boundaries at n≤7)")
print("  is STRONGER than what the LES implies.")
print()
print("  The LES only tells us: either i_* is injective (and [tv-only z]")
print("  can be nonzero), or i_* is zero (and [z] could be the generator).")
print("  But the Ghost Cycle Theorem says: at n≤7, tv-only cycles are")
print("  ALWAYS boundaries, regardless of rank(i_*).")
print()
print("  This is an additional structural property of the path homology")
print("  chain complex that goes beyond the LES.")
print()
print("DONE.")
