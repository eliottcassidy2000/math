# Complete Proof: β₂ = 0 for Tournaments on n ≤ 5 Vertices

## Theorem
For every tournament T on n ≤ 5 vertices, β₂(T) = 0 in GLMY path homology.

## Proof

### Base cases (n ≤ 4)
Verified exhaustively:
- n = 3: 4 tournaments, all have β₂ = 0.
- n = 4: 64 tournaments, all have β₂ = 0.

### The n = 5 case

By strong induction, assume β₂(S) = 0 for all tournaments S on ≤ 4 vertices.

**Step 1.** We show: if β₁(T) = 0, then ∃v with β₁(T\v) = 0.

At n = 4: β₁(S) > 0 iff t₃(S) = 2 (the unique β₁ = 1 case at n=4).

Suppose β₁(T) = 0 but β₁(T\v) > 0 for all v ∈ V(T).
Then t₃(T\v) = 2 for all v.
Since Σ_v t₃(T\v) = 2 · t₃(T) (each 3-cycle is in n-3 = 2 subtournaments):
  5 · 2 = 2 · t₃(T), so t₃(T) = 5.

Now t₃ = C(5,3) - Σ_v C(d⁺(v), 2). With t₃ = 5:
  Σ C(d⁺(v), 2) = 10 - 5 = 5 and Σ d⁺(v) = C(5,2) = 10.
The unique solution is d⁺(v) = 2 for all v (regular tournament).

At n = 5, every regular tournament satisfies β₁ = 1 (verified: 24 regular
tournaments, all have the 5-cycle as a non-boundary 1-cycle).

This contradicts β₁(T) = 0. So ∃v₀ with β₁(T\v₀) = 0. ∎ (Step 1)

**Step 2.** For the vertex v₀ from Step 1, the LES of (T, T\v₀) gives:

  H₂(T\v₀) → H₂(T) →^{j*} H₂(T, T\v₀) →^{δ} H₁(T\v₀) →^{i*} H₁(T) → ...

By induction: β₂(T\v₀) = 0, so j* is injective.

Since β₁(T) = 0: the map i*: H₁(T\v₀) → H₁(T) = 0 is zero.
So ker(i*) = H₁(T\v₀). By exactness: im(δ) = ker(i*) = H₁(T\v₀).
Since β₁(T\v₀) = 0: im(δ) = 0, so δ = 0.
Then ker(δ) = H₂(T, T\v₀) = im(j*) = H₂(T).
So h₂_rel = β₂(T). But h₂_rel = β₂(T) + β₁(T\v₀) = β₂(T) + 0 = β₂(T).

Wait — we need h₂_rel = 0. This follows from β₁(T\v₀) = 0 and the
direct computation of the relative homology. But actually, the formula
h₂_rel = β₂(T) + β₁(T\v₀) holds, and we WANT β₂(T) = 0. With β₁(T\v₀) = 0:
h₂_rel = β₂(T).

So we need h₂_rel(v₀) = 0 to conclude β₂(T) = 0. But h₂_rel = β₂(T)!

This is circular. Let me fix the argument.

**Corrected Step 2 (Case β₁(T) = 0):**

We handle the case β₁(T) = 0 separately.

Every 1-cycle z ∈ Z₁(T) is a boundary (β₁ = 0). In particular, every 1-cycle
in T\v₀ is a boundary in T. Since β₁(T\v₀) = 0, every 1-cycle in T\v₀ is
ALSO a boundary in T\v₀. So no "filling from v₀" is needed.

The relative chain complex Ω_*(T)/Ω_*(T\v₀) has:
  h₂_rel = ker(∂₂^rel)/im(∂₃^rel)

We compute h₂_rel directly from the relative complex.
Since β₁(T\v₀) = 0: the LES gives h₂_rel = β₂(T) + 0 = β₂(T).
So h₂_rel = β₂(T), and we haven't gained anything.

Hmm, this approach is genuinely circular for β₁ = 0.

**Step 2 (Case β₁(T) > 0):**

When β₁(T) > 0: verified (HYP-263) that h₂_rel(v) = 0 for ALL v.
This gives β₂(T) ≤ h₂_rel(v) = 0, so β₂(T) = 0.

But HYP-263 at n=5 is verified computationally, not proved algebraically.

**Status:** The n = 5 proof is INCOMPLETE. Step 1 is proved. Step 2 needs
either (a) a direct proof of h₂_rel(v₀) = 0 from the relative complex,
or (b) a proof of HYP-263 at n=5. Both remain open.

## What's actually proved

The algebraic argument proves: β₁(T) = 0 ⟹ ∃v with β₁(T\v) = 0 at n=5.

Combined with exhaustive verification of β₂ = 0 at n ≤ 6, this gives a
partially algebraic, partially computational proof.

## Key gap: proving h₂_rel(v) = 0 for some v

The LES gives h₂_rel(v) = β₂(T) + ker(i*_v). For the proof, we need to
show h₂_rel(v) = 0 for some v. But h₂_rel = β₂(T) + ker(i*_v), and
we can't assume β₂(T) = 0.

UNLESS we can compute h₂_rel directly from the relative complex and show
it's 0 for some v. This is what HYP-258 asserts (verified computationally).

The WEAKER approach: if h₂_rel ∈ {0,1} (HYP-257) and Σ h₂_rel ≤ 3 (HYP-262),
then for n ≥ 5: Σ ≤ 3 < 5, so some v has h₂_rel = 0.

Author: opus-2026-03-08-S49
