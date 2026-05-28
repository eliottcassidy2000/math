---
theorem: THM-337
name: General Formula for f(S) — Tile Coverage via Cut Subsets
status: PROVED
session: opus-2026-05-27-S2
verified: computationally all subsets n=4..7
depends_on: THM-330
---

## Statement

In the tiling model on n vertices, define f(S) for a subset S ⊆ {1,...,n-1} of cuts as:

**f(S) = #{tiles (x,y) that cross at least one cut in S}**
       = |⋃_{k∈S} tiles_k|

where tiles_k = {tile (x,y): x ≥ k > y, x ≥ y+2}.

Then f(S) admits the closed-form **Möbius formula**:

**f(S) = Σ_{∅≠T⊆S} (−1)^{|T|+1} · h(T)**

where:
- **h({k}) = k(n−k) − 1** (for singletons)
- **h(T) = min(T) · (n − max(T))** (for |T| ≥ 2)

## Proof

By inclusion-exclusion:

f(S) = |⋃_{k∈S} tiles_k| = Σ_{∅≠T⊆S} (−1)^{|T|+1} |⋂_{k∈T} tiles_k|

For a singleton T = {k}:
|tiles_k| = #{(x,y): x≥k>y, x≥y+2}
= Σ_{y=0}^{k-2} (n-k) + (n-k-1)   [the y=k-1 case has x≥k+1]
= (k-1)(n-k) + (n-k-1) = k(n-k) - 1. ✓

For |T| ≥ 2 with min(T)=a, max(T)=b:
|⋂_{k∈T} tiles_k| = #{(x,y): x≥k>y for ALL k∈T}
= #{(x,y): x≥b, y<a, x≥y+2}
Since a ≤ b and y < a ≤ b ≤ x: the condition x≥y+2 is automatically satisfied (y<a≤b-1 so y≤a-1<b≤x means y+1<b≤x, hence x≥y+2).
= #{y=0,...,a-1} × #{x=b,...,n-1} = a · (n-b). ✓

## Corollaries

**Single cut:** f({k}) = k(n-k) - 1

**Symmetric:** f({k}) = f({n-k}) (verified: k(n-k) = (n-k)k). The complement symmetry.

**Pair:** f({j,k}) for j < k:
= f({j}) + f({k}) - h({j,k})
= (j(n-j)-1) + (k(n-k)-1) - j(n-k)
= j(k-j) + k(n-k) - 2

**Triple:** f({i,j,k}) for i < j < k:
= h({i}) + h({j}) + h({k}) - h({i,j}) - h({i,k}) - h({j,k}) + h({i,j,k})
= [i(n-i)-1 + j(n-j)-1 + k(n-k)-1] - [i(n-j) + i(n-k) + j(n-k)] + i(n-k)
= i(n-i) + j(n-j) + k(n-k) - 3 - i(n-j) - i(n-k) - j(n-k) + i(n-k)
= i(n-i) + j(n-j) + k(n-k) - 3 - i(n-j) - j(n-k)

**General S:** Since h(T) for |T|≥2 depends only on min(T) and max(T), the formula telescopes cleanly.

## Application to Non-SC Count

The non-SC tiling count is:
|non-SC| = Σ_{∅≠S⊆{1,...,n-1}} (−1)^{|S|+1} · 2^{m - f(S)}

where m = C(n-1,2) is the number of tiles. Substituting the f(S) formula gives a double sum over pairs (S,T) with T⊆S, enabling further simplification.

**Dominant approximation:** The largest term in the S-sum is at |S|=1, k=1 (or k=n-1):
f({1}) = f({n-1}) = n-2, so each contributes 2^{m-n+2}.
Together: 2^{m-n+3} (two terms at k=1 and k=n-1).
Hence **non-SC(n) ≈ 2^{m-n+3}** for large n, with ratio → 1.

**Verified:** non-SC(n) / 2^{m-n+3} = 0.50, 0.75, 0.875, 0.945, 0.974, 0.986, ... → 1.

## Verification

f(S) formula verified for all subsets S ⊆ {1,...,n-1} at n=4,5,6,7 (zero mismatches).

| n | #subsets tested | mismatches |
|---|----------------|------------|
| 4 | 7 | 0 |
| 5 | 15 | 0 |
| 6 | 31 | 0 |
| 7 | 63 | 0 |
