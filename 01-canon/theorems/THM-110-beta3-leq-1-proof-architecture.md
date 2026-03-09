# THM-110: beta_3 <= 1 for All Tournaments — Proof Architecture

**Status:** CONJECTURE (computational evidence complete; algebraic ingredients identified but not proved)
**Filed by:** kind-pasteur-2026-03-09-S46
**Depends on:** THM-098, THM-108 (beta_2=0), LES of pair (T, T\v)

## Main Theorem (THM-098 for k=1)

For every tournament T on n >= 3 vertices, beta_3(T) <= 1 in GLMY path homology.

Equivalently: rank(d_4) >= ker(d_3) - 1 (rank near-saturation).

## Computational Evidence

| n | Method | beta_3 range | Confirmed |
|---|--------|-------------|-----------|
| 3-5 | exhaustive | beta_3 = 0 always | YES |
| 6 | exhaustive (32768) | beta_3 in {0,1}; 320 have beta_3=1 | YES |
| 7 | sampled (2000) | beta_3 in {0,1}; ~8.3% have beta_3=1 | YES |
| 8 | sampled (200) | beta_3 in {0,1} | YES |

Zero violations found. Gap = rank(d_4) - (ker(d_3) - 1) is always 0 or 1.

## Proof Strategy: LES Induction

**Base case:** n <= 5, beta_3 = 0 always. DONE (exhaustive).

**Induction step:** Assume beta_3 <= 1 for all tournaments on < n vertices.
For an n-vertex tournament T, use the long exact sequence of the pair (T, T\v):

    ... -> H_3(T\v) -> H_3(T) -> H_3(T, T\v) -> H_2(T\v) = 0

Since H_2(T\v) = 0 (THM-108), the map H_3(T) -> H_3(T, T\v) is surjective.
From LES exactness:

    beta_3(T) <= beta_3(T\v) + dim H_3(T, T\v)

**Strategy:** Find v with beta_3(T\v) = 0 ("good vertex for beta_3").
Then: beta_3(T) = dim H_3(T, T\v).
If dim H_3(T, T\v) <= 1, we're done.

## Key Ingredients (Both Computationally Verified)

### Ingredient 1: Good Vertex Existence for beta_3

**Claim:** For every tournament T with beta_3(T) > 0, there exists v with beta_3(T\v) = 0.

| n | Method | beta_3>0 tournaments | All have good vertex |
|---|--------|---------------------|---------------------|
| 6 | exhaustive | 320 | 320/320 = 100% |
| 7 | sampled (500) | 34 | 34/34 = 100% |
| 8 | sampled (200) | 31 | 31/31 = 100% |

At n=6: beta_3 is COMPLETELY fragile — ALL 6 deletions give beta_3=0.
At n=7: typically 5-7 out of 7 deletions give beta_3=0.
At n=8: typically 4-7 out of 8 deletions give beta_3=0.

### Ingredient 2: Relative Homology Bound

**Claim:** dim H_3(T, T\v) <= 1 for all tournaments T and all vertices v.

| n | Method | Pairs tested | Max dim H_3(T,T\v) |
|---|--------|-------------|---------------------|
| 6 | exhaustive (beta_3>0) | 1920 | 1 |
| 7 | sampled | 147 | 1 |

At n=6: ALL 1920 values of dim H_3(T, T\v) are EXACTLY 1 (for beta_3>0 tournaments).
At n=7: dim H_3(T, T\v) in {0, 1}.

### LES Isomorphism Check

When beta_3(T\v) = 0 and H_2(T\v) = 0:
    H_3(T) ≅ H_3(T, T\v) (by LES exactness)

Verified at n=6: 1920/1920 = 100% match.

## Structural Observations

### beta_3 > 0 Tournament Types (n=6)

Two types:
- **Type A** (80 tournaments): Score (1,1,1,4,4,4). Omega_4 = 0. ker(d_3) = 1, beta_3 = 1.
  H_3 generator uses 9 paths on 9 vertex 4-subsets (not all C(6,4)=15).
- **Type B** (240 tournaments): Score (2,2,2,3,3,3). ker(d_3) = 7, rank(d_4) = 6, beta_3 = 1.
  H_3 generator uses 36 paths on all 15=C(6,4) vertex 4-subsets.

### ker(d_3) Range at n=7

ker(d_3) ranges from 10 to 46 at n=7 (for beta_3>0 tournaments).
When beta_3=1, rank(d_4) = ker(d_3) - 1 always.

### Quotient Proportionality

ALL ker(d_3) basis vectors project to PROPORTIONAL elements in H_3.
This is equivalent to dim(H_3) <= 1 and directly implies rank near-saturation.
Verified exhaustively at n=6 (240/240 Type B).

### Cokernel Direction

The cokernel direction (H_3 generator, expressed in ker(d_3) basis)
varies across tournaments — it is NOT a universal direction.
The algebraic mechanism for proportionality is tournament-dependent.

## Comparison with beta_2 = 0 Proof (THM-108)

| Feature | beta_2 = 0 | beta_3 <= 1 |
|---------|-----------|-------------|
| LES uses | H_2(T\v)=0 by induction | H_2(T\v)=0 (THM-108) |
| Good vertex | b_1(T\v) <= b_1(T) | beta_3(T\v) = 0 |
| Relative bound | H_2(T,T\v) = 0 iff injectivity | H_3(T,T\v) <= 1 |
| Fragility | Complete at n=6 | Complete at n=6 |
| Algebraic proof | YES (THM-108+109) | NOT YET |

## Open Problems

1. **Prove good vertex existence algebraically.** What property of v
   ensures beta_3(T\v) = 0? Analogy: for beta_2, "good" means b_1(T\v) <= b_1(T).

2. **Prove relative H_3 bound algebraically.** Why is dim H_3(T, T\v) <= 1?
   The relative complex C_*^rel consists of Omega chains through vertex v.

3. **Prove quotient proportionality directly.** Why do all ker(d_3) elements
   project proportionally onto the quotient ker(d_3)/im(d_4)?

4. **Extend to beta_5.** The seesaw mechanism (THM-098) predicts beta_5 in {0,1} too.

## See Also
- THM-098 (Boolean odd Betti conjecture)
- THM-108 (beta_2 = 0 proof)
- THM-109 (good vertex existence for beta_2)
- THM-103 (beta_1 <= 1)
- HYP-349+ (various beta_3 hypotheses from S46)
