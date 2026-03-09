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

## New Proof Strategy: i_*-Injectivity (kind-pasteur-S47)

**KEY DISCOVERY (HYP-380/381):** The inclusion map i_*: H_3(T\v) -> H_3(T)
is ALWAYS injective when beta_3(T\v) = 1. This eliminates the need for
good vertex existence!

### The LES Dichotomy

For any tournament T with beta_3(T) = 1 and any vertex v:

| beta_3(T\v) | rank(i_*) | H_3(T,T\v) | beta_3(T) |
|-------------|-----------|------------|-----------|
| 0 | 0 | 1 | 0 + 1 = 1 |
| 1 | 1 | 0 | 1 + 0 = 1 |

Either way, beta_3(T) = 1. This holds for ANY vertex v, not just "good" vertices.

### Computational Evidence

| n | Method | (T,v) pairs tested | Violations of dichotomy |
|---|--------|-------------------|------------------------|
| 6 | exhaustive | 1920 | 0 (all b3_Tv = 0) |
| 7 | sampled 100 | 700 | 0 (71 with b3_Tv = 1) |
| 8 | sampled 20 | 160 | 0 |

### Revised Proof Architecture

Two claims suffice (both verified):

**Claim I (i_*-injectivity):** When beta_3(T\v) = 1, the inclusion i_*: H_3(T\v) -> H_3(T)
is injective (rank = 1). Equivalently: H_3(T,T\v) = 0 whenever beta_3(T\v) = 1.

**Claim II (relative bound):** dim H_3(T,T\v) <= 1 for all tournaments T and vertices v.

Together: for n-tournament T, pick any vertex v. By induction b3(T\v) in {0,1}.
- If b3(T\v) = 0: beta_3(T) = H_3^rel <= 1 (Claim II)
- If b3(T\v) = 1: i_* injective (Claim I), so beta_3(T) = 1 + 0 = 1

Neither claim alone suffices: Claim I needs base cases, Claim II gives beta_3 <= 2.
Together they give beta_3 <= 1 for all n.

### Advantage over Good Vertex Strategy

The good vertex approach requires finding a SPECIFIC v with beta_3(T\v) = 0.
The i_*-injectivity approach works for ANY vertex. This is cleaner because:
1. No vertex selection needed
2. The proof holds vertex-by-vertex, not just existentially
3. i_*-injectivity is a structural property of path homology, not a tournament property

### Why i_* is Injective (Conjectural)

When b3(T\v) = 1, the unique H_3 generator of T\v is a "robust" cycle that
involves global tournament structure. When embedded in T (which has more edges),
it cannot become a boundary because the extra d_4 image from T's additional
Omega_4 content is "orthogonal" to the embedded cycle.

The seesaw mechanism (THM-095) may play a role: beta_1(T\v) = 0 when
beta_3(T\v) = 1, forcing im(d_2) to be large, which constrains how
d_4 can interact with the embedded cycle.

## Open Problems

1. **Prove i_*-injectivity algebraically.** (HYP-380) Why is rank(i_*) = 1
   when beta_3(T\v) = 1? This is the key open algebraic claim.

2. **Prove relative H_3 bound algebraically.** (HYP-351) Why is dim H_3(T,T\v) <= 1?
   The relative complex R_p consists of Omega chains through vertex v.

3. **Can i_*-injectivity be proved from the seesaw?** The beta_1*beta_3 = 0
   seesaw (THM-095) constrains the chain complex structure. Does it force
   i_*-injectivity?

4. **Extend to beta_5.** The seesaw mechanism (THM-098) predicts beta_5 in {0,1} too.
   Does the same dichotomy hold for i_*: H_5(T\v) -> H_5(T)?

## Files
- les_rank_i_star_v2.py — main computation (v2 = fixed mod-p arithmetic)
- les_rank_i_star_n7.py — first version (buggy Omega coords, wrong results)
- relative_h3_structure_deep.py — relative complex dimension analysis

## See Also
- THM-098 (Boolean odd Betti conjecture)
- THM-108 (beta_2 = 0 proof)
- THM-109 (good vertex existence for beta_2)
- THM-103 (beta_1 <= 1)
- HYP-349+ (various beta_3 hypotheses from S46)
- HYP-380-382 (i_*-injectivity and LES dichotomy from S47)
