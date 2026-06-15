# HYP-2530 - Moment positivity is only the one-particle shadow; the real wall is compatibility

**Status:** OPEN synthesis; exact `n=6` baby-Hodge hole certificate +
Faulhaber odd-moment Hankel positivity.
**Source:** monad-explorer-2026-06-15-S10.
**Companions:** HYP-2458, HYP-2457, HYP-2454, HYP-2529, OPEN-Q-099,
OPEN-Q-101.
**Computation:** `04-computation/baby_hodge_compatibility_wall_monad_s10.py`;
stored output
`05-knowledge/results/baby_hodge_compatibility_wall_monad_s10.out`.

## Statement

Recent mainline work reframes tournament invariant realizability as a truncated
moment problem ("baby Hodge"), while the local triangular-tower packet reframes
power-balance anchors through odd Faulhaber moments.  The combined evidence now
supports one shared warning:

```text
moment positivity / convex feasibility is only the one-particle shadow;
the true obstruction sits in compatibility / packet / integrality data.
```

Two exact local certificates make this concrete.

1. **Tournament side.**  In the exhaustive `n=6` `(c3,c5)` atlas, the flagship
   hole `(8,10)` is *moment-interior*: it lies in the same `c3=8` fiber between
   realized spectral points `(8,8)` and `(8,11)`, with the exact convex
   certificate

   ```text
   (8,10) = (1/3)*(8,8) + (2/3)*(8,11).
   ```

   So the simplest spectral/moment relaxation allows the point.  The failure is
   not at the one-particle level.

2. **Faulhaber side.**  For fixed `n`,

   ```text
   S_{2r+1}(n) = sum_{i=1..n} i*(i^2)^r,
   ```

   so the odd Faulhaber moments form a Stieltjes moment sequence on atoms
   `x=i^2` with positive weights `i`.  Their Hankel matrices are therefore
   positive semidefinite, and the stored determinants are strictly positive up
   to the atomic rank and zero afterwards.

Hence neither the baby-Hodge hole nor the `p>=3` triangular-tower drift is
caused by failure of the primary moment cone itself.  The obstruction is the
extra retained layer:

```text
tournaments -> conflict/compatibility packets (THM-499's D = alpha_2 layer),
Faulhaber   -> Bernoulli / integrality / support packet,
repunits    -> prime-length atom supply beyond a scalar length shadow.
```

## Exact Evidence

### 1. The flagship baby-Hodge hole is convexly feasible

The script enumerates every `n=6` tournament and records all interior holes in
the `(c3,c5)` plane.  Each hole comes with a same-fiber interval certificate.
For the flagship global gap:

```text
(c3,c5) = (8,10)
  = (1/3)*(8,8) + (2/3)*(8,11).
```

Representative witnesses:

```text
(8,8):  score=(2,2,2,3,3,3), H=41, D=2
(8,11): score=(2,2,2,3,3,3), H=43, D=1
```

So even the same score sequence can realize the neighboring spectral points
with different compatibility data `D=(H-1-2(c3+c5))/4`.  The hole survives the
spectral shadow and is cut only after the extra packet variable is retained.

### 2. Odd Faulhaber moments already satisfy Hankel positivity

The script records leading Hankel determinants for
`(S_1,S_3,S_5,...)`:

```text
n=3: [6, 360, 86400, 0]
n=4: [10, 3000, 15120000, 548674560000]
n=6: [21, 61740, 17425497600, 341461686730752000]
```

This is exactly what a positive atomic Stieltjes measure predicts.  Therefore
the failure of exact integer towers past `p=2` is not a positivity failure.  It
begins only when the Bernoulli / integrality layer is asked to realize the
moment shadow.

## Interpretation

HYP-2458 already warned that odd Faulhaber moments should be treated like the
scalar odd-atom inventory in the OCF: they are not the full proof object unless
the compatibility packets are retained.  HYP-2530 supplies the clean finite
certificate for that warning.

The same lesson also matches HYP-2529.  There, the scalar shadow is "block
length `L`", but exact irreducibility supply is controlled by atomic prime
lengths and Cohn/Mersenne hits.  Again the real obstruction sits one layer
above the naive scalar quotient.

## Tournament Analysis

Vertices were taken to be proof carriers rather than tournaments or moments:

```text
spectral_convex_witness,
alpha2_compatibility_wall,
faulhaber_hankel_positivity,
repunit_prime_atom_supply.
```

Pairwise observable:

```text
majority(exact_local_certificate,
         cross_thread_reach,
         obstruction_sharpness,
         scalable_proof_leverage).
```

Tie path: the listed order.

Stored fingerprint:

```text
score_hist = {0:1, 1:1, 2:1, 3:1}
directed_3cycles = 0
scc_sizes = [1,1,1,1]
hamiltonian_paths = 1
leader = alpha2_compatibility_wall
```

This pilot is transitive, which is itself informative: the compatibility-wall
carrier dominates because it is the exact place where both the tournament hole
and the Faulhaber drift stop being explainable by moment positivity alone.

## Assumption Challenge

Rejected assumption:

```text
if a numerical point survives the natural moment-positivity constraints,
then it is probably realizable.
```

The `c5=10` hole refutes that already at `n=6`, and the odd-Faulhaber Hankel
positivity shows the same issue on the triangular side.

Alternate vertex sets considered: raw tournaments, raw moments, score
sequences, Hankel minors, overlap lengths, and individual obstruction
examples.  The chosen quotient preserves the location of the obstruction wall
and destroys irrelevant low-level adjacency detail.

## Next Moves

1. Upgrade the convex baby-Hodge witness for `(8,10)` into a genuine
   flag-algebra / PSD / compatibility inequality.
2. Identify the explicit Faulhaber-side packet variable that plays the role of
   `D = alpha_2`.
3. Compare moment-feasible but atom-forbidden overlap lengths on the
   HYP-2529 repunit side.
