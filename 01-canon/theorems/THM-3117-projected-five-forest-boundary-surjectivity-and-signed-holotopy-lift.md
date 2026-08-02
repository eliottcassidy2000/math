---
id: THM-3117
title: "Projected five-forest boundary surjectivity and signed holotopy lift"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; independent audit pending.  On the
  actual K8 and K9 macro spaces of THM-3110, the projected boundary from
  five-edge forest chains onto rank-four component partitions is surjective
  over Q.  Hence every rational rank-four current, including both labelled
  product-Gamma Ewens currents, has an exact rational signed forest-cycle
  lift.  Conversely, exact one-sign face propagation proves that no lift
  supported only over the nonzero macro fibres and preserving their signs
  can be a cycle.  Signed fibre null-pairs are therefore sufficient and some
  violation of literal fibrewise positivity is necessary.  This is a
  chain-level existence/boundary theorem, not product-Gamma positivity.
source: root/gmc3000-audit-2026-08-02
audit: >
  The primary exact companion constructs the full projected boundary matrix
  implicitly and proves full row rank modulo two by deterministic bitset
  Gaussian elimination: 1701/1701 for K8 and 6951/6951 for K9, with K6/K7
  controls.  The secondary exact companion reconstructs both live Ewens
  currents and kills every literal forest variable by simultaneous one-sign
  face propagation.  Normal, optimized, stored-output, hash, and independent
  mathematical audits are still required before promotion.
depends_on:
  - THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction
related:
  - THM-3112-cycle-weighted-young-subgroup-gap-and-uniform-octopus-boundary
script: 04-computation/gmc_projected_five_forest_boundary_thm3117.py
output: 05-knowledge/results/gmc_projected_five_forest_boundary_thm3117.out
script_sha256: 689b54dde79d65c017dbceb927ed7f678a340161aa7c0bee3590d09a5e8965dd
output_sha256: 21c4e3312987775e30cd69ab65e2082b414a0d4f26cd9329bce9ed67a6c7e2d9
secondary_script: 04-computation/gmc_product_gamma_rank4_forest_holotopy.py
secondary_output: 05-knowledge/results/gmc_product_gamma_rank4_forest_holotopy.out
secondary_script_sha256: 5b40ba16ee89709bbc62212e576a192a37303940a5d6ed6bee3a8a0c3345f175
secondary_output_sha256: f50cc7852b694f8066d31a3efea4ab26bcdb91455454f791d367d2d31240a20f
hash_basis: LF-normalized bytes
---

# THM-3117 -- projected five-forest boundary surjectivity and signed holotopy lift

**PROVED CANDIDATE + VERIFIED-EXACT; independent audit pending.**

THM-3110 reduces the open arbitrary-anchored width-three product-Gamma
coefficient problem to two rational currents on rank-four set partitions of
eight and nine labelled macro packets.  The obvious uniform lift of those
currents to spanning forests is not a cycle.  There are two logically
different possible conclusions: either the macro current has a genuine
homological obstruction, or the forest fibres need internal signed
cancellation.

This theorem decides between them.  Exact signed cycle lifts always exist on
the two relevant macro spaces.  What fails is precisely the literal
fibrewise-positive gauge.

## 1. Forest chains and the projected boundary

For `N in {8,9}`, let `C_r(N;Z)` be the free abelian group on the `r`-edge
forests of the complete graph `K_N`.  Every edge set is written in
lexicographic order.  Define the ordinary simplicial deletion boundary

```text
partial_r[e_1,...,e_r]
 =sum_(j=1)^r (-1)^(j-1)[e_1,...,omit e_j,...,e_r].    (1)
```

Deletion preserves acyclicity, so `(C_*(N),partial)` is a chain complex and
`partial_(r-1) partial_r=0`.

Let `V_4(N;Z)` be the free abelian group on set partitions of `[N]` having
atomic rank four, equivalently `N-4` blocks.  There are

```text
dim V_4(8)=S(8,4)=1701,
dim V_4(9)=S(9,5)=6951.                                (2)
```

For a four-forest `F`, let `pi(F)` be its component partition and define the
orientation-forgetting projection

```text
P_4:C_4(N;Z)->V_4(N;Z),          P_4[F]=[pi(F)].       (3)
```

The load-bearing map is

```text
T_N=P_4 partial_5:C_5(N;Z)->V_4(N;Z).                 (4)
```

It is not asserted to be a chain map.  It records the five rank-four flats
obtained by deleting one edge from a five-forest, with their alternating
signs.

## 2. Exact surjectivity on the two product-Gamma spaces

Reduce the integer matrix of `T_N` modulo two.  All five alternating signs
become one.  Enumerating five-forests lexicographically and performing exact
bitset Gaussian elimination gives

```text
N        target rows       rank over F_2      columns read to full rank
6             31                31                       433
7            301               301                      3187
8           1701              1701                     13683
9           6951              6951                     44524.             (5)
```

The `N=6,7` lines are stable lower controls; the `N=8,9` lines are the live
THM-3110 spaces.  Since the reduction modulo two has full row rank, some
maximal minor of the original signed integer matrix is odd.  Its determinant
is nonzero over `Q`.  Therefore

```text
T_8 and T_9 are surjective over Q.                     (6)
```

This inference is exact: full rank over one finite field is a certificate of
full rank over characteristic zero, not a probabilistic modular heuristic.

## 3. Every macro current has a signed cycle lift

Let `W in V_4(N;Q)` be arbitrary.  By `(6)`, choose

```text
Y in C_5(N;Q)                  with T_N Y=W.           (7)
```

Put `C=partial_5 Y`.  Then

```text
partial_4 C=0,                 P_4 C=W.                (8)
```

Thus every rational rank-four macro current is the projection of an exact
rational forest cycle.  In particular `(8)` applies to both labelled Ewens
currents `W_1,W_2` in THM-3110.  Clearing denominators gives the equivalent
integral statement for an integer multiple of each current.

The proof is existential.  It gives no useful denominator, sparsity,
locality, symmetry, or positivity bound for `Y` or `C`.

## 4. Literal same-sign lifts are impossible

Surjectivity does not rescue the tempting positive construction.  Write
`F_4^*(N)` for four-forests whose component partition has nonzero THM-3110
weight.  A literal same-sign distribution would be

```text
C=sum_(F in F_4^*) sign(W(pi(F))) x_F [F],
x_F>=0,
sum_(F:pi(F)=pi)x_F=|W(pi)|.                           (9)
```

It uses no forest over a zero macro fibre and does not reverse the sign of a
nonzero fibre.  At a three-face `G`, the cycle equation is

```text
sum_(F covers G) sign(W(pi(F))) epsilon(G,F)x_F=0.     (10)
```

Whenever all live oriented signs in `(10)` agree, nonnegativity forces every
incident variable to zero.  Simultaneously iterating this exact implication
gives

```text
                         W_1 / K8                 W_2 / K9
nonzero macro fibres       480                      1620
literal variables          1440                     4860
variables killed by round  1368,72                  4038,551,244,27
remaining                     0                         0.                 (11)
```

Equation `(9)` is therefore incompatible with `partial C=0`.  Combined with
`(8)`, this gives the exact boundary:

```text
signed rational forest lift                    always exists;
literal sign-preserving nonnegative lift       never exists.             (12)
```

Consequently every lift of `W_1` or `W_2` must use at least one of the
following:

1. an opposite-sign forest coefficient inside a nonzero macro fibre; or
2. a nontrivial zero-sum packet of forests over a zero macro fibre.

Those are the signed fibre null-pairs hidden by the collected rank-four
current.

## 5. Relation to the Young-subgroup quotient

The theorem separates homological existence from order positivity.  The
signed forest complex has enough room: `(6)` says there is no remaining
linear-topological obstruction.  But the sign constraint needed for a
positive coefficient proof is destroyed by `(11)`.

THM-3112's Young-subgroup projection gap keeps a different quotient.  It
forgets the alternating tree presentation while retaining monotonicity under
partition coarsening.  A positive proof may therefore still exist after
grouping by block-size type and transporting mass upward in the refinement
order.  The present theorem neither proves nor obstructs that order-theoretic
transport.  It says only that a forest-chain proof must carry signed
intra-fibre data before taking its positive quotient.

## 6. Exact companions and scope

Run

```text
python 04-computation/gmc_projected_five_forest_boundary_thm3117.py
python -O 04-computation/gmc_projected_five_forest_boundary_thm3117.py

python 04-computation/gmc_product_gamma_rank4_forest_holotopy.py
python -O 04-computation/gmc_product_gamma_rank4_forest_holotopy.py
```

and compare with the two stored outputs named in the frontmatter.  The
primary companion proves `(5)` using exact bit arithmetic.  The secondary
companion reconstructs the live rational currents, verifies every count in
`(11)`, and checks an explicit three-face forcing witness in each bank.

This theorem does not prove the all-degree row extremum, arbitrary anchored
product-Gamma goodness, SFC in width three, GMC(2), NC2, LRC(14), JC(2), or
DC(2).  It is a sharp existence/no-positive-gauge theorem for the rank-four
forest lift only.
