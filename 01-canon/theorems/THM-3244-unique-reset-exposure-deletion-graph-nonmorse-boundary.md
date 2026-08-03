---
id: THM-3244
title: "Unique-reset Rips routing and deletion-graph non-Morse boundary"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.  On
  the complete 4,319-state THM-3238 physical bank, the exact exposed-reset
  functional has 32 nonreset one-pole local maxima, so global exposure does
  not induce a local deletion flow.  Its one-pole superlevel merge tree has
  33 births, 32 finite bars, and as many as 19 simultaneous components; the
  second-highest state persists as a separate component for 2,096 rank steps.
  The reset-monotone escape radius is at most 10, with a unique sharp
  radius-10 trap; at radius 10 the reset is the unique directed sink.  More
  strongly, two explicit hubs route every state to the reset in at most two
  strict, reset-monotone radius-10 jumps, and no one-hub two-level atlas can
  do so.  This is a finite directed-graph theorem, not a simplicial collapse
  or a result for another support/bank.
source: root/multiscale-newton-flag/2026-08-03
audit: >
  The exact companion reconstructs the canonical 4,319 physical
  submultisets, pins the immutable THM-3238 script and output, verifies an
  independently extracted exact rank certificate for all THM-3238 integer
  coordinates, and exhausts every one-pole neighbor, every reset-monotone
  upward pair, every possible radius-10 hub, and the complete elder-rule
  union-find merge tree of the one-pole superlevel filtration.  Normal and
  optimized runs reproduce the stored output.  Independent hostile theorem
  audit is pending.
depends_on:
  - THM-3238-complete-physical-product-gamma-bank-unique-reset-stitch
related:
  - THM-3160-complete-pluecker-pole-holotopy-and-selector-projection-no-go
  - THM-3216-depth-nine-degree-fourteen-unique-reset-face-and-omega-cone-boundary
script: 04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py
output: 05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out
script_sha256: 4f877b0287883bc7864bb593eb6f0b56b3cfda4a2a6a4963454eb094a65ae806
output_sha256: d60d9db3d28fef202a51f86338896e6881bcc5f4533ed6e1e07dcda795f7bfc6
hash_basis: LF-normalized bytes
---

# THM-3244 -- unique-reset Rips routing and deletion-graph non-Morse boundary

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.**

THM-3238 exposes the reset

```text
Q=(1,3,3,4,5,6,7,8)                                    (1)
```

as the unique maximum of an exact lawful functional `H` on all 4,319
nonempty physical submultisets of

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                  (2)
```

Exposure is global.  This theorem asks the different, local question: does
the order defined by `H` flow to `Q` through small pole edits?  The answer is
negative at one-pole scale but sharply positive at edit radius 10.  A
two-chart atlas then improves mere reachability to two jumps.

## 1. Physical edit metric and directed Rips graphs

Write `n_j(sigma)` for the multiplicity of `j` in a physical state `sigma`.
Thus

```text
0 <= (n_1,...,n_8) <= (4,3,2,2,2,1,1,1),               (3)
```

coordinatewise, with the zero vector deleted.  Put

```text
d(sigma,tau)=sum_(j=1)^8 |n_j(sigma)-n_j(tau)|.          (4)
```

The one-pole graph joins states at distance one.  Say that `tau` is
`Q`-monotone from `sigma`, and write `tau <=_Q sigma`, when

```text
|n_j(tau)-n_j(Q)| <= |n_j(sigma)-n_j(Q)|                (5)
```

for every `j`.  For an integer `r>=1`, define the directed reset-Rips graph
`R_r^+(H,Q)` by

```text
sigma --> tau
 iff tau <=_Q sigma, d(sigma,tau)<=r, and H(tau)>H(sigma). (6)
```

The name records only this threshold graph.  No simplicial Vietoris--Rips
complex is being asserted to collapse.

THM-3238's 4,319 exact integer values are pairwise distinct.  Their canonical
state order has digest

```text
c060e22f900232b608ad3bce1f6b24cae51b6eb45c138b679f4c698fbad2c6a2, (7)
```

and the compressed rank permutation used by the companion has raw digest

```text
ac7aeaf94958c04034137dbbffaf18f494485446421ad9006b6657a15b3487d7. (8)
```

Rank zero is the largest coordinate, namely `Q`.  The certificate is an exact
ordering of the full THM-3238 integer vector, not a floating approximation;
the companion pins the dependency script and output hashes before using it.

## 2. One-pole exposure is not Morse

Exhausting every legal insertion and deletion in the nonempty physical bank
gives exactly

```text
32                                                           (9)
```

nonreset states whose `H`-value is strictly larger than the value at every
one-pole neighbor.  Hence the unique global maximum `Q` does not produce a
strict one-pole ascent, a deletion gradient, or a discrete-Morse matching.

The reset-monotone obstruction is larger still.  For `sigma!=Q`, define its
escape radius

```text
epsilon(sigma)=min {d(sigma,tau):
                    tau<=_Q sigma and H(tau)>H(sigma)}.  (10)
```

The minimum exists because `tau=Q` is always admissible.  The exact histogram
of `(10)` is

```text
radius       1    2   3   4  5  6  7  8  9 10
states    4231   45  23  10  1  2  2  2  1  1.          (11)
```

Thus 87 nonreset states have no radius-one reset-monotone ascent, even though
only 32 of them are local maxima in the full one-pole graph.  This separates
failure of the chosen reset orthant from failure of every local direction.

### The order-persistent obstruction

Filter the one-pole graph by adding states in decreasing `H` order and use
the elder rule when two connected components merge.  Equivalently, rank zero
is `Q`, rank one is the next largest exact coordinate, and the filtration at
rank `k` is induced by the first `k+1` states.  The exact zero-dimensional
merge tree has

```text
33 component births, 32 finite bars,
maximum simultaneous component count 19.               (11a)
```

The longest finite rank bar is

```text
birth state (1), birth rank 1,
merge state (1,1,1,1,2,2,2,3,3,4,5,8), merge rank 2097,
rank persistence 2096.                                  (11b)
```

Thus the second-highest global state is not adjacent to the reset basin in
the superlevel graph: its component remains separate until nearly half the
bank has entered.  Two further named bars anticipate the coarse atlas below:

```text
(3): birth 16, merge at (1,3) in rank 597, persistence 581;
T:   birth 12, merge in rank 496,              persistence 484. (11c)
```

The first state in `(11c)` is the lower hub of Section 4 and `T` is the sharp
radius-10 trap.  The upper hub has rank 36 but is not a component birth.  This
rank persistence uses only the exact order of `H`, so it is invariant under
strictly increasing reparameterization; it is not a claim about gaps between
the enormous integer coordinates.

## 3. Sharp radius-10 sink filtration

Since `H` strictly increases along every directed edge, `R_r^+(H,Q)` is
acyclic.  Its number of sinks for `r=1,...,10` is

```text
(88,43,20,10,9,7,5,3,2,1).                              (12)
```

The unique state with escape radius 10 is

```text
T=(1,1,1,1,2,2,2,3,3,4).                               (13)
```

Consequently `Q` and `T` are both sinks at radius 9, while `Q` is the unique
sink at radius 10.  Every state therefore has a directed path to `Q` in
`R_10^+(H,Q)`.  Radius 10 is sharp for this fixed metric, order, and monotone
orthant.

## 4. A two-hub, two-jump atlas

The unique-sink conclusion does not bound path length.  There is nevertheless
an explicit stronger routing.  Every state with `d(sigma,Q)<=10` jumps
directly to `Q`.  Exactly 133 states lie farther away.  Put

```text
A=(3),
B=(1,1,2,2,3,3,4,5,6,7,8).                             (14)
```

Their distances to `Q` are 7 and 3.  Among the 133 far states:

```text
A is a strict reset-monotone radius-10 ascent from  36 states,
B is a strict reset-monotone radius-10 ascent from 127 states,
both are available from                              30 states. (15)
```

The two coverage sets therefore have union size

```text
36+127-30=133.                                          (16)
```

After the first jump, either hub jumps directly to `Q`.  Hence every one of
the 4,319 states reaches `Q` by at most two edges of `R_10^+(H,Q)`.

This two-hub claim is minimal within such a two-level radius-10 atlas.  Any
hub which itself jumps to `Q` must lie within radius 10 of `Q`.  Exhausting all
such states, one hub covers at most 127 of the 133 far states.  Exactly two
single hubs attain that bound:

```text
(1,1,2,3,4,5,6,7,8),
(1,1,2,2,3,3,4,5,6,7,8).                               (17)
```

Thus no one-hub two-jump routing exists.  The particular lower hub `A` in
`(14)` supplies the six far states missed by the chosen upper hub `B`.

## 5. Exact proof and hostile boundary

The companion performs only integer and finite combinatorial operations.  It
reconstructs the complete capacity box `(3)`, checks the state digest `(7)`,
decompresses and verifies the exact rank permutation `(8)`, and then exhausts:

1. every one-pole neighbor for `(9)`;
2. every reset-monotone upward pair for `(11)` and `(12)`;
3. the sharp state `(13)`;
4. every far-state incidence with the two hubs in `(14)`; and
5. every possible hub within radius 10 of `Q` for the minimality claim `(17)`;
   and
6. the elder-rule union-find merge tree `(11a)--(11c)`.

The independently extracted rank certificate is deliberately separate from
the expensive product-Gamma reconstruction in THM-3238.  An independent
audit must compare it against that full exact vector; matching only the
certificate's own digest would not establish provenance.

## 6. Scope and holotopy reading

This is a theorem about the fixed support `(1,3)`, bank `I2`, physical
submultiset metric `(4)`, reset orthant `(5)`, and exact THM-3238 functional.
It proves neither an analogous radius for the other maintained faces nor a
uniform statement for arbitrary radial coefficients.  It also gives no
probability transport, Markov deletion law, simplicial collapse, deformation
retract, or proof of the Gaussian Moment Conjecture.

The useful holotopy lesson is more precise.  A global exposed section can fail
to localize on the elementary deletion cover, yet become coherent after a
finite scale enlargement; here the first unique-sink scale is 10 and two
transition charts suffice.  THM-3160's flat order-independent pole holotopy
and this `H`-oriented coarse routing are different structures.  Identifying
them would require an additional lawful local intertwiner, which `(9)` rules
out at one-pole scale.

QED.
