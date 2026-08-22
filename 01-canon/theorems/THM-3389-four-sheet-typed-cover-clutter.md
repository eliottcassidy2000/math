---
id: THM-3389
title: "Four-sheet typed cover clutter"
status: >
  PROVED analytic complete-cochain criterion + FINITE-EXACT literal q=4
  clutter/atlas + INDEPENDENTLY HOSTILE-AUDITED.  On a four-sheet fibre, odd
  transverse speeds block singleton sheets and speeds congruent to two modulo
  four block antipodal pairs.  Inclusion-minimal full covers therefore have
  block partitions 2+2, 2+1+1, or 1+1+1+1.  An exact complete affine gap
  cochain with zero triangle circulation decides common phase.  The literal
  clutter has 36 edges of ranks (2:4,3:15,4:17), independence profile
  (1,11,51,118,123,44,3,0,...), and classifies all 619 q=4 rows of THM-3387,
  with no core rescue.  This gives no refined-ledger decrement or LRC(14).
source: codex-2026-08-14-q4-typed-cover-clutter
audit: independent blocker/CRT/Helly proof, 4950 rank-two pairs through 400, 6525 rank-three cases through 60, 5985 rank-four cases through 41, exact atlas replay, endpoint hostile, dilation and harmonic audit
depends_on:
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
related:
  - THM-3388-three-sheet-phase-triangle-cover-clutter
  - THM-3366-all-sector-complement-clock-completion
script: 04-computation/lrc14_q4_typed_cover_clutter_thm3389.py
output: 05-knowledge/results/lrc14_q4_typed_cover_clutter_thm3389.out
script_sha256: cd963c20ff47c9840222c6bfd95088e3d649ac90edd518353c0f64f4f5ec9bfd
output_sha256: 45459e60a69bea9d7e99746fb1b1ad8dc8f79506bb27eb7f287f887c38270adc
semantic_sha256: 5e19b6e083e8b1faf6c1a082c6e321f148b2e5ff85bf93aa1d15c92542bf87b2
hash_basis: LF-normalized bytes
---

# THM-3389 -- four-sheet cover is a typed clutter, not a tournament

**PROVED analytic complete-cochain criterion + FINITE-EXACT literal `q=4`
clutter/atlas + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and connection contract

THM-3387 identifies full sheet cover as the exact transverse obstruction.
THM-3388's proved `q=3` analysis shows why pairwise compatibility needs
an integral phase cochain.  At `q=4`, the new feature is not a larger
tournament: blocker vertices themselves have two different sheet capacities.

| field | connection |
|---|---|
| source | transverse speeds on a degree-four fibre |
| target | a typed set-cover clutter on `Z/4Z` |
| map | send a speed to its blocked residue class and its affine centre interval |
| preserved | block size, strict overlap, cyclic sheet labels, common phase |
| destroyed by an ordinary graph | simultaneous sheet block, phase magnitude, higher circulation |
| sidecar | blocker type plus a complete affine gap cochain |
| decisive hostiles | `(2,7,11)` at rank three and order `(1,3,11,5)` at rank four pass every pair test but do not glue |

## 2. Two blocker species

Let `u` be transverse to `q=4`, and put

```text
g=gcd(u,4),                  m=4/g.                     (1)
```

Across sheets `t+k/4`, the phases `u(t+k/4)` have `m` distinct values,
each repeated `g` times.  Their spacing is `1/m`, while a danger arc has
length `1/7<1/m` for `m=2,4`.  Thus at a fixed source time, a firing speed
blocks exactly one residue class modulo `m`, of size `g`:

```text
u odd       (g=1): one singleton sheet;
u=2 mod 4   (g=2): one antipodal parity pair.           (2)
```

In an inclusion-minimal full cover, a singleton cannot lie inside a selected
pair, two pair blocks cannot coincide, and singleton blocks cannot repeat.
The selected blocks therefore partition `Z/4Z`.  The only possibilities are

```text
4=2+2,              4=2+1+1,              4=1+1+1+1.  (3)
```

Consequently minimal cover edges have ranks two, three, or four with exactly
the species patterns displayed in `(3)`.  This typing conclusion is analytic
and holds for arbitrary positive transverse speeds, not only the literal
pool.

## 3. Complete affine cochain criterion

Choose one of the partitions `(3)`.  Give each owner `i` a speed `u_i` and a
representative sheet label `k_i`: pair owners use parity labels `0,1`, and
singleton owners use labels in `Z/4Z`.  The canonical representatives are

```text
2+2:       k=(0,1),
2+1+1:     k=(0,1,3),
1+1+1+1:   k=(0,1,2,3),                              (4)
```

up to cyclic rotation and reflection.  A lifted centre has the form

```text
x_i=a_i/u_i-k_i/4,                  a_i in Z,           (5)
```

with interval radius `1/(14u_i)`.  For each oriented pair define

```text
p_ij=4u_i u_j(x_i-x_j)
    =4(u_j a_i-u_i a_j)+(k_j-k_i)u_i u_j.               (6)
```

The exact affine and strict-overlap constraints are

```text
p_ij == (k_j-k_i)u_i u_j (mod 4 gcd(u_i,u_j)),
7|p_ij|<2(u_i+u_j).                                    (7)
```

Take `p_ji=-p_ij`.  The normalized edge values

```text
delta_ij=p_ij/(4u_i u_j)                               (8)
```

are centre differences exactly when their circulation vanishes on every
triangle.  Clearing denominators gives

```text
u_h p_ij+u_i p_jh+u_j p_hi=0                           (9)
```

for all distinct `i,j,h`.

Necessity of `(7)`--`(9)` follows immediately from `(5)`--`(6)`.  Conversely,
triangle closure makes the complete rational `1`-cochain `(8)` the
coboundary of rational vertex potentials `z_i`.  The first congruence in
`(7)` says that for every pair, the two shift conditions

```text
z_i+s in (1/u_i)Z-k_i/4,
z_j+s in (1/u_j)Z-k_j/4                                (10)
```

are compatible.  Choose a common denominator `M` for the potentials and the
candidate shift.  Multiplying `(10)` by `M` produces ordinary integer
congruences modulo `M/u_i`; their pairwise compatibility is necessary and
sufficient by the generalized CRT.  Thus one shift `s` puts every `z_i+s` in
its required centre lattice.

The second inequality in `(7)` makes every pair of selected owner intervals
overlap.  For the real lifts supplied by the CRT, ordinary one-dimensional
Helly gives a common open interval.  Equivalently, the circular-arc argument
works because a pairwise-intersecting family with empty total intersection
would cover the circle, whereas these selected intervals have total length at
most `4/7<1`.  Hence all owner intervals share one source time.  We have
proved:

```text
a typed partition in (3) is a full-cover edge
iff it admits a complete affine cochain satisfying (7) and (9). (11)
```

This is an explicit `H^1` closure test on `K_r`, `r=2,3,4`: the edge fibres
are affine integral sets, and all triangle classes must vanish.

## 4. The rank-two gcd law

For a `2+2` edge write `u=2a`, `v=2b` with `a,b` odd.  The smallest absolute
gap in `(7)` is `2 gcd(u,v)`.  Thus the gap set is nonempty exactly when

```text
u+v>7 gcd(u,v).                                        (12)
```

Equivalently, after collapsing each antipodal pair to one sheet, this is
THM-3387's `q=2` gcd graph for reduced speeds `a,b`.  In the literal pool the
four pair edges are

```text
(2,14), (6,10), (6,14), (10,14).                       (13)
```

The same-looking inequality is lawful here because it is transported through
an exact parity-pair quotient, not by analogy.

## 5. Literal clutter and exact q=4 atlas

The literal transverse vertex set is

```text
V={1,2,3,5,6,7,9,10,11,13,14}.                         (14)
```

Applying `(11)` gives `36` inclusion-minimal cover edges:

```text
rank 2: 4,                rank 3: 15,                rank 4: 17. (15)
```

Their complete list is frozen in the exact output.  The independent-set
profile by cardinality is

```text
(I_0,...,I_11)=(1,11,51,118,123,44,3,0,0,0,0,0).      (16)
```

There are three independent six-sets, also frozen in the output.  A literal
six-speed `q=4` body chooses `c=1,2,3` core clocks from `{1,2,3}` and `6-c`
transverse speeds.  Therefore the globally safe row count is

```text
3I_5+3I_4+I_3=3*44+3*123+118=619.                      (17)
```

An independent exact event sweep checks all `2,541` candidate rows.  Every
row containing a clutter edge leaks outside its core danger union, so there
is no q=4 core rescue.  Thus the pointwise-exact total is also `619`, exactly
the q=4 entry in THM-3387.

## 6. Pair tests still do not suffice

For the typed rank-three assignment `(2,7,11)` with labels `(0,1,3)`, all
three affine pair-gap sets are `{-2,2}`, but no choice has zero triangle
circulation.  At rank four, assign `(1,3,11,5)` to labels `(0,1,2,3)`; all six
pair-gap sets are nonempty, yet no complete cochain satisfies `(9)`.

Thus neither a graph clique nor a tournament recovers the cover clutter.
Orienting the pairs adds a gauge but not the integral gap value or its cycle
class.  A simultaneous antipodal block is also literal data: replacing it by
two cosmetic arcs forgets their common owner.

Strictness is essential.  The smallest closed-boundary-only typed cover is
`(2,5,7)` with cochain `(-2,-2,2)`: at `t=15/28`, speeds `2` and `5` touch
their danger boundaries while speed `7` fires centrally.  Replacing `<` by
`<=` would therefore add a false rank-three edge.

Positive cochain controls are

```text
(2,14):       p_01=-4;
(10,1,9):     (p_01,p_02,p_12)=(-2,2,2);
(1,3,7,5):    (p_01,p_02,p_03,p_12,p_13,p_23)
               =(-1,-2,-1,1,2,3).                     (18)
```

## 7. Typed ternary ancestry and harmonic support

Multiplying every speed by a positive odd integer preserves blocker species.  Under
the source-time rescaling, sheet label `k` is permuted to `sk mod 4`; hence
common odd dilation preserves full-cover edges.

Start from the typed `2+1+1` edge `{1,9,10}` and multipliers `7,11,13`.  At
depth `d`, the free ternary ancestry has `3^d` words, while integer support
has `binom(d+2,2)` scales with multinomial multiplicities.  The root orbits
are disjoint by their `2`-, `3`-, and `5`-adic valuations.  With root harmonic
mass

```text
1+1/9+1/10=109/90,                                     (19)
```

the distinct-support and word-weighted masses obey

```text
1001H_d=311H_(d-1)-31H_(d-2)+H_(d-3),
1001W_(d+1)=311W_d.                                    (20)
```

Here the first recurrence holds for `d>=3`, and the second for `d>=0`.

The full orbit is a structured subset of the harmonic series with mass

```text
(109/90)(1-1/7)^(-1)(1-1/11)^(-1)(1-1/13)^(-1)
=109109/64800.                                          (21)
```

As in THM-3388's proved ternary orbit, ancestry words, exponent-lattice
support, collision multiplicity, and harmonic weight remain four distinct
representations.

## 8. Verification and scope

The standard-library companion:

- checks all `1,225` pair-type rows below `200` against `(12)` and an exact
  event sweep;
- checks `1,900` typed rank-three and `715` rank-four rows beyond the literal
  pool by independent event and cochain routes;
- independently enumerates every literal minimal edge, the entire
  independence profile, all `2,541` body rows, and absence of core rescues;
- freezes pairwise-positive closure hostiles and positive cochains of every
  rank; and
- checks sixteen ternary dilation shells, collision counts, recurrences, and
  exact harmonic mass.

An independent hostile audit then proved the blocker classification and the
CRT/Helly converse afresh; checked all `4,950` rank-two pairs through `400`,
all `6,525` typed rank-three cases through `60`, and all `5,985` odd
rank-four cases through `41`; reconstructed the `36`-edge clutter and every
one of the `2,541` atlas rows from `1,042` exact event samples; and recovered
`619` exact rows with no core rescue.  New pairwise-only hostiles
`(10,1,41)` and `(1,3,27,5)` survived beyond the companion's cutoffs.

It has no floating literal or optimization-dependent `assert`.  Reproduce
with

```text
python 04-computation/lrc14_q4_typed_cover_clutter_thm3389.py
python -O 04-computation/lrc14_q4_typed_cover_clutter_thm3389.py
```

Ordinary and optimized runs LF-normalized-byte-match the stored output.

This theorem classifies the q=4 slice of THM-3387.  It does not transport
arbitrary reflected phase, close a new refined-ledger row, physically realize
a drift tail, or prove LRC(14).
