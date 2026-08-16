---
id: THM-3464
title: "U-spine q=123 rank-eight coexistence and exact ZMC prefix through q=227"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Exact two-layer and
  divisor-minimum analysis gives rho_ZMC(41)=rho_ZMC(123)=8 and
  rho_ZMC(227)=9.  At q=123 the minimum is attained both by the threefold
  pullback of a primitive Q=41 packet and by a primitive mixed-order Q=123
  packet.  The first seven U-spine ranks are (6,4,7,9,8,4,9).  No endpoint
  current, spectral closure, density change, or LRC(14) consequence follows.
source: codex-2026-08-15-u-spine-q123-rank-eight
audit: >
  exact integer/Fraction two-layer banks, complete branch-and-memo set-cover
  searches, q=41 witness and threefold pullback, divisor ancestry, 271816 direct/fraction mask cells, residue
  normalization, mode-width, open-cell, multiplicity, dependency, semantic,
  AST/security, and normal/optimized replay gates; independent clean-room
  meet-in-the-middle lower bounds, geometry reconstruction, Q41 pullback,
  immutable-hash, security, documentation, and normal/optimized/stored replay
  audit PASS
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum
  - THM-3461-literal-half-twist-common-centre-lifts-and-q83-rank-nine-boundary
related:
  - THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order
script: 04-computation/lrc_u_spine_q123_q227_zmc_rank_thm3464.py
output: 05-knowledge/results/lrc_u_spine_q123_q227_zmc_rank_thm3464.out
script_sha256: 74dcb9931f64de4f0c194aa2d38a70e393c78e6525cb5f4eb439ca350764bba5
output_sha256: 122029ec19709dd4e34d89d06fcde3b3dc129523de17b88092e9f49ae6e9a8e8
semantic_sha256: 992ebd92709f39f2250bb9a21a10b6e0bb3c6e21860a3bf388d826e23b500b50
hash_basis: LF-normalized bytes
---

# THM-3464 -- U-spine q=123 rank-eight coexistence and exact ZMC prefix through q=227

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The exact companion and all internal controls pass.  A clean-room
meet-in-the-middle implementation independently reproduces every lower bound,
the Q41 pullback, geometry, immutable package, and all replay gates.

## 1. Inheritance and typed rank problem

On the distinguished parabolic Berggren branch, put

```text
c_t=2t^2+2t+1,
q_t=2c_t+1=(2t+1)^2+2.                               (1)
```

THM-3416 defines `rho_ZMC(q)` as the least number of distinct positive
**transverse** owners whose strict danger blocks cover all `q` cyclic sheets
at one THM-3398 mode centre with complete cochain zero.  THM-3405 reduces
every such family, after its active gcd is recomputed, to exactly two possible
Boolean centre layers.  THM-3416 supplies the divisor-minimum formula, while
THM-3455 classifies the U-spine only through rank seven.  THM-3461 gives the
first exact value beyond that cap, `rho_ZMC(83)=9`.

The canonical near miss is the numerical coincidence

```text
rho_ZMC(51)=7=sqrt(51-2),
rho_ZMC(83)=9=sqrt(83-2).                              (2)
```

The next U-spine label is `123`, where the same extrapolation predicts `11`.
The exact answer below is `8`.

## 2. The complete two-layer finite model

For `epsilon in {0,1}`, `0<r<q`, and `ell in Z/qZ`, write

```text
w_(r,epsilon)(ell)=r(2ell+epsilon) mod 2q,
B_q^epsilon(r)={ell:
  14 min(w_(r,epsilon)(ell),2q-w_(r,epsilon)(ell))<2q}. (3)
```

These are exactly the fixed-zero and half-twist danger masks from THM-3405.
The companion checks `(3)` against the direct rational inequality

```text
||r(ell+epsilon/2)/q||<1/14                           (4)
```

on `271,816` sheet cells for `q=3,41,123,227`.  Every noncore residue modulo
`2q` has a mask already represented by `1<=r<q`; equality deduplication keeps
the least representative.  The omitted residue `r=q` is the inadmissible
core direction: it covers every sheet in the zero layer and none in the half
layer, which is precisely why the transverse qualifier is load-bearing.

For a fixed layer, discard a mask only when it is a strict subset of another
available mask.  This preserves existence of a cover with at most `k` masks:
replace the smaller mask by its superset, deleting a duplicate if necessary.
The remaining exact set-cover search stores a state `(covered sheets, slots)`.
At each node it chooses an uncovered sheet with the fewest live coverers and
branches over all of them.  It prunes only if the union of all live masks
misses a sheet or if the sum of the `slots` largest individual gains is less
than the number of missing sheets.  Both are necessary conditions, so the
branch-and-memo search is complete in both directions.

## 3. Divisor-layer coexistence at q=123

The nontrivial divisors of `123` are `3,41,123`.  The complete banks and
cap-seven searches are:

| `Q` | layer | unique/maximal masks | states / branches / memo hits | cover through 7 |
|---:|---:|---:|---:|---:|
| 3 | 0 | 1/1 | 1 / 0 / 0 | no |
| 3 | 1 | 1/1 | 1 / 0 / 0 | no |
| 41 | 0 | 20/20 | 1 / 0 / 0 | no |
| 41 | 1 | 40/40 | 307 / 325 / 19 | no |
| 123 | 0 | 61/61 | 246 / 245 / 0 | no |
| 123 | 1 | 121/121 | 16,949 / 17,104 / 156 | no |

These all-bank exclusions are stronger than the primitive exclusions needed
in the divisor-minimum formula.  Hence both `rho_ZMC(41)>7` and
`rho_ZMC(123)>7`.  At `Q=41`, the half-layer packet

```text
(3,5,11,19,28,33,37,39)                              (4a)
```

covers all `41` sheets.  It has widths `(6,6,6,6,13,6,6,6)` and occupancy
profile `1^35 2^6`, so `rho_ZMC(41)=8`.  Its threefold pullback

```text
(9,15,33,57,84,99,111,117)                           (4b)
```

covers at `q=123` with active gcd three.  Independently, the primitive
full-modulus half-layer packet

```text
(1,40,42,81,82,83,117,122)                            (5)
```

covers every sheet, so

```text
rho_ZMC(123)=8.                                        (6)
```

The search itself rediscovers `(5)` after `82` states and `81` branches.  Its
mode widths are

```text
(4,11,39,18,123,4,18,11),                             (7)
```

and the sheet multiplicity profile is

```text
1^106 2^6 4^11.                                       (8)
```

Every owner has private sheets.  The active gcd is one; all selected modes
have common physical source centre `1/(2*123)`, and an open base-coordinate
radius `2/581` keeps the packet strict while the core remains safe.  At this
centre `A_i=r_i`, so every THM-3410 projective wedge is zero.

The quotient-order profile of `(5)`, with
`m(r)=123/gcd(123,r)`, is

```text
order 3:   1 owner,
order 41:  3 owners,
order 123: 4 owners.                                  (8a)
```

Thus the minimum rank at `q=123` is attained in two distinct divisor layers:
the inherited pullback `(4b)` from primitive `Q=41`, and the primitive
full-modulus `Q=123` family `(5)`, whose owners mix all three quotient orders.
Primitivity of one realization does not imply noninheritance of the numerical
grade.  The grade therefore does not determine ancestry; THM-3416's divisor-
minimum decomposition itself remains a disjoint minimum over primitive
quotient layers.  The fourfold overlaps in `(8)` also show why XOR or a
partition model would lose the primitive positive witness.

## 4. The prime q=227 boundary

Because `227` is prime, there is no proper divisor layer.  The zero and half
banks have `113` and `226` unique maximal masks.  Complete cap-eight searches
give:

| layer | states | branches | memo hits | cover through 8 |
|---:|---:|---:|---:|---:|
| 0 | 3,681 | 3,792 | 112 | no |
| 1 | 254,223 | 257,248 | 3,026 | no |

The stronger zero-layer cap-nine search also fails after `356,355` states,
`365,568` branches, and `9,214` memo hits.  The half layer instead finds

```text
(2,6,10,215,217,219,221,223,225)                       (9)
```

after `25,892` states.  Direct union independently checks the frozen tuple,
and therefore

```text
rho_ZMC(227)=9.                                        (10)
```

The mode widths are `(3,3,3,10,10,10,10,10,10)`, the open base radius is
`1/315`, and the multiplicity profile is

```text
1^176 2^38 3^13.                                      (11)
```

Again the family is primitive, every owner has a private sheet, the source
centre is `1/(2*227)`, the complete mode cochain vanishes, and all projective
wedges vanish.

## 5. Exact first-seven U-spine rank word

Equation `(1)` gives

```text
(q_1,...,q_7)=(11,27,51,83,123,171,227).              (12)
```

THM-3461 supplies the first four ranks `(6,4,7,9)`.  Equation `(6)` supplies
the fifth.  Since `9|171`, THM-3416 gives `rho_ZMC(171)=4`.  Equation `(10)`
supplies the seventh.  Thus the exact prefix is

```text
(rho_ZMC(q_1),...,rho_ZMC(q_7))=(6,4,7,9,8,4,9).       (13)
```

In particular `(2)` is coincidence, not a recurrence law:

```text
rho_ZMC(123)=8 != 11=sqrt(123-2).                      (14)
```

This does not alter THM-3455's period-`1683` cap-seven word or its natural
and harmonic densities.  The new ranks refine two labels outside that subset;
they do not change membership in `{t:rho_ZMC(q_t)<=7}`.

## 6. Tournament, ternary-tree, and harmonic-series boundary

At each `q`, the honest finite object is a set-cover hypergraph: vertices are
sheets and labelled hyperedges are owner masks.  A pairwise tournament on
owners records neither a three-way union nor an eight- or nine-edge cover, so
it cannot recover `(6)` or `(10)`.  A ternary word can encode each sheet's
membership state (`absent`, `private`, `multiple`), but that is a lossless
serialization only while owner labels and the full incidence matrix remain
attached.

Likewise every rank-defined subset of U-spine indices is a subset of the
natural numbers and may be sampled in the harmonic series.  A finite exact
prefix such as `(13)` has no density consequence.  Periodicity, an automaton,
or another all-index theorem is required before a logarithmic coefficient may
be inferred.

## 7. Exact companion and scope

Run from the repository root:

```bash
python 04-computation/lrc_u_spine_q123_q227_zmc_rank_thm3464.py
python -O 04-computation/lrc_u_spine_q123_q227_zmc_rank_thm3464.py
```

The standard-library companion pins four theorem dependencies, uses exact
integers and `Fraction`, contains no `assert`, file write, subprocess, dynamic
evaluation, or network path, and freezes both full mask banks by SHA-256.

The candidate concerns the exact fixed-centre/ZMC ranks at two moduli and the
derived U-spine prefix.  It does not produce an endpoint coefficient, a
relation-residue current, a bispectrum, a physical LRC row, a decrement, or a
solution of LRC(14).  LRC(14) remains open.
