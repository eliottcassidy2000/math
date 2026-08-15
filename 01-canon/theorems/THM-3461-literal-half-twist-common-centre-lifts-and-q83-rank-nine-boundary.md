---
id: THM-3461
title: "Literal half-twist common-centre lifts and the q=83 rank-nine boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Literal half-twist masks
  retain a fixed common THM-3398
  centre and zero cochain.  The q=11,27,51 packets have exact ZMC ranks
  6,4,7, while an exact two-layer search gives rho_ZMC(83)=9.  Their
  projective current collapses to one ray.  No endpoint current or LRC(14)
  consequence.
source: codex-2026-08-15-common-centre-zmc-lift-q83
audit: >
  exact integer/Fraction normal and optimized replay; direct rational mask,
  mode-width, gcd-gauge, wedge, open physical-cell, two-layer set-cover,
  ancestry, affine-lift, same-mask hostile, and CRT controls; independent
  clean-room meet-in-the-middle mathematics audit PASS; independent immutable
  package replay, dependency, hash, security, routing, and scope audit PASS
depends_on:
  - THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph
  - THM-3398-general-finite-mode-sheet-cover-cochain
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3410-projective-cochain-wedge-ray-tree-tariff-and-residue-scalar-hubs
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3429-prime-fibre-activity-descent-for-mixed-order-half-twist-seven-covers
  - THM-3453-global-literal-half-twist-cap-seven-support-classification
  - THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum
related:
  - THM-3459-rule30-ternary-intersection-factorial-truth-lift-and-keller-boundaries
script: 04-computation/lrc_literal_half_twist_common_centre_q83_rank9_thm3461.py
output: 05-knowledge/results/lrc_literal_half_twist_common_centre_q83_rank9_thm3461.out
script_sha256: 8e98c0f89a4035eeec6c4e3688d51ab546928a14fe275a32ab48ec65bcb473a4
output_sha256: 1734f1c4f38fb4cf61097cffc80eecfc6f1335816507842cfbca9a315d62ce42
semantic_sha256: 97647a6cf82b23de0ddb962b6b2b714d6ca9904408af87e001fa37d0c9b37e90
hash_basis: LF-normalized bytes
---

# THM-3461 -- literal half-twist common-centre lifts and the q=83 rank-nine boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The elementary lift, exact companion, independent clean-room `q=83`
lower-bound mathematics, and immutable packaged files have all passed audit.

## 1. The fixed centre that the rank word forgot

For `q>=2`, a transverse owner `r`, and `ell in Z/qZ`, THM-3453 uses

```text
B_(q,r)={ell: ||r(2ell+1)/(2q)||<1/14}.                (1)
```

THM-3398 uses

```text
D_(q,r)(c)={ell: ||r(c+ell/q)||<1/14}.                 (2)
```

Substitution, not a new search, gives

```text
B_(q,r)=D_(q,r)(1/(2q)).                               (3)
```

Retain the unique maximal selected block containing `(3)`.  Its centre
residue satisfies

```text
h_r == r  (mod 2q).
```

Choose the integer `n_r` with `r=2qn_r+h_r`.  The THM-3398 centre lift is

```text
x_r=(n_r+h_r/(2q))/r=1/(2q).                          (4)
```

It is the same for every owner.  Hence the complete affine mode cochain is
zero.  This is the direct repair recorded as MISTAKE-404: the compressed rank
word forgets the mode sidecar, but the literal mask itself was evaluated at a
fixed common centre.

For an active owner family `U`, put `d=gcd(U)` and `g=gcd(q,d)`.  At `(4)` the
THM-3405 scalar is

```text
a=2q d x=d,
```

so the MISTAKE-389 divisibility gate `g|a` always passes.  This does not turn
an arbitrary synchronized physical time into a mode centre; it verifies the
specific half-twist centre already present in `(1)`.

## 2. Exact packets at the first three parabolic labels

The following labelled packets cover all sheets at `(4)`:

| `q` | owners `U` | `d,g` | mode widths `w_i` | common base-`y=qx` radius | multiplicity profile |
|---:|---|---|---|---|---|
| 11 | `(1,2,3,5,7,9)` | `(1,1)` | `(4,11,4,4,4,4)` | `2/63` | `1^11` |
| 27 | `(3,15,18,21)` | `(3,3)` | `(6,6,27,6)` | `1/49` | `1^27` |
| 51 | `(1,11,12,18,23,34,35)` | `(1,1)` | `(2,2,9,9,2,51,2)` | `1/245` | `1^42 2^4 3^3 4^2` |

The first two are exact partitions.  The third is an OR cover, not an XOR
cover, although every owner has private sheets.  The companion prints every
literal mask and derives each width from

```text
w=q-7 gcd(q,r)(s-1),                                  (5)
```

where `s` is the unique short cyclic phase-block length.

Their exact zero-mode-cochain ranks are

```text
rho_ZMC(11)=6,       rho_ZMC(27)=4,       rho_ZMC(51)=7. (6)
```

For `11`, this is THM-3416's rank-six atom.  For `27`, division by three gives
THM-3416's inherited `q=9` rank-four atom `(1,5,6,7)`.  Its active gcd is three, so a
primitive ambient physical row needs the breaker owner `1`.  For `51`,
THM-3416 excludes every rank at most six and THM-3453 supplies the seven-owner
witness.  The companion also exhausts both literal Boolean layers below the
displayed ranks.

## 3. These are open physical cells, not endpoint artifacts

Use THM-3387 with core `C={1}` and the primitive rows

```text
q=11: F=(1,2,3,5,7,9,11),
q=27: F=(1,3,15,18,21,27),
q=51: F=(1,11,12,18,23,34,35,51).                    (7)
```

At base point `y=1/2`, the core is safe and the transverse owners cover the
whole `q`-sheet fibre.  The radii in the table give open intervals contained
in `B_q(U)` and disjoint from the core danger arc.  With
`D=14 lcm(F)/q`, exact non-`D`-grid controls are respectively

```text
y=193/378,       149/294,       1227/2450.             (8)
```

Thus the packets are genuine local physical obstructions to one quotient
step.  They are not thirteen-speed counterexamples: a local transverse fibre
cover does not say that a full LRC row has no globally safe time.

## 4. The fourth parabolic label has exact rank nine

The first four THM-3455 labels are

```text
(q_1,q_2,q_3,q_4)=(11,27,51,83),
```

previously typed only as ranks `(6,4,7,>7)`.  The last entry sharpens to

```text
rho_ZMC(83)=9.                                         (9)
```

Here is why the finite search is exhaustive.  Let an arbitrary active family
at prime modulus `83` have gcd `d`, and write its owners as `u_i=dv_i`.
Transversality makes `d` invertible modulo `83`.  At the two THM-3405
canonical centres, relabelling sheets by `ell -> d ell` sends the family to
the literal zero and half layers for the residues `v_i`.  It is therefore
enough to inspect the two complete banks.

After equality deduplication they contain `41` zero-layer and `82` half-layer
masks; all are inclusion-maximal.  The companion recursively stores states
`(covered sheets, slots left)`, branches on an uncovered sheet with the
fewest live owners, and applies exact remaining-mass and remaining-union
bounds.  The zero-layer cap-eight search visits `6` states and `5` branches.
The half-layer cap-eight search visits `27,992` states and `28,809` branches,
with `818` memo hits.  Neither finds a cover.  The half-layer tuple

```text
(1,13,14,27,41,42,55,69,70)                           (10)
```

covers with multiplicity profile

```text
1^66 2^12 3^5.                                        (11)
```

An independent implementation does not use this recursive solver.  In the
zero layer, all `41` masks have size `11` and share sheet zero, so eight masks
cover at most

```text
1+8*10=81<83.                                          (11a)
```

In the half layer, `41` odd masks have size `12` and omit sheet `41`, while
`41` even masks have size `11` and contain it.  Pairing the other sheets
antipodally reduces them to translates in `C_41` of

```text
A={0,8,21,24,27,31},       B={0,1,2,27,31}.            (11b)
```

Every cover needs an even mask; translation fixes it as `B_0`.  The
meet-in-the-middle audit enumerates exact union levels

```text
(1,81,3240,85320,1657445)                              (11c)
```

through four candidates.  It compares `88,603` unions of at most three with
`1,738,296` unions of at most four, performing `6,878,449` exact superset
probes; no compatible pair exists.  Independently written Python and
JavaScript routes agree.  A further immutable-package audit replayed both
companion modes, checked every dependency and frozen hash, and reconstructed
the masks, centres, widths, witnesses, physical cells, CRT lifts, and scope.

## 5. Why the current still vanishes

THM-3410 associates

```text
A_i=2q u_i x_i,          P_ij=A_i u_j-u_i A_j.         (12)
```

At `(4)`, `A_i=u_i`, so every column is a positive multiple of `(1,1)` and

```text
P_ij=0                                                   (13)
```

for every pair.  The primitive-ray quotient has one vertex; both ray-tree
tariffs are zero.  Owner labels, masks, private sheets, widths, and the
physical clock survive, but no nonzero projective current survives.

The minimal information-loss hostile is already at `q=11`.  Owners `1` and
`23` have the same half-twist mask and the same centre residue modulo `22`,
but their source interval radii are

```text
2/77          and          2/1771.                     (14)
```

Thus a mask does not determine endpoint scale or Fourier support.  This is
the same boundary seen in THM-3459: symmetric mask overlap is only the first
compiler stage; a labelled star, physical clock, and endpoint/current
sidecar remain separate.

## 6. The `7 x 13` CRT entry gate is open, but only the entry gate

Because `gcd(2q,91)=1` for `q=11,27,51`, an owner can be changed inside its
class modulo `2q` while independently imposing `7`- and `13`-adic data.  The
companion verifies the exact examples

| `q` | preserved residue | lift | `v_13` | residue modulo 7 |
|---:|---:|---:|---:|---:|
| 11 | 1 | 5083 | 1 | 1 |
| 11 | 2 | 138580 | 2 | 1 |
| 27 | 3 | 57135 | 1 | 1 |
| 27 | 15 | 215475 | 2 | 1 |
| 51 | 1 | 19279 | 1 | 1 |
| 51 | 11 | 1538069 | 2 | 1 |

Each lift preserves its literal mask, selected half-centre, and zero wedge.
This opens a valuation-decoration gate; it does not construct the missing row
relation, endpoint coefficients, delayed word, relation address, or
nonconstant `169`-twist current.  In particular it is not a `7 tensor 13`
bispectrum nonvanishing theorem.

The `q=51` affine-lift hostile must also remain attached: owners `1` and `35`
agree modulo `34`, while their lift characters modulo `3` differ.  Their
literal masks differ.  Reducing to quotient order or modulo `34` destroys the
cover information needed by THM-3429.

## 7. Exact companion and scope

Run from the repository root:

```bash
python 04-computation/lrc_literal_half_twist_common_centre_q83_rank9_thm3461.py
python -O 04-computation/lrc_literal_half_twist_common_centre_q83_rank9_thm3461.py
```

Both transcripts match the stored output after LF normalization.  The script
uses only the Python standard library, has no optimization-sensitive
assertions, pins all eight dependencies, checks `16,742` rational mask cells,
and freezes the complete two-layer banks by SHA-256.

The theorem proves a fixed-half-twist common-centre lift and an exact
`q=83` ZMC rank boundary.  It does not construct an endpoint current, a live
LRC row, a decrement, a spectral closure, a Jacobian map, or an LRC(14)
solution.  LRC(14) remains open.
