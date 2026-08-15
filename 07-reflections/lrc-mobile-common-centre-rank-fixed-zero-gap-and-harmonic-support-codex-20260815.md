# Mobile common-centre rank versus the fixed-zero rank

**Date:** 2026-08-15  
**Status:** THM-3401 PROVED for the fixed source centre `t=0`;
PROVED-ELEMENTARY general zero-centre CRT upper construction;
FINITE-EXACT mobile common-centre minima for `q=15,...,28` on the literal
owner pool `{1,...,14}`; PROVED-ELEMENTARY capacity-sharp mobile rank six at
`q=23`.  The unique physical but non-mobile-common-centre rank-six edge at
`q=15` is FINITE-EXACT.  This is an unnumbered complement to THM-3401, not a
replacement for its true fixed-zero theorem.  No LRC(14) ledger decrement
follows.

## 1. The gauge that THM-3401 fixes

[THM-3398](../01-canon/theorems/THM-3398-general-finite-mode-sheet-cover-cochain.md)
proves that a physical cyclic sheet cover is exactly a selected-mode block
cover carrying one complete affine cochain.  Its centre lattice for owner
`u_i` and mode residue `h_i` is

```text
L_i=(Z+h_i/(2q))/u_i,                                 (1)
```

and its cochain is

```text
p_ij=2q u_i u_j(x_i-x_j).                            (2)
```

[THM-3401](../01-canon/theorems/THM-3401-centered-transverse-sheet-cover-rank-fifteen-through-twenty-eight.md)
proves the minimum rank when the physical source time is fixed at zero:

```text
r_0(q)=phi(q)/2+omega(q)-1_(q prime),  15<=q<=28.    (3)
```

There is a second, strictly larger zero-cochain carrier.  Define a **mobile
common-centre certificate** to be one selected mode per chosen owner such
that

```text
c belongs to every L_i for some rational c,
the selected sheet blocks cover Z/qZ.                 (4)
```

Taking every `x_i=c` makes every `p_ij=0`.  Let
`r_mob^[14](q)` be the minimum number of owners in `(4)` with owners
restricted to `{1,...,14}`.  The superscript is load-bearing.  Fixed zero is
one allowed centre, so

```text
r_mob^[14](q)<=r_0(q).                               (5)
```

The converse identification is false.  Translating a physical time `c` to
zero would translate a sheet label by `qc`, because

```text
u(c+ell/q)=u(ell+qc)/q.                              (6)
```

Only when `qc` is integral is this a common permutation of `Z/qZ`.  At the
half-sheet centres responsible for the new compression, `(6)` is not a legal
sheet relabelling.  Equivalently, `p=0` remembers equality of centre lifts
but does not force their common additive gauge `c` to be zero.  Thus
THM-3401's precise statement remains true, while any prose identifying its
fixed-zero slice with the whole zero-cochain slice needs this correction.

MISTAKE-383 is the closest corrected near miss: bounded-rank edge lists do
not determine a full physical clutter.  The canonical hostile remains the
`q=8` `4+2+2` block partition with nonempty pair fibres but no closed
cochain.  The least-used relevant sidecar is now the additive centre gauge.

The live concept board is:

1. THM-3398 selected modes and complete affine cochains;
2. THM-3401's fixed-zero quotient blocks;
3. the mobile common-centre gauge in `(4)`;
4. prime kernels versus longer symmetric blocks;
5. dilation as a degree-graded branch transplant;
6. all-tie block clutters of owner size four and six;
7. harmonic support of the strict gap `r_0-r_mob^[14]`.

## 2. The zero-centre CRT construction

Let `p` range over the distinct prime divisors of a composite `q`.  The
owner

```text
u_p=q/p                                                (7)
```

has at centre zero the kernel block

```text
K_p={ell:p divides ell}.                              (8)
```

Their union is exactly the nonunit residues.  For each unit sign class
`{+-a}`, choose a unit owner with inverse `+-a`.  Its centre-zero trimode is

```text
B_u={0,+-u^(-1)}.                                     (9)
```

The nonzero parts of `(9)` partition the units.  All selected modes have
`h=0`, every gap is zero, and `(8)`--`(9)` cover all sheets.  The rank is

```text
omega(q)+phi(q)/2             if q is composite,
phi(q)/2                      if q is prime.          (10)
```

This proves an abstract fixed-zero upper construction for every `q>=15`.
For `15<=q<=28`, all chosen owners are at most 14 and THM-3401 proves that
`(10)` is minimal.  Beyond 28 it remains only an upper construction: longer
unit blocks can beat it, as THM-3401's `q=29` cycle already shows.

## 3. Exact fixed-zero/mobile comparison

The finite atlas enumerates every mode centre of every owner, retains the
largest block at each owner-centre pair (the blocks form an inclusion chain),
and exhausts owner subsets in increasing size at every common centre.  It
therefore computes `(4)` exactly, rather than sampling source times.

| `q` | fixed-zero `r_0` | capped mobile `r_mob^[14]` | gap | canonical mobile owners |
|---:|---:|---:|---:|---|
| 15 | 6 | 6 | 0 | `(1,2,3,4,5,7)` |
| 16 | 5 | 4 | 1 | `(2,6,10,14)` |
| 17 | 8 | 8 | 0 | `(1,2,3,4,5,6,7,8)` |
| 18 | 5 | 4 | 1 | `(2,10,12,14)` |
| 19 | 9 | 9 | 0 | `(1,2,3,4,10,11,12,13,14)` |
| 20 | 6 | 6 | 0 | `(1,3,10,11,12,13)` |
| 21 | 8 | 8 | 0 | `(1,2,3,4,5,10,13,14)` |
| 22 | 7 | 6 | 1 | `(2,4,9,10,13,14)` |
| 23 | 11 | 6 | 5 | `(1,4,5,7,9,11)` |
| 24 | 6 | 6 | 0 | `(1,5,7,8,11,12)` |
| 25 | 11 | 7 | 4 | `(1,4,7,8,9,10,11)` |
| 26 | 8 | 8 | 0 | `(1,2,3,5,7,9,11,13)` |
| 27 | 10 | 9 | 1 | `(1,3,4,5,6,7,8,11,13)` |
| 28 | 8 | 8 | 0 | `(1,3,4,5,9,11,13,14)` |

Hence the mobile sequence is

```text
6,4,8,4,9,6,8,6,6,6,7,8,9,8,                       (11)
```

and mobility helps exactly at

```text
C={16,18,22,23,25,27}.                               (12)
```

This is the decisive separation between fixed-zero and zero-cochain rank.
An independent formula-level atlas also classifies the surviving sheet-gauge
twist

```text
theta=qc mod 1.                                      (12a)
```

For every cover-capable exact centre in the capped `q=15,...,28` atlas,
`theta` is either zero or `1/2`; every strict improvement in `(12)` occurs
only at `theta=1/2`.  Thus this range has an exact Boolean twist carrier, but
the block/mode sidecar is still needed to recover rank.

The owner cap cannot be dropped.  At `q=25,c=1/50`, owners
`(1,9,10,11,19,21)` give a six-owner partition, below the capped value seven;
at `q=27,c=1/54`, owners `(3,15,18,21)` give a four-owner partition, below
the capped value nine.  Both are exact hostile controls.  At fixed zero,
speed types reduce modulo `q`; at half twist they reduce only modulo `2q`, so
the representatives `1,...,14` cease to be exhaustive.

### 3.1 Dilation branch transplants

The rank-four minima are pure double lifts of earlier perfect partitions:

```text
q8 -> q16: (1,3,5,7) -> (2,6,10,14),  c=1/32,
q9 -> q18: (1,5,6,7) -> (2,10,12,14), c=1/12.       (13)
```

Their values `qc=1/2` and `qc=3/2` exhibit exactly the half-sheet obstruction
in `(6)`.  This is an exact degree-graded branch transplant under THM-3398's
dilation action.  The owner cap is load-bearing: further dilation eventually
leaves `{1,...,14}`.

### 3.2 The q=23 long-block boundary

At `q=23` every literal owner is a unit and every selected block has at most

```text
ceil(23/7)=4                                           (14)
```

sheets, so at least six owners are necessary.  At centre `c=1/2`, owners

```text
(1,4,5,7,9,11)                                       (15)
```

carry disjoint blocks of sizes `(4,3,4,4,4,4)`.  They partition all 23
sheets.  Thus, independently of the exhaustive search,

```text
r_mob(23)=6.                                         (16)
```

This is five below `r_0(23)=11`.  The canonical mobile certificates are
perfect partitions exactly for `q in {16,18,22,23}`.

### 3.3 The q=15 positive-drift boundary

The full physical `q=15` rank-six clutter has 157 edges.  Exactly 156 admit
a mobile common-centre certificate, and all 156 have gcd type

```text
(5,3,1,1,1,1).                                       (17)
```

The unique physical edge outside the mobile slice is

```text
E_*=(1,7,8,9,11,13),       gcd type (3,1,1,1,1,1).  (18)
```

One exact THM-3398 witness orders it as `(9,1,7,8,11,13)`, occurs at source
time `32939/36960`, and has a nonzero complete cochain with displayed norms

```text
sum |p_ij|=50,                 max |p_ij|=6.          (19)
```

Those are witness norms, not minima.  The exact statement is the zero versus
nonzero boundary: exhaustive centre enumeration rules out every zero
cochain, while the displayed cochain realizes the physical edge.

## 4. Tournaments, trees, and harmonic subsets

For every mobile certificate the intrinsic pair observable `(2)` is zero.
The four-owner and six-owner examples are therefore complete all-tie objects,
not ordinary tournaments.  Perfect certificates are cyclic block partitions;
forcing an orientation introduces gauge and forgets the cover predicate.
The faithful finite carrier is a block clutter decorated by owner, centre,
and mode length.

The trimodes `(9)` have ternary shape and the dilations `(13)` form rooted
branches, but `(15)` uses four-phase blocks.  A ternary ancestry word alone
cannot classify the rank sequence: mode length is a required sidecar.  This
is the same general lesson as a tournament with missing or bidirected edges:
the relation alphabet must retain ties and overlap, not collapse them into
arbitrary arcs.

Regarding every subset of the naturals as a subset of the harmonic series,
the strict-gap support `(12)` has exact masses

```text
sum_(q in C) 1/q                 = 776071/2732400,
sum_(15<=q<=28) (r_0(q)-r_mob^[14](q))/q
                                    = 1579159/2732400. (20)
```

The first mass remembers where mobility matters; the second remembers its
rank saving.  Neither scalar remembers centres, owner words, block order, or
the dilation/long-block distinction.  Those sidecars are indispensable.

## 5. Typed analogies and next frontier

[THM-3397](../01-canon/theorems/THM-3397-torsor-killing-versus-effective-boundary-valuations.md)
shows that killing a class or finite torsor need not make a specified divisor
effective.  Here killing the cochain differences does not fix the additive
centre gauge and does not by itself cover the sheet fibre.  This is a grammar
match, not a map from LRC blocks to JC divisors.

[THM-3400](../01-canon/theorems/THM-3400-discounted-norming-orbit-commutator-flux-tariff.md)
quantifies leakage when exact commutation fails.  Edge `(18)` is a first
positive-cochain control, but `(19)` is unoptimized.  The next exact invariant
is

```text
lambda_q(U)=minimum cochain norm among physical realizations of U,       (21)
```

with `lambda_q(U)=0` exactly on the mobile common-centre slice.  Optimizing
`(21)` across ranks would turn the fixed/mobile/physical trichotomy into a
genuine leakage tariff.

## 6. Reproduction and scope

Run

```text
python 04-computation/lrc_mobile_common_centre_crt_rank_probe_20260815.py
python -O 04-computation/lrc_mobile_common_centre_crt_rank_probe_20260815.py
```

The artifact is exact-rational, contains no `assert`, pins its dependencies,
checks an independently hashed full physical `q=15` clutter, and requires
ordinary and optimized transcripts to agree.  Integrity hashes are recorded
in the results index and stored output after replay.

The cheapest next decisive tests are:

1. an independent implementation of the mobile centre atlas;
2. physical rank versus `r_mob^[14]` for all `q=16,...,28`;
3. classification of perfect mobile partitions beyond `q=23`;
4. exact optimization of `(21)` on the exceptional `q=15` edge;
5. an all-owner continuation beyond 28, separated from the literal cap.

Nothing here classifies core rescues, transports a current into the refined
LRC ledger, changes THM-3401's fixed-zero ranks, or proves LRC(14).
