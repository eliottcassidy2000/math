---
id: THM-3469
title: "Three-times-p half-twist eight-owner cover boundary and periodic exact-rank family"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
  AUDIT PENDING.  For p>=13 and q=3p, one explicit eight-owner half-twist
  template covers exactly when p is not 7 or 35 modulo 42.  On p=14k-1 it
  gives rank at most eight, and THM-3455 sharpens this to an exact periodic
  rank-4/6/7/8 word with minimal period 24322155 and exact natural/harmonic
  densities.  No endpoint current or LRC(14) consequence follows.
source: codex-2026-08-15-q123-q291-affine-template
audit: >
  self-contained nearest-multiple and strict septimal-gap proof; exact
  symbolic channel, threshold, 4325328-sheet boundary, rational-mask,
  mode-centre, hostile/repair, CRT-period, rank-count, U-spine, dependency,
  semantic, security, and normal/optimized replay gates; independent audit
  pending
depends_on:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3455-berggren-q-spine-cap-seven-atom-sieve-and-fibonacci-rank-spectrum
related:
  - THM-3461-literal-half-twist-common-centre-lifts-and-q83-rank-nine-boundary
  - THM-3464-u-spine-q123-rank-eight-break-and-divisor-layer-certificate
script: 04-computation/lrc_three_p_half_twist_eight_owner_template_thm3469.py
output: 05-knowledge/results/lrc_three_p_half_twist_eight_owner_template_thm3469.out
script_sha256: 436a9a97d275ad250cca333afb187fd4c2bd8ef2bd2b6e9641b938acdc6d6c17
output_sha256: 64c968ec6e1694f4dd45ad47420036eb949c63759d8c5c9de0e26fd3586e3dc4
semantic_sha256: 42c2ecaa420914e22e57a9415efd5dce73ec6de026dfec5bd06e19d24a2582dc
hash_basis: LF-normalized bytes
---

# THM-3469 -- three-times-p half-twist eight-owner cover boundary

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT / INDEPENDENT
AUDIT PENDING.**

The proof and deterministic companion pass their internal gates.  This file
remains outside the proved dependency graph until an independent audit checks
the universal strict-boundary argument and immutable package.

## 1. The affine owner template

Let `p>=13`, put `q=3p`, and choose the eight distinct transverse owners

```text
U_p=(1,p-1,p+1,2p-1,2p,2p+1,3p-6,3p-1).             (1)
```

At the half-twist centre `c=1/(2q)`, the literal danger mask of owner `r` is

```text
B_q(r)={ell in Z/qZ:
  14 dist_(2q)(r(2ell+1),0)<2q}.                     (2)
```

Every mask in `(1)` is nonempty.  Since owner one makes the active gcd one,
THM-3405 identifies `(2)` with selected THM-3398 modes at one common source
centre and complete cochain zero.  The statement below is about this fixed
zero-cochain slice, not arbitrary synchronized physical times.

## 2. Nearest-multiple proof of the exact cover boundary

As `ell` runs through `Z/qZ`,

```text
x=2ell+1                                                (3)
```

runs through the odd residues modulo `2q=6p`.  The greatest integer numerator
distance allowed by the strict inequality in `(2)` is

```text
h=floor((3p-1)/7),
7h<3p<=7(h+1).                                        (4)
```

If `3|x`, owner `2p` kills the phase exactly:

```text
2p x == 0  (mod 6p).                                  (5)
```

Assume now `x==sigma (mod 6)` with `sigma in {+1,-1}`.  Choose a nearest
multiple `mp` of `p` (either choice works at a tie) and write

```text
x=mp+y,        m in Z/6Z,        |y|<=floor(p/2).     (6)
```

For an owner `r=ap+b`, `b in {+1,-1}`, one has

```text
r x == (a sigma+b m)p+b y  (mod 6p).                  (7)
```

The following table makes the coefficient of `p` vanish modulo six.  Thus the
chosen owner's cyclic numerator distance is exactly `|y|`.

| nearest multiple `m` | `sigma=+1` | `sigma=-1` |
|---:|---|---|
| 0 | `1` | `1` |
| 1 | `p-1` | `p+1` |
| 2 | `2p-1` | `2p+1` |
| 3 | `3p-1` | `3p-1` |
| 4 | `2p+1` | `2p-1` |
| 5 | `p+1` | `p-1` |

Consequently the first six nonbackbone channels cover every phase with
`|y|<=h`.

It remains to treat a gap phase with `|y|>=h+1`.  Since `x` is odd,

```text
(3p-6)x == 3p-6y  (mod 6p).                           (8)
```

For `0<=|y|<=p/2`, its cyclic distance is `3p-6|y|`.  At the first possible
gap this is at most `h` exactly when

```text
3p-6<=7h.                                             (9)
```

If `7` does not divide `p`, `(4)` makes `(9)` automatic.  If `p=7s`, then

```text
h=3s-1,                                               (10)
```

and only the single remainder `|y|=3s` can fail: both the nearest-multiple
channel and `(8)` have distance `h+1`.

- If `s` is even, `p`, every multiple `mp`, and `3s` are even, so this
  remainder cannot occur for the odd phase `x`.
- If `3|s`, every phase with this remainder is divisible by three and is
  already covered by `(5)`.
- If `s` is odd and `3` does not divide `s`, the phase
  `x=2p+3s=17s` is odd and not divisible by three.  It is missed by every
  owner in `(1)`.

Therefore the union in `(2)` is all of `Z/qZ` exactly when

```text
p mod 42 not in {7,35}.                               (11)
```

The boundary is strict and sharp.  At `p=35`, `x=85` (sheet `42`) is the
minimal displayed hostile.  At `p=14`, parity removes the exceptional phase;
at `p=21`, the order-three backbone absorbs it.

For completeness, the eighth mask is genuinely active for every `p>=13`:
choose an odd `x` nearest `p/2`.  Its distance under `(8)` is at most six,
which is at most `h` from `p=16` onward; the cases `p=13,14,15` are direct.
The other masks have a channel in the table with a phase within three of its
target multiple.

## 3. An exact rank family on p=14k-1

Set

```text
p=14k-1,      q=3p=42k-3,      k>=1.                 (12)
```

Then `p mod 42` is one of `13,27,41`, so `(11)` always gives

```text
rho_ZMC(42k-3)<=8.                                    (13)
```

THM-3455's complete cap-seven atom sieve makes `(13)` exact.  After rank
priority is applied,

```text
rho_ZMC(42k-3)=4
  iff k==2 (mod 3);

rho_ZMC(42k-3)=6
  iff not rank 4 and
      [k==4 (mod 5) or k==4 (mod 11) or k==5 (mod 23)];

rho_ZMC(42k-3)=7
  iff neither lower rank and
      [k==1 (mod 13) or k==11 (mod 17) or k==27 (mod 29)];

rho_ZMC(42k-3)=8
  otherwise.                                          (14)
```

Rank five never occurs.  The rank-six atom `25` adds no class: its condition
`k==9 (mod 25)` is already contained in `k==4 (mod 5)`.  The odd rank-seven
atom `51` reduces to the modulus-17 clause in `(14)`.

The exact rank word in `k` has minimal period

```text
P=3*5*11*13*17*23*29=24,322,155.                     (15)
```

Its census is

| exact rank | count in one period | density / harmonic coefficient |
|---:|---:|---:|
| 4 | `8,107,385` | `1/3` |
| 5 | `0` | `0` |
| 6 | `4,934,930` | `14/69` |
| 7 | `1,818,080` | `33056/442221` |
| 8 | `9,461,760` | `57344/147407` |

The counts follow directly by independent CRT coordinates.  For example the
rank-eight count is

```text
2*4*10*12*16*22*28=9,461,760.                         (16)
```

Every prime in `(15)` is essential: impose its one triggering residue and
zero at the other six primes.  Adding `P` divided by that prime preserves all
other coordinates and removes the trigger, changing rank `4`, `6`, or `7` to
rank `8`.  Hence no maximal proper divisor, and therefore no proper divisor,
is a period.

Periodicity gives both natural and logarithmic densities.  For each rank `r`,

```text
#{k<=N:rho_ZMC(42k-3)=r}=delta_r N+O(1),
sum_(k<=N,rho_ZMC(42k-3)=r) 1/k=delta_r log N+O(1).   (17)
```

This is an honest subset-of-the-harmonic-series theorem in the family index
`k`; it is not the convergent reciprocal sum over the quadratic labels.

## 4. Intersection with the parabolic Berggren U-spine

For

```text
q_t=(2t+1)^2+2,                                       (18)
```

one has `q_t=42k-3` for an integer `k` exactly when

```text
t mod 21 in {5,8,12,15}.                              (19)
```

Thus the template supplies rank at most eight on a periodic set of U-spine
indices with natural and harmonic coefficient `4/21`.  The first four
representatives are

```text
(t,k,q_t,rho_ZMC)=(5,3,123,8),
                    (8,7,291,8),
                    (12,15,627,6),
                    (15,23,963,4).                    (20)
```

In particular THM-3464's `q=123` witness was the first instance of a general
affine owner template, and the next unresolved U-spine label `q=291` also has
exact rank eight.  The periodic membership `(19)` does not say that every one
of its labels has rank eight; the lower-rank clauses in `(14)` remain active.

## 5. Representation boundary

The proof is not a tournament argument.  Its exact carrier is a labelled
cover hypergraph with one order-three residue-class backbone, six
nearest-multiple channels, and one septimal gap repair.  Pairwise owner arcs
do not determine their union.  A ternary encoding can record
`uncovered/private/multiple`, but owner labels and the full incidence matrix
are required to reconstruct the certificate.

At the fixed centre, the complete mode cochain vanishes.  This does not create
an endpoint coefficient, relation-residue current, bispectrum, physical LRC
row, or decrement.  LRC(14) remains open.

## 6. Exact companion

Run from the repository root:

```bash
python 04-computation/lrc_three_p_half_twist_eight_owner_template_thm3469.py
python -O 04-computation/lrc_three_p_half_twist_eight_owner_template_thm3469.py
```

The standard-library companion checks the symbolic twelve-channel table and
all `42` strict residue classes, `4,325,328` direct sheet incidences for
`13<=p<=600`, `42,048` independent rational-mask cells, all mode centres and
widths for a complete residue window, the sharp hostile and both repairs,
`100,000` rank-classifier values, exact CRT counts/minimality witnesses, and
the U-spine intersection.  It uses explicit exceptions under `-O`, pins both
dependencies, freezes a semantic digest, and performs no file write, dynamic
evaluation, subprocess, or network action.
