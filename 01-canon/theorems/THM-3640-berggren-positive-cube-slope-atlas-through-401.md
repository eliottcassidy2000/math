---
id: THM-3640
title: "Berggren positive two-cube slope atlas through denominator 401"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + FINITE-EXACT + VERIFIED-EXACT /
  PENDING INDEPENDENT HOSTILE AUDIT.  Among 8,195 primitive parity-correct
  slopes through denominator 401, exact local screening through prime 997
  leaves 104 rows.  Complete generalized-Pell class and compiler-address
  orbit decisions leave exactly 42 slopes, each with an explicit invariant
  positive cone and hence an infinite positive two-cube ray.  Twenty-seven
  slopes are new beyond THM-3594.  This is bounded; no parametrization,
  density, or asymptotic is claimed.
source: kps-s189 / THM-3594 bounded slope-atlas continuation, 2026-08-21
audit: >
  PENDING -- author-side normal and optimized replays hash-pin the audited
  parent, all candidates, prime obstructions, 104 complete LMM/unit-orbit
  ledgers, 42 positive seeds, and the semantic transcript.  Independent
  clean-room reconstruction remains the promotion gate.
depends_on:
  - THM-3594-berggren-positive-cube-slope-atlas-through-201
related:
  - THM-3370-berggren-two-cube-biquadratic-norm-collision
  - THM-3375-berggren-positive-two-cube-pell-ray
  - THM-3580-berggren-positive-cube-slope-atlas-completion-through-101
script: 04-computation/berggren_positive_cube_slope_atlas_401_thm3640.py
output: 05-knowledge/results/berggren_positive_cube_slope_atlas_401_thm3640.out
script_sha256: 596ba45542c18b971fa73cdeaad8d48410b04be45961cf092bfea70468318d5d
output_sha256: 3ec90d9f1fdef199d780273e9b75c381967aee62ae520e441b5abb35e8b12cfa
semantic_sha256: 879a912b16752dc8b97393affeff61387eedc45d23c490a7c9ce5f54c1836374
hash_basis: raw LF bytes for files; canonical JSON for semantic ledger
---

# THM-3640 -- positive two-cube slope atlas through 401

**RESERVED / PROVISIONAL PROOF CANDIDATE + FINITE-EXACT + VERIFIED-EXACT /
PENDING INDEPENDENT HOSTILE AUDIT.**  This doubles THM-3594's denominator
range while retaining both halves of its proof: local primes may exclude a
slope, but only complete generalized-Pell classes and their finite compiler
orbits may certify the surviving slopes.

## 1. Universe and exact compiler

For

```text
3<=n<=401,       n odd,       m even,       n/2<m<n,
gcd(m,n)=1,                                               (1)
```

put

```text
d=n^2u^2+2,
e=uW,
a=mn^2u^3+(2m+n)u,
q=m^2n^2u^4+2m(m+n)u^2+1.                              (2)
```

There are exactly `8,195` slopes in `(1)`.  The load-bearing conic is

```text
3W^2=n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2).                (3)
```

Every local survivor below has `3` dividing neither `m` nor `n`.  Define

```text
K=(4m^2-n^2)/3,
C=4(2m^2+2mn-n^2)/3,
h=nu.                                                   (4)
```

Then `(3)` becomes

```text
W^2-Kh^2=C,                                             (5)
```

and the original odd compiler address is exactly

```text
W odd,                     h=n mod 2n.                 (6)
```

Whenever `(5)--(6)` hold with `0<W<nh`, formulas `(2)` give positive
distinct integers

```text
x=(d+e)/2,                 y=(d-e)/2,
r=(a-1)/2,
x^3+y^3=a^2+2=(2r+1)^2+2.                              (7)
```

The companion rechecks both polynomial identities behind `(7)` on the
selected seeds and on one full compiler-period iterate.

## 2. Local exclusion is only a first gate

For every prime `p<=997`, the verifier exhausts `(3)` modulo `p`.  The `168`
primes exclude `8,091` slopes and leave exactly `104`.  The exact ledgers have
digests

```text
obstruction histogram:
c430037d87fe8b94adb4a1b73ec1e3adf813215d50a7dcd41579581c9acbc990,

ordered survivor list:
2b1e24ba10e59deac190e1b4d20a0afe3792c138ecc84e06fa22e1671448ff8e. (8)
```

No existence conclusion is drawn from passing this screen.

## 3. Complete class and address-orbit decision

For every row surviving `(8)`, the complete Lagrange--Mollin--Matthews
algorithm enumerates all integer-solution classes of `(5)`.  If

```text
epsilon=P+Q sqrt(K),                 P^2-KQ^2=1,        (9)
```

is the fundamental positive unit, its action on `(W,h) mod 2n` is finite.
Every listed class is therefore followed through its full orbit and tested
against `(6)`, including signs and conjugation.

The `104` survivor rows have LMM-class-count histogram

```text
0:38, 1:4, 2:10, 3:9, 4:2, 5:1, 6:27, 9:3, 12:10.    (10)
```

Their numbers of parity-admissible classes have histogram

```text
0:62,                 2:22, 4:13, 6:2, 8:5.           (11)
```

Thus precisely `42` slopes survive globally:

```text
(14,23), (26,29), (26,47), (98,101),
(86,125), (86,149), (110,149), (122,149),
(86,167), (98,167), (110,167), (146,167),
(170,191), (110,197), (182,197),
(122,215), (182,215), (146,239), (182,239),
(158,263), (218,263), (230,263),
(146,269), (230,269), (242,269),
(170,293), (182,293), (194,293), (254,293), (278,293),
(158,311), (278,317), (218,335), (254,335),
(182,359), (194,359), (206,359), (278,365),
(194,383), (338,383), (350,383), (206,389).            (12)
```

The first fifteen are exactly THM-3594's atlas.  The remaining twenty-seven
are new.  The full representatives, units, finite orbit periods and good
phases have digest

```text
5390c4b137d3fdd90d90212a7b3c26cac2b819ff1fd1b7e0950149cf6e593265. (13)
```

## 4. A useful slope-coordinate normalization

Put

```text
u_0=n-m,                     v_0=2m-n.                 (14)
```

Then

```text
m=u_0+v_0,                  n=2u_0+v_0,

K=v_0(4u_0+3v_0)/3,
C=4(2u_0^2+6u_0v_0+3v_0^2)/3.                         (15)
```

Every one of the `104` local survivors satisfies

```text
3 divides u_0,                    3 does not divide v_0. (16)
```

This also follows structurally from integrality of both quantities in `(15)`:
with `3` dividing neither `m` nor `n`, the second numerator forces `mn=1`
modulo `3`, hence `m=n mod 3`.  The transformed list for `(12)` is pinned in
the transcript.  It makes the next question precise—Vieta moves in the
`(u_0,v_0)` slope plane—but `(16)` is not itself a parametrization.

The even good-class counts in `(11)` are likewise a bounded exact fact.  They
are compatible with sign/conjugation pairing, but no all-denominator pairing
theorem is asserted here.

## 5. Every admissible row gives an infinite positive ray

Choose any row in `(12)`, a good class, and a phase satisfying `(6)`.
Independent sign and conjugation make `W,h>0`.  Since

```text
sqrt(K)<n                                               (17)
```

is equivalent to `m<n`, positive unit iteration eventually enters `W<nh`.
For a complete unit-orbit period define

```text
L=nh-W,                         H=nW-Kh.               (18)
```

At the selected seed `L,H>0`; under the period unit
`P_T+Q_T sqrt(K)` they obey

```text
L'=P_T L+Q_T H,
H'=P_T H+KQ_T L.                                      (19)
```

Thus the positive cone is invariant, `(6)` recurs, and `h` grows strictly.
Equation `(7)` compiles infinitely many distinct positive two-cube points.
One exact seed and its next period iterate are checked for every slope in
`(12)`; the seed ledger digest is

```text
30af440ae6817b71045a3351d8480aef55c5ce62ac46580659738266f74fd6e5. (20)
```

## 6. Verification and stopping boundary

The assertion-free companion pins the promoted THM-3594 theorem, script and
output, reconstructs every candidate and local obstruction, enumerates every
class and complete address orbit, and verifies the cone compiler.  Normal and
optimized runs execute the same gates and must agree with the stored output.

Reproduce with

```bash
python3 04-computation/berggren_positive_cube_slope_atlas_401_thm3640.py
python3 -O 04-computation/berggren_positive_cube_slope_atlas_401_thm3640.py
```

The result is bounded by `(1)`.  It proves no unbounded slope family,
parametrization, natural density, asymptotic count, or complete census of
integers representable as sums of two distinct cubes.  Independent hostile
reconstruction remains required for promotion.
