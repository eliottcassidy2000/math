---
id: THM-2994
title: "First-gap norm-core Hurwitz and Hermite--Biehler prefix through width sixteen"
status: PROVED + VERIFIED-EXACT + CERTIFIED-ARB + INDEPENDENTLY HOSTILE-AUDITED
source: codex-gmc-first-gap-hurwitz-prefix-2026-07-30
audit: >
  Independent immutable replay and proof audit accepted the canonical
  THM-2994 package after all three required repairs. Normal, optimized, and
  stored transcripts are byte-identical; reciprocal Routh scaling, strict
  Hurwitz inference, Hermite--Biehler typing, Arb certification, hashes,
  record digest, and the repaired THM-2991 global-return control all pass.
depends_on:
  - THM-2969-first-gap-wall-stripped-resultant-norm-core-atlas
related:
  - THM-2982-first-gap-wall-stripped-norm-core-strict-ulc-through-thirty-four
  - THM-2989-first-gap-wall-stripped-all-width-leading-edge-positivity
  - THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return
script: 04-computation/gmc_first_gap_hurwitz_hermite_biehler_prefix_thm2994.py
output: 05-knowledge/results/gmc_first_gap_hurwitz_hermite_biehler_prefix_thm2994.out
script_sha256: ff6c934bf179aa6fbf59d3e2fddfac8fa70e97b2741c1bf5ac2e9d705417a8fa
output_sha256: a2c354eafd4f3e6cbaee84b2936b5f3b2f152ea31eb81237e8cfa5618d94aebb
hash_basis: LF-normalized bytes
---

# THM-2994 -- first-gap Hurwitz/Hermite--Biehler prefix

**PROVED + VERIFIED-EXACT + CERTIFIED-ARB + INDEPENDENTLY
HOSTILE-AUDITED.**

## 1. Statement

For every integer width

```text
6 <= M <= 16,                                             (1)
```

let `C_M(n)` be the genuine first-gap wall-stripped norm core of PROVED
THM-2969 for support `(0,1,2,M)`.  Every zero `z` of `C_M` satisfies

```text
Re z < 0.                                                 (2)
```

Thus all eleven cores in `(1)` are Hurwitz stable.

Write uniquely

```text
C_M(x)=E_M(x^2)+x O_M(x^2).                              (3)
```

The classical Hermite--Biehler criterion then says that the zeros of
`E_M(-t)` and `O_M(-t)` are real, positive, and strictly alternate.  The
orientation begins with an `E` zero.  It ends with `O` when `deg C_M` is odd
and with `E` when `deg C_M` is even.

The companion independently isolates the Hermite--Biehler legs for
`M=6,7,8`.  For `M=9,...,16`, interlacing in this candidate is the theorem
corollary of the exact Hurwitz result, not a separately executed root census.

This is a finite root-location theorem.  It proves neither PF-infinity nor a
Newton-circuit no-return property.

## 2. Exact reciprocal Routh certificate

The native coefficient list of `C_M` is constant-first.  If `d=deg C_M`,
read that list as the descending coefficient list of

```text
C_M^vee(s)=s^d C_M(1/s).                                (4)
```

THM-2969 gives `C_M(0)>0`, so neither core has a zero at zero.  Inversion
preserves the open left half-plane because

```text
Re(1/z)=Re(z)/|z|^2.                                    (5)
```

Consequently `C_M` is Hurwitz if and only if `C_M^vee` is Hurwitz.

Construct the Routh array of `(4)`.  The exact primitive first-column records
are:

| `M` | degree | positive pivots | max first-column pivot height (bits) | first-column SHA256 |
|---:|---:|---:|---:|:---|
| 6 | 121 | 122 | 38376 | `bf5e565214f4cf3257d1a6a1055b8bf6b58b9ae68bd55c1aa32236dbd814f3bb` |
| 7 | 144 | 145 | 55189 | `7e1ed405d784ce7c07664867ae9c3fc32547bd512057357881a18c41aa2e5f61` |
| 8 | 164 | 165 | 74026 | `60844d75371190769d5db5435ebe5da049c693afd29ea48b1ddf261d11b436ab` |
| 9 | 184 | 185 | 97503 | `c9bea1dc9b8ac06c50d46d928e34bc1caaeb26fd2fd24debdcb94ffcc2863cb5` |
| 10 | 205 | 206 | 119890 | `7151ba9b535619b9006d6f72abeeef6a8ff979c75a5476b3b04024d968dfcd9f` |
| 11 | 226 | 227 | 148547 | `ce7eddf2ecee0dda531b1bfc7cc8123a1402c62c08cffddaaecdf49238103fa3` |
| 12 | 244 | 245 | 174313 | `0082b706509af1d93a5afe5f0882d0395e19def94b80105a50a5da3cc7f53ae8` |
| 13 | 268 | 269 | 209343 | `eb765ac99562c13db36d8a47dfd4142f04bb196164b153b2d1faa4ad4d28bf68` |
| 14 | 288 | 289 | 242649 | `474e601046f515353c6a6aa2f10cdb67328c9c9c39cf61bbec9d53249c1d7339` |
| 15 | 308 | 309 | 283228 | `e70693e4413446420592c976921aeddd7d6c73529fa0380c427b3a95de226b68` |
| 16 | 329 | 330 | 321582 | `5d3065df61595d91a34815b597d702fba8ca2c1e8e13ee2448369fb3634d48fc` |

There are no zero or negative pivots.  The Routh--Hurwitz criterion proves
that every zero of `(4)` lies strictly in the open left half-plane.  Equations
`(4)--(5)` prove `(2)`.  Width `13` is a useful internal hostile: its pure
quadratic and cubic resonance walls occur at the distinct roots `9` and `10`,
yet the certificate remains strict.

## 3. Why the fraction-free recursion is lawful

The executed recursion never forms giant rational rows.  Suppose stored
consecutive rows `U,L` are respectively positive multiples `alpha,beta` of
the ordinary rational Routh rows.  Since the already certified lower pivot is
positive, the raw next row

```text
N_j=L_0 U_(j+1)-U_0 L_(j+1)                            (6)
```

is `alpha beta` times the positive ordinary lower pivot times the ordinary
next row.  Hence it is a positive multiple of that row.  Dividing `N` by its
positive integer content changes no sign.  The initial rows are positive
multiples of themselves.  Induction proves that every printed pivot has the
same sign as the corresponding ordinary Routh pivot.  This establishes the
exact input to Routh--Hurwitz without a floating sign decision.

## 4. Hermite--Biehler typing

On the imaginary axis, `(3)` gives

```text
C_M(i sqrt(t))=E_M(-t)+i sqrt(t) O_M(-t).               (7)
```

Because `(2)` is strict, `C_M` has no imaginary-axis zero.  The
Hermite--Biehler theorem applied to `(7)` gives real simple alternating zeros
for the two displayed legs.  Positive coefficients orient the first crossing
as `E`; the degrees orient the final crossing.

As an independent interval check, the companion isolates all leg zeros at
the first three widths:

```text
M=6: deg(E),deg(O)=60,60;  E<O<...<O,
M=7: deg(E),deg(O)=72,71;  E<O<...<E,
M=8: deg(E),deg(O)=82,81;  E<O<...<E.                  (8)
```

Every isolating ball is strictly positive, every imaginary part is exactly
zero, and every adjacent pair is interval-separated.  Equation `(8)` is
corroboration, not a hidden premise for `(2)`.

## 5. Rigorous nonreal-root corroboration

Arb isolation of the full core through width `12` gives exact certified
`(real,nonreal)` counts

```text
M=6,...,12:
(37,84), (44,100), (46,118), (54,130),
(59,146), (66,160), (64,180).                          (9)
```

All real-part balls are strictly negative and every alleged nonreal ball has
imaginary part separated from zero.  In particular the first core already
has `84` nonreal zeros.  The result is therefore genuinely Hurwitz and is not
a disguised real-negative-root/PF-infinity theorem.

## 6. Sharp no-return boundary

Hurwitz stability alone does not imply even ULC: `x^2+x+1` is Hurwitz, while
its binomially normalized middle Newton ratio is `1/4`.

Much more strongly, PROVED THM-2991 constructs PF-infinity, Hurwitz, strictly
ULC polynomials with an arbitrarily long improving leading Newton-ratio prefix
and a later global return below the leading ratio.  Its smallest frozen strong
control is

```text
P_(2,20)=(x+1)^2(x+3)(x+20)^2
        =x^5+45x^4+607x^3+2283x^2+2920x+1200,
(R_1,R_2,R_3,R_4)
 =(810/607,368449/205470,5212089/3544880,42632/34245),
R_1<R_2,  R_4<R_1.                                      (10)
```

Thus a fixed root class, even PF-infinity, cannot prove the missing
family-specific no-return theorem after THM-2989.  The old two-cluster
degree-four control only turns downward and does not cross `R_1`; it is not
used here.  The exact connection is

```text
source:     the wall-stripped norm core;
target:     a zero-free right half-plane / HB pair;
preserved:  all coefficients and the full reciprocal Routh chain;
destroyed:  none inside the finite root theorem;
missing for no-return:
            a family-specific wall/response continuation law controlling
            every unobserved Newton circuit.                           (11)
```

## 7. Scope

This theorem proves only the finite Hurwitz statement `(1)--(2)` and its
Hermite--Biehler corollary.  It
does not prove:

- Hurwitz stability for `M>=17`;
- PF-infinity or real-rootedness at any displayed width;
- strict ULC beyond any already proved finite atlas;
- the THM-2989 no-return obligation;
- arbitrary-radial NC2, a new Gaussian moment cutoff, or GMC(2).

The companion binds the LF hash of the immutable THM-2969 exact engine,
reconstructs every core from that engine, verifies all eleven primitive Routh
records, checks `(8)--(10)`, and reports a deterministic transcript.  Normal
and optimized replays are byte-identical (19 lines, 4666 bytes).  Reproduce
with

```text
python 04-computation/gmc_first_gap_hurwitz_hermite_biehler_prefix_thm2994.py --output thm2994.normal.out
python -O 04-computation/gmc_first_gap_hurwitz_hermite_biehler_prefix_thm2994.py --output thm2994.opt.out
```

The frozen LF-normalized hashes are

```text
script  ff6c934bf179aa6fbf59d3e2fddfac8fa70e97b2741c1bf5ac2e9d705417a8fa
output  a2c354eafd4f3e6cbaee84b2936b5f3b2f152ea31eb81237e8cfa5618d94aebb
```

The independent hostile audit rederived reciprocal orientation, fraction-free
Routh scaling, strict Hurwitz inference, Hermite--Biehler interlacing, Arb
certification, and scope.  It required the repaired three-cluster control,
precise first-column height label, and canonical THM-2994 namespace installed
here.  The final immutable normal and optimized replays both equal the stored
transcript byte-for-byte: `19` lines and `4666` LF bytes.

**QED.**
