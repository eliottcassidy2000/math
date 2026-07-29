---
id: THM-2863
title: "Endpoint Prony splitter and carry character-three intertwiner"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The
  frequency-cleared numerators of four consecutive 91-unit endpoint
  multiplier probes form a nonsingular two-node Prony system.  Its fixed
  and relative-character-three summands split exactly, and the latter has
  a canonical normalized equivariant match with the character-three
  projection of the THM-2851 carry derivative.  Simultaneous full-current
  survival, signed physical projection, and E3 transport remain open.
source: root/lrc-endpoint-prony-character-three-2026-07-28
depends_on:
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2857-endpoint-galois-carry-torsor-and-phase-alignment-sidecar
related:
  - THM-2286-endpoint-prony-lift-bank-and-sharp-owner-multiplier-landing
  - THM-2408-endpoint-prony-resultant-clock-separation-and-shared-node-boundary
  - THM-2502-endpoint-boolean-newton-carry-tournament-and-dipole-boundary
  - THM-2857-endpoint-galois-carry-torsor-and-phase-alignment-sidecar
  - THM-2861-endpoint-hermitian-edge-holonomy-and-semilinear-clutch-test
script: 04-computation/lrc14_endpoint_prony_carry_character3_thm2863.py
output: 05-knowledge/results/lrc14_endpoint_prony_carry_character3_thm2863.out
script_sha256: 3d09e2f06cd17c38a34b00fa133f67963287173e5a69092e4e2d2b190f4459b9
output_sha256: 2cf1902b19fd6ae9c6780521da0e431c8c7e1887235ea94ef4a8a5ca7cf9c7e4
hash_basis: LF-normalized bytes
---

# THM-2863 -- endpoint Prony splitter and carry character-three intertwiner

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The apparent two-endpoint scalar in THM-2847 can be split without choosing
one of the thirteen carry sections.  Four consecutive nonzero
multiplier probes recover its two Prony nodes and their signed weights.
Relative Galois descent then separates a fixed line from one faithful
`C_13` character.  That charged line is the same representation carried
by the character-three projection of THM-2851's oriented carry derivative.

This is an exact algebraic interface and a physical endpoint-factor replay,
not an LRC(14) closure.  Canon does not yet provide the four coefficients
as simultaneously surviving coefficients of one fully typed current, a
lawful signed endpoint projection, or the q11-to-q7 E3 transporter.

## 1. The four-sample physical endpoint bank

Put

```text
xi=zeta_2366,             F=Q(zeta_91)=Q(xi^13).
```

On the literal q3 interval of THM-2847, probe the right endpoint at

```text
Y(m)=12+26m,              m=1,2,3,4.                     (1)
```

Thus `Y=(38,64,90,116)`, and every entry is a unit modulo `91`.  Let
`C_m` be the builder's frequency-cleared endpoint numerator.  The actual
interval Fourier coefficient differs from `C_m` by an explicit known
nonzero scalar proportional to `Y(m)^(-1)` (with the fixed Fourier-sign
normalization).  Clearing this scalar preserves vanishing and is necessary
before applying the constant-coefficient recurrence below.

Exact replay of the endpoint builder at the three live rows

```text
(q,origin)=(3,00),(11,00),(11,12)
```

gives the same four nonzero pairs in both certified finite fields.  The
four reduced cyclotomic exponent pairs are

```text
(2203,65), (2216,234), (2229,403), (2242,572).           (2)
```

The `(3,12)` row stays empty.  Formula `(2)` is exactly

```text
C_m=alpha_L lambda_L^m-alpha_R lambda_R^m,               (3)

alpha_L=xi^2190,       lambda_L=xi^13,
alpha_R=xi^2262,       lambda_R=xi^169.
```

At `m=1`,

```text
C_1=xi^2203-xi^65
   =zeta_1183^624-zeta_1183^510,                         (4)
```

the canonical THM-2847/2857 endpoint scalar.  The replay in this section
certifies these frequency-cleared endpoint numerators only.  It does not
assert that the other factors of one marked current survive simultaneously
at all four values of `Y(m)`.

## 2. Exact two-node Prony splitting

Both nodes lie in `F`, are distinct, and have respective orders `182` and
`14`.  No integer-index numerator vanishes: `C_m=0` would force

```text
72+156m=0 mod2366,
```

which is impossible because `gcd(156,2366)=26` does not divide `72`.
Hence

```text
C_(m+2)
 =(lambda_L+lambda_R)C_(m+1)
 -lambda_L lambda_R C_m.                                (5)
```

More generally every consecutive Prony window is nonsingular:

```text
Delta_j=C_j C_(j+2)-C_(j+1)^2
 =-alpha_L alpha_R (lambda_L lambda_R)^j
   (lambda_L-lambda_R)^2 !=0.                           (6a)
```

In particular, the determinant for recovering the two recurrence
coefficients from `C_1,C_2,C_3,C_4` is

```text
Delta=C_1 C_3-C_2^2
 =-alpha_L alpha_R lambda_L lambda_R
   (lambda_L-lambda_R)^2 !=0.                           (6)
```

Thus the four samples recover the unordered nodes as the roots of

```text
T^2-(lambda_L+lambda_R)T+lambda_L lambda_R.
```

Their different orders canonically orient them.  Once the nodes are known,
two adjacent samples split the signed endpoint terms:

```text
L_j=alpha_L lambda_L^j
   =(C_(j+1)-lambda_R C_j)/(lambda_L-lambda_R),

-R_j=-alpha_R lambda_R^j
   =(lambda_L C_j-C_(j+1))/(lambda_L-lambda_R).          (7)
```

The signs and weight-to-node assignment in `(3)` are load-bearing.
Swapping the two endpoint weights preserves the recurrence polynomial but
changes which line is charged.  This is the precise node-only boundary
behind THM-2502: Prony works here because the oriented samples, not merely
the node set, are retained.

## 3. The relative character-three line

The generator of `Gal(Q(xi)/F)` is

```text
sigma(xi)=xi^1275,       1275=1+14*91.                  (8)
```

It fixes both nodes and the right weight, while

```text
sigma(alpha_R)=alpha_R,
sigma(alpha_L)=omega^3 alpha_L,     omega=xi^182.        (9)
```

Consequently, for every section `r in F_13`,

```text
sigma^r(C_m)=omega^(3r)L_m-R_m.                         (10)
```

With Fourier convention

```text
fhat(k)=sum_(r=0)^12 f(r)omega^(-kr),
```

equation `(10)` has support exactly `{0,3}`:

```text
widehat{sigma^r(C_m)}(0)=-13R_m,
widehat{sigma^r(C_m)}(3)= 13L_m.                        (11)
```

Thus the endpoint orbit function lies in

```text
F*(trivial character) direct-sum F*(character 3),        (12)
```

and its reduced summand is one-dimensional.  At `m=1`,
`-R_1=zeta_1183^624` is the Galois mean and
`L_1=-zeta_1183^510` is the centered eigenline.

## 4. The normalized carry intertwiner

THM-2851 gives the oriented first-carry derivative

```text
d h=449(delta_1-delta_0).                               (13)
```

Its character-three Fourier coefficient is

```text
(d h)hat(3)=449(omega^(-3)-1) !=0.                      (14)
```

Let `Pi_3` denote the character-three projector.  There is a unique
`C_13`-equivariant isomorphism normalized to carry the displayed generator
`Pi_3(d h)` to the centered endpoint orbit.  It is multiplication by

```text
iota_m=
  13L_m/[449(omega^(-3)-1)].                            (15)
```

Both spaces in `(15)` are one-dimensional; the word `unique` refers to
this normalization, not to an unbased isomorphism between character
lines.  This is the smallest possible nonzero linear repair: the reduced
representation has dimension one.  It is categorically different from a
deterministic carry section, whose sharp transitive-set cost remains
thirteen states by THM-2851.

## 5. Exact verification

The companion pins THM-2847, THM-2851, and THM-2857.  It:

1. rebuilds the literal endpoint interval and all four multiplier samples
   on the four q/origin rows;
2. verifies `(2)--(7)`, including all-index numerator nonvanishing,
   every-window nonsingularity, and the two node orders, in exact integer
   and finite-field arithmetic;
3. checks the relative Galois action, Fourier support `{0,3}`, and both
   coefficients in `(11)`;
4. constructs `(15)` on the actual projected carry generator; and
5. checks a swapped-weight node-only hostile.

Normal, optimized, and stored-output replay must be byte-identical after
LF normalization.  LF-normalized SHA-256:

```text
script  3d09e2f06cd17c38a34b00fa133f67963287173e5a69092e4e2d2b190f4459b9
output  2cf1902b19fd6ae9c6780521da0e431c8c7e1887235ea94ef4a8a5ca7cf9c7e4
```

An independent immutable hostile audit replayed the final blobs and
separately rederived the cleared exponent progression, all-index
nonvanishing, every shifted Hankel minor, both node orders and split signs,
the relative Galois decomposition, Fourier normalization, and the
normalized intertwiner `(15)`.  It accepted the endpoint-numerator-only
scope and found no defect.

## 6. Connection contract

```text
source:
  THM-2847 oriented endpoint scalar and
  THM-2851 first-carry derivative;

target:
  a one-dimensional endpoint/carry character-three interface;

map:
  four-sample Prony recovery, signed two-sample splitting, relative
  Galois centering, then the normalized intertwiner (15);

preserved:
  endpoint orientation, multiplier progression, relative Galois
  character, and carry orientation;

destroyed / still missing:
  positivity of the separated boundary term, simultaneous survival of
  the rest of the fully typed current, lawful access to the signed carry
  derivative, E3 transport, owner phase, and a deterministic section;

needed sidecar:
  one common-gauge full-current multiplier bank plus a lawful signed
  projection and q11-to-q7 E3 co-support;

cheapest decisive test:
  exhibit two (nodes known) or four (nodes unknown) consecutive 91-unit
  multiplier coefficients of one fully typed q3/q11 current, apply the
  known frequency clearing, and verify that (7) survives every retained
  factor.
```

No row is excluded and the LRC(14) ledger remains `165`.
