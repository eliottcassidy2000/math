---
id: THM-3711
title: "Complete all-scale W008 arm-square nonentry in the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT.  The full W008 support ray
  (n,n;2n,n,n) is Darboux-empty in the y=0 collision ring for every positive
  integral scale.  THM-3613 leaves only one placement at n=2 and two at every
  n>=3.  In the first family one collision row forces the negative arm
  coefficient to have an H^2 factor; in the second, two consecutive rows do.
  The n=3 extra arm addresses have the same square factor.  This closes the
  whole named word W008, not the full 3x4 cell or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion reconstructs the W008
  fibre word, deduplicates actual supports, reproduces the inherited parity
  rejection, and verifies both all-parameter differential integrations with
  their exact constants and signs.  An independent derivation checked the
  post-parity census, UFD bases, both normalized ODEs, every live-scale and
  exceptional arm address, coefficient-module divisibility, transcript, and
  hashes.  Normal and optimized runs byte-match the frozen transcript.
depends_on:
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3708-y0-w002-complete-all-scale-nonentry
script: 04-computation/jacobian_y0_w008_complete_arm_square_thm3711.py
output: 05-knowledge/results/jacobian_y0_w008_complete_arm_square_thm3711.out
script_sha256: ca000de4d2ae53abb4289aa3742c83ad54fe1d5bd67657fdc40f74fc76824772
output_sha256: 4011081487ae029e1886b3f534ab763ce059a30419b5d2e9bdd678b07c3e53c6
hash_basis: LF-normalized bytes
---

# THM-3711 -- the entire W008 ray is empty

**PROVED + VERIFIED-EXACT.**  Reversing the long gap in W002 produces W008,
but reversal is not a regularity symmetry of the collision ring.  The proof
below is therefore intrinsic: it uses the two low collision rows to propagate
the negative singleton base into the scalar arm.

Work over `C` in the THM-3696 collision ring.  Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

Thus the bracket of homogeneous terms of weights `r,s` has coefficient
`W_(r,s)(F,G)`.

## 1. Exact post-parity census

The W008 ray and its ordered fibres are

```text
(X,Y;U,V,W)=(n,n;2n,n,n),                     n>=1,

00; 10; 01+20; 02+11; 03+12+21; 13+22; 23.            (2)
```

Anchor a nonsingleton fibre at `(-2,1)` or `(1,-2)`, impose the singleton
same-sign/zero gate, deduplicate actual support pairs, and apply the exact
THM-3613 common-base parity obstruction.  The inherited candidate and
survivor counts are

```text
scale                  n=1       n=2       n=3       n>=4
inherited candidates     0         2         4          6
post-parity survivors    0         1         2          2.            (3)
```

The surviving actual supports are exactly

| family | scalar fibre | `wt(P)` | `wt(Q)` | live scales |
|---|---|---|---|---|
| X | `03+12+21` | `(1-2n,1-n,1)` | `(-2n-2,-2,n-2,2n-2)` | `n>=2` |
| Y | `13+22` | `(1-2n,1-n,1)` | `(-3n-2,-n-2,-2,n-2)` | `n>=3` |

The companion reconstructs `(2)` and `(3)` directly rather than counting
different arm anchors on the same support as different placements.

## 2. The common negative singleton base

Write the coefficient polynomials of `P,Q` as `f_i,g_j`.  In both surviving
families, singleton rows `00` and `10` connect `f_0,f_1,g_0` in one negative
Wronskian-zero component.  In family X the absolute weights in this component
are

```text
2n-1, n-1, 2n+2,
```

and in family Y they are

```text
2n-1, n-1, 3n+2.
```

In both cases the gcd is one because `gcd(2n-1,n-1)=1`.  UFD
common-power rigidity therefore gives, for nonzero scalars,

```text
f_0=aH^(2n-1),                  f_1=c_0H^(n-1),        (4)

family X: g_0=b_0H^(2n+2),
family Y: g_0=b_0H^(3n+2).
```

Write the active weight-one coefficient as

```text
f_2=d_0K.                                               (5)
```

The negative coefficient-module law in THM-3696 gives `h|f_0`; since `h` is
squarefree, `(4)` implies

```text
h|H.                                                    (6)
```

## 3. Family X: one row creates the arm square

The zero row `01+20` is

```text
W_(1-2n,-2)(f_0,g_1)+W_(1,-2n-2)(f_2,g_0)=0.          (7)
```

Substitute `(4),(5)` and work in `C(b)`.  After division by the common
nonzero factor, `(7)` is precisely

```text
(g_1/H^2)'=gamma(HK)',
gamma=2(n+1)b_0d_0/[a(2n-1)].                         (8)
```

Hence, with an arbitrary scalar `kappa`,

```text
g_1=H^2(kappa+gamma HK).                               (9)
```

The scalar fibre is `03+12+21`.  Address `21` has weights `(1,-2)` at every
scale, and its negative coefficient is the square-bearing `g_1` in `(9)`.
At `n=3` there is one additional arm address, `12=(-2,1)`; its negative
coefficient is `f_1=c_0H^2`.  There are no other arm addresses.  By `(6)`,
each negative arm coefficient and its derivative vanish at `b=+-1`.
THM-3696 says every non-arm scalar address also vanishes there.  The scalar
row therefore cannot equal one.  Family X is empty for every `n>=2`.

## 4. Family Y: two rows transport the square

The first low row `01+20` now reads

```text
W_(1-2n,-n-2)(f_0,g_1)+W_(1,-3n-2)(f_2,g_0)=0.       (10)
```

It integrates completely to

```text
g_1=H^(n+2)(kappa+gamma HK),
gamma=(3n+2)b_0d_0/[a(2n-1)].                         (11)
```

Substitute `(11)` in the next zero row `02+11`:

```text
W_(1-2n,-2)(f_0,g_2)+W_(1-n,-n-2)(f_1,g_1)=0.        (12)
```

After the same exact rational-derivative normalization,

```text
(g_2/H^2)'=-eta(HK)',
eta=c_0(n-1)gamma/[a(2n-1)],                          (13)
```

and therefore

```text
g_2=H^2(lambda-eta HK)                                (14)
```

for an arbitrary scalar `lambda`.

The scalar fibre is `13+22`.  Address `22=(1,-2)` is always an arm address,
and `(14)` makes its negative coefficient square-bearing.  At the exceptional
scale `n=3`, `13=(-2,1)` is the only additional arm address and its negative
coefficient is again `f_1=c_0H^2`.  The same evaluation at either arm gives
zero instead of one.  Family Y is empty for every live scale `n>=3`.

## 5. Exhaustion and scope

THM-3613 removes every inherited placement outside X and Y.  Sections 3 and
4 remove X and Y uniformly, including the only exceptional multiple-arm
scale.  Consequently

```text
no retained W008 3 x 4 Darboux pair exists in R at any n>=1.             (15)
```

Output exchange gives the transposed conclusion.  This is a complete
nonentry theorem for the named word W008.  It does not close W003--W007,
arbitrary `3 x 4` supports, general quartic C3 data, or `JC(2)`.

Run

```bash
python3 -B 04-computation/jacobian_y0_w008_complete_arm_square_thm3711.py
python3 -B -O 04-computation/jacobian_y0_w008_complete_arm_square_thm3711.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
