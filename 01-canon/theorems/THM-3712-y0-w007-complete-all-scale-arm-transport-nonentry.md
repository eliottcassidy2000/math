---
id: THM-3712
title: "Complete all-scale W007 arm-transport nonentry in the y=0 collision ring"
status: >
  PROVED + VERIFIED-EXACT.  The full W007 support ray
  (2n,n;n,n,n) is Darboux-empty in the y=0 collision ring at every positive
  integral scale.  THM-3613 leaves no placement through n=2, three actual
  supports at n=3, and four thereafter.  One or two zero collision rows force
  every surviving negative arm coefficient to carry H^2.  The C/D supports
  merge at n=3, where the same calculation kills both eligible arm addresses.
  This closes the named word W007, not the full 3x4 cell or JC(2).
source: jc-quartic-c3-construction / 2026-08-22
audit: >
  INDEPENDENTLY HOSTILE-AUDITED.  The exact companion reconstructs the W007
  fibre word, deduplicates actual supports, reproduces the inherited parity
  rejection and n=3 merger, and checks all four all-parameter differential
  integrations.  An independent derivation checked the C/D triple signs,
  quotient polynomiality, common-base gcds, exact C/D identification at n=3,
  and both arm coefficients on the merged support.  Normal and optimized runs
  byte-match the frozen transcript.
depends_on:
  - THM-3603-three-by-four-support-collision-cone-and-fibre-cut-atlas
  - THM-3606-exponent-two-three-by-four-scalar-singleton-gate-atlas
  - THM-3613-three-by-four-size-seven-ray-parity-gate
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3708-y0-w002-complete-all-scale-nonentry
  - THM-3711-y0-w008-complete-all-scale-arm-square-nonentry
script: 04-computation/jacobian_y0_w007_complete_arm_transport_thm3712.py
output: 05-knowledge/results/jacobian_y0_w007_complete_arm_transport_thm3712.out
script_sha256: 52fabfdf0fbf862724c3ba8c26efdb26654c16ec0df2aa087578554bef1229d6
output_sha256: 67deba4374ab48351e82762941b3191fa922c1fd6ecd73883a6279ea7ba09786
hash_basis: LF-normalized bytes
---

# THM-3712 -- the entire W007 ray is empty

**PROVED + VERIFIED-EXACT.**  Work over `C` in the THM-3696 collision ring.
Put

```text
h=1-b^2,
W_(r,s)(F,G)=sF'G-rFG'.                                (1)
```

The mechanism is an arm-transport law: zero collision rows carry the common
negative singleton base into whichever weight `-2` coefficient was not
visible to the singleton graph.

## 1. Exact post-parity census

The W007 ray and its ordered fibres are

```text
(X,Y;U,V,W)=(2n,n;n,n,n),                     n>=1,

00; 01; 02+10; 03+11+20; 12+21; 13+22; 23.            (2)
```

After the inherited scalar/singleton gates, actual-support deduplication, and
the THM-3613 parity obstruction, the exact counts are

```text
scale                  n=1       n=2       n=3       n>=4
inherited candidates     2         4         6          8
post-parity survivors    0         0         3          4.            (3)
```

The four tail families are

| family | scalar fibre | `wt(P)` | `wt(Q)` |
|---|---|---|---|
| A | `12+21` | `(-2n-2,-2,n-2)` | `(1-2n,1-n,1,n+1)` |
| B | `12+21` | `(1-2n,1,n+1)` | `(-2n-2,-n-2,-2,n-2)` |
| C | `13+22` | `(1-3n,1-n,1)` | `(-2n-2,-n-2,-2,n-2)` |
| D | `13+22` | `(-2n-2,-2,n-2)` | `(1-3n,1-2n,1-n,1)` |

All four formulas are live for `n>=3`, but at `n=3` C and D are the same
actual support pair

```text
P=(-8,-2,1),                    Q=(-8,-5,-2,1).        (4)
```

It retains both arm anchors.  This is one coefficient system, not two.

## 2. Common negative bases and the arm divisor

Singletons `00,01` connect `f_0,g_0,g_1` in one negative Wronskian-zero
component.  The absolute weight triples in A--D are respectively

```text
A: (2n+2,2n-1,n-1),       B: (2n-1,2n+2,n+2),
C: (3n-1,2n+2,n+2),       D: (2n+2,3n-1,2n-1).       (5)
```

Every triple has gcd one.  For A or B, a possible factor three in the first
two large entries does not divide `n-1` or `n+2`, respectively.  For C, the
first entry removes the only possible factor two in the last two; D already has
`gcd(3n-1,2n-1)=1`.  UFD common-power rigidity therefore supplies one base
`H` in each family with exponents equal to `(5)`.  Explicitly,

```text
A: f_0=aH^(2n+2), g_0=b_0H^(2n-1), g_1=c_0H^(n-1),
B: f_0=aH^(2n-1), g_0=b_0H^(2n+2), g_1=c_0H^(n+2),
C: f_0=aH^(3n-1), g_0=b_0H^(2n+2), g_1=c_0H^(n+2),
D: f_0=aH^(2n+2), g_0=b_0H^(3n-1), g_1=c_0H^(2n-1). (6)
```

All displayed scalars are nonzero.  Negative coefficient membership and the
squarefreeness of `h` imply

```text
h|H                                                     (7)
```

in every family.

## 3. Families A and B: direct one-row transport

In A write the weight-one coefficient `g_2=d_0K`.  The zero row `02+10` is

```text
W_(-2n-2,1)(f_0,g_2)+W_(-2,1-2n)(f_1,g_0)=0.          (8)
```

After substituting `(6)`, it is exactly

```text
(f_1/H^2)'=gamma_A(HK)',
gamma_A=2(n+1)ad_0/[(2n-1)b_0].                       (9)
```

Thus

```text
f_1=H^2(kappa+gamma_A HK).                            (10)
```

The persistent scalar arm `12=(-2,1)` has negative coefficient `f_1`.  At
`n=3`, the extra arm `21=(1,-2)` has negative coefficient
`g_1=c_0H^2`; so all A arms die to second order on `h=0`.

In B write the weight-one coefficient `f_1=d_0K`.  The same row, with the
roles reversed, integrates to

```text
g_2=H^2(kappa+gamma_B HK),
gamma_B=2(n+1)b_0d_0/[a(2n-1)].                       (11)
```

Here the persistent arm is `12=(1,-2)`, with negative coefficient `g_2`.
There is no additional arm address.  Thus A and B are empty.

## 4. Family C: triple polynomiality followed by arm transport

The positive singleton `23` gives

```text
f_2=d_0K,                         g_3=t_0K^(n-2).      (12)
```

Put `Z=f_1/H^(n-1)` in `C(b)`.  The zero triple `03+11+20` integrates to

```text
Z=kappa+alpha(HK)^(n-2)-beta HK,
alpha=a t_0(3n-1)/[c_0(n+2)],
beta=2(n+1)b_0d_0/[c_0(n+2)].                         (13)
```

In particular, `(13)` proves rather than assumes that `Z` is polynomial.
The lower double `02+10` then becomes

```text
(g_2/H^2)'=theta Z',
theta=2(n+1)b_0/[a(3n-1)],                            (14)
```

so

```text
g_2=H^2(lambda+theta Z).                              (15)
```

The persistent arm `22=(1,-2)` now dies to second order.  At `n=3`, the
additional arm `13=(-2,1)` has

```text
f_1=H^(n-1)Z=H^2Z,                                   (16)
```

so it dies as well.

## 5. Family D and the n=3 merger

The positive singleton now gives

```text
f_2=d_0K^(n-2),                    g_3=t_0K.           (17)
```

The zero triple `03+11+20` directly integrates the arm coefficient:

```text
f_1=H^2[
 kappa+alpha_D HK-beta_D(HK)^(n-2)
],

alpha_D=2(n+1)a t_0/[(2n-1)c_0],
beta_D=(3n-1)b_0d_0/[(2n-1)c_0].                     (18)
```

For `n>=4`, `13=(-2,1)` is the only arm address, so `(18)` closes D.  At
`n=3`, address `22=(1,-2)` is also eligible.  But `(4)` says that this is
exactly the already treated C support, and `(15)` gives its second negative
arm coefficient the factor `H^2`.  Thus the merger is closed without
discarding either arm orientation.

## 6. Arm evaluation, exhaustion, and scope

By `(7)`, every negative arm coefficient displayed in `(10),(11),(15),
(16),(18)` and its derivative vanish at `b=+-1`.  THM-3696 says that every
other complementary-weight address vanishes there.  Hence no scalar row can
equal one.  Together with the exact census, this proves

```text
no retained W007 3 x 4 Darboux pair exists in R at any n>=1.             (19)
```

Output exchange gives the transposed conclusion.  This closes the named word
W007, not W003--W006, arbitrary `3 x 4` supports, general quartic C3 data, or
`JC(2)`.

Run

```bash
python3 -B 04-computation/jacobian_y0_w007_complete_arm_transport_thm3712.py
python3 -B -O 04-computation/jacobian_y0_w007_complete_arm_transport_thm3712.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
