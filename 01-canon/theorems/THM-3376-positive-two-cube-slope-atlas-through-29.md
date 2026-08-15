---
id: THM-3376
title: "Positive Berggren two-cube slope atlas through denominator 29 and a second Pell ray"
status: PROVED + FINITE-EXACT + VERIFIED-EXACT / INDEPENDENT AUDIT PENDING
source: kps-s181
depends_on:
  - THM-3375
related:
  - THM-3370
companion: 04-computation/berggren_positive_cube_slope_atlas_29_kps_s182.py
output: 05-knowledge/results/berggren_positive_cube_slope_atlas_29_kps_s182.out
script_sha256: 7ec10b47139f10117a9f80155ac040ac864bb329ea5687ac1c890bec27529d30
output_sha256: 0e5660b462b6f656b8cf09e26e5f6826e74fe626535cb8188fa3fd0f16fa4a44
hash_basis: LF-normalized bytes
---

# THM-3376 -- the first slope atlas has exactly two Pell rays

**PROVED + FINITE-EXACT + VERIFIED-EXACT / INDEPENDENT AUDIT PENDING.**

## 1. Primitive parity-correct slope universe

THM-3375's reduced-seed compiler starts from

```text
d=n^2u^2+2,
a=mn^2u^3+(2m+n)u,
e=uW,                                                    (1)
```

and reduces the cube condition to

```text
3W^2=n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2).                (2)
```

For an odd target `Q_r`, the positive cube pair has odd `d=x+y`.  In this
ansatz that forces `n,u` odd; then `a=2r+1` forces `m` even.  A common odd
factor of `m,n` can be absorbed into `u`, so take `gcd(m,n)=1`.  Positivity of
the limiting chamber asks `n/2<m<n`.

There are exactly `47` primitive parity-correct slopes with

```text
3<=n<=29,       n odd,       m even,       n/2<m<n.     (3)
```

Among these, equation `(2)` has an integer solution with odd `u,W` exactly
for

```text
                         (m,n)=(14,23),(26,29).         (4)
```

## 2. Finite modular obstruction atlas

For each slope, the companion enumerates all `(u,W) mod p` in `(2)`.  Modulo
`3`, `35` slopes have no point.  The twelve survivors of that first screen are

```text
(4,7),(8,11),(10,13),(14,17),(10,19),(16,19),
(14,23),(20,23),(16,25),(22,25),(20,29),(26,29).       (5)
```

The remaining exact obstructions are

| prime | slopes with no point modulo that prime |
|---:|---|
| `5` | `(4,7),(14,17)` |
| `7` | `(10,13),(20,23),(16,25)` |
| `11` | `(8,11),(20,29)` |
| `13` | `(10,19),(16,19)` |
| `23` | `(22,25)` |

This excludes `45` of the `47` slopes over the integers.  The two slopes in
`(4)` are not merely locally unobstructed: the next sections give explicit
integer seeds and norm-one units.  Thus `(4)` is an iff inside the finite
universe `(3)`, not a bounded-height search.

## 3. The new `(14,23)` ray

For `(m,n)=(14,23)`, equation `(2)` becomes

```text
W^2-44965u^2=676.                                      (6)
```

It has the positive odd seed

```text
(W_0,u_0)
=(16715091502370452792679,
  78826357654129385469).                               (7)
```

The positive norm-one unit

```text
P+Q sqrt(44965),                                       (8)
P=1062835917426709745982924462536293407897885896493741491943228169,
Q=5012206133771469409285758349474043581500293093464114075203972
```

has `P` odd and `Q` even.  Define

```text
W_(j+1)=P W_j+44965 Q u_j,
u_(j+1)=Q W_j+P u_j.                                  (9)
```

Then `(6)` and odd parity persist, while `W_j,u_j` grow strictly.  Put

```text
d_j=529u_j^2+2,
e_j=u_jW_j,
a_j=7406u_j^3+51u_j,
q_j=103684u_j^4+1036u_j^2+1.                         (10)
```

Coefficient comparison gives

```text
a_j^2+2=d_jq_j,
4q_j-d_j^2=3e_j^2.                                   (11)
```

The seed norm is positive, so `W_0/u_0>sqrt(44965)`, while direct integer
comparison gives `W_0/u_0<529`.  Under `(9)` the ratio decreases toward
`sqrt(44965)`.  Therefore

```text
0<e_j<529u_j^2<d_j                                    (12)
```

for every `j`.  Hence

```text
x_j=(d_j+e_j)/2 > y_j=(d_j-e_j)/2 > 0,
x_j^3+y_j^3=a_j^2+2=Q_((a_j-1)/2).                    (13)
```

This is a second infinite positive ray.  Its first displayed point is

```text
x_0=2302290678332599157988107845344339684655411,
y_0= 984700897345246967397057762882950552473960,
r_0=1813711014853445365774251238977913323291598894943013194651502886.
                                                               (14)
```

The huge seed in `(7)` explains why the earlier search through `u<=5000`
found only THM-3375's `(26,29)` ray.  No minimality claim for `(7)` is needed
or made.

## 4. Two genuinely different interior rays

For `(26,29)`, THM-3375 gives

```text
W^2-522261u^2=2692,
(W_0,u_0)=(5059,7).                                   (15)
```

The limiting positive-chamber ratios of `e/d` on the two rays are

```text
sqrt(44965)/529,
sqrt(522261)/841=sqrt(621)/29.                        (16)
```

They are unequal by exact cross-squaring.  Thus `(14)` and THM-3375 are
distinct asymptotic rays on the positive cubic surface, not two seeds on the
same Pell orbit.

## 5. Audit and scope

The standard-library companion reconstructs all `47` slopes, every modular
point set, the first-obstruction partition `35+2+3+2+2+1=45`, both seed norms,
both unit norms, six iterates of each ray, all parity/factor/square/cube and
positivity identities, strict depth growth, `(14)`, and the distinct limiting
ratios.  Normal and optimized outputs are byte-identical.  LF-normalized
source/output hashes are

```text
7ec10b47139f10117a9f80155ac040ac864bb329ea5687ac1c890bec27529d30
0e5660b462b6f656b8cf09e26e5f6826e74fe626535cb8188fa3fd0f16fa4a44.
```

The atlas is exact only for the reduced-seed ansatz `(1)` and denominators
`n<=29`.  It does not classify all positive Berggren/two-cube points, prove
either displayed seed minimal, count the rays, or give an asymptotic.  It has
no LRC, FC, JC, or AMM consequence.
