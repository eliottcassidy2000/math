---
id: THM-3547
title: "Positive two-cube slope atlas through denominator 101 and four Pell rays"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT / INDEPENDENTLY AUDITED WITH WORDING
  REPAIRS.  Among
  the 528 primitive parity-correct reduced slopes through denominator 101,
  an exhaustive screen by every prime at most 199 excludes 517 and leaves 11
  explicitly listed screen survivors.  Four survivors carry directly
  certified positive odd Pell orbits, giving four infinite families of
  positive ordered cube pairs with x^3+y^3=(2r+1)^2+2.  No global insolubility or
  everywhere-local-solubility statement is made for the other seven.
source: kps-s185
depends_on:
  - THM-3375-berggren-positive-two-cube-pell-ray
related:
  - THM-3375-berggren-positive-two-cube-pell-ray
  - THM-3376-positive-two-cube-slope-atlas-through-29
companion: 04-computation/berggren_positive_cube_slope_atlas_101_kps_s183.py
output: 05-knowledge/results/berggren_positive_cube_slope_atlas_101_kps_s183.out
script_sha256: bdb8cd4fbd14235ee144d80c4766aed5117097321cae8ec87b1aad1ffff1212c
output_sha256: fb0b50e2384b8834b381c3ded37358d636c7ab475b34d1c5091ff7f3f6162f49
hash_basis: LF-normalized bytes
---

# THM-3547 -- four positive Pell rays survive the slope atlas through 101

**PROVED + FINITE-EXACT + VERIFIED-EXACT / INDEPENDENTLY AUDITED WITH WORDING
REPAIRS.**
THM-3375 found one positive ray and THM-3376 found a second.  The larger exact
atlas shows that the mechanism is neither unique nor visible at small height:
two more rays have enormous first certified seeds.

## 1. Reduced-slope compiler

For coprime integers `(m,n)` with

```text
n odd,        m even,        n/2<m<n,                  (1)
```

put

```text
d=n^2u^2+2,
e=uW,
a=mn^2u^3+(2m+n)u,
q=m^2n^2u^4+2m(m+n)u^2+1.                              (2)
```

The nontrivial identity needed to compile a sum of two cubes is controlled by

```text
3W^2=n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2).                (3)
```

Indeed the first identity below holds identically, while direct expansion gives

```text
a^2+2=dq,
4q-d^2-3e^2
 =u^2[n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2)-3W^2].        (4)
```

Thus for `u!=0`, equation `(3)` is equivalent to `4q-d^2=3e^2`.  Every
certified ray below has positive odd `u`.

For odd `u,W`, define

```text
x=(d+e)/2,        y=(d-e)/2,        r=(a-1)/2.         (5)
```

Then `(4)` implies

```text
x^3+y^3=d(d^2+3e^2)/4=dq=a^2+2=(2r+1)^2+2.            (6)
```

Positivity is the remaining chamber condition `0<e<d`.

## 2. Exact finite prime screen

The companion enumerates the finite universe

```text
3<=n<=101, n odd, m even, n/2<m<n, gcd(m,n)=1.         (7)
```

It contains exactly `528` slopes.  For every prime `p<=199`, and for every
still-live slope, the code exhausts all residue pairs `(u,W) mod p` in `(3)`.
The first-obstruction counts are

| prime | excluded at this first obstruction |
|---:|---:|
| 3 | 395 |
| 5 | 21 |
| 7 | 44 |
| 11 | 11 |
| 13 | 7 |
| 17 | 8 |
| 23 | 2 |
| 29 | 5 |
| 31 | 4 |
| 37 | 3 |
| 41 | 2 |
| 47 | 1 |
| 53 | 2 |
| 59 | 3 |
| 61 | 2 |
| 71 | 2 |
| 79 | 1 |
| 83 | 2 |
| 89 | 2 |

These counts sum to `517`.  The exact eleven screen survivors are

```text
(14,23), (26,29), (26,47), (38,47), (38,53), (50,53),
(50,71), (62,95), (74,95), (74,101), (98,101).        (8)
```

Here and below `(m,n)` is the ordering.  Surviving `(8)` means only that no
obstruction was found modulo a prime at most `199`; it is not an assertion of
solubility over every completion.

## 3. Four direct Pell certificates

For four slopes, dividing `(3)` by three gives

```text
W^2-Du^2=C.                                             (9)
```

The following table supplies an odd positive seed `(W_0,u_0)` and a positive
norm-one unit `(P,Q)`, meaning `P^2-DQ^2=1`.

| `(m,n)` | `(D,C)` | `(W_0,u_0)` | `(P,Q)` |
|---|---|---|---|
| `(14,23)` | `(44965,676)` | `(16715091502370452792679,78826357654129385469)` | `(1062835917426709745982924462536293407897885896493741491943228169,5012206133771469409285758349474043581500293093464114075203972)` |
| `(26,29)` | `(522261,2692)` | `(5059,7)` | `(454592181623521601375,629039857366305528)` |
| `(26,47)` | `(364485,2116)` | `(1576376215602550757032908052542957159292644701,2611079190314729909962215655566812890419809)` | `(24107434053023296515322082382772446705687987537902965213862639312622666158519,39931089269621661093589841956801775898467910099449466227642200335186151476)` |
| `(98,101)` | `(95940405,38404)` | `(220788520986775462919614657060547064697173954066965079893,22541131703851573370870656292410483656470170855181183)` | `(105189902320645145211266433397204047145119454659522482848985261748433399672601631439105713727025690111,10739233323941537083546286758788229602860630133788513945467635000056068525534702555933010597208288)` |

Every displayed norm identity is checked directly with integers.  Define the
recurrence

```text
W'=PW+DQu,                    u'=QW+Pu.                 (10)
```

It preserves `(9)`.  In all four rows `D,P,W_0,u_0` are odd and `Q` is even,
so `(10)` preserves odd parity.

## 4. The positive invariant cone

Set

```text
L=n^2u-W,                    H=n^2W-Du.                 (11)
```

For every seed in the table, `L,H>0`.  Under `(10)` they satisfy the exact
positive recurrence

```text
L'=PL+QH,                    H'=PH+DQ L.                (12)
```

Thus the cone `L,H>0` is invariant.  In particular,

```text
0<W<n^2u,
0<e=uW<n^2u^2<d,                                      (13)
```

so `(5)` gives `x>y>0` at every depth.  Also `u'>u` for positive data, hence
the orbit supplies infinitely many distinct positive integer solutions of
`(6)`.

Multiplication by the positive unit `P+Q sqrt(D)>1` gives
`W_j/u_j -> sqrt(D)`.  Since `e_j=u_jW_j` and
`d_j=n^2u_j^2+2`, the four limiting chamber ratios are

```text
e/d -> sqrt(85)/23, sqrt(621)/29, sqrt(165)/47,
       sqrt(9405)/101,                                  (14)
```

respectively.  Their squares are distinct rational numbers, so the four rays
are asymptotically distinct.  In fact the reduced slope can be recovered from
any compiled positive ordered pair: with

```text
N=sqrt(x+y-2)=nu,              a=sqrt(x^3+y^3-2),
m/n=(a-N)/(N(x+y)).                                     (15)
```

Thus different displayed slopes cannot produce the same ordered pair, though
no claim is made that their target-value sets are disjoint.  The companion prints and verifies the exact
first `(x,y,r)` member of each ray; the huge members for `(26,47)` and
`(98,101)` explain why bounded-height searches missed them.

## 5. Verification and boundary

The companion is standard-library Python.  It checks the complete universe,
all `46` primes through `199`, every norm and parity identity, six recurrence
steps on every ray, the invariant cone, strict growth, `(4)`, and the exact
cube identity `(6)`.  Normal and optimized executions are byte-identical and
equal the stored output.

This theorem does not classify the seven remaining screen survivors, claim
the displayed seeds are minimal, or prove an asymptotic count of positive
representations.  The sharp next arithmetic target is to separate, for those
seven slopes, three different possibilities that a finite prime screen cannot
distinguish: a later local obstruction, a global norm-class obstruction, or a
positive odd Pell orbit with another enormous regulator.

The independent audit reconstructed all `528` slopes, all first-obstruction
counts, the eleven survivors, every seed/unit norm, the cone recurrence,
parity, growth, cube identities, replay hashes, and limits.  It supplied the
`u!=0` precision in `(4)` and the limit justification above.
