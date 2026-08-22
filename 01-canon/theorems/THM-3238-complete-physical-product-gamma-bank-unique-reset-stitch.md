---
id: THM-3238
title: "Complete physical product-Gamma bank unique-reset stitch"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the
  support-(1,3), bank-I2 product-Gamma selector, one explicit positive
  twenty-two-row combination of lawful coarsening-upset responses vanishes
  at Q=(1,3,3,4,5,6,7,8) and is strictly negative on every other nonempty
  physical pole submultiset.  The complete 4,319-state selector cone is
  therefore the singleton delta_Q at every horizon D>=14 and in the
  all-degree intersection.  This is one fixed support/bank face, not arbitrary
  radial coefficients or the Gaussian Moment Conjecture itself.
source: root/multiscale-newton-flag/product-gamma-width3/2026-08-03
audit: >
  The cache-free exact companion reconstructs all 4,319 physical states from
  the pinned product-Gamma formula, verifies all twenty-two upsets and degree
  normalizers, clears the positive rational multipliers primitively, and
  evaluates all 4,319 integer coordinates.  It also replays the 173-state
  THM-3216 tail resurrection, proves the exact degree-five two-endpoint no-go,
  and exhausts all 486 principal upsets (28 point correctly at both endpoint
  hostiles; none repairs the full bank alone).  Independent proof/typing and
  separate-cache audits accept every implication and exact object.  A further
  independent hostile audit rebuilt the 4,319-state bank, every upset/minimal
  antichain, the full sign census and singleton-cone inference, and confirmed
  that the first depth-ten hostile Q+{2,4} changes from positive under the old
  functional to strictly negative under the stitched cocircuit.  Fresh
  cache-free normal and optimized replays both reproduce the immutable stored
  transcript and all declared hashes.
depends_on:
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3216-depth-nine-degree-fourteen-unique-reset-face-and-omega-cone-boundary
related:
  - THM-3219-complete-reset-upper-filter-principal-upset-exclusion
  - THM-3222-universal-product-gamma-reset-upper-filter-collar
  - THM-3244-unique-reset-exposure-deletion-graph-nonmorse-boundary
script: 04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py
output: 05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out
script_sha256: 201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa
output_sha256: 77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b
hash_basis: LF-normalized bytes
---

# THM-3238 -- complete physical product-Gamma bank unique-reset stitch

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3216 exposed the reset through depth nine but its eighteen-row functional
turns positive on 173 deeper states.  THM-3222 supplies universal principal
collars, but principal collars do not glue by themselves.  The missing object
is an overlap cocircuit: eight nonprincipal near-fine upsets, combined with
eleven inherited templates and three principal rows, orient the entire
physical submultiset bank coherently.

This is a global exposed-face theorem, not a local gradient theorem.  A unique
supporting maximum does not by itself give a deletion flow, discrete Morse
contraction, or transport to another support pair.

## 1. Complete physical bank

Use THM-3216's support `(1,3)` and signed bank `I2`:

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1),
Q=(1,3,3,4,5,6,7,8).                                    (1)
```

Let `S(P)` be the set of all nonempty submultisets of `P`.  The pole
multiplicities are `4,3,2,2,2,1,1,1`, so

```text
|S(P)|=(4+1)(3+1)(2+1)^3(1+1)^3-1=4319.                 (2)
```

The exact depth census is

```text
(8,33,93,200,348,507,631,678,
 631,507,348,200,93,33,8,1).                             (3)
```

Thus `(2)` is the complete physical bank, not a depth truncation.  For a
partition `mu` of `N`, retain the response

```text
G_N^sigma(mu)
 =Phi^sigma(h_N)m_mu[Q-sigma]
  -Phi^sigma(m_mu)h_N[Q-sigma].                          (4)
```

For every nonempty proper coarsening upset `U`, write

```text
r_(N,U)(sigma)=sum_(mu in U)G_N^sigma(mu).               (5)
```

Every selector-feasible probability law has nonnegative expectation against
each row `(5)` by THM-3127.

## 2. Twenty-two lawful rows

Write `1^k` for `k` trailing ones.  The last column gives the exact positive
factor `a_i`; the response multiplier is

```text
theta_i=a_i/M_(N_i),                                     (6)
```

where the degree normalizers are those in `(8)` below.  The first eleven row
templates occur in THM-3216, rows 12--14 are principal, and rows 15--22 are
the nonprincipal stitch.

| `i` | `N` | `|U_i|` | minimal generators | `a_i` |
|---:|---:|---:|---|---:|
|1|5|1|`(5)`|`7/916609201`|
|2|14|130|`(3,2,1^9);(2,2,2,1^8)`|`7177018/925091523`|
|3|10|40|`(3,1^7);(2,2,1^6)`|`113613/103885675`|
|4|12|74|`(2,2,1^8)`|`14821346/514459311`|
|5|14|128|`(4,1^10);(3,2,1^9);(2^6,1,1)`|`1800950/915068561`|
|6|14|132|`(2,2,1^10)`|`234911529/973260547`|
|7|14|129|`(5,1^9);(3,3,1^8);(2,2,2,1^8)`|`12760257/713572348`|
|8|14|127|`(5,1^9);(4,2,1^8);(3,3,1^8);(2^4,1^6)`|`1171444/999246305`|
|9|10|41|`(2,1^8)`|`1806523/275215201`|
|10|11|55|`(2,1^9)`|`5625173/775067384`|
|11|11|54|`(3,1^8);(2,2,1^7)`|`3559607/824858825`|
|12|13|98|`(2,2,1^9)`|`30666900/415449457`|
|13|13|100|`(2,1^11)`|`7423537/97767608`|
|14|14|134|`(2,1^12)`|`386487777/737119657`|
|15|5|5|`(3,1,1);(2,2,1)`|`2/599801583`|
|16|9|28|`(3,1^6);(2,2,1^5)`|`5467/645985968`|
|17|11|51|`(3,2,1^6);(2,2,2,1^5)`|`81670/325987229`|
|18|13|96|`(4,1^9);(3,2,1^8);(2^4,1^5)`|`1928501/649376836`|
|19|13|97|`(4,1^9);(3,2,1^8);(2,2,2,1^7)`|`1131782/357891057`|
|20|14|122|`(5,1^9);(4,2,1^8);(3,3,3,1^5);(2^5,1^4)`|`9291/416817778`|
|21|11|53|`(3,1^8);(2,2,2,1^5)`|`768511/977040433`|
|22|14|126|`(5,1^9);(4,2,1^8);(3,3,1^8);(3,2,2,2,1^5);(2^5,1^4)`|`452405/829943856`|

Every displayed set is checked directly to be a nonempty proper coarsening
upset with the displayed minimal antichain.  The unchanged exact degree
normalizers on the complete bank are

```text
M_5 =1300713120,             M_6 =591017995680,
M_7 =67528413334656,         M_8 =9558239814808320,
M_9 =1141773245941342464,    M_10=98177407096144199040,
M_11=6402822693369065077104,
M_12=623312546828577253020720,
M_13=39394150378793250375693600,
M_14=2528571670236939601303479360.                       (8)
```

## 3. Primitive exact dual and exhaustive signs

Clear the denominators of the twenty-two `theta_i` in `(6)` and divide their
integer tuple by its gcd.  This gives a positive primitive vector `c` with
coefficient bit lengths from 887 through 936.  Its full decimal tuple is
printed in the immutable companion output; its digest is

```text
SHA256(c_1,...,c_22)
=bf09da9f72279700fed46afd02e88b9b4798ab228d349ddf665b944bfb50fb0b. (9)
```

The rational table `(6)` is itself a complete exact specification of that
primitive tuple; the output records every decimal coordinate as an additional
copy/paste audit surface.  Define

```text
H(sigma)=sum_(i=1)^22 c_i r_(N_i,U_i)(sigma).             (10)
```

Exact integer evaluation on every state in `S(P)` gives

```text
H(Q)=0,
H(sigma)<0                         for all sigma!=Q.       (11)
```

The negative/zero/positive census by depth is

```text
(8/0/0,33/0/0,93/0/0,200/0/0,348/0/0,507/0/0,
 631/0/0,677/1/0,631/0/0,507/0/0,348/0/0,200/0/0,
 93/0/0,33/0/0,8/0/0,1/0/0).                             (12)
```

The closest negative state is the singleton `(1)` and the farthest is `(8)`.
Their exact integer coordinates are printed in the output; in particular the
strict margin is not inferred from floating point.  The ordered state-bank
digest is

```text
c060e22f900232b608ad3bce1f6b24cae51b6eb45c138b679f4c698fbad2c6a2. (13)
```

## 4. Singleton consequence at every horizon

Let `lambda` be feasible at any horizon `D>=14`.  Each row used in `(10)` is
then lawful, so

```text
sum_sigma lambda_sigma H(sigma)>=0.                       (14)
```

But `(11)` and `lambda_sigma>=0` make the left side nonpositive, with equality
only if `lambda` is supported at `Q`.  Since the quotient alphabet `Q-Q` is
empty, `delta_Q` satisfies every degree.  Therefore

```text
C_D^(physical)={delta_Q}             for every D>=14,
intersection_(D>=14) C_D^(physical)={delta_Q}.            (15)
```

This strictly supersedes THM-3216's depth-at-most-nine conclusion for this
one support/bank face.  THM-3222 remains the reusable local-collar theorem
across all 230 maintained faces.

## 5. Why one collar cannot do this

Let `H_9` be THM-3216's promoted eighteen-row functional.  On depths 10--16
its exact negative/zero/positive census is

```text
(488/0/19,313/0/35,155/0/45,52/0/41,
 9/0/24,0/0/8,0/0/1).                                   (16)
```

Thus precisely 173 tail states resurrect.  The universal degree-five
principal row has the global identity

```text
r_5(sigma)=1440(p_5[Q]-p_5[sigma]).                       (17)
```

Both `P` and `P\{8}` have `H_9>0`, while

```text
r_5(P)=-6117120,
r_5(P\{8})=41068800.                                    (18)
```

Hence no real multiple of `r_5` repairs both endpoint hostiles.  More
strongly, the companion exhausts all 486 principal upsets in degrees 5--14.
Exactly 28 point negatively at both endpoints, but none admits a nonnegative
one-row multiplier making `H_9+lambda r` strict on all 4,318 nonreset states.

The eight nonprincipal rows are therefore not cosmetic numerical slack.  They
are overlap data needed to glue the principal collars.  In holotopy language,
the exposed endpoint section globalizes, while a one-chart collar atlas does
not; the nonprincipal antichains are the finite transition cocycle.

## 6. Scope

This theorem proves the complete physical selector bank only for support
`(1,3)` and product-Gamma bank `I2`.  It does **not** prove the analogous
statement for the other 229 maintained faces, arbitrary radial polynomial
coefficients, NC2 in full generality, or the Gaussian Moment Conjecture.  It
also proves no monotone one-pole deletion route: a globally exposed singleton
need not be a local discrete-Morse maximum.

QED.
