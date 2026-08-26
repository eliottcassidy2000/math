---
id: THM-4228
title: "Common-gcd two-outsider periodic-observable Haar ray"
status: >
  PROVED ANALYTIC LEMMAS + FINITE-EXACT + INDEPENDENTLY AUDITED FIXED-POOL
  RESULT. For the displayed thirty-label pool, every nine-body is Haar-safe
  after adjoining any two distinct positive outsiders whose gcd is at least
  3467. The analytic proof treats the two outsiders as one periodic
  observable and proves the sharp universal density floor 66/91. The number
  3467 is sharp only for the sufficient activation-deck certificate, not for
  literal Haar safety. THM-4150 transfers the body result to infinite
  odd-tail LRC(14) families. Arbitrary finite pair entry and LRC(14) remain
  OPEN.
source: codex-lrc14-common-gcd-ray-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
related:
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
primary_script: 04-computation/lrc14_common_gcd_two_outsider_periodic_observable_ray_thm4228.cpp
primary_output: 05-knowledge/results/lrc14_common_gcd_two_outsider_periodic_observable_ray_thm4228.out
independent_audit_script: 04-computation/lrc14_common_gcd_two_outsider_periodic_observable_ray_independent_audit_thm4228.cpp
independent_audit_output: 05-knowledge/results/lrc14_common_gcd_two_outsider_periodic_observable_ray_independent_audit_thm4228.out
primary_script_sha256: 917184e996dffc038ef00a551625da7055522d993be8ac5bdaeb69e87cd7d340
primary_output_sha256: 78399448d2ba9bcb3fef31bbcdb19ea2e55dd6cf5129abd71aaefa1643f34d44
independent_audit_script_sha256: 55ed8cb6fee4e24f6a0365c0cf4f73eb5ae59c04f5e7e996179d2097c2697ded
independent_audit_output_sha256: 2d75e5c1ef79b458888156da2e8adc2e81020fbf1e8a5ff361c2ad541c3ee6ee
primary_include_sha256: 25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d
independent_include_sha256: 58817d5f5e1a8cc07384f3ea82a1feb221af37ab0907afde890ab4fbdd949137
hash_basis: raw LF bytes
primary_audit: >
  PASS. The canonical pool-cell implementation groups 2,721 failure and
  adjacency atoms, performs 152,170,990 colex zeta contributions, enumerates
  all 5,852,925 depth-eight repairs and all 14,307,150 nine-bodies, and freezes
  the complete strict-edge semantic ledger. It also directly rescans the
  activating witness and checks 79,851 primitive-pair overlap controls.
independent_audit: >
  ACCEPT. A separately implemented fixed-wall geometry uses reverse-colex
  storage, recursive subset enumeration, and 2,493,648,782 ungrouped cell and
  adjacency contributions. It reproduces the semantic ledger, activation
  filtration, every body verdict, and the exact boundary witness. Literal
  joint-wall controls independently test 5,021 primitive pairs, including the
  sharp 1:13 density equality. Optimized and Clang UBSan output streams match.
---

# THM-4228 -- common-gcd two-outsider periodic-observable Haar ray

**PROVED ANALYTIC LEMMAS + FINITE-EXACT + INDEPENDENTLY AUDITED FIXED-POOL
RESULT; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Retain the thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive label set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

> **Common-gcd two-outsider theorem.** Let `q,r` be distinct positive
> integers with
>
> ```text
> gcd(q,r)>=3467.                                       (3)
> ```
>
> Then for every `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {q,r}))>=4/63.                          (4)
> ```

Equivalently, for every integer `g>=3467`, every two distinct positive
integers `u,v`, and every `B in binom(P,9)`,

```text
mu(G_(B union {gu,gv}))>=4/63.                         (5)
```

The labels `q,r` in `(3)` are automatically outside `P`, since both are at
least their gcd and `max(P)=290`. By THM-4150, for every positive integer `c`
and every two distinct positive odd integers `a,b`, some `x in R/Z` satisfies

```text
min_(w in 2c(B union {q,r}) union {a,b})||wx||>=1/14. (6)
```

The eleven even body labels and two odd tails in `(6)` are distinct, so these
are infinite thirteen-speed LRC(14) families. THM-4231 separately covers the
entire northeast quadrant of distinct pairs `q,r>=17548`; the present theorem
reaches below that quadrant along large-gcd rays. It proves neither arbitrary
finite pair entry nor full LRC(14).

The closest proved mechanism is THM-4170's endpoint discrepancy bound for one
comb. The canonical hostile is THM-4207's comparable pair `(50,51)`, where
separate marginal decks fail to compose already at depth four. The corrected
near miss is that a sufficient-deck cover signals certificate failure, not an
unsafe body; THM-4203/4207/4211 give literal positive countercontrols. The
least-used sidecar is the common divisor retained when the pair is factored as

```text
(q,r)=g(u,v).                                           (7)
```

The live concept board was

```text
depth-eight base repair | gcd quotient | joint pair observable
endpoint cocycle | activation filtration | predecessor cover. (8)
```

Sequentially adjoining `q` and `r` creates order-`q` interval complexity.
Factoring `(7)` instead turns `G_q intersect G_r` into a single pullback at
frequency `g`; the primitive pair shape stays inside one periodic observable.

## 2. Sharp overlap floor for two integer combs

Write, up to null endpoints,

```text
D_n={y in R/Z:||ny||<1/14}=T\G_n.                     (9)
```

Every `D_n` has measure `1/7`. We first prove the extra overlap forced by the
integer-comb geometry.

> **Two-comb overlap lemma.** If `u,v` are distinct positive integers, then
>
> ```text
> mu(D_u intersect D_v)>=1/91,                         (10)
> mu(G_u intersect G_v)>=66/91.                        (11)
> ```
>
> The constant in `(10)` is sharp: equality holds exactly when the reduced
> ratio is `u:v=1:13` or its swap.

To prove the lemma, divide by `gcd(u,v)`. Indeed, multiplication by the gcd is
Haar preserving, so the two intersection measures do not change. We may
therefore write `1<=u<v` and `gcd(u,v)=1`.

The components of `D_u` and `D_v` are centered respectively at `m/u` and
`n/v`, with half-widths `1/(14u)` and `1/(14v)`. The map

```text
(m,n) |-> vm-un  mod uv                                (12)
```

is a bijection from their `uv` ordered component pairs to `Z/(uv)`: reducing
an equality in `(12)` modulo `u` and modulo `v` recovers `m` and `n`
separately. Choose the centered representative in the half-open interval
`[-uv/2,uv/2)`. If `uv` is even, the alternative antipodal representative
has the same absolute value and in any case contributes zero, because
`u+v<7uv`. The same inequality says that the sum of the two component
half-widths is below `1/2`, so no second circular wrap can contribute. For
this representative `k`, the overlap length is

```text
ell_k=[min(2u,u+v-14|k|)]_+/(14uv),                   (13)
```

where `[z]_+=max(z,0)`. Thus, with

```text
w(t)=[min(2u,u+v-14t)]_+,
gamma=mu(D_u intersect D_v),                           (14)
```

we have the exact formula

```text
14uv gamma=w(0)+2 sum_(k>=1) w(k).                    (15)
```

The function `w` is nonnegative and decreasing, and its trapezoidal area is

```text
integral_0^infty w(t)dt=uv/7.                         (16)
```

If `v>=16`, the left Riemann sum and `(16)` give

```text
w(0)+2 sum_(k>=1)w(k)
 =2 sum_(k>=0)w(k)-2u
 >=2uv/7-2u
 >=2uv/13.                                             (17)
```

For `v<=13`, the term `w(0)=2u` already gives the last bound in `(17)`. If
`v=14`, the three terms at `k=0,+-1` sum to `4u`; if `v=15`, they sum to
`4u+2`. These also exceed `2uv/13`. The last inequality in `(17)` is strict
for `v>=16`, and the `v=14,15` estimates are strict. In the remaining range,
equality forces `v=13`; then `w(1)=[u-1]_+`, so equality further forces
`u=1`. Conversely, at `(u,v)=(1,13)`, only `k=0` has positive length and
`(15)` equals `2`, so `gamma=1/91`. Dividing by `14uv` proves `(10)` and its
equality statement.

Finally, inclusion-exclusion gives

```text
mu(G_u intersect G_v)
 =1-1/7-1/7+gamma
 >=66/91,                                              (18)
```

which proves `(11)`.

## 3. The periodic-observable endpoint cocycle

Let `A` be a measurable subset of the circle with `mu(A)=beta`, and define
the centered primitive on the real line by

```text
H_A(t)=integral_0^t (1_A(s)-beta)ds.                   (19)
```

The zero mean makes `H_A` one-periodic. For an oriented arc `I` of length
`ell in [0,1]`,

```text
integral_I(1_A-beta)=mu(A intersect I)-beta ell.       (20)
```

If the right side is positive, it is at most

```text
min(beta,ell)-beta ell<=beta(1-beta).                  (21)
```

Choose the orientation from a point where `H_A` is minimal to one where it is
maximal. Equations `(20)--(21)` prove the sharp general estimate

```text
osc(H_A)<=beta(1-beta).                                (22)
```

Sharpness is attained when `A` is a single arc: taking `I=A` in `(20)` gives
exactly `beta(1-beta)`.

Now let `U` be a union of `c` positive-length circular intervals and let
`m_g(y)=gy` on the circle. For each interval component, choose an oriented
real lift `[x,z]` with `x<z<=x+1` (equivalently, split a component at zero).
Periodicity gives

```text
integral_x^z (1_A(gy)-beta)dy
 =(H_A(gz)-H_A(gx))/g.                                 (23)
```

Summing `(23)` over the components, while ignoring null isolated points,
proves the periodic-observable discrepancy lemma

```text
mu(U intersect m_g^(-1)(A))
 >=beta mu(U)-c osc(H_A)/g.                            (24)
```

Take `A=G_u intersect G_v`. By `(18)`, `beta>=66/91`; also `beta<=6/7`, and
`beta(1-beta)` is decreasing for `beta>=1/2`. Therefore

```text
osc(H_A)<=1650/8281.                                   (25)
```

The factorization `(7)` is now exact:

```text
G_(gu) intersect G_(gv)=m_g^(-1)(G_u intersect G_v).  (26)
```

Combining `(18)` and `(24)--(26)` yields, for every such `U`,

```text
mu(U intersect G_(gu) intersect G_(gv))
 >=(66/91)mu(U)-1650c/(8281g).                         (27)
```

No independence or coprimality assumption is hidden in `(27)`. Coprimality
only chooses a unique representative in `(7)`.

## 4. The depth-eight activation deck

The common fixed-pool wall denominator is

```text
D=18,241,159,416,480.                                  (28)
```

For `R in binom(P,8)`, put

```text
U_R=G_(P\R),
mu(U_R)=m_R/D,                                         (29)
```

and let `c_R` be the number of positive-length circular interval components
of `U_R`. Define

```text
T_R=297m_R-26D.                                        (30)
```

This normalization records exactly the worst-density slack:

```text
(66/91)(m_R/D)-4/63=2T_R/(819D).                       (31)
```

For every repair with `T_R>0`, define its sufficient activation time

```text
kappa(R)=ceil(7425 c_R D/(91T_R)),                     (32)
```

and the nested common-scale deck

```text
E_g={R in binom(P,8):T_R>0 and kappa(R)<=g}.           (33)
```

The constants in `(32)` are the reduced form of

```text
(1650/8281)(c_R/g)<=2T_R/(819D).                       (34)
```

Because `g` is an integer and `T_R>0`, `kappa(R)<=g` is equivalent to the
unrounded ratio in `(32)` being at most `g`. Equality is admissible: the
target Haar bound is the non-strict inequality `mu>=4/63`.

Consequently, `(27)` and `(31)--(34)` give the load-bearing implication

```text
R in E_g
 =>mu(G_((P\R) union {gu,gv}))>=4/63                  (35)
```

for every pair of distinct positive `u,v`. The deck `E_g` is sufficient. It
can omit literally lawful repairs because `(18)` and `(22)` deliberately
discard pair-specific density and endpoint phase.

## 5. Exact transition at 3467

Two independent exact programs prove

```text
|binom(P,8)|=5,852,925,
#{R:T_R>0}=5,052,990,
#{R:T_R=0}=0,                                          (36)

|E_3466|=821,120,
|E_3467|=821,923,                                      (37)

tau(E_3466)<=9,
tau(E_3467)>9.                                         (38)
```

For `(38)`, both programs inspect every one of the

```text
binom(30,9)=14,307,150                                 (39)
```

labelled nine-bodies. At `g=3467`, every body misses an active repair. Any
transversal of size at most nine could be extended inside `P` to a nine-body,
so the all-nine-body scan proves `tau(E_3467)>9`, not merely the absence of an
exactly-nine cover. At the predecessor, the explicit body

```text
B_*={88,95,170,193,240,252,264,286,290}               (40)
```

meets every edge of `E_3466`. An activating repair disjoint from it is

```text
R_*={8,16,85,132,143,145,168,176}.                    (41)
```

The exact geometry and slack of `(41)` are

```text
m_(R_*)=1,920,678,988,176,
c_(R_*)=224,
T_(R_*)=96,171,514,659,792,                            (42)

7425 c_(R_*)D/(91T_(R_*))
 =210,474,916,344,000/60,714,340,063
 =3466.642577... .                                     (43)
```

The sufficient lower-bound excess over `4/63` is exactly

```text
-19,506,842,821/8,172,402,168,912,345 at g=3466,
 21,700,654,421/16,349,520,092,105,655 at g=3467.      (44)
```

Thus `R_*` activates at 3467. The cover `(40)` and the exhaustive no-cover
scan prove that 3467 is the exact transition for the nested certificate
`E_g`. Equation `(44)` is not a witness of literal danger at 3466.

The complete strict-edge ledger is

```text
b270307777887d42.                                      (45)
```

Both implementations freeze `(36)--(45)` and agree on 44,319,858,349 body-to-
edge incidence checks, with maximum first-disjoint position 821,129.

## 6. Proof of the theorem

Fix distinct `q,r` satisfying `(3)`, write

```text
g=gcd(q,r), q=gu, r=gv, gcd(u,v)=1,                   (46)
```

and fix `B in binom(P,9)`. Since the activation decks are nested and
`g>=3467`, we have `E_3467 subset E_g`. By `(38)`, choose
`R in E_3467 subset E_g` disjoint from `B`. Equation `(35)` gives

```text
mu(G_((P\R) union {q,r}))>=4/63.                      (47)
```

Since `R intersect B=empty`,

```text
B union {q,r} subset (P\R) union {q,r}.               (48)
```

Safe-set inclusion reverses label inclusion, so `(47)--(48)` prove `(4)`.
The multiplication map on the circle preserves Haar measure; THM-4150 then
proves `(6)`. **QED.**

## 7. Connection contract, audit, and failure boundary

The out-of-the-box bridge is

```text
source:       a pair (q,r) with retained common divisor g, and base U_R
target:       a lawful joint depth-eight repair
map:          (q,r)=g(u,v), A=G_u intersect G_v,
              G_q intersect G_r=m_g^(-1)(A), then endpoint cocycle (24)
preserved:    exact pair correlation inside A, labelled repair, base mass,
              component count, common scale, target threshold
destroyed:    literal endpoint addresses and pair-specific beta/oscillation
sidecar:      gcd g and c_R; literal joint walls below the certificate boundary
decisive test: tau(E_g)>9.                              (49)
```

The primary computation groups fixed-pool failure and adjacency atoms and
uses an ordinary-colex zeta transform plus Gosper subset enumeration. The
independent audit rebuilds the walls separately, uses reverse-colex storage
and recursive enumeration, and scatters every cell and adjacency contribution
without the primary grouping. It directly rescans `(41)`. A separate literal
joint-wall bank checks 5,021 primitive pairs through denominator 128; the
observed maximum primitive oscillation there is `990/8281` at ratio `1:13`.
That finite observation is only a hostile control: the proof uses the global
analytic ceiling `1650/8281` from `(22)`.

Reproduction from `04-computation/` is:

```bash
clang++ -std=c++20 -O3 -DNDEBUG \
  lrc14_common_gcd_two_outsider_periodic_observable_ray_thm4228.cpp \
  -o /tmp/thm4228-primary
/tmp/thm4228-primary

clang++ -std=c++20 -O3 -DNDEBUG -pthread \
  lrc14_common_gcd_two_outsider_periodic_observable_ray_independent_audit_thm4228.cpp \
  -o /tmp/thm4228-independent
/tmp/thm4228-independent
```

The boundary is strict:

1. `3467` is sharp only for the sufficient universal-density/primitive-
   oscillation activation deck. Exact pair geometry can succeed earlier.
2. The theorem says nothing about distinct pairs with gcd below `3467`.
   THM-4231 separately handles such pairs when both labels are at least
   `17548`; smaller-coordinate entry remains open in general.
3. The diagonal `u=v` analytically collapses to one comb and does not provide
   two distinct outsider labels.
4. No arbitrary LRC(14) row is proved to enter this fixed-pool chart.
5. LRC(14) remains open.
