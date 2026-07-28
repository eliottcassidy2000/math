---
id: THM-2791
title: "Full arm-orbit transfer and lower-central chord"
status: >
  RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
  INDEPENDENT AUDIT.  The complete clock-one THM-2782 address orbit has,
  in each raw target column 3 through 11, exactly two positive marked
  cylinders.  Their gap is the graded chord
  Z1 Z2^10 Z3 Z4^11 Z5, with first lower-central transgression 10 in
  Z2/Z3.  The coefficient has full period and every address conductor;
  its canonical pushforward to the 169-address quotient is the positive
  two-point chain c z^6(1+z), whose permitted scalar or coefficient-
  Bockstein normalization is a central group-algebra unit.  This is
  transfer along a genuine partial graded-chord germ, not pointwise
  descent, pure-Z1 packet covariance, or an endpoint-origin allocation.
  Until independent promotion, no proved result may depend on this
  candidate.
source: root/graded-semantic-arm-transfer-2026-07-28
depends_on:
  - THM-2782-semantic-arm-right-wing-local-unit-and-endpoint-deck-boundary
  - THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate
related:
  - THM-2657-odometer-carry-root-lift-nonsplit-extension-and-cech-cocycle
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2712-semantic-following-congruence-lock-and-address-coboundary-descent
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
script: 04-computation/lrc14_full_arm_orbit_lower_central_chord_thm2791.py
output: 05-knowledge/results/lrc14_full_arm_orbit_lower_central_chord_thm2791.out
script_sha256: 2824d62c237fd9ac831d23236e6987ecabe96bebd68ba37a9abd0bb685ad0716
output_sha256: 9f2b8e69b9de430f201adb7758f98fff7bf505c5bd792b03f40b3ee7c9f46edd
hash_basis: LF-normalized bytes
---

# THM-2791 -- full arm-orbit transfer and lower-central chord

**RESERVED / PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT / AWAITING
INDEPENDENT AUDIT.**

THM-2782 finds a three-point physical central-direction segment and one
same-`Omega` descent failure.  The complete orbit is substantially more
rigid.  Its main raw target columns are two-point coefficient chains with
full cyclotomic spectrum, and the remote second point differs from the
canonical central successor by one explicit lower-central word.  After
the permitted scalar or coefficient-Bockstein normalization, quotient
pushforward turns that chord into a unit on the THM-2779 central cycle.

The positive result is a **transfer**, not a descent.  The two coarse
neighbours come from remote high-digit representatives.  There is an exact
partial physical translation germ between them, but it has the full graded
gain, not the pure central gain.

## 1. Complete physical coefficient universe

Put

```text
p=13,                      R=p^6=4,826,809,
N=p^5=371,293,             n_0=3,454,614,
c=790,161,473,087,466,480.                              (1)
```

Fix physical clock one and lawful source-target label `sigma=0`.  For
`tau in F_13` and every address offset `k in Z/R`, define `F_tau(k)` by the
same order of physical operations as THM-2782:

1. construct the full `Rwing` packet with both relative-present
   complements;
2. push it to the target coordinate;
3. retain the `c3` root-one half;
4. restrict to the open target cylinder based at address `n_0+k`; and
5. take the actual carry-six/root-one coefficient of
   `D^6 Q_(3,{1,2})`.

No semantic or address value is assigned after integration.  The companion
constructs the interval carrier first and then scans every one of the
`R` address cylinders.

For each

```text
tau in {3,4,...,11},
```

the complete answer is

```text
F_tau(k)=c   for k in {0,689364},
F_tau(k)=0   otherwise.                                 (2)
```

The two physical addresses are

```text
n_0=3,454,614,              n_1=4,143,978.              (3)
```

Both open cylinders retain:

```text
rail eight;
clock one;
E3 -> Q_(3,{1,2});
both relative-present complements;
all lawful target labels tau=3,...,11 at sigma=0;
target root one;
actual predecessor carry six.                            (4)
```

The semantic stability radius equals the chosen open-cylinder radius at
both selected centers, so the record is constant on each open cylinder.
At each cylinder the right boundary point changes the semantic record and
is load-bearing; both endpoints remain omitted under the strict-open
convention.  Each positive cell has exact weighted mass

```text
c/R=60,781,651,775,958,960/371,293.                      (5)
```

## 2. The physical chord and its first transgression

On the full address group let

```text
O:n -> n+1,                   X:n ->14n,
Z_r=O^(13^r).                                             (6)
```

Direct substitution gives

```text
X O X^(-1)=O^14,             [X,Z_r]=Z_(r+1),            (7)
```

where the commutator convention is
`[A,B]=A B A^(-1) B^(-1)`.

The gap in `(2)` is

```text
689364=13a,                 a=53028=1+13*4079.           (8)
```

Its base-thirteen digits are

```text
a=(1,10,1,11,1)_13.                                      (9)
```

Therefore the positive physical chord has the exact lower-central word

```text
Z_1^a=Z_1 Z_2^10 Z_3 Z_4^11 Z_5.                       (10)
```

Let `u` denote `Z_1` in `Z[C_(13^5)]`.  Relative to the canonical coarse
central step, the discrepancy factors as

```text
u^a-u
 =u(u^13-1)(1+u^13+...+u^(13*4078)).                    (11)
```

It is divisible by `u^13-1` and not by `u^169-1`.  After division by the
first factor, reduction by `u^13=1`, and coefficient reduction modulo
thirteen, the quotient is

```text
4079u=10u.                                               (12)
```

Thus the first nonzero lower-central transgression is exactly

```text
10 in Z_2/Z_3.                                           (13)
```

THM-2782's `j=13` hostile is the action of `Z_2`: the failure occurs on the
first commutator grade, not in an unspecified high digit.

## 3. An exact partial physical germ

The chord is not merely a coincidence of two integrals.  Its physical
circle translation is

```text
delta=7*689364/13^6=371196/371293.                       (14)
```

Translation by `delta` maps the entire first open target cylinder onto the
second.  Their restricted weighted carriers are single whole-cylinder
pieces with the common rail weight

```text
27,581,135,604.                                          (15)
```

Moreover

```text
13^6 delta=4,825,548,        13^5 delta=371,196          (16)
```

are integers.  Here the delayed phase uses `R=13^6`, the predecessor carry
uses `N=13^5`, and the selected deep/root speed is
`C_3=742586=2N`.  Hence the delayed base, predecessor carry, and root half
are pointwise invariant on the restricted translation.  Both delayed
carry-six pairs are exactly

```text
(0,c).                                                    (17)
```

The selected marked integrand therefore transports along the local
`Z_1^a` germ.  This does not say that any present factor, the whole packet,
or the full circle is covariant under `delta`; both restricted carrier
products happen to equal one on their respective cylinders.  In particular,
it proves no pure-`Z_1` action.

## 4. Exact period and every address conductor

On the central coset write `k=13j`.  Equation `(2)` becomes

```text
f_tau/c=delta_0+delta_53028
        in Z[C_(13^5)].                                  (18)
```

Because `53028` is a thirteen-unit and the group order is odd, this
two-point support has trivial translation stabilizer.  Hence

```text
minimal central-arm period=13^5,
minimal full-address period=13^6.                        (19)
```

For every `1<=r<=5`, the finite difference

```text
(Z_r-1)F_tau
```

has exactly four signed cells and

```text
||.||_1=4c,                   ||.||_2^2=4c^2.            (20)
```

The next element `Z_6` is the identity and its difference is zero.

For a character `chi` of `C_(13^6)`,

```text
Fhat_tau(chi)
 =c(1+chi(689364)^(-1)).                                (21)
```

Here the transform uses the offset coordinate `k` from Section 1.  If one
instead regards the same row as an absolute-address function supported at
`n_0` and `n_0+689364`, its transform acquires the harmless common phase
`chi(n_0)^(-1)`.  The last character value has odd prime-power order and
cannot equal `-1`.
Thus **every one of the `13^6` full-address characters survives**.  The
same argument on `(18)` gives every central-arm character.  In particular,
all primitive conductors

```text
13,13^2,13^3,13^4,13^5                                 (22)
```

survive.

## 5. Positive quotient transfer normalizes to a central unit

Reduction modulo `169` sends the addresses in `(3)` to

```text
85=7+13*6,                    98=7+13*7.                 (23)
```

The address-chain pushforward is therefore

```text
q_!F_tau=c(delta_(7,6)+delta_(7,7))
        =c z^6(1+z),                                      (24)
```

where `z` is the central THM-2779 coordinate at fixed low digit seven.
All twelve nontrivial central characters are nonzero.

The two-point chain is stronger than a nonzero spectrum.  In every
coefficient field `K` of characteristic different from two,

```text
(1+z)(1-z+z^2-z^3+...+z^12)=1+z^13=2                  (25)
```

in `K[C_13]`.  Hence

```text
(1+z)^(-1)=1/2(1-z+z^2-z^3+...+z^12).                  (26)
```

After dividing by the nonzero scalar `c`, or after the inherited
coefficient-Bockstein normalization `(c/13) mod13=2` in characteristic
thirteen, convolution by the normalized version of `(24)` is an
automorphism of the whole thirteen-coordinate central coefficient module.
It is a cyclic module generator.  The raw integral chain in `(24)` is not
itself a unit of `Z[C_13]`.

This is a **pushforward**, not a pullback descent.  Descent would require
constancy on the `13^4`-point fibres and already fails at

```text
F_tau(0)=c,                     F_tau(169)=0.             (27)
```

The adjacent coarse vertices in `(24)` come from physical central lifts

```text
j=0,                         j=53028,
```

whereas the physical adjacent lift `j=1` is empty.  The inverse `(26)` is
signed and generally nonintegral.  Equations `(24)--(26)` do not produce
same-high-digit ancestry, positivity of the inverse, pure-`Z_1` covariance,
or a THM-2779 endpoint-origin allocation.

## 6. Exceptional target and full decoder

The raw `tau=12` column has exactly `121` positive central-arm cells.  Its
pushforward to `C_13`, in the normalized `j mod13` coordinate, is

```text
9N_13+(1+z+z^2+z^3).                                    (28)
```

It has minimal period `13^5`, and all five primitive conductors in `(22)`
survive.

After dividing the complete raw target table by its common positive content
`c`, let

```text
A=delta_0+delta_53028
```

be the raw chord, and let `B` be the exceptional 121-cell profile.  Applying
the displayed integral lift of THM-2771's `K_beta` to this normalized target
table gives

```text
G_q=A_q A+B_q B,                                         (29)

(A_q,B_q)_(q=0)^12 =
((58,5),(64,5),(59,5),(55,7),(48,12),(51,2),(48,8),
 (53,2),(56,9),(50,8),(47,11),(49,0),(55,3)).            (30)
```

Exact cyclotomic block tests show that every one of the thirteen decoded
columns retains every conductor in `(22)`.  Each quotient pushforward is
nonconstant and hence has all twelve nontrivial central characters.

This is a useful target-decoder survivor.  THM-2782's signed `Q3`
projection does not canonically extend from its three-arm augmentation
role to the complete central orbit, so no such extension is claimed.

## 7. Exact invoice and remaining square

The result has the following typed invoice:

```text
source:
  the complete physical clock-one THM-2782 coefficient array on Z/13^6;

target:
  one positive two-vertex coefficient chain on the central Omega cycle
  whose scalar/Bockstein normalization is a group-algebra unit;

map:
  exact address-chain pushforward modulo 169;

preserved:
  raw coefficient and mass, positivity, target semantic word, target root
  one, carry six, lawful labels, relative-present cuts, low digit, and
  every coarse central character;

lost:
  four higher central digits, pure-Z_1 physical adjacency, global packet
  covariance, and endpoint-origin allocation;

positive sidecar:
  one exact same-integrand partial germ with the full graded chord (10);

needed sidecar:
  a natural map from that graded physical germ to the THM-2779
  endpoint-origin central edge, retaining ancestry and current phase.     (31)
```

The cheapest next test is not another local support scan.  It is to ask
whether one marked endpoint/cospan bank transports the address gain

```text
13*(1,10,1,11,1)_13
```

to the endpoint-origin central action, or whether every endpoint edge is
forced to have the pure digit word `(1,0,0,0,0)`.

## 8. Exact companion and scope

Run

```bash
python 04-computation/lrc14_full_arm_orbit_lower_central_chord_thm2791.py
python -O 04-computation/lrc14_full_arm_orbit_lower_central_chord_thm2791.py
```

Both modes byte-match

```text
05-knowledge/results/lrc14_full_arm_orbit_lower_central_chord_thm2791.out.
```

The companion uses exact integers and `Fraction`s and contains no Python
`assert`.  It reconstructs the physical carriers, scans all `13^6`
addresses and the full `13^5` central coset, verifies both positive
cylinders and their semantic/root/carry/mass data, exhausts the exceptional
121-cell column, checks each cyclotomic conductor and all thirteen decoded
columns, proves the group-ring unit identity, and verifies all
lower-central differences and the partial physical germ.

No full physical Heisenberg action, pure central packet covariance,
endpoint-origin allocation, THM-2542 root transition, row exclusion, or
LRC(14) conclusion is proved.

**Awaiting independent audit; not QED.**
