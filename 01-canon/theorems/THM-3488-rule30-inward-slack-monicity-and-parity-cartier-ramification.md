---
id: THM-3488
title: "Rule 30 inward slack monicity, parity-Cartier ramification, and shallow transverse atlas"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  The complete physical inward tail has a unique lowest ballistic Newton
  face, so every target slack polynomial is monic and every global Hasse
  layer has an exact causal onset.  A two-channel Cartier transform gives
  the exact repeated-root law after summing all source depths, and every
  translated-root source coordinate has a shallow nonzero witness.  No
  Rule 30 prize, density result, or unrestricted complexity lower bound is
  claimed.
source: root-rule30-next-targets-20260816
audit: >
  An independent hostile audit rederived the ballistic Newton faces and
  targetwise monicity, the full parity-Cartier inversion and valuation law,
  and every unique owner in the source-Taylor residue atlas.  Expanded exact
  controls covered 1,024 Taylor orders and 4,096 abstract channel pairs;
  ordinary and optimized runs equal the stored transcript byte-for-byte:
  ACCEPT.
depends_on:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas
related:
  - THM-3480-rule30-staircase-transducer-entropy-and-nonrectangular-macroblock-compiler
  - THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum
  - THM-3466-factorial-face-stokes-and-keller-boundary-current
script: 04-computation/rule30_inward_slack_cartier_thm3488.py
output: 05-knowledge/results/rule30_inward_slack_cartier_thm3488.out
script_sha256: be219b6d26a5cdc472e03bab6c282eebdfcdfdd674699d6087b6447099b4df0e
output_sha256: c5783f48ef5992de3209939f0c457e466db72ffab058e121f7131942f8bb0310
hash_basis: raw bytes
---

# THM-3488 -- Rule 30 inward slack monicity, parity-Cartier ramification, and shallow transverse atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The first three radial source strips form THM-3471's exact unmarked boundary
circuit.  This theorem starts immediately inside that circuit, at source
depth `u=3`, and keeps the Green slack marker.  Three exact structures survive
the sum over **all** inward depths:

1. a unique lowest ballistic Newton face;
2. an invertible even/odd Cartier atlas for the repeated-root filtration; and
3. a uniformly shallow leading witness in every translated-root source
   coordinate.

The structures are marked and filtered.  Specializing the marker to one can
and does destroy them.

## 1. Inheritance and conventions

The closest mechanism is THM-3476's fixed-strip ramification law.  The
canonical hostile is THM-3471's nonzero three-strip current which becomes zero
at `q=1`.  The corrected near miss is therefore to sum unbounded source depth
only while retaining the native slack and its two parity channels.  The
least-used sidecars are the selected-distance packet and source-depth parity.

All algebra is over `F_2`.  Retain THM-3471's notation

```text
alpha_u(d)=F_(u+d)(d),
K(n,r)=[x^r](1+x+x^2)^n,
t=u+1+2d+v.                                          (1)
```

The complete marked inward tail, with the boundary circuit removed, is

```text
R^in(z,q)
 =sum_(u>=3) R_u(z,q)
 =sum_(u>=3,d,v>=0) alpha_u(d)K(d+v,v)
       z^(u+1+2d+v)q^v.                              (2)
```

Every target coefficient

```text
C_t(q)=[z^t]R^in(z,q)                                (3)
```

is a polynomial: only finitely many triples in (2) reach a fixed `t`.

## 2. The unique ballistic Newton face

Give a monomial its ballistic weight

```text
omega(z^t q^v)=t-v.                                  (4)
```

On a physical event (1),

```text
omega=u+1+2d.                                        (5)
```

Thus changing Green slack moves along one Newton face, while source depth and
distance choose the face.

### Theorem 2.1 (lowest face and the first gap)

The inward tail has

```text
boxed:
 Face_4 R^in = z^4/(1+zq),
 Face_5 R^in = 0,
 Face_6 R^in = z^6/(1+zq).                           (6)
```

#### Proof

Weight four forces `(u,d)=(3,0)`.  Directly from the third Rule 30 row,

```text
a_3(-3),...,a_3(3)=1101111,
alpha_3(0)=a_3(0)a_3(1)=1.                           (7)
```

The binary Green recursion gives

```text
K(2n,2n)=K(n,n),
K(2n+1,2n+1)=K(n,n),                                (8)
```

so `K(v,v)=1` for every `v>=0`.  The complete face is therefore

```text
sum_(v>=0) z^(v+4)q^v=z^4/(1+zq).                   (9)
```

Weight five has the only candidate `(u,d)=(4,0)`, and
`alpha_4(0)=0`.  Weight six has candidates `(5,0)` and `(3,1)`; the exact
edge words give `alpha_5(0)=1` and `alpha_3(1)=0`.  The same central Green
ray proves the last identity in (6).  This proves the theorem.  `square`

More generally, with `x=zq` and

```text
D_d(x)=sum_(v>=0)K(d+v,v)x^v,                       (10)
```

the weight-`n` face is the finite expression

```text
Face_n R^in
 =z^n sum_(u+1+2d=n) alpha_u(d)D_d(x).              (11)
```

The face filtration is therefore exact and nonrectangular, not an
approximation by finitely many source depths.

## 3. Full-tail monicity and the global jet staircase

### Theorem 3.1 (monic slack polynomial)

For every `t>=4`,

```text
boxed: deg_q C_t=t-4 and [q^(t-4)]C_t=1.             (12)
```

For `t>=6`, its top three possible coefficients are

```text
([q^(t-4)]C_t,[q^(t-5)]C_t,[q^(t-6)]C_t)=(1,0,1).  (13)
```

#### Proof

At target `t`, every event has

```text
v=t-u-1-2d <= t-4.                                  (14)
```

Equality uniquely forces `(u,d,v)=(3,0,t-4)`.  Equations (7)--(8) give its
coefficient one.  Equations (13) are the target slices of the next two faces
in (6).  `square`

Put `q=1+epsilon` and define the complete inward Hasse layers

```text
J_j(z)=[epsilon^j]R^in(z,1+epsilon).                 (15)
```

Since a polynomial of degree below `j` has zero `j`th Hasse coefficient,
(12) immediately gives the exact causal staircase

```text
boxed:
 ord_z J_j=j+4,
 [z^(j+4)]J_j=1                    for every j>=0.   (16)
```

Thus no Hasse layer is emptied by the unbounded source-depth sum.  The
detecting order in (16) moves with the target; this is not a bounded-jet
theorem.

## 4. The full-sum parity-Cartier atlas

THM-3476's ramification-two law was stripwise.  It becomes an exact
two-channel law after summing all inward strips.

Fix `t`.  For each `3<=u<t`, put

```text
R_u=t-u-1=2L_u+delta_u,             delta_u in {0,1},
beta_(u,t,d)=alpha_u(d)K(R_u-d,R_u),
                                      0<=d<=L_u.     (17)
```

Green palindromy turns `beta` into the coefficient of the event with slack
`R_u-2d`.  Aggregate the reversed selected-distance packets by their parity:

```text
H_(delta,t)(X)
 =sum_(3<=u<t, delta_u=delta) sum_(0<=d<=L_u)
      beta_(u,t,d) X^(L_u-d).                        (18)
```

Because `delta_u=t-u-1 mod 2`, these are equivalently the two source-depth
parity channels at fixed `t`.

### Theorem 4.1 (Cartier inversion and valuation law)

For every target,

```text
boxed: C_t(q)=H_(0,t)(q^2)+q H_(1,t)(q^2).           (19)
```

Write

```text
h_(delta,j)=[eta^j]H_(delta,t)(1+eta),
J_i=[epsilon^i]C_t(1+epsilon).                       (20)
```

Then, for every `j>=0`,

```text
boxed:
 J_(2j)=h_(0,j)+h_(1,j),
 J_(2j+1)=h_(1,j),                                  (21)

 h_(1,j)=J_(2j+1),
 h_(0,j)=J_(2j)+J_(2j+1).                           (22)
```

Let

```text
r_delta=ord_(X+1) H_(delta,t) in N union {infinity}, (23)
```

where the zero polynomial has order infinity.  Then

```text
ord_(q+1) C_t =
  infinity,                       r_0=r_1=infinity;
  2 min(r_0,r_1),                 r_0 != r_1;
  2r+1,                           r_0=r_1=r<infinity. (24)
```

#### Proof

An event in (18) has slack

```text
R_u-2d=delta_u+2(L_u-d),                            (25)
```

which proves (19).  In characteristic two,
`(1+epsilon)^2=1+epsilon^2`; substituting this in (19) gives (21), and (22)
is immediate.  If the two channel orders differ, the earlier even coordinate
in (21) is one.  If both first occur at `r`, their even coordinates cancel
and the odd coordinate is one.  This proves (24).  `square`

Equation (24) is the exact full-sum repair of the fixed-strip doubling law:
ramification still doubles unequal channel orders, while a tie creates one
additional odd layer.  No individual source owner survives the aggregation
in (18).

### 4.1 Two minimal hostile controls

At target six only the even-slack channel survives:

```text
C_6(q)=1+q^2=(1+q)^2.                               (26)
```

This realizes unequal-channel doubling.  The earliest cross-parity scalar
cancellation is at target nine:

```text
H_(0,9)=1,
H_(1,9)=1+X+X^2,
C_9(q)=1+q+q^3+q^5.                                 (27)
```

Both channel orders are zero, so `C_9(1)=0` but its first Hasse coefficient
is one.  Monicity therefore does not survive scalarization even at a small
physical target.

## 5. A shallow atlas for every source-transverse coordinate

The preceding results filter after Green transport.  There is also a sharp
source-local statement before target extraction.  Put

```text
T(z,X)=sum_(u>=3,d>=0) alpha_u(d)z^uX^d.             (28)
```

Let `W=z^2+O(z^3)` be THM-3476's small root of

```text
P(z,X)=X^2+(1+z)X+z^2.                              (29)
```

With `Y=X+W`, one has `P=Y(Y+1+z)`, and the second factor is a unit.  Thus
the `Y`-adic Hasse coordinates are a fixed coordinate chart compatible with
the `P`-adic filtration.  Define

```text
G_k(z)=[Y^k]T(z,W+Y).                                (30)
```

Hasse expansion gives

```text
G_k(z)=sum_(u>=3,r>=0)
 alpha_u(k+r) binom(k+r,k) z^u W^r.                 (31)
```

Lucas--Kummer says `binom(k+r,k)=1 mod 2` exactly when the bitwise AND
`k AND r` is zero; moreover `ord_z W=2`.  Hence a surviving pair `(u,r)`
starts at `z^(u+2r)`.

### Theorem 5.1 (exact 64-class shallow-owner atlas)

One has

```text
ord_z G_0=3,             unique leading owner (u,r)=(3,0). (32)
```

For `k>=1`, the following disjoint table is exact.  The source distance of
the owner is `d=k+r`.

| condition on `k` | `ord_z G_k` | unique `(u,r)` |
|---|---:|---:|
| `k=2 mod 4` | 3 | `(3,0)` |
| `k=1,5,7 mod 8` | 4 | `(4,0)` |
| `k=4 mod 8` | 5 | `(5,0)` |
| `k=0 mod 8` | 6 | `(4,1)` |
| `k mod 64` in `{3,35}` | 7 | `(7,0)` |
| `k mod 64` in `{19,51}` | 8 | `(8,0)` |
| `k mod 64` in `{43,59}` | 9 | `(9,0)` |
| `k mod 64=11` | 11 | `(11,0)` |
| `k mod 64=27` | 12 | `(4,4)` |               (33)

Consequently

```text
boxed:
 G_k !=0 for every k>=0,
 3<=ord_z G_k<=12,
 sup_k ord_z G_k=12.                                 (34)
```

In particular,

```text
ord_P T=0.                                           (35)
```

Since THM-3476's boundary source `S_[0,2]` is a multiple of `P`, the complete
physical source `S_[0,2]+T` also has `P`-adic order zero.  The boundary
`P`-factor does not persist one step into the physical inward completion.

#### Proof of the atlas

The closed edge recurrences from THM-3476 give exact state returns, not fitted
prefixes.  For `u=3,...,12`, the preperiod/period pairs of `alpha_u` are

```text
(1,4),(0,8),(1,8),(0,16),(0,32),
(1,32),(0,64),(0,64),(0,64),(0,64).                 (36)
```

The right prefix through `b_u` and left prefix through `ell_(u+1)` are closed
under the triangular recurrences; equality of the complete prefix state at
the endpoints in (36) proves the words for all `d`.

Every proposed leading order in (33) is at most twelve.  A competing term in
(31) must therefore have `u<=12` and `r<=4`; every `u>=13` or `r>=5` starts
later.  The words (36), the no-carry condition `k AND r=0`, and the 64 residue
classes leave exactly the unique winners displayed in (33).  Each `W^r` has
leading coefficient one, so the unique winner cannot cancel.  Residue 27
attains order twelve infinitely often, proving sharpness.  The finite table
is checked from complete recurrent states by the companion.  `square`

The theorem controls only the first `z`-coefficient of each coordinate
`G_k`; deeper source depths still contribute to every full coordinate.

## 6. Preservation and loss ledger

| map | preserves | destroys | required sidecar / boundary |
|---|---|---|---|
| event -> ballistic face | exact `t-v=u+1+2d` and the face current | individual ternary path ancestry | slack `v` and source intercept |
| inward tail -> `C_t(q)` | exact target and slack polynomial; monic top owner | individual `u,d` below the top face | full marker or selected-distance packet |
| `q -> 1` | scalar inward contribution | top slack, face, and parity owner | Hasse tower; (26)--(27) are hostile |
| two Cartier channels -> scalar | XOR of the two channel values | separate `u` parity and repeated-root order | first jet for values; paired tower for all coordinates |
| `T(z,X) -> G_k(z)` | exact translated-root source coordinate | decomposition of nonleading coefficients by depth | `u,d` if owner identity is needed |
| source completion modulo `P` | exact unmarked Green evaluation | every higher transverse coordinate | full `Y`- or slack-jet tower |

The unique face is an endpoint-grade parity, not a unique ternary Green path:
`K(v,v)=1` can summarize many paths before reduction modulo two.  No
tournament or Berggren ancestry statement follows.

## 7. Exact companion and validity boundary

The companion is intended to run identically with ordinary and optimized
Python:

```text
python3 04-computation/rule30_inward_slack_cartier_thm3488.py
python3 -O 04-computation/rule30_inward_slack_cartier_thm3488.py
```

It uses explicit failure gates and checks:

1. centered Rule 30 rows against independent closed edge recurrences;
2. complete state return for every source cycle in (36);
3. polynomial and binary-digit Green engines, including `K(v,v)=1`;
4. native-event target polynomials, faces (6), monicity, and the global jet
   staircase on declared finite universes;
5. the full Cartier atlas and valuation law, including (26)--(27);
6. all 64 source-Taylor competitor classes and a direct truncated
   characteristic-root substitution.

The inequalities (5), (14), and (31), the Green recursion (8), Frobenius in
(21), and the exact recurrent state returns are the universal proofs.  The
finite universes audit the implementation rather than replace those proofs.

## 8. No-prize scope and next exact targets

Monicity proves `C_t(q)` is nonzero as a marked polynomial for every `t>=4`.
It does not prove `C_t(1)` is nonzero; (26)--(27) already show the opposite.
The staircase uses jet order `t-4`, so it supplies no bounded carrier.  The
Cartier atlas retains source-depth parity but loses individual depth, source
time, and THM-3481's terminal phase owner.  The shallow Taylor table controls
only initial `z`-terms, not long-time coefficients or computational cost.

Accordingly this theorem proves no Rule 30 nonperiodicity, balance, or
complexity prize.  The next exact targets are:

1. combine the two Cartier channels with the calibrated terminal phase;
2. classify cancellation between ballistic faces after `q=1`; and
3. decide whether the full-tail first-live Hasse order is unbounded, rather
   than inferring it from one fixed strip.
