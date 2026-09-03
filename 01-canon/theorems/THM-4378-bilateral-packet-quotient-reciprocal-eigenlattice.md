---
id: THM-4378
title: "Bilateral packet quotient and reciprocal eigenlattice"
status: >
  PROVED ELEMENTARY RELATIVE TO THM-4368, THM-4369, THM-4375, AND THM-4377 +
  VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT; JC(2) OPEN. For every finite
  bilateral packet prefix at ell>=3, the complete integral trace quotient is
  a full truncated polynomial
  lattice carried by a width-one central spine. The parity-signed reciprocal
  involution descends and has explicit saturated invariant and anti-invariant
  lattices. Their sum has an exact elementary two-primary Smith defect; the
  unsigned coordinate swap does not descend. The ell=2 diagonal boundary is
  classified separately, and the gluing class of the natural fibre-weighted
  aggregate is exactly the parity of THM-4377's source-fibre imbalance. No
  bracket, seam, JC(2), or DC(2) consequence is asserted.
source: root + reciprocal_kernel / JC2 continuation session, 2026-09-03
depends_on:
  - THM-4368-diagonal-boundary-valuation-triangular-address-and-simplex-stream-rank
  - THM-4369-source-packet-pascal-circuit-kernel-and-boundary-basis
  - THM-4375-bilateral-source-cone-reciprocal-orbits-and-fibre-imbalance
  - THM-4377-reciprocal-source-mutation-and-boundary-jet-cokernel
related:
  - THM-4364-sharp-binomial-diagonal-annihilator-hierarchy
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
mistake_firewall:
  - MISTAKE-222
  - MISTAKE-434
primary_script: 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_thm4378.py
primary_output: 05-knowledge/results/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_thm4378.out
primary_script_sha256: b93dec50de6e75b571c3f1b1cdc1649472504c19b62dfb0c09baf88c7ddaa51b
primary_output_sha256: c23b9afc7682476e1ce0f1e8d9fb1e5e9939395821811871df971386a8d58b6a
independent_referee_script: 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_independent_referee_thm4378.py
independent_referee_output: 05-knowledge/results/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_independent_referee_thm4378.out
independent_referee_script_sha256: 70801c2c943c48258dc6b57781edc4f33248c36b09000662eb3387d93346c31e
independent_referee_output_sha256: 60127767ac254b6069d45e80acdb0466c6927dfd9074a88619f6dc5aca0ac73f
hash_basis: raw LF bytes
audit: >
  PRIMARY AND CHECKED-IN INDEPENDENT EXACT AUDITS PASS + INDEPENDENT
  SCRATCH-AUDIT SIGNAL. The primary rebuilds
  the bilateral band directly, verifies its width-one integral spine,
  reconstructs every packet in that basis, checks the signed reciprocal
  intertwiner, the eigenlattice Smith blocks, norm/difference images, ell=2
  boundary, smallest unsigned-reflection hostile, and the source-fibre parity
  bridge. The recurrence referee independently checks the same bridge.
  Independent broadcast
  audit c99e9177b5 separately found R_0 C=2e-C^R, the checkerboard repair,
  elementary two-gluing, and no hidden four-torsion. Normal, optimized, and
  fixed-hash-seed primary and referee replays agree with their frozen raw-LF
  outputs.
---

# THM-4378 -- Bilateral packet quotient and reciprocal eigenlattice

**PROVED ELEMENTARY RELATIVE TO THM-4368, THM-4369, THM-4375, AND THM-4377 +
VERIFIED-EXACT + INDEPENDENTLY VERIFIED-EXACT. THIS IS AN INTERNAL FINITE
PACKET-QUOTIENT AND SOURCE-PARITY THEOREM. IT DOES NOT ASSERT A BRACKET
OBSTRUCTION, CHART ENTRY, SEAM CONTROL, `JC(2)`, OR `DC(2)`.**

## 1. Setup and theorem

Fix `ell>=2` and put

```text
s=ceil(ell/2),             rho=ceil(ell/3),
q=z(1-z),                  x=2z-1,
tau0=2rho-1.                                                    (1)
```

THM-4368 associates to a positive boundary pair `(u,v)` the normalized
complete diagonal trace

```text
F_(u,v)(z)=(-1)^(u-1) z^(u-1)(1-z)^(v-1).                     (2)
```

For `R>=tau0`, let

```text
B_(ell,R)={(u,v) in B_ell : u+v-1<=R},
D=R-tau0,                                                        (3)
```

where `B_ell` is THM-4375's bilateral source cone. Let `M_(ell,R)` be
the free abelian group on symbols `e_(u,v)` for `(u,v) in B_(ell,R)`, and
define

```text
Phi:M_(ell,R) -> Z[z],       Phi(e_(u,v))=F_(u,v),
K_(ell,R)=ker(Phi),          Q_(ell,R)=M_(ell,R)/K_(ell,R).     (4)
```

By THM-4369, `K_(ell,R)` is equivalently the complete expanded-packet or
complete fixed-order stream kernel after restriction to this bilateral
packet module. The theorem is as follows.

### Theorem A (a width-one spine is the full integral quotient)

If `ell>=3`, then

```text
Q_(ell,R) = q^(rho-1) Z[z]_(<=D)                              (5)
```

as an integral trace lattice, not merely after tensoring with `Q`. An
explicit packet basis is, up to the harmless signs in `(2)`,

```text
q^(rho-1) q^k       from (u,v)=(rho+k,rho+k),       2k<=D,
q^(rho-1) zq^k      from (u,v)=(rho+k+1,rho+k),     2k+1<=D.   (6)
```

Thus the loops and one chosen orientation of the nearest-neighbour arcs in
THM-4375's path-power band already carry the entire packet quotient. Wider
arcs add packet symbols and relations, but no new quotient direction. In
particular,

```text
rank Q_(ell,R)=D+1,
rank K_(ell,R)=|B_(ell,R)|-(D+1),                              (7)
```

and `K_(ell,R)` is saturated in `M_(ell,R)`; the quotient `Q_(ell,R)` is
torsion-free.

### Theorem B (only the signed reciprocal action descends)

The formula

```text
iota(e_(u,v))=(-1)^(u-v)e_(v,u)                               (8)
```

defines an involution of `M_(ell,R)`, preserves `K_(ell,R)`, and descends on
`Q_(ell,R)` to

```text
J(f)(z)=f(1-z).                                                (9)
```

The unsigned coordinate swap `e_(u,v) -> e_(v,u)` does **not** in general
preserve the packet kernel, even though the bilateral type set itself is
swap-stable.

### Theorem C (integral eigenlattices and exact two-glue)

For `ell>=3`, the invariant and anti-invariant lattices are

```text
Q^+ = q^(rho-1) Z[q]_(<=floor(D/2)),                           (10)

Q^- = q^(rho-1) x Z[q]_(<=floor((D-1)/2)).                    (11)
```

When `D=0`, `(11)` denotes the zero lattice. Each of `Q^+` and `Q^-` is
individually saturated in `Q`. Their direct sum need not be saturated:

```text
Q/(Q^+ direct_sum Q^-) ~= (Z/2Z)^(ceil(D/2)),                 (12)

Sat_Q(Q^+ direct_sum Q^-)=Q.                                  (13)
```

Equivalently, reciprocal eigensplitting is integral after inverting `2`, and
the complete obstruction before inversion is the elementary Smith defect in
`(12)`; there is no hidden odd-primary defect.

### Theorem D (difference, norm, and Tate quotients)

For `ell>=3`,

```text
(1-J)Q=Q^-,                                                   (14)

Q^+/(1+J)Q ~= Z/2Z,     D even,
             0,          D odd.                              (15)
```

Thus the anti-invariant Tate quotient is zero and the invariant Tate quotient
is the single even-top class in `(15)`. The ordinary quotients are

```text
Q/(1-J)Q ~= Z^(floor(D/2)+1),                                 (16)

Q/(1+J)Q ~= Z^(ceil(D/2)) direct_sum
             { Z/2Z,  D even;
               0,     D odd. }                               (17)
```

### Theorem E (the `ell=2` boundary)

At `ell=2`, the bilateral cone contains only `(u,v)=(w,w)`. Hence

```text
Q_(2,R)=Z[q]_(<=floor(D/2)),          J=id,                   (18)
K_(2,R)=0,
Q^+=Q,                               Q^-=0.                   (19)
```

There is no eigenlattice gluing defect. Instead,

```text
Q/(1-J)Q=Q,
Q/(1+J)Q ~= (Z/2Z)^(floor(D/2)+1).                            (20)
```

The inheritance pass was:

- closest proved mechanism: THM-4369's integral Pascal-circuit kernel and
  boundary basis;
- canonical hostile: THM-4375's first nonfixed bilateral fibre pair;
- corrected near miss: swap-stability of packet **types** does not imply that
  unsigned swap preserves signed packet relations;
- least-used sidecar: the integral, rather than rational, eigensplitting of
  the complete packet quotient.

The live board was

```text
bilateral band | Pascal kernel | central spine | signed reciprocity
integral eigenlattice | Smith two-glue | norm/difference quotient.       (21)
```

## 2. Proof of the width-one quotient

By THM-4375,

```text
B_ell={(u,v):u,v>=rho, |u-v|<=s-1}.                          (22)
```

Write

```text
u=rho+i,                    v=rho+j.                         (23)
```

Then `(3)` and `(22)` become

```text
i,j>=0,                     |i-j|<=s-1,
i+j<=D.                                                         (24)
```

With `r=rho-1`, equation `(2)` factors as

```text
F_(u,v)=(-1)^(rho+i-1) q^r z^i(1-z)^j.                      (25)
```

Every trace in the bilateral prefix therefore belongs to the right-hand side
of `(5)`.

For `ell>=3`, one has `s-1>=1`, so all the width-zero and width-one pairs

```text
(i,j)=(k,k),                  2k<=D,
(i,j)=(k+1,k),                2k+1<=D                       (26)
```

belong to `(24)`. Their residual polynomials are exactly

```text
q^k,                          zq^k.                          (27)
```

Ordered by residual degrees `0,1,...,D`, the list

```text
1,z,q,zq,q^2,zq^2,...                                         (28)
```

has coefficient matrix triangular with diagonal

```text
1,1,-1,-1,1,1,...
```

in the monomial basis `1,z,...,z^D`: `q^k` and `zq^k` have respective
degrees `2k` and `2k+1`, both with leading coefficient `(-1)^k`. Hence `(28)`
is a unimodular basis of `Z[z]_(<=D)`. This proves `(5)` and `(6)`.
Equation `(7)` follows from the first isomorphism theorem. Since the image is
free, the kernel is saturated; `(5)` itself shows that the quotient is free.

The proof also explains the width collapse. It is a quotient statement, not
a claim that wider packet types or their source fibres disappear before
passing through `Phi`.

## 3. The signed involution and the first hostile

The bilateral band and every finite block prefix are stable under coordinate
swap. Formula `(8)` squares to the identity. Directly from `(2)`,

```text
Phi(iota(e_(u,v)))
 =(-1)^(u-v)F_(v,u)(z)
 =F_(u,v)(1-z)
 =J(Phi(e_(u,v))).                                           (29)
```

Thus `Phi iota=J Phi`; in particular, `iota` preserves `K` and induces `(9)`.
The parity sign is forced by the trace normalization and is not cosmetic.

The smallest hostile occurs at

```text
ell=3,       (s,rho)=(2,1),       R=2,       D=1.           (30)
```

All three types `(1,1),(1,2),(2,1)` are bilateral. THM-4369's first Pascal
circuit is

```text
C=e_(1,1)-e_(1,2)+e_(2,1).                                  (31)
```

Its trace is

```text
1-(1-z)-z=0.                                                  (32)
```

Unsigned reflection sends it to

```text
R_0(C)=e_(1,1)-e_(2,1)+e_(1,2),
Phi(R_0(C))=1-(-z)+(1-z)=2.                                  (33)
```

So unsigned swap does not descend. In contrast, `(8)` gives

```text
iota(e_(1,2))=-e_(2,1),       iota(e_(2,1))=-e_(1,2),
iota(C)=C.                                                     (34)
```

There is no smaller hostile: `ell=2` has only fixed pairs, and the `R=1`
prefix at `ell=3` also has only the fixed type `(1,1)`.

## 4. Integral reciprocal eigenlattices

Factor out the fixed nonzerodivisor `q^(rho-1)`. The unimodular basis `(28)`
gives every `f in Z[z]_(<=D)` a unique expression

```text
f=A(q)+zB(q),                                                 (35)
```

where

```text
deg_q A<=floor(D/2),
deg_q B<=floor((D-1)/2).                                      (36)
```

Since `J(q)=q` and `J(z)=1-z`,

```text
J(A+zB)=(A+B)-zB.                                             (37)
```

Equation `J(f)=f` forces `2B=0`, hence `B=0` in the torsion-free polynomial
ring. Equation `J(f)=-f` forces `B=-2A`, hence

```text
f=(1-2z)A=-xA.                                                (38)
```

This proves `(10)` and `(11)`. Each is the kernel of one of the endomorphisms
`J-1` and `J+1` of the torsion-free lattice `Q`, so each is saturated.

Now order `(28)` in two-element blocks `(q^k,zq^k)`. The corresponding
eigenvectors are

```text
q^k,                         xq^k=-q^k+2zq^k.                (39)
```

Their inclusion matrix in each paired block is

```text
[ 1 -1 ]
[ 0  2 ].                                                     (40)
```

Adding the first column to the second gives the Smith block `diag(1,2)`.
If `D` is even, the unpaired terminal vector `q^(D/2)` contributes one unit
Smith block. There are exactly `ceil(D/2)` paired blocks. This proves
`(12)`--`(13)`. Concretely, the quotient map in `(12)` records each
coefficient of `zq^k` modulo two. Its kernel is the eigenlattice sum because

```text
2zq^k=q^k+xq^k.                                               (41)
```

At the hostile `(30)`, `Q=Z{1,z}`, while

```text
Q^+ direct_sum Q^-=Z{1,2z-1}.                                (42)
```

The determinant is two, so this is also the first eigen-gluing defect.

## 5. Difference and norm

On each paired block, equations `(37)` and `(39)` give

```text
(1-J)(zq^k)=xq^k,
(1+J)(zq^k)=q^k,
(1+J)(q^k)=2q^k.                                             (43)
```

The first identity proves `(14)`. If `D` is odd, every invariant basis vector
`q^k` is paired with `zq^k` and is therefore a norm. If `D=2m` is even, this
is true for `k<m`, while the unpaired top vector `q^m` has norm lattice
`2Zq^m`. This proves `(15)`.

Quotienting the central basis blockwise by `(43)` gives `(16)` and `(17)`.
In cohomological language,

```text
Q^-/((1-J)Q)=0,
Q^+/((1+J)Q)=1_[D even] Z/2Z.                               (44)
```

This two-torsion is distinct from the `ceil(D/2)` eigen-sum gluing classes in
`(12)`: `(12)` measures failure of integral eigensplitting, whereas `(44)`
measures failure of the norm to hit the even top invariant.

## 6. The `ell=2` endpoint

For `ell=2`, `(s,rho)=(1,1)`, and `(22)` forces `u=v=w`. Put `w=1+k`.
The row bound is `2k<=D`, and

```text
F_(w,w)=(-1)^kq^k.                                           (45)
```

These distinct-degree polynomials are integrally independent, proving
`(18)`--`(19)`. The reciprocal action is the identity, so its difference is
zero and its norm is multiplication by two. This proves `(20)`. Notice that
the exceptional endpoint has norm two-torsion but no failure of integral
eigensplitting, because there is no anti-invariant direction at all.

## 7. The source-fibre parity bridge

There is one exact connection to the source sidecar isolated by THM-4377.
Fix `d>0`, a bilateral orbit

```text
e_+=e_(w+d,w),              e_-=e_(w,w+d),                 (46)
```

and a prefix containing its block. Let `mu_+,mu_-` be the two source-fibre
sizes from THM-4375/4377, let `bar(e)_+,bar(e)_-` be the classes of these
packet symbols in `Q`, and let

```text
gamma:Q -> Q/(Q^+ direct_sum Q^-)
```

be the two-glue quotient from `(12)`. Then

```text
gamma(mu_+ bar(e)_+ + mu_- bar(e)_-) !=0
        iff mu_+-mu_- is odd.                                (47)
```

Indeed, up to unit signs,

```text
Phi(e_+)=q^(w-1)z^d,
Phi(e_-)=q^(w-1)(1-z)^d.                                     (48)
```

Write uniquely

```text
z^d=A_d(q)+zB_d(q).                                          (49)
```

Modulo `q`, the relation `z^2=z-q` gives `z^d=z` for every `d>0`, so
`B_d(0)=1`. Therefore `e_+` has a nonzero gluing vector, and reciprocal
substitution changes `B_d` to `-B_d`, which is the same vector modulo two.
Under the identification `(5)`, the weighted class has gluing vector

```text
(mu_++mu_-) gamma(bar(e)_+)
 =(mu_+-mu_-) gamma(bar(e)_+),                               (50)
```

proving `(47)`.

In THM-4377's postclock range `w>=s`,

```text
mu_+-mu_-=beta=min(w-s+d+1,2d).                              (51)
```

Thus the gluing class toggles whenever the growing boundary-jet rank grows by
one, and vanishes permanently at the full-jet threshold `w>=s+d-1`, where
`beta=2d`. This is an exact coupling between source-fibre imbalance and the
packet quotient. It does not say that a bracket or seam functional detects
`gamma`.

## 8. Consequence, audit, and scope

Any additive consumer on a finite bilateral prefix which factors through the
complete THM-4368 packet trace is determined integrally by the width-one
spine `(6)`. Consequently a new obstruction cannot arise merely from adding
longer reciprocal arcs to that same trace quotient. It must either constrain
the spine coefficients or read information discarded by `Phi`, such as the
individual source-monomial fibres, which THM-4375 shows can be imbalanced, or
a bracket/seam/entry sidecar.

The primary rebuilds the bilateral band directly, verifies the width-one
integral spine, reconstructs every packet in that basis, checks the signed
reciprocal intertwiner, the eigenlattice Smith blocks, norm/difference images,
the `ell=2` boundary, and the smallest unsigned-reflection hostile. Its finite
universe is

```text
ell=2..80,             quotient excess D=0..30,
reciprocal/Smith blocks D=0..60.                              (52)
```

These checks verify consequences of the proved formulas, not an open claim.
The checked-in independent referee imports no primary code and instead works
in the quadratic free module

```text
Z[q] direct_sum z Z[q],                 z^2=z-q.              (53)
```

Its two multiplication recurrences rebuild every band packet for `ell=2..60`
and `D=0..24`; separate coefficient boxes through `D=32` recover the
eigenconditions, Smith parity map, difference, norm, and both hostiles.
An independent scratch audit recorded in broadcast commit `c99e9177b5`
separately obtained, for every local Pascal cell,

```text
R_0 C_(u,v)=2e_(v,u)-C_(v,u),
iota C_(u,v)=(-1)^(u-v) C_(v,u),                            (54)
```

and independently confirmed that the gluing has exponent two with no hidden
four-torsion. Its checkerboard sign `(-1)^(i+j)` in residual coordinates is
the same as `(-1)^(u-v)` because these exponents have the same parity. That
audit is corroborating evidence rather than a checked-in proof dependency.

Reproduce from the repository root:

```text
python3 -B 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_thm4378.py
python3 -B -O 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_thm4378.py
python3 -B 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_independent_referee_thm4378.py
python3 -B -O 04-computation/jc2_bilateral_packet_quotient_reciprocal_eigenlattice_independent_referee_thm4378.py
```

The theorem concerns finite linear combinations of THM-4368 bilateral packet
types. It does not classify arbitrary source polynomials, retain individual
source exponents, construct a bracket obstruction, prove that a hypothetical
counterexample enters this chart, cross the seam, prove termination, produce
a Keller pair, or prove `JC(2)` or `DC(2)`.
