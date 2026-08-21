---
id: THM-3574
title: "Universal reducible target-graph component unit no-go"
status: >
  PROVED + VERIFIED-EXACT.  Every irreducible component of every reducible
  nonzero polynomial target-graph pullback for the fixed THM-1300 map has a
  nonconstant unit and hence is not A2.  THM-3573 writes the graph as
  phi_H=4H(1+bH+4aH^2).  Its source pullback factors universally as S_H Q_H.
  On every factor of S_H, R(F1,F2) has inverse 1+xy; on every factor of Q_H,
  M(F1,F2) has an exact polynomial inverse.  Both R and M are nonconstant,
  and quasi-finiteness makes every component dominate the target graph.
source: kps-s188
depends_on:
  - THM-3573-polynomial-target-graph-pell-parameter-descent-classification
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
related:
  - THM-3568-reducible-target-graph-component-euler-no-go
  - THM-3571-quadratic-target-graph-euler-no-go
companion: 04-computation/jacobian_universal_reducible_graph_component_unit_no_go_kps_s188.py
output: 05-knowledge/results/jacobian_universal_reducible_graph_component_unit_no_go_kps_s188.out
script_sha256: a0f53bdc4576bf586a982f291590b75e47bab07f37530d6055431f783f7ef333
output_sha256: 7188e9f7b4a6e9e781abd2cd72e528dee6fc43c09cb18cc19758dd5fbb758571
hash_basis: LF-normalized bytes
---

# THM-3574 -- universal reducible target-graph component unit no-go

**PROVED + VERIFIED-EXACT.**  The proper-factor route through triangular
target graphs is closed in every degree.  Every component of every reducible
nonzero graph pullback carries an explicit nonconstant unit, so no component
can be a source coordinate plane.

All rings and varieties are over `C`.  Use the fixed THM-1300 source map
`F=(F1,F2,F3)` and put

```text
u=1+xy.                                                  (1)
```

## 1. Statement

Let `0!=phi in C[a,b]`, and suppose

```text
P_phi=F3+phi(F1,F2)                                    (2)
```

is reducible.  For every irreducible factor `f` of `P_phi`, the coordinate
ring

```text
C[x,y,z]/(f)                                            (3)
```

has a nonconstant unit.  Consequently `V(f)` is not isomorphic to `A2`, and
`f` is not a coordinate polynomial.  Thus `(2)` has no source-coordinate
factor.

The exclusion `phi!=0` is necessary: `F3=x(2-3xy-x^2z)`, whose factor `x` is
a source coordinate.

## 2. The universal source factorization

By THM-3573, reducibility of `(2)` is equivalent to the existence of
`0!=H in C[a,b]` such that

```text
phi=phi_H=4H U,

U=1+bH+4aH^2.                                          (4)
```

Introduce three more target polynomials

```text
R=1+2bH+12aH^2,

M=48a^2H^2+8abH+16a-b^2,

V=1-bH+12aH^2.                                         (5)
```

First treat `H` as an independent indeterminate.  Direct expansion in
`C[x,y,z,H]`, followed by the lawful substitution `H=H(F1,F2)`, gives

```text
F3+phi_H(F1,F2)=S_H Q_H,

S_H=x+2u H(F1,F2).                                     (6)
```

Before the substitution, the second factor is the compact polynomial

```text
Q_H =
 8H^2x^2y^2z+24H^2xy^3+16H^2xyz+32H^2y^2+8H^2z
 +2Hx^2yz+6Hxy^2+2Hxz+2Hy-x^2z-3xy+2.                (7)
```

Equation `(7)` is then evaluated at `H(F1,F2)`.  This is a polynomial
factorization; no localization or generic specialization is used.

## 3. The unit on the linear packet

The first exact Bezout identity is

```text
R(F1,F2)u-1=S_H A_H,                                   (8)
```

where, before substituting `H(F1,F2)`,

```text
A_H=
 6Hx^3y^3z+18Hx^2y^4+18Hx^2y^2z+42Hxy^3
 +18Hxyz+24Hy^2+6Hz+y.                                 (9)
```

Hence on every irreducible component whose equation divides `S_H`, the
target function `R(F1,F2)` is a unit with inverse `u`.

It is not constant on such a component.  If `m=deg_a(H)`, the unique top
`a`-term in `(5)` gives

```text
deg_a(R)=2m+1.                                         (10)
```

Moreover, THM-2473 makes `F` quasi-finite.  Every hypersurface component of
`(2)` therefore maps dominantly to the graph `T_phi~=A2`: its two-dimensional
image closure cannot be a proper closed subset of the irreducible target
surface.  Dominance injects `C[a,b]` into the component function ring, so the
nonconstant target polynomial `R` remains nonconstant there.

## 4. The unit on the quadratic packet

The short target identity driving the second component is

```text
H^2M+1=UV.                                              (11)
```

THM-3573's core packet is

```text
E_phi(X)=(RX+2H)(RMX^2-2HMX+4U).                       (12)
```

On the source, the second factor in `(12)` is exactly divisible by `Q_H`:

```text
M(F)(R(F)x^2-2H(F)x)+4U(F)=Q_H B_H                    (13)
```

for a polynomial `B_H in C[x,y,z]`.  Here and below `F` in parentheses means
`(F1,F2)`.  The exact quotient in `(13)` is not printed because its expanded
form has 152 terms; polynomial division over `Q` is an active companion gate.

Define the polynomial

```text
K_H=-H(F)^2
    -V(F)(R(F)x^2-2H(F)x)/4.                           (14)
```

Combining `(11)` and `(13)` gives a literal Bezout identity

```text
1-M(F)K_H=Q_H C_H                                     (15)
```

with `C_H in C[x,y,z]`.  Thus `M(F1,F2)` is a unit on every irreducible
component whose equation divides `Q_H`.

It is again nonconstant.  The unique top `a`-term gives

```text
deg_a(M)=2m+2,                                         (16)
```

and the same quasi-finite dominance argument applies.  Therefore no
irreducible component from `Q_H`, including any further splitting of `Q_H`,
is `A2`.

Every irreducible factor of `P_phi=S_HQ_H` divides at least one of the two
packets.  Equations `(8)` and `(15)` therefore supply a nonconstant unit on
every component, proving the theorem.

## 5. Geometry and boundary checks

The proof does not require either packet to be irreducible.  It also avoids
Euler-characteristic cancellation: the obstruction lives separately on
each component and survives every further factorization.

For `H=-h(b)/2`, equations `(5)--(6)` become exactly THM-3565/3568's first
resonance packet:

```text
S_H=x-u h(F2),

R=3ah^2-bh+1,

M=12a^2h^2-4abh+16a-b^2.                              (17)
```

Thus THM-3568's degree-one Euler obstruction is the first slice of the
universal unit theorem.  The present result is stronger in degree range but
weaker in topology: it does not compute the Euler characteristic of the
components.

At the excluded boundary `H=0`, `(6)` reduces to

```text
F3=x(2-3xy-x^2z),                 R=1.                 (18)
```

The purported unit on the linear packet becomes constant and that packet is
the genuine coordinate plane `x=0`.  This is the canonical hostile showing
that `phi!=0` is load-bearing.

## 6. Exact verification

Run

```bash
python3 04-computation/jacobian_universal_reducible_graph_component_unit_no_go_kps_s188.py
python3 -O 04-computation/jacobian_universal_reducible_graph_component_unit_no_go_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion verifies the
target identities `(11)--(12)`, source factorization `(6)--(7)`, both Bezout
certificates `(8)` and `(15)`, the degree rows `(10)` and `(16)` through six
independent controls, the first-resonance regression `(17)`, and the zero-row
hostile `(18)`.  The two certificate quotients have respectively 8 and 152
terms; no numerical inference is used.

**QED.**
