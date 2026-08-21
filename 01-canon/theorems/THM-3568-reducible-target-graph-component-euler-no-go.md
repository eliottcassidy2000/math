---
id: THM-3568
title: "Reducible target-graph component Euler no-go"
status: >
  PROVED + VERIFIED-EXACT.  For every nonzero h in C[b], neither irreducible
  component of the THM-3565 pullback associated to
  phi=-2h^3a+b h^2-2h is A2.  The degree-one component is exactly
  A2 minus V(3ah^2-bh+1), with a nonconstant unit.  If rho,mu,nu count the
  distinct roots of h, (bh)^2-2bh+4, and bh-2, respectively, the degree-two
  component has compactly supported Euler characteristic
  -3+3rho+2mu+nu: it is 2 for constant h and at least 5 otherwise.  Thus no
  reducible deg_a(phi)=1 target graph has a source-coordinate factor.
source: kps-s188
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-3560-jelonek-euler-gate-monomial-target-shear-no-go
  - THM-3565-resonant-linear-a-target-graph-factor-classification
related:
  - THM-3559-affine-target-coordinate-pullback-no-go
companion: 04-computation/jacobian_reducible_graph_component_euler_no_go_kps_s188.py
output: 05-knowledge/results/jacobian_reducible_graph_component_euler_no_go_kps_s188.out
script_sha256: e281372c6a6c3ee1a7689288847d44a032b152aebf5e20def393cb7289e72b52
output_sha256: 77845935c7a5414e9d90d77aef7858670f34518607d7b55e6fc78910019499d7
hash_basis: LF-normalized bytes
---

# THM-3568 -- reducible target-graph component Euler no-go

**PROVED + VERIFIED-EXACT.**  THM-3565 classifies every reducible target
graph with `deg_a(phi)=1`.  Here both of its components are excluded from
being coordinate planes, uniformly in the degree and singularities of the
parameter polynomial.

All varieties are over `C`; `chi` denotes compactly supported Euler
characteristic.  Target coordinates are `(a,b,c)` and source coordinates are
`(x,y,z)`.

## 1. The reducible family and its two components

Fix `0!=h in C[b]` and put

```text
phi=-2h^3a+b h^2-2h.                                      (1)
```

THM-3565 proves that the pullback of `T_h=V(c+phi)` has exactly two
irreducible components: one of generic degree one and one of generic degree
two over `T_h~=A2`.  The degree-one component has equation

```text
R_h=x-(1+xy)h(F2)=0.                                      (2)
```

Introduce the three target polynomials

```text
D=3ah^2-bh+1,
C=12a^2h^2-4abh+16a-b^2,
P=(bh)^2-2bh+4,                 T=bh-2.                  (3)
```

The Jelonek equation on the target graph has the exact factorization

```text
L_phi=D^2 C.                                               (4)
```

This square is not a multiplicity to be counted by Euler integration; the
reduced Jelonek section is `V(D) union V(C)`.

## 2. The degree-one component is a punctured plane

On the target open set `D!=0`, define

```text
x=h/D,
y=b-3ah,
z=aD^3-y^2D(D+3).                                        (5)
```

Direct substitution gives

```text
1+xy=1/D,             F1=a,             F2=b,
F3=-phi.                                                    (6)
```

Conversely, modulo `(2)`, exact polynomial identities give

```text
D(F1,F2)(1+xy)=1,           y=F2-3F1h(F2).               (7)
```

The remaining inverse coordinate is then `(5)`.  Hence the restriction of
`F` is an isomorphism

```text
R_h ~= A2_(a,b) minus V(D).                                (8)
```

The polynomial `D` is nonconstant and becomes a unit on `(8)`.  Thus the
coordinate ring of `R_h` has a nonconstant unit, whereas `C[A2]^*=C*`.
Therefore

```text
R_h is not A2.                                              (9)
```

This upgrades the constant-`h` Kummer cylinder from THM-3559 to the full
polynomial family.

## 3. Three reduced root supports

Let

```text
rho = # distinct roots of h,
mu  = # distinct roots of P=(bh)^2-2bh+4,
nu  = # distinct roots of T=bh-2.                          (10)
```

All counts are reduced supports; repeated roots count once.

The curve `A_h=V(D)` projects isomorphically to the complement of the roots
of `h` in the `b`-line: `D=0` has no point when `h=0`, and otherwise solves
uniquely for `a`.  Hence

```text
chi(A_h)=1-rho.                                            (11)
```

For `B_h=V(C)`, projection to the `b`-line has two points unless either
`h=0` or `P=0`, and one point on each of those disjoint supports.  Indeed

```text
disc_a(C)=64P,                    C|_(h=0)=16a-b^2.        (12)
```

Constructible Euler integration along this finite-fibre projection gives

```text
chi(B_h)=2-rho-mu.                                         (13)
```

The reduced intersection of the two curves is exactly the `T`-support.  The
identity

```text
h^2C+T^2=4(ah^2+1)D                                      (14)
```

shows that on `D=0`, where `h!=0`, the condition `C=0` is equivalent to
`T=0`.  Every root of `T` gives one point of `A_h intersect B_h`, so

```text
chi(A_h intersect B_h)=nu.                                (15)
```

This same support is precisely the omitted-curve intersection.  On the
omitted curve, parameterized by `a=b^2/12`, `c=4/(3b)`, the graph equation
restricts to

```text
c+phi(a,b)=-T^3/(6b).                                     (16)
```

Thus multiplicity three in `(16)` still contributes only the reduced count
`nu`, as required by THM-3560's constructible strata.

## 4. Fibre strata of the degree-two component

Write `Q_h` for the irreducible degree-two component from THM-3565.  Combine
the exact `3/1/0` fibre census of THM-2473 with the isomorphism `(8)`:

```text
target stratum                         # points of Q_h

A2 minus (A_h union B_h)                       2,
A_h minus B_h                                  1,
B_h minus A_h                                  0,
A_h intersect B_h                              0.          (17)
```

For example, over `B_h minus A_h` the complete fibre has one point and it
lies on `R_h`; the quadratic pair has escaped.  Over `A_h minus B_h`, the
degree-one component is absent and the unique complete-fibre point lies on
`Q_h`.  This typing prevents the common mistake of assigning the full
`3/1/0` census separately to each component.

Euler integration of `(17)` gives

```text
chi(Q_h)
 =2[1-chi(A_h)-chi(B_h)+nu]+[chi(A_h)-nu]
 =2-chi(A_h)-2chi(B_h)+nu
 =-3+3rho+2mu+nu.                                         (18)
```

## 5. The formula never equals one

If `h` is constant and nonzero, then

```text
(rho,mu,nu)=(0,2,1),              chi(Q_h)=2.             (19)
```

If `h` is nonconstant, then `rho>=1`.  The polynomial `bh` is nonconstant;
the two equations

```text
bh=1+i sqrt(3),                    bh=1-i sqrt(3)          (20)
```

have nonempty disjoint root supports, so `mu>=2`.  Likewise `bh=2` has a
root, so `nu>=1`.  Therefore

```text
chi(Q_h)>=-3+3+4+1=5.                                    (21)
```

In every case `chi(Q_h)!=1`; hence `Q_h` is not isomorphic to `A2`.
Together with `(9)`, neither irreducible factor is a coordinate polynomial.

Consequently:

> No reducible target graph with `deg_a(phi)=1` has a source-coordinate
> factor under the fixed THM-1300 map.

This closes the entire factor-bearing first-resonance family in all degrees.
It does not rule out an irreducible `deg_a(phi)=1` complete pullback being a
coordinate plane, nor any resonant degrees `4,7,...`; those remain separate
Euler/irreducibility problems.

## 6. Exact verification

Run

```bash
python3 04-computation/jacobian_reducible_graph_component_euler_no_go_kps_s188.py
python3 -O 04-computation/jacobian_reducible_graph_component_euler_no_go_kps_s188.py
```

The ordinary and optimized transcripts agree.  The companion verifies
`(4)--(7)`, `(12)`, `(14)`, and `(16)` symbolically.  It then recomputes all
three reduced root counts and `(18)` on constant, squarefree, multiple-root,
and repeated-omitted-root controls; `h=b^2-3` has
`bh-2=(b+1)^2(b-2)` and checks that `nu=2`, not three.

**QED.**
