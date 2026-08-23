---
id: THM-3836
title: "The nonlinear cubic chart is a factor-cofactor Darboux packet"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every polynomial
  Keller atlas of the THM-3811 nonlinear cubic
  surface satisfies an exact cubic factorization P=C*S and cofactor equation
  delta(C)=lambda*P.  The two factors are comaximal and allocate all
  components of the three cubic pencil members.  A differentiated-determinant
  argument excludes every allocation by whole pencil members, so at least
  one cubic fibre has nonempty comaximal factors mapping to both intrinsic
  G_m arms.  Their ideals are the exact etale base changes of those arms;
  every irreducible factor has degree at least two, so the pencil member has
  degree at least four.  Equality is two disjoint smooth G_m conics.  This
  removes THM-3831's irreducible-h hypothesis in all degree.
source: root + jc_quartic_c3_construct / homogeneous root-ratio and component-selection lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn-boundary, 2026-08-23).
  The audit rederived chart clearing and the Keller cofactor, checked both
  comaximality arguments scheme-theoretically, retained whole irreducible
  multiplicities, verified every j=1,2,3 monochromatic endpoint, grounded
  the final nonunit contradiction in THM-3822's actual h-arm, and checked
  the two intrinsic arm labels under etale base change.  Normal and optimized
  runs byte-match the frozen transcript and both raw hashes agree.  The
  companion verifies P=C*S from the
  triangular chart, the cleared Keller cofactor, all spectral resultants and
  arm values, the three Euler/cofactor identities, boundary monomials, and
  characteristic-zero multipliers.  A second independent hostile audit
  (root / audit_thm3836) proved the exact base-change ideals, smoothness,
  reducedness, component-degree floor, and equality classification.
depends_on:
  - THM-3811-ramification-class-unit-criterion-and-nonlinear-cubic-packet
  - THM-3831-intrinsic-spectral-pencil-fibre-atlas-and-forced-cubic-two-arm-hit
  - THM-3832-nonlinear-cubic-root-ratio-triangular-birational-chart
related:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
  - THM-3835-polynomial-marked-root-ratio-nonentry
  - THM-3830-coordinate-cross-bichromatic-split-nonentry
  - THM-3837-comaximal-line-hyperbola-affine-bichromatic-contact-nonentry
  - THM-3843-nonlinear-cubic-keller-atlas-total-degree-contradiction
script: 04-computation/jc2_cubic_factor_cofactor_darboux_packet_thm3836.py
output: 05-knowledge/results/jc2_cubic_factor_cofactor_darboux_packet_thm3836.out
script_sha256: 692f6caec2b3398da6ebb7cfd33f475178557375db25b29aa75d7aa3ad0d423e
output_sha256: 212befa0196033325c42f1f11eaa2fd94b79bc8b54be8f173ed7861a896ec0b4
semantic_sha256: 2aaac6405fb5979bb64afd6dc6548c66bf8e457e87728a7d9282e11cbd85071a
hash_basis: raw LF bytes
---

# THM-3836 -- the cubic factor/cofactor packet forces a two-arm fibre

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `K` of characteristic zero.
Let

```text
psi:A2_(x,y) -> U                                                (1)
```

be a dominant etale plane atlas of the THM-3811 surface, and suppose its
composite with `U -> A2_(A,C)` has

```text
Jac(A,C)=lambda in K*.                                          (2)
```

Write `R=K[x,y]` and use the same letters for the pulled-back intrinsic
functions.  In particular,

```text
Ck-mh=1.                                                         (3)
```

Define the homogeneous binary forms

```text
P=3h^3+7h^2k+k^3,
Q=7h^2+3k^2,
B=6h^3+7h^2k-k^3,                                               (4)

S=h^2 Q C-kB.                                                    (5)
```

Then the atlas satisfies the exact factor/cofactor packet

```text
P=C S,                                                           (6)
delta(C):=k Jac(h,C)-h Jac(k,C)=lambda P.                        (7)
```

Moreover `(C,S)=R`.  The three distinct roots `a_1,a_2,a_3` of

```text
r(Z)=3Z^3+7Z^2+1                                                (8)
```

give pairwise comaximal pencil members

```text
p_i=h-a_i k,                    P=3 p_1p_2p_3.                   (9)
```

The factorization `(6)` allocates every irreducible component, with its full
multiplicity, to exactly one of `C` and `S`.  At least one `p_i` is allocated
to **both** sides.  Consequently there are nonconstant `u_i,v_i in R` with

```text
p_i=epsilon u_i v_i,            (u_i,v_i)=R,                    (10)
```

such that the two nonempty fibre schemes map respectively to

```text
u_i=0  ->  U_i^-: C=0,
v_i=0  ->  U_i^+: C=-1/a_i^2.                                  (11)
```

These are the two intrinsic `G_m` arms of THM-3831.  Thus every Keller atlas,
including the reducible-`h` branch, must meet both arms of some cubic
spectral fibre.  The two schemes are smooth and reduced, every irreducible
component has degree at least two, and

```text
deg(h-a_i k)>=4.                                              (11a)
```

Equality in `(11a)` consists of two disjoint smooth affine conics, each
isomorphic to `G_m`.  No degree or support bound is assumed.

## 1. Clearing the triangular chart

On the dense function-field chart of THM-3832 put `z=h/k` and

```text
r=3z^3+7z^2+1,          q=7z^2+3,
b=6z^3+7z^2-1,          s=z^2qC-b.                              (12)
```

The chart formula `k=r/(Cs)` and homogeneity give

```text
P=k^3r,                    S=k^4s.                               (13)
```

Hence `CS=P` in the function field, and therefore in the domain `R`.  The
same chart turns the Keller equation into

```text
Jac(z,C)=lambda k r.                                            (14)
```

But

```text
Jac(z,C)=[k Jac(h,C)-h Jac(k,C)]/k^2,
kr=P/k^2.                                                        (15)
```

Clearing `k^2` proves `(7)`.  These identities use the rational ratio only;
THM-3835 shows that treating `z` as a source polynomial would lose the live
case.

## 2. The two global factors are comaximal

Suppose a maximal ideal of `R` contained both `C` and `S`.  Equation `(6)`
would put `P` in that ideal, so `h=a_i k` at the corresponding geometric
point for some root `a_i` of `(8)`.  The determinant `(3)` makes `k` nonzero
there.  With

```text
b(Z)=6Z^3+7Z^2-1,             q(Z)=7Z^2+3,                     (16)
```

exact spectral arithmetic gives

```text
b+q=2r,       Res(r,Zqb)!=0.                                    (17)
```

Thus `b(a_i)=-q(a_i)!=0`, whereas `(5)` at `C=0,h=a_i k` gives

```text
S=-b(a_i)k^4!=0,                                                 (18)
```

a contradiction.  Nullstellensatz proves `(C,S)=R`.

Similarly, if `p_i=p_j=0` for `i!=j`, then `(a_i-a_j)k=0`, hence
`h=k=0`, contradicting `(3)`.  Therefore the three ideals `(p_i,p_j)` are
the unit ideal.  Unique factorization applied to `(6),(9)` now gives the
claimed component allocation.  This argument retains reducible and
nonreduced pencil members: coprimality forces the complete multiplicity of
each irreducible factor onto one side.

## 3. Whole-member allocation is impossible

Call `p_i` **monochromatic** when all its irreducible factors go to the same
side.  Suppose all three members were monochromatic.  Since `(2)` makes `C`
nonconstant, unique factorization would give

```text
C=gamma product_(i in I) p_i,
gamma in K*,                  1<=j:=|I|<=3.                     (19)
```

The right side is homogeneous of degree `j` in `(h,k)`.  For every
homogeneous `F(h,k)` of degree `j`, the chain rule and Euler identity give

```text
delta(F)=j F Jac(h,k).                                          (20)
```

Combining `(7),(19),(20)` with `P=CS` and cancelling the nonzero `C` yields

```text
j Jac(h,k)=lambda S.                                            (21)
```

All roots `a_i` are nonzero by `(17)`.  Reducing `(19)` modulo `h` therefore
gives, for a nonzero scalar `c_I`,

```text
C=gamma c_I k^j mod h.                                         (22)
```

The determinant `(3)` says that `h` divides `Ck-1`, hence

```text
h divides gamma c_I k^(j+1)-1.                                 (23)
```

Take the Jacobian bracket with `h`.  Because `j+1` is nonzero in
characteristic zero, `(23)` implies

```text
h divides k^j Jac(h,k),
h divides k^j S                                                   (24)
```

by `(21)`.  But `(3)` gives `gcd(h,k)=1`, while `(5)` gives the exact boundary

```text
S mod h = k^4.                                                   (25)
```

Thus `gcd(h,S)=1`, and `(24)` would make `h` a unit.  This is impossible:
THM-3822 gives the proper intrinsic arm `B/(h)=K[k,k^-1]`, so `h` is a
nonunit and nonconstant in `B`; dominance of `(1)` makes `B -> R` injective,
so its pullback cannot be a scalar unit.  The contradiction proves that some
`p_i` is not monochromatic.

## 4. Exact base-change ideals and the component-degree floor

For that `i`, let `u_i` be the product of the factors assigned to `C` and
`v_i` the product assigned to `S`, retaining multiplicities.  Both are
nonconstant, and a common zero would give `C=S=0`; hence `(u_i,v_i)=R` and
the CRT gives the scheme decomposition `(10)`.

Absorb scalar units and write

```text
p_i=u_i v_i,              C=u_i C_0,             S=v_i S_0.   (25a)
```

If a maximal ideal contained both `v_i,C_0`, it would contain both `C,S`,
contrary to `(C,S)=R`; similarly for `u_i,S_0`.  Hence

```text
(v_i,C_0)=R,              (u_i,S_0)=R,                        (25b)
(p_i,C)=(u_i),            (p_i,S)=(v_i).                      (25c)
```

The ideals `(25c)` say scheme-theoretically that `V(u_i)` and `V(v_i)` are
the exact base changes of the two intrinsic arms.  Etale base change makes
them smooth and reduced; in particular `u_i,v_i,p_i` are squarefree.  Every
irreducible component maps nonconstantly to its target `G_m`, since an etale
map is quasi-finite.  The pulled-back Laurent coordinate `k` is consequently
a nonconstant unit on every component.

No component can be a line: its coordinate ring would be `K[t]`, whose only
units are scalars.  Thus every irreducible factor of both `u_i` and `v_i` has
degree at least two.  Both sides are nonempty, which proves `(11a)`.

If equality holds, `u_i` and `v_i` are irreducible quadratics.  Their
projective closures are smooth conics.  Each affine conic has a nonconstant
unit, so its boundary cannot be a single point; Bezout with the line at
infinity gives exactly two distinct boundary points.  Over `K` it is
therefore `P^1` minus two points, hence `G_m`.  Comaximality makes the two
conics disjoint.  The pair `xy-1,xy-2` shows that this geometric degree floor
is sharp before the global factor/cofactor equations are imposed.

## 5. Exact arm labels and scope

On `p_i=0`, the determinant reads

```text
k(C-a_i m)=1,                                                   (26)
```

so `k` is a unit.  The `C`-side has `C=0`.  On the `S`-side, use
`h=a_i k`, `b(a_i)=-q(a_i)`, and `q(a_i)!=0` in `(5)` to obtain

```text
S=k^4 q(a_i)(a_i^2 C+1).                                       (27)
```

Thus `S=0` gives `C=-1/a_i^2`, proving `(11)` with the exact signs.  Since
`psi` is etale, each nonempty source component maps nontrivially into the
corresponding intrinsic Laurent arm.

This closes the monochromatic, reducible-`h`, and linear-component escape
routes.  THM-3837 independently closes its line--hyperbola selector model;
the sharp uneliminated geometry before the remaining equations is
conic--conic.  THM-3843 subsequently combines this theorem's global packet
with total-degree arithmetic to exclude the full polynomial plane atlas.
No general Jacobian counterexample or proof is claimed.
**QED.**
