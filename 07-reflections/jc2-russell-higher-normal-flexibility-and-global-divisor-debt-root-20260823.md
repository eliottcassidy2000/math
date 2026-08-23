# Russell higher-normal flexibility and global divisor debt

**Status (2026-08-23): THM-3856, THM-3860, THM-3861, THM-3867, and THM-3868
are PROVED + VERIFIED-EXACT + independently hostile-audited; THM-3862 is
PROVED + CITED + independently proof-audited.  JC(2) remains OPEN.**

This session began from the apparent algebraization wall in THM-3846: every
unimodular arm jet lifts formally, but the canonical common-normal
resummation acquires a nonsquare.  The decisive change of viewpoint was to
stop treating that common normal coordinate as intrinsic.  The free tangent
coefficient at every higher order is a real coordinate, not disposable
gauge.  It changes the rationality type of the resummation.

The result is a sharper three-level frontier:

```text
finite formal neighborhoods: unobstructed;
canonical linear-normal resummation: nonrational for a colliding arm;
higher-normal rational resummation: possible, but globally polar;
global polynomial Darboux pair on the Russell surface: still open.          (1)
```

## 1. Inheritance pass and live concept board

The closest proved mechanism is THM-3846's canonical strip
`Bhat=k[s][[z]]`, `{z,s}=1`.  The canonical hostile is THM-3843: every
hypothetical Russell arm is a noninjective immersed normalization.  The
corrected near miss is to read THM-3846's square obstruction as invariant;
the theorem itself never made that overclaim, and THM-3860 supplies the
explicit higher-normal hostile.  The least-used sidecar is the finite
normalization between the Darboux target and the Russell ring; it remembers
deleted ramification that the arm conductor forgets.

The session board was:

| Object | Representation | Preserved invariant | Missing coordinate |
|---|---|---|---|
| Russell surface `Y` | `B` and `k(s,z)` | Poisson bracket, `Pic(Y)=Z/3`, constant units | effectivity at every height-one divisor |
| arm `L` | immersion `(a,b)` plus Bezout normal row | normalization and first jet | higher tangent gauges |
| formal Darboux lift | coefficient sequence in `k[s][[z]]` | bracket to every order | convergence/rationality/global regularity |
| conductor `kappa` | weighted collision graph | finite branch contacts | deleted ramification and infinity address |
| finite completion | `R subset S subset B` | degree and divisorial valuations | visible/deleted sheet decomposition unless `(e_i,f_i)` is retained |
| cubic factor packet | tower divisor `v` and factor-cofactor law | exact support equations | intrinsic identification with `kappa` |

Every useful pull in the session changed at least two board entries.  The
quartic normal classification closed the first three nontrivial
representation cells; the rational lift exposed the missing divisor
coordinate; the
completion contract showed why `kappa` alone cannot pay that divisor.

## 2. The formal lift space is an affine tangent torsor

For

```text
A=sum a_n(s)z^n,                     C=sum b_n(s)z^n,
```

the exact coefficient law is

```text
[z^m]{A,C}
 =sum_(i+j=m+1)(i a_i b_j'-j a_i'b_j).                       (2)
```

If `a_1b_0'-a_0'b_1=1`, the order-`m` equation has all solutions

```text
(a_(m+1),b_(m+1))
 =-R_m/(m+1)(a_1,b_1)+tau_m(s)(a_0',b_0').                    (3)
```

The forced part is normal to the arm; `tau_m` is tangent to it.  At the
first new order,

```text
(a_2,b_2)=-(W/2)(a_1,b_1)+tau_1(a_0',b_0'),
W=a_1b_1'-a_1'b_1.                                           (4)
```

THM-3846 chooses `tau_m=0` in a coherently resummed common-normal coordinate.
That choice is canonical and useful, but it is only one point of the affine
torsor `(3)`.  Thus a finite-order search which sets tangent terms to zero is
not testing the full formal problem.

The same formal-flexibility/global-algebraization split has now appeared in
the independent inverse-discriminant lane.  Provisional THM-3855 finds a
surjective formal coefficient-gradient recursion for one-place cubic
discriminants, while provisional THM-3853 says the first homogeneous
quadratic truncation cannot realize its target.  This is a structural
analogy, not a map between the two moduli problems: the Russell recursion is
symplectic and normal-adic, while the other is a discriminant inverse problem
in the target maximal ideal.

## 3. Polynomial normal depth through four closes completely

THM-3856 classifies every pair `A,C in k[s,z]` with

```text
deg_z A,deg_z C<=2,                  J_(z,s)(A,C)=lambda!=0.    (5)
```

The top Wronskian makes the quadratic coefficient row a constant target
direction.  After a constant target change, every genuinely quadratic pair
is exactly

```text
C=b(s)+beta z,
A=rho C^2+dC-(lambda/beta)s+a0,        rho,beta in k*.          (6)
```

It has the polynomial inverse

```text
s=(beta/lambda)(rho C^2+dC+a0-A),
z=(C-b(s))/beta.                                               (7)
```

All zero-top and zero-linear cases reduce to the same triangular endpoint.
Therefore a polynomial normal strip of depth at most two cannot restrict to
a self-identifying arm.  This result is uniform in the degree of `b(s)`;
the closed coordinate is transverse depth, not arm complexity.

THM-3861 continues the classification through `deg_z<=3`.  Writing the two
cubic coefficient rows and equating the six Jacobian buckets first makes the
top cubic row a constant target direction.  After a constant target change,
the only potentially new branch has transverse bidegree `(3,2)`.  Its two top
rows have the Kummer form

```text
(p,v)=(P h^3,V h^2),                    P,V in k*.             (C1)
```

The lower buckets integrate to a single coordinate `X=beta/h` and a shifted
arm parameter `T=b-X^2/(4V)`, with final equation

```text
h(MT+e)T'=lambda.                                            (C2)
```

If `h` is nonconstant, polynomiality forces `X` to have a pole over every
prime of `h`, and the unique cubic pole in the arm coefficient cannot cancel.
If `h` is constant, `(C2)` says a nonconstant polynomial and its derivative
multiply to a unit, again impossible.  Thus the `(3,2)` branch is empty; the
remaining `(3,1)` branch is explicitly triangular, and every degree drop is
THM-3856.  Consequently no polynomial normal strip through transverse degree
three can realize the self-identifying Russell arm.

THM-3867 closes the next depth.  Its `(4,1)` and `(4,2)` rows reduce by the
target shears `A-rho C^4` and `A-rho C^2` to THM-3861.  In the new `(4,3)`
row the top coefficients have Kummer scale `(r,q)=(Rh^4,Qh^3)`.  After a
moving scale and cubic depression, all lower buckets integrate to one rational
coordinate `W` with autonomous final equation

```text
lambda/h=-(3Q/(16R^2))
          (W^2-2EW+4RS^2/W^3)W'.                           (C3)
```

At a prime of nonconstant `h`, `W` can have a pole, a zero, or vanish
identically.  Canceling the corresponding pole of the cubic arm in those
three channels leaves respectively the nonzero quartic-arm residues
`-RB^2/(9Q^2)` or `-Rg^4/3`.  Constant `h` makes the final polynomial product
a unit and is also impossible.  Thus the genuine `(4,3)` row is empty.

The boundary is exact.  These theorems do not treat:

```text
one coordinate of z-degree >4;
both coordinates of z-degree 5;
rational functions in (s,z);
infinite z-adic expansions of global elements of B.                           (8)
```

Thus normal degree five is the first remaining polynomial cell.

A second hostile audit reaches the same closure by a genuinely different
route: after a determinant-one target shear, irreducible-prime valuations in
the lower buckets exclude nonconstant `(3,2)` Kummer factors, while the
constant factor gives `J=(KB+d)B'`, never a nonzero scalar.  Agreement between
that local-DVR proof and the global `X=beta/h` pole proof is strong evidence
that the degree-three boundary is structural rather than an artifact of one
normal form.

## 4. The square gate moves, but the fixed nodal seed is globally rigid

For a constant nonzero Wronskian `W=w`, take a rational function `phi(z)`
with

```text
phi(0)=0,             phi'(0)=1,             phi''(0)=-w.       (9)
```

Put

```text
Z=phi(z),
S=s/[(1+w phi(z))phi'(z)]+f(z),              f(0)=f'(0)=0.     (10)
```

Then

```text
J_(Z,S)(a(S)+alpha(S)Z,b(S)+beta(S)Z)=1+wZ,
J_(z,s)(Z,S)=1/(1+wZ),                                        (11)
```

so the composite bracket is exactly one.  The conditions `(9)` make
`S=s+O(z^2)` and `Z=z+O(z^2)`, preserving the arm and first normal row.

For the minimal nodal packet, the Mobius choice

```text
phi=z/(1+wz/2),
S=s(1+wz/2)^3/(1+3wz/2)                                      (12)
```

gives a pair in `Frac(B)^2`.  This is the decisive hostile against promoting
the THM-3846 squareclass to an invariant obstruction.  The new pair still
has a pole: in the explicit nodal formula, `A` has order `-1` on the Laurent
prime `z+4c^3`.

More generally, every rational `phi` satisfying `(9)` has a finite nonzero
regular point where

```text
(1+w phi)phi'=0.                                               (13)
```

Otherwise `phi+1/w` would be `gamma/(1+az)^d`; its first two jets would force
`(d+1)/d=-1`.  At the divisor `z=z_0`, the coefficient of `s` in `(10)` has a
pole.  Since `s=1/(3r)-c^3/(3z_0^3)` is nonconstant on that Laurent divisor,
no scalar `f(z)` cancels the leading coefficient.  The quadratic nodal
coordinate then has twice the pole order.

Allowing `Z_s!=0` really does move the pole, but THM-3868 proves that it cannot
remove it for this fixed nodal seed.  Eliminating `Z` from

```text
A_0=9c^6S^2-Z/(3c^3),
C_0=27c^9S^3-3c^3S-(3/2)SZ
```

gives the monic recovery cubic

```text
S^3+[(2-3A)/(9c^6)]S+2C/(27c^9)=0,                         (R1)
f'(S)=[2/(9c^6)](1+Z/(2c^3)).                              (R2)
```

If a rational composite `A,C` were regular in the normal Russell ring `B`,
`(R1)` would force `S in B` and then `Z in B`.  The density identity would
make `1+Z/(2c^3)` a unit; since `B*=k*`, `Z` would be constant, contradicting
the bracket.  This excludes every rational precomposition of the fixed seed,
not just the vertical subclass.

The formal freedom remains genuine.  An exact rational family with
`D=1-eta s z^2` has `Z_s` first at order `z^3`, satisfies the density exactly,
and moves the old divisor `z+4c^3` to a mixed `(r,z)` divisor; its nodal output
still has a pole there.  The next live operation is therefore a **different
seed whose hidden-control recovery is not monic**, or a nonrational control
with a separate normalization audit.  Merely adding mixed tangent gauge to
the same seed is now closed.

## 5. Proved finite-completion contract

THM-3862 records the following conditional consequences of a hypothetical
Russell Darboux map.  The theorem has received a separate proof audit; the
conditional object remains hypothetical because JC(2) remains open.

Assume a Darboux map `Psi=(P,Q):Y->A2` exists.  Put

```text
R=k[P,Q],                 K=Frac(B),
S=normalization of R in K,                  Xbar=Spec S.        (15)
```

Because every element of `S` is integral over `R subset B` and `B` is
normal, `S subset B`.  The induced birational quasi-finite map
`Y->Xbar` is an open immersion.  Normal surface Cohen--Macaulayness and
miracle flatness make `Xbar->A2` finite locally free of rank
`d=[K:Frac(R)]`.

Two independent necessities follow.

First, `S/R` cannot be monogenic.  If `S=R[theta]=R[T]/(f)` with `deg f=d>1`,
then etaleness on `Y` makes `f'(theta)` a unit in `B`.  THM-3785 gives
`B*=k*`.  But the unique power-basis expression of `f'(theta)` has nonzero
`theta^(d-1)` coefficient `d`, so it cannot be scalar in characteristic
zero.

Second, every irreducible branch component `Gamma` of the finite completion
must have affine normalization `A1`.  A ramification prime above `Gamma` is
omitted from the etale open `Y`.  The valuation lemma of THM-3841 sends
`Gamma` into the nonproperness set of the polynomial plane composite
`Psi o phi`, where `phi:A2->Y` is the Russell atlas.  Jelonek--Lason's
Theorem 3.2 covers that set by polynomial curves.  A dominant `A1` curve
lifts through the normalization of `Gamma`; extending to projective
normalizations forces genus zero and exactly one missing point.  Over `C`,
Nguyen Van Chau's one-point-at-infinity theorem further forces all branch
components, and the arm component already in the same nonproperness set, to
share one target infinity address.

Primary sources for the two external steps are
[Jelonek--Lason](https://arxiv.org/html/1411.5011v2) and
[Nguyen Van Chau](https://arxiv.org/abs/math/0305088).  The
[Dubouloz--Palka](https://arxiv.org/abs/1701.01425) nonproper etale maps are
`Y->Y` hostiles to naive pseudo-plane properness, not maps `Y->A2` and not a
contradiction to this contract.

The arm conductor and deleted ramification must not be identified.  If the
arm image itself were a branch component, THM-3849 gives two visible etale
contributions over its generic point: `L` contributes one and the nonempty
residual divisor `D` contributes at least one.  An omitted ramified prime
contributes at least two.  Hence

```text
Gamma_arm a branch component            ==>            d>=4.  (16)
```

In particular, at surface degree three the arm remains a Jelonek component
caused by the escaping source atlas, but it cannot be a branch component of
the finite completion.  The polynomial `kappa` records visible `L`--`D`
contact; it does not record the missing ramification prime, the different, or
the common infinity address.

## 6. Exact connection ledger

| Source | Target | Map | Preserved | Destroyed | Needed sidecar |
|---|---|---|---|---|---|
| arm jet | formal lift | recursion `(3)` | bracket to each order | rationality and poles | tangent gauges plus height-one valuations |
| canonical lift | rational higher-normal lift | symplectic precomposition `(10)` | arm and first jet | global regularity | divisor of every denominator |
| arm curve | collision graph | implicit gradient and `kappa` | finite contacts and delta weights | branch status/infinity/deleted sheets | finite completion and `(e_i,f_i)` |
| deleted prime | target branch | `E -> pi(E)` | valuation and nonproperness component | sheet decomposition | ramification/residue degrees |
| discriminant profile | formal cubic order | inverse-gradient recursion (provisional THM-3855) | formal discriminant | polynomial termination, units, atlas | index/different plus global divisor audit |
| bichromatic tower | arm polynomial | specialization in THM-3839 | exact support equations | intrinsic collision meaning | compare `div(v)` with `div(kappa)`, not just degrees |

The repeated mechanism is not “formal methods fail.”  Formal methods expose
the complete local freedom.  The failure occurs when a local solution is
asked to be effective simultaneously at every global divisor.

## 7. Generated task portfolio

### Anchor -- break fixed-seed monic recovery

THM-3868 shows that arbitrary mixed tangent gauge does not escape the fixed
nodal seed: regular outputs recover the controls integrally, and the
power-basis derivative is the forbidden density unit.  Parameterize the
smallest alternative seeds with the same arm and first-normal packet, eliminate
one hidden control, and classify the first point at which either

```text
recovery ceases to be monic, or seed Jacobian != power-basis derivative.   (17)
```

The cheap hostile is a seed differing only by a target automorphism, which
preserves the trap.  A positive signal must exhibit a genuine loss of
integrality while retaining an exact density control; a changed formula alone
is not enough.

### Niche -- classify quintic polynomial normal depth

THM-3867 closes `deg_z<=4`, so `deg_z=5` is the first polynomial cell.  Derive
all ten Jacobian buckets without bounding the `s`-degrees, remove every
target-shear row, and isolate the genuine `(5,4)` Kummer packet.  The quartic
proof needed simultaneous regularity of two arm coefficients; at depth five
three leading pole channels can tie, so the decisive object is their complete
initial-form ideal rather than another one-coordinate valuation.

### Geometry -- build a nonmonogenic all-`A1` branch completion

The completion must have:

```text
every reduced branch component normalized by A1;
one common target point at infinity;
no unit represented by its index/different;
constant units on the etale open.                              (19)
```

THM-3852 warns that one affine-line component is insufficient: its affine
profile always carries a nonpolynomial companion.  Provisional THM-3853/3855
instead provide a one-place formal inverse-discriminant laboratory.  The
next decisive test is polynomial termination together with the index-form
unit test, not another puncture count alone.

### Contact -- compare tower zeros with the conductor

THM-3839's nonconstant Kummer tower `v` is the cheapest remaining sparse
profile.  For the fixed nodal arm, compare the valuation packet of `v` with
`kappa=-(t^2-1)`.  A degree match is meaningless; the test must preserve each
collision address and the `8/7` tower multiplicities.  The strongest possible
outcome is an impossibility from incompatible local orders.  The honest
stopping outcome is that `v` controls support algebra while `kappa` controls
target-curve contact, requiring an additional map.

### Incoming signal -- retain the ordered reduction word and its cocycle

The incoming continuant audit gives a useful warning: a scalar digit mean or
product loses the sign holonomy, while the ordered Euclidean word survives
only when paired with its native cocycle.  On the JC side, proved
THM-3736/3740 already close constant Cohn tails and one polynomial right
factor by source automorphisms.  Thus the first live Cohn grammar needs at
least two interacting polynomial factors.  This is not a transfer to the
Russell problem, but it suggests a cheap comparison: compose two elementary
higher-normal factors in both orders and record

```text
(ordered factor word, determinant cocycle, complete pole divisor).        (20)
```

If swapping the factors changes the pole divisor while preserving scalar
degree/product data, then an ordinal such as transverse depth is only a task
scheduler, not a structural rank.  The conductor `kappa` gives the matching
hostile: it retains total incident collision weight but loses edge allocation
and deleted-sheet data.  Any word-based construction must therefore carry a
target-specific valuation cocycle rather than another scalar compression.

### Wildcard -- formal flexibility as an algebraization detector

Compare the Russell recursion `(3)` with the provisional inverse-
discriminant recursion of THM-3855.  Ask for a common obstruction language in
terms of height-one valuations of the first finite truncation whose local
formal correction is surjective.  Do not identify their moduli spaces or
transfer a theorem; retain the normal-coordinate versus coefficient-gradient
sidecar.

## 8. Stopping certificate

The session did not construct a global Darboux pair and did not prove or
disprove JC(2).  It did prove the complete polynomial classification through
transverse degree four, a rational higher-normal lift, the all-rational
precomposition barrier for its fixed nodal seed, and the nonmonogenic
finite-completion branch contract.  The old square obstruction survived only
in its exact canonical class; the stronger invariant is the fixed seed's
monic recovery/different-unit trap.  The next honest frontiers are an
alternative seed or nonrational control, while the orthogonal bounded
polynomial frontier begins at transverse degree five.
