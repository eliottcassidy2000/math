# Russell higher-normal flexibility and global divisor debt

**Status (2026-08-23): THM-3856, THM-3860, THM-3861, THM-3867, THM-3868,
and THM-3871 are PROVED + VERIFIED-EXACT + independently hostile-audited;
THM-3862 is PROVED + CITED + independently proof-audited.  The alternative
quadratic-seed computation is an independently verified FINITE-EXACT scout,
not a theorem or Keller construction.  JC(2) remains OPEN.**

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
bounded normal classifications through quintic depth closed four successive
nontrivial representation cells; the rational lift exposed the missing
divisor coordinate; the completion contract showed why `kappa` alone cannot
pay that divisor.

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

## 3. Polynomial normal depth through five closes completely

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

THM-3871 closes quintic depth.  Its genuinely new rows begin from the
root-free UFD packet

```text
(w,c_j)=(R h^5,Q h^j),                  j=2,3,4.             (C4)
```

The `(5,2)` row has an unavoidable local arm residue `3R/8`.  The `(5,3)`
row acquires one conserved polynomial and a constant one-form; after the only
possible balanced scaling, arm cancellation and the conserved equation have
resultant `441`.  The `(5,4)` row acquires two conserved polynomials.  Its
`W`-pole, unit, finite-zero, identically-zero, and constant-scale channels are
all exhaustive, and the balanced arm/conserved system has resultant

```text
3171942400000 = 2^23*5^5*11^2.                             (C5)
```

This extra conserved data is load-bearing.  The point
`(X,Y,Z)=(1,3/10,-23/10)` cancels both leading arm equations, so an arms-only
analysis has a genuine false positive; the two conserved residuals are
`476/25` and `321/50`.  Thus degree five is the first depth at which the
leading divisor picture alone loses the theorem.

The boundary is exact.  These theorems do not treat:

```text
one or both coordinates of z-degree >=6;
rational functions in (s,z);
infinite z-adic expansions of global elements of B.                           (8)
```

Thus normal degree six is the first remaining bounded polynomial cell.

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

The first alternative-seed scout now identifies exactly how monic recovery
can fail without pretending that failure is Keller.  In normalized controls
`x=Z/(3c^3)`, `u=3c^3S`, consider the complete affine-`p,q` quadratic cell

```text
A=u^2-x+p(u)x^2,
C=u^3-u-(3/2)ux+q(u)x^2.                                  (Q1)
```

If `R(u;A,C)` is the Sylvester eliminant and

```text
L=(3/2)p(u)u-q(u),
M=q(u)(u^2-A)-p(u)(u^3-u-C),                              (Q2)
```

then exact elimination gives

```text
pR=pM^2+LM+(u^2-A)L^2,             R_u|_(Q1)=L D,          (Q3)
```

where `D=J_(x,u)(A,C)`.  Every affine stratum has scalar leading coefficient
except one: for `e!=0`,

```text
p=-4e^2,                  q=-4e^2u+e,                       (Q4)
```

the degree-three recovery coefficient is
`e(16e^2(A-1)-1)/2`.  This is genuine nonintegrality, not merely a bad
primitive element.  The symplectic pole chart

```text
u=1/t,
x=1/(2et)-1/(8e^2)-t/(4e)+v t^2
```

makes `A,C` polynomial while `u,x` both have valuation `-1`.  But on `t=0`
the seed Jacobian is `16e^2v`, not a unit constant.  Thus the scout finds the
first exact escape from THM-3868's monic-recovery trap and simultaneously
explains why it is not a Keller construction.

A distinct family
`p=alpha(u-u0)`, `q=(3/2)alpha u0(u-u0)` keeps `u` and
`w=alpha(u-u0)x` integral while `x` may pole.  Its two boundary differents are
`(3X+2)/2` and `(3X-2)/2`.  The conditional Keller identity
`D{u,X}=2w-1` becomes a useful divisor/unit test only after proving that a
source component realizes the full `X`-line; that missing geometric step is
why no contradiction is claimed.  Finally, the statement that quadratic
corrections are “gauge” is only a unique formal normal-coordinate 2-jet: it
need not be symplectic or an automorphism, and exact equality already leaves
a cubic `C` residual and a quartic `A` residual.

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

### Anchor -- impose constant density on the first nonintegral seed

The affine-quadratic scout exhausts the first alternative-seed cell.  Its
unique nonintegral locus `(Q4)` pays for recovery failure with the nonconstant
boundary Jacobian `16e^2v`.  The next task is therefore not another quadratic
coefficient search.  Quotient by the unique formal normal-coordinate 2-jet,
retain the first genuine cubic residuals, and solve simultaneously for

```text
hidden-control nonintegrality,          D in k*,
and effectivity at every pole divisor.                                (17)
```

The cheapest bounded cell adds affine-in-`u` coefficients of `x^3` to `(Q1)`
and recomputes both the eliminant leading divisor and the full Jacobian.  A
positive signal must keep the pole chart while making `D` constant.  If that
cell is empty, move to a seed with two integral sidecars and nonmonogenic
recovery rather than changing primitive element again.

### Niche -- classify sextic polynomial normal depth

THM-3871 closes `deg_z<=5`, so `deg_z=6` is the first bounded polynomial
cell.  Derive all twelve Jacobian buckets without bounding the `s`-degrees,
remove target-shear rows, and classify the genuine `(6,j)` packets.  The
quintic false positive proves that arms and top Kummer equations are no longer
a sufficient state: every conserved penultimate bucket and both original arm
coefficients must stay in the packet.  The cheapest hostile deletes one such
bucket and should be required to manufacture a spurious cancellation family.

### Geometry -- build a nonmonogenic all-`A1` branch completion

The completion must have:

```text
every reduced branch component normalized by A1;
one common target point at infinity;
no unit represented by its index/different;
constant units on the etale open.                              (19)
```

THM-3852 warns that one affine-line component is insufficient: its affine
profile always carries a nonpolynomial companion.  Incoming proved THM-3873
extends that failure to the first fixed-coordinate non-graph triangular
parabola: every polynomial transverse quotient forces a distinct reduced
companion whose normalization omits at least two projective points.  Hence
that marked-root family cannot satisfy the one-infinity condition in `(19)`.
This is a map from one explicit branch profile to its discriminant companion,
not a classification of all finite completions.  Provisional THM-3853/3855
instead provide a one-place formal inverse-discriminant laboratory.  The next
decisive test is polynomial termination together with the index-form unit
test in a higher-triangular or non-profile branch, not another puncture count
alone.

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
sidecar.  At the time of this reflection THM-3872 was still audit-pending; it
has since been independently hostile-audited and promoted.  Its useful
procedural signal remains: quotient cosmetic additions by the entire ideal
vanishing at the marked addresses before choosing a minimal affine section,
then record the residual left by that section.  Applied to `(Q1)`, this says
the cubic-seed task must carry the full zero-jet addition ideal alongside the
chosen normal 2-jet; otherwise a “unique” residual may be an artifact of the
section.

## 8. Stopping certificate

The session did not construct a global Darboux pair and did not prove or
disprove JC(2).  It did prove the complete polynomial classification through
transverse degree five, a rational higher-normal lift, the all-rational
precomposition barrier for its fixed nodal seed, and the nonmonogenic
finite-completion branch contract.  It also isolated, at FINITE-EXACT scout
scope, the unique affine-quadratic locus where hidden-control integrality
genuinely fails; the same pole chart shows its Jacobian is nonconstant, so it
is not a Keller construction.  The old square obstruction survived only in
its exact canonical class; the stronger invariant is the interaction among
recovery monicity, the different, and every boundary divisor.  The next
honest anchor begins with genuine cubic seed residuals or nonmonogenic hidden
coordinates, while the orthogonal bounded polynomial frontier begins at
transverse degree six.
