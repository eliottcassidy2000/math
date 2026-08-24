---
id: THM-3928
title: "Split affine conic cannot support an irreducible one-place torus sextic"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Let an
  irreducible nonlinear degree-N torus branch have affine normalization A1,
  and suppose its quadratic coefficient is the product of two distinct
  affine lines. Parallel lines are impossible. For nonparallel lines, the
  induced map from the branch normalization to the Cardano cusp has degree
  at least ceil((N+1)/2). In degree six its degree is 4, 5, or 6, with exact
  line-degree packets (6,2), (6,4), or (6,6). For a general chosen component
  these are genuine high-fold candidates. When the branch is the full
  irreducible sextic (so `deg q<=3`), the first has multiple infinity support,
  the second fails its weighted degree equation, and the third requires an
  impossible binary-form identity. Hence no distinct-affine-line conic
  supports a full irreducible one-place torus sextic. On the quadratic
  resolvent, splitting the conic supplies only one intrinsic three-torsion
  direction: the universal divisor-relation lattice is Z plus Z/3, not two
  copies of Z/3. In the same classical full-sextic grammar a double-line
  conic factors the discriminant into degree-at-most-three pieces. Thus all
  affine singular conics are closed for the full irreducible sextic, while
  arbitrary-`q` high-fold and double-line components and a conic containing
  the infinity line remain open.
source: jc_degree6_one_place / post-THM-3925 singular splitting-conic lane, 2026-08-23
audit: >
  TWO INDEPENDENT HOSTILE AUDITS PASS (jc-cohn3709 and
  incoming_thm3928_audit/root, 2026-08-23). They reconstructed the Cardano
  integrality and parallel-line UFD steps, projective degree/contact ledger,
  all three full-sextic row contradictions, and every divisor/Smith-form
  implication. The first audit repaired the distinct-line statement by
  separating the general componentwise fold bound from the full irreducible
  sextic closure. The second found the remaining double-line overstatement:
  degree-at-most-three also requires the classical `deg q<=3` full-
  discriminant grammar. MISTAKE-469 and the explicit `y=x^6` hostile freeze
  the surviving arbitrary-`q` boundary. LF-normalized normal and optimized
  streams match the frozen LF output in all 820 active gates; raw script,
  output, and semantic hashes agree.
related:
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3925-fivefold-conic-contact-torus-sextic-one-place-fold
  - THM-3926-unit-ideal-cubic-primitive-class-genus-two-boundary
script: 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
output: 05-knowledge/results/jc2_split_affine_conic_one_place_fold_degree_thm3928.out
script_sha256: b26b18721aa67ecd5212c1edc0bf39f28505db6f631b3dc50de998e12c07a2b6
output_sha256: fbed1442286d6b3ffa1756d9efa1e606a1774f984e740d1b34afc8164ec2c978
semantic_sha256: b8f9d6a7fe151c72db445042197047033803066df0adfb566445b0fb27e65690
hash_basis: raw LF bytes
---

# THM-3928 -- splitting the conic multiplies fold degree, not C3 directions

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Let

```text
p=ell_1 ell_2 in k[x,y],                 Delta=4p^3-27q^2,   (1)
```

where `ell_1,ell_2` are nonassociate nonconstant affine-linear forms and
`q in k[x,y]`. Let `Gamma` be an irreducible nonlinear component of
`Delta=0` such that `p` is nonzero in `k(Gamma)`. Assume that the projective
normalization of `Gamma` is `P1` and that exactly one normalization point
lies above the line at infinity. Equivalently, the normalization of the
affine curve is

```text
nu:A1_t -> Gamma.                                             (2)
```

Assume also that the coefficient map `(p,q)` is nonconstant on `Gamma`.
Write `N=deg Gamma`.

Then the zero lines of `ell_1,ell_2` cannot be parallel. If they are
nonparallel, the induced finite map from `Gamma^nu=A1_t` to the normalization
of the Cardano cusp has degree `d` satisfying

```text
d >= ceil((N+1)/2).                                          (3)
```

In particular, for a sextic `N=6`,

```text
d in {4,5,6}.                                                 (4)
```

After exchanging the two lines, the degrees of their pullbacks are exactly

```text
d=4: (6,2),             d=5: (6,4),             d=6: (6,6). (5)
```

This first gives a fold-degree barrier. Section 4 then eliminates all three
rows in `(5)` under the additional assumption that `Gamma=V(Delta)` is the
full irreducible sextic, proving that a distinct-affine-line conic cannot
support any irreducible nonlinear one-place torus sextic. Equivalently in
this quadratic-`p` setting, the closure assumes `deg q<=3` and that `Delta`
has no additional component. The componentwise fold bound `(3)` does not
need this extra hypothesis.

The quadratic-resolvent ledger is equally rigid. Under the normality and
generic-splitting hypotheses stated in Section 5, the divisor relations
intrinsic to the Cardano norm have presentation

```text
Z^4 / <(1,1,0,0),(0,0,1,1),(3,0,3,0)> = Z direct_sum Z/3.   (6)
```

Thus the two components supply one diagonal three-torsion candidate and one
free anti-diagonal direction. They do not by themselves supply two
independent three-classes. A second three-class can arise only if additional
global divisor data imposes a new relation on the free direction.

## 1. The Cardano parameter is polynomial

Pull `(1)` back to `A1_t`. In `k(t)` one has

```text
4(p o nu)^3=27(q o nu)^2.
```

Define

```text
h=3(q o nu)/(2(p o nu)).                                    (7)
```

Direct cancellation gives

```text
p o nu=3h^2,                         q o nu=2h^3.            (8)
```

The first identity says that `h` is integral over `k[t]`: it satisfies the
monic equation

```text
H^2-(p o nu)/3=0.
```

Since `k[t]` is integrally closed in `k(t)`,

```text
h in k[t].                                                   (9)
```

The nonconstancy assumption makes `d=deg h` positive. Notice the exact
scope: one place proves `(9)`, not `d=1`.

The image of `(p,q)|Gamma` is dense in the standard cusp
`4P^3-27Q^2=0`, whose normalization is

```text
h |-> (3h^2,2h^3).
```

Therefore

```text
[k(t):k(h)]=deg h=d.                                        (10)
```

This is the precise fold degree. In particular, the coefficient map is
birational on `Gamma` if and only if `d=1`.

## 2. Parallel components are impossible

Put

```text
a_i=ell_i o nu in k[t].                                     (11)
```

Equation `(8)` gives

```text
a_1 a_2=3h^2.                                               (12)
```

Suppose the two affine lines are parallel. After rescaling their equations,

```text
a_2=lambda a_1+c,                         lambda*c!=0.       (13)
```

Thus `gcd(a_1,a_2)=1`. The UFD `k[t]` and `(12)` force each `a_i` to be a
scalar times a square. Since `k` is algebraically closed, absorb the scalar
square roots and rewrite `(13)` as

```text
g^2-f^2=c.                                                   (14)
```

But `(g-f)(g+f)=c` is a nonzero scalar, so both factors, and hence `f,g`,
are constant. Then `h` is constant by `(12)`, a contradiction. No
nonconstant torus branch in this scope exists for two distinct parallel
components.

There is also an ambient warning. For intersecting lines, `dp=0` at their
intersection. For parallel lines, after an affine change `p=u(u-c)`, so
`dp=0` on `u=c/2`. Hence the coefficient map `(p,q)` can never have
everywhere nonzero Jacobian when `p` is a distinct-line affine conic. This
does not assume or invoke the planar Jacobian conjecture.

## 3. The projective degree forces a large fold

Now suppose the lines are nonparallel. Then

```text
(x,y) |-> (ell_1,ell_2)                                     (15)
```

is an affine automorphism. Thus `(a_1(t),a_2(t))` is still a polynomial
birational parametrization of a degree-`N` curve with one place at infinity.
If

```text
e_i=deg a_i,
```

homogenizing this parametrization to degree `max(e_1,e_2)` gives a
basepoint-free map `P1 -> P2`. It is birational onto `Gamma`, so pulling back
a generic line proves

```text
N=max(e_1,e_2).                                             (16)
```

Neither `a_i` can be constant: otherwise `Gamma` lies in an affine line and,
being irreducible of dimension one, is that line. Hence `e_1,e_2>=1`.
Taking degrees in `(12)` now gives

```text
e_1+e_2=2d.                                                 (17)
```

Equations `(16)-(17)` immediately imply

```text
2d>=N+1,
```

which is `(3)`.

The same ledger records the unique infinity point geometrically. A line
with `e_i<N` passes through that infinity address with contact `N-e_i`.
Since two nonparallel lines have different projective infinity points, at
most one `e_i` is below `N`; equivalently, one of them equals `N`, as already
forced by `(16)`.

## 4. Every sextic fold row fails

Assume now that `Gamma=V(Delta)` is the full irreducible degree-six curve;
in particular `deg q<=3`. Set `N=6` and exchange the lines so that `e_1=6`.
Equation `(17)` becomes

```text
e_2=2d-6.                                                    (18)
```

The inequalities `1<=e_2<=6` give precisely the three rows `(4)-(5)`:

| fold degree `d` | line degrees `(e_1,e_2)` | contact of the second line at infinity |
|---:|---:|---:|
| 4 | `(6,2)` | `4` |
| 5 | `(6,4)` | `2` |
| 6 | `(6,6)` | `0` |

In particular `d=1` is not merely absent; `d=2,3` are absent as well. It
remains to test the three high-fold rows. Write `q_3` for the homogeneous
cubic part of `q` in the affine coordinates `(ell_1,ell_2)=(x,y)`.

### 4.1 Fold degree four has more than one infinity point

For `(e_1,e_2,d)=(6,2,4)`, the weighted degrees of `x^3` and `x^2y` on the
normalization are respectively `18` and `14`, whereas

```text
deg(2h^3)=12.                                                (18a)
```

No other cubic monomial has either forbidden degree, so their coefficients
in `q` vanish. Consequently

```text
q_3=y^2(alpha x+beta y).
```

The degree-six part of the branch equation is

```text
Delta_6=y^3[4x^3-27y(alpha x+beta y)^2].                    (18b)
```

It vanishes at `[1:0]`, while the bracket is nonzero there and, being a
nonconstant binary cubic, has a zero elsewhere. Thus `Delta_6` has at least
two projective support points. An irreducible degree-six curve with only one
normalization place at infinity has only one projective infinity point, so
this row is impossible.

### 4.2 Fold degree five misses degree fifteen

For `(e_1,e_2,d)=(6,4,5)`, the evaluated degrees of the four cubic monomials
are

```text
deg(x^3,x^2y,xy^2,y^3)=(18,16,14,12).                      (18c)
```

The degree-`18` and degree-`16` terms are unique, so their coefficients must
vanish in an identity with `2h^3`, whose degree is `15`. Every remaining
monomial of total degree at most three has evaluated degree at most `14`.
It cannot produce degree `15`. This row is empty.

### 4.3 Fold degree six has no one-support binary top form

For `(e_1,e_2,d)=(6,6,6)`, one-place infinity would force

```text
Kx^3y^3-27q_3(x,y)^2=L(x,y)^6,                  K!=0,       (18d)
```

after allowing the nonzero rescaling induced by affine coordinate changes.
There is no such binary-form identity. If `L` is one of the coordinate
forms, say `L=x`, write

```text
q_3=ax^3+bx^2y+cxy^2+dy^3.
```

The `y^6` and `x^2y^4` coefficients successively give `d=c=0`; the
`x^3y^3` coefficient then gives `K=0`, a contradiction.

If both coefficients of `L` are nonzero, rescale the coordinates to
`L=x+y`. The `x^6,x^5y,xy^5,y^6` coefficients give

```text
-27a^2=-27d^2=1,                  b=3a,       c=3d.         (18e)
```

The `x^4y^2` coefficient gives `ad=-1/27`, hence `d=a`. Therefore
`q_3=a(x+y)^3` and `-27q_3^2=(x+y)^6`; equation `(18d)` again forces
`K=0`. This exhausts all linear forms `L` and kills the last row.

Hence no irreducible one-place sextic in the theorem's scope exists. The
fold-degree table is not a list of survivors; it is the complete finite
pre-obstruction list dispatched by `(18a)-(18e)`.

The bound is genuinely about a nonlinear branch. For example,

```text
p=3xy,                              q=2x^3
```

has the line `y=x` as a component of `Delta=0`; on that line
`h=x` has degree one. This is exactly the excluded linear case and
shows where the proof would fail: one of the pullback line forms becomes
proportional to the other.

## 5. Splitting creates one torsion direction, not two

Let `rho^2=27`, and let `B` be the normalization of the integral quadratic
resolvent

```text
k[x,y,w]/(w^2-(27q^2-4p^3)).                                (19)
```

Assume that `q mod ell_i` is nonzero for both `i`. At the generic point of
each line the quadratic cover is split and unramified, with two height-one
primes

```text
D_i^+=(ell_i,w+rho q),                 D_i^-=(ell_i,w-rho q). (20)
```

Put `gamma=rho q+w`. Its norm under `w |-> -w` is

```text
gamma gamma_bar=4p^3.                                      (21)
```

At `D_i^+`, the conjugate factor is a unit, so `(21)` gives valuation three;
at `D_i^-`, `gamma` is a unit. There are no other height-one zeros. Therefore

```text
div(ell_1)=D_1^++D_1^-,
div(ell_2)=D_2^++D_2^-,
div(gamma)=3D_1^++3D_2^+.                                  (22)
```

The relation matrix in the ordered basis
`(D_1^+,D_1^-,D_2^+,D_2^-)` has Smith form

```text
diag(1,1,3,0),                                              (23)
```

proving `(6)`. After using the first two relations, write

```text
x=[D_1^+],                 y=[D_2^+].
```

The norm relation is only

```text
3(x+y)=0.                                                   (24)
```

Thus `x+y` is the diagonal Kummer candidate while `x-y` remains a free
direction in the universal presentation. The map from this presentation to
the actual `Cl(B)` may have further kernel. In particular `(24)` alone does
not prove that its torsion candidate is nonzero. Conversely, if the actual
class group contains two independent three-classes generated by these
components, some additional global principal-divisor or saturation relation
must turn the free direction into three-torsion. That is a separate invoice,
not a benefit obtained merely by splitting `p`.

More generally, for `r` squarefree components the analogous norm packet is

```text
Z^(2r) / <D_i^++D_i^- (1<=i<=r), 3 sum_i D_i^+>
       = Z^(r-1) direct_sum Z/3.                            (25)
```

There is still only one intrinsic three-torsion direction.

## 6. Double lines: the classical full sextic factors, but components survive

The other affine singular conic is even more rigid. Choose `sigma in k` with
`sigma^2=3`. If `p=ell^2`, then identically

```text
4p^3-27q^2=(2ell^3-3sigma q)(2ell^3+3sigma q).              (26)
```

The identity alone gives no degree bound because the standing component
hypotheses allow arbitrary `q`. Under the additional classical torus-sextic
assumptions

```text
deg q<=3,                 Delta itself irreducible of degree six,           (27)
```

both factors in `(26)` have degree at most three. Even if one is scalar,
zero, or repeated, the full discriminant is reducible or nonreduced; it is
never an irreducible sextic. Thus double lines are empty in the classical
full-discriminant grammar `(27)`.

They are not empty for the broader component grammar. Put

```text
ell=x,        f=y-x^6,        q=(2x^3-f)/(3sigma).                         (28)
```

Then

```text
4x^6-27q^2=f(4x^3-f).                                                       (29)
```

The component `Gamma=V(f)` is an irreducible degree-six curve with
normalization `t |-> (t,t^6)` and one place at infinity. On it, `p=t^2` is
nonzero in the function field and `q=2t^3/(3sigma)` is nonconstant. This is
the hostile witness separating full-discriminant reducibility from exclusion
of a chosen component. Arbitrary-`q` double-line components therefore remain
open; see MISTAKE-469.

Their resolvent lattice remains a useful reducible hostile. If
`q mod ell` is nonzero, then

```text
div(ell)=D^++D^-,                         div(gamma)=6D^+.
```

The corresponding presentation has Smith form

```text
Z^2/<(1,1),(6,0)> = Z/6.                                  (30)
```

Its three-primary Kummer class is `2[D^+]`; again there is only one
three-direction. This is a class-lattice control, not an exclusion of the
surviving arbitrary-`q` component `(28)`.

Likewise, if one projective component of the conic is the chosen line at
infinity, the affine polynomial `p` is linear rather than a product of two
affine line equations. Its missing projective divisor must be retained in
the boundary ledger, and `(3)` is not asserted. Thus the arbitrary-`q`
distinct-line high folds and double-line components are the honest remaining
affine-component design space, while for the full classical sextic only this
infinity-component case remains.

Reproduce the exact arithmetic packet with

```bash
python3 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
python3 -O 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
```

After platform newlines are normalized to LF, both streams must byte-match
the frozen LF output named in the metadata.
**QED.**
