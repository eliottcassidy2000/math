---
id: THM-3928
title: "Split affine conic forces a high-degree Cardano fold on a one-place branch"
status: >
  PROVED + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING. Let an
  irreducible nonlinear degree-N torus branch have affine normalization A1,
  and suppose its quadratic coefficient is the product of two distinct
  affine lines. Parallel lines are impossible. For nonparallel lines, the
  induced map from the branch normalization to the Cardano cusp has degree
  at least ceil((N+1)/2). In degree six its degree is 4, 5, or 6, with exact
  line-degree packets (6,2), (6,4), or (6,6). Thus a split affine conic can
  never give a birational/no-fold coefficient presentation of a one-place
  sextic. On the quadratic resolvent, splitting the conic supplies only one
  intrinsic three-torsion direction: the universal divisor-relation lattice
  is Z plus Z/3, not two copies of Z/3. Extra independent three-classes need
  an additional global relation. A double-line conic makes the torus
  discriminant a product of two cubics, so cannot underlie an irreducible
  sextic. This closes the double-line and birational distinct-line lanes,
  but degree-4/5/6 folds and conics with an infinity component remain open.
source: jc_degree6_one_place / post-THM-3925 singular splitting-conic lane, 2026-08-23
audit: >
  PROVISIONAL PROOF CANDIDATE AWAITING INDEPENDENT HOSTILE AUDIT. The proof
  separates the normalization-integrality step, parallel-line UFD
  obstruction, projective degree ledger, fold-degree interpretation, and
  quadratic-resolvent divisor presentation. The assertion-free exact
  companion verifies the Cardano identities, all degree tables, the split
  and double-line Smith forms, the double-line discriminant factorization,
  the arbitrary-component extension, and sharp excluded controls. Normal and
  optimized replay must byte-match the frozen output before promotion.
related:
  - THM-3844-two-cusp-polynomial-branch-quadratic-resolvent-design-gate
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-3920-affine-plane-boundary-unibranch-depressed-cubic-chart-obstruction
  - THM-3925-fivefold-conic-contact-torus-sextic-one-place-fold
script: 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
output: 05-knowledge/results/jc2_split_affine_conic_one_place_fold_degree_thm3928.out
script_sha256: ab66d5f711b4291861b8354c0d533223ba27fc3456ecc79999ddcd4d92c46aa7
output_sha256: 2b1d8bec1d41f132c39bb01e76f34c140f8c6845ed64f60bdab7f5b15d1f48db
semantic_sha256: 5d4acb8186b838b177e2a6951f9a96e1fb580a903f8c4fe99ca394fce2c608cd
hash_basis: raw LF bytes
---

# THM-3928 -- splitting the conic multiplies fold degree, not C3 directions

**PROVED + VERIFIED-EXACT; INDEPENDENT HOSTILE AUDIT PENDING.** Work over an
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

Consequently a distinct-line splitting conic cannot produce a birational
coefficient presentation of a nonlinear one-place torus sextic. The fold is
not a small defect: it has degree at least four.

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

## 4. Exact sextic table and the no-fold obstruction

Set `N=6` and exchange the lines so that `e_1=6`. Equation `(17)` becomes

```text
e_2=2d-6.                                                    (18)
```

The inequalities `1<=e_2<=6` give precisely the three rows `(4)-(5)`:

| fold degree `d` | line degrees `(e_1,e_2)` | contact of the second line at infinity |
|---:|---:|---:|
| 4 | `(6,2)` | `4` |
| 5 | `(6,4)` | `2` |
| 6 | `(6,6)` | `0` |

In particular `d=1` is not merely absent; `d=2,3` are absent as well. Any
split-affine-conic sextic design must carry a fourfold-or-higher polynomial
coefficient fold. If the intended construction requires the coefficient
map to be birational along its branch, this lane is empty.

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

## 6. Double lines factor; high folds and the infinity component remain

The other affine singular conic is even more rigid. Choose `sigma in k` with
`sigma^2=3`. If `p=ell^2`, then identically

```text
4p^3-27q^2=(2ell^3-3sigma q)(2ell^3+3sigma q).              (26)
```

Both factors have degree at most three. Even if one is scalar, zero, or
repeated, the result is reducible or nonreduced; it is never an irreducible
sextic. Thus double lines are empty for the branch geometry sought here.

Their resolvent lattice remains a useful reducible hostile. If
`q mod ell` is nonzero, then

```text
div(ell)=D^++D^-,                         div(gamma)=6D^+.
```

The corresponding presentation has Smith form

```text
Z^2/<(1,1),(6,0)> = Z/6.                                  (27)
```

Its three-primary Kummer class is `2[D^+]`; again there is only one
three-direction. This is a class-lattice control, not a surviving sextic.

Likewise, if one projective component of the conic is the chosen line at
infinity, the affine polynomial `p` is linear rather than a product of two
affine line equations. Its missing projective divisor must be retained in
the boundary ledger, and `(3)` is not asserted. Together with the
degree-`4,5,6` distinct-line folds, this is the honest remaining
singular-conic design space.

Reproduce the exact arithmetic packet with

```bash
python3 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
python3 -O 04-computation/jc2_split_affine_conic_one_place_fold_degree_thm3928.py
```

Both streams must byte-match the frozen output named in the metadata.
**QED.**
