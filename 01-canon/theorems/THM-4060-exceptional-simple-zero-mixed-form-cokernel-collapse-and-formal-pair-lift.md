---
id: THM-4060
title: "Exceptional simple-zero mixed-form cokernel collapse and formal pair lift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with the strict
  formal/local boundary stated below.  On the exceptional affine
  three-branch packet, the mixed closed-form cokernel is the four-dimensional
  filtered space K[[A]]/A^4, rather than the one-class-per-degree fixed-a
  cokernel.  Every normalized simple-zero displacement H=t+O(t^2) therefore
  has an all-order formal local mixed pair with source Jacobian 12.  This
  proves no convergence, algebraization, global polynomial pair, or
  consequence for JC(2).
source: jc2-double-zero-rebuild-20260824 / mixed-pair sidecar, 2026-08-25
audit: >
  PASS.  The characteristic-zero companion replays the exact Q(alpha)
  engine of THM-4054, verifies the moving-endpoint transgression and carrier
  factors, reconstructs every closed coefficient triple through degree four,
  and obtains the complete coupled ranks for cutoffs 0,...,8.  A standalone
  finite-field implementation reconstructs the exceptional branches,
  triangle, fixed-output, mixed, and polynomial-exact banks at all four roots
  modulo 137; for every cutoff 4,...,12 the mixed and exact ranks are
  rows-4, the fixed rank is rows-(N+1), and all visible pure stable powers
  lie in the mixed image.  Both normal/optimized pairs match their stored
  transcripts byte for byte; both scripts have zero Python assert nodes and
  zero float literals.  The split-prime calculation is an independent
  hostile audit, not the characteristic-zero all-order proof.
depends_on:
  - THM-4058-exceptional-affine-triangle-period-and-simple-zero-monomial-ladder
  - THM-4054-exceptional-affine-simple-zero-retained-packet
related:
  - THM-3623-russell-cylinder-even-general-vertical-fold-all-order-closure
  - THM-3630-russell-cylinder-noneven-formal-survivor-no-finite-jet-bound
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
script: 04-computation/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060.py
independent_script: 04-computation/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060_independent_audit.py
output: 05-knowledge/results/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060.out
independent_output: 05-knowledge/results/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060_independent_audit.out
script_sha256: 5faa1db5d9aa127e71f22f6bbb4ceee56b3dfe7addd5eb26b2416d04bfe8098a
independent_script_sha256: f168c3d8a5d6cd4f92a6d2a18ddd4708c1f29da142e00eab2b294426f16bce2d
output_sha256: 803ee72dcfc77a60020bcafc85c67b2c5f04836315dc7c8963700ffee5927d63
independent_output_sha256: 03faaf8109e2d3469550888c3bb6fc2f9339b1c3e0742cf2d693a872a455ba79
hash_basis: raw LF bytes
---

# THM-4060 -- the fixed-output ladder transgresses away in the mixed complex

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with a strict
formal/local boundary.**  THM-4058's growing fixed-output obstruction is real,
but it is the derivative of a first-output circulation carrier in the mixed
closed-form complex.  Four lower classes survive; no normalized simple-zero
reparametrization hits them.

All notation is inherited from THM-4058.  In particular, work over the
exceptional field `K=Q(alpha)`, let `A=a+3/4`, and rename the affine stable
coordinate `w` there as `u`.  The three affine source images are

```text
u=R_i(A,c),                   i=0,1,2.                 (1)
```

The pairwise intersection vertices, their orientation, the nonzero field
element `delta=(26/15)rho`, and the triangular operators `Pi,Lambda` are
exactly those of THM-4058.

## 1. The mixed exact-form operator

For a target germ `f in K[[A,c,u]]`, write

```text
f_i(A,c)=f(A,c,R_i(A,c)).                              (2)
```

The affine Darboux pair is `(A,-4c)`.  Its mixed linearized bracket on the
three branches is

```text
D_i(f,g)=partial_c g_i-4 partial_A f_i.                (3)
```

Here both derivatives are total derivatives of the restrictions in `(2)`,
with the other branch coordinate held fixed.  Since
`Jac_(x,u)(A,c)=-3`, the corresponding source-density operator is `-3D`;
thus `(3)` is the bracket normalization, not a change in the source
Jacobian convention of THM-4054.

The operator `(3)` is also exactly the restriction of the closed target
two-form

```text
-4 df wedge dc+dA wedge dg.                            (4)
```

In fact `(4)` presents every closed formal target two-form.  Write an
arbitrary one as

```text
omega=P dA wedge dc+Q dA wedge du+S dc wedge du.       (5)
```

It is closed precisely when

```text
P_u-Q_c+S_A=0.                                         (6)
```

Choose `f_u=S/4` and `g_u=Q` by formal integration in `u`.  Then
`E=P+4f_A-g_c` is independent of `u` by `(6)`.  Adding to `g` a formal
`u`-independent primitive whose `c` derivative is `E` makes `E=0`, and
gives `(4)=(5)`.  The integrations raise target degree by at most one, so
this representation respects the complete degree-at-most-`N+1` domain for
a retained cutoff `N`.

## 2. Circulation around the infinitesimal triangle

For a common target germ `f`, put

```text
I_f(A)=Pi((f_0,f_1,f_2)).                              (7)
```

All three vertices agree modulo `A^5`.  At fixed `A`, every edge therefore
has both `(c-c_*)` and `(u-u_*)` of order at least five.  Taylor-expand `f`
about the common vertex jet.  The constant and `c`-linear line integrals
telescope, the `u`-linear term is the signed triangle area, and every
higher term has still larger valuation.  Hence

```text
I_f in A^10 K[[A]].                                    (8)
```

THM-4058 gives

```text
I_u=Pi(u)=-(15/52)delta^2 A^10+O(A^11),                (9)
```

whose leading coefficient is nonzero.  Multiplication by arbitrary series
in `A` now proves the exact ideal identities

```text
{I_f:f in K[[A,c,u]]}=A^10 K[[A]],
{Lambda(f):f in K[[A,c,u]]}=A^5 K[[A]].               (10)
```

Differentiate `(7)` with respect to `A`.  Every moving-endpoint term at a
vertex cancels with the adjacent edge because the two branch restrictions
are values of the same target germ.  Therefore

```text
Pi(partial_A f)=dI_f/dA,
Lambda(partial_A f)=(d/dA+5/A)Lambda(f).               (11)
```

The apparent pole in `(11)` is harmless by `(10)`.

## 3. The growing fixed-output cokernel collapses to four classes

The `g` part of `(3)` is the fixed-`a` derivative of THM-4058, so its
complete cokernel is measured by `Lambda`.  Equations `(10)--(11)` show
that the additional first-output part has response ideal

```text
Lambda(D(f,0))
 =-4(d/dA+5/A) A^5K[[A]]
 =A^4K[[A]].                                           (12)
```

Indeed `A^(k+5)` is sent to the nonzero multiple `(k+10)A^(k+4)`.

Consequently

```text
coker D is isomorphic as a filtered K-vector space to
K[[A]]/A^4K[[A]],                                      (13)

r lies in im D
 iff [A^0]Lambda(r)=...=[A^3]Lambda(r)=0.              (14)
```

The filtered graded proof of THM-4058 is recursive: at degree `n`, its
rank-`3n+2` diagonal block supplies the next target correction without
changing lower degrees.  Thus cutoff solutions in `(14)` may be chosen
compatibly degree by degree; no bare inverse-limit surjectivity is assumed.

The finite-cutoff form is equally sharp.  If

```text
M_N=direct_sum_(i=0)^2 K[A,c]/(A,c)^(N+1),             (15)
```

then the complete target-jet bank in `(3)` has

```text
dim coker(D_N)=min(N+1,4),
rank(D_N)=3 binom(N+2,2)-min(N+1,4).                   (16)
```

This is a characteristic-zero all-`N` statement.  It follows from the
complete fixed-output criterion of THM-4058 and the triangularly invertible
map in `(12)`, not from extrapolation of a rank table.

The first-output carriers can be made completely explicit.  For every
integer `s>=1`, use `f=u^s`.  From THM-4058,

```text
Lambda(u^s)=-(s/2)rho A^(s+4)+O(A^(s+5)).              (17)
```

and `(11)`,

```text
Lambda(D(u^s,0))
 =2s(s+9)rho A^(s+3)+O(A^(s+4)).                       (18)
```

The coefficient never vanishes in characteristic zero because `rho!=0`.
Thus `(18)` kills successively every fixed-output period class of degree at
least four and explains why precisely four lower classes survive.

## 4. Every normalized simple zero has a formal mixed pair

Let

```text
H(t)=t+O(t^2) in K[[t]],
u=H(t),                    t=h(u).                     (19)
```

This is a formal target-coordinate change.  On the affine branch packet
`(1)`, a pair with constant source Jacobian `12` must have bracket density

```text
-4h'(R_i) dA wedge dc.                                 (20)
```

The required correction to the affine base form is the common target
scalar restriction

```text
r_i=-4(h'(R_i)-1).                                     (21)
```

By `(10)`, `Lambda(r)` belongs to `A^5K[[A]]`; in particular its first four
coefficients vanish.  Criterion `(14)` therefore supplies formal target
germs `f,g` satisfying

```text
D(f,g)=r.                                              (22)
```

The unit at the next step can be tracked explicitly.  First solve

```text
-4(d/dA+5/A)lambda=Lambda(r).                          (22a)
```

Here `lambda` belongs to `A^6K[[A]]`.  Since `Lambda(u)` is `A^5` times a
unit, take

```text
f=(lambda/Lambda(u))u in (A,c,u)^2.                    (22b)
```

The residual has zero complete period, and the recursive fixed-output lift
may take `g in (A,c,u)^2`.  Define

```text
Omega=-4dA wedge dc-4df wedge dc+dA wedge dg.          (23)
```

This form is closed, its `dA wedge dc` coefficient is a unit, and `(22)`
gives on every branch

```text
Omega|_i=-4h'(R_i)dA wedge dc.                         (24)
```

The auxiliary germs in `(22)` are not asserted to make
`(A+f,-4c+g)` a pair: that wedge would contain the additional term
`df wedge dg`.  The actual, generally different pair comes from factoring
the closed form `(23)`.

For completeness, the formal Darboux step is elementary here.  The kernel
of the nonvanishing closed two-form `(23)` has a generator
`-4partial_u+O(A,c,u)`.  Recursive formal straightening over `K` makes it
`partial_U`; only units and positive integers are divided out.  Contraction
then removes all `dU` terms, while closedness makes the remaining coefficient
independent of `U`; hence `Omega=b(P,Q)dP wedge dQ`.  Taking

```text
F=P,                 G=integral_0^Q b(P,s) ds          (25)
```

gives `Omega=dF wedge dG`.

Finally substitute `u=H(w)` in these target germs.  On the source,

```text
dA wedge dc=-3 dx wedge du,
du=H'(t)dt,                 h'(H(t))H'(t)=1.           (26)
```

Equations `(24)--(26)` give

```text
Jac_(x,t)(F(A,c,H(w)),G(A,c,H(w)))=12                 (27)
```

simultaneously in the three completed source germs.  This is an all-order
formal mixed pair for every normalized simple-zero displacement `(19)`.
After an embedding `K->L`, the same construction works for coefficients in
any characteristic-zero extension `L`; a complex-coefficient `H` produces a
pair after that scalar extension, not one uniform pair over `K`.

## 5. The leading monomial repair

For

```text
H(t)=t+gamma t^m,             m>=2,                   (28)
```

formal inversion gives

```text
h'(u)=1-m gamma u^(m-1)+O(u^(2m-2)).                  (29)
```

The leading response of the right side `(21)` is

```text
-2m(m-1)gamma rho A^(m+3).                            (30)
```

Equation `(18)` cancels it with the first-output term

```text
f_lead=-((m-1)/(m+9))gamma u^m.                        (31)
```

Thus THM-4058's sharp fixed-`a` obstruction is not merely absent from the
full mixed tangent bank at its first occurrence: every member of its
monomial ladder lies in the transgressed range, and `(23)--(27)` integrate
the repair to a possibly different formal pair.

## 6. Audit gate

The characteristic-zero companion replays the frozen THM-4054 exact field
engine and checks `(8)--(18)` and `(30)--(31)`.  Its complete coupled matrix
table is

| `N` | rows | columns | rank | cokernel |
|---:|---:|---:|---:|---:|
| 0 | 3 | 6 | 2 | 1 |
| 1 | 9 | 18 | 7 | 2 |
| 2 | 18 | 38 | 15 | 3 |
| 3 | 30 | 68 | 26 | 4 |
| 4 | 45 | 110 | 41 | 4 |
| 5 | 63 | 166 | 59 | 4 |
| 6 | 84 | 238 | 80 | 4 |
| 7 | 108 | 328 | 104 | 4 |
| 8 | 135 | 438 | 131 | 4 |

It also builds the complete degree-at-most-four kernel of `(6)`, of
dimension `85`, and reverse-reconstructs every basis vector by `(4)`.  Its
serialized table hash is

```text
1fefca39e8602eaae4b2e1095a603e0457e475665ffa28b0bd90faa39bf38fde. (32)
```

The standalone hostile companion uses no repository import.  Over

```text
F_137,              alpha=44,82,92,134,                (33)
```

it independently reconstructs the triangle, with

```text
rho=(8,85,12,135),       delta=(23,56,103,97).          (34)
```

For every root and every `N=4,...,12`, it obtains fixed-output cokernel
`N+1`, mixed and polynomial-exact cokernel `4`, equality of the mixed and
exact images, and membership of every visible `1,u,...,u^N`.  It also checks
the carrier in the equivalent basis `A^k u`, `k=0,...,12`.  Its table hash is

```text
6833a2c709ed77ab45f3e02cd3bdda198e649d5c9f546740b9ccc602962ff0fe. (35)
```

This finite-field calculation is an independent hostile control.  The
characteristic-zero all-order proof is Sections 1--5.

## 7. Reproduction

```bash
python3 04-computation/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060.py
python3 -O 04-computation/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060.py
python3 04-computation/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060_independent_audit.py
python3 -O 04-computation/jc2_exceptional_simple_zero_mixed_form_cokernel_formal_pair_thm4060_independent_audit.py
```

Both normal/optimized pairs are byte-identical to their stored transcripts.
The four raw-byte script/output hashes are pinned in the frontmatter.

## 8. Connection contract and strict boundary

| field | contract |
|---|---|
| source | closed formal target two-forms, equivalently the mixed exact-form presentation `(4)` about `(A,-4c)` |
| target | three retained branch-density germs, or their cutoff-`N` quotients |
| map | restriction to the affine graphs `(1)`, followed by the normalized triangle period of THM-4058 |
| preserved predicate | exact finite-cutoff membership, the all-order formal cokernel, and constant source Jacobian after a formal simple-zero target reparametrization |
| destroyed information | convergence, algebraicity, polynomial degree, regularity outside the collision chart, behavior on other source fibres, injectivity, properness, and control at infinity |
| needed sidecar | convergence/algebraization and removal of the denominator in `a=e/(b+4)`, followed by a global collision and Keller-map analysis |
| cheapest decisive tests | the area coefficient `(9)`, carrier response `(18)`, rank formula `(16)`, closed-form identity `(6)`, and direct all-root split-prime membership |

This theorem concerns the formal completion at the one exceptional
three-point collision.  It does not construct polynomial or convergent
target functions, a global pair on the Russell cylinder, a polynomial
endomorphism of the plane, a Keller counterexample, or any consequence for
`JC(2)`.  The global problem remains **OPEN**.
