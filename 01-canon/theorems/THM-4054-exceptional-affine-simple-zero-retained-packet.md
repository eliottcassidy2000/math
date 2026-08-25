---
id: THM-4054
title: "Exceptional simple-zero multigerm exact-form saturation and fixed-a fifth-order obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, in the finite/local scope
  stated below.  Over the exceptional quartic field, for both H=t and
  H=t+t^2, the complete cutoff-five retained image of arbitrary target
  two-forms equals the retained image of polynomial exact target two-forms;
  both images have rank 59 in 63 coordinates and contain the constant
  density.  The affine fold H=t has the explicit local Darboux pair
  (a,-4c).  For H=t+t^2, keeping F=a works through cutoff four and fails
  sharply at cutoff five, with a displayed nonzero quartic-field response.
  Nevertheless the full mixed tangent bank about (a,-4c) contains the
  required first-order correction -24t through cutoff five.  These are
  three-branch finite multigerm statements.  They construct no global
  polynomial pair, no coherent all-order lift, no Keller map, and no
  counterexample to JC(2).
source: jc2-double-zero-rebuild-20260824 / affine simple-zero side lane, 2026-08-24
audit: >
  PASS IN THE NARROW FINITE/LOCAL SCOPE.  Independent array arithmetic
  reconstructed the exceptional polynomial and the three source germs,
  checked Q_alpha'(-1,0,1)=(-9/2,1,9/2), verified the parity-free local
  identities a=e/(b+4)=q/D^2 and Jac_(x,q)(a,c)=-3, and recovered at the
  split prime 137 the cutoff-five ranks 59,59,57,59 for the arbitrary,
  polynomial-exact, fixed-a, and mixed-tangent banks.  A separate four-root
  modular companion matched all four split reductions; the exact K-result
  then specializes to all complex embeddings.  The production companion
  computes the promoted ranks, response, and norm exactly over
  K=Q(alpha), with byte-identical normal and optimized transcripts.  The
  audit explicitly rejected global, all-order, exact-form-implies-Darboux,
  and wholesale-THM-3629 readings.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
  - THM-3687-russell-cylinder-exceptional-quartic-actual-j0-lift
related:
  - THM-4058-exceptional-affine-triangle-period-and-simple-zero-monomial-ladder
  - THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction
  - THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary
  - THM-3631-russell-cylinder-noneven-closed-form-order-five-survival
  - THM-3737-russell-cylinder-exceptional-quartic-jacobian-image-hyperplane
script: 04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054.py
independent_script: 04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054_independent_audit.py
output: 05-knowledge/results/jc2_exceptional_affine_simple_zero_retained_packet_thm4054.out
independent_output: 05-knowledge/results/jc2_exceptional_affine_simple_zero_retained_packet_thm4054_independent_audit.out
script_sha256: 1e84ce5674c74c8b2d504fe8f57341be86ab20e1eab5c8fafe7494bb52d488a2
independent_script_sha256: cc12c81d8e4d614fbdc0146035cabc0c82b8a64c25ba53f0d9a27f969210b41d
output_sha256: 55a78d9bdff23d85058acff9a20e2903c801512bc68f90b2bdc3b1dfdb731005
independent_output_sha256: 96b2cf4ad92128c377876c193d16a2bc118edb741a13d4e9df01fe0937c26fed
hash_basis: raw LF bytes
---

# THM-4054 -- exceptional simple-zero exact-form saturation and the fixed-`a` fifth-order obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, with a strict
finite/local boundary.**  The exceptional collision polynomial is not even.
After repairing that typing error, its affine simple-zero boundary has a
clean three-part anatomy: an explicit local Darboux pair for `H=t`, complete
polynomial-exact saturation through cutoff five for both `H=t` and
`H=t+t^2`, and a sharp cutoff-five obstruction only when the first output is
frozen to the local coordinate `a`.  The latter obstruction disappears in
the full first-order mixed tangent space, but no nonlinear pair is thereby
integrated.

All rings and retained calculations are in characteristic zero unless a
finite-field hostile control is explicitly labelled.

## 1. Exceptional field, collision polynomial, and retained multigerm

Let

```text
K=Q[alpha]/(72783360 alpha^4-77822208 alpha^3-28419741 alpha^2
                                      +7849770 alpha-1276420).       (1)
```

Put

```text
q_0(x)=-3/4+x-(27/4)x^2-2x^3+(9/2)x^4+x^5,
P(x)=x^2-2x^4+x^6=x^2(1-x^2)^2,
Q_6=q_0-(259/36)P,
R_1=P-x^2P=x^2(1-x^2)^3,
R_2=4P-9xP=(4-9x)P,
theta=(520/9)alpha^2-(1688/81)alpha-5717/729,
Q_alpha=Q_6+theta R_1+alpha R_2.                       (2)
```

At the three retained source points `p=-1,0,1`, exact evaluation gives

```text
(Q_alpha(-1),Q_alpha(0),Q_alpha(1))=(-3,-3/4,-3),
(Q_alpha'(-1),Q_alpha'(0),Q_alpha'(1))=(-9/2,1,9/2).   (3)
```

In particular,

```text
Q_alpha'(0)=1, so Q_alpha is not even.                  (4)
```

Use the exponent-two Russell compiler

```text
D=1+x^2q,
b=(D-1)(D+2)^2,
c=xD(D+2),
e=q(D+3),                                               (5)
```

and the two simple-zero stable displacements

```text
H_0(t)=t,                    H_1(t)=t+t^2.              (6)
```

For either displacement set `q=Q_alpha(x)+H_i(t)` and `w=t`.  At `t=0`
the three points in `(3)` all map to

```text
(b,c,e,w)=(0,0,-3,0).                                  (7)
```

The regular target parameters used for the arbitrary-form packet are

```text
y=c/3,                    z=e+3,                    w. (8)
```

For `p in {-1,0,1}`, write `x=p+xi`.  Let

```text
P_N=direct_sum_(p=-1,0,1) K[xi,t]/(xi,t)^(N+1),        (9)
```

where a target two-form is recorded by its three pulled density classes
relative to `d xi wedge dt`.  Thus

```text
dim_K P_N=3 binom(N+2,2),             dim_K P_5=63.    (10)
```

For

```text
Omega=A dy wedge dz+B dy wedge dw+C dz wedge dw,       (11)
```

define `T_(H,N)(Omega)` to be this three-branch retained packet.  Because
`y,z,w` all vanish at the common target point, coefficient monomials of
degree greater than `N` are invisible in `P_N`.  Hence the complete
arbitrary-form domain at cutoff `N` is supplied by

```text
A,B,C in K[y,z,w]_(degree<=N).                          (12)
```

At `N=5`, `(12)` has

```text
3 binom(8,3)=168                                        (13)
```

columns.  This is a complete local target two-form jet universe, not a
sampled ansatz.

## 2. The evenness repair and the parity-free local Darboux pair

The hypothesis `Q` even in THM-3629 is load-bearing for its shifted `D=0`
sign collision and hence for its conclusion that, for nonlinear `H`, both
members of every hypothetical global pair must mix the surface and stable
coordinates.  Equation `(4)` shows that this conclusion cannot be imported
for `Q_alpha` and `H=t+t^2`.

The local affine pair survives for a different, parity-free reason.  On the
chart `b+4!=0`, which contains the common target point `(7)`, put

```text
a=e/(b+4).                                               (14)
```

Directly from `(5)`,

```text
b+4=D^2(D+3),             a=q/D^2,
Jac_(x,q)(a,c)=-3.                                      (15)
```

These identities make no parity assumption on `Q`.  Along
`q=Q_alpha(x)+H(t)` they give

```text
Jac_(x,t)(a,c)=-3H'(t).                                (16)
```

Therefore the affine displacement has the exact local Darboux pair

```text
F=a,                    G=-4c,
Jac_(x,t)(F,G)=12                    when H=t.          (17)
```

The associated two-form extends globally even though `a` does not.  Namely,
on the surface,

```text
eta=da wedge dc
   =(e/2)db wedge dc+(3c/8)db wedge de
                    -(b+2)/8 dc wedge de               (18)
```

is a global polynomial exact form.  An explicit primitive is

```text
alpha_surf=-(ce/4)db+(e(b+2)/4)dc+(c(b+2)/8)de,
d alpha_surf=eta.                                      (18a)
```

The compiler pulls `eta` back to `-3 dx wedge dq`.  Thus `-4eta` pulls back
to the density `12` for `H=t`.
This parity-free fragment is proved by `(14)--(18)` directly; it is not an
application of the even-`Q` theorem as a whole.

The distinction is essential: `(17)` is a pair only on the local chart,
while `(18)` is a global two-form with no asserted global polynomial
decomposition `dF wedge dG`.

## 3. Complete cutoff-five polynomial-exact saturation

Let `E_N` be the target two-forms

```text
d beta,       beta=U dy+V dz+W dw,
U,V,W in K[y,z,w]_(degree<=N+1).                       (19)
```

Potential monomials depending only on the coordinate of their differential
have zero derivative.  After these zero columns are removed, `(19)` gives

```text
231 nonzero exact-form generators at N=5.              (20)
```

Every element of `E_N` belongs to the arbitrary universe `(12)`, and
potential terms of degree greater than `N+1` cannot affect `P_N`.  Exact
elimination over the quartic field `(1)` gives the following complete table:

| displacement | rows | arbitrary columns | exact columns | arbitrary rank | exact rank | constant in both images |
|---|---:|---:|---:|---:|---:|---:|
| `H=t` | 63 | 168 | 231 | 59 | 59 | yes |
| `H=t+t^2` | 63 | 168 | 231 | 59 | 59 | yes |

Since the exact-form image is contained in the arbitrary-form image and the
ranks agree, this proves the equality

```text
T_(H,5)(E_5)=T_(H,5)(all target two-form jets)          (21)
```

for each `H` in `(6)`.  Both sides have codimension four in `P_5`.  In the row order by `t`-degree, then `xi`-degree, then
`p=-1,0,1`, the four retained relations have
nonzero support counts

```text
(3,7,14,23).                                            (22)
```

Their canonical serializations have SHA-256

```text
H=t:       bb0a2ba7f46b6bf721c875a1a7a19000e1b3b8845ba11bc4707bd299b9740852,
H=t+t^2:   bbaf81804c69200a9214711f4f3865acb39831f22af02d623436e788dd6b633e. (23)
```

The two displacement-specific images in `(21)` are not the same subspace of
`P_5`: the span of their two four-dimensional relation spaces has dimension
six.  Thus polynomial-exact saturation holds separately on each source, but
does not identify the affine and quadratic-displacement packets.

All four relations are zero on the terminal `t^5` value block.  If
`mathbf 1` denotes the packet whose constant coefficient is one on each
branch and whose other coefficients vanish, then

```text
mathbf 1 in im T_(H,5)|_(E_5)                          (24)
```

for both displacements.  Thus a polynomial exact target two-form realizes a
common nonzero constant through total source degree five.  This is an exact-form image statement; no pair decomposition is inferred.

The equality `(21)` concerns a finite retained pullback image.  It does not
identify the exact-form bank with the image of the nonlinear map
`(F,G)|->dF wedge dG` in global target functions.

## 4. The fixed-`a` lane fails first at cutoff five for `H=t+t^2`

Use the vanishing local coordinates

```text
A=a+3/4,                    c,                    w.   (25)
```

For fixed first output `F=a`, define

```text
D_(a,H,N): K[A,c,w]_(degree<=N+1)/K -> P_N,
G |-> retained density of da wedge dG.                 (26)
```

This domain is complete at cutoff `N`, since `A,c,w` vanish at `(7)` and a
term of degree greater than `N+1` cannot enter the retained density.  Along
the source, the exact formula behind `(26)` is

```text
Jac_(x,t)(a,G(A,c,w))=-3H'(t)G_c+a_xG_w.               (27)
```

For `H=t`, the choice `G=-4c` in `(17)` solves the problem exactly.  For
`H=t+t^2`, exact quartic-field elimination gives

| cutoff | packet rows | nonzero `G` columns | rank | relation dimension | constant in image |
|---:|---:|---:|---:|---:|---:|
| 4 | 45 | 55 | 40 | 5 | yes |
| 5 | 63 | 83 | 57 | 6 | no |

A cutoff-four solution truncates to every lower cutoff, while completeness
of `(26)` prevents any higher target Taylor term from changing cutoff five.
Thus five is the first failure order in this fixed-`a`, three-multigerm lane.

In the canonical RREF gauge for the row order used in `(23)`, exactly one of
the six displayed cutoff-five relation-basis vectors has nonzero response on
`mathbf 1`.  That response is

```text
rho=
 -2073506706944/1678822119
 +(372679949312/62178597) alpha
 -(184159683584/6908733) alpha^2
 -(73442787328/2302911) alpha^3.                       (28)
```

Its field norm is

```text
N_(K/Q)(rho)=
 28278059768285603108255604733711150481408
 /392425657272606224710564875 !=0.                    (29)
```

Consequently the obstruction survives every complex embedding of `(1)`.
Equations `(28)--(29)` obstruct only a mate for the fixed local first output
`a`; they do not obstruct arbitrary mixed pairs.

## 5. The full mixed tangent bank absorbs the quadratic perturbation

Work over the dual numbers `K[epsilon]/(epsilon^2)` and set

```text
H_epsilon=t+epsilon t^2,
F_epsilon=a+epsilon f,
G_epsilon=-4c+epsilon g,                               (30)
```

where `f,g` are polynomial germs in `A,c,w`.  On the affine base source
`H=t`, the first-order pair correction is

```text
L(f,g)=12f_A+4c_xf_w-3g_c+a_xg_w.                     (31)
```

The uncorrected change of `H'` contributes `24t`, so preserving density `12`
through first order requires

```text
L(f,g)=-24t.                                           (32)
```

At cutoff five the complete degree-at-most-six bank for `(f,g)` has

```text
166 columns, rank 59, relation dimension 4.            (33)
```

Its image equals the complete arbitrary retained image for the **affine base
source `H=t`** in Section 3, and the packet `-24t` belongs to that image.
It is not the distinct arbitrary image for `H=t+t^2`.  Hence polynomial
local germs `f,g` exist for which `(30)` has constant density through cutoff
five modulo `epsilon^2`.

This is a tangent-space statement.  It neither integrates the deformation
to `epsilon=1` nor supplies corrections at order `epsilon^2`, at source
cutoff six, or at any later order.  It also does not prove that both outputs
of a hypothetical global pair must mix; that unavailable conclusion was the
evenness error repaired in Section 2.

## 6. Connection contract

The exact connection proved here is the following.

| field | contract |
|---|---|
| source | arbitrary two-form jets in the regular parameters `(y,z,w)`; polynomial exact forms `d beta`; and, separately, pair tangents `(f,g)` about `(a,-4c)` |
| target | the 63-dimensional three-branch packet `P_5` of total `(xi,t)` degree at most five |
| map | pull back through `(x,t)->(x,Q_alpha(x)+H(t),t)`, then retain the three Taylor densities |
| preserved predicate | equality of source densities through cutoff five; exactness of a chosen target form; the common constant-density target |
| destroyed information | target terms above the cutoff, coherence between cutoffs, convergence, global regularity of `a`, global decomposability as `dF wedge dG`, injectivity/properness, and polynomial degree at infinity |
| needed sidecar | a coherent all-order pair lift plus a global denominator/decomposition argument; for a Keller consequence, global polynomial functions and the usual noninjective-collision verification |
| cheapest decisive tests | exact rank equality for arbitrary versus exact images; augmentation by `mathbf 1`; the fixed-`a` left response `(28)` and norm `(29)`; augmentation of the tangent image by `-24t` |

This contract prevents two tempting but false upgrades: exact two-form
saturation is not global Darboux decomposition, and a full finite tangent
image is not nonlinear integrability.

## 7. Inheritance and hostile controls

The inheritance pass is as follows.

- **Closest proved mechanism.**  The parity-free algebraic fragment of
  THM-3629 supplies `(14)--(18)`.  THM-3631 is the closest finite/local
  analogue: for a different non-even hostile, polynomial closed forms and
  local pair germs survive through cutoff five before an arbitrary-form
  failure at six.
- **Canonical hostile.**  THM-4046 proves that on this same exceptional
  quartic, every nonzero critical displacement `H in t^2 C[t]` has an
  arbitrary-form constant-pullback obstruction at retained order eight.
  The present simple-zero displacements have `H'(0)=1` and are disjoint from
  that theorem's scope.
- **Corrected near miss.**  Applying THM-3629's nonlinear mixing conclusion
  to `Q_alpha` was invalid because `(4)` contradicts its evenness hypothesis.
  Only the directly rechecked parity-free identities survive.
- **Least-used sidecar.**  THM-3737's actual restriction image is the
  hyperplane `ker Lambda` and controls stagewise actual target coefficients.
  The present domain is instead the complete local two-form jet bank.  The
  two images must not be identified without an actual-representative and
  globality sidecar.

Useful hostile controls for the companions are:

1. `H=t` with the exact pair `(a,-4c)`;
2. `H=t+t^2` at cutoff four, where fixed `a` still succeeds;
3. the same fold at cutoff five, where arbitrary and exact forms succeed but
   fixed `a` fails;
4. all four split roots `44,82,92,134` modulo `137`;
5. the dual-number target `-24t`, which the tangent bank must contain.

## 8. Reproduction and artifact gate

The frozen characteristic-zero production companion and independently
implemented modular audit satisfy the following gates.

The production companion, in its default run:

1. construct `K` and `Q_alpha` exactly and print `(3)--(4)`;
2. build the complete arbitrary, polynomial-exact, fixed-`a`, and mixed
   tangent universes without importing stored matrices;
3. run both `H=t` and `H=t+t^2` at cutoff five and the sharp fixed-`a`
   positive control at cutoff four;
4. compute the exact `K`-ranks directly, rather than infer the
   polynomial-exact rank from a modular sandwich;
5. verify image containment, constant membership/nonmembership, terminal
   projection, relation residuals, `(28)`, and `(29)`;
6. print the canonical relation hashes; and
7. contain no Python `assert` statement.

The independent audit must reconstruct the local arrays without importing
production code, check the local identities in `(15)`, verify a good split
prime and all four exceptional roots, and reproduce the decisive rank and
membership table.  Normal and optimized runs of both companions must be
byte-identical to their frozen LF transcripts.

Required final commands are

```bash
python3 -B 04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054.py
python3 -O -B 04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054.py
python3 -B 04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054_independent_audit.py
python3 -O -B 04-computation/jc2_exceptional_affine_simple_zero_retained_packet_thm4054_independent_audit.py
```

The four raw-LF SHA-256 values are recorded in the frontmatter.  Both
normal/optimized pairs are byte-identical to the frozen outputs.

## 9. Strict boundary

This theorem proves only statements about the three completed source germs
at `x=-1,0,1` through total source degree five, plus the exact affine local
pair `(17)` and the dual-number tangent statement `(30)--(33)`.  It does not
prove any of the following:

- a global regular or polynomial decomposition of the exact forms in
  Section 3 as `dF wedge dG`;
- a global polynomial mate for the local rational coordinate `a`;
- failure of arbitrary mixed pairs for `H=t+t^2`;
- integration of the first-order correction or a coherent all-order lift;
- convergence, algebraization, or polynomial degree control;
- the even-`Q` shifted-sign collision or nonlinear mixing conclusion of
  THM-3629;
- a new statement for critical displacements `H in t^2 C[t]`, already in the
  separate scope of THM-4046;
- an arbitrary Darboux pair on the Danielewski surface, whose global support
  problem remains open beyond the proved cells including THM-3592; or
- a global Keller map or a counterexample to `JC(2)`, which remains **OPEN**.

The positive result is exact-form saturation of one finite simple-zero
packet.  The negative result is a sharp obstruction in one fixed-coordinate
lane.  The tangent escape explains why the latter does not close the full
mixed problem.  **QED.**
