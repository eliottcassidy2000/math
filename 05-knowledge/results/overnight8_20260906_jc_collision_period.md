# Collision-preserving ninth source jets cannot change the retained wedge period

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED.** The
[independent proof and exact audit](overnight8_20260906_jc_collision_period_audit.md)
passes without a mathematical correction. This
strengthens the [seventh fixed-prefix source theorem](overnight7_20260906_jc_fifth_module.md)
by one source order when the labelled collision is preserved. It applies to
arbitrary smooth-target map perturbations in the stated fixed source
coordinates, hence in particular to actual polynomial source constructions.
It does not settle unrestricted source charts, formal unit densities,
polynomial termination, or planar JC(2).

## 1. Inheritance and the apparent escape

The closest proved mechanism is
[THM-4046 / exceptional-quartic-j7-lift-and-j8-obstruction](../../01-canon/theorems/THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction.md):
every arbitrary target two-form pulled back by the exceptional quadratic
source has a complete retained relation

```text
0=Lambda(J8)+R4(J6)+R2(J4)+R0(J2)+S0(J0),
S0(1)=kappa!=0,                  Lambda=L/18,
L(v)=5v_- -18v_0+13v_+.
```

The seventh theorem proves that every source-map change O(t10) preserves
all these entries literally. Its canonical hostile is t9*x: this changes
J8 by `36(-1,0,1)` for `dy wedge dz`, with Lambda value16. It shows
the unrestricted order estimate is sharp. It does not preserve the
three-branch collision at ninth order.

The corrected near miss is to treat an exposed earlier coefficient as an
admissible escape. The least-used coordinate is the common target motion
at that same order. Tangent reparametrizations disappear inside a wedge,
and the weighted tangent relation kills the remaining common motion.

The live board is: fixed source clocks; simultaneous branch contact;
tangent-relation cokernels; actual target two-forms; source reparametrization;
and low-order versus late correction. Targeted searches of THM-4046,
THM-4424, the seventh source report and the correction ledger retain the
distinction between THM-4046's ninth **target-pair** coefficients and the
ninth **source-map** perturbations studied here.

## 2. A general leading wedge-period lemma

Let K have characteristic zero, let Y be a smooth target, and use fixed
source coordinates `(x,t)`. Let f and ftilde be formal maps into one
completed smooth target neighborhood. Suppose N>=2 and

```text
ftilde(x,t)-f(x,t)=t^N v(x)+O(t^(N+1)).                    (1)
```

Let distinct source points x_i map at t=0 to the same target point y0.
Put `T_i=partial_x f(x_i,0)`. Suppose all T_i are nonzero and at least
two of their one-dimensional spans are distinct. Suppose sections z_i(t)
through x_i have common f-image modulo t^(N+1), and sections ztilde_i(t)
through those same labelled x_i have common ftilde-image modulo t^(N+1).
Let constants ell_i satisfy

```text
sum_i ell_i T_i=0.                                       (2)
```

For any regular target two-form omega, let J_r and Jtilde_r be the
coefficients of the corresponding `dx wedge dt` pullback densities.
Then

```text
Jtilde_r=J_r as functions of x, r<N-1,
sum_i ell_i (Jtilde_(N-1)(x_i)-J_(N-1)(x_i))=0.            (3)
```

The second conclusion is a weighted retained equality, not pointwise
equality of the last coefficient.

**Proof.** The two families of sections agree modulo t^N. Induct on their
first potentially different order j<N. Equation (1) contributes nothing
at that order, and subtraction of their common-image equations gives
`T_i r_i=b` for one common target vector b. The intersection of two
distinct tangent lines is zero, so b=0, and all r_i=0 because T_i!=0.
This proves the induction, including equality of the common image prefix.

At order N the difference instead gives

```text
v(x_i)+r_i T_i=b,                                        (4)
```

for one common vector b. In fixed smooth coordinates, differentiating (1)
in x preserves its order; differentiating in t gives
`N t^(N-1) v+O(t^N)`. Coefficients of omega change only at order N.
Therefore the first possible density difference is exactly

```text
Jtilde_(N-1)(x_i)-J_(N-1)(x_i)
   =N omega_(y0)(T_i,v(x_i)).                              (5)
```

The term with the changed x derivative is order N, so it cannot enter
(5), regardless of the old t derivative. By (4), alternation kills r_i T_i.
Summing (5) with (2) gives
`N omega_(y0)(sum ell_i T_i,b)=0`, proving (3). This proof is invariant
under a change of smooth target coordinates. It requires a single coherent
target germ, not separately assigned branch forms. No closedness or
decomposability of omega is assumed. This proves the lemma.

## 3. Apply the lemma to the exceptional triple

Use the characteristic-zero exceptional field and polynomial Q of
[THM-4424 / russell-constant-normal-debt-discriminant-contact-correspondence](../../01-canon/theorems/THM-4424-russell-constant-normal-debt-discriminant-contact-correspondence.md).
Let Phi be the actual Russell compiler

```text
D=1+x^2q, C=xD(D+2), E=q(D+3),
f0(x,t)=(y=C(x,Q+t^2)/3,z=E(x,Q+t^2)+3,w=t).
```

The local target is smooth at the common image of x=-1,0,1, and `(y,z,w)`
are regular coordinates there. The three source tangents are

```text
T_-=(1,-9,0), T_0=(1,4,0), T_+=(1,9,0),
5T_- -18T_0+13T_+=0.                                    (6)
```

The uncorrected constant pencil's branch sections through fourth s-order
have common image modulo s5=t10. Thus N=9 meets precisely the old
collision premise in the lemma. Let ftilde be **any** map into this smooth
target with

```text
ftilde=f0+O(t9),
its three labelled branches have common image modulo t10. (7)
```

The lemma keeps every J0 through J7 identical, including all retained
x-jets. It keeps Lambda(J8) identical by (6), which is the only retained
order-eight combination used in THM-4046. Hence the complete old relation
remains valid for every target two-form on every map (7). A nonzero
constant source density c would again give `c*kappa=0`, impossible.

**Conclusion.** No arbitrary regular target two-form has nonzero constant
pullback density on any map (7). In particular no pair of actual target
functions yields a nonzero constant source Jacobian there. This includes
polynomial source planes and more general actual constructions meeting (7);
the proof does not claim every arbitrary map in that class is supplied by
a polynomial source automorphism. It is stronger than an exclusion limited
to target-descended Hamiltonian variations.

## 4. Exact source conditions and sharp hostiles

For graph perturbations

```text
q=Q(x)+t2+t9 h(x)+O(t10), w=t,                            (8)
```

the ninth collision equation is soluble if and only if `L(h)=0`.
Indeed the lower coefficients are fixed, the five section/common-surface
unknowns have rank five, and the sole residual is the tangent-relation
cokernel. At each base point, writing `e_i=(-9,4,9)`, one has

```text
T_i=(1,e_i,0),
U_i=partial_q(y,z,w)=(2/3,-2,0),(0,4,0),(-2/3,-2,0),
det_(y,z)(T_i,U_i)=4.
```

The cokernel applied to the forcing U_i h_i is `4L(h)`. Its vanishing
is also sufficient because the rank-five system has only that one
cokernel coordinate. Higher terms in (8) cannot affect this order.

For `omega=A dy wedge dz+B dy wedge dw+Cform dz wedge dw`, let A0,B0,C0
be its three coefficients at the common base point. Direct differentiation
gives

```text
delta J8(x_i)=36 A0 h(x_i),
delta Lambda(J8)=2 A0 L(h).                               (9)
```

The other two form slots first change at order nine. Thus every legal
ninth collision coefficient kills the apparent period payment. The
polynomials 1, `4x^2-9x`, and `x(x^2-1)` are explicit kernel controls;
`h=x` has L=8 and is the precise hostile to dropping collision preservation.
It changes the retained matrix but cannot serve as a legal ninth repair.

The same condition can be read before choosing a graph. For an actual
ambient source perturbation beginning

```text
(x,q,w)=(x+t9 r(x), Q(x)+t2+t9 s(x), t+t9 k(x))+O(t10),
h_i=s(x_i)-Q'(x_i)r(x_i),
```

its target forcing is `r_i T_i+h_i U_i+k_i e_w`. The full ninth
collision equation is soluble exactly when

```text
L(h)=0,                 k_-=k_0=k_+.                     (10)
```

There are three independent cokernel conditions: one surface relation
and two stable-coordinate differences. Tangential r_i is absorbed by
section changes. The first period variation is

```text
2 A0 L(h) + (B0 L(k)+C0 L(e*k))/2.                        (11)
```

Equation (10) makes (11) zero. This is a criterion on coefficients of
actual source changes, not a claim that arbitrary triples r,s,k define
globally invertible polynomial automorphisms. The graph shear subset is
actual, and already supplies nontrivial admissible examples.

There is an exact converse for this retained triple: vanishing of (11)
for **every** base two-form is equivalent to (10). Its A coefficient
forces L(h)=0. The B and C coefficients force L(k)=L(e*k)=0, whose joint
kernel is precisely the constant line: the two rows are independent and
both kill `(1,1,1)`. Thus the complete three-coordinate wedge-period
response detects every ninth collision incompatibility. One selected
two-form is insufficient. For h=x and k=x the period is
`16A0+4B0+81C0`; choosing `(A0,B0,C0)=(1,-4,0)` hides the collision
split even though L(h)=8 and the stable values are unequal.

## 5. Scope and procedural consequence

The source-to-target map is differentiated pullback followed by the
tangent-relation covector. It preserves the old retained obstruction
because common branch motion survives and tangential section motion is
annihilated. Forgetting collision congruence permits h=x and loses this
conclusion. Forgetting the fixed source density would wrongly exclude
formal units; that claim is not made.

For an even source graph, every change O(t9) is already O(t10), so the
seventh theorem suffices. The new content is the odd ninth jet and arbitrary
target-map perturbations preserving the labelled contact. Any escape from
this combined fixed-prefix class must change a source coefficient at
order at most eight, or give up the specified ninth collision congruence.
Neither possibility is ruled out here.

This is not a consequence of the old target-pair ninth correction alone:
the source map now changes. The new preserved quantity is the leading
wedge period under simultaneous collision motion. The next decisive
object is the earlier source-jet space with both its collision conditions
and its complete retained two-form response, not a single changed J8 entry.

## 6. Reproduction and independent audit

The primary imports no repository mathematics. It reconstructs the actual
compiler tangent and normal vectors, the complete nine-equation collision
matrix, the universal wedge identity, both directions of the retained
collision criterion, and the one-form hostile. Fifteen literal local
compiler perturbations verify all three arbitrary two-form slots, every
unchanged J0..J7 coefficient, and the exact J8 variation. All **59** gates
pass normally and optimized with identical LF output.

The separate referee uses fresh exceptional embeddings `(449,120)` and
`(467,169)`, overdetermined section solving, and a dual-number compiler.
It checks actual collision congruences, graph/tangential/stable controls,
and 864 literal density/form cases including nonconstant target coefficients.
Its **1,871** gates pass normally and optimized with identical output.
It audits the general lemma and all quantifiers analytically; the finite
controls are not its proof. The parallel-tangent hostile in that audit
explains the distinct-tangent hypothesis in the section-uniqueness step.

```text
python -B 04-computation/overnight8_20260906_jc_collision_period.py
python -B -O 04-computation/overnight8_20260906_jc_collision_period.py
python -B 04-computation/overnight8_20260906_jc_collision_period_audit.py
python -B -O 04-computation/overnight8_20260906_jc_collision_period_audit.py

primary source a7295b4fe6b4053b5cd71678acd0bc73d8df53107bba0bc699c5c78b0d4cac42
primary output b267f206667e1bc7104a864eabf3c4f02dd1828e1c9f2919096180bcfd755d47
audit source 01600a2c46790da5716a71017202e05181add970fdab598eae76e354c0bc96c4
audit output 24893fd9f63d3788dab4ed677e8846a0593a423773fee2d324a8fd24a05fec18
```

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
