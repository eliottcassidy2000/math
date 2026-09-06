# The fifth collision coefficient is source-admissible, but misses the target-form obstruction

**Status: PROVED in the polynomial-source and filtered target-Hamiltonian
categories defined below; INDEPENDENTLY AUDITED.**
[Independent proof and exact audit: PASS](overnight7_20260906_jc_independent_audit.md). This resolves the
maintained fifth-order test after separating three distinct admissibility
questions. It supplies neither a Keller pair nor a JC(2) result.

1. The actual polynomial source-plane coefficient space is all `K[x]`.
   A triangular polynomial source automorphism realizes the required fifth
   coefficient with every lower coefficient fixed.
2. For the full space of target-descended Hamiltonian generators, including
   arbitrary cancellations among their lower pullback jets, the retained
   normal-value image is exactly `ker(5,-18,13)`. It cannot pay that debt.
3. Regardless of descent, **every** source perturbation starting at `t^10`
   leaves the known `J_8` constant-Jacobian obstruction unchanged. The fifth
   `s=t^2` collision correction arrives one density order too late.

## 1. Inheritance and the question made precise

The current frontier and active guardrails keep JC(2), chart entry, and
termination open. The closest proved result is
`THM-4424-russell-constant-normal-debt-discriminant-contact-correspondence.md`:
the exceptional constant pencil has a unique transverse correction
`chi5=-9*kappa/20!=0`, and exact contact order five. The earlier
`synthesis_20260905_transgression.md` supplies the determinant `-288`
formal lifting system. `THM-4411-first-order-collision-transgression-seminormal-tradeoff.md`
supplies its normal functional

```text
L(h)=5h(-1)-18h(0)+13h(1),                 L(x)=8.
```

The canonical hostile is `h=x`: it pays the seminormal response by splitting
the original triple at first order. The corrected near miss here is to
equate formal collision repair with payment of an actual target-form debt.
The least-used sidecars are the old **polynomial source-plane definition**
and the order lost when differentiating a source-map jet.

The lawful source charts are recovered from
`THM-3612-russell-cylinder-even-fold-nongraph-collision-jet-rigidity.md`,
equations (4)-(6), rather than inferred from ambient coefficient freedom.
The target category is the actual Russell cylinder ring of
`THM-3605-russell-cylinder-graph-slice-puncture-no-filling.md` and
`THM-3629-russell-cylinder-linear-vertical-fold-global-form-boundary.md`.
The retained target-form obstruction is the complete, arbitrary-form
`THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction.md`.

No claim identifies arbitrary polynomial source embeddings with
target-descended Hamiltonian moves. They are explicitly different spaces.

## 2. Fixed objects and clocks

Work over the characteristic-zero exceptional field

```text
K=Q[alpha]/(72783360alpha^4-77822208alpha^3-28419741alpha^2
                                  +7849770alpha-1276420).
```

Let Q be the exceptional polynomial of THM-4424, with the Russell compiler

```text
D=1+x^2q,  B=(D-1)(D+2)^2,  C=xD(D+2),  E=q(D+3),
C^2 E=B(B+4).
```

Write `Phi(x,q)=(B,C,E)` and `q0(x,t)=Q(x)+t^2`. The stable target
coordinate is `w=t`. At `x=-1,0,1`, the undeformed surface image is
`(0,0,-3)`. Local surface coordinates are `(C,Z=E+3)`; `(y=C/3,z=E+3,w)`
are the normalized coordinates used by the retained J-matrices.

The ordinary tangent rows in `(C,E)` are

```text
t_-=(3,-9), t_0=(3,4), t_+=(3,9),
5t_- -18t_0+13t_+=0.                                           (1)
```

Keep the two clocks distinct: the collision pencil uses `s`, while the
actual quadratic source fold uses `s=t^2`. A fifth s-coefficient is a
tenth t-coefficient of the source map.

THM-4424 gives, with kappa nonzero at all four exceptional embeddings,

```text
chi5=-9*kappa/20,
q_s=Q+s+chi5*s^5*x+O(s^6).                                     (2)
```

For a prescribed `q_s=Q+s+s^5 h(x)+O(s^6)`, the same rank-five repair
system at the first unpaid coefficient gives

```text
labelled triple agrees modulo s^6  iff
L(h)=8chi5=-18*kappa/5.                                        (3)
```

All lower section and target coefficients are already determined. At this
degree, subtract the known solution (2); the remaining condition is exactly
`L(h-chi5*x)=0`. This proves (3), with no linearization of the lower debt.

## 3. Actual polynomial source construction

For every `h in K[x]`, define

```text
E_h(x,t)=(x,Q(x)+t^2+t^10 h(x),t).                             (4)
```

Its image is the closed polynomial plane

```text
q-Q(x)-w^2-w^10 h(x)=0,
K[x,q,w]/(q-Q(x)-w^2-w^10 h(x)) = K[x,w].                       (5)
```

Thus these are genuine source `A^2` charts of the old non-graph type, not
merely coefficient vectors. More strongly, the ambient polynomial shear

```text
T_h(x,q,w)=(x,q+w^10 h(x),w),                                 (6)
```

has determinant one and polynomial inverse obtained by subtraction.
It takes the old source plane to (5), fixes x and w, and fixes every source
coefficient below t-order ten. With relative source form `dx wedge dq`,
it is the time-one shear generated by the polynomial Hamiltonian
`-w^10 integral h(x) dx`. In characteristic zero every h has that polynomial
primitive. This is a **source** primitive; it is not asserted to descend
through the compiler.

Consequently the actual polynomial source normal space at this order is
`K[x]`, and `L` maps it onto K. Taking `h=chi5*x` realizes (3) by an actual
polynomial plane and source automorphism. The three labelled sections have
equal complete target images modulo `t^12`, since (3) is modulo `s^6`
and all three stable coordinates are t. This makes no all-order claim for
the truncated polynomial (4); higher collision coefficients require their
own corrections.

## 4. The full filtered target-Hamiltonian response

Let `R=K[B,C,E]/(C^2E-B(B+4))`. For an **arbitrary actual target polynomial**
`F in R[w]`, put

```text
H(x,q,w)=F(Phi(x,q),w),
A_F(x,t)=H(x,q0(x,t),t).
```

A Hamiltonian first variation with generator H sends `(x,q)` with velocity
`(H_q,-H_x)`. Rewriting the moved graph at fixed x gives the effective
normal velocity

```text
-H_x-q0_x H_q = -partial_x A_F.                                (7)
```

The residue form pulled back from the Russell surface is `3 dx wedge dq`;
rescaling target Hamiltonians by the nonzero constant three gives exactly
the normalization in (7). Hence this is an actual descended Hamiltonian
map, not a guessed coefficient rule.

Define the K-linear space

```text
H_10={F in R[w]: partial_x A_F is divisible by t^10},
M_10={-[t^10]partial_x A_F : F in H_10} subset K[x].              (8)
```

It contains **all** lower-jet cancellations in F. It is not defined by
requiring each target monomial to be individually divisible by w^10.
Nor is M_10 claimed to be a K[x]-submodule; multiplication by x need not
preserve its defining condition or its retained image.

**Theorem (exact retained image).**

```text
ev_(-1,0,1)(M_10)
 = span_K{(3,3,3),(-9,4,9)}
 = ker(5,-18,13).                                               (9)
```

In particular `L(M_10)=0`, so the affine debt hyperplane (3) has empty
intersection with this full filtered descended normal space.
The leading primitive has an exact stronger six-coordinate description:
its three values are equal, and its three derivatives lie in the tangent
plane in (9). Every such value/derivative packet occurs. This classifies
the retained first jets, not the entire global polynomial space M_10.

Here is an all-degree proof of the upper bound, including cancelled target
prefixes. From (8), the coefficients of A_F below t-order ten are constant
in x. Subtract their polynomial in the stable target coordinate w. This
changes neither (7) nor membership in (8), and gives

```text
A_F(x,t)=t^10 a(x)+O(t^11).                                    (10)
```

Write F in completed local target coordinates `(C,Z,w)`. Its w=0 term F0
vanishes on each of the three base-curve germs, by (10). Differentiating
the three restrictions once gives `dF0(t_i)=0`; two independent tangent
rows force `dF0=0` at the triple. Differentiating twice then gives
`Hess(F0)(t_i,t_i)=0`. The three directions in (1) force the complete
Hessian to vanish: in the symmetric-coordinate basis its matrix is

```text
[9,-54,81; 9,24,16; 9,54,81],        determinant=63180!=0.        (11)
```

Choose the old uncorrected branch sections from THM-4424 through fourth
s-order. Their surface images `gamma_i(t)` agree modulo t^10. Since the
constant term of the surface Hessian of F is zero by (11), Taylor
expansion shows that their pulled-back surface gradients agree one order
further:

```text
grad_(C,Z)F(gamma_i(t),t)
 - grad_(C,Z)F(gamma_j(t),t) = O(t^11).                          (12)
```

Indeed the surface displacement is O(t^10), its Hessian coefficient has
zero constant term, and quadratic displacement terms are still higher.
The same Taylor argument using `dF0=0` gives
`F(gamma_i(t),t)-F(gamma_j(t),t)=O(t^11)`, so (10) forces the three
values of a to be equal.
Let `t_i(t)=partial_x gamma(x,t)|_(x=z_i(t))`. There is a formal relation
`sum ell_i(t)t_i(t)=0`, normalized with `ell(0)=(5,-18,13)`; a nonzero
two-tangent minor is a unit. Using (12) in the chain rule gives

```text
sum ell_i(t) partial_x A_F(z_i(t),t) = O(t^11).
```

By (10), its leading coefficient is
`5a'(-1)-18a'(0)+13a'(1)=0`, proving the upper bound in (9).
The proof also applies to one coherent completed target germ whose three
pullbacks have the normalized prefix (10). It assumes no degree bound on F.

The lower bound is actual and polynomial: `F=w^10 C` and
`F=w^10(E+3)` lie in H_10 and give the negatives of the two displayed
tangent-coordinate rows. Those rows span the two-dimensional kernel in
(9). Thus equality is proved, not merely a necessary period condition.
Adding `F=w^10` attains the common primitive-value coordinate independently,
proving the asserted complete three-dimensional primitive first-jet image.

This is a general local move: at an ordinary triple point, vanishing of a
descended zeroth coefficient kills its gradient and Hessian. An Nth-order
collision congruence then gives an `(N+1)`st-order gradient congruence,
which forces the leading derivative packet into the tangent-relation
kernel. N=10 is the present application.

## 5. A sharper obstruction independent of Hamiltonian descent

**Filtered pullback lemma.** If two source maps into a smooth target agree
modulo `t^N`, then their pullbacks of any regular target two-form have
`dx wedge dt` densities agreeing modulo `t^(N-1)`.

Proof: target coefficient functions still agree modulo t^N. Differentiating
in x preserves that order, while differentiating in t lowers it by at most
one. Expand each coordinate wedge and its coefficient; every difference
has order at least N-1. This holds equally in completed regular coordinates
and for polynomial maps. It uses no parity assumption on the target form.

Apply it with N=10 to (4), or to **any** graph

```text
q(x,t)=Q(x)+t^2+O(t^10),              w=t.                       (13)
```

The complete densities `J_0,...,J_8`, including all their retained x-jets,
are unchanged for every target form. The THM-4046 relation therefore
remains literally the same:

```text
0=Lambda(J8)+R4(J6)+R2(J4)+R0(J2)+S0(J0),
S0(1)=kappa!=0.                                                (14)
```

A nonzero constant source density would have `J0=c!=0` and `J1=...=J8=0`,
so (14) would read `0=c*kappa`, impossible. Hence **no arbitrary target
two-form**, and in particular no actual target pair, has nonzero constant
source Jacobian on any graph (13). This excludes every fifth-or-later
source repair in this frozen-prefix class, whether its primitive descends
or not. It also applies to a fully formal compensated graph, for this
fixed source density and these target germs.

Do not replace "nonzero constant" by "formal unit": densities such as
`1+t^8 a(x)+...` can be units in a completed ring and are not excluded.
An arbitrary formal source reparametrization can change the density by a
nonconstant unit; no such change is silently treated as a polynomial
constant-Jacobian coordinate change.

The order bound is sharp. Direct computation gives

```text
Jac_(x,q)(C,E)=6(D+1)(D^2+2D-2),
Jac_(x,q)(C,E)|_(retained points)=12.
```

For the normalized target form `d(C/3) wedge d(E+3)`, a perturbation
`t^10 h(x)` first changes the retained density at J9, by `40h(i)`.
For the actual payer `h=chi5*x`, its Lambda response there is `-8*kappa`.
It cannot reach J8. By contrast a `t^9 x` perturbation changes J8 by the
retained vector `36(-1,0,1)`, whose Lambda value is 16. This only shows
that the earlier-order matrix can change; it constructs no constant form.
If one retains even source dependence, an escape must change at least the
fourth s-jet; allowing odd dependence first exposes the ninth t-jet.

## 6. Exact controls and remaining scope

The standalone companion imports no repository mathematical code. Symbolic
identities check the field element, source shear, gradient/Hessian linear
systems, and the full compiler wedge. Independent finite-field local
series at `(p,alpha)=(421,126),(443,112)` solve all six collision equations
through s^5 for zero, paying, constant-kernel, quadratic-kernel, and unpaid
constant controls. They recover `(kappa,chi5)=(180,340),(361,347)`.

A separate complete retained primitive bank assigns weight two to a local
normalization displacement and weight one to t. All target monomials
`C^a Z^b w^c` with `2a+2b+c<=12` give 140 columns. The 135 constraints
kill all coefficients with t-degree<10 in that weighted window across
the three branches. This is a finite-jet **relaxation** of full prefix
vanishing, and contains every actual admissible primitive. Its ranks are
105 before and 107 after the three normal outputs are added; adjoining
the L-row alone leaves rank105. Adding both the primitive values and normal
outputs gives rank108; both branch-value differences annihilate the kernel.
The polynomial witnesses above attain
the two response dimensions. These modular controls challenge the analytic
proof; they are not characteristic-zero upper-bound arguments.

Normal and optimized executions pass **226** live exact gates and have
identical LF output. Reproduce with

```text
python -B 04-computation/overnight7_20260906_jc_fifth_module.py
python -B -O 04-computation/overnight7_20260906_jc_fifth_module.py
source SHA-256 0a478c1e94ff0b9ad69de2727801f6d0413f7df4a1275aa0a75f0885a4f5cdc1
output SHA-256 c73b297befb08788e6d5ad12472df9477d93f7cac8161bd25c56d6964ea63b5f
```

## 7. Portfolio comparison and connection contract

| Lane | Exact comparison | Information that must stay |
|---|---|---|
| Anchor: LRC body/event completion | A local admissible object can fail the consumer's actual entry condition | Source-plane admission here is separate from its target Jacobian predicate |
| Niche: Laurent carried root response | Complete carriers and both clocks fix the sign/consumer map | The s^5 correction is t^10 in the source and first visible at density t^9 |
| Wildcard: Smith unit jets / grid palettes | Restoring a forgotten coordinate repairs a particular quotient, not every later operation | Source Hamiltonian existence does not supply a descended primitive or erase a lower retained obstruction |

The concrete maps are polynomial source shears (6), target-Hamiltonian
normal restriction (7), and filtered two-form pullback (14). They preserve,
respectively, a closed source plane and finite collision repair; the actual
descended normal response; and every retained lower target-form coefficient.
Passing directly from the first map to a Keller conclusion loses the latter
two predicates. The cheapest decisive hostile is the order-ten payer itself:
it fixes the collision coefficient while the old order-eight obstruction
remains untouched.

The maintained question can therefore be replaced by a sharper one:
which **earlier** source jets, odd ninth-order terms, or other actual source
coordinates can change the J8 relation while retaining the intended lower
collision data? This report does not answer that next question. It also
does not classify unrestricted JC(2) source charts, all target primitives,
or polynomial termination. The fixed-prefix fifth-order lane is resolved
at the exact categories above, and the global problems remain open.
