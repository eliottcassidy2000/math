---
id: THM-4063
title: "Finite-graph period connection and higher-ramification no-go"
status: >
  PROVED ABSTRACT PERIOD/CONNECTION THEOREM + REFUTED COMPLETENESS ANTECEDENT
  FOR THE DISPLAYED FIGURE-EIGHT REALIZATION + CONDITIONAL ABSTRACT
  TWO-CYCLE COKERNEL + UNCONDITIONAL RAMIFICATION NO-GO + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. Under explicit period-completeness and full
  opening-lattice hypotheses, a finite branch graph has mixed cokernel
  R^beta/nabla L; aligned carrier valuations nu_i give length sum(nu_i-1).
  The displayed two-cycle carrier has relative exponents (q,2q), but THM-4067
  proves that its common ambient target is not period-complete: the fixed
  defect surjects R[[c]]/(A^q), and for q>=2 the actual mixed cokernel
  surjects R[[c]]/(A^(q-1)). Thus the conditional length 3q-2 is inapplicable
  to that realization. Separately, every pullback through H(t)=a t^e+...
  has Jacobian factor H'(t); e>=2 forbids a nonzero constant Jacobian,
  sharply. This proves no convergence, algebraization, global pair, or JC(2).
source: codex-frontier-synthesis-breakthrough-20260825 / mixed-form conductor lane
audit: >
  PASS after replacing a crossing polygon, typing endpoint evaluation,
  demoting the figure-eight cokernel to its exact conditional scope, and
  naming the source-ring Jacobian ideal. The SymPy companion checks the
  embedded moments, all pure-c periods through degree twelve, the triangle
  specialization, ramification valuations, and a characteristic-five
  hostile. A Fraction-only audit independently checks embedding, moment
  reconstruction, carrier ledgers through eight openings, and ramification
  leads through e=19. Both normal/optimized pairs byte-match; both scripts
  have zero assert nodes and zero float literals. THM-4067 subsequently
  hostile-audited the conditional antecedent itself and refuted it for the
  displayed T=R[[c,u]], M=all branchwise densities; this changes no abstract
  implication, carrier moment, or ramification result here.
depends_on:
  - THM-4058-exceptional-affine-triangle-period-and-simple-zero-monomial-ladder
  - THM-4060-exceptional-simple-zero-mixed-form-cokernel-collapse-and-formal-pair-lift
related:
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-4054-exceptional-affine-simple-zero-retained-packet
  - THM-4067-seminormal-period-kernel-and-figure-eight-completeness-obstruction
script: 04-computation/jc2_finite_graph_period_connection_ramification_thm4063.py
output: 05-knowledge/results/jc2_finite_graph_period_connection_ramification_thm4063.out
independent_audit_script: 04-computation/jc2_finite_graph_period_connection_ramification_thm4063_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_finite_graph_period_connection_ramification_thm4063_independent_audit.out
script_sha256: 3be9075dfca04a0389645ef980cb0c47b958ffc6fb8a74654c44a2fdf574a24f
output_sha256: b5326dd6bcaf4536e132a008530c20191b7a06534f0ccf8d9739302f3dae60e2
independent_audit_script_sha256: be239e482359bcb2920c7417ef1139ff568d7dcca8f43cbcc254712d8b203a8a
independent_audit_output_sha256: 08ae928e55c7b69283ba687f8250062200af5c2a87f263798ef6c3ac1b04d318
hash_basis: raw LF bytes
---

# THM-4063 -- mixed graph periods carry a connection, not a Betti scalar

**PROVED in the typed scopes below.** THM-4058/4060 show that a normalized
triangle period measures the fixed-output cokernel and that its first-output
variation is an Euler connection. The first part of this theorem identifies
the exact finite-graph mechanism under explicit completeness hypotheses. The
second part gives an independent sharp obstruction at higher ramification.

## 1. The formal graph-period setup

Let `k` be a characteristic-zero field, `R=k[[A]]`, and `F=k((A))`. Let
`Gamma` be a finite connected oriented graph. Choose a `k`-basis

```text
z_1,...,z_beta of H_1(Gamma;k),
beta=dim_k H_1(Gamma;k),                              (1)
```

and represent each `z_i` by its oriented edge coefficients.

Work in one centered branch coordinate `c`. Every vertex `v` has a section

```text
c_v(A) in A R.                                       (2)
```

For every oriented edge `e`, put `B_e=R[[c]]`. Condition `(2)` makes endpoint
evaluation and coefficientwise formal integration in `c` well defined. Let
`M` be an `R`-submodule of `direct_sum_e B_e`, and let `T` be an `R`-module
of common target germs with restriction map

```text
j:T -> direct_sum_e B_e.                             (3)
```

Assume:

1. `j(T) subset M`;
2. the restrictions of every `f in T` have one common value at each incident
   vertex;
3. the fixed-`c` derivatives

   ```text
   delta_c f=(partial_c j_e(f))_e,
   delta_A f=(partial_A j_e(f))_e                    (4)
   ```

   belong to `M`.

More general edge algebras are allowed if the same endpoint evaluation,
primitive, and differentiation maps are explicitly supplied; they are not
automatic from graph incidence.

For a cycle `z`, define

```text
P_z(m)=sum_e z_e integral_(c_s(e)(A))^(c_t(e)(A))
                         m_e(A,c) dc.                 (5)
```

## 2. Incidence cancellation and moving-endpoint transgression

For every `f,g in T`,

```text
P_z(delta_c g)=0,
P_z(delta_A f)=d P_z(j(f))/dA.                       (6)
```

The first identity is the endpoint sum pairing the common vertex values with
`partial z=0`. Differentiate `(5)` for the second. At a vertex `v`, the
moving-endpoint terms contribute

```text
j(f)(v)c_v'(A)
 [sum_(t(e)=v)z_e-sum_(s(e)=v)z_e]=0.                (7)
```

Thus all endpoint motion cancels by incidence, leaving the fixed-`c`
derivative in `(6)`. No completeness assumption is used in this section.

## 3. The period connection theorem

Put

```text
P=(P_z1,...,P_zbeta):M->F^beta.                      (8)
```

Now impose two load-bearing hypotheses:

1. **period completeness**

   ```text
   ker P=delta_c(T);                                  (9)
   ```

2. **full opening lattice:** `J=P(M)` is a free full-rank `R`-lattice in
   `F^beta`, so for some `B in GL_beta(F)`,

   ```text
   J=B R^beta.                                       (10)
   ```

Define

```text
Lambda=B^(-1)P,
L=Lambda(j(T)) subset R^beta,
nabla=d/dA+B^(-1)B',                                 (11)
```

where `L` is the normalized common-germ carrier submodule. For
`kappa in k^*`, let

```text
D:T direct_sum T -> M,
D(f,g)=delta_c g-kappa delta_A f.                    (12)
```

Then `nabla L subset R^beta`, and there is a canonical `k`-vector-space
isomorphism, after choosing `B`,

```text
boxed: coker D ~= R^beta/(nabla L).                  (13)
```

### Proof

By `(9)`, the period map identifies

```text
M/delta_c(T) ~= J.                                   (14)
```

Normalization by `B` identifies the right side with `R^beta`. If
`lambda=Lambda(j(f))`, then `(6)` gives

```text
Lambda(delta_A f)
 =B^(-1)(B lambda)'
 =lambda'+B^(-1)B'lambda
 =nabla lambda.                                      (15)
```

The remaining image of `(12)` in `(14)` is precisely
`kappa*nabla L`; the nonzero scalar is harmless. This proves `(13)` and also
shows `nabla L subset R^beta`.

The quotient need not be an `R`-module because `nabla` is a connection, not
an `R`-linear map. A filtration statement requires either strictness
hypotheses or the quotient filtration transported through `(14)`; none is
silently asserted here.

If `B` is replaced by `BU`, `U in GL_beta(R)`, then `L` changes to
`U^(-1)L` and `(11)` undergoes the usual gauge transformation. Multiplication
by `U^(-1)` identifies the two quotients, so `(13)` is independent of the
chosen opening frame.

## 4. Aligned carrier valuations

Suppose that after one constant basis change,

```text
B=diag(A^q_i),
L=direct_sum_i A^nu_i R e_i,
q_i>=0,                  nu_i>=1.                    (16)
```

Then

```text
nabla(A^(nu_i+n)e_i)
 =(nu_i+n+q_i)A^(nu_i+n-1)e_i.                       (17)
```

Every displayed coefficient is nonzero in characteristic zero. Hence

```text
nabla L=direct_sum_i A^(nu_i-1)R e_i,
coker D ~= direct_sum_i R/A^(nu_i-1),
dim_k coker D=sum_i(nu_i-1).                          (18)
```

Thus the relative carrier valuations `nu_i`, not `beta` or the opening
orders alone, determine the aligned length. In the unaligned case the exact
invariant is `(L,nabla)` itself. If `nabla L` is an `R`-lattice with Smith
exponents, their sum gives the length; if `L` is rank deficient or the image
is not a lattice, the cokernel may be infinite-dimensional.

For THM-4058/4060, after a constant normalization,

```text
Gamma=triangle, beta=1,
B=A^5, L=A^5R,
nabla=d/dA+5/A.                                      (19)
```

Equation `(18)` is exactly `R/A^4`, of dimension four. This is a theorem
specialization, not a new proof of THM-4060's period completeness or
closed-form presentation.

## 5. An embedded two-cycle carrier hostile

Let `epsilon=A^q`, `q>=1`. Take two oriented triangles in the `(c,u)` plane,

```text
C_1: (0,0),(1,0),(2,1),
C_2: (0,0),(-1,2),(-3,1),                            (20)
```

and scale every coordinate by `epsilon`. The closed triangles lie in
opposite half-planes `c>=0` and `c<=0`, share only the origin, and every edge
has nonzero `c`-increment. Their union is an embedded figure-eight graph with
`beta=2`.

Let `T=R[[c,u]]` be the common ambient-germ module, let `M` be all branchwise
densities, and use the two raw cycle periods `integral_(C_i) f dc`. Exact
integration gives

```text
v_u =(-1/2,-5/2),
v_cu=(-1/2,10/3),
det[v_u,v_cu]=-35/12.                                (21)
```

The exact common-germ period carrier is

```text
C=epsilon^2 R v_u+epsilon^3 R v_cu.                 (22)
```

To prove all degrees in `(22)`, expand

```text
f=sum a_(n,r,s) A^n c^r u^s.                        (23)
```

Every pure `c` period vanishes. The unique degree-one survivor is `u`, with
vector `v_u`. A monomial of total `(c,u)` degree `d>=2` contributes
`A^n epsilon^(d+1)m_(r,s)`. Since the two vectors in `(21)` form a basis,
write `m_(r,s)=alpha v_u+beta v_cu`; then its contribution belongs to the
right side of `(22)`. Reverse containment comes from `f=u` and `f=cu`, with
arbitrary coefficients in `R`.

Every branchwise period is divisible by `epsilon`. A constant density on one
nonshared edge of each loop gives `epsilon` times the corresponding cycle
coordinate, so

```text
J=epsilon R^2.                                       (24)
```

Let `V=[v_u,v_cu]`. Since `V in GL_2(k)`, choose the opening frame
`B=epsilon V`. Equations `(22)--(24)` give

```text
L=epsilon R e_1 direct_sum epsilon^2 R e_2,
(nu_1,nu_2)=(q,2q),
nabla=d/dA+(q/A)I.                                   (25)
```

Consequently, **for any separately proved period-complete realization with
this same carrier data**, `(18)` gives

```text
coker D ~= R/A^(q-1) direct_sum R/A^(2q-1),
dim_k coker D=3q-2.                                  (26)
```

This differs from the Betti-only guess `2(q-1)` by `q`. Equation `(26)` is an
abstract conditional implication, and no such realization is constructed
here.

THM-4067 subsequently resolves the antecedent **negatively for the displayed
choice `T=R[[c,u]]`, `M=direct_sum_e R[[c]]` above**. The nonincident edges

```text
e_3:u=c/2,                  e_5:u=c/2+(5/2)epsilon
```

force every common-target derivative to satisfy
`(delta_c f)_3=(delta_c f)_5 mod epsilon`. Nevertheless, with the orientations
in `(20)`,

```text
m=(0,2,1,0,0,0)
```

has periods `0*epsilon+2*epsilon-2*epsilon=0` and `0`, while
`m_3-m_5=1`. More strongly, THM-4067 proves

```text
ker P/delta_c(T) ->> R[[c]]/(A^q),                    (26a)

coker D ->> R[[c]]/(A^(q-1))       for q>=2.          (26b)
```

Both targets are infinite-dimensional over `k` when nonzero. Thus `(26)` is
inapplicable to this displayed common-ambient realization. Its unconditional
survivor is the exact moment carrier `(22)`, opening lattice `(24)`, and the
warning that graph incidence and opening order do not determine either the
common-germ carrier or the hidden graph-gluing defect.

## 6. Completeness and characteristic firewalls

- THM-4067 identifies the missing invariant exactly. With edgewise primitives,
  `ker P/delta_c(T)` is the quotient of the graph value-equalizer by the
  common-target restriction algebra. In reduced curve settings this is the
  seminormal defect; the mixed cokernel has an additional exact left term
  before the connection quotient `(13)`.
- Incidence does not imply period completeness. The curves
  `y(y-x^m)=0` have the same two-branch graph but normalization image
  `{(f,g):f=g mod x^m}` and conductor
  `x^m k[[x]] direct_sum x^m k[[x]]`; contact thickness is invisible to the
  graph. THM-3696 is the canonical three-branch version, with an extra
  derivative law and conductor `b^2(1-b^2)^2 k[b]`.
- A tree has `beta=0`, but period completeness/no-hidden-conductor-jets is
  still needed for surjectivity.
- Zero area, dependent moment vectors, nonuniform opening orders, or an
  `A`-dependent Smith basis can change `L` and the connection.
- In characteristic `p`, coefficients in `(17)` vanish at infinitely many
  degrees. The aligned image need not be a lattice and the cokernel can
  become infinite-dimensional. Characteristic zero is load-bearing.
- Equation `(13)` computes the displayed mixed operator. Identifying it with
  every closed target two-form also requires THM-4060's formal closed-form
  presentation and compatibility; a period solution is not automatically a
  pair.

## 7. Higher ramification has a sharp Jacobian divisor

Let `k` have characteristic zero and

```text
H(t)=a t^e+O(t^(e+1)),       a!=0, e>=1.             (27)
```

Suppose target pullbacks factor through

```text
Phi_H:(x,t)->(x,u)=(x,H(t)).                          (28)
```

For formal target germs `Fbar,Gbar in k[[x,u]]`, the chain rule gives

```text
Jac_(x,t)(Fbar(x,H(t)),Gbar(x,H(t)))
 =H'(t)Jac_(x,u)(Fbar,Gbar)(x,H(t)).                 (29)
```

Since `ord_t H'=e-1`, every such Jacobian lies in

```text
(H'(t))=(t^(e-1)) inside k[[x,t]].                   (30)
```

The same calculation for any target two-form gives

```text
Phi_H^*(Omega)=H'(t)K_Omega(x,H(t)) dx wedge dt.     (31)
```

Therefore

```text
boxed: e>=2 implies no target pair has a nonzero
       constant source Jacobian.                    (32)
```

This no-go is independent of graph periods and already holds for the whole
two-form image. It cannot be repaired by Darboux factorization.

The divisor is sharp whenever the pre-ramification target algebra contains a
unit-Jacobian pair. In the exceptional affine packet, `(A,-4c)` has
pre-substitution Jacobian `12`, so its pullback has Jacobian `12H'(t)`.
Consequently the ideal generated in `k[[x,t]]` by all target-pair Jacobians is
exactly `(H')`.

For `e=1`, `H'(0)` is a unit and `H` has a formal inverse. This removes the
ramification obstruction but does not prove existence: the graph-side
criterion and closed-form/pair gates remain. A sufficient carrier condition
for common-germ corrections is `L subset nabla L`; THM-4060 has
`A^5R subset A^4R`.

In three formal target variables, a closed two-form with a unit coefficient
has rank two; its kernel can be formally straightened and the form factors as
`dF wedge dG`. Without a unit coefficient that Darboux step is not licensed.
The need to divide by nonunit `H'` at `e>=2` is useful intuition, but no
Puiseux equivalence is asserted.

## 8. Audit and scope

The primary companion verifies the embedded moment matrix, pure-`c`
telescoping, THM-4060 specialization, ramification valuation, and the first
positive-characteristic resonance. The independent companion uses only exact
fractions and elementary polynomial arithmetic, checks the embedding and
moment basis separately, and replays the carrier/ramification ledgers.

Reproduce both normal/optimized pairs with

```text
python3 -B 04-computation/jc2_finite_graph_period_connection_ramification_thm4063.py
python3 -B -O 04-computation/jc2_finite_graph_period_connection_ramification_thm4063.py
python3 -B 04-computation/jc2_finite_graph_period_connection_ramification_thm4063_independent_audit.py
python3 -B -O 04-computation/jc2_finite_graph_period_connection_ramification_thm4063_independent_audit.py
```

The displayed figure eight is now proved period-incomplete by THM-4067; for
`q>=2` its actual mixed cokernel already has an infinite contact quotient.
This theorem and its repair give no convergence or algebraization of a formal
pair, no polynomial Keller map, and no progress from the local packet to
global `JC(2)` or `DC(2)`. **QED.**
