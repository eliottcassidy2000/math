---
id: THM-4067
title: "Seminormal period kernel and figure-eight completeness obstruction"
status: >
  PROVED ABSTRACT EXACT SEQUENCES + PROVED ORDINARY-NODE SUFFICIENT
  CRITERION + PROVED CONDUCTOR-LENGTH COROLLARIES + REFUTED DISPLAYED
  THM-4063 FIGURE-EIGHT COMPLETENESS ANTECEDENT + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. With exact edgewise primitives, zero graph
  periods reconstruct precisely the graph value-equalizer S_Gamma, so
  ker(P)/delta_c(T) is S_Gamma/Acal. The mixed cokernel has an exact left
  correction by this defect before the THM-4063 connection quotient. Genuine
  ordinary-node realizations are period-complete; y(y-x^m)=0 has defect
  length m-1. In THM-4063's displayed figure eight, two nonincident parallel
  edges force a hidden congruence modulo A^q. Its fixed period defect
  surjects k[[A,c]]/(A^q) for every q>=1, and for q>=2 its mixed cokernel
  surjects k[[A,c]]/(A^(q-1)); both are infinite-dimensional over k. The
  abstract THM-4063 connection, carrier moments, and ramification no-go
  survive. No global pair or consequence for JC(2) follows.
source: codex-frontier-synthesis-creative-20260825b / graph-period lane
audit: >
  PASS. The SymPy companion checks abstract graph dimension ledgers, contact
  and ordinary-multiple-point lengths, the THM-3696 missing class, all six
  figure-eight lines and orientations, the minimal witness, 504 fixed and
  504 mixed monomial contact valuations, and 126 explicit quotient lifts. A
  no-import Fraction-only audit independently reconstructs incidence/cycle
  ranks, exact contact derivative matrices, branch equations, orientations,
  the witness, and the same 1,008 valuations using a binomial dictionary
  engine. Normal/optimized outputs byte-match; both scripts have zero assert
  nodes and zero float literals.
depends_on:
  - THM-4063-finite-graph-period-connection-and-ramification-no-go
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
related:
  - THM-3955-node-cotangent-normalization-kernel-and-conductor-torsion
  - THM-3957-triple-normal-crossing-cotangent-conductor-kernel-and-normalization-cokernel
  - THM-4058-exceptional-affine-triangle-period-and-simple-zero-monomial-ladder
  - THM-4060-exceptional-simple-zero-mixed-form-cokernel-collapse-and-formal-pair-lift
script: 04-computation/graph_seminormal_period_kernel_thm4067.py
independent_script: 04-computation/graph_seminormal_period_kernel_thm4067_independent_audit.py
output: 05-knowledge/results/graph_seminormal_period_kernel_thm4067.out
independent_output: 05-knowledge/results/graph_seminormal_period_kernel_thm4067_independent_audit.out
script_sha256: 09f790d639b29f9174dfdc6c0c538ac1f8a57ca768cb6fe0e79137d81d13fab6
independent_script_sha256: f76f5149a799c06365c52377582d2914c9158acee0ff3add625ece5f9199f9ea
output_sha256: 5f65180cc121d2840dbdd2c516842c27428228bc7df82800d5c8a87c7095745f
independent_output_sha256: 7a81491f1ceecf6367827cf7a2c39985495ce54d7add97f69d8d989a62b7293d
hash_basis: raw LF bytes
---

# THM-4067 -- graph periods end at the seminormal boundary

**PROVED in the typed scopes below.** THM-4063 correctly made period
completeness a load-bearing hypothesis. This theorem identifies that
hypothesis exactly, adds the missing term when it fails, and resolves the
displayed figure-eight antecedent negatively. The failure is not Betti data:
it comes from a contact congruence between nonincident supporting branches.

## 1. Exact edge-Poincare data and the graph equalizer

Let `Lambda` be a commutative `Q`-algebra and let `Gamma` be a finite
connected oriented graph. For every edge `e`, suppose that there are
`Lambda`-modules and maps

```text
Lambda -> B_e --d_e--> M_e,
ev_(e,s),ev_(e,t):B_e->Lambda,
I_e:M_e->Lambda                                      (1)
```

such that

```text
0 -> Lambda -> B_e --d_e--> M_e -> 0,
I_e(d_e b)=ev_(e,t)(b)-ev_(e,s)(b).                  (2)
```

The first map in `(2)` is the constant inclusion. Thus every edge density
has a primitive, unique up to a constant, and `I_e` obeys the fundamental
theorem. For THM-4063,

```text
Lambda=R=k[[A]],       B_e=M_e=R[[c]],
d_e=partial_c,                                         (3)
```

and `I_e` is coefficientwise formal integration between endpoint sections
in `A R`. Characteristic zero is load-bearing for primitive division.

Put

```text
Btilde=direct_sum_e B_e,       M=direct_sum_e M_e,

S_Gamma={b in Btilde:
 all incident endpoint values ev_(e,v)(b_e)
 agree at every vertex v}.                            (4)
```

Let `T` be a module of common target functions with restriction map
`j:T->Btilde`, assume `j(T) subset S_Gamma`, and write

```text
Acal=j(T).                                             (5)
```

Assume diagonal constants lie in `Acal`. For `m in M`, let `a(m)` be the
edge one-cochain with entries `I_e(m_e)`, and define

```text
P(m)=[a(m)] in H^1(Gamma;Lambda),       J=P(M).        (6)
```

Equivalently, choose a cycle basis and record every oriented cycle integral.

## 2. The period-kernel exact sequence

There is a canonical exact sequence of additive `Lambda`-modules

```text
boxed:
0 -> S_Gamma/Acal --dbar--> M/d(Acal) --Pbar--> J -> 0. (7)
```

Consequently

```text
ker P/d(Acal) ~= S_Gamma/Acal,                         (8)

ker P=d(Acal) iff Acal=S_Gamma.                       (9)
```

Thus, under exact edgewise primitives, period completeness is exactly the
claim that common target functions realize every graph-continuous primitive.

### Proof

If `b in S_Gamma`, equation `(2)` makes every cycle integral of `db`
the oriented sum of endpoint values. These telescope at each vertex, so

```text
d(S_Gamma) subset ker P.                              (10)
```

Conversely, take `m in ker P` and choose an edgewise primitive `b_e` by
`(2)`. The increments `I_e(m_e)` form a graph one-cochain with zero sum on
every cycle. Choose a root vertex and integrate the increments along a
spanning tree. The cycle conditions make the resulting vertex potentials
`phi_v` independent of path and give

```text
I_e(m_e)=phi_(t(e))-phi_(s(e)).                       (11)
```

Shift each `b_e` by its unique edgewise constant so that its source value is
`phi_s`. Equation `(2)` then gives target value `phi_t`; the shifted tuple is
in `S_Gamma` and differentiates to `m`. Hence

```text
ker P=d(S_Gamma).                                     (12)
```

The kernel of `d` on `S_Gamma` consists of edgewise constants. Connectivity
and the vertex equalities make them one diagonal constant, which belongs to
`Acal`. Therefore differentiation induces an injection
`S_Gamma/Acal->M/d(Acal)`. Its image is the kernel of `Pbar` by `(12)`, and
`Pbar` is onto `J` by definition. This proves `(7)--(9)`.

## 3. The conductor correction to the mixed connection

Now specialize to `Lambda=R=k[[A]]`, with `k` of characteristic zero. Assume
the total fixed-`c` base derivative of restricted targets descends to a map

```text
delta_A:Acal->M.                                      (13)
```

For `kappa in k^*`, put

```text
D:Acal direct_sum Acal -> M,
D(f,g)=d g-kappa delta_A f.                           (14)
```

Let

```text
Q=S_Gamma/Acal ~= ker P/d(Acal),
E=M/d(Acal),
U=image(delta_A Acal -> E).                           (15)
```

There is an exact sequence of `k`-vector spaces

```text
boxed:
0 -> Q/(Q intersect U) -> coker D
  -> J/P(delta_A Acal) -> 0.                          (16)
```

The intersection is inside `E`, and no splitting is asserted. Indeed `(7)`
is `0->Q->E->J->0`; quotienting by `U` gives `(16)` because the new kernel
is `(Q+U)/U ~= Q/(Q intersect U)`.

Under THM-4063's full opening-lattice hypothesis, choose its frame `B` and
carrier `L`. Its transgression computation identifies the rightmost term as

```text
J/P(delta_A Acal) ~= R^beta/(nabla L).                (17)
```

Thus THM-4063's connection quotient is the exact right quotient in `(16)`.
Period completeness is the sufficient condition `Q=0` that deletes the left
correction. When `Q!=0`, first-output motion may kill part of it, but that
payment is the explicit intersection `Q intersect U`, not a Betti count.

## 4. Seminormal and conductor interpretation

Suppose `Acal` is a reduced curve-function ring, `Btilde` is its finite branch
normalization, and the graph equalizer `(4)` is its seminormalization
`Acal^sn`. Then

```text
Q=Acal^sn/Acal.                                       (18)
```

This is the seminormal defect, not the whole normalization defect. If

```text
conductor c=(Acal:Btilde),                             (19)
```

then `cQ=0`, since `c Acal^sn subset c Btilde subset Acal`, and

```text
0 -> Acal^sn/Acal -> Btilde/Acal
  -> Btilde/Acal^sn -> 0.                             (20)
```

For a completed reduced curve singularity over an algebraically closed field,
with `r` smooth branches and delta invariant `delta`, the last quotient in
`(20)` records only the `r-1` independent residue-value mismatches. Hence

```text
length(Q)=delta-(r-1).                                (21)
```

### 4.1 Ordinary nodes are sufficient

At an ordinary node,

```text
k[[x,y]]/(xy)
 ~= {(f,g) in k[[x]] direct_sum k[[y]]:f(0)=g(0)}.    (22)
```

There is no defect in `(18)`. More intrinsically, for the two local branch
ideals `I,J`, the exact Chinese-remainder sequence

```text
0 -> T/(I intersect J) -> T/I direct_sum T/J
  -> T/(I+J) -> 0                                    (23)
```

says that equality in the node residue ring is the complete gluing condition.
Since the normalization quotient is supported at the nodes, the local
equalities globalize. Therefore a connected realization with exact edge
primitives, only ordinary two-branch nodes, and no unrecorded ambient
intersections or contacts satisfies `Acal=S_Gamma` and is period-complete.

The last two exclusions are load-bearing: an abstractly embedded compact
edge graph need not record intersections or contacts of its full supporting
formal branches.

### 4.2 Exact contact and multiple-point lengths

For the canonical hostile

```text
A_m=k[[x,y]]/(y(y-x^m))
   ={(f,g) in k[[x]]^2:f-g in x^m k[[x]]},            (24)

A_m^sn={(f,g):f(0)=g(0)},
conductor_(Btilde/A_m)=x^m k[[x]] direct_sum x^m k[[x]],
```

one obtains

```text
A_m^sn/A_m ~= x k[[x]]/x^m k[[x]],
length=m-1,                                           (25)

dA_m={(p,q):p-q in x^(m-1)k[[x]]}.                   (26)
```

The converse in `(26)` follows by integrating `p-q` and choosing equal
constants. Thus even a tree, with no cycle periods, can have exact
fixed-output cokernel length `m-1`.

For a plane ordinary `r`-fold point with distinct smooth tangents,
`delta=binom(r,2)`, so `(21)` gives

```text
length(Q)=binom(r-1,2).                               (27)
```

Nodes have zero defect, plane triple points have one, and plane fourfolds
have three. “Ordinary” does not replace the two-branch nodal hypothesis.

### 4.3 THM-3696 is the one-class triple specialization

For THM-3696, put

```text
S={f in C[b]:f(-1)=f(0)=f(1)}.                        (28)
```

Its exact reconstruction theorem says

```text
R_0={f in S:f'(-1)+4f'(0)+f'(1)=0}.                  (29)
```

Therefore `S/R_0 ~= C`, as predicted by `(27)` for `r=3`. The explicit
element

```text
b^2X=b^3(1-b^2)                                      (30)
```

has value zero at all three addresses and derivative response `-4`, so it
generates the missing class. THM-3696's conductor `X^2C[b]` annihilates it.

This function defect is not the cotangent defect of THM-3955/3957. A node
has `Q=0` here even though its Kahler differentials contain conductor torsion;
the coordinate triple crossing is seminormal on functions even though
normalized differentials have a conductor-axis cokernel. Those distinct
types must not be merged.

## 5. The displayed THM-4063 figure eight is not period-complete

Return to

```text
k characteristic zero,       R=k[[A]],
epsilon=A^q,                  q in Z, q>=1,
T=R[[c,u]],                   M=direct_sum_(e=1)^6 R[[c]]. (31)
```

Use THM-4063's two oriented triangles, now writing their scaled vertices as

```text
C_1: O=(0,0), B_1=(epsilon,0), B_2=(2epsilon,epsilon),
C_2: O=(0,0), D_1=(-epsilon,2epsilon),
                     D_2=(-3epsilon,epsilon).         (32)
```

Orient them in the displayed order and close each cycle. The six edge
orientations and `c` bounds are

```text
e_1: O->B_1,       0->epsilon,
e_2: B_1->B_2,     epsilon->2epsilon,
e_3: B_2->O,       2epsilon->0,
e_4: O->D_1,       0->-epsilon,
e_5: D_1->D_2,     -epsilon->-3epsilon,
e_6: D_2->O,       -3epsilon->0.                     (33)
```

Their supporting line equations are

```text
e_1: u=0,                       e_4: u=-2c,
e_2: u=c-epsilon,               e_5: u=c/2+(5/2)epsilon,
e_3: u=c/2,                     e_6: u=-c/3.          (34)
```

The nonincident edges `e_3,e_5` are parallel and separated by
`(5/2)epsilon`. For a common target `f`, write

```text
r_3=f(A,c,c/2),
r_5=f(A,c,c/2+(5/2)epsilon).                           (35)
```

Taylor expansion in `u` gives

```text
r_3-r_5 in epsilon R[[c]],
partial_c r_3-partial_c r_5 in epsilon R[[c]].        (36)
```

Both derivatives in `(36)` use the same total branch operator
`partial_c+(1/2)partial_u`; this is exactly THM-4063's
`delta_c f=partial_c(j_e(f))`, not ambient partial differentiation.

It follows that

```text
chi_fix:ker P/delta_c(T) -> R[[c]]/(epsilon),
[m] |-> [m_3-m_5]                                      (37)
```

is well-defined. It is surjective. Given `phi(c) in R[[c]]`, set

```text
H(c)=integral_0^c phi(s) ds,       v=H(2epsilon),

h_1=0,
h_2=(v/epsilon)(c-epsilon),
h_3=H(c),
h_4=h_5=h_6=0.                                       (38)
```

Since `H(0)=0`, one has `v in epsilon R`, so `(38)` is regular. Its endpoint
values agree at every graph vertex: on the first triangle they are
`0,0,v,0`, and all values on the second are zero. Hence `h in S_Gamma`, so
`m=dh` has zero graph periods, while `(37)` sends it to `[phi]`. Therefore

```text
boxed:
ker P/delta_c(T) ->> R[[c]]/(A^q).                    (39)
```

The quotient is infinite-dimensional over `k` for every `q>=1`. In
particular, THM-4063's period-completeness hypothesis is false for its
displayed common-ambient realization.

The minimal witness takes `phi=1`:

```text
h=(0,2(c-epsilon),c,0,0,0),
m=delta_c h=(0,2,1,0,0,0).                           (40)
```

Using the oriented bounds `(33)`, its periods are

```text
P_(C_1)(m)=0*epsilon+2*epsilon+1*(-2epsilon)=0,
P_(C_2)(m)=0,                                         (41)
```

but `m_3-m_5=1` violates `(36)`. This is an all-`q` proof, not a finite
truncation.

## 6. The actual mixed cokernel is already infinite for q>=2

The same parallel pair gives a stronger consequence for THM-4063's displayed
mixed operator. Differentiate `(35)` with respect to `A` at fixed `c`.
Because the branch separation has order `q`, one loses at most one order:

```text
delta_A(T)_3-delta_A(T)_5 in A^(q-1)R[[c]].           (42)
```

Equation `(36)` gives the stronger `A^q` divisibility for `delta_c(T)`.
Therefore, for

```text
D(f,g)=delta_c g-kappa delta_A f,
```

the branch-difference map descends to a surjection

```text
boxed:
coker D ->> R[[c]]/(A^(q-1)).                         (43)
```

Surjectivity follows by choosing an arbitrary branchwise density on `e_3`
and zero on `e_5`. The order in `(42)` is sharp: for the target `f=u`,

```text
delta_A(f|e_3)-delta_A(f|e_5)
 =-(5q/2)A^(q-1).                                     (44)
```

Consequently the displayed mixed cokernel is infinite-dimensional over `k`
for every `q>=2`. At `q=1`, `(43)` has zero target and makes no mixed-cokernel
claim; the fixed period defect `(39)` remains nonzero and infinite.

## 7. Exact repair of THM-4063

No abstract claim in THM-4063 is retracted. Its theorem explicitly assumed
period completeness, and its figure-eight length was explicitly conditional.
The exact disposition is:

1. incidence cancellation and moving-endpoint transgression survive;
2. the abstract connection theorem under period completeness and a full
   opening lattice survives and is strengthened by `(16)`;
3. the exact figure-eight moments, opening lattice, and carrier exponents
   `(q,2q)` survive;
4. the conditional length `3q-2` is inapplicable to the displayed
   `T=R[[c,u]]`, because `(39)` refutes its antecedent and `(43)` gives an
   infinite actual mixed quotient for `q>=2`;
5. `3q-2` remains valid for a separately proved period-complete realization
   with the same carrier data; none is constructed here. Simply enlarging
   the target to `S_Gamma` can change `L` and is not such a construction; and
6. THM-4063's higher-ramification Jacobian divisor and no-go survive
   unchanged.

The first failed implication was not in the conditional theorem statement.
It was the geometric search heuristic that an embedded union of compact edge
segments could be treated as independent supporting formal branches. The
parallel pair `(34)` is disjoint as segments but congruent modulo `A^q` as
formal target restrictions.

## 8. Audit, connection contract, and strict boundary

The primary companion uses exact SymPy arithmetic. It checks graph dimension
ledgers, conductor lengths, `(30)`, all lines and orientations, `(40)--(41)`,
504 fixed and 504 mixed monomial valuations through total degree six for
every `q=1,...,6`, 126 explicit lifts from `(38)`, and sharpness `(44)`.

The independent companion imports no production code and uses only
`Fraction`, binomial expansion, and a standalone row reducer. It reconstructs
the incidence and cycle ranks, complete contact derivative matrices for
`m=1,...,8`, all branch equations and signs, and the same 1,008 valuations by
a coefficient-dictionary path. These finite windows are hostile controls;
the arbitrary-degree proofs are `(2)`, `(36)`, `(38)`, and `(42)`.

| field | contract |
|---|---|
| source | common target restrictions `Acal`, or the exact displayed ambient ring `R[[c,u]]` in Sections 5-6 |
| target | all branchwise densities `M`, with oriented graph periods |
| map | edge derivative `d`; in the mixed case `D=delta_c-kappa delta_A` |
| preserved predicate | exact period-kernel membership, graph-continuous primitive descent, conductor/contact quotient, and the displayed fixed/mixed nonmembership |
| destroyed information | convergence, algebraicity, polynomial degree, behavior outside the branch packet, injectivity, properness, and global target geometry |
| needed sidecar | for a positive construction, proof that `S_Gamma/Acal` is killed in the mixed complex, followed by closed-form factorization, convergence/algebraization, and global control |
| cheapest decisive tests | compare every pair of formal branch restrictions modulo its contact ideal; test `(40)` and the sharp target `f=u` before computing a carrier Smith form |

Reproduce the two independent paths in normal and optimized modes:

```bash
python3 -B 04-computation/graph_seminormal_period_kernel_thm4067.py
python3 -B -O 04-computation/graph_seminormal_period_kernel_thm4067.py
python3 -B 04-computation/graph_seminormal_period_kernel_thm4067_independent_audit.py
python3 -B -O 04-computation/graph_seminormal_period_kernel_thm4067_independent_audit.py
```

Each stream byte-matches its stored output. Both scripts have zero Python
assert nodes and zero float literals; the four hashes are pinned in the
frontmatter.

This theorem proves no convergence, algebraization, global polynomial pair,
Keller map, or consequence for `JC(2)` or `DC(2)`. It supplies the exact
missing local sidecar and closes only the displayed figure-eight antecedent.
**QED.**
