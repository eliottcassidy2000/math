---
id: THM-4167
title: "Tournament exposure-capacity deletion support moment and parity holonomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Every exposure-capacity
  coordinate is a positive tagged OCF sum. Vertex deletion removes exactly
  the atoms using the deleted vertex, giving coordinatewise monotonicity and
  an exact outside-support moment. The quadratic Johnson numerator and
  denominator have support sizes three and four, so their restriction deck
  obeys an exact parity-dependent normalized transport law. Odd parent order
  contracts the weighted child tilt by one half; even parent order amplifies
  it by two.
source: codex-frontier-synthesis-creative-20260826at
depends_on:
  - THM-002-ocf
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4123-balanced-cardinality-ear-average-growth-and-johnson-layer-lattice
  - THM-4128-johnson-slice-support-envelope-and-exposure-centrality-criterion
  - THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets
related:
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4137-strong-tournament-centrality-complete-order-ten
  - THM-4163-order-eleven-homogeneous-pair-johnson-centrality
  - THM-4177-root-split-johnson-current-and-source-sink-boundary-census
script: 04-computation/tournament_capacity_deletion_holonomy_thm4167.py
output: 05-knowledge/results/tournament_capacity_deletion_holonomy_thm4167.out
independent_audit_script: 04-computation/tournament_capacity_deletion_holonomy_thm4167_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_capacity_deletion_holonomy_thm4167_independent_audit.out
script_sha256: a44adacdd1b7cdccccb9a81225f99073a8e468136b306c7369c1e8770d3e1f8a
output_sha256: 34cd665df7d7b3c00293d8f50a3280c240841c22e1746efdd5a0d81baec0a99c
independent_audit_script_sha256: 236e807a7bc4120a2fdd1fe803ac7adb35d9bd1bd05a49452041adc50038c1b1
independent_audit_output_sha256: 9cc785a952b98373d5cae3894ab00d789847d6cd27057a0d7bd75defff708ce1
hash_basis: raw LF bytes
primary_audit: >
  PASS. A self-contained exact implementation checks 410 literal tagged-OCF
  capacity/support-moment gates through order four, 31,512 marked-word
  deletion identities through order five, all 1,966,080 capacity-deletion
  comparisons at order six, four arbitrary-tensor restriction controls, and
  the named order-eleven and order-twelve hostiles. Normal, optimized, and
  hash-seeded streams byte-match.
independent_audit: >
  ACCEPT. A clean-room C++ referee uses literal exposed-word enumeration for
  the small universes and four-ear Hamilton-path differences for the named
  large rows. It independently repeats all 31,512 marked identities, all
  1,966,080 order-six monotonicity gates, the restriction algebra, and both
  hostile packets. Clang O0/O3 and GCC O3 streams byte-match.
---

# THM-4167 -- capacity deletion support and parity holonomy

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The order-eleven Johnson problem cannot be inherited from order ten by
deleting a vertex and then using the actual card capacity. The parent
capacity restricted to a card contains additional odd-cycle mass whose
support used the deleted vertex. This theorem identifies that missing tensor
exactly and then shows that the normalized Johnson tilt nevertheless has a
simple restriction-deck average. The normalization alternates with parity:
odd order gains a factor `1/2`, while even order pays a factor `2`.

## 1. Tagged OCF formula

Let `T` be a tournament on vertex set `V`, let `e={i,j}`, and adjoin a new
vertex `x`. For `S subseteq V`, orient its incident arcs by

```text
x -> u  iff  u in S,
F_T(S)=H(T+x_S).                                           (1)
```

Retain the integer exposure capacity of THM-4128/4162,

```text
c_e(T)=F_T({i})+F_T({j})-F_T(empty)-F_T({i,j}).            (2)
```

Thus `c_e=Q(i,j)+Q(j,i)=2w_e`. Define `P_e(T)` as a **tagged** family of
odd-cycle packings:

- in ambient state `T+x_{ {i} }`, the unique cycle using `x` has local
  segment `j -> x -> i`;
- in ambient state `T+x_{ {j} }`, the unique cycle using `x` has local
  segment `i -> x -> j`.

The ambient state is part of the tag, and `U(Gamma)` denotes the union of the
cycle vertices in the packing.

> **Theorem 1 (positive capacity atoms).** For every tournament and edge,
>
> ```text
> c_e(T)=sum_(Gamma in P_e(T)) 2^|Gamma|.                  (3)
> ```

### Proof

Apply THM-002 to the four terms of `(2)`. Every packing avoiding `x` occurs
in all four states and cancels. A packing using `x` has one distinguished
cycle. In a singleton state, if the predecessor of `x` is outside `e`, the
canonically corresponding cycle collection occurs in the two-element state
with opposite sign. Those terms cancel. The only survivors have the predecessor and successor of
`x` equal to the two endpoints of `e`; they are precisely the two tagged
families above. Each retains its ordinary OCF weight `2^|Gamma|`. **QED.**

The tags are load-bearing: a directed cycle is an oriented cyclic word, and
the two ambient ear states must not be identified merely because their
supports agree.

## 2. Deletion support moment

Fix `v notin e`. A tagged packing belongs to `P_e(T-v)` exactly when the
corresponding packing of `T` avoids `v`. Therefore

```text
delta_(v,e)
 :=c_e(T)-c_e(T-v)
 =sum_(Gamma in P_e(T), v in U(Gamma))2^|Gamma|
 >=0.                                                      (4)
```

This proves coordinatewise exposure-capacity monotonicity under deletion in
every order. Put

```text
A_e(T)=sum_(v notin e)c_e(T-v),
Omega_e(T)=sum_(Gamma in P_e(T))
             |U(Gamma) intersect (V minus e)| 2^|Gamma|.   (5)
```

Double-counting pairs `(v,Gamma)` gives the exact support moment

```text
boxed: (n-2)c_e(T)=A_e(T)+Omega_e(T),
       Omega_e(T)>=0.                                     (6)
```

Thus `Omega` is not a normalization error. It is exactly the first moment of
the outside support of the positive tagged OCF atoms. The labelled card
tensors lose precisely this mass.

## 3. Marked-word refinement and compensation

Let `E_e(T)` be the vertex words in which the endpoints of `e` are adjacent
and every other adjacent step follows an arc of `T`. The exposed-block
interpretation of THM-4128 gives

```text
|E_e(T)|=c_e(T).                                          (7)
```

For `P=(P_1,...,P_(n-1)) in E_e(T-v)`, put `s_k=1` when
`v -> P_k`, and `s_k=0` otherwise. Define

```text
r_v(P)=#{k:s_k=1, s_(k+1)=0},
m_(v,e)(P)=1  iff the marked e-gap has signature 0 -> 1.  (8)
```

The number of legal insertions of `v` that do not split the marked gap is

```text
a_(v,e)(P)=1+r_v(P)-m_(v,e)(P).                           (9)
```

Indeed the endpoint insertions plus all internal `0 -> 1` slots total
`1+#(1 -> 0)` by binary transition balance; the marked legal slot must then
be removed. Let `O_(v,e)` count full exposed words whose deletion of `v`
creates a reversed shortcut and hence does not land in `E_e(T-v)`. Partition
the full words into deletion fibres and orphans. Equations `(7)--(9)` give

```text
boxed:
c_e(T)-c_e(T-v)
 =sum_(P in E_e(T-v))(r_v(P)-m_(v,e)(P))+O_(v,e).          (10)
```

Combining `(4)` and `(10)` yields the exact compensation inequality

```text
sum_P r_v(P)+O_(v,e) >= sum_P m_(v,e)(P).                 (11)
```

Marked-slot theft can occur pointwise, but redundant insertions and orphans
pay its total invoice. This is the marked-edge refinement of THM-4094's full
deletion-incidence deficit.

## 4. Quadratic restriction deck

Now let `z=(z_e)` be any symmetric edge vector on the fixed oriented complete
graph underlying `T`. Define

```text
d_i(z)=sum_(j!=i)z_ij,
h_i(z)=sum_(i->j)z_ij-sum_(j->i)z_ij,
C(z)=sum_i h_i(z)d_i(z),
D(z)=sum_(e<f, e intersect f=empty)z_e z_f.               (12)
```

Direct expansion cancels every edge square and every mixed in/out pair at a
vertex:

```text
C(z)=2(sum_(e<f with common tail)z_ez_f
       -sum_(e<f with common head)z_ez_f).                 (13)
```

Every monomial of `C` is supported on three vertices, while every monomial
of `D` is supported on four. Hence, writing `z|_(V-v)` for restriction,

```text
boxed:
sum_v C(z|_(V-v))=(n-3)C(z),
sum_v D(z|_(V-v))=(n-4)D(z).                              (14)
```

No tournament-specific positivity is used in `(13)--(14)`.

## 5. Parity holonomy law

Assume `n>=5`, `D(z)>0`, and

```text
D(z|_(V-v))>0                    for every v.              (15)
```

The denominator hypothesis is necessary for the normalized statement. For
example, an order-five tensor supported only on two disjoint edges has
positive parent `D` but several zero restriction denominators.

For `k>=4`, define the signed normalized tilt

```text
tau_k(z)= (k-3)C(z)/(2D(z)),       k even,
tau_k(z)= (k-3)C(z)/(4D(z)),       k odd.                  (16)
```

Let

```text
lambda_v=D(z|_(V-v))/sum_u D(z|_(V-u)).                   (17)
```

Substitution of `(14)` into `(16)--(17)` proves

```text
boxed:
tau_n(z)=1/2 sum_v lambda_v tau_(n-1)(z|_(V-v)), n odd;
tau_n(z)=2   sum_v lambda_v tau_(n-1)(z|_(V-v)), n even.  (18)
```

The weights are positive and sum to one. Condition `(15)` holds for the
actual tournament capacity tensor `z=c(T)`: THM-4128 gives positive
disjoint-edge energy. It also holds for every restricted parent tensor when
`n>=5`. Indeed `(4)` gives

```text
c(T)|_(T-v)=c(T-v)+delta^(v),       delta^(v)>=0,          (19)
```

and the left side coordinatewise dominates the actual card tensor, whose
disjoint-edge energy is positive by THM-4128. In particular,
`D(c(T)|_(T-v))>0` for every `v`.

## 6. Exact order-eleven reframe

At `11 -> 10`, set

```text
b^(v)=c(T)|_(T-v)=c(T-v)+delta^(v).                       (20)
```

These are the **holonomy-corrected card tensors**, not the actual card
capacities. Equation `(18)` becomes

```text
tau_11(c)=1/2 sum_v lambda_v tau_10(b^(v)).                (21)
```

THM-4128 says that every rational support-floor optimizer is central under
the strict gate `|tau|<1`. Equality can tie a central layer with an outer
layer; the non-strict inequality asserts only the existence of a central
optimizer. Consequently the exact strict order-eleven criterion is

```text
boxed: |sum_v lambda_v tau_10(b^(v))|<2.                  (22)
```

The triangle inequality gives the stronger sufficient target

```text
sum_v lambda_v |tau_10(b^(v))|<2.                         (23)
```

Equation `(23)` is sufficient, not equivalent to centrality. It is still far
weaker than asking every corrected restriction to satisfy the order-ten wall
`|tau_10|<1`.

The support moment also reconstructs the complete parent tensor:

```text
A+Omega=9c.                                               (24)
```

Thus homogeneous Johnson ratios and exchange directions can be recovered
from the labelled deletion tensor plus `Omega`. Absolute floors are not
claimed to scale unless their base term is scaled as well.

## 7. Boundaries, hostiles, and replay

Primeness does not remove the correction and does not make every card strong.
The exact prime order-eleven code `3169369058263173` has ten strong cards,
but deletion `v=10` is nonstrong and fails the order-ten rational gate:

```text
parent: H=23685, (C,D)=(4220068008,88725253576),
        rho=1055017002/11090656697;
card:   H=2037,  (C,D)=(755197384,1016215996),
        7|C|=5286381688 > 2D=2032431992.                  (25)
```

Conversely, the THM-4133 order-twelve counterexample shows why the parity
factor cannot be discarded. Every restricted parent tensor is individually
central at the odd scale, with

```text
max_v |tau_11|=9073595176/12026131621<1,
sum_v lambda_v tau_11=-53092739331/80871049732.            (26)
```

The even-order factor `2` in `(18)` gives the failing parent load

```text
tau_12=-53092739331/40435524866.                           (27)
```

The primary exact audit is

```text
python3 -B 04-computation/tournament_capacity_deletion_holonomy_thm4167.py
python3 -B -O 04-computation/tournament_capacity_deletion_holonomy_thm4167.py
PYTHONHASHSEED=271828 python3 -B \
  04-computation/tournament_capacity_deletion_holonomy_thm4167.py
```

Compile the independent exposed-word/ear-response referee with

```text
clang++ -O3 -std=c++17 \
  04-computation/tournament_capacity_deletion_holonomy_thm4167_independent_audit.cpp \
  -o /tmp/thm4167_independent_audit
/tmp/thm4167_independent_audit
```

An `-O0` Clang build and an `-O3` GCC build produce the same frozen output.
The theorem is all-order; the finite audits are hostile controls rather than
the source of the proof. **QED.**
