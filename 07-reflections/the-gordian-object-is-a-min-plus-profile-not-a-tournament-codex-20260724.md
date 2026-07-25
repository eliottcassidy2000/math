---
source: codex-2026-07-24-knot-relations
status: >
  SYNTHESIS. The proved mathematical outputs are THM-2176, THM-2182, and
  THM-2183. The external
  transfers below are separated into cited facts, exact repo analogues,
  executable proposals, and no-go statements. No knot-to-LRC or
  knot-to-tournament equivalence is asserted.
related:
  - THM-082
  - THM-840
  - THM-853
  - THM-1862
  - THM-1975
  - THM-2162
  - THM-2174
  - THM-2176
  - THM-2182
  - THM-2183
---

# The Gordian object is a min-plus profile, not a tournament

## 1. Inheritance pass

The closest proved mechanism was not in knot theory. It was THM-840's exact
criterion

```text
ker(O) subset ker(O after T)
```

for an observation to update under an operation, followed by THM-853's finite
Myhill--Nerode refinement of a recursive state. The closest already-realized
knot transfer was THM-2162: retain a whole safe core, integrate its signed
boundary, and only then charge the incoming comb.

The canonical hostile example was the old knot--tournament tangent T225.
Reversal of a directed 3-cycle preserves the Hamiltonian-path count at `n=4`
but fails at `n=5`. A skein-shaped identity or a Reidemeister-shaped move is
not enough to create a knot invariant.

The corrected near miss was MISTAKE-245. A bounded LRC bank missed the
structured replacement `10->20`; THM-2162 repaired the claim with an infinite
analytic tail and the exact one-swap gap. This is the same failure genus as
connected-sum nonadditivity: separately optimized or bounded local pieces can
miss a thin common-intermediate shortcut.

The least-used sidecar was THM-1975's path-cover profile. Its finite verified
universe shows the right shape: scalar `H` does not compose under cyclic
substitution, while the full path-system profile does.

## 2. Concept board

| object / representation | predicate | operation | lost coordinate | cheapest test |
|---|---|---|---|---|
| Gordian graph / min-plus kernel | crossing-change distance | connected sum | intermediate knot if evaluated only at the unknot | recover the `7_1` shortcut from a landmark |
| continuation profile `Pi_K(J)` | future unknotting cost | append `J` | diagrams and isotopy witnesses | split `7_1` from its mirror using continuation `7_1` |
| defect cocycle `sigma(K,J)` | strict subadditivity | merge sums | additive one-body gauges | enforce the triple cocycle identity |
| rooted Gordian order | membership in a minimal unknotting path | descend toward the unknot | incomparable branches | test product-cone monotonicity |
| tournament order-join | edit response and compositional invariants | `A▷-` | ambient quotient shortcuts | apply THM-2183's block-mixing uncrossing |
| LRC safe-core profile | future safe measure/truth | insert combs | joint normalized-tail overlap | compare THM-2182's `{3,12}` and `{4,6}` tails |
| inherited-color amplitude | whole-fiber interference | sum over perfect matchings | matching identity after scalarization | compare full amplitude with a tournament projection |

## 3. What THM-2176 changes

For knots, the scalar

```text
u(K)=d_G(K,U)
```

is one matrix entry of the faithful min-plus kernel

```text
P_K(A,B)=d_G(K#A,B),
P_K tensor P_L=P_(K#L).
```

Equivalently,

```text
u(K#L)=inf_A [d_G(K,A)+u(L#A)].
```

The separate-summand bound uses only the landmark `A=U`. A common-intermediate
search is therefore a restricted tropical matrix product.

The universal operation-ready refinement is

```text
Pi_K(J)=u(K#J).
```

Equality of these profiles is the largest connected-sum congruence contained
in equality of `u`. This statement is stronger and more precise than saying
that the profile "contains more information": it is the coarsest possible
refinement on which every future connected sum updates deterministically.

The first splitter is forced. For `K=T(2,7)`,

```text
u(K)=u(mirror K)=3,
Pi_K(K)=6,
Pi_(mirror K)(K)<=5.
```

Thus mirroring is invisible to the one-body scalar here and visible to a
single continuation query.

The interaction

```text
sigma(K,L)=u(K)+u(L)-u(K#L)
```

is a nonnegative symmetric normalized 2-coboundary. Its cocycle identity

```text
sigma(K,L)+sigma(K#L,M)
 =sigma(L,M)+sigma(K,L#M)
```

is the exact associativity law for all many-body savings. The binary
"symbiont" relation from Brittenham--Hermiller's 2026 follow-up is only the
support `sigma>0`.

The rooted relation

```text
A <=_G K  iff  u(K)=d_G(K,A)+u(A)
```

is a graded partial order. The defect is monotone on its product cones, so
strict pairs form an up-set. This proves the abstract poison-pair propagation
used in the torus-knot families without treating it as a special diagram
trick.

Finally, the exact split

```text
sigma(K,L)
 = [u(K)-d_G(K#L,L)]
 + [d_G(K#L,L)+u(L)-u(K#L)]
```

separates translation catalysis from bypass of the separated-summand
intermediate. Signature forces the first bracket to vanish for
`T(2,7)#mirror(T(2,7))`. The known counterexample is a folded Gordian square:
all four separated sides retain length three, but the root-to-sum diagonal
has length at most five instead of six.

Connected-sum homogenization does not repair the scalar. Fekete gives a
homogeneous limit, but signature fixes the two one-body limits at three while
the repeated five-change certificate keeps the mirror-pair limit at most
five.

The wider published atlas sharpens the mechanism. The `5_1,8_2` seed and all
`tau`-calibrated torus/Baader families are pure bypass in both directions.
The `4_1,9_10` seed is pure bypass in the figure-eight direction, while its
reverse translation term is unresolved. The candidate `3_1,10_6` seed is
conditional on `u(10_6)=3`. Across the two audited papers, no example proves
positive translation catalysis. That is now the genuinely new knot target.

There is also a clean quotient interpretation. Give ordered pairs of knots
the `l1` product Gordian metric, then forget the two summand labels by
`(K,L)|->K#L`. The contraction of root distance is exactly `sigma(K,L)`.
Thus the interaction defect measures information lost by forgetting the
decomposition marker, not merely failure of one numerical formula.

## 4. Four relation algebras, only one of them a tournament

The supplied sources become coherent after naming their coefficient algebra.

| world | carrier | aggregation | composition | faithful information |
|---|---|---|---|---|
| knots | Gordian cost kernel | `min` | `+` | optimal intermediates and costs |
| inherited colors | perfect-matching amplitude | complex `sum` | complex `product` | destructive/constructive fiber interference |
| LRC | signed endpoint/Fourier current | signed `sum` or integral | insertion/convolution | owner, phase, scale, magnitude |
| tournaments | complete asymmetric Boolean relation | `or` | relational composition | one orientation per pair |

Krenn--Gu--Soltesz assign each perfect matching an inherited vertex coloring
and sum matching weights over the whole coloring fiber. Projecting this to a
pairwise orientation discards edge absence, multiedges, endpoint colors,
matching multiplicity, and complex cancellation. The object is a weighted
relation tensor, not a tournament.

Zakharov's word-overlap relation likewise allows both directions, neither
direction, and loops. Its theorem concerns a forbidden directed bipartite
rectangle. The proof splits first-hit layers, derives a convolution inequality,
chooses the positive root of a generating polynomial, and optimizes the root.
Tournament completion would erase the forbidden-rectangle predicate.

The general lesson is not that tournaments are too weak. It is that a
tournament is the correct terminal quotient only when the source supplies an
intrinsic complete antisymmetric comparator.

## 5. What the other supplied papers contribute

### Structured multiplication: search in the correct quotient

Rybin--Zhang--Luo compute `XX^t` more cheaply than treating the two factors as
generic. Their discovery loop is:

```text
generate rank-one candidates
 -> enumerate exact linear relations
 -> solve a minimum cover of every target expression.
```

This is a concrete search architecture for two repo problems.

For knots, generate a certified landmark dictionary, compute or bound
`d_G(K,A)` and `u(L#A)`, and cover target pairs by short common-intermediate
certificates. The mirror involution is structural input, analogous to
transpose symmetry, not an after-the-fact label.

For LRC, generate low-support signed relation packets, verify their exact
rational span, and solve a minimum cover of all deletion/owner obligations.
The output is only algebraic until phase owner, scale, and analytic-uniformity
sidecars pass hostile controls. THM-2174 is the canary: extensive algebraic
and finite phase state can survive while the target measure changes.

### Word overlap: first-hit layers rather than pairwise orientation

Zakharov proves the sharp asymptotic product bound

```text
mu(A)mu(B) <= (1/n)(n/(n+1))^(n+1) ~ 1/(en)
```

for two length-`n` word families with no suffix-prefix overlap. A legitimate
LRC transfer would require a cyclic owner word whose first-hit layers are
actually disjoint and whose shift densities obey the same convolution.
Cutting a cyclic word at an arbitrary point loses rotation, reflection, and
endpoint owner. The cheapest test is therefore finite and exact: construct the
shift layers on one named denominator packet and check the convolution before
borrowing the inequality.

### Inherited colors: aggregate the full fiber

The matching amplitude

```text
A_G(c)=sum_(matching M inducing c) product_(e in M) w_e
```

is the sum-product analogue of THM-2022's entire balanced face and THM-2162's
entire signed endpoint boundary. Atomwise positivity or an unsigned norm is
orthogonal to cancellation over a fiber. The transfer is methodological:
choose the fiber before summing. There is no proved matching model of Gordian
distance or LRC owner currents here.

### The anonymous shared conversation

The supplied [public ChatGPT conversation](https://chatgpt.com/share/6a60b2eb-0b64-83ee-9c76-7931ca1de063)
was recovered in a signed-in browser after the first source pass. It concerns
Goemans' cost-enhanced single-source unsplittable-flow conjecture, not the
classical Dinitz list-coloring conjecture. Its seven-vertex construction was
then rebuilt from the graph, checked against Conjecture 1.3 of
Traub--Vargas Koch--Zenklusen, generalized to a nonempty open
three-parameter region, and independently exhausted at graph level in
[THM-2177](../01-canon/theorems/THM-2177-planar-counterexample-to-goemans-unsplittable-cost-flow-conjecture.md).

The native binary carrier is again not a tournament. Three cheap path choices
have the symmetric conflict graph `K_3`; their fractional marginals sum to
`16/15`, violating the triangle stable-set facet. An orientation would erase
the independent-set predicate. The faithful quotient is the conflict graph
together with its fractional point and complementary cost separator.

## 6. Tournament consequences

Under tournament order-join,

```text
H(A▷B)=H(A)H(B)
```

by THM-1862. Thus `log H` has zero interaction cocycle and its continuation
profile collapses to the scalar. This is the opposite of unknotting number.

Under cyclic substitution, THM-1975 gives the relevant warning on its verified
finite universe: scalar `H` is not enough, while the path-cover profile is.
The open general transfer kernel should be viewed as a tournament continuation
profile, not merely as "one more invariant."

The metric question has a sharper exact answer. THM-2183 proves, in the
unmerged isomorphism quotient,

```text
d_iso(A▷X,B▷Y)=d_iso(A,B)+d_iso(X,Y).                 (1)
```

A two-vertex image-swap uncrossing strictly improves every block-mixing
bijection. Thus fixed order-join patches are isometric Cartesian products and
have no analogue of knot translation catalysis, even after the cut marker is
forgotten.

The ambient quotient still has shortcuts: its four-vertex tournament layer
contains a triangle and is not a partial cube. The distinction is therefore
operation-local. Relabelling creates global bypasses, but the constant one-way
cut prevents them from entering a fixed join patch. Knots lack precisely this
Monge uncrossing once the connected-sum sphere is forgotten.

## 7. LRC consequence and next targets

The common-intermediate transfer has already been realized by THM-2162; it
should not be advertised as new. The new contribution is the continuation
criterion. For a finite speed core `E`, define

```text
Pi_E(F)=target(E union F)
```

for the chosen target: exact safe measure, weak nonemptiness, or a named
certificate. Equality of full profiles is the universal insertion congruence.
Any proposed finite carrier must prove that its kernel lies inside this
future equivalence.

THM-2174 is the first hostile splitter. Endpoint phase labels, even together
with substantial relation and owner data, do not determine the future measure
because the current retains a `1/W` scale.

THM-2182 now gives the positive whole-profile counterpart. If `L` is a
positive multiple of every core endpoint denominator, then

```text
measure(G_(E union LC))=measure(G_E)measure(G_C).     (2)
```

The map `t=(j+x)/L` preserves the entire normalized tail product, not its
marginals. It closes every aligned seven-core/six-tail row with the explicit
floor `1/(49e_2)`. Conversely, normalized tails `{3,12}` and `{4,6}` have the
same phase-zero word, zero one-comb currents, equal one-comb masses, and equal
reciprocal sum, yet their joint masses differ by `measure(G_E)/84`. This
isolates joint ratio/gcd overlap as an indispensable interaction sidecar.

The next concrete portfolio is:

1. **Knot anchor:** search specifically for positive translation catalysis
   `C_J(K)>0`; all currently calibrated published families are pure bypass.
   Use the certified dictionary `{U,7_1,5_1,8_2}` as the first continuation
   table and keep interval values where exact unknotting numbers are unknown.
2. **Tournament niche:** classify which other modular/substitution cuts admit
   THM-2183-style Monge uncrossing; the ambient quotient triangle is the
   mandatory hostile control.
3. **LRC wildcard:** combine THM-2178's mod-14 rank fork and THM-2179's
   body-touching packet with THM-2182's exact product. The decisive supplier
   would force or manufacture an endpoint-aligned `7+6` split without losing
   distinctness.
4. **Word-overlap hostile probe:** test the first-hit convolution on one cyclic
   endpoint-owner packet; stop immediately if rotations make the layers
   non-disjoint.
5. **Fiber-amplitude probe:** encode one existing finite signed LRC packet as a
   symbolic whole-fiber sum before attempting any matching interpretation.

The cross-domain rule is:

```text
compose in the native relation algebra first;
evaluate a scalar or tournament shadow last.
```

## Primary sources used

- Brittenham--Hermiller,
  [*Unknotting number is not additive under connected sum*](https://arxiv.org/pdf/2506.24088).
- Brittenham--Hermiller,
  [*Unknotting number and connected sums: The knots 4_1 and 5_1*](https://arxiv.org/pdf/2601.18757).
- Zakharov,
  [*An isoperimetric inequality for word overlap*](https://arxiv.org/pdf/2602.20143).
- Rybin--Zhang--Luo,
  [*XX^t Can Be Faster*](https://arxiv.org/pdf/2505.09814).
- Krenn--Gu--Soltesz,
  [*Questions on the structure of perfect matchings inspired by quantum physics*](https://arxiv.org/pdf/1902.06023).
- Traub--Vargas Koch--Zenklusen,
  [*Single-Source Unsplittable Flows in Planar Graphs*](https://arxiv.org/pdf/2308.02651).
