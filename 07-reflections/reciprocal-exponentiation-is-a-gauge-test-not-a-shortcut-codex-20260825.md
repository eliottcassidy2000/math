---
source: codex-lrc14-abc-exponent-reciprocity-20260825
status: >
  RESEARCH SYNTHESIS. THM-4105--4107 are unconditional proved results;
  THM-4108 is unconditional for its separation and gauge obstruction and
  conditional on ABC for radical saturation. LRC(14), ABC, and the disputed
  IUT-to-ABC implication remain OPEN / CONTESTED at their recorded scopes.
external:
  - https://arxiv.org/abs/2604.23906
  - https://link.springer.com/article/10.1007/s40315-021-00369-6
  - https://arxiv.org/abs/2111.07791
  - https://www.math.uni-bonn.de/people/scholze/WhyABCisStillaConjecture.pdf
---

# Reciprocal exponentiation is a gauge test, not a shortcut

## Verdict

The useful hidden key in `a^b` versus `b^a` is not the direction of the
inequality. Raw direction is only the scalar potential `log(x)/x`. The useful
structure appears when exponentiation is asked one of three better-typed
questions:

1. do reciprocal phase equations descend to one physical time?
2. what edgewise normalization created a tournament cycle?
3. after primitive normalization, what radical information survives an
   additive packet?

Those questions produced four durable results:

- THM-4105: reciprocal commutators exactly characterize a primitive row's
  physical phase orbit and quantitatively control arrival;
- THM-4106: pair optimum plus directed owner mass recovers every mixed-parity
  primitive ratio, and a cross-`v_2` tree recovers a whole row;
- THM-4107: pairwise gcd normalization creates a sparse exponent tournament
  whose curl is exactly divisibility/`2:3` gauge incompatibility, but the
  complete labelled tournament is LRC-blind;
- THM-4108: ABC conditionally forces logarithmic radical saturation of
  coprime `a^b+-b^a`, while the natural LRC straddle packet fails an exact
  lift-gauge test and cannot force the loneliness margin.

The net LRC advance is a sharper representation and a smaller named
obligation: **simultaneous phase/arrival synchronization across a lossless
reciprocal-owner tree**. It is not a proof of LRC(14).

## Inheritance pass and live board

The closest proved mechanism was THM-4095's exact two-speed optimum and order
blindness. The canonical hostile was AP13, supplemented by V26 and the
nonprimitive phase packet `(2,4);(0,1/2)`. The corrected near miss was the old
hope that radical support or tournament orientation could stand in for owner
and arrival. The least-used sidecars were directed Haar owner mass and the
full physical one-parameter orbit.

The board stabilized as follows:

| object | representation | invariant | operation | destroyed information |
|---|---|---|---|---|
| primitive speed row | reciprocal phase commutators | exact orbit membership | Bezout descent | nothing at defect zero |
| two-speed observer | sum/difference clocks + owner cells | primitive ratio in mixed parity | Haar integration | common gcd sheet |
| reciprocal powers | signed prime-valuation defect | primitive powered ratio | positive/negative part split | radical projection loses depths/slots |
| exponent tournament | edgewise gcd chart | dilation class + sparse reversal type | contextual normalization | residue, margin, phase, arrival |
| LRC straddle | `(q,p,D)` determinant packet | local signed margin | integer lift change | ABC height/radical are not gauge invariant |

Every new idea was compared against all five rows. That prevented three
tempting but false transfers: scalar power order to tournament curl, pair
reconstruction to common arrival, and ABC height control to an LRC margin.

## Surprise 1: reciprocal exponent equations are physical descent

For a primitive integer vector `v`, the equations

```text
z_i^(v_j)=z_j^(v_i)
```

look pairwise but are actually global. A Bezout word reconstructs the unique
common base `x`, so on the circle they say exactly that

```text
theta_i=v_i t mod 1
```

for one physical time. This is the right exponentiation analogue of a
rank-one Pluecker condition: all two-coordinate commutators vanish, and
primitivity removes the residual torsion sheet.

The quantitative version is equally revealing. Maximum commutator defect is
bi-Lipschitz to distance from the orbit up to explicit constants
`1/(2 max v_i)` and the minimum Bezout `l1` tariff. The next LRC compiler can
therefore distinguish an algebraically compatible packet from a physically
arriving one without pretending that arrival already implies loneliness.

## Surprise 2: the owner imbalance is a reciprocal Pythagorean leg

For primitive `m<n`, the owner switch is governed by the two clocks

```text
n-m,                    n+m.
```

Their square-wave correlation is zero in the odd shell and
`1/(n^2-m^2)` in mixed parity. Hence the smaller-speed owner excess is

```text
1/[2(n^2-m^2)].
```

Meanwhile the optimized pair defect is `1/[2(n+m)]`. Dividing the two
recovers `n-m`; the optimum recovers `n+m`. Addition and subtraction decode
the pair, while their product is the first leg of

```text
(n^2-m^2,2mn,n^2+m^2).
```

This is the cleanest bridge found in the session to the repo's Berggren and
Pythagorean carriers. Its boundary is equally useful: the decoder collapses
every odd pair to `(1/2,1/2)`, but an all-odd whole row is already safe at the
antipode. In a mixed row, the unequal-`v_2` graph is connected, so only twelve
informative pair ratios are needed to reconstruct thirteen speeds.

The forgotten coordinate is not arithmetic. It is **when the twelve edges
can be realized together**.

## Surprise 3: cycles measure normalization holonomy

Raw `a^b ? b^a` cannot cycle because it sorts `log(a)/a`. Pairwise gcd
normalization changes the object. For `a<b`, the normalized smaller speed
loses exactly when

```text
a divides b                 or                 3a=2b.
```

On a sorted triple, cyclicity occurs exactly when the two adjacent exception
bits agree and the long-edge bit differs. This is a genuine tournament, but
its curl comes from charts that cannot be glued to one vertex potential.

That explains several older tournament lessons in one line:

```text
global vertex normalization -> gradient / transitive;
edge-dependent normalization -> possible holonomy / cycles.              (1)
```

The LRC hostile is decisive. AP13 and every `{1,...,12,p}` with prime `p>13`
have the same complete labelled normalized exponent tournament. AP13 has
margin `1/14`; every companion row has the explicit `1/13` witness. Thus
even genuine curl is not enough unless the normalization preserves the target
predicate.

The out-of-box survivor is the `2<->3` exchange graph. Moving from classic
power level one to every higher integer level flips exactly the edges
`{2g,3g}`. Its components are valuation paths preserving the prime-to-six
part and total `v_2+v_3`. This deserves comparison with the existing mod-six
LRC carriers, but only after residues and owners are attached.

The incoming THM-4099 gap-transfer theorem sharpens this tournament-only lane.
Because every edge of the normalized exponent tournament is an explicit
divisibility or `2:3` test, its squarefree local insertion factors are
arithmetic predicates rather than black-box orientations. That is promising
for exact recurrences of `H` on initial segments and connects directly to
THM-4102's selected-ear bank. It cannot improve LRC by itself: the complete
gap polynomial is a function of the labelled tournament, so it is identical
on the AP13/prime hostile fibre.

## Surprise 4: ABC is strongest where the LRC bridge is weakest

For coprime distinct bases, the unconditional separation

```text
|a^b-b^a| >= max(a^b,b^a)/9
```

makes the sum and difference comparable to the powered height. Assuming ABC,
their radicals consequently have logarithmic exponent tending to one. This
is a real conditional theorem, stronger in exponent than fixed-signature
power estimates because the bases grow only logarithmically in the height.

But the natural LRC straddle packet fails before ABC can help. Replacing a
circle representative `t` by `t+N` adds the same `Nuw` to both additive
terms while leaving the determinant, margin, and physical phases fixed.
Hence ABC height is a coordinate of the chosen integer lift, not an invariant
of the LRC state. AP13 is the excluded boundary `0+1=1`, and the exact family

```text
(u,w,alpha,beta)=(n+1,n,1,1)
```

has determinant one but unbounded `q/D=2n+1`. ABC is perfectly compatible
with these packets.

The correct use of ABC is therefore downstream of an actual positive
three-term relation and upstream of a separate owner selector. Radical size
is a budget, not a residue witness.

## Challenged assumptions and repaired forms

| tempting assumption | minimal hostile | strongest survivor |
|---|---|---|
| ordered `a^b,b^a` creates tournament structure | scalar potential `log x/x` | edgewise gcd normalization creates typed holonomy |
| pairwise compatible phases are physical | nonprimitive `(2,4);(0,1/2)` | primitive commutators iff one physical time |
| lossless pair ratios solve the row | `{1,3}`/`{1,5}` with different third-runner optima | cross-`v_2` tree reconstructs arithmetic; locations remain |
| a rich exponent tournament predicts loneliness | AP13 versus `{1,...,12,p}` | tournament is a diagnostic only with residue/arrival sidecars |
| ABC bounds the LRC straddle margin | lift gauge and `(n)+(1)=(n+1)` | conditional radical gate on a fixed primitive packet |
| large powered height means large primitive ABC height | rational equality Stern ray | signed valuation defect records what survives gcd removal |

## Procedural generator for the next session

The best next tasks arise by crossing objects, representations, operations,
and missing sidecars:

1. **Object:** live `11+2` row. **Representation:** cross-`v_2` owner tree.
   **Operation:** choose a low-clock tree. **Target:** common safe cell.
   **Hostile:** incompatible edge maximizers.
2. **Object:** normalized exponent tournament. **Representation:** two-colour
   reversal deck. **Operation:** add one modulus. **Target:** a signed blocker
   card. **Hostile:** AP13/prime labelled twins.
3. **Object:** actual short Graver relation. **Representation:** primitive
   ABC packet plus valuation slots. **Operation:** split low/rich radical.
   **Target:** owner-resolved prime selection. **Hostile:** lift gauge.
4. **Object:** rational equality Stern ray. **Representation:** signed
   valuation defect. **Operation:** one-coordinate perturbation. **Target:**
   a nontrivial primitive packet. **Stopping reason:** common gcd still absorbs
   exponential height.
5. **Object:** reciprocal commutator matrix. **Representation:** sparse
   Bezout anchors. **Operation:** minimize `B(v)`. **Target:** quantitative
   arrival within the current 165 valuation orbits. **Hostile:** rows needing
   three anchors such as `(6,10,15)`.

The session's reusable research move is not “try exponentiation.” It is:

```text
first test whether the comparator is a global potential;
then identify the contextual gauge that creates curl;
finally audit whether that gauge preserves the target predicate.           (2)
```

That move produced both the new carrier and its exact stopping theorem.

## Scope

LRC(14) is **OPEN**. ABC is **OPEN**. The IUT-to-ABC implication remains
**CONTESTED** in the repo ledger. No theorem here proves irrationality,
transcendence, Mahler's `3/2` problem, the Jacobian conjecture, or a global
tournament inequality. The advance is a set of proved coordinates, exact
firewalls, and a more precise synchronization frontier.
