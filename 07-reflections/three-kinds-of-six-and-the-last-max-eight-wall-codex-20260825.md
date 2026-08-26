# Three kinds of six and the last maximum-eight wall

**Session reflection / provenance, 2026-08-25.** This is not a truth source.
For current statements use
[THM-4139](../01-canon/theorems/THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier.md),
[THM-4143](../01-canon/theorems/THM-4143-two-term-collision-wall-critical-boundary-monodromy-exclusion.md),
and the maintained navigation files. External-source status and transfer
firewalls are recorded in
[the source ledger](../05-knowledge/reference/CORE-PAPERS-PETERSEN-INSTANTON-EXACT-SIX-2026-08-25.md).

## Portfolio and inheritance

The session used the required three-lane portfolio.

- **Anchor:** squeeze the last open collision wall in THM-4053's exact-`M=8`
  planar-Jacobian trichotomy.
- **Niche:** completely classify rational preperiodic points of
  `x^2-29/16` and look honestly for period six.
- **Wildcard:** separate the Mersenne exception at `63`, Petersen-minor flow
  structure, and instanton passage asymptotics into typed mechanisms that
  might suggest lawful new observers.

The closest inherited mechanism was the fixed-sheet/orbit-merger obstruction
of THM-4138/4141. The hostile example was THM-4134's repeated edge, where
naive specialization loses roots. The corrected near miss was an arbitrary
boundary-response allocation that ignored Galois labels. The least-used
sidecar was the identity `P^4-PY^2=TP^3` on the two-term wall.

The live concept board was:

1. affine critical-scheme length;
2. normalized boundary packets with residue fields;
3. quotient loss under a central sign or a collided branch;
4. exact-period assembly by prime-power depth versus CRT;
5. obstruction structure under graph minors;
6. rank-one leading terms plus controlled defects.

## What “six” meant in each lane

The same integer appeared through three genuinely different carriers.

### Rational dynamics: no period six

For

```text
f(x)=x^2-29/16,
```

the complete rational preperiodic set is

```text
{+-1,+-3,+-5,+-7}/4.
```

Its only rational cycle is

```text
-7/4 -> 5/4 -> -1/4 -> -7/4.
```

Thus the requested rational period-six search terminates negatively for this
specific parameter, by an all-height denominator-and-escape proof rather than
a bounded search. The unique projective map cyclically interpolating these
three points has an `SL_2` lift `B` with

```text
B^3=-I,                         B^6=I.
```

This is order six in a central lift, not a six-cycle of the quadratic map.
The quotient `SL_2 -> PGL_2` destroys the sign and reduces the visible order
to three.

### The `3-4-5` coordinate is exact but local

The marked three-cycle chart has arithmetic-progression representatives
`t=1,-2,-1/2`. At `t=1`, Euclid's pair `(t+1,t)=(2,1)` produces

```text
(3,4,5),        D(1)=4,        denominator(c)=16=4^2,
29=5^2+2^2,     7=4^2-3^2.
```

This explains every requested number in one coordinate specialization. It
does not turn Pythagorean triples into a general cycle-construction theorem.
The preserved data are the marked AP and its relabeling action; the destroyed
data under forgetting the marking are the chosen Euclid pair and orientation.

### Mersenne `63`: no primitive prime does not mean no new period

At exponent six,

```text
2^6-1=63=3^2*7,
ord_3(2)=2,                 ord_7(2)=3,
```

so there is no primitive prime divisor. Nevertheless

```text
ord_9(2)=6,                 ord_21(2)=6.
```

These are two different mechanisms: prime-power depth at `9`, and CRT order
assembly at `21`. For squaring on `G_m`, the exact-period-six polynomial is

```text
Phi_9 Phi_21 Phi_63,
```

with `54` exact-period-six points in nine cycles. The modulus `63` is their
join, not their unique primitive carrier. A Boolean “period six” observer
forgets the conductor vector and therefore conflates `9`, `21`, and `63`.

For `x^2-29/16`, exact period six does occur after reduction: modulo `43^2`
there are `126` exact-period-six points, or `21` cycles, in the tube above the
rational three-cycle. Modulo `63` there are only three-cycles and no
six-cycle. This makes the local/global distinction concrete: a modular
period-six lift is neither a rational period-six orbit nor a proof that one
exists.

## The useful transfer: collision labels, not numerical coincidence

The dynamics did not supply a Jacobian counterexample. Its useful lesson was
typed: a quotient can merge distinct objects while a small sidecar restores
them. In the dynamical chart the missing sidecar is the central sign or the
marked-cycle parameter. On the Jacobian wall it is the residue-field label of
each normalized puncture.

That distinction repaired the Anchor. On `Theta=-Delta`, the top row
contracts exactly, and the final wall splits into:

```text
Phi=0: impossible by the critical-value cubic;
Phi!=0,K!=0: critical length 19, packet (6,6,3,2,2,1);
Phi!=0,K=0:  critical length 17, packet (6,6,5,1).
```

The generic packet contains rational entries `(6,6,3,1)` and one irreducible
quadratic closed point `(2,2)`. That label forces the exact response degrees
`20` or `16`; the `K=0` packet is rational and forces degree `18`. The full
responses contradict the one-pivot commutator shape. The finite quadratic
response contributes exactly two transpositions, one orbit-merger unit too
few to overcome the nineteenth critical fixed sheet. This closes the wall
relative to the inherited target and transport theorems.

The comparison that survived is therefore:

```text
quotient/normalization  +  lost label  +  exact lift/response sidecar.
```

The superficial comparison “both involve six” did not survive.

## What the two new preprints contributed—and did not

Pintér's Petersen-contraction theorem supplies a legitimate obstruction
pattern: every bridgeless `P/e`-minor-free multigraph has a nowhere-zero
four-flow. Its vertices and relation are graph branch sets and minor
adjacency, not Jacobian punctures or a tournament. Minor containment forgets
the embedding and boundary multiplicities, so no planar-Jacobian implication
was imported. A potentially lawful future test is to ask whether contraction
of a stable-fibre dual graph preserves the **monodromy-index deficit** together
with a labelled boundary sidecar. The cheapest hostile is `K_10`, which
contains the minor but still has a four-flow.

Shapiro's double-well theorem isolates another useful proof interface: a
rank-one leading kernel plus defect estimates and a summable majorant yields
a Poisson passage law at fixed mean. The required objects are a positive path
measure, spectral gap, kernel decomposition, and uniform domination. A
deterministic cycle and a Jacobian covering have none of these by default.
The only retained inspiration is methodological: never pass from fixed-sector
asymptotics to a normalized global law without an escaping-sector bound.

Neither preprint proves or even directly advances `JC(2)`. Their value here
was to sharpen the sidecar and domination gates.

## Precise next frontiers

The exact-`M=8` result changes the next session's geometry. Re-running another
collision wall is now saturated. The following objects are sharper targets.

1. **Exact weight nine on the same seam.** The new top monomials solving
   `2i+3j=9` are `p^3y` and `y^3`. Build the complete forced polynomial,
   enumerate its support/collision strata, and for each stratum compare
   critical length with normalized boundary defect and residue-degree labels.
   The cheapest decisive test is a symbolic support-and-resultant compiler
   with a hostile leading-row specialization.
2. **A response theorem for residue degree at least three.** The quadratic
   carrier was rigid because one base change led to a rank-one rational
   elliptic surface. Classify what survives for cubic and quartic boundary
   closed points: source point, residue extension, horizontal image, lost
   Galois data, and polynomial-section sidecar must all be explicit.
3. **Entry into the reduced `2:3` seam.** Closing a cell is not an entry
   theorem. Identify the exact operation sending a hypothetical planar
   Keller pair to the cell, what degree/width it preserves, and which
   counterexamples evade it. This remains more globally valuable than one
   more internal coefficient case.
4. **Dual-graph contraction with index labels.** Test whether a labelled
   stable dual graph has a finite forbidden-minor obstruction for failure of
   the orbit-merger budget. Vertices should be components or puncture blocks,
   edges should be actual node attachments, and labels must retain residue
   degree and ramification index. An unlabelled minor statement is too lossy.
5. **Arithmetic dynamics beyond the fixed parameter.** For rational exact
   six-cycles of general `x^2+c`, keep Stoll's conditional global conclusion
   separate from unconditional finite and local results. A feasible exact
   lane is to classify prime-power lifts of the `-29/16` three-cycle by the
   multiplier of `f^3`, explaining why `43^2` doubles the period and which
   primes can do so.

The session's strongest meta-observation is not promoted to canon yet: a
support collision can create new depth without a new prime or a new monomial.
At `63`, depth and CRT assembly restore order six; on the Jacobian wall,
normalization restores two branches after a top-row contraction. These are
evidence for a “retain valuation/label depth after support collapse” move,
but a third independent theorem-level instance is needed before adding it to
`META-PATTERNS.md`.
