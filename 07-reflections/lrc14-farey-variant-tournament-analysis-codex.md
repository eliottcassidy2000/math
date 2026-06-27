# LRC14 Farey-Variant Tournament Analysis

Source: codex, 2026-06-24  
Status: exploratory proof work; new computation added at
`math/04-computation/lrc14_farey_variant_tournament_codex.py`.

## Executive Takeaway

The four Farey mutations do not become a proof by themselves.  The existing
Farey note already gives the warning: optimal witness denominators can grow, so
checking a bounded Farey order cannot close LRC14.

The useful thing is different:

```text
Farey order = perturbation grammar.
mod 7       = sector / Paley-Fano clock.
mod 27      = signed shell clock C = 2n-1 at n=14.
```

After these projections, the four transforms become packet classifiers.  The
new computational signal is that the transform

```text
a/b -> a*b
```

on strict `F_14` creates an induced octahedral shell subgraph in the mod-27
transition graph, on shells

```text
{3,4,6,9,12,13}.
```

That is not decorative: it is the same finite graph type already isolated by
the octahedral current work, namely `L(K4) = K6 minus a perfect matching`.  The
best new hook is therefore:

```text
the numerator-denominator product is seeing the low-height relation product,
and its mod-27 Farey transitions land on the octahedral current carrier.
```

## Current Proof Map

Local notes put the LRC14 proof state roughly here.

1. Non-covering sets are controlled.  If the 13-speed set omits a multiple of
   some `q <= 14`, then `t=1/q` gives a witness.

2. A counterexample must therefore be a covering set: it must meet every
   divisibility class `q=2,...,14`.  That is over-determined.  Local exact
   searches found covering examples far above the `1/14` floor, with minima
   around `0.097`.

3. The row range has been reduced in the cap program to finite rows `k=8..13`.
   The old consecutive-max wall is dissolved.  The live analytic gaps are the
   multi-carrier joint-decorrelation leg and the razor-thin `k=9` Delsarte
   certificate/formalization.

4. Pure tournament `H`, score variance, and Jensen-style energy are too coarse.
   The object must be a metric winding tournament: circular order plus gap
   lengths, anchored at the observer.

This session's additions sit inside item 4.  They are a way to manufacture
finite packet graphs that remember more than scalar energy.

## The Four Farey Variants

For each reduced Farey fraction `a/b` in strict `F_14`, I tested:

```text
P(a,b) = a*b
S(a,b) = a+b
D(a,b) = b^a
N(a,b) = a^b
```

The raw sequence behavior is not decisive.  All four are highly nonmonotone in
Farey order; raw descents and inversions organize intuition but do not provide a
bound.  The meaningful readout appears after sector/shell projection.

Key output from `lrc14_farey_variant_tournament_codex.py`:

| transform | adjacent descents | pair inversions | mod-7 Paley blowup triangles | mod-27 shell edges | induced octahedra |
|---|---:|---:|---:|---:|---:|
| `a*b` | 31 | 496 | 9408 | 39 | 1 |
| `a+b` | 31 | 608 | 10134 | 40 | 0 |
| `b^a` | 31 | 438 | 9234 | 39 | 0 |
| `a^b` | 24 | 441 | 8223 | 26 | 0 |

The distinction is sharp.  `a+b` is the most inversion-heavy and has the largest
Paley blowup triangle count, but it does not produce the octahedral mod-27
carrier.  The exponential variants create more residue collapse: `a^b` has many
mod-27 zero/shell-1 hits and eight adjacent zero differences.  That matches the
older hyperoperation lesson: exponential growth is too lacunary and loses the
dense additive obstruction.

The product transform is the useful one because it is neither purely additive
nor explosively lacunary.  It is the height/product carrier: exactly the kind of
quantity the affine relation-lattice note said Fourier tails see through
coefficient products.

## Why The Octahedron Matters

The existing octahedral current work says the repeated-packet graph is not a
free tournament.  It is:

```text
K6 minus the affine-zero perfect matching
= octahedron graph
= L(K4).
```

This session found the same graph as an induced shell-transition subgraph of
the `a*b` Farey transform:

```text
shells {3,4,6,9,12,13}.
```

That makes a concrete proof target:

```text
Show that all low-height dangerous Farey perturbations factor through this
octahedral shell carrier, and that its lifted current has controlled divergence
and curl once the known wall packets are included.
```

Equivalently: make the octahedral Hodge problem do real work.  The local
Kirchhoff balance lives on `L(K4)`; the residual curl must be routed to
low-height walls or bounded by the high-height reciprocal/Fourier tail.

## Clebsch, Half-Cube, Paley

These graphs should be used as carriers, not as analogies.

### Octahedron

Role: local current carrier for repeated-residue/fold packets.

The product-Farey signal supports this: `a*b` is the only one of the four
variants that produces an induced octahedron in the LRC shell graph.

### Clebsch

Role: cut-space quotient.

The local S105 computation verified Clebsch as the folded 5-cube, the
`SRG(16,5,0,2)` carrier.  It preserves one-flip signed covariance but destroys
some scalar missed-depth data.  That is exactly the right warning: use Clebsch
for cut labels, not as a scalar proof.

### Half-Cube

Role: covariance/exclusion complement of the Clebsch carrier.

In the new script, the half-cube complement on the folded 5-cube has:

```text
16 vertices, 80 edges, degree 10, 160 triangles, 40 induced octahedra.
```

This makes it a rich reservoir of octahedral subcarriers.  It is a good place
to organize many local tail terms, but it is too dense to be the final metric
object by itself.

### Paley / Fano

Role: sector clock and flat design model.

Paley `T7` has score histogram `{3:7}` and 14 directed triangles, matching the
Fano/BIBD sector story.  It is excellent for the mod-7 sector design, but it
does not know the mod-27 signed shell by itself.  For LRC14 the Paley/Fano piece
must be paired with the `C=27` shell clock.

## The Route Tournament

I ranked eight proof routes by a simple binary-relation tournament on:

```text
preserves LRC predicate,
keeps metric data,
finite exactness,
actionability,
risk.
```

The resulting order was transitive:

```text
1. metric winding tournament: order + gaps
2. covering-set split + equidistribution
3. k=9 Delsarte razor certificate
4. octahedral Hodge tail on L(K4)
5. Clebsch / half-cube cut covariance
6. Paley / Fano sector design
7. four Farey variants as classifiers
8. scalar H / Jensen / additive-energy only
```

This is the right hierarchy.  The Farey variants are not the proof engine.
They are a discovery lens for finite packet carriers, and the product variant
has now pointed at the octahedral carrier.

## Proposed Proof Hook

The next proof-shaped object should be:

```text
Boundary/Farey packet P:
  interval endpoint or local perturbation fraction a/b
  sector label mod 7
  shell label mod 27
  product-height label a*b
  metric gap bracket around the observer
```

Define binary relations on these packets:

```text
R_sector(P,Q): Paley/Fano relation of sector labels mod 7.
R_shell(P,Q): signed shell relation mod 27.
R_octa(P,Q): product-shell transition lies in the octahedral carrier.
R_cut(P,Q): Clebsch/half-cube cut relation after tangent anchoring.
R_metric(P,Q): observer bracket compatibility.
```

Then the proof target becomes:

```text
Any LRC14 counterexample induces a packet complex satisfying all five
relations.  Non-covering packets are already eliminated.  Covering packets
must either:

  (A) fall into the octahedral low-height carrier, where Hodge balance plus
      wall incidence gives a nonnegative reservoir; or

  (B) escape to high height / multi-carrier span, where equidistribution and
      reciprocal-height tails give decorrelation.
```

The missing rigorous lemmas are concrete:

1. Product-Farey carrier lemma: low-height dangerous perturbations project to
   the `a*b` shell-transition graph, and the only induced six-shell obstruction
   is octahedral up to the known symmetries.

2. Octahedral Hodge lemma: the `L(K4)` lifted current has no negative residual
   after the wall packets are added; curl is confined to the eight triangular
   faces and paid by known low-height relations.

3. Clebsch cut-stability lemma: tangent-refined folded-cube classes preserve
   the signed covariance needed by the metric packet, even though scalar depth
   is lost.

4. Paley/Fano sector lemma: the mod-7 sector design supplies the uniform part
   of the Delsarte/decorrelation bound; deviations must enter through mod-27
   shell packets.

5. Metric winding lemma: every packet certificate must include the actual
   observer gap bracket.  This blocks the previously refuted scalar routes.

## Bottom Line

The strongest new lead is:

```text
Farey product a*b -> mod-27 shell transitions -> induced octahedron
{3,4,6,9,12,13}.
```

That is a real finite carrier, aligned with the existing octahedral current
program.  It suggests the right next move is not more raw Farey search, but a
formal packet theorem showing that the low-height part of any LRC14 obstruction
must pass through this octahedral Hodge carrier, while the high-height part is
handled by the already-live decorrelation/equidistribution machinery.
