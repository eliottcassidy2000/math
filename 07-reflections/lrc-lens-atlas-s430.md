---
source: codex-2026-05-31-S430
status: exploratory synthesis
tags:
  - lonely-runner
  - tournament-recursion
  - incidence-hypergraph
  - view-obstruction
  - zonotopes
  - distance-graphs
---

# LRC Lens Atlas

The question was: how many different lenses can LRC be seen through, and which
of them are really isomorphic problems rather than just analogies?

The short answer from S430 is:

```text
at least 20 useful lenses,
11 in one formulation/transport component,
10 connected as proof technologies around the same obstruction core.
```

Some overlap because a lens can be both a formulation and a proof machine.  The
central object is not "runners."  It is a labelled incidence problem.

```text
rows       = constraints/endpoints/cells/cuts that must be protected
columns    = speeds/protectors/arcs/gates
witness    = uncovered row or positive gap
bad object = full cover with every boundary protected
```

## External Map

The web literature confirms that LRC is unusually polyglot.

- Perarnau and Serra's survey presents it as a problem with many partial
  results, tools, and applications.
- The modern finite-checking papers formulate the runner statement as
  `||v_i t|| >= 1/(n+1)`, reduce to integer velocities, and then use zonotope
  geometry and lattice-point bounds.
- Henze and Malikiosis identify a covering-radius / lattice-zonotope
  reformulation connected to view-obstruction and billiard motion.
- Bienia, Goddyn, Gvozdjak, Sebo, and Tarsi connect the problem to
  nowhere-zero flows and view obstructions.
- Giri and Kravitz define Lonely Runner spectra by distances from proper
  subtori of `(R/Z)^n` to the center point.
- Recent computer-assisted work is moving fast: Rosenfeld proved eight
  runners; Trakulthongchai proved nine and ten; Sungkawichai and
  Trakulthongchai's April 2026 preprint proves the `k=10,11,12` nonzero-speed
  cases, i.e. up to 13 total runners in the `k+1` convention.

That last point matters for repo language.  Our "n=14" frontier is still beyond
the currently reported `k=12` / 13-total-runner result, but it is now adjacent
to the public frontier rather than far beyond it.

## Exact Lenses

The formulation/transport component in `lrc_lens_atlas_s430.py` is:

```text
runner frame
Diophantine norm
forbidden interval cover
endpoint-protection hypergraph
coarse denominator sieve
torus line / proper subtorus
view obstruction
distance graph / circular coloring
nowhere-zero flow implication family
lattice zonotope covering radius
Lonely Runner spectrum
```

I would call these "close enough to transport proof ideas," but not all in the
same way.  Some are strict restatements.  Others, like the flow theorem and
coarse sieve, are implication or necessary-constraint lenses: proving the LRC
instance supplies the external combinatorial statement, or a counterexample
must satisfy the extracted row condition.

The repo's THM-357 is the key local normal form.  Once LRC becomes:

```text
forbidden endpoints must all be strictly protected
```

most of the other formulations are ways of choosing a basis for the same rows.

## Proof Lenses

The second component is not exact in the same sense, but it is where proofs
live:

```text
finite checking
prime-product sieve
integer program / set cover
Bruhat-Tits frontier
adelic gap-debt product
natural x+2 / x*2 modes
product-sum and Egyptian ledgers
mixed thresholds
shifted LRC and coloopless zonotopes
```

The variant lenses are especially useful as warning systems.  The shifted LRC
now has counterexamples in broader zonotope classes, which is the LRC analogue
of a recurring tournament lesson in this repo: a nearby operation can look
identical until a complement, quotient, or shift changes the category.

## Tournament Translation

The tournament bridge is not that LRC "is a tournament."  It is that both are
complete labelled incidence systems.

```text
LRC forbidden endpoint       <-> Hamiltonian-path cut
LRC protecting interval      <-> backward arc crossing the cut
unprotected endpoint         <-> bad cut / non-SC witness
all endpoints protected      <-> every cut protected inside an SCC
p-adic endpoint depth        <-> projection defect after recursive motion
finite speed enumeration     <-> tournament enumeration modulo S_n/complement
distance/circular coloring   <-> Omega(T) chromatic/fugacity axis
```

This is why the natural-number modes keep reappearing.

```text
addition shadow:       x -> z iff x<z      transitive tournament
multiplication shadow: x -> z iff x|z      divisibility DAG
```

Tournament recursion has:

```text
n -> n+2       add-two / top-row motion
n -> 2n        blowup / row motion
```

LRC denominator recursion has:

```text
odd root b -> b+2       column payload changes
n -> 2n                 row endpoint debt transfer
```

So the strongest new statement from S430 is HYP-1900: LRC and tournament
recursion share a labelled incidence core.  The proof technologies differ, but
the failure mode is the same:

```text
no uncovered row,
no positive gap,
no private leaf,
no nonzero projected defect.
```

## New Questions

1. Can the recent prime-product sieve be rewritten as a dual certificate in
   the S420 row/column IP matrix?
2. Does the Giri-Kravitz spectrum recursion `S(n)` accumulating on `S(n-1)`
   correspond to tournament vertex deletion, cut quotient, or root-spectrum
   projection?
3. Can the view-obstruction cube arrangement be turned into a tournament
   tiling polytope whose facets are good cuts and odd-cycle conflict rows?
4. Can disproof searches be inverted: generate leafless labelled protection
   hypergraphs first, then ask whether any primitive speed set realizes them?
5. Is the row/column natural-number split the common skeleton behind zonotope
   paving, p-adic endpoint depth, finite-checking sieves, and product-sum
   factor packing?

## Web Sources

- Perarnau and Serra, `The Lonely Runner Conjecture turns 60`:
  https://arxiv.org/abs/2409.20160
- Malikiosis, Santos, and Schymura, finite checking and zonotopes:
  https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/linearly-exponential-checking-is-enough-for-the-lonely-runner-conjecture-and-some-of-its-variants/A51A991DE89B8C9C2E2FF13FBD4501DA
- Henze and Malikiosis, covering radii / view obstruction / zonotopes:
  https://arxiv.org/abs/1609.01939
- Bienia, Goddyn, Gvozdjak, Sebo, and Tarsi, flows and view obstructions:
  https://pagesperso.g-scop.grenoble-inp.fr/~seboa/sebo_files/papers/jctb98runner.pdf
- Giri and Kravitz, Lonely Runner spectra:
  https://arxiv.org/abs/2304.01462
- Rosenfeld, eight runners:
  https://arxiv.org/abs/2509.14111
- Trakulthongchai, nine and ten runners:
  https://arxiv.org/abs/2511.22427
- Sungkawichai and Trakulthongchai, eleven, twelve, and thirteen runners:
  https://arxiv.org/abs/2604.23906
- Jensen, mixed thresholds:
  https://arxiv.org/abs/2605.27941
- Blanco, Criado, and Santos, shifted LRC counterexamples / coloopless
  zonotopes:
  https://arxiv.org/abs/2603.24784
