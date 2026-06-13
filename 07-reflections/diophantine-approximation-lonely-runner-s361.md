---
source: codex-2026-05-30-S361
status: exploratory synthesis after broad repo search
tags:
  - diophantine-approximation
  - lonely-runner
  - bohr-sets
  - quotient-gaps
  - endpoint-protection
  - residue-calculus
---

# Diophantine Approximation Threads Toward Lonely Runner

## Search Pass

I searched broadly for literal and adjacent Diophantine approximation language:

```text
diophantine, approximation, Kronecker, Dirichlet, Bohr, Weyl,
continued fractions, Farey, Minkowski, geometry of numbers,
lattice, torus, mod 1, nearest integer, fractional part,
covering radius, quotient gaps, endpoints.
```

The useful repo hits fell into five strata.

## Stratum 1: Direct Lonely Runner Diophantine Approximation

The strongest direct spine is:

```text
HYP-1794  Lonely Runner as a quotient-gap certificate
HYP-1802  endpoint-protection obstruction
HYP-1810  unit-boundary skeleton
HYP-1811  no all-protected endpoint core
HYP-1812  endpoint kernel pressure
THM-357   endpoint-protection trichotomy
```

The common object is the open forbidden union

```text
F(V) = union_{v in V} {t in R/Z : ||v t|| < 1/(k+1)}.
```

This is a finite inhomogeneous Bohr-cover problem.  The LRC witness is not an
ordinary good rational approximation; it is an **anti-approximation** time:

```text
||v t|| >= 1/(k+1) for every v.
```

For integer speeds, every endpoint lives in

```text
Q(V) = (k+1) * lcm(V).
```

THM-357 says the continuous circle contributes no hidden fourth case:

```text
positive open cell        -> witness interval
full measure + endpoint   -> boundary witness
full open cover           -> all endpoints protected
```

So the Diophantine approximation problem becomes finite endpoint incidence.

## Stratum 2: Distance-Graph Colorings

`07-reflections/lonely-runner-distance-graph-colorings-s359.md` translates LRC
into regular circular coloring of the integer distance graph:

```text
c_t(x) = t x mod 1.
```

The usual Diophantine approximation map `t -> (v_1 t,...,v_k t)` becomes a
restricted family of homomorphic circle colorings.  This matters because an
arbitrary circular coloring is too flexible; LRC asks for a single multiplier.

The multiplier language makes the exact obstruction:

```text
bad multipliers cover the multiplier circle
boundary multipliers are rational endpoints
endpoint protection is finite arithmetic incidence
```

This is a cleaner bridge to graph theory than the vague phrase
"simultaneous approximation."

## Stratum 3: Quotient Covers Outside LRC

The base-42 Erdos-Straus thread is the closest non-LRC analogue.

In `07-reflections/base42-erdos-straus.md`, identity families cover residue
classes, hard classes are quotient gaps, and multi-`r` fallback repairs
boundary classes.  The grammar is:

```text
choose finite quotient
cover residue fibers by parametric identities
study the classes that survive the first cover
```

LRC has the same grammar, but the fibers are geometric endpoint cells on the
time circle rather than congruence classes of primes:

```text
Erdos-Straus: residue class survives identity cover
LRC: endpoint survives forbidden-arc cover
```

This is why the repo's residue calculus is relevant.  It is not replacing
Diophantine approximation; it is discretizing the surviving boundary.

## Stratum 4: Lattice/Torus/Rapidity Echoes

There are many lattice and torus hits that are not directly LRC but are useful
as analogies:

- THM-252: rapidity values lie in a logarithmic prime lattice, with Baker's
  theorem giving independence.
- `lattice_tricks_s91i.out`: denominator-killing primes and roots of unity
  show how a modulus can trivialize rational coordinates.
- `torus_helix_mobius_s90bj.out`: integers as points in a prime-exponent
  lattice, with squarefree points as a boundary skeleton.
- continued-fraction/golden-ratio threads: useful as "bad approximation"
  metaphors, but mostly not operational for LRC yet.

The most transferable lesson is the denominator-killing one.  In LRC, the
denominator

```text
Q(V) = (k+1)lcm(V)
```

is not just a nuisance.  It is the finite phase space where endpoint exposure
is exact.  The right primes are not "small denominators" in the classical
continued-fraction sense; they are primes that make quotient fibers visible or
force divisibility.

## Stratum 5: False Friends

Some hits are important to filter out.

`04-computation/tournament_dirichlet_s291b.py` is about a Dirichlet series,
not Dirichlet approximation.  It is still useful as an analogy for dominance
by smallest terms, but it should not be cited as a Diophantine approximation
anchor.

Many continued-fraction hits concern rational tournament ratios, golden
shadows, or generating functions.  These are number-theoretic, but not
currently evidence for LRC.

The strong LRC route is not:

```text
find a good rational approximant by Dirichlet/Kronecker.
```

The strong route is:

```text
prove a finite anti-Bohr set cannot be all-protected at its boundary.
```

## External Constraint

The literature check supports this distinction.

Horvat and Stoffregen prove an almost-all velocity result and explicitly warn
that a traditional Diophantine approximation approach cannot settle the
problem.  Tao's 2017 paper improves the general lower bound and gives a finite
checking reduction for bounded integer speeds.  Chow's Bohr-set work shows
that inhomogeneous Bohr sets can contain structured generalized arithmetic
progressions, extending a Tao construction.  Jensen's mixed-threshold work
adds Fourier formulas for safe/unsafe indicators.

Taken together, these suggest:

```text
Metric/density DA explains why generic speeds are easy.
Finite endpoint incidence explains why exceptional speeds are hard.
Bohr/GAP/Fourier kernels may bridge the two.
```

References checked:

- Horvat-Stoffregen, `arXiv:1103.1662`.
- Tao, `arXiv:1701.02048`.
- Chow, `arXiv:1703.07016`.
- Jensen, `arXiv:2605.27941`.

## Creative Connection To Lonely Runner

The repo's current LRC program can be reframed as:

```text
Lonely Runner = finite anti-Bohr boundary theorem.
```

For a speed set `V`, each forbidden arc is a Bohr neighborhood of an integer:

```text
B_v = {t : ||v t|| < 1/(k+1)}.
```

The safe set is the complement of the union of these Bohr neighborhoods.
Classical simultaneous approximation wants to put `tV` near an integer vector.
LRC wants the opposite: keep every coordinate away from integers by the
critical distance.

The Diophantine difficulty is not that the orbit is hard to describe.  For
integer speeds, the boundary is rational and finite.  The difficulty is that
the bad Bohr neighborhoods can overlap so efficiently that the open cells
disappear.  THM-357 then says the only possible counterexample has all
endpoints protected.

That gives a three-layer attack:

1. **Unit skeleton.**  HYP-1810 says known tight examples leave exactly the
   nonzero unit residues `(Z/(k+1)Z)^* /(k+1)` unprotected as forbidden
   endpoints.  This is the first residue quotient.  The point `0` is part of
   the Dirichlet-pigeonhole orbit but is a forbidden center, not an endpoint.

2. **Boundary peeling.**  HYP-1811 asks whether every full-measure endpoint
   graph has an unprotected endpoint after peeling.  This is the incidence
   proof layer.

3. **Kernel pressure.**  HYP-1812 asks whether a nonnegative Fejer/Riesz
   kernel on `Z/QZ` must see positive closed safe mass.  This is the analytic
   proof layer.

The new Diophantine insight is that these are not separate routes.  They are
the same route at different resolutions:

```text
unit quotient        = first denominator
endpoint graph       = exact finite inhomogeneous Bohr boundary
Fejer/Riesz kernel   = Fourier certificate of surviving safe mass
GAP/Bohr structure   = possible descent mechanism inside the boundary graph
```

## New Hypothesis

This search suggests HYP-1813:

```text
Every hypothetical all-protected LRC endpoint core admits a smaller
Bohr-boundary quotient, obtained either by dividing a forced common
divisibility class or by extracting a structured progression inside the
endpoint set.  Iterating this descent produces an unprotected endpoint.
```

In plainer language: if a counterexample protects the obvious unit residues,
the protection edges should create a new quotient where some endpoint is less
protected, not more.  The obstruction cannot stay all-protected through every
Bohr/GAP refinement.

## Concrete Next Experiments

1. Add endpoint-core peeling to the S359/S360 scripts.
2. For each peeled layer, record the quotient denominator and the residue
   subgroup generated by unprotected endpoints.
3. Search for generalized arithmetic progression structure inside protected
   endpoint sets, inspired by Bohr-set/GAP results.
4. Build the closed safe product on `Z/QZ` and test shifted Fejer/Riesz
   kernels from HYP-1812.
5. Compare near-tight positive-gap examples against the unit-boundary skeleton:
   do they fail at the unit quotient, or only at higher endpoint denominators?

The key psychological shift: stop asking whether Dirichlet approximation can
find a time.  Ask why the finite anti-Dirichlet boundary cannot be completely
protected.

## S362 Formalization Addendum

The next session turned two parts of this note into formal objects.

THM-358 proves the initial-segment case:

```text
V={1,...,n-1} has safe set {a/n : gcd(a,n)=1}.
```

This is exactly the equality case of Dirichlet's pigeonhole approximation
argument.  The unit-boundary skeleton is therefore not just a pattern in the
data; for the standard tight family it is the rigid equality case of the basic
Diophantine lemma.

THM-359 proves that endpoint/interval peeling computes the largest finite
protection core.  The S362 computation implements this and finds empty cores in
all inherited full-measure primitive-box examples.  In tight cases, the first
peeled quotient layer is `unit_mod_n`; in near-tight positive-gap cases, the
first layer lives at higher endpoint denominators.

THM-360 adds the first divisibility filter inside the unit layer: a unit
endpoint `a/n` can only be strictly protected by a speed divisible by `n`.
Thus any full-open-cover counterexample must contain a speed divisible by
`k+1` before higher Bohr-boundary structure can even enter.
