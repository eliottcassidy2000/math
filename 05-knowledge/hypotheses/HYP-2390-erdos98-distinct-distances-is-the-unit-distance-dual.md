---
id: HYP-2390
name: erdos98-distinct-distances-is-the-unit-distance-dual
status: SURVEY + bridge (mechanism & CM tension VERIFIED) + a creative attack direction
date: 2026-06-08
session: claudebox-2026-06-08-S722
external:
  - "Erdos problem 98 (erdosproblems.com/98, OPEN, no prize, tags geometry/distances)"
  - "Pach: h(n) < n^(log2 3); Erdos-Furedi-Pach: h(n) < n exp(c sqrt(log n))"
depends_on:
  - THM-440   # u(22) unit-cocyclic extension: max extension degree 3 (the radius-1 slice of 98)
  - HYP-2376  # unit-distance momentum twistor (angle = log U(1)); CM rotation orbits (S719)
  - HYP-2350  # unit-distance product anatomy (S712)
  - THM-441   # convolution-correlation-adjoint duality / autocorrelation operator (S714)
provisional_id: true
---

# HYP-2390: Erdos problem 98 is the dual of our unit-distance work; the obstruction is the shared "3"

## The problem (OPEN)

`h(n)` = the largest number such that any `n` points in `R^2` with **no 3 collinear and no 4 concyclic**
determine at least `h(n)` distinct distances. **Does `h(n)/n -> infinity`?** Erdos could not prove even
`h(n) >= n`. Known UPPER bounds (constructions with few distinct distances): Pach `h(n) < n^(log2 3)`;
Erdos-Furedi-Pach `h(n) < n exp(c sqrt(log n))`.

## The bridge to the repo (the dual)

Problem 98 and our unit-distance work (S711-S719) are the two opposite extremizations of the SAME radial
autocorrelation (Patterson) of a point set:
- **unit distances** MAXIMIZE one distance's multiplicity = the autocorrelation PEAK `u(n)`;
- **problem 98** MINIMIZES the number of distinct distances = the autocorrelation SUPPORT, under a
  general-position constraint.

**THE SHARED "3" (mechanism, VERIFIED).** No 4 concyclic => no config point has 4 others at a single
distance (4 equidistant points lie on a circle centered at that point => 4 concyclic) => **every
distance-graph has max degree <= 3** => each distance multiplicity `m(r) <= 3n/2` =>
```
   h(n) >= ceil((n-1)/3).
```
This "3" is *exactly* THM-440/S719's "max unit-cocyclic extension degree = 3" that caps `u(22)` at 60:
the local obstruction "you cannot put 4 of these points on a unit circle" becomes, globally over all
radii, "no 4 concyclic." **Problem 98 is the all-radii form of the u(22) ceiling.** (VERIFIED: random
no-4-concyclic 7-point configs all have distance-graph max degree <= 3.)

## The CM tension (VERIFIED) -- why 98 is hard for our methods

Our few-distinct-distance optimum is a CM lattice (S719: few norms = peaked support). But a CM lattice
point has `r_Q(D)` others on each norm-`D` circle: `r_Q(1)=r_Q(3)=6`, `r_Q(7)=r_Q(13)=12`. Every one is
`>= 4` => **massively 4-concyclic => FORBIDDEN by problem 98.** So problem 98 forbids exactly the CM
structure that minimizes distinct distances; its extremal configs must be "CM-broken" (the
Erdos-Furedi-Pach perturbation keeps few distinct distances while destroying all 4-concyclicity). The
unit-distance OPTIMUM (peaked, concyclic, CM) and the problem-98 OPTIMUM (spread, general position) sit at
opposite poles of one autocorrelation.

## Three reframes (exploring threads)

1. **Coloring (repo's coloring-unification):** minimizing distinct distances = decomposing `K_n` into the
   FEWEST distance-classes each `<= 3`-regular AND simultaneously realizable in the plane with no 3
   collinear. `D >= (n-1)/3` is the `3`-regular packing bound; the planar realizability is the hard part.
2. **Sidon / additive energy (our S718 twistor / half-system thread):** "no 4 concyclic" is a Sidon-type
   (`B_2`) ban on coincidences. A Sidon set has ALL distances distinct (`D = C(n,2)`, max support); few
   distinct distances needs many coincidences; no-4-concyclic caps the coincidences. The gap between
   `(n-1)/3` and Erdos's hoped `>= n` is the slack between the WEAK consequence (per-point multiplicity
   `<=3`) and the FULL ban (no 4 on ANY circle, including non-config centers) -- which the trivial bound
   does NOT use.
3. **Temperature (S720/S721):** distinct distances = support of the radial measure; concyclic/CM = the
   COLD crystalline (peaked, minimal-support) configs. No-4-concyclic forbids the coldest configs.
   `h(n)/n -> infinity` asks: **does forbidding the cold structure force the support to WARM up
   (superlinear)?**

## The creative attack direction

The trivial bound `h(n) >= (n-1)/3` uses only the per-point multiplicity `<= 3`. The UNUSED strength of
"no 4 concyclic" bounds the number of concyclic 4-tuples to ZERO -- a global ban on circle-incidences,
i.e. on isosceles/perpendicular-bisector coincidences. **Conjecture/route:** bound the additive energy of
the distance structure (the count of isosceles 4-tuples) using the zero-concyclic ban, then convert energy
-> support via Cauchy-Schwarz (Elekes-style) to push `h(n)` past `(n-1)/3` toward superlinear. The repo's
autocorrelation operator `MM*` (THM-441) is the right bookkeeping: distinct distances = `rank`/support of
`MM*`'s radial spectrum; the concyclic ban = a bound on the off-diagonal "circle" incidences of `MM*`.

## Next
- improve constructions (the Erdos-Furedi-Pach `n exp(c sqrt log n)`; a CM-broken / perturbed-lattice
  family) and compute the small-`n` constrained minimum distinct-distance table;
- make the additive-energy route precise: does zero-concyclic => `o(n^2)` isosceles 4-tuples => `D = omega(n)`?
- the radius-1 slice is THM-440/u(22): use the u(22) census as a microscope on the n=22 case of 98.
