# Tournament boundaries seen from every prime

2026-11-03. **PROVED scoped / FINITE-EXACT / OPEN**, separated below.
Same-agent proof review with independent exact algorithms. No global
Hamilton-spectrum, Collatz, planar-Jacobian or plane-colouring claim.

This continues `digit_memory_all_bases_20261103` from clean origin/main
`d53992645e` in `/tmp/math-wt-digit_memory_all_bases_20261103`.
The clean integration at `a55a46e2a9` adds the S174 norm bridge and the
return18979 carrier obstruction; see the final section of CONNECTIONS.

The main result is an all-prime extension of the two-boundary-incidences
principle. With M old uniform modules of fixed path-cover contents c_i,
adding q vertices without changing old edges forces the divisor

    gcd_{|U|=max(M-2q,0)} product_{i in U} c_i.

At each prime it retains the smallest max(M-2q,0) valuations. It is also
a determinantal divisor / exterior-power content, and has a Fitting-ideal
version over non-UFD coefficient rings. This is an exact bridge from
literal Hamilton paths to arithmetic and geometry, not an analogy based
only on repeated integers. See AP1--AP2 in [PROOFS](PROOFS.md).

Further results and their boundaries:

* **PROVED, with a complete finite minimality component:** two explicitly
  labelled tournaments of orders4 and5 have contents1,1 but their ordinal
  sum has content5. Their slot polynomials cover complementary residue
  pairs modulo5. Total order9 is minimal for a new prime to arise this way.
  Individual minima lose the shared argument; the full p-adic valuation
  functions preserve it. AP3.
* **PROVED exact counting rule:** deleting two named exceptional vertices
  gives a three-segment, seven-choice-per-module compiler. It preserves
  actual modules and roots. A28-vertex heterogeneous example has
  H=1,133,904,592,664,865. AP4; its displayed count is FINITE-EXACT.
* **PROVED elementary synthesis:** coherent residues at every prime do
  not by themselves select an ordinary integer. A uniform bound on their
  balanced representatives does. An explicit non-diagonal idempotent and
  a unit-valued product-formula hostile locate the missing information.
  AP5--AP6 connect this to every-base digit tricks and Collatz clocks.
* **PROVED corollary of the inherited cut quadratic:** the gcd of every
  one-birth response is exactly the gcd of its integral quadratic
  coefficients. **FINITE-EXACT hostile:** the9-vertex H55 source has this
  gcd1 yet none of its512 one-birth cuts reaches111,167 or169. AP7. A
  primitive ideal is not an actual compatible continuation.

All numbers here count raw directed Hamilton paths, not H/Aut, isomorphism
classes, undirected paths or presentation counts.

## Inheritance and concept board

Closest proved mechanism: TM2--TM3 in
[the previous packet](../digit_memory_all_bases_20261103/TOURNAMENT.md)
give sharp boundary budgets for common-content3 modules and hidden-bit
interaction. Fixed slot divisors and wrapper valuations already have the
stronger complete treatment in
[G3--G5](../coverage_machine_induction_20261103/gluing/NOTE.md).

Canonical hostile: equal full profiles of the two marked5-vertex seeds
do not determine nonuniform reader responses. Corrected near miss: the
H15 neutral update increases the uniform divisor3->15 while decreasing
the named continuation divisor3->1; these are distinct observables
([ordered response](../ordered_response_ten_20261103/INCOMING.md)).
Least-used sidecar: the index of a prime ideal above a rational prime,
as exposed by the norm19 direction in
[cubic adelic geometry](../cubic_adelic_20261016e.md).

|Lane|Live object|Retained information and decisive test|
|---|---|---|
|Anchor|Actual-source tournament returns|Old edges, roots and full cut response; primitive ideal must not be mistaken for a return|
|Anchor|Module cover contents|All prime valuations at once; AP1 against arbitrary one/two-birth cuts|
|Niche|Shared-slot p-adic functions|Argument x, not only minima; order4+5 complementary zero sets|
|Niche|Exterior powers and ideals|Actual minors/ideal generators; six pair coordinates for four modules|
|Wildcard|All-prime residue trees|Common source and real height; non-diagonal CRT idempotent|
|Wildcard|Square-sum and MRB probe banks|Use existing quadratic reconstruction, retaining names and future context|

The board changed after reading concurrent work: the new exact response
rank theorem concerns degree<=2q Walsh functions, whereas our exterior
map concerns arithmetic determinants. The common boundary incidence is
real; the two exterior-looking representations are not the same algebra.
See [CONNECTIONS](CONNECTIONS.md) for maps, losses and next tests.

## Reproduction and next target

From repository root:

```sh
python3 -B 05-knowledge/results/all_prime_tournament_geometry_20261103/reproduce.py
```

Uses the Python standard library and the previous packet's literal graph
DP. Runs normal and optimized modes and checks byte-identical receipts.
See [AUDIT](AUDIT.md) for the complete finite universes and independent
paths. No subagent was used.

The next decisive experiment is a **marked prime-phase return atlas** with
two feedback ports. Incoming `return18979_ports_20261103` closes the uniform
one-apex +2 route and identifies the constant-mode debt1182. Retain the
actual18979 root43/spare6 and the separately supplied18977 source; compute
ideals and feasible values separately. Seek the exact target coefficients
(3729,-3181), rather than accepting the one-apex substitute(3729,-1999).
Extending the new compiler to arbitrary two readers requires joint endpoint
profiles; replacing them by marginals is not licensed. Full ternary source
closure remains OPEN.
