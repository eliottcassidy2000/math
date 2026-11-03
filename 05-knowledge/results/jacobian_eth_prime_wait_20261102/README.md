# Jacobian classification, ETH scope, and the 290/227 arrival clocks

2026-11-02. Worktree `codex/session-jacobian_eth_prime_wait_20261102`, inherited
from `5e1cedf452`. This packet separates proofs, exact finite checks,
external statements and open classification problems.

The mathematical advance is [W1–W4](WEIGHTED-CLASSIFICATION.md): all target
fibres in the inherited weighted Keller family, arbitrary disjoint-cycle
inertia at infinity, and a partition lower bound for inequivalent maps.
Many remain inequivalent while sharing the exact same positive component.
The number-theory advance is the explicit CRT delay family and the MRB
orientation calculation in [PRIME-MRB.md](PRIME-MRB.md).

## Inheritance and five-concept board

| lane | live concept | inherited mechanism / hostile / next operation |
|---|---|---|
| Anchor | Keller map versus its degree | THM-1330 is explicitly an atlas, not a classification; degree nine already has both atoms and composites. Retain the nonproper-component ledger. |
| Anchor | missing inverse chart | THM-3441's three-point boundary defeats extending a simple-root rule through C=0. Split x=0 from gamma=0 before counting. |
| Niche | deterministic prime arrival | The prime-prefix correction in `duck_primes_20260925` separates retained rank from original cutoff. Retain the zero-hazard prefix and CRT candidate divisors. |
| Wildcard | MRB observation kernel | The contour audit's h_K changes the sum while preserving finite samples, jumps and asymptotics. Test a proposed formula on the whole functional. |
| Niche | ETH and exact counting | The October-27 intake preserved clause width, input size and counting semantics. Re-read present formal statements before transporting an algorithm. |

The least-used relevant sidecar was the hypothesis that **all nonzero roots
of R are simple** in THM-3448. Removing it, while retaining the source
chart, produced the new partition theorem. The canonical corrected near
miss is MISTAKE-236's rejection of a “complete classification” reading of
THM-1330. No chronological log was used as current truth.

## What is actually classified for Jacobian counterexamples

The current canon and primary literature support the following different
levels of statement:

| question | present answer and scope |
|---|---|
| Which dimensions have known complex polynomial constant-Jacobian counterexamples? | Every dimension >=3, by the explicit constructions and stabilization. JC(2) remains open. |
| Which generic degrees occur in each dimension >=3? | Exactly {1,3,4,5,...}; THM-3438 plus the inherited degree-two exclusion. This is geometric/function-field degree, not coordinate degree. |
| Which degrees contain composition atoms? | Every degree >=3; the weighted examples have global monodromy S_n. |
| Which degrees contain composites of two nonunits? | Exactly ab with a,b>=3. Degree nine is the first mixed grade. |
| Are all maps in each degree classified up to polynomial automorphisms? | No. THM-1330 explicitly leaves this open, including uniqueness of the degree-three atom outside tested design families. |
| What new classification is proved here? | W1 classifies every target fibre of the specified weighted family. W2 describes its generic boundary inertia without the simple-nonzero-root restriction. W3 gives many inequivalent maps, not a count of all maps. |

The primary [Gao exposition](https://arxiv.org/html/2608.00222) explains the
tangent-sweep construction and higher-dimensional direction fields. Its
Appendix A already classifies the displayed cubic and quartic fibres;
these were recovered rather than republished as new. The repository's
quartic and cyclic-inertia results are stronger dependencies for the
specific extension pursued here.

**New separation.** At degree five, W3 constructs four inequivalent
determinant-one maps with the same global group S_5. Their exceptional
component signatures are respectively

    (3 finite, identity), (3 finite, transposition),
    (2 finite, transposition), (2 finite, 3-cycle).

Every map also has the ordinary component (3 finite, transposition);
the *unordered multiset* of both components proves inequivalence even
when a coordinate change swaps them. Thus neither generic degree nor
global monodromy classifies the map. The extra information is concrete
boundary behaviour, rather than an unspecified increase in abstraction.

The degree-seven seed R=w³(w-1)(w-2)³ exposes the mechanism especially
clearly: two triple roots produce a 2-cycle and a 3-cycle at infinity.
The difference is their position. The perturbation term -BCw vanishes at
zero but not at 2, so the two Newton equations have different exponents.
The marked zero retains one finite branch; root multiplicity alone would
lose that information. This is an exact instance of an offset cycle law.

## ETH: a current source audit, not a new proof claim

The old [October-27 intake](../external_preprints_mub6_local_smoothing_eth_merge_20261027.md)
treated the fresh manuscript as an unverified external claim. This session
read the current primary formalization interfaces and scope notes. Their
content hashes and retrieval time are in [external_source_hashes.json](external_source_hashes.json).
No local Lean build, comparator run or independent kernel replay was
performed; the local Jacobian/prime proofs do not depend on ETH being false.

There is positive external replay evidence beyond the original source
announcement: the [community replication thread](https://ethdiscussion.blog/2026/10/30/lean-replications/)
and [Carlos Toledo's report](https://carlostoledo.co/eth/) describe successful
kernel checks and distinguish formal theorem acceptance from interpretation
of the machine model. These are cited reports, not executions by this session.

* [OpenAI scope note](https://github.com/openai/math/blob/main/lean/docs/378.md):
  for fixed k, the documented formal result counts canonically encoded
  k-CNF assignments in polynomial space with time
  `(L+1)^O_k(1) 2^O_k(n(log log n)^2/sqrt(log n))`.
  The accompanying manuscript states a sharper exponent with one log-log
  factor, and a sparse-clause improvement. Keep those scopes separate.
* [Anthropic statement file](https://github.com/anthropics/formal-math/blob/main/eth/MainResults.lean):
  one finite-alphabet multi-stack machine and one constant C precede k;
  it states exact counting within
  `2^[C min(km/sqrt(log m), k² n log log n/sqrt(log n))]` for its defined
  duplicate-free encoding, n>=4, m>=2, all variables occurring.
  The [README](https://github.com/anthropics/formal-math/blob/main/eth/README.md)
  supplies build, standard-axiom and independent-comparator instructions.
  Its challenge-file `sorry` placeholders are the specification, not the
  claimed proof; this distinction was checked before interpreting them.

These sources concern subexponential exact counting, not polynomial-time
SAT. Thus they do not establish P=NP. Their asymptotic bounds also say
nothing by themselves about practical crossover sizes.

There is a precise cost check for our tournament program. Substituting an
N=Theta(v²)-variable position encoding into the sharper fixed-width
exponent gives order `v² log log v/sqrt(log v)`, which exceeds v.
Even an N=Theta(v log v) encoding yields order
`v sqrt(log v) log log v`, again exceeding v. Neither bound improves the
existing O(v²2^v) Hamilton-path subset DP merely by invoking ETH failure.
A useful reduction would need, among other things, a much smaller N,
fixed clause width, polynomial encoding length, and exact control of
witness multiplicity. For the displayed sharper rate, a sufficient scale
is `N=o(v sqrt(log v)/log log v)` with N growing polynomially in v.
This is a reduction-size calculation, not an impossibility theorem for
other algorithms or a consequence for the tournament spectrum conjecture.

The compressed tensor-contraction description in the external README is
an interesting research move: a rank tradeoff is iterated while the
auxiliary size changes. For our work the test would be an exact contraction
identity preserving endpoint conditions and source multiplicities. No such
identity is established by an analogy to a tournament grammar.

## Other imported advances: which radical picture is supported?

The [trusted advances ledger](../../reference/TRUSTED-ADVANCES-2026-10-27.md)
already distinguishes named-expert or independently re-proved results
from unreviewed manuscripts. This session does not promote the entire
external collection to proved canon.

* The unit-distance disproof and its [human exposition](https://arxiv.org/abs/2605.20695) use arithmetic
  configurations that defeat a broad asymptotic expectation. They do not
  decide the chromatic number of the plane: many edges and a high chromatic
  obstruction are different predicates. The finite exact unit-distance
  results in the repo remain separate from the asymptotic conjecture.
* The cited zeta progress uses a positive/rank inequality on a finite
  compression; [Lamzouri's independent proof](https://arxiv.org/abs/2609.02882)
  replaces that framework by a Hilbert-space inequality. It does not prove RH, odd-zeta irrationality, or a Collatz
  cancellation bound. Its reusable question is whether an actual positive
  form and a trace estimate can be constructed in the intended source space.
* The local-smoothing and other fresh manuscripts in the ledger remain
  external claims at their recorded tier unless individually verified.

The shared picture is not that finite patterns dictate every global
answer. It is that the **allowed operations and the information they
preserve** matter more than a summary statistic. Here this has four
different, precise meanings: a Keller inverse needs its boundary chart;
a prime first hit needs its starting clock; an MRB filter needs its sign
sector; and a SAT reduction needs the actual encoding and multiplicity.

## Connection ledger: maps, preserved predicates, and failures

| source -> target | map / what is preserved | lost information and decisive test |
|---|---|---|
| R's root multiplicities -> Keller infinity cycles | W2's Newton expansions; generic sheet count and inertia | positions of roots and all higher gluing data; W3 gives a lower bound, not completeness |
| integer partitions -> inequivalent polynomial maps | formula (6) in W3, followed by the inherited lift | no inverse classification; compare the unordered component ledger |
| affine prime candidates -> CRT progression | preserve specified proper divisors | primality of later slots; 5m+1 need not be prime |
| 227's Stern obstruction -> N=227 mod165 | preserve the first ten quadratic candidate failures | growing allowed a range; N=557, a=12 is the hostile |
| units -> sign quotient | U(290)/±1 has 56 classes and preserves the ±1 target event | orientation characters; the signed MRB sum converges while the merged filter diverges |
| MRB profile -> jump plus finite data | inherited h_K preserves those observations | the alternating functional changes explicitly |
| SAT instance -> counted tournament object | only a proposed encoding at present | clause width, variable inflation, auxiliary witness multiplicity; no automatic runtime gain |

## Incoming work integrated after the first checkpoint

The clean-tree fetch brought in `b687cd28ca`, `7b76044ca1` and
`aaf2c12555`. The first two are the independent
[Session144 synthesis](../session144_jacobian_eth_mrb_20261103/SYNTHESIS.md).
They agree on the distinction between degree spectra and full map classes,
and add stronger results which this packet inherits rather than duplicates:

* The corrected equivariant form allows bivariate coefficients. The
  primitive nonlinear-scale question cannot be phrased as raw degree:
  an equivariant shear makes that degree arbitrarily large within one
  equivalence class. Our W3 uses intrinsic target nonproper-component
  data under unrestricted equivalence instead of raw displayed degree.
* Its universal positive-coordinate theorem and disjoint iteration levels
  concern source orbit-zero divisors in an equivariant category. Our
  target nonproper components are different objects. W4 supplies a strong
  comparison: p_part(n-3) inequivalent primitive maps can share the
  **identical** positive coordinate x(1-xy/2+x²z). Complete information
  about that single coordinate still cannot classify their joint map.
* The incoming MRB theorem strengthens h_K to smooth-Laplace perturbations
  preserving all algebraic-order asymptotics uniformly on the entire
  right half-plane. Its complete periodic-mean criterion includes our
  even-modulus specialization. It also proves a certified precision
  algorithm and distinguishes cutoff from Abel finite parts; none
  recognizes a modular closed form. Our contribution on this lane is
  the explicit 290/227 CRT and sign-sector synthesis, not a claim to
  improve that stronger analytic obstruction.
* Its bounded-integer Keller/SAT interface counts declared finite witnesses
  and includes verification gates. The W1 family admits a different,
  unconditional exact reduction: **two univariate gcds** count every
  complex fibre at a rational target. This preserves the full target
  predicate within the specified family without Boolean encoding.

The latest [native module-exchange packet](../uniform_prime_transport_20261103/README.md)
separates fixed-attachment means from legality after native re-marking;
it does not claim spectrum completeness. This complements the earlier
[endpoint-block recurrence](../primegame_port_spectrum_20261103/FOREST_RESPONSE_RECURRENCE.md):
compression is faithful for its specified operation class, and a new
guard can demand information outside it. The MRB orientation loss and
the missing Keller chart provide concrete mathematical analogues, not
unproved reductions between these problems.

## Audit, searches and handoff

Repository searches used exact theorem slugs, fibre/degree/monodromy
synonyms, their correction lineage, and numerical 290/227 occurrences.
Confirmed older hits are linked in PRIME-MRB. Hash substrings, theorem
numbers and session numbers were excluded as mathematical connections.
The search was broad but is not a claim to have read every historical
occurrence. Recovered MRB packets from October 22 and 24 prevented both
re-proving their weaker statements and reviving the disproved q^4 formula.

The checker outputs record exact universes and hostile controls. The
general arguments appear in the proof files. All checks use explicit
exceptions rather than Python assertions, so optimized execution preserves
them. Normal and optimized outputs must match before close-out.

Next decisive work:

1. At fixed weighted degree, test whether root-position cross-ratios yield
   a positive-dimensional family of inequivalent maps after allowing all
   polynomial source/target automorphisms. The component ledger alone will
   usually be too coarse; a branch-curve invariant may survive.
2. Seek a global-cover obstruction beyond S_n and determine which proposed
   covering data admit polynomial gluing. This is the actual gap behind a
   complete classification, not the already solved degree-value spectrum.
3. For prime arrivals, compare explicitly specified random *starting-window*
   models with the fixed-origin clock; retain the local CRT sieve state.
4. For MRB, evaluate the oriented residue functional without replacing it
   by a density multiple. Any proposed modular identity must overcome the
   existing h_K obstruction using additional global information.

The primary requested universal classifications remain open where stated.
The partition theorem and the exact CRT/MRB identities are completed,
reviewable steps rather than declarations that those larger problems closed.
