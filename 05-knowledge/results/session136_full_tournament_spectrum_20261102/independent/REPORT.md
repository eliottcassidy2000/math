# Session136 independent hostile audit

**Verdict: PASS for the stated framework theorems and finite certificates.**
The unrestricted conjecture that all positive odd H except 7 and 21 occur
remains **OPEN**. No inherited finite bank or exclusion is presented here as
a newly proved or independently regenerated result.

## Accepted results and boundaries

1. **Finite-bank witness compression.** The attained-only bound
   `mu(h)<=ceil(22 log_80267 h)` follows from the stated batching argument,
   inherited positive banks, exact strong floor, and six strict residual
   inequalities. All new integer inequalities pass independent checks.
   It improves the ordinary deletion-core coefficient, not the smaller
   strong-only coefficient. Replacement preserves H, not induced inclusion,
   strongness, or continuation states. The min-plus recursion and strict
   cost-atom framework are sound. The separate cofinality note is an
   elementary multiplicative consequence, not integer completeness; the
   explicit 65/67 threshold is exactly 138.

2. **Context tomography and capped substitution algebra.** The exact
   homogeneous-context detector depth is `2c-2` for a first profile difference
   at c. The active-module barrier `H>=2^n+1` and the tagged finite capped
   state are sound. The tag and small profiles cannot be erased. Arbitrary
   substitution composes through those states; this justifies full-depth
   fixed points, not merely bounded-depth experiments. Independent fixed
   points reproduce 58 and 65 states for the two J/C3 seed grammars and
   their complete supports through 127. Allowing all quotients of order at
   most five gives 86 states and 49 values through 127. The exact missing
   list is independently reproduced. In particular the verified width
   results are `w(13)=w(19)=5` and `w(31)=6`, not an unbounded width theorem.
   The width-five strong-output cutoff is safe: large active inputs
   saturate by the module barrier; all-small strong outputs of total order
   at least eleven have H at least the inherited `f(11)=225>127`.
   Nonprime five-vertex quotients reduce through a proper module to
   quotients/factors of order at most four.

3. **Obstruction lifetimes and basis growth.** The interval `[d(F),H(F))`,
   endpoint conservation, and SCC maximum-ratio rule have the asserted
   exact directions. Empty/unit cases are separated. Pigeonholing fixed-cycle
   strong tournaments proves the superpolynomial subsequence lower bound
   for explicit full bases only; it says nothing about compressed
   descriptions or recognition complexity. The uniform two-atom lower
   bound is valid, with both H15,d5 atoms checked literally and distinguished
   by score multisets. Independent enumeration recovers all strong classes
   through order seven and the 104 strong seven-vertex births at d=37.
   Scalar paths in the lifetime DAG are not compatible embedding chains.

4. **Positive even-path cut geometry.** The OCF-derived R kernel, strict
   uncrossing, terminal chains, compatible low splits and forced-hole bounds
   have correct factors and punctured-domain conditions. The circulation
   kernel is exactly the information lost by the cut function. The n6
   code83 balanced minimum is verified literally; no smaller strong parent
   fails the singleton/complement minimum test. The vertex-transitive block
   argument and primitive specialization are sound: block sizes form a
   divisibility chain, and oddness gives the logarithm base three. Circulant
   minima reduce to proper-subgroup representatives. The F1 normalization
   in the primitive mean/minimum formula matches inherited THM-4111.
   These are one-parent holes, not global spectral holes. Independent cyclic
   controls recover 52 missing bottom-band levels at n9 and 1,062 at n11.

5. **Strong one-apex local universality.** See `LOCAL_CLOSURE.md` for the
   full proof audit. The positive-word realization, odd commutator,
   synchronized CRT reset, polynomial numeric word bound and unbounded
   residue fibers are valid. The result is profinite odd-density, not
   integer completeness, natural density, or a short-witness theorem.

6. **Inherited-lens hostile.** The root GF(2) example is valid: the
   transitive six-tournament has H1 and code280 has H41, while both adjacency
   matrices have rank powers `(6,5,4,3,2,1,0)` over GF(2), hence the same
   single nilpotent Jordan block. This corroborates loss of H under that
   lens, not a new classification theorem.

## Independent executed verification

No independent script imports producer code or a producer data catalogue.
`kernel.py` uses literal vertex permutations, backwards endpoint recursion,
literal permutation-and-cut covers, and literal odd-cycle collections. The
three H engines agree on every labelled tournament through order five.

| script | independent universe/consequence |
|---|---|
| `audit.py` | 1,100 labelled tournaments at orders 0..5; 1,098 literal ear extensions; 2,850 literal one-hole context cases; gluing depth, interaction and module bounds; H23 collision; compression base/step arithmetic |
| `local_check.py` | 255 actual binary-word graphs; 63 strong words; positive A/B orbits for all moduli 1..80 (156,760 states); 43,679 odd commutator/inverse states; dyadic exponents 0..9 |
| `grammar_check.py` | complete finite tagged fixed points: 58/65 states, exact supports <=127, exact seed-bank delta {13,39,117} |
| `width5_check.py` | all labelled order-five prime-quotient classification; separate forward quotient-word kernel; 86-state fixed point; 49 exact supported values; literal H31 witness |
| `obstruction_check.py` | independent fixed-cycle class enumeration through order seven; class counts 1,1,6,35,353 and full birth histograms; literal H15 atoms; every scaling threshold 5..100000 |
| `ear_check.py` | own even-path network; 948 literal child checks; all strong parents through order five; cyclic5/7/9/11 cut tables, with five direct larger-child checks per parent |
| `tiny_extensions.py` | all 64 order-four quotient isomorphism checks; literal H19; GF(2) image-space rank hostile; six residual inequalities; exact 65/67 amplifier threshold |

Every listed executable was run normally and with `python -O`; matching
receipts are byte-identical. `normal.json` and `optimized.json` correspond
to `audit.py`; other names use their script stem, except `local_check.py`
uses `local.*`, `obstruction_check.py` uses `obstruction.*`,
`grammar_check.py` uses `grammar.*`, `ear_check.py` uses `ear.*`, and
`width5_check.py` uses `width5.*`.

## Corrections, inheritance and scope discipline

- The compression coefficient comparison was clarified before acceptance:
  approximately 1.94809 is below ordinary-core 2.73072 but above strong-only
  1.86400. The finite positive bank is inherited, not rerun here.
- A producer-only min-plus checker initially compared an uncapped product
  closure against an order-six census. It was repaired by truncating costs
  above six before comparison; no theorem statement depended on that error.
- The ordinary path-cover substitution law, same-H profile hostile,
  exclusions 7/21, strong floors and earlier capped methods are explicitly
  inherited. A repository search is not a literature-priority certificate.
- Arbitrary positive profiles, cut networks or formal affine states are not
  assumed tournament-realizable. Actual witnesses or constructive state
  induction supply realization where needed.
- Capping is allowed only for nonnegative integer formulas and preserves
  the threshold predicate, not exact values above the cap. Modulo equivalence
  is not exact-H equality. Reencoding loses context and embeddings.

The first attempted independent scratch write failed because `/tmp` was
full. The parent restored storage; the empty partial file was replaced and
all tests rerun successfully. No blocked action remains. No git mutation,
fetch or push was performed by this auditor.
