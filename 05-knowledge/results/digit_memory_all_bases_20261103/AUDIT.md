# Audit and handoff

**Same-agent proof review, with separate exact implementations.** No external
referee or subagent was used. The statements are elementary all-size proofs;
the finite runs authenticate formulas, conventions and adversarial controls.
No empirical rate is promoted to universal convergence.

## Validity gates

* DM1 splits the B-supported divisor from its coprime complement before
  using an inverse; digital root zero is separated from positive sources.
* DM2 includes the empty suffix and the short-suffix empty-acceptance case.
  Its K satisfies BOTH B^K>=m and s|B^K. A logarithmic size bound alone
  would not ensure the valuation threshold when the base has mixed primes.
* DM3 excludes zero sources and allows leading-zero suffixes. The pure
  p-power-base case has a valuation quotient; the mixed-base construction
  cannot be reused when its coprime factor C equals one.
* DM4 distinguishes an empty legality guard, and proves the least translation
  period rather than calling it a minimal state count.
* DM5 observes the full ordered word, not just its odd weight; it uses the
  shortcut map with ordinary clock k+o. Formal continuation at ROOT is
  allowed. All source residues have positive representatives.
* DM6's minimality target is the source-tail residue t mod M0, equivalently
  n mod2^kM0 given the word. Preserving only n mod M0 could be cheaper.
  Ambient terminal modulus and number of attainable phases are distinct.
* TM1 uses an actual uniform module partition. The run skeleton keeps all
  disjoint visited paths, and its stabilizer quotient gives integer positive
  coefficients. The one-run contribution proving strict read separation is
  disjoint from the remaining skeletons; no cycle reserve is reused.
* The feedback cut is unique because BOTH old sides have more than r vertices.
  Strongness alone would not prove this reconstruction. Read contexts retain
  their selected block and local vertex-0 port, not just an unmarked count.
* Code distance measures adjacency errors; it is not a bound on the change
  of a Hamilton count after an arc error. No such noise claim is made.
* The same skeletons prove 3^M divides every storage H, so increased
  connectivity does not supply missing prime Hamilton counts.
* TM2 counts distinct modules adjacent to new vertices in each Hamilton
  path, not every module to which those vertices have nonuniform arcs.
  There are at most 2q such incidences. All old arcs are unchanged; arbitrary
  old-edge edits are explicitly outside the budget theorem.
* TM3 fixes attachment arcs independently of hidden bits. Its degree bound
  concerns the raw multilinear response, not arbitrary nonlinear decoding
  of an integer count. Matching small marginals gives equal means only.

## Exact finite universes

`check.py`:

* all 1,500 base/modulus pairs B=2..16, m=1..100: DM2 compared with full
  deterministic remainder-machine partition refinement;
* 153,000 divisibility sources in those pairs;
* all 8,424 affine-quotient parameters a=-6..6, b=-4..4, d=1..8, M=1..9,
  including empty guards and nonunit coefficients, with period found by
  literal exhaustive translation tests;
* 137,088 mixed-base positive pairs and 7,499 pure-base equal-valuation
  pairs, B=2..20, p in {2,3,5,7}, 1<=n<n'<=64;
* 150 complete Collatz periods, k=0..9, M in
  {1,2,3,4,5,6,7,8,9,11,18,25,27,54,81}: 267,003 source residues;
* 4,590 CRT precision cells and 4,080 adaptive source cells through k=8,
  with actual ordinary-clock replay and all-odd hostile controls.

`tournament_check.py`:

* exact seed and query unordered profiles by subset path DP and least-vertex
  cover recursion;
* 48 literal storage/read graphs: all four two-bit words, r=1..4 and both
  addressed readers; Hamilton counts compared with formula (3);
* all 7,760 deletion sets of size <=r for these storage graphs;
* full profile equality for the four literal r=2 twelve-vertex graphs;
* larger r=5,9,10,20 connectivity controls on explicit graphs, not exhaustive
  deletion censuses at those orders. Their general guarantee comes from proof.
* twelve coupled one-vertex extensions: all four two-block bit words and all
  eight three-block words. The prime 30431 and exact v3=1 at H=12573525
  give sharp boundary-budget controls. The quadratic mixed difference is
  312 and the cubic mixed difference is zero, as predicted by TM3.

`audit_independent.py` imports neither main checker:

* 15,300 complete modular multiplication orbits, m=1..50, multiplier1..12,
  every starting residue, compared with prime-power transients and periods;
* 168 smaller base/modulus machines checked by literal suffix-response vectors;
* 7,650 word/phase cases using rational affine composition and direct image sets;
* 30 adaptive-modulus minima, trying every smaller endpoint modulus and
  rejecting it by a literal phase collision;
* eight profiles by permutations and cut masks, deduplicating actual unordered
  path covers; all 120 seed relabellings reject an isomorphism;
* twelve wrapper/read counts by independently enumerating collapsed cyclic
  words, compared with the literal graph receipts.
* nonmemoized directed-path enumeration of all four twelve-vertex coupled
  readers; deterministic trial division proves 30431 prime; the separate
  seven-vertex sharpness component has exactly 93 paths by literal DFS.

Normal and optimized replay, receipt equality and documentation-check outcome
are recorded in the close-out section below. Reproduction commands are in
README.md. A fresh implementation is evidence against shared code mistakes;
it is not an independent mathematical referee.

## Next decisive steps and unchanged gaps

The most productive surprise was the positive full-profile read difference:
it allowed an all-size connectivity upgrade without searching large graphs.
The subsequent coupled-reader probe did leave the 3^M divisor ideal, giving
the prime 30431 and a sharp two-module-factor budget per added vertex.
Next, compatible output sources are needed to repeat such escapes.
For Collatz the adaptive terminal modulus is a correct refinement interface;
the decisive remaining issue is still paid, source-specific height descent.
Use the exact all-odd hostile at every attempted bounded-depth closure.

No intentional changes were made to other sessions or their theorem statuses.
The initial full checkout hit disk capacity; the session used a sparse
worktree and materialized the documentation checker's required tracked files
without changing their contents. The first documentation run failed on those
absent sparse paths, before they were materialized; that was an environment
failure, not a theorem or documentation correction.

## Close-out checks

All three normal/optimized pairs agree byte for byte. The documentation
checker passed after its required files were materialized. Relative packet
links were checked against this worktree or tracked HEAD; no broken link
was found. `git diff --check` is part of the final Git close-out.

Receipt SHA-256 values:

- `certificate.json`: `4e403b74aed901ae9c882a55038b92966db793f0785b43dae94107a2b6f5537f`
- `tournament_certificate.json`: `0443eac2864300c29791a435da7353c82139e88c35723360130867157f0740f8`
- `audit_certificate.json`: `4602fadc47b0fdaf6bf3a8a31999616f72e43e0d44f66b8d853d2769c53ff0f3`
