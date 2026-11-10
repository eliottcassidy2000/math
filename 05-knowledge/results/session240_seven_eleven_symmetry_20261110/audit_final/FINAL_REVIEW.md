# Final adversarial synthesis review: Session240

2026-11-10. **PASS for the mathematical synthesis and stated finite scope.**
No mathematical blocker was found. Verification-mode wording identified below
must track the final scripts; this review is not a replacement for either
independent census audit or the final replay.

## Reviewed artifacts and precise conclusions

Reviewed the packet README, gaps/PROOF.md, paley/REPORT.md, mixed/PROOF.md,
provenance/README.md and SOURCES.json, both independent AUDIT.md reports,
REPRODUCE.md, replay.py, THM-5100/5101/5102, and navigation S240.

1. **Elevenfold gap:** the prime-action reduction, sevenfold threshold2163,
   core minimum92411, path-cover kernel, and free ordered-path action are
   correctly inherited. The new exact elevenfold F2 table supplies the sharp
   strong non-core threshold1848759. Multiple moving orbits use the inherited
   product bound's equality characterization, not the inequality alone;
   oddness and divisibility give the additive22. One moving orbit is a module;
   a quotient triangle lifts to an induced inflation. Vertex insertion is only
   non-strict monotonicity. The proof correctly leaves ambient equality open.
2. **Matching restoration:** the objects are point/block bijections satisfying
   skew incidence, acted on by the side-preserving design group. Distinct
   blocks make the point action determine the block action. Full matching
   enumeration, rather than a quotient of group orders, proves the claimed
   complete fibres. Ordered paths must remain tagged by their matching.
   The induced free action gives the same quotient9/1729 and no stronger
   raw-H denominator. PSL isomorphism identifications remain inherited;
   group orders alone are not treated as their proof.
3. **Joint order77 carrier:** every cycle length divides77. Through order19,
   the only possible faithful profiles are (7)(11) at18 and (7)(11)(1) at19.
   The generator acts transitively on the77 cross pairs, forcing uniform
   arcs. At18 the quotient is non-strong; at19 strongness forces a triangle.
   Conversely every stated inflation has the actual order77 action. The
   six allowed run words give 3ab+xb+ay+xy. Simultaneously attainable positive
   minima yield3049061631 and exactly6*10*2=120 labelled fixed-action minima.
   This proves the minimum-order stratum only. No assumption that group
   order divisible77 implies an element of order77 is made.
4. **Recovered domains:** the synthesis distinguishes affine and projective
   actions, fixed tournaments and incidence designs, rational Collatz labels
   and integer/positive predicates, binary image codes and original matrices,
   and full norm circles versus marked knight subsets. It explicitly retains
   the physical-packing cited premise and the all-period theorem's inherited
   status. No global Collatz, packing, Paley extremality, or Chan–Pak result
   is inferred from the repeated numbers7 and11.

## Source/access and inexpensive independent checks

`source_check.json` records read-only verification of all13 listed source
blobs against accessible commit
`bdaec8654d7b225cb5b270331c19c2e2b4f1332b`. Every path resolves to exactly its
recorded blob, and every blob is accessible. The pin does not depend on the
old local5806 commit. This is an object-integrity check, not an external
literature-priority search or an independent proof of every cited source.

`checks.json` records additional independent checks: exhaustive run words
through length8 (there can be at most six runs because z occurs once), all
order77 cycle profiles through19, the77-element cross-pair orbit, mixed
minimum arithmetic, and existence of every relative Markdown target in the
packet at review time. Larger enumerations are supported by the existing
independent design and gap audits; this final review did not duplicate them.

## Corrections raised before finalization

- REPRODUCE.md initially stated that explicit failures survive -O, while
  audit_design/verify_independent.py, paley/verify.py and paley/audit.py used
  assertions for mathematical gates. Their ordinary runs remain valid.
  An optimized run with disabled assertions plus saved-output comparison is
  not a rerun of those gates. Either convert them to explicit checks and
  replay, or limit the wording. This was reported to the coordinator.
- THM-5101 should explicitly say **side-preserving** design group, as its
  source proof already does. This excludes incidence correlations swapping
  the two sides and makes the general action statement fully unambiguous.
- The README and S240 still said final replay in progress when reviewed.
  Promote this only after the coordinator confirms the final replay; the
  reviewed saved normal and optimized receipts both report PASS.

## Publication-safe boundaries

The mathematical packet and wrappers consistently distinguish proof with
finite-exact inputs, inherited claims, labelled orientations, and actual
subgroup actions. They make no literature novelty claim. These are ordinary
exact programs plus human proofs, not formal proof-assistant certificates.
The independent design audit covers Paley7/11/19; the producer's separate
brute/all32 control audit covers its broader order11 census. Do not blur those
scopes. Remaining open problems include arbitrary odd target realization,
full ambient elevenfold equality, the all-order77 minimum, general matching
fibre classification, and global Paley Hamiltonian extremality.

## Closure after corrections

**All review corrections closed; publication-safe mathematical PASS.**
Re-inspected the current files after correction. THM-5101 now explicitly names
its side-preserving design group. AST inspection confirms no remaining assert
statements in paley/verify.py, paley/audit.py, or the independent design verifier;
their check helpers raise RuntimeError. Independently reran the design
verifier's inexpensive deliberate-failure gate in normal and optimized modes:
both exited1 with RuntimeError. The saved independent design normal/optimized
outputs agree after removing only timing fields. Its audit now accurately
states optimization survival and links its receipts.

The coordinator reports fresh successful combined normal and optimized replays.
Inspected both updated REPLAY receipts: each records PASS for all seven Python
jobs and the separate nineteen-vertex C++ witness job. The full suites were not
needlessly duplicated by this final reviewer. Earlier reproduction concerns
are therefore resolved, not outstanding mathematical qualifications. Final
session/publication status is appropriately reserved for the coordinator's
publication completion; this review does not certify a push or remote state.
