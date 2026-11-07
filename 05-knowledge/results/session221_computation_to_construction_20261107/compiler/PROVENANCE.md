# Machine and Rule 30 provenance for the Chan–Pak constructor

Pinned truth surface: `eliottcassidy2000/math@38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817`. All paths below refer to that commit, not the dirty shared checkout. `SOURCES.json` supplies exact paths, SHA-256 hashes, and immutable GitHub URLs. Repository theorem IDs below always mean their full slug in that manifest. Startup read AGENTS, START-HERE, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL, the Rule 30/tournament frontier and relevant CORE-PAPERS entries. No shared-repository file was edited by this research worker.

## Findings first

The strongest inherited construction is **not universality**: Session133 already has a 32-fraction program computing the ordinary Hamiltonian-path count of every pointed one-apex tournament, including the binary address needed for future updates. Session135 supplies actual priority-correct Turing translations and an exact positive-polynomial cap lemma. THM-4880 supplies the target nonnegative ribbon transfer. Session203 `saturated-profile-membership/PROOF.md` and Session216 `inverse/PROOF.md` already implement the entire finite cap-BFS/shortest-witness method for other tournament grammars. Their specialization yields the all-budget finite selector and marked local-edit observer in `THEOREM.md`.

Existing count-construction automata, fixed-budget certificates and spectral exclusions must remain credited. The contribution here is an explicit portable composition with honest termination, witness and edit interfaces, not a newly discovered spectrum range. No Rule 30 physical orbit has been mapped to a prescribed odd tournament count. No universality result proves Chan–Pak.

## Seven live source objects and transfer contracts

### 1. The marked tournament tape, including its two obstructions

**Source:** Session65 `session65-cryptids/REPORT.md` §3 (lines 141–218), within `05-knowledge/results/session65_cryptid_tournament_interfaces_20261026/`; accompanying `session65-k-write/REPORT.md` §§1–3 and `session65-palindrome-guard/REPORT.md` §§1–4.

- **Object/map:** encode a bit by `G0=TT3` or `G1=C3`, differing in one arc; encode a word by their ordered join `E(w)`. Canonical SCC order decodes the word, including zero runs. The exact isomorphism-edit metric is Hamming distance.
- **Preserved:** full finite tape with named window, head, control and program; a changed bit is one controlled marked arc rewrite.
- **Lost:** `H(E(w))=3^wt(w)` and triangle count retain only weight, not word order or machine time. `01` and `10` have the same count but different macros. A one-arc edit can take `2k+1` native machine steps.
- **Sidecar:** ordered window/origin, head/control/program, actual write guard. Not the autonomous full-order `K` dynamics. Literal full-order K cannot perform a one-bit write at tape lengths `m=2,3 mod4`. Pure-blue endpoints give a separate black-K obstruction; `E(w)` is pure blue exactly for palindromic `w`. Guards/padding remove these stated obstructions only, not the simulation gap.
- **Target/productive transfer:** replace the independent TT3/C3 blocks by the **ordered ribbon alphabet**. An alphabet bit at site `j` is now the orientation of `{j,j-s}`; its write is a real local arc flip and its observed H is given by the target profile transfer. This is the marked local-write map in THEOREM §4.
- **Cheapest decisive test:** equal-count words `000,001` at span one; toggle the middle bit. Counts become `1,3`. It is a literal same-H, same-write hostile against scalar dynamics, checked in `verify.py`.

### 2. The exact one-apex machine and binary sidecar

**Source:** `05-knowledge/results/session133_fractran_spectrum_refinement_20261101/path_counter/THEORY.md` §§1–3 (lines 10–190). It credits Session124's digital inverse and the October20 feedback-word formula, rather than claiming them anew.

- **Object/map:** pointed `(T,x)` with transitive `T-x` corresponds bijectively to its binary attachment word. For inversion weights `w(1)=1`, `w(d>=2)=2^(d-2)`, `c(b)=sum b_i(1-b_j)w(j-i)` and `H=1+2c`. Append maps are `0:(c,t,l)->(c+t+l,2t+l,0)` and `1:(c,t,l)->(c,2t+l,1)`.
- **Preserved:** ordinary H, odd-cycle count, chronological word address `Q=2t+l`, actual graph constructor. A **fixed 32-fraction FRACTRAN list** computes these from a sentinel word encoding; outputs are designated prime exponents, not the whole integer.
- **Lost:** H alone loses the address and next response; forgetting the marked apex also merges different pointed fibres. Runtime is unary/exponential in word length, not a fast binary algorithm.
- **Sidecar:** Q, previous letter and input bit order, apex mark and source word, control-prime/first-eligible invariant.
- **Target/productive transfer:** use this as the concrete machine-to-H positive control before attempting multiroot/ribbon constructions. Compile its controller into the actual TM/CA languages through item 3, but retain the graph word separately. The still closer `05-knowledge/results/feedback_word_hamilton_spectrum_20261020/proofs.md` §3 already proves the exact adjacent `10->01` jump formula, strict rearrangement, and the no-unit-gap theorem for strong interior swaps at order>=6. Novel claims must lie in new fixed-length gap spectra, constructive prescribed differences, optimal-size statements or physical arc-flip interpretation, not that local formula or single-swap Lyapunov law.
- **Cheapest decisive test:** `f(4)=f(6)=5`, but `f(8)=9` and `f(12)=11`. The scalar state fails under one appended zero. The fixed-H inverse already has a proved finite search `n in [(h-1)/2,h-1]`; bounded search here is not a new discovery.

### 3. The certified machine compiler and positive cap

**Source:** `05-knowledge/results/session135_fractran_tm_structure_20261101/machine_translation/PROOF.md` §§3–8; `acceleration/THEORY.md` §§2,4. THM-4858 supplies the exact TM-to-radius-one-CA wrapper (53 symbols for Leng's 25-state machine), using standard constructions.

- **Object/map:** ordered FRACTRAN -> valuation registers with first-eligible guards/rollback -> literal TM transitions; general finite register programs also have cancellation-safe two-fraction instruction encodings. Separately, nonnegative polynomial checkpoint maps descend exactly through `z->min(z,B+1)`.
- **Preserved:** configurations and halting at declared boundaries, exact native costs, bounded-output values under positive checkpoint maps.
- **Lost:** prime renaming can lose natural density; readout without the inert-factor sidecar loses integer identity. Microstep saturation loses countdown time. Compilation changes cost and machine size, not the underlying conjecture.
- **Sidecar:** actual boundary decoder, control/program, observer, primitive cost; positive rather than signed state coordinates.
- **Target/productive transfer:** specialize the already implemented Session203/Session216 capped inverse method to THM-4880's 27/119-dimensional nonnegative profile matrices. Session216 `05-knowledge/results/session216_overnight_constructions_20261107/inverse/PROOF.md` supplies the closest exact theorem, including shortest parent-pointer words, the enormous all-target bound and output pruning; it explicitly credits `session203/saturated-profile-membership/PROOF.md`. This guarantees a finite all-length target selector at every finite B. Retain predecessor words to emit actual tournaments. `THEOREM.md` gives the proof and precise complexity boundary. The two-sided capped context product also observes supplied local edits.
- **Cheapest decisive tests:** cap commutation against exact profile evolution and independent graph DP; reject signed `10-9` and decrement `12,18 -> 4,10` hostiles. The three finite fixed points in the release are executable, not inferred from large finite prefixes.

### 4. The Rule 30 observer quotient and exact transition test

**Source:** THM-3456, *left-permutive-trace-bijection-and-rule30-seed-boundary*, §1 and §2; THM-4013, *rule30-depth-fifteen-history-failure-and-adaptive-routed-repair*, §§2–3.

- **Object/map:** a left-permutive CA configuration maps to `(center trace, positive initial half-line)` bijectively; bare trace language discards the inverse left boundary. Separately, physical phase records are grouped by finite observations and tested for equal next observations.
- **Preserved:** full free-input trace cylinders or the explicitly retained finite observation, not the distinguished seed or next transition automatically.
- **Lost:** named seed and operation-compatible future. Rule30/Rule60 have identical free trace languages but radically different single-seed centers. Three histories plus carry still collide at physical records `1597/32967`; six histories retain ten finite mismatches.
- **Sidecar:** inverse boundary for seed membership; physical phase owner and routed chains for tested transitions. Adaptive closure through depth15 is finite-only.
- **Target/productive transfer:** transfer the **quotient test**, not Rule30 evolution: declare `q(u)`, then require every letter to preserve equality of q. For target construction, equality of current H is insufficient, while the capped profile vector is provably a right congruence. Marked edits additionally need right context.
- **Cheapest decisive test:** the actual one-apex `Q=4,6` collision, and the span-one edit hostile in item1. No theorem identifies Rule30 configurations with tournament constructors here.

### 5. The finite inverse-boundary monoid and the reset trap

**Source:** THM-4204, *rule30-debruijn-reset-and-dyadic-prefix-saturation*, §§2,6 (lines 83–167,512–571).

- **Object/map:** a binary **spatial** output word maps to a product of two four-boundary-state matrices. The resulting transformation monoid has41 elements; shortest resets are `0010,1010`. Cyclic predecessor count is the trace.
- **Preserved:** that inverse cyclic spatial count and boundary transformation.
- **Lost:** the marked temporal center after the unpointed product saturates at prefix length15. Rotation preserves trace, not product rank: `0101` has rank2, `1010` rank1.
- **Sidecar:** marked boundary state, operation/reading direction, and temporal/physical admissibility. A finite monoid alone is not a solver for a different observable.
- **Target/productive transfer:** capped **full block matrices** form a finite compositional monoid for ribbon counts. Only a reset proven to preserve the target readout under all relevant left/right contexts can be used as a constructor gadget. Current scalar or rank-one collapse is not enough.
- **Cheapest decisive test:** compare all bounded left/right context evaluations before declaring two block actions interchangeable. Signed/eigenvalue resets require separate exact checks; Session220's rank-two block contains a doubling eigenmode but its physical H also has a nonzero265-mode.

### 6. The lossless Cartier current tree

**Source:** THM-4210, *rule30-lossless-dyadic-block-current-cartier-tree*, §§2–5,8–9 (lines142–423,625–703); closest earlier source THM-3463, *rule30-mealy-section-suffix-parity-current-and-complexity-boundary*, §§1–2.

- **Object/map:** retain the exact physical block current `J`; form `F_n=qJ_n+J_(n+1)` and both spatial Cartier channels `E=C0 F`, `O=C1 F`. The two physical temporal children are the two temporal parities of **E**, not E and O.
- **Preserved:** every current coefficient, the marked 2-kernel node and physical child transitions when the full state is retained.
- **Lost:** E-only descent drops the transverse O channel. It is center-null at the marked decimated observation, yet re-enters nonlinear physical current generation. Infinite-state tree presentation is not finite-state automaticity.
- **Sidecar:** O plus initial current for reconstruction, physical-admissibility relation and bit-reading convention. Finite profile counts are not an infinite lower bound or a Rule30 prize proof.
- **Target/productive transfer:** transfer the **future-observable state principle**: a coordinate invisible to current H can matter after a local edit or continuation. Use exact left/right context responses, or the full nonnegative profile, before quotienting. A deterministic all-scale digit inverse would need a proved closed residual state, not only fitted small-bit outputs.
- **Cheapest decisive test:** compute suffix-response distinguishers for any proposed compressed constructor state. THM-4880 already supplies a14-observable span-two basis and irreducible13-dimensional nonempty-word space, so simply rediscovering the need for extra state is inherited.

### 7. The existing target transfer and corrected limits

**Source:** THM-4878, *graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses*; THM-4880, *ribbon-tournaments-bounded-span-back-arcs-give-a-finite-state-transfer-for-hamiltonian-paths-strong-connectivity-is-regular-the-span-2-series-has-a-rank-14-integer-representation-irreducible-on-nonempty-words-47-needs-span-3*; THM-4881, *reachable-states-of-the-span-2-ribbon-transfer-congruence-coset-of-index-2-16-exact-tripling-block-sparse-response-sets-and-residue-bias*.

- **Object/map:** bounded-cutwidth/spine or bounded-span ribbon words -> positive finite profile matrix products -> exact H; regular strong-connectivity sidecar; quotient by future H observables gives signed linear representations.
- **Preserved:** H under every append continuation; full profile retains fragment-boundary data. Product with the actual condensation DFA retains strongness.
- **Lost:** graph identity and word multiplicities after state merging; positivity in the minimal signed basis; arbitrary graph coverage at fixed span. Span2 excludes47, while span3 realizes it.
- **Sidecar:** actual witness word, graph construction, positive basis for saturation, root responses or condensation state when those predicates are wanted.
- **Target/productive transfer:** the finite-cap compiler and marked edit observer above. It gives a certified testbed for state-dependent control laws, not such a law automatically.
- **Cheapest decisive test:** reproduce THM-4881's bounded spectrum via a distinct positive cap representation, then use exact shortest witness extraction as a regression oracle for proposed local-control constructions. All budget99 coverage is inherited.

## Corrections and originality boundary

- `01-canon/MISTAKES.md`, **MISTAKE-701**: TM/FRACTRAN checker alignment overstated300 steps as covering300 when it covered66; the corrected literal audit covers300. Tournament GF(2) automata classify merged nodes only through5 and fail at every order>=6. Paley/QR-code facts were already THM-481. The GF(2) map is not an H-preserving replacement for the positive profile transfer.
- **MISTAKE-719**: rank(M-I)>=25 does not exclude unipotence; the repaired short-block no-go uses exact eigenvalue1 multiplicity. The14-dimensional representation is reducible; irreducibility applies to the13-dimensional nonempty-word space. These bounds do not exclude all longer or state-dependent blocks.
- **MISTAKE-720**: cut-state methods and saturation were already known; empty-word conventions changed published Hankel ranks; two-to-the32 coverage did not belong to the claimed spine family. No new finite-state method or large coverage range is asserted here.
- **MISTAKE-721**: invariant-lattice containment was not equality without the recursion. Tripling starts sharply after the third repetition; it does not supply additive digit control. The all-length `H<=99` fixed point is explicitly inherited.
- No new external literature search was performed. Originality means “not located in the bounded cited repository sources,” never “first in the literature.” Generic TM/CA simulation, nonnegative semiring saturation, matrix-product segment trees, and frontier DP are standard. The useful deliverable is their exact combination with this repository's physical tournament words, preserved predicates, lost coordinates and executable checks.
