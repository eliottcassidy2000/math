# Independent adversarial audit: native descent and positive interfaces

**Verdict: PASS within the stated mathematical scope. The endpoint-domain clarification is resolved in the current producer text. No Collatz or LRC(14) proof results.** No producer files were modified; no publication, external fetch, or raw git-blob fetch was performed.

## Sources read directly

- `../descent/POSITIVE_ATOMS_AND_INTEGER_DESCENT.md`, including the added §4a bounded-reset compiler.
- `../positive/POSITIVE_INCIDENCE_AND_ATOMS.md`.
- `../../session100_reframes/controllers/RETURN_CONTROLLER_TRAP.md` and its actual `check_return_controller_trap.py` source.
- `../../math/05-knowledge/results/collatz_atomic_rank_interface_20261011.md`, especially A1–A4.
- `../../session103_words_necklaces/triangularity/AVOIDING_WORD_TRIANGULARITY.md`, including its actual all-k argument.

The independent checker imports only Python standard-library modules. Reading producer/controller scripts is not importing or reusing their routines. Historical controller receipts were not represented as independently rerun.

## Mathematical review

### 1. Literal native rank

For active states, P>=3 implies B(P)>=1. Weak contraction makes B nonincreasing. Branch contraction P'<=rho P makes B(P')<=B(P)-1: B(P)-1 is a nonnegative admissible exponent for P'. A leaf subtracts one queue credit; a branch gains two credits but loses at least three energy-rank units. Thus every active-to-active native edge lowers 3B+Q by at least one. Direct terminal edges lower the positive active rank to zero. The strict threshold in B causes no boundary/off-by-one error; the displayed 163,162,38,4,0 example checks exactly.

The potential hypothesis is authenticated by the actual (3,3) controller, not by a statement about arbitrary prescribed valuation words. Small native inputs 1,3,5,7,9,11 are not branch inputs, while the extremal branch 13 with valuation 3 attains 27/28. Branches with larger valuations are covered by the inherited potential argument. The (5,5) 83-cycle repeats its entire native/queue state and forbids a well-founded strict rank, including along an infinite history-augmented run.

### 2. Atomic extraction and aggregate budget

The strict inequality A_K<w(x) is necessary; equality allows x to survive. The proposed uniform truncated budget yields tau(x)<=floor(M/w(x)) by retaining the original source atom in every nonnegative cost. Summing the displayed aggregate drift telescopes to the truncated cost identity. These implications do not require conditional drift; conversely, the hypothesis is not supplied by termination alone for a fixed prior. The countdown example with masses 2^-j and times 2^j correctly separates the two claims.

Existence of the hitting-time rank is a deterministic termination equivalence, not a recipe proving termination. The effective compiler is explicitly conditional on a proved effective bound. Source reweighting/conditioning and signed cancellation are correctly excluded. A fixed positive atom does not import Haar almost-everywhere results onto integer sources. The randomized caveat is appropriate.

### 3. Carry/source preservation and reset compiler

Direct expansion verifies both carry recurrences and L'=bL+2^D[(2^a-b)n0-1]. L is exactly 2^D(n0-n_J), hence its sign tests the immutable source rather than the preceding endpoint. For positive prescribed valuations the modulo 2^(D+1) cylinder enforces exact chronology, not just final divisibility. Independent replay tests both cylinder representatives and neighboring odd nonmembers.

The contracting 5-word (1,1,5) fixes 13 and pays zero source debt; the 3-word has rational fixed point 19/101, not an integer counterexample. No denominator rank is asserted. The auxiliary 3-return from 3 to 5 terminates but has L=-4; it is a decisive hostile to blindly paid resets.

The added §4a compiler is valid: F(n)-F(n') includes B(n')+1, so it pays the entire refreshed rank at any strictly smaller label. Internal edges retain their source label. Fresh-rank bounds, nonterminal edge availability, and genuine terminals are expressly assumed. This is conditional compilation, not evidence that raw controller returns satisfy its false source-descent premise. The independent checker also tests resets from nonzero residual ranks, a stronger finite inequality than needed.

### 4. Positive cone, graph alternative, and endpoint scope

The compact convex-hull minimum-norm argument proves the finite-dimensional strict separator theorem as written. Rational approximation is legitimate because all finitely many inequalities are strict; integer scaling then works. This replaces only a zero/nonzero predicate, not state reconstruction.

A nonnegative incidence circulation contains a directed cycle. A topological ordering gives a strictly positive edge functional. K4 has incidence rank three independently of orientation, yet exactly 24 of its 64 orientations are transitive and positive-detecting. Signed kernel dimension three is not in conflict with cone injectivity. The wedge 1∧f faithfully records potential differences; this is no identification of six arbitrary edge coordinates with three retained coordinates.

The Session103 telescoping formula agrees with the source's all-k theorem. Column matching by canonical cycle cuts must accompany endpoint-weight ordering (the inherited theorem supplies this). Cyclicity alone is correctly not promoted to observer singularity for arbitrary digraphs.

**Domain clarification resolved (current §5 re-read):** the endpoint section now explicitly assumes a finite nonempty set of distinct positive integer speeds. This correctly excludes empty V, for which Safe(V) would be the whole circle with no endpoint tags, and zero or negative speeds outside its enumeration. Under this explicit domain, the endpoint equivalence is correct: feasibility is constant within endpoint cells, and closed feasible arcs or isolated feasible points have a feasible endpoint. Birth/filter is a disjoint union by owner tags, not positions. Rejected old tags cannot revive; the new-speed birth term means accepted counts need not decrease. The 1,...,13 fixture checks exactly: 182 candidates, 12 accepted tags at six isolated locations.

**Publication-link clarification resolved:** the Session103 relative link uses `session103_words_necklaces_20261029`, the intended exact publication target confirmed by the coordinating agent. The source actually inspected locally is `session103_words_necklaces`. No local-path discrepancy is being treated as a mathematical objection.

## Novelty assessment

October11 A3 already gives the energy/credit rank 3b(E)+k and A1/A2/A4 already handle atomic integrability and controller-work expectations. Session104's 3B(P)+Q is a correct specialization, not a new abstract integer-ranking theorem. Its substantive local interface improvements are an explicit rank on **each native accelerated step**, source-marked carry debt, and a clear budget-to-bounded-search statement. The §4a finite-prefix-sum reset compiler is an elementary conditional assembly, not a missing source-payment proof. Positive-cone separation and DAG potentials are classical; the K4/observer/atom connection is useful packaging with its information-loss boundaries retained.

## Reproduction and evidence limits

```sh
python3 -B session104_descent_exterior/audit_descent/check_independent.py > session104_descent_exterior/audit_descent/receipt.json
python3 -B -O session104_descent_exterior/audit_descent/check_independent.py > session104_descent_exterior/audit_descent/receipt.optimized.json
cmp session104_descent_exterior/audit_descent/receipt.json session104_descent_exterior/audit_descent/receipt.optimized.json
```

Both modes PASS and receipts agree byte for byte. All guards use explicit exceptions, never `assert`. Scope includes 3,072 exact native rank edges, all 64 K4 orientations, 1,020 valuation words (length 1–4, valuations 1–4, multipliers 3,5,7), exact source/carry/debt replay, 186 endpoint cell decompositions, 558 endpoint insertions, aggregate budget controls, bounded-reset controls, and ten malformed-input rejections. Receipt counts distinguish individual identities from fixtures. These finite checks support the audit but do not prove the infinite claims; those depend on the arguments reviewed above.

The endpoint-domain correction changes no formula or numerical claim; the checker already enforced the corrected domain. No numerical rerun was needed for this report clarification.
