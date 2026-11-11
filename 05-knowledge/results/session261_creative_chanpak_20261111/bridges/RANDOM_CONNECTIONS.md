# Random-source connection ledger

**Status: provenance and scoped transfer tests, not inherited-result audits.** Pin: f92e070445e43c59f1252f2979db564df88d7b60. Python's `random.Random(seed).sample(pool,n)` was applied to paths in their `git ls-tree -r --name-only` order. Sources were read with `git show`; their claims are accepted only at their stated scope.

## Wave 1: five genuinely uncontrolled result pulls

Seed **2611703**. Pool: all `.md` paths under `05-knowledge/`, excluding paths containing the literal substring `session26`; size **12,374**. Five sampled paths:

1. `05-knowledge/results/nonzero_preparation_20261103/README.md`
2. `05-knowledge/results/session187_named_response_20261104/gate_audit/REPORT.md`
3. `05-knowledge/results/session56_rmq_lca_connections_20261025/audit-forest-mirror/corollary/REPORT.md`
4. `05-knowledge/results/session178_prime_flow_honeycomb_20261104/primeflow/README.md`
5. `05-knowledge/results/routefree_suffix_gate_20261021/README.md`

### Transfers tested

- **Native preparation → apex cloning.** Map saved root to twin-apex module, retain actual graph and root rather than scalar H. Preserved: old graph induced on deletion. Lost: readiness and continuation sidecars. Cheap test: C3 gives 2^m+1. Outcome: amplification, not a small correction. The stronger all-order bound (3) in `PATH_COVER_INTERFACES.md` measures the entire newborn contribution.
- **Named gate occupied runs → alternating two-apex word.** Map an occupied gap to a path-cover bin; preserve the first/last colour constraints and exact number of runs. Lost: vertex names after quotient to word, harmless only because U is transitive. Test: actual DP at ABAB produces 29. Outcome: successful literal 27-state automaton and infinite closed family, not a scalar analogy.
- **Forest complement involution → domination/odd-cycle pairing.** Source uses an involution whose fixed point controls parity; target would need enough paired objects to reserve old paths. Map attempted: complement vertex choices inside a biclique. Preserved: parity balance of subsets. Lost: arbitrary-edge legality and Hamilton seams. Cheap test: an arbitrary bipartite incidence does not make complement pairs legal directed paths. Outcome: rejected as a universal bijection, repaired to the even-subset count E_i and factorial component-block theorem. Here the source's parity pattern becomes an explicit subset identity, and all joins are checked.
- **Prime-flow multiplicative cuts → graph components.** Source combines independent factors before comparing totals. Target map: expand ∏d_i−1 over the set of active domination components. Preserved: component product count. Missing sidecar: a legal order of corresponding path blocks. Cheapest test: two independent edges provide both orders, giving exactly 8 newborn weight; three provide 6 orders instead of the 4 needed. Outcome: the factorial surplus is the proof mechanism. No prime valuation theorem is imported into tournament counting.
- **Route-free source ownership → disjoint OCF packing classes.** Map an owned source to a uniquely recoverable rooted odd cycle; retain chosen vertex subset and component permutation. Preserved: no double counting. Lost if one simply sums overlapping biclique lower bounds. Test: overlapping selected subsets can generate the same cycle. Outcome: only the induced-component family is summed globally; the general biclique bound is stated individually.

Board update: the original domination injection is not the only useful map. Working with owned odd-cycle objects bypasses its global Hall obstruction on a non-chain class. Extra apexes can escape a restricted image even while amplifying the positive defect.

## Wave 2: unrelated canonical theorem pulls

Seed **2611704**. Pool: all `.md` paths under `01-canon/theorems/`; size **4,130**. Three sampled paths:

1. `01-canon/theorems/THM-840-hamming-five-continuation-congruence-boundary.md`
2. `01-canon/theorems/THM-3821-cubic-pseudoplane-rz2-odd-ladder-terminal-riccati-gate.md`
3. `01-canon/theorems/THM-2684-three-tooth-rail-envelope-diagonal-arrival-law-and-full-dilation-nilpotence.md`

### Transfers tested

- **Operation-kernel criterion → two-apex continuation.** Source observation O has a deterministic update only if its fibres remain fibres after the named operation. Target O=(n,H,Hdelete,d), operation=replace the root by TT₂. Cheap test: search a small exact collision, not a large count range. Result: rows (16,1,3,7,14) and (20,1,2,7,14), both root0, have O=(5,9,1,9) but outputs23 and29. Thus the actual theorem applies with no analogy gap. The repaired state for this operation is (a₁,a₂); an arbitrary clone-size query uses the whole skeleton profile.
- **Kummer odd ladders → parity-filtered clone runs.** Only the representation idea is transferred: constrain allowed exponents before summing. Map odd root multiplicity to odd clone-run length, preserving oddness of each directed cycle. Lost: all differential, valuation and Darboux structure; none is asserted relevant. Cheap test C3 reaches equality. Result: summing odd subsets gives exactly 2^{m−1}, producing the sharp amplification theorem. This is an elementary combinatorial theorem, not an implication of the Kummer result.
- **Rail support versus projection → multi-apex bin automaton.** Source warns that nonempty projected transitions do not imply legal joint paths. Target must remember middle-bin nonemptiness plus first/last colours simultaneously. Cheap test allowing an empty middle incorrectly adds illegal adjacent apex symbols to a two-run skeleton. Result: retain three-state bins and an explicit accepting set; the matrix recurrence is certified on that full carrier. No nilpotence or rail theorem is claimed for tournaments.

Board update: full boundary support, rather than a count alone, determines continuation. The two-apex family realizes an escape but remains sparse. The most productive cross-theme transfer was the method of preserving an operation's necessary memory; the most productive counting transfer was product expansion plus owned block permutations.

## Reproducibility and scope

The pools and random seeds record actual choices, including unhelpful draws; they were not retrospectively chosen to make a connection look random. Reading stopped at the relevant initial material (roughly the first 8–9KB for long files). No source was used beyond that inspected scope. `path_cover_probe.py` uses a separate seed **2611705** for new component-block hostile controls, not source selection.
