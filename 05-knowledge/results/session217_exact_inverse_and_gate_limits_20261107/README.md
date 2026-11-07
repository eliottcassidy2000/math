# Session217: exact inverses and construction limits

**Session complete; Chan–Pak completeness remains OPEN.** The target is every odd tournament Hamiltonian-path count at least 23. This packet separates all-order proved results from finite controls. No new unrestricted finite-range record, eventual threshold, or full-spectrum proof is claimed.

## 1. Q/Z words are two counters, with a sharp logarithmic-square spectrum

[Proof and replay](paired/README.md), [independent review](paired/audit/REVIEW.md).

The left-reset Q and sink-side Z moves commute up to rooted graph isomorphism for every source, including transient sources. Every word is determined by its two letter counts. A length-k word therefore has at most k+1 rooted graph outcomes, not 2^k.

The exact closed forms prove that any fixed finite source bank, followed by arbitrary Q/Z words and the inherited width-eight compiler, supplies only O((log N)^2) distinct outputs at most N. A single rooted directed triangle supplies a matching Omega((log N)^2) family. Thus the bound is sharp, but this mechanism cannot cover an eventual interval of odd integers. Singleton and constant-axis cases are handled explicitly.

This changes the interpretation of Session216's still-open paired-source criterion: a synchronized Q/Z seed would yield a useful two-counter family, but would not by itself give binary interval amplification. The prior synchronization identities and explicit paired ray remain valid.

## 2. Adaptive uniform-module gates still reduce to multiplication

[Theorem and replay](bridges/README.md), [independent review](bridges/audit/REVIEW.md).

Consider any menu of tournament templates, finite or infinite, with each arbitrary input used exactly once as a uniform module. The selector may inspect the entire input graphs. If the output depends only on their Hamiltonian-path counts h_i, it must be a value C times their product supplied by a template whose input holes are singleton strong components.

The mechanism is simultaneous ordinal padding: it preserves each h_i but drives every unsafe template above any fixed target by a lower bound uniform in template size. Consequently even an unbounded adaptive menu cannot implement arbitrary-input addition such as x+y+1 or a nonzero affine shift. This is a scoped obstruction to uniform-module scalar gates, not to nonuniform cut constructions, stateful interfaces, or the three-core conjecture.

## 3. Complete fixed-core inverse in square-root time

[Algorithm, implementation and replay](inverse/README.md), [independent audit](inverse/audit/REVIEW.md).

For any fixed finite bank of nonempty strong cores, exact all-length membership in the two-sided-padding/star-closure template can be decided in

    O(sqrt(N) (log N)^2) arithmetic operations,

with O(log N) live node records using depth-first traversal. Short suffixes use low-prefix-bit recursion; long suffixes use a new reversed-suffix recursion. Both leave at most one exact candidate per branch. This is complete over all legal lengths, unlike a bounded balanced-band search. The complexity is an arithmetic-operation bound, not a unit-cost bit bound.

The integer-only C++ implementation uses the 80 Session216 states and explicitly accepts only 1<=N<=10^15. Independent checks compare complete recipe sets through 1,500 and 106 additional targets through 10,000. The final DFS executable, compiled with undefined-behavior sanitizer, agrees with exhaustive all-length support through 3,000 and passes seven additional controls. Sixty small literal strong graphs were counted by subset DP; larger counts use the independently checked exact construction formula.

The prefix-only recipe for 10^12+39 has 40 vertices, one fewer than its previously retained balanced-band witness. This is an individual construction, not a minimum-order theorem or an unrestricted coverage record. The permanent grammar omissions 47 and 63 remain absent at all lengths and remain realizable by other tournaments.

## 4. Low-run padding is a zero-density shortcut

[Proof and replay](run-entropy/README.md), [independent audit](run-entropy/audit/REVIEW.md).

For any fixed finite strong-core bank and every fixed 0<=rho<1/4, the targets having a recipe with at most rho log2(N) one-runs in the two padding words satisfy

    B_rho(X) <= X^(H2(2rho)+o(1)) = o(X).

A fixed run budget yields only polynomial-logarithmic growth. Thus bounded-run or sufficiently low-run restrictions cannot establish eventual coverage. General finite-state languages can allow many runs and are not excluded. The exact three-core reversal-product formula is retained with its logarithmic defect and explicit provenance to earlier one-apex work, not advertised as a new general mechanism.

## 5. Precise remaining proof obligation

The sufficient [Session216 three-core occupancy conjecture](../session216_overnight_constructions_20261107/three-core/PROOF.md) remains OPEN: every odd N>=65 should have a legal recipe. Together with the already verified small cases and separate 47/63 witnesses, this would prove the requested Chan–Pak target.

Faster membership does not prove universal occupancy. The current work rules out Q/Z-only finite-bank amplification, arbitrary-input scalar uniform-module gates, and low-run shortcuts in their stated scopes. A successful construction still needs target-dependent high-complexity recoding, another nonuniform returning interface, or a different all-target argument. No such argument was obtained in this session.

## Checkpoints and verification

1. [Q/Z commutation, sharp counting limits and adaptive module gates](https://github.com/eliottcassidy2000/math/commit/69d163a8680781b1fe042d76783c47ad1ceb783b).
2. The current final packet adds the complete inverse, run-entropy theorem, independent reviews and this synthesis.

Both checkpoints preserve concurrent main-branch work and use non-forced expected-head updates. All promoted results carry independent mathematical reviews, explicit computational scopes, replay commands and input hashes. Historical freezes retain their original labels and paths; current release reviews control the promoted claims. No full-spectrum proof, eventual threshold, positive-density theorem or unrestricted finite-range extension is claimed.
