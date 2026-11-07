# Session217: exact inverses and construction limits

**Research ACTIVE; Chan–Pak completeness remains OPEN.** The target is every odd tournament Hamiltonian-path count at least 23. This packet separates all-order proved results from finite controls. No new unrestricted finite-range record, eventual threshold, or full-spectrum proof is claimed.

## 1. Q/Z words are two counters, with a sharp logarithmic-square spectrum

[Proof and replay](paired/README.md), [independent review](paired/audit/REVIEW.md).

The left-reset Q and sink-side Z moves commute up to rooted graph isomorphism for every source, including transient sources. Every word is determined by its two letter counts. A length-k word therefore has at most k+1 rooted graph outcomes, not 2^k.

The exact closed forms prove that any fixed finite source bank, followed by arbitrary Q/Z words and the inherited width-eight compiler, supplies only O((log N)^2) distinct outputs at most N. A single rooted directed triangle supplies a matching Omega((log N)^2) family. Thus the bound is sharp, but this mechanism cannot cover an eventual interval of odd integers. Singleton and constant-axis cases are handled explicitly.

This changes the interpretation of Session216's still-open paired-source criterion: a synchronized Q/Z seed would yield a useful two-counter family, but would not by itself give binary interval amplification. The prior synchronization identities and explicit paired ray remain valid.

## 2. Adaptive uniform-module gates still reduce to multiplication

[Theorem and replay](bridges/README.md), [independent review](bridges/audit/REVIEW.md).

Consider any menu of tournament templates, finite or infinite, with each arbitrary input used exactly once as a uniform module. The selector may inspect the entire input graphs. If the output depends only on their Hamiltonian-path counts h_i, it must be a value C times their product supplied by a template whose input holes are singleton strong components.

The mechanism is simultaneous ordinal padding: it preserves each h_i but drives every unsafe template above any fixed target by a lower bound uniform in template size. Consequently even an unbounded adaptive menu cannot implement arbitrary-input addition such as x+y+1 or a nonzero affine shift. This is a scoped obstruction to uniform-module scalar gates, not to nonuniform cut constructions, stateful interfaces, or the three-core conjecture.

## Continuing positive route

The concrete sufficient target remains the [Session216 three-core arithmetic conjecture](../session216_overnight_constructions_20261107/three-core/PROOF.md): every odd N>=65 has a legal recipe in its 80-state grammar. A new complete inverse algorithm and its independent audit are being finalized. An inverse decides membership for each specified target; universal occupancy remains a separate proof obligation.

Both released packets carry independent mathematical reviews, self-contained controls and exact input freezes. Historical frozen manuscripts preserve their original status labels; the current reviews govern the results above.
