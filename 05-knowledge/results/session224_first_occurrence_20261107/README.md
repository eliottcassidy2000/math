# Session 224: first-occurrence tournament structures

**Status:** PROVED refutation; FINITE-EXACT complete birth atlas through order eight; recovered all-order H=33 formulas. Universal odd-value realizability remains OPEN.

For the actual number H(T) of directed Hamiltonian paths, define b(h) as the least nonempty tournament order attaining h. Track its complete birth fiber by

`(b(h), s(h), p(h))`,

where s counts self-converse isomorphism classes and p counts pairs of distinct converse classes. The fiber contains s+2p classes and s+p converse orbits. Converse reversal explains why individual orbits have size one or two; it does **not** force the entire fiber to be one orbit.

## Main findings

- **Smallest counterexample:** h=9 is born at order5 with **two self-converse classes**, not a converse pair. Both are strong, rigid, and have score sequence (1,1,2,3,3). [Elementary proof and verifier](proof/PROOF.md).
- **Prior 33 merged:** h=33 is born at order6 with **five classes = one self-converse class + two converse pairs**. Independent labelled enumeration recovers all five. At every order n>=6, prior work gives n²−5n−1 ordinary classes, exactly one self-converse class, and n(n−5)/2 converse orbits. These are recovered results, not new priority claims. [Provenance](PROVENANCE.md), [six-vertex replay](33/check33.py).
- **Complete finite ledger:** [birth_fibers.csv](census/birth_fibers.csv) has all320 counts attained through order8. Of these,285 have multiple converse orbits at birth. The exhaustive class totals are 1,1,2,4,12,56,456,6880; orbit mass independently proves completeness. [Census proof](census/PROOF.md), [receipt](census/receipt.json).
- **Constructive replacement:** retain a fiber of structures rather than a single representative or type bit. The ordered strong-component grammar gives exact multiplicative construction, minimum-order optimization, and converse symmetry. [Structural proof](proof/PROOF.md).

The320 observed birth orders are globally minimal because all smaller orders were exhausted. Values absent through eight vertices are not declared impossible. Globally7 and21 are excluded by prior theorems; the assertion that every other positive odd value occurs remains the Chan–Pak completeness problem. Existing audited finite coverage and stronger single-run claims are distinguished in [provenance](PROVENANCE.md).

## Reproduction and audit

From repository root:

```
python3 05-knowledge/results/session224_first_occurrence_20261107/census/verify.py
python3 05-knowledge/results/session224_first_occurrence_20261107/proof/check.py
python3 05-knowledge/results/session224_first_occurrence_20261107/33/check33.py
```

The census reads the repository's existing Session24 atlas; its exact input hashes are in the receipt. Repeat with `python3 -O`; checks remain active. The copied C++ augmentation provides a separate byte-for-byte order-eight rebuild. No strong-only filter or automorphism normalization is used for H.

## Research posture and scope

Anchor: first-birth fibers and the proposed dichotomy. Niche: reconcile the full33 multiplicity theory. Wildcard: test whether first-born structures alone retain enough information for future construction. The relevant concept board is scalar H, converse orbits, strong-component words, minimal order, and vertex-extension ancestry. The minimal hostile shows loss of structure even after retaining H, scores, automorphism order, and converse type. This packet does not claim a resolution of global realizability or LRC(14).
