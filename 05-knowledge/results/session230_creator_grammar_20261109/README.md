# Session 230: first creators, prime kernels and positive cut growth

**Status: COMPLETE.** Independently audited proofs and exact finite certificates; normal and optimized replays agree. Chan–Pak
completeness and the all-order five-operation grammar remain OPEN.

[Read the synthesis](SYNTHESIS.md).

## Results

- **63 is an operation-relative orphan.** Its two first creators are W5
  substitutions with three transitive-pair modules. They need no new prime
  seed. The entire previously known Mersenne fan family has two prime
  kernels: W5 for even exponents, W7 for odd exponents.
- **Birth minima are not structural primitives.** The complete birth library
  is redundant under substitution (already at H=5), yet insufficient under
  substitution alone: prime W7 has the old count31. An independently
  certified census first finds eight prime noncreators at order seven and
  214 at order eight.
- **Controlled prime growth has an exact law.** The critical alternating-hub
  family satisfies H(W_(n+2))=4H(W_n)-n. The cited Schmerl–Trotter
  classification gives an all-order creator-rooted grammar with proper
  substitution, prime-preserving one-vertex growth, and only three explicit
  critical-family two-vertex lifts. This is not a minimal basis or a proof
  that every odd count occurs.
- **The omission proof object computes generation.** An explicit positive
  odd-cycle/even-path kernel gives every attachment as H+2 times a directed
  cut, with a minimum-cut birth certificate. A 31-word certificate excludes
  7 and21 for all apex-transitive tournaments without using their global
  exclusion. Prior directed-cut calculus is explicitly credited.
- **Recovered sharp connectivity boundary (THM-4811).** Minimum semidegree>=k implies
  connectivity>=k exactly up to the universal order threshold4k+1, sharp at
  4k+2. The floor-based near-ordinal-sum consequence itself holds at every
  order when its H-floor hypothesis holds.

## Proofs and audits

- [Prime kernels, fan decomposition and grammar](primes/REPORT.md) ·
  [independent orbit-mass and proof audit](primes/audit/REPORT.md).
- [Positive cut kernel and subclass omissions](cut/PROOF.md) ·
  [independent audit](cut/audit/AUDIT.md).
- [Review of inherited order-eleven grammar and sharp connectivity proof](audit/REVIEW.md).

The incoming S17 full independent audit applied its orientation and scope repairs and settled class-level creator genealogy through order eleven. We preserve that stronger result and its source. Session230 did **not** rerun the huge census; it independently checked all 97 exceptional witnesses and the complete order-eight class certificate. The earlier missing-audit finding is retained only as baseline history in our scoped review.

## Reproduce

From repository root, with Python3 and gcc:

    python3 05-knowledge/results/session230_creator_grammar_20261109/verify.py
    python3 05-knowledge/results/session230_creator_grammar_20261109/verify.py --optimized

The frozen [verification receipt](VERIFICATION.json) records all nine passing checks. The runner builds the inherited candidate-class generator in a temporary
copy, then executes nine checks against frozen receipts. The independent
canonical-form/orbit-mass certificate proves completeness through order
eight, with class totals 1,1,2,4,12,56,456,6880. All Python checks remain active
under optimization. No existing research files are overwritten by replay.

## Scope

No single representative per count is treated as the whole birth fiber.
No arbitrary insertion or unrestricted same-size arc rewriting is passed
off as a structural theorem. Count existence, class generation and minimal
basis are separate problems. First-birth strongness and creator ancestry
outside their verified ranges remain open; the all-order grammar has a
cited classification dependency, not a new independent classification.
