# Audit and reproducibility

2026-10-25. **PASS in the displayed scopes.** Mathematical proofs and finite
controls have separate roles. No computation here is a proof of an open
headline conjecture, and no Lean formalization is claimed.

## Independent reviews

- The planar proof was independently reconstructed in a polynomial jet ring,
  with denominator-cleared elimination and rational charts having finite
  poles. See [the detailed audit](anchors_audit.md): **28** independent
  checks, in addition to the author's **661** symbolic/finite controls.
  The actual scalar Jacobian has degree `12m-1`, below the individual
  degree-`15m-1` products but above every error term of degree at most
  `11m-1`. The constant and zero branches are included.
- The Collatz claims have [an independent native audit](collatz_audit.md):
  **917** narrow exact gates. The phase checker has **161,420** exact gates;
  the refill checker replays **256** positive sources, **64** residue
  controls, and **641** independent phase fixtures. These corroborate the
  all-source/all-h proofs; they do not extrapolate from a scan.
- The root reviewer checked the weighted deficit identity, positivity,
  full-gap clone obstruction, intrinsic weight-k anchor order, exclusion of
  other good-order types by vertex count, and both directions of the
  span-two matching bijection. The proof's strongness sentence was clarified:
  only the construction omitting exterior gaps can be strong; general
  records retain their source and sink. Counts and receipts are unchanged.
  Canonicalization does not use the theorem's anchor decoder. There are
  **178** independent feedback-DP controls and complete new-branch quotients
  through k=6; the matching bijection is exhaustive through k=7.
- The detour checker enumerates every labelled tournament through order six,
  without the inherited mixedness filter. All eleven success certificates
  are independently replayed with directed arc sets. The three failures
  have no black path at the target; invariance of path colour under its own
  K move rules out incoming black moves. A separate order-seven endpoint-four
  hostile is checked over all 5,040 permutations of that fixed target.
- The root reviewer checked all shape equivalences, including covariance
  positivity, semisimplicity, finite vertex actions and the infinite-order
  rational rotation hostile. The standard-library checker passes **80**
  exact controls. The Newton map is a semiring map only in its stated
  nonnegative-coefficient domain; a determinant still needs coefficients.
  Bounded proof search is polynomial in the **unary bound**, conditionally
  on P=NP, and says nothing about existence of conjecture proofs.

The new local normalization and scope clarifications do not alter the
hash-pinned predecessor packets. The older genus-six classification is
preserved; CURRENT-FRONTIER routes its new continuation. THM-4734 and
HYP-9348 were directly repaired because their omitted filters were
demonstrably false, with lineage in MISTAKES.md. No unresolved dispute was
used to keep those false claims current.

## Commands and frozen artifacts

From the repository root, with Python 3 and SymPy (the symbolic checks used
SymPy 1.14.0):

```bash
python3 -B 04-computation/shape_ports_20261025/verify_packet.py
python3 -B 04-computation/shape_ports_20261025/verify_packet.py --optimized
```

`--both` runs both modes. Individual commands appear in each proof. All gates
raise explicit exceptions, so optimized Python cannot disable them.
`MANIFEST.json` freezes every packet proof, checker and receipt except itself;
paths are repository-relative. The packet verifier compares parsed JSON and
checks every recorded file hash before running any named checker.

Git baseline: `2b96df316ec324bb8785c45d672f5562f5f2abfa`. All inherited mathematical references
in the proofs were read at that baseline. Primary literature links in
shapes.md were separately checked; they support their stated background,
not an RH/JC/Collatz/P-vs-NP equivalence.
