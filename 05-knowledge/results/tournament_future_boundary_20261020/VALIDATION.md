# Validation and scope

The final normal and `-O` runs of `verify.py --max-order 100 --large-order 200`
both completed successfully and produced byte-identical `receipt.json`.
There are **34,705 always-active checks**. The source uses explicit
exceptions; no proof check depends on a removable Python assertion.

Independent implementations, implemented by the same author:

- An unmerged odd-partition Burnside sum through60.
- Direct formal clone substitution through18.
- Ascending versus descending elimination through24.
- Literal pair-orbit orientation enumeration and homogeneous-pair tests
  for every odd cycle type through order6.
- Full permutation enumeration of the Hamiltonian-path sets through order8.
- A subset DP for minimum feedback on all 88 tested no-singleton run
  codes through order10.

Additional external/inherited controls:

- All 77 entries of the retrieved OEIS A000568 table, orders0–76.
- The independent recorded 1,333-digit A100 hash, also used in incoming
  `602286dc6f`, agrees exactly.
- The previously generated pair-free core bank through order8.
- The inherited tribonacci path counts and half-size blue counts through14.

The 4,840 path codes at orders1–14 each pass exact labelled child
reconstruction, endpoint restoration, reflected-code restoration,
canonical-core preservation, and the sharp `3r-2` bound.
Every one of the 5,273 admissible three-letter words with1–8 runs is used
in the extremal core-size census. The feedback certificate separately
checks 1,676 triangle witnesses, disjointness of their edge sets, and
609 explicit attaining orders.

The table of maximum stored states is an exact implementation statistic.
It is not an asymptotic complexity proof. The full A200 output has an
all-order proof by the recurrence but does not have an independently
enumerated order200 object census. The fixed-point-free series is rational
symmetry mass; its coefficients are not claimed to count a subclass.

The source and proofs were self-audited against the stated hostile examples.
No independent-agent audit, Lean formalization, or external-priority claim
is made for this packet. Incoming work's independent reviews belong to its
own results and are not silently transferred to these new propositions.

`manifest.json` records hashes of the local proof, source, reference and
receipt files. It deliberately excludes itself. The two starting/incoming
commits locate the research context without treating historical reflections
as current truth.
