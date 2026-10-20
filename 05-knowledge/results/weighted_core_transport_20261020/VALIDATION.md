# Validation and finite scope

The new general statements have proofs in the three companion notes and
README. Computation checks their stated finite consequences; it is not
substituted for the proofs. No Lean formalization or novelty claim is made.

## Reproduction

From any working directory:

```text
python <repo>/05-knowledge/results/weighted_core_transport_20261020/audit.py
```

The audit runs all three scripts normally and with `-O`, rejects assertions
in their source, checks successful exits and empty stderr, and compares
literal stdout with the retained `.out` files. `--refresh` updates outputs
only after the two runs agree. The scripts use Python's standard library.

Two frozen utilities from `metagraph_growth_20261020` supply independent
small tournament canonicalization, class BFS, path cubes and homogeneous
chain decomposition. Those utilities were committed before this packet.
Their imports do not execute their main routines. New core generation is
by alternating attachments; the weight implementation independently starts
from the prior arc-flip census.

## Exact universes

- **Growth:** every one of 33,866 attachments to every labelled parent of
  orders 1 through 5; clone dichotomy, surviving old pairs, component
  equation, alternating gate, and full labelled transition histogram.
- **Structural recurrence:** complete core and weighted class banks through
  order 8, including 6,880 ordinary classes. Every individual core receives
  symmetry-corrected root mass exactly one. Labelled orbit masses sum to
  2^binom(n,2). Complete fixed-path cubes through order 7 independently
  recover exactly the same class sets.
- **Labelled mass:** core-size kernel and signed-binomial formal-series
  inversion agree through order 20. Each core-size stratum agrees with the
  ordered-chain substitution count in the generated bank.
- **Weight allocations:** all pair-free cores through order 7, all positive
  weight orbits through total weight 12: 208,126 orbit cases. Coordinate
  increments, direct compositions and cycle-type coefficient formulas
  agree. The primitive C3 number-theory formula is separately checked
  through total weight 60.
- **Converse:** anti-automorphism coset counts reconstruct self-converse and
  converse-merged counts through order 7, independently of the ordinary
  weight counts.
- **Cycle index:** chain substitution and inverse through weighted degree
  7; rooted differential identities through degree 6. These are formal
  series, not analytic convergence claims.
- **Path allocations:** 2,785 multiset shuffle rank/unrank words of length
  at most 7 over at most 3 consecutive labels; all 1,100 fixed-path cube
  states through order 6; 2,200 dilation controls; 1,109 path-root controls;
  5,164 legal endpoint clone insertions; 1,098 closed/open controls.
- **Hostiles:** the explicit equal-inversion allocation pair, equal-group
  but different-cycle weight response, equal-full-cycle-index but different
  growth response, and transitive-complement core growth through order 100
  in the allocation script. The last family also has a general proof.

The weight totals above 7 are restricted to cores of order at most 7;
they are never presented as full A000568 values above that order.
The growth script's complete order-eight census is separate.

## Independent review

Three parallel research lanes reviewed each other's claims. The allocation
proofs were reviewed by the root and the weight-counting lane. The growth
dichotomy, transition kernel, symmetry weights and cycle-index recurrence
were reviewed by both other lanes. A notation collision was repaired:
unlabelled core count is written |C_N| in the weight recurrence, while
the inheritance note uses c_n for labelled core counts.

The root additionally generated the order-eight bank independently and
found the full-cycle-index growth hostile. This is stronger than merely
testing core order or automorphism group order.

An incoming pre-audit theorem's affine-palindrome test gap is recorded in
`inheritance.md`, section 10. Its all-order alphabet claim is not used as a
dependency of this packet.

## Artifact integrity

`manifest.json` records SHA-256 after CRLF-to-LF normalization for each
packet source, note and retained output, excluding the manifest itself.
No cache files are retained. The normal and optimized replay output is
saved as `audit.out`.
