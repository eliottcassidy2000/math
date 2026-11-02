# Tournament stability, residue tensors, and exact port recurrences

**New proved range: weighted tournament stability is sharp for every n>=42.**
The incoming frontier was48; this session resolves42..47, preserving the
fair-pair equality classification. Orders10..41 remain open. The extension
has a same-agent proof audit and exact rational certificates.

- [Stability proof](STABILITY_42.md): unused matching capacity replaces a
  heavy/light cutoff, then concavity reduces the remaining mass interval
  to endpoint checks. [Exact receipt](verification.json).
- [Tensor residue theorem](TENSOR_COSETS.md): a sharp full-jet inequality
  for every modulus m>=2, with two proofs. The Eisenstein proof needs a
  unit phase and a coset norm floor; it also works at ramified/composite
  moduli. Tournament Pfaffian data still guarantee only parity.
- [Recurrences](RECURRENCES.md): exact collision-block transitions for
  prime wheels, the guarded differential generating-function recurrence,
  and a cycle-marked Burnside functional equation with fixed ports.
  A000568 is reproduced through55 against the old packed-response engine;
  fixed one- and two-port counts are computed through55 body vertices
  (total orders56 and57 respectively). [Exact receipt](recurrences.json).
- [Recovery and source boundaries](RECOVERY.md): older5/11,56,57,
  primorials and the supplied arithmetic patterns, with explicit maps and
  failed transfers. [Concept board](BOARD.md), [audit](AUDIT.md),
  [next test](HANDOFF.md).

The81/27 control survives unchanged:67 shared labelled cyclic triples,
but1620 versus1728 regular binary attachments. Counts do not reconstruct
metagraph incidence. The prime-gap champion4 intervals are exactly
[131,139) and[491,547) in the checked range; the endpoints are ties.

Reproduce from the repository root:

```powershell
python -X utf8 -B 05-knowledge/results/tensor_primorial_ports_20261102/verify.py
python -X utf8 -B 05-knowledge/results/tensor_primorial_ports_20261102/recurrences.py
python -X utf8 -B 05-knowledge/results/tensor_primorial_ports_20261102/seal.py --check
```

Both checkers use explicit exceptions and give identical output under -O.
scout.py is floating exploration, not a proof or a finite exhaustive census.
The mathematical proof, not the number of checks, supplies each unbounded
quantifier. No resolution of a major open conjecture is claimed.
