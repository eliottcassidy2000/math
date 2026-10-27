# Session 74 — Leveraging eleven without conflating its roles

Source pin: `4959b4933cba7e6561b9f19b1d68278ecf6cb2f5`.
**Checkpoint 1: independently reviewed proof plus exact certificates.**
The user's `{2,3,11}` is recovered with its original integer-bracket mechanism;
a separate Paley11 thread now has an exact six-move forest-frame bridge.
Repeated appearances of 11 are not treated as a common mathematical object.

## The new result

For every tournament on `n>=2` vertices, let
`L(T)=sum_i |s_(i)-i|/2` for sorted outdegrees, and let `beta(T)` be the minimum
cycle rank of its actual backward graph over all rank orders. Then

    beta(T) >= max(0, L(T)-n+2).

Thus every forest-framed tournament has `L<=n-2`. A legal native path
complement reverses precisely the arcs outside its named directed Hamilton
path; the inherited score update changes L by at most one. Paley11 starts at
L=15, while an eleven-vertex forest frame has L<=9. Six explicit legal moves
attain that bound:

**The exact native distance from Paley11 to the forest-framed class is 6.**

The certificate attains intrinsic cycle ranks `5,4,3,2,1,0` at depths 1–6,
with L descending from 14 to 9. These are exact global minimum cycle ranks at
each stated depth. At depth zero Paley11 has inherited beta=10, not 6.
The particular witness's feedback numbers are not claimed to be the global
feedback minima at every depth. Free relabelling and converse preserve the
lower bound. The distance to the transitive class remains unresolved here.

- [Proof, scope, six paths and sharpness](session74-paley/REPORT.md)
- [Independent proof audit and controls](session74-score-audit/REPORT.md)
- [Labelled certificate](session74-paley/certificate.json)
- [Producer receipt](session74-paley/normal.json) and
  [independent receipt](session74-score-audit/independent-normal.json)

The bound applies to actual backward graphs, not merely an equinumerous family
of forests. Edge count `<=n-1` does not imply acyclicity. Nor does the score
bound imply `f<=n-2`: the final forest witness has exact feedback number 10.

## Recovering the earlier elevens

[The inheritance map](session74-inheritance/INHERITANCE.md) recovers and scopes:
odd-square bracket escape `{2,3,11}`; the marked square-sum port 11 in Q25;
the golden period-five quotient of size 11; Terras period eleven; the finite
valuation-potential obstruction; Mills' `(2,11)` pair; Camion residue tests;
and the unrelated count of eleven seven-vertex spectra.

The same-bracket integer multiplication classification and its real-multiplier
extension were already proved. They are inherited, not rebranded as new.
The original decoder checker is explicitly replayed from its exact source
bytes. All other historical source programs, if any, remain inert snapshots.

## Verification and limits

Run from this directory:

```bash
python3 -B verify_packet.py --replay
python3 -O -B verify_packet.py --replay
python3 -B verify_provenance.py
python3 -O -B verify_provenance.py
python3 -B test_packet_corruption.py .
python3 -B test_provenance_corruption.py .
```

Python replays use the standard library. The independent forest enumeration
also requires `g++` (C++17); its wrapper compiles only in a temporary directory,
runs optimized and undefined-behavior-sanitized variants, and preserves frozen
receipts. Native `reproduce.sh` is retained as provenance but is not the
immutable replay entry point, because it writes outputs in its own directory.

[Verification scope](VERIFICATION.md), [dependency readback](DEPENDENCY_ORIGINS.checkpoint1.json),
and the checkpoint inventories distinguish proved arguments from bounded exact
checks. No formal Lean proof, literature-priority claim, LRC(14) closure,
Collatz integrality result or universal decoder is asserted.
[The original plan](session74-root/SESSION_PLAN.md) is explicitly a research
plan; its uncompleted candidates are not checkpoint results.
