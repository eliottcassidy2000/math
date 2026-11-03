# Return at 6325, observer completion and what halving forgets

2026-11-03. **PROVED scoped / FINITE-EXACT.** Same-agent proof review;
independent checks and exact universes are recorded in [AUDIT](AUDIT.md).

* [RETURN6325](RETURN6325.md): twenty native moves return the actual source
  at its old root, without deleting old vertices; the next child is 18979.
  Names, root/spare responses, all five Q4 phases and operation clocks survive.
* [O1--O3](PROOFS.md): observer completion is an exact context-selection
  problem. The hostile eight-source fibre has completion infinity/2/2/1
  at incoming widths 1/2/3/4. Its response partitions refine 1->4->8 while
  invisible-symmetry orders are 54->9->9. Faithful symmetry is weaker than
  faithful source observation.
* The persistent-spare legal macro language has H-only response ranks
  1->4->10, with seventeen linear directions forever hidden. Keeping four
  named scalar readouts gives ranks 4->14->27 and complete response recovery
  within two macros; twenty-three extra scalar queries are necessary and
  sufficient on that rational span.
* [O4--O5](PROOFS.md): a one-vertex read head strictly orders adjacent twins.
  After structural halving, B fibres of size 2^a need B(2^a-1) single-port
  queries if only fibre sets survive, B*2^(a-1) if the unordered contraction
  hierarchy survives, and zero if its orientations are retained.
* [O6](PROOFS.md): the new neutral return has exact clone carry
  1342(2^m-2). Its first cloned root has defect 27 at every m>=2; later
  clone marks have larger defects. This particular family cannot remain
  1-ready or 2-ready by cloning/re-marking within that fibre alone.

Full ternary source closure, the complete Hamilton spectrum/A000568 and
Collatz convergence remain OPEN. The saved next source and decisive tests
are in [HANDOFF](HANDOFF.md); prior mechanisms and incoming work are credited
in [RECOVERY](RECOVERY.md).

Reproduce from repository root with Python, SymPy and g++:

```text
python -X utf8 -B 05-knowledge/results/return6325_observer_20261103/reproduce.py
python -X utf8 -B -O 05-knowledge/results/return6325_observer_20261103/reproduce.py
```
