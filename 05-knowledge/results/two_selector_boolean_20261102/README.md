# Boolean repairs and the cost of adding a marked operation

**PROVED all-size restricted-carrier theorem; complete FINITE-EXACT seven-vertex
observers; independently computed algorithms, same-agent audit.** The universal
Hamilton-path spectrum conjecture and full metagraph classification remain OPEN.

The requested two-edge trap test was already completed by incoming
[prime-inverse work](../prime_inverse_ports_20261102/MOVES.md):480 of768 traps
reach H17 in two flips, while288 require three. This session independently
replays the two-edge test and extends its actual repair cubes.

## Main advances

- **A minimum closed carrier at every size for one involutive selector.**
  With complement C and one specified legal involutive re-marking R, retain
  H, H after C, a three-way legality phase, and one additional conditional
  count after re-marking and complement. These data are sufficient and
  necessary for this command language. [All-size proof](CARRIER.md).
- **The17/71/33 example now has exact closure counts.** At order seven the
  observer grows from3008 states for C alone to10455 for C,R. Adding the
  second specified selector S gives15874 states. All32768 actual tiles are
  checked independently by endpoint DP and bounded distinguishing words.
- **An executable Boolean inverse, with its boundary measured.** Shortest
  legal command words reach H17 from148 tiles under C,R and152 under C,R,S.
  The latter needs at most four commands. The constructor records actual
  masks, selected paths and cumulative label transport. It does not assert
  that every source lies in this small operation basin.
- **Pairwise source laws miss actual target support.** Two positive laws on
  actual six-vertex tournaments have every one-edge and two-edge marginal
  equal and both mean H=17, but target-hit probabilities0 and1/5. Their
  squared errors are144/5 and136/5. Three-edge dependence remains observable.
  [Exact laws and shared-path certificates](JOINT_REPAIRS.md).
- **The supplied inspiration is separated into valid maps and failed claims.**
  The nine-cycle is a translated modulo-nine clock; magic-square and Lucas
  interfaces retain their arithmetic guards. The pasted spectra, projective
  pentad/anchor rule and MRB decay formula are refuted by exact tests and an
  all-size counting bound. [Proofs and recovered connections](INSPIRATION.md).

The shared mechanism is concrete: adding operations may expose formerly
invisible incidence, and correlating operations requires joint source data.
Closure under each operation separately does not imply closure under their
alternating words. A more accurate carrier improves exact reasoning; source
coverage is still an additional theorem needed for a universal inverse.

## Reproduce and construct

From the repository root:

```text
python -B -X utf8 05-knowledge/results/two_selector_boolean_20261102/experiment.py
python -B -X utf8 05-knowledge/results/two_selector_boolean_20261102/independent.py
python -B -X utf8 05-knowledge/results/two_selector_boolean_20261102/inspiration_audit.py
python -B -X utf8 05-knowledge/results/two_selector_boolean_20261102/constructor.py 868 --selectors 2
python -B -X utf8 05-knowledge/results/two_selector_boolean_20261102/replay.py
```

The last command returns the shortest word C,S,R,C with raw counts
93,31,31,31,17 and all source/target label orders. The same tile is unreachable
under C,R alone. Each verifier also passes with `-O`, with byte-identical
receipts. [AUDIT.md](AUDIT.md) records universes, source pins and scope.
