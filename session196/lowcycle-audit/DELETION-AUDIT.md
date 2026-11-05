# Independent audit: cycle-preserving deletion extension

**VERIFIED.** The five frozen extension files are preserved in `deletion-source/`, with SHA-256 hashes in `deletion_hashes.json`. The audited `DELETION.md` hash is `bd3980acc35daab753e6d79d365983eae48654c17f5cc202c308464ee8550073`.

`deletion_audit.py` uses the independently reconstructed 32-state Moore machine from the main audit, encodes deletion outcomes as a bitset, and explores the complete reachable lifted automaton. It reproduces **1,066 reachable states** and verifies that every failure to preserve cycle count has capped length at most 2 for `C=1`, or at most 3 for `C=2`.

The update invariant is exact: deleting the last letter leaves the previous current state, while deleting an earlier letter follows that deletion-state's transition on the new letter. Multiplicities of identical residual states are irrelevant because the assertion concerns existence of a deletion. Saturating length at 4 loses no information needed for either threshold. Complete finite closure, rather than a length cutoff, proves the all-order assertion.

Direct independent adjacency counting checks every listed irreducible word and each possible single-base-vertex deletion. The irreducibles are exactly:

- `C=1`: `2`, `30`.
- `C=2`: `120`, `212`, `230`, `302`, `310`, `320`, `321`.

A separate exhaustive pass through the independent cycle machine selects all **1,817** low-cycle words of base lengths 1 through 7. Independent adjacency/subset counting and sparse transfers agree on every selected word. The entire inventory, minimum values, and minimizing witnesses agree with `lowcycle_minima.json`. In particular, base-length-seven minima are **75** for one cycle and **115** for two cycles, and no selected word of base length at most six has `H=67`.

The second obstruction proof is valid: words of greater length reduce by cycle-preserving base deletions to length seven; induced-extension monotonicity then forces `H>=75`, or `H>=115` when `C=2`. Together with the finite census, this excludes 67 at every length in the defined low-cycle two-apex class, independently of the path-vector saturation argument. It does not exclude 67 for `C>=3` or outside this class.

The constructive reverse-deletion conclusion is also valid for rooted orientation words: repeatedly delete until irreducible, then reverse the deletions. It is not an unlabeled uniqueness statement or an inverse theorem for every admissible scalar path count.

Both producer extension scripts were rerun against the audited main source snapshot. All five output/source hashes still match the producer after reproduction. No unresolved issue was found.

Reproduction from this audit directory:

```
python3 deletion_audit.py
PYTHONPATH=source python3 deletion-source/deletion_certificate.py
PYTHONPATH=source python3 deletion-source/lowcycle_minima.py
```

Artifacts: `deletion_results.json`, `deletion.log`, `producer_deletion.log`, and `producer_minima.log`.
