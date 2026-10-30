# Ramsey boundaries with retained ports and symmetry

The main proved implication is that any hypothetical `(3,3,5)`-avoiding
coloring on56 vertices has a color-preserving automorphism group whose
order divides **224**. Restriction to one13-vertex neighborhood is faithful
after fixing its root. Every29-vertex blue neighborhood also has exactly
**104 red,104 green,198 blue edges**. These constrain a complete search;
they do not yet exclude all56-vertex colorings.

The bounds in the checked 2026 revision of
[Radziszowski's survey](https://www.cs.rit.edu/~spr/ElJC/ejcram18.pdf) are

    45 <= R(3,3,5) <= 57,
    61 <= R(3,3,6) <= 91.

Both exact values remain OPEN. The survey upper91 is credited to
[Boza](https://arxiv.org/html/2603.10851v1). **Update2026-10-30:** Dauth's
new preprint gives R(3,3,6)<=86; its exact certificate was replayed in
[N1 of the follow-up source audit](../ramsey_neighbourhood_lifts_20261030/NEIGHBOURHOODS.md).
The old91 boundary is superseded; conditional structures near that
boundary are not current existence targets. No Ramsey improvement is
claimed on behalf of this repository.

## Results and certificates

| Status | Result | Where |
|---|---|---|
| PROVED with CITED inputs | Single-port faithful restriction and group order dividing224 | [P1-P5](PROOFS.md) |
| PROVED | General edge-budget transfer and residual fixed-point obstruction | [P2/P5](PROOFS.md) |
| FINITE-EXACT | Seven17-vertex ports give67 rooted and1065 ordered two-root types | [verification.json](verification.json) |
| FINITE-EXACT | Two27-vertex local witnesses, including one with nonuniform margins | [AUDIT.md](AUDIT.md) |
| PROVED identities; FINITE-EXACT controls | Labelled, weighted and rooted attachment recurrences | [RECURRENCES.md](RECURRENCES.md) |
| FINITE-EXACT | GF(16) colored control, Paley17, and the81/27 hostile pair | [verify.py](verify.py) |

## Inheritance and concept board

The anchor is the requested Ramsey boundary; LRC(14) remains the repo's
open long-term anchor. The niche is the earlier Singer57/E7 action
obstruction. The wildcard is the Clebsch five-frame over GF(16).

Five concepts were kept live: saturated degree budgets; marked joint
neighborhoods; faithful restriction of symmetries; fixed-point residues;
and exact incidence totals versus averages. Their concrete interactions
are the formulas and proofs above.

- **Closest proved mechanism:** the classical saturated-neighborhood
  argument and the repo's exact attachment orbit recurrences.
- **Canonical hostile:** identical labelled cyclic triples with
  automorphism orders81 and27. Recomputed, not replaced by a triangle count.
- **Corrected near miss:** average four was temporarily treated as every
  row or edge having count four. The repaired model retains the variable.
- **Least-used relevant sidecar:** the actual fixed-point character in the
  Singer57/E7 obstruction. It suggests checking which local symmetries
  extend; it does not identify the two56-dimensional objects.

The new local witness and the full action theorem pull in different
directions: the first shows two ports can coexist, while the second
limits how their automorphisms can extend through the third port. This
is why a complete search must retain both adjacency and root positions.

## Reproduction

From the repository root, run:

```text
python -X utf8 -B 05-knowledge/results/ramsey_saturated_ports_20261030/verify.py
python -X utf8 -B 05-knowledge/results/ramsey_saturated_ports_20261030/recurrences.py
python -X utf8 -B 05-knowledge/results/ramsey_saturated_ports_20261030/independent_audit.py
python -X utf8 -B 05-knowledge/results/ramsey_saturated_ports_20261030/seal.py
```

The native verifier and recurrence checker use the standard library.
The separate independent audit uses NetworkX. Generating fresh local
witnesses uses python-sat; checking the saved matrices does not.
Normal and optimized runs must agree. The source manifest pins inherited
Git blobs and exact artifact bytes. This is a same-agent mathematical
audit with separate algorithms, not a Lean build or an external review.

See [SOURCES.md](SOURCES.md) for primary claims and recovered repo routes,
[AUDIT.md](AUDIT.md) for the corrected implication, and
[HANDOFF.md](HANDOFF.md) for the next decisive search.
