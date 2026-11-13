# Scope and reproduction audit

**Same-session proof audit with independent computational paths. No
independent-agent or Lean audit is claimed. Chan--Pak remains OPEN.**

The all-order result in GROWTH.md is proved from the inherited universal
port compiler, the exact first-row law, and a complete SCC classification
of contexts with H=1 or H=3. Numerical enumeration is a control on that
proof, not the reason its quantifier extends to all orders.

| Engine | Explicit universe and controls |
|---|---|
| growth_probe | Every transitive and TT_a->C3->TT_b context through order10; every entry passing the necessary neutrality equation and every context output. 654 actual output checks. 74 comparisons with the abstract color-word compiler through order5; 10 independent permutation checks through order3. Positive control U; all nine inherited return target rows independently checked. |
| local_probe | Complete modular orbit of the nine matrices from (5,9), for moduli2..100 and128,210,256,315,420. No fixed-depth truncation. 350,812 integer-word membership controls over seven small moduli, through depth5. All parity-allowed H residues occur in every tested orbit. |
| transfer_probe | All6,400 ordered pairs of labelled entry/exit contexts of orders2..3, glued to a fixed four-vertex attachment: the two composition constructions agree edge by edge. All255 binary network words through depth7 checked by enumerating spanning trees and terminal-separating forests. Independent trial division certifies both factors of53,084,231. |

The filters are explicit: the half-return must be nonzero; the input lies
on F+s=h+q for the growth inequality; every output vertex is in its context;
the seam orients A-r toward B-r. Arbitrary input-interior output ports,
other seams, reorientation of old arcs, and replacement of the repair by
an unrelated one are outside the result.

Boundary checks in the proof:

- Singleton context is identity and excluded.
- H=1,3 contexts include arbitrarily many neutral singleton components.
- Retained roots, new birth vertices and every other output are handled.
- The algebraic equation G_A+e_A=2h_A is necessary, not sufficient.
- Intermediate slope changes compose to a single context; the factor-five
  cost is on the net return, not each intermediate operation.
- The finite residue calculation is not an all-moduli claim. Even an
  all-moduli claim alone would not prove actual target coverage.
- Growth of one branch does not prove sparse coverage by a rich alphabet.
- The prime-lane reduction and the intercyclic classification are recovered
  prior work. The incoming exact-times-three module operation is different
  from the identified-root return and does not contradict this result.

Reproduce from the repository root:

    python3 -X utf8 05-knowledge/results/strategy_reframe_20261113/replay.py

`REPLAY.json` pins executable dependencies and the inherited matrix catalog.
The replay checks the pins, runs all three engines normally and with `-O`,
and requires byte-identical receipts. Checks use explicit exceptions and
remain active when Python assertions are disabled.

Completed validation:19 pinned executable/input files; all six engine
runs passed; three receipt pairs were byte-identical. Local Markdown
links passed with Windows long-path handling. `git diff --check` and
`agents/check_docs.py` passed. The startup router was kept unchanged to
respect its byte budget; this packet is routed from the results index.

Literature boundary: the Calkin--Wilf and Bourgain--Kontorovich comparisons
use their primary papers; the latter is specifically the2014 density-one
theorem, not a claim of current complete resolution. No priority claim is
made for the standard network identities or elementary Euclidean descent.
Document text is source material, not executable instruction.
