# Independent audit of the Session181 multiport +4 witness

The candidate passes an independently written endpoint subset dynamic program
(`audit.py`), including direct full-tournament counts for source and target,
all named incidences, five changed pairs, strong source/target structure,
strong suffix, physical inports, joint suffix field and root two-slot defect.
All checks remain active with Python `-O`.

The source snapshot `source_interface.json` is copied unchanged from Session180
`returns/source_interface.json`. `audit.py` reads this local snapshot and
`../multiport/candidate.json`; it is portable within this result packet.

The native compiler `native_audit.py` adapts the independently verified Session180
compiler. It finds explicit directed Hamiltonian guards, replays two native
complements per edge reversal on a fresh temporary helper, deletes only that
helper, and checks the literal target adjacency. Its separately compiled C++
subset counter checks every completed-state count. Normal and optimized Python
runs gave identical results:

170807, 217551, 1215017, 1412269, 1294011, 512425.

All six completed states are strong. The schedule uses ten complements and
five temporary helpers, preserving all seventeen inherited names. Full guards,
cuts, and intermediate rows are in `native_audit.json`.

The constructive worker's general identity checker was additionally rerun under
`python -O`; `formula_controls.json` records its 3069 passing cases. This audit
does not claim an independent complete rerun of the larger finite search menu.

The established target has H=512425, root E2=516115 and defect3690. Thus it is a
literal legal native +4 transformation, not a zero-ready return and not a
recursive closure theorem. Neither scalar equality nor the finite search proves
an infinite invariant family.
