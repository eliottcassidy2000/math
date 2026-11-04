# Independent audit of the incoming integration note

**PASS: source fidelity, scope, arithmetic, pinned links, and independent small literal controls.** No mathematical correction is required. This audits NOTE.md, not a fresh execution of the incoming native compiler or its complete receipts.

## Sources checked independently

The GitHub connector fetched the three primary files directly at commit `cc49aaf93bb89592e11fbcaf2ce86385ff1d302c` under `05-knowledge/results/multiport_grammar_20261103/component/`:

- `NATIVE_SIX.md`, Git blob `57c83e1940ec9ab92db06c6b13dae48d5d8709de`;
- `PREPARED_FOUR.md`, Git blob `0a8249e316f4a2bcd6f82cb82be95b0c18e0c746`;
- `OBSERVER_HOSTILE.md`, Git blob `037c8d4a168422df7ad880d613ebe78a2a9f103a`.

The commit parent was independently checked as `3ec7c475dbdcda7c6f3224813d93ccbdceeb3b10`, confirming the one-commit-after-S186 statement. Each local copy agrees with its pinned source except for one added terminal blank line. The two linked audit sidecars were also fetched at the same commit: `NATIVE_SIX_AUDIT.md` (blob `fccc08d6c1a24cca30b2f46d2a97bfb804bdf103`) and `PREPARED_FOUR_AUDIT.md` (blob `594405248a6c9fc3d1336f60fb237ed7e96d12ff`). Their independent analytic/source/replay PASS attributions are accurate.

## Accepted claims and boundaries

1. The native repair is genuinely all-clock t>=1, preserves all old names/root, and uses only temporary helpers, legal native complements, and deletion of those helpers. The source provides a structural bulk strongness argument and distinct finite patches at t=1,2. The note correctly avoids extending the unchanged bulk graph to an empty prefix: its t=2 unpatched count would be219, not223.
2. The displayed seven-coordinate hostile and both post-flip signatures match the pinned source. The common named action is available in both graphs. The statement concerns an assembly signature, not equality of complete named two-port response functions.
3. The terminal preparation applies to the actual bulk output, hence t>=3. Its (sigma,tau)=(2,0) gives exactly (H,F)->(3H+4,3F). The note does not claim it for the t=1,2 patches. The nine-vertex475 and eleven-vertex475 distinction is accurately attributed to the source; those graphs are not identified by H.
4. The all-auxiliary obstruction is restricted to subsequent C2 digits2/4 in the C2 image. The transported auxiliary supports indefinite zero steps. The note correctly avoids a reusable nonzero or two-branch claim, and does not rule out other recipes or genuine re-preparation.
5. The connection to the session's rank theorem is limited to preserving relevant named information. It does not transfer reconstruction between the two different families or claim closure of a fixed two-port observer under arbitrary-depth operations.

## Independent exact controls

`audit_note.py` is a new standalone subset-endpoint checker, without imports from the incoming producer. It reconstructs the physical clock-one hostile attachments and checks423->611 versus423->471, the corresponding root defects112->114/84, and strongness. It separately reconstructs the actual clock-three source, bulk target, and two permanent C2 births, checking

    H: 157 -> 475 -> 1429,
    F at the retained root: 169 -> 521 -> 1563.

The old induced graph is unchanged under both births. These are small literal controls; the all-clock statements remain attributed to the source proofs and their existing audits. Normal and optimized outputs agree byte for byte in `audit_note.json` and `audit_note.optimized.json`.

Optional clarity only: writing “actual bulk output (t>=3)” in NOTE.md would make its already-correct domain explicit at the point of use. No status promotion of the incoming packet depends on this audit.
