# Complete dyadic guards for shared-source rectangle circuits

**PROVED + FINITE-EXACT, 23 October 2026.** Read `THEOREM.md` for proof,
source scope, realization hypotheses and the preserved/lost-coordinate ledger.

Main increment: paired perturbations of the original 3,528-node native
atlas need only the complete marked label `(Y0,Y1,Y2,Y3 mod2; X mod128)`.
Its group is `(Z/2)^4 + Z/128`; the 1,695 occupied vertex labels are in
`node_labels.json`. Exact rowwise neutrality retains EVERY original dyadic
harmonic bank, not just the saved primitive modulus-256 vector. Odd primary
characters and native resources remain separate obligations.

The compiler retains every physical shared row. Four parity-cut tests and
seven integer carry layers diagnose the precise obstruction. Two individually
unsafe rectangles can cancel at shared target rows through a balanced
multi-source charge. Explicit path resources, full shared records and a
reserved unchanged complement remain necessary for native injectivity.

## Reproduce

    python3 -B check.py --repo /path/to/math --out normal.json
    python3 -O -B check.py --repo /path/to/math --out optimized.json
    cmp normal.json optimized.json

Both runs pass **60,205 checks** and have identical output. The generated
finite-ring controls include overlapping source/target roles, negative
amplitudes, cut parity, nonzero carries, the old fixed-source tariff, and the
weighted-C4 scalar-phase/current hostile.

    python3 -B recheck_filtration.py --repo /path/to/math
    python3 -O -B recheck_filtration.py --repo /path/to/math

Both reruns match the entire pinned original order-eight filtration, with
696 original implementation checks. This is a rerun, not an independent
implementation of that filtration. The absent historical order<=7 input is
not silently replaced or claimed rerun.

## Compile a supplied menu

    python3 -B compile_menu.py --menu algebraic_circuit_menu.json

The example uses original marked vertices and cancels sources `{0,55}`
against `{1,46}` at common targets398/528. It is explicitly an algebraic
certificate, NOT an assertion of available native paths. Each component
rectangle fails; their combined owned row debt vanishes. The normal and
optimized compiler receipts agree exactly.

Optional resource-ledger and shared-record sections check only the supplied
keys for capacity/conflicts. The compiler does not authenticate their native
meaning. The separate native producer must supply actual matrices, paths,
owners, kernel lifts, converse squares and disjoint F/R resource orbits.

`SOURCE_PINS.json` pins all original mathematical inputs. `MANIFEST.json`
and `SHA256SUMS` freeze this packet. No files in the repository checkout were
modified by this lane, and no Git actions were performed.
