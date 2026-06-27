# S239 Reflection: Renormalized Polymer / Dirichlet Bridge

The useful surprise in this pass was the namespace collision.  I first claimed
HYP-3072 locally for this polymer/Dirichlet bridge, then the checkpoint rebase
revealed incoming S238 had already taken HYP-3072 for the cross-carrier
pullback portfolio.  Renumbering to HYP-3073 was not just bookkeeping: it made
the structure clearer.  HYP-3072 is the portfolio of possible carriers;
HYP-3073 tests two old carriers for theorem-facing payload.

The signed-polymer lesson is that the old failure was not simply "polymer gas
bad."  It was "untyped absolute activity bad."  HYP-2645's relation-packet
figures show why: raw R6 density does not predict signed correction.  `odd_AP`
has fewer relations than `near_AP` but larger correction, while wide/Sidon
packets have tiny signed activity.  A plausible next theorem would type
activities by packet class, isolate AP-like slow sectors, and route wide/Sidon
or repeated-residue sectors through finite-cell or character channels.

The Dirichlet lesson is complementary.  Sidecars should be treated as boundary
conditions, not as labels pasted onto a route after the fact.  The toy network
is intentionally modest, but it already separates raw route conductance `1/2`
from legal sidecar conductance `9`; the phantom F7 class remains a one-unit
boundary exit rather than an anonymous hard bucket.

Tournament Analysis used proof carriers and renormalization/energy obligations,
not runners, raw routes, or median centers.  The challenged assumption is that
a final LRC proof route must be judged at the level of route labels.  This pass
suggests an alternate quotient predicate: a sidecar-completed residual cannot
cover or discharge illegally if either its typed signed activity is controlled
or its boundary current survives Schur complementation to a named exit.

Next work: build the actual HYP-2963 typed-polymer ledger and residual sidecar
graph.  The proof-facing outputs should be `signed_polymer_packet_type`,
`signed_activity_budget`, `finite_cell_route`,
`schur_complement_conductance`, `sidecar_energy_exit`, and
`phantom_f7_boundary_atom`.

