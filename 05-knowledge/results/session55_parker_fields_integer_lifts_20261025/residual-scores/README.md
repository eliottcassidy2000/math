# Residual-score proof and independent test packet

Read THEOREM.md for R1–R4 and the odd-order supplement. The packet separates
the true even-cycle half rule from the path's extra residual edge, proves a
necessary-and-sufficient cut criterion, and gives all-size prescribed-high-set
existence theorems with their exact small exceptions. It retains native path
endpoints and recovers the Session52 odd-order Eulerian lift.

Reproduce with Python 3, no packages:

    python3 check.py > normal.json
    python3 -O check.py > optimized.json
    cmp normal.json optimized.json
    python3 verify_manifest.py

The checker exhausts all labelled residual orientations for the FIXED path
and FIXED cycle at orders 3 through 6. Independently it enumerates every
nonnegative candidate score vector bounded by the graph's degrees and with
the exact total. It compares actual support, all subset inequalities, and a
constructive bipartite-matching algorithm, rather than trusting just one route.
It constructs every prescribed high-set pattern at even orders 4 through 14,
and every odd-path single-high pattern and odd-cycle Eulerian prescription
at orders 3 through 15. Those larger-order tests are exhaustive over the
stated score-pattern universes, NOT all orientations or isomorphism classes.
At orders 4 and 6 it exhausts the fixed-path cube, checks endpoint criteria,
the internal high-set action, the regularizing extension, and its Eulerian
cycle residual. Both small graph exceptions and the incorrect path half
count are explicit hostile controls. Checks use exceptions, not asserts.

The pinned Session52 source is copied solely for the inherited endpoint and
regular-cycle bridge, with no dependence on its other enumeration or topology
claims. Manifest checks cover every packet input and output. The checker is
self-contained and does not require the source checkout or a network.

Scope boundaries: score support is not coefficient multiplicity, quotient
class count, or wiring; no map to square-magic arithmetic has been proved.
