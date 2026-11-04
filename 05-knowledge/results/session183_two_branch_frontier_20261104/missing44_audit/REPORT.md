# Independent missing44 audit

**FINITE-EXACT / VERIFIED.** Starting from the actual Session182 digit4 state on the same eighteen names, the Session183 candidate reaches Hamiltonian path count **4,611,841 = 3(1,537,279)+4**. This completes the missing finite word44 return, not reusable two-digit closure.

The independent predecessor subset DP in `count.cpp` counts the actual full tournament, separately from the suffix response formula. The complete count trace is

1537279, 2014907, 1495091, 1559067, 2144877, 2910413, 4610117, 5056377, 5311829, 4216865, 4611841.

`native_audit.py 44` compiles ten edge reversals into twenty literal native path complements, using ten distinct temporary names100 through109 and deleting only those helpers. It checks every directed Hamiltonian path, both complement results, the target rows, and strongness of the source-sized completed states and all three helper-containing states per reversal. Existing names are retained. The schedule and cuts are fully recorded in `native_audit44.json`.

The independent `ports_audit.py` reconstructs the suffix and obtains h=40085, F=633897, J=499420, root E2=5111261, and root defect499420. It verifies the exact ten changed pairs, transitive unchanged prefix50,51,55 and incoming ports1,4,14,21,32,40,54. The suffix changes are exactly {1,4},{14,21},{33,54}. `ports_audit.json` includes all fourteen responses and the full fourteen-by-fourteen ordered joint matrix.

Both auditors pass under normal Python and `python -O` with identical stdout. Conditions use explicit exceptions rather than assertions. No startup helper or prohibited repository file was accessed; no packages installed.

Reproduction from this directory:

- `python native_audit.py 44` (requires existing g++ and sibling missing44 input files).
- `python ports_audit.py`.
- Repeat both with `python -O`.

The native compiler is adapted from Session182's independent auditor; the suffix-response auditor is a separate implementation. These are exact finite certificates, not proof of an indefinitely iterable branching family, spectrum classification, or LRC(14).

## Arbitrary-prefix isolation and search scope

The proof's coefficient identity was independently checked: source intercept14863 and output intercept−459335 give −503924, while633897−3·190302=62991. For a nonempty prefix A, digit4 therefore requires62991(p−8a)=4(1−a). The Hamilton-path cutting injection proves p≥(n+1)a; the equation and a≥1 give p≤8a, hence n≤7. The coprimality gcd(62991,4)=1 forces a≡1 modulo62991, and a≤7!=5040 forces a=1. Any directed triangle has three distinct directed Hamilton paths; sequential vertex insertion preserves each triangle order and produces three distinct full Hamilton paths. Thus a=1 implies transitivity. For a transitive prefix p=2^n, so p=8 forces n=3. This establishes the claimed all-order isolation without assuming a bound on prefix order in advance.

The actual source JSON is exactly equal to Session182 closure/candidate4.json. Lexicographic enumeration independently places triple(1,63,89) at index7569 among the91 choose3 suffix-edge triples. The search stops on that hit; the remaining triples are untested. Its terminal-SCC test and positive monotone clock scan have the stated scope. The inequality J_W≤F_W follows by restricting the second path's first vertex to W in the definition of F_W, so the clock cutoff and pruning are valid. Full-order14 path-cover totals are bounded by15·14!, well within signed64-bit range; port sums and the pruned clock products also fit. No full rerun of all7569 search candidates was needed or performed for this certificate audit.
