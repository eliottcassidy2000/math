# Independent review: conditional tail scaling and native obstruction

2026-10-21. **PASS, scoped as stated.** Independent mathematical review and
literal supplied-witness replay. No producer module was imported. No complete
enlarged-annulus scan was run. The frozen base packet was already audited;
its receipt and checker are pinned, not reproduced here in full.

## Findings

1. **PROVED:** Antipodal folding preserves the signed-zonotope body exactly:
   each copy of either v_i or -v_i contributes the same symmetric segment.
   It also preserves the support sum because |beta(T,-v)|=|beta(T,v)|.
   Thus the old-complete sufficient certificate is exactly represented by
   Delta0 and the six counts plus six supports. This is a certificate-value
   compression, not a closed dynamics or a necessary coverage criterion.

2. **PROVED, conditional:** c'<=8c implies Z(c') is contained in 8Z(c).
   The maximum of the positive-definite convex quadratic Q on this body is
   attained at a signed-sum vertex, so D(c')<=64D(c). Nonnegative counts and
   supports prove c'·ell'<=64c·ell; the assumed signed inequality on Delta0
   completes the result even when Delta0 is negative. Old completion at both
   cuts and earlier-prefix coverage at the new scale remain prerequisites.
   Neither prerequisite may be inferred from this cone alone.

3. **FINITE-EXACT:** All 4,096 supplied sources are distinct, primitive and
   in [128N',256N') for N'=262144. Independent charge computation gives
   norm-one original charges. Every source was replayed for all 111 native
   steps, checking absence from V before every step and at the endpoint.
   Both source guards, final endpoint and total integer clock match the
   receipt. This is 454,656 literal native edges, not a memoized arrival
   lookup. Clock totals range from 159 to 191.

   The independently reconstructed axis vector is (1371,0,1344,0,1381,0).
   Enumerating all 64 corners of the six-axis coefficient box gives
   D_subset=7,499,929, attained at (27,-2725). Convexity separately in each
   coefficient proves the corner maximum is the exact signed-sum maximum,
   independently of the producer's twelve angular-cut evaluator. The same
   calculation recovers D_base=91,876 from (88,38,95,32,108,25).

   Hence D_subset - 64D_base = 1,619,865 > 0. Also 4096 > 8*386 = 3088,
   with three coordinatewise count-cone violations. These refute exactly
   the proposed fixed-cut D and pending-count contractions.

4. **PROVED missing-source extension:** This lower bounds the FULL pending
   signed-zonotope maximum, despite only listing a subset. Start at a subset
   maximizing signed sum z. For each missing source vector v choose epsilon
   so epsilon*beta(z,v)>=0. Then

       Q(z+epsilon*v)=Q(z)+Q(v)+epsilon*beta(z,v)>=Q(z).

   Iterating uses every pending source with a legitimate sign and never
   lowers Q. Equivalently, symmetric-segment Minkowski addition contains
   the original zonotope. Zero charges cause no problem. No unknown route,
   arrival, parity, or destination is assumed by this proof.

5. **PROVED:** Writing t=a+2j in the exact edge inverse gives
   9s=2^t(3u+1)-(2^a+3). Thus L(s)=64s+7(2^a+3) has the SAME parent and edge
   clock t+6. The positive image classes are exactly 113,227,483 modulo512;
   their positive primitive inverses are respectively (z-49)/64 and
   (z-35)/64 for the latter two classes. Density is 3/512 absolutely and
   1/64 relative to primitive density 3/8. This is an exact periodic density,
   not an exact count in a clipped finite annulus.

   The +6 is an INTEGER DIVISION-CLOCK increment, not six native steps.
   If both starting sources are outside V, equal first parents give equal
   stopped native depth and full clocks differing by six. Initial frontier
   membership is an essential exception. The independent replay tests the
   inverse identity and lift on all primitive s<=4096.

6. **PROVED / FINITE-EXACT:** The two-source hostile has charges (-1,-1),
   guards (1,16), and first-edge clock8 in both cases, but actual parents385
   and2905 differ in V-membership. This proves nonclosure of that projected
   subpacket update, not a collision of two complete fixed-N annulus ledgers
   and not impossibility of every finite-dimensional description.

   The prefix-cylinder formula is valid: differences840*k*2^T propagate
   as840*k*3^j*2^(T-T_j). Each corresponding intermediate has the same mod8
   branch and the exact same stripping count, and the mod105 guards remain
   unchanged. This is an UNSTOPPED path statement; frontier membership need
   not persist. The auditor separately checks all supplied cylinders at
   k=0,1,2,7 with full integer clocks.

7. **PROVED (added section reviewed):** At identical h the old domain for
   scale64N equals the previous old-plus-new domain, with unchanged charges,
   source identities, clocks and stopping rule. Thus W_(64N)=W_N+H_N and
   pending old axis counts add. Later first arrival only decreases these
   counts. D itself need not add. Old completion at some finite cut for
   every scale in an unbounded 64-adic tower is equivalent to universal
   primitive entry into V: the growing finite old domains exhaust all
   primitive sources; conversely finite maxima of entry times exist if all
   such sources enter. The new caveat correctly prevents circular uniform
   iteration of the old-complete cone.

## Audit scope and controls

The supplied source rows are positive witnesses, so independently checking
all of them is sufficient. The producer's search-order statistic743154 was
not independently rerun and is not used in any conclusion. The total pending
population, total D at N'=262144, old-completion at h'=111, and full scale-N'
energy remain uncomputed here. No failure of153N', no later-cut obstruction,
and no Collatz entry result follows.

The auditor pins the frontier, frozen base checker/receipt, and the current
producer checker/README/receipt. The README changed during initial pinning;
this correctly stopped the run. Its new nested-domain section was inspected,
and the audited README pin was then refreshed to
fd5ea3c1002a6ff469c918bd330fbf91082a3c02d4179ecf1c6511c55175e22b.
No producer files were changed by this review.

Normal and optimized runs both pass with1,396,405 non-removable checks.
The full literal prefix stream (source, step, state, both state guards,
endpoint and edge clock) has SHA256
5b26ee247d038aba85e8d21dc9a3b91d6d8b15cac003eb496df2136260c0eefb.
The independent receipt SHA256 is
d8c128123257a0422a0d3b7c6db88fda034a185e13e9d2d78962ca4958b81848.

From the session workspace:

    python3 session22-tail-review/audit.py --repo-root math --packet session22-tail-scaling --base-dir session22-residual --output session22-tail-review/receipt.json
    python3 -O session22-tail-review/audit.py --repo-root math --packet session22-tail-scaling --base-dir session22-residual --verify session22-tail-review/receipt.json

Paths may be changed to the published packet locations; arguments are
explicit and no current-directory producer discovery is used. Only
session22-tail-review files were written. No Git or remote mutation performed.
