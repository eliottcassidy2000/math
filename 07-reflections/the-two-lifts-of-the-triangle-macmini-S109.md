# The two lifts of the triangle — and where parity laws live

*mac-mini-2026-07-15-S109. Session artifacts: THM-865 (refutation), THM-866 (proof),
figurate-two-axes draft.*

The triangular numbers sit at the corner (s, d) = (3, 2) of the figurate array
N(s,d,m) = (s−2)C(m+d−2,d) + C(m+d−2,d−1), and the two ways out of the corner are
different in KIND, not just direction. The shape axis is **affine** — one more side glues
one more triangle (P(s+1,m) = P(s,m) + T_{m−1}; n² = T_{n−1} + T_n is its first step, the
Mode-A staircase gluing). The dimension axis is **integral** — one more dimension is a
partial sum (coning). Affine steps preserve degree; integral steps raise it. That single
distinction generates everything the owner noticed: polygonal columns stay quadratic,
Pascal columns grow; rows sum to Moser (a degree-4 truncation) vs 2^r; diagonals to G(n)
(a degree-4 quasi-polynomial) vs Fibonacci.

The bridge between the two lifts is a *filtration, not a formula*: grade the simplex
count by support (= runs, = open faces of the dilated simplex), and the polygonal world
is the 1-skeleton — the part of every dimension that is still made of edges, i.e., still
made of glued triangles' boundaries. The deficit scheme is then Pascal itself with its
first three columns deleted: **the difference between the two generalizations of the
triangle is the part of Pascal beyond what the generalizations share.** Cutting the shared
frame off Pascal leaves Pascal. The triangle keeps pointing at itself.

Two completeness proofs landed today with the same skeleton and opposite fates, and the
contrast is the real lesson. The axis level lattice has NO HOLES because a conserved exact
step exists: F3's tie-splitting flip always advances x by exactly +8, and the only
stuck configuration is the transitive order (THM-866 — greedy walk + conserved step =
interval, the same shape as Brown's criterion for complete sequences, Zeckendorf's
exchange, and the binary numeration at the tight corner). The locker tournament's H mod 4
has NO LAW because no conserved pairing exists: the divisor-pairing involution that tames
τ(n) refuses to act on odd cycles — its exit-cell shadow looks real through m = 18, then
m = 20 breaks the compensation and primes never had it (THM-865). Parity statements
survive exactly as far as an explicit involution or an explicit conserved step carries
them; where the mechanism ends, the law ends — Rédei's digit is constant because the
reversal involution is total; digit 1 decays (THM-466(iv), opus-S316 probe (ii)) because
nothing total replaces it, and D_n is now the cleanest witness: a single canonical
arithmetic family where digit 1 is already lawless.

Smallest good fact of the day: Moser's missing 32nd region and G's missing 34th
(33 = 34 − 1 at n = 8) are the SAME object one filtration level up — the unique minimal
3-run pattern, cell (m,k) = (3,3), deficit exactly 1, appearing at r = 2j+1 in rows and
n = 3j+2 in diagonals. The famous circle-problem "hole" is not an accident of geometry;
it is the first interior point of the first 2-face of the dilated simplex.
