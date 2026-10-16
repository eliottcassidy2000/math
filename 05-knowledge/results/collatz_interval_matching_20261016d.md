# Interval certificates for native Collatz cancellation

2026-10-16. **PROVED:** every supplied finite native controller admits an
optimal source-atom matching for all annular rows in an octave, encoded by
fewer than `16X/3` interval records. Valid records survive common native
extension and every physical-bin refinement. A subpower bound for their
coverage deficit suffices for the dyadic square-root energy hypothesis.
**FINITE-EXACT:** the declared five-scale, eight-reset controls.
**OPEN:** the uniform coverage-deficit bound, paired half-energy, independent
ROOT positivity, and universal convergence. Self-audited; no independent
referee or historical priority claim is made.

## 1. Inheritance and the change of representation

The closest mechanism is the actual target-fibre matching in
[native survivor packets](collatz_survivor_packets_20261015e.md).
[Refinement debt](collatz_refinement_debt_20261016a.md) already supplies a
nonaccumulating Gram budget; it does not pay the coarse mean.
The one-reset `7X` failure at `k=23` and the angular-debt correction in
[paired energy](collatz_paired_energy_20261015g.md) remain in force.
The incoming [marked elliptic lift](collatz_curve_bridge_20261016c.md)
preserves a supplied word, but identifies cancellation only after projection
to the shared integer endpoint and its weighted clock sign.

The unused coordinate exploited here is the *interval of physical rows in
which an original source occurs*. A fixed pairing of sources can miss
cancellation because different rows need different partners. We retain
each source and allow only the matching to switch at authenticated integer
boundaries. This gives a certificate format with universal finite coverage;
its quantitative arithmetic score is still unpaid.

The concept board is original source atoms / literal path clocks / common
terminal / interval capacity / refinement / the cubic geometric seam.
The companion [cubic packet](cubic_crossedge_packets_20261016d.md) keeps
the geometric ports just as literally. The connection is preservation of
the next operation's interface, not an identification of unit walks with
directed Collatz words. The inherited META-PATTERNS cards on shared labels
and operation closure apply; no new general card is claimed.

## 2. Source atoms and their exact physical intervals

Use the shortcut map `T(n)=n/2` on evens and `(3n+1)/2` on odds, stopped
at ROOT 1. Supply a finite common observed graph whose paths from the
original sources terminate. Write its authenticated endpoint and literal
clock as `e(n),d(n)`. Every path edge is checked against T. A reset is not
renamed as a single T step.

Fix an integer `X>=1`. For integer `X<=N<2X`, set

\[
 v_N=\sum_{N<n\le4N} w(n)(-1)^{d(n)}[e(n)],
 \qquad w(n)=3[n\equiv2\pmod3]-1.
 \tag{1}
\]

Replace source n by `|w(n)|` distinct unit atoms `(n,c)`, with
`0<=c<|w(n)|`. Each atom retains n, c, e(n), d(n), and sign
`sign(w(n))(-1)^d(n)`. Its set of physical rows is the integer interval

\[
 I_n=[\max(X,\lceil n/4\rceil),\min(2X,n))\cap\mathbb Z.
 \tag{2}
\]

An interval record consists of two distinct source atoms, an interval
`[l,h)` contained in both sets (2), and their supplied native paths.
Their endpoints must agree and their **weighted** signs must be opposite.
For each atom, the intervals on all its records must be disjoint. These
are capacity conditions at each physical row, not a license to reuse a
source's mass twice in the same row.

Let `C(N)` count records covering N. There are exactly `4N` active atoms:
the `3N` consecutive sources contain N weights of absolute value two and
2N weights of absolute value one. Cancelling each recorded pair gives

\[
 \|v_N\|_1\le B(N):=4N-2C(N).
 \tag{3}
\]

This implication needs no optimality claim, ROOT oracle, eventual sign,
or density model. A proposed family of intervals can be checked by its
actual paths and capacity events alone.

## 3. An optimal certificate with a linear number of records

**I1 (PROVED).** For every supplied controller above there is an interval
certificate with at most

\[
 P_X=4X-2+\sum_{n=X+1}^{2X-1}|w(n)|<\frac{16X}{3}
 \tag{4a}
\]

records and

\[
                 B(N)=\|v_N\|_1\quad(X\le N<2X).
 \tag{4}
\]

**Proof.** At N=X, match opposite atoms separately in each common terminal
fibre until only one sign remains in each fibre. At most `2X` pairs are
created. Move from N to N+1 by removing source N+1 and adding the four
sources `4N+1,...,4N+4`. When an atom leaves, close its current record;
if it was paired, return its partner to the unmatched pool. Keep all
other pairs. Match the available opposite atoms in each fibre again.

Count records at their other end. Every record closed before 2X is closed
by one departing atom, which can be charged only once. There are
`sum_(X<n<2X)|w(n)|` such atoms. At the final row, at most `4X-2` pairs
remain, since its raw mass is `4(2X-1)`. This proves (4a). Explicitly,
P_X is `(16X-9)/3`, `(16X-10)/3`, or `(16X-11)/3` as X is 0, 1, or 2
modulo three. Keeping the final norm also saves `||v_(2X-1)||_1/2` from
this bound.

On every row, every fibre's unmatched atoms have one sign. Their total is
exactly that fibre coefficient's absolute value, proving (4). Closing
the remaining records at 2X completes the certificate. QED.

This budget is sharp for unrestricted terminal assignments with the
same source weights. Send every source to one formal terminal with even
clock. Every row is then zero because w has period three and zero sum.
Every departing atom must have been paired on the preceding row, so at
least one distinct record ends for each such atom; two atoms departing
together belong to the same source and have the same sign, so cannot
be paired together. The last row needs exactly `4X-2` further records.
Thus every exact interval certificate needs at least P_X records. This
formal assignment is **not a native Collatz controller**; it proves
sharpness of the universal representation bound, not arithmetic sharpness.

The bound counts interval records. It does **not** bound native path
lengths, endpoint bit lengths, the size of the supplied controller, or
the time to discover a controller with a small score. In the implementation,
lazy priority queues select the same least-labelled partners as a direct
ordered-set implementation; queue operations do not change the mathematics.

### Why a single fixed matching is too restrictive

Take the restricted subledger consisting of the three original sources
`3,4,9`, each with its original weight -1. Their actual ROOT clocks are
`5,2,13`, respectively, giving signs `+,-,+` at endpoint 1.
For N=2, its active atoms are 3 and 4; for N=3 they are 4 and 9. Both
restricted row norms are zero. Attaining both with a fixed matching would
use source 4's unit capacity once against 3 and once against 9.
Even fractional fixed capacities cannot do this: each row needs mass one.

Interval records `(3,4,[2,3))` and `(4,9,[3,4))` repair the obstruction
exactly. This is a hostile in a specified source subledger, not a claim
that those are the complete annular rows. The missing coordinate was
the physical row in which capacity is consumed.

The sharper hostile to overreading I1 is already known: the halving-only
controller has a linear certificate but energy asymptotic to
`(700/27)X^3`; see [odd-ray payment](collatz_odd_ray_overlap_20261015d.md),
equation (12). A small certificate format does not imply a small score.

## 4. Source-preserving refinement and changing terminals

**I2 (PROVED).** Every valid interval certificate remains valid after a
common native extension of its terminals. If e is sent to f in exactly
tau further steps, both atoms of a pair acquire the same added clock tau.
Their endpoints still agree, and their weighted signs remain opposite.
The source capacities and physical intervals do not change.

The full integer clocks must still be retained and replayed. Replacing a
clock by `d+2` leaves all sign tests unchanged but generally falsifies the
path receipt; the checker includes this hostile.

For any physical bin J, a record supplies the integer pair capacity
`|J intersection [l,h)|`. Each source's total used capacity is at most
`|J intersection I_n|`, exactly its multiplicity in the sum of the original
rows. Splitting J simply intersects the same intervals with the child
bins. Thus every refinement retains the source coefficients and yields

\[
 \frac{\|\sum_{N\in J}v_N\|_1^2}{|J|}
 \le \frac{(\sum_{N\in J}B(N))^2}{|J|}
 \le \sum_{N\in J}B(N)^2.
 \tag{5}
\]

Coarse rows can have additional cancellation across physical rows, so no
equality with their optimum is asserted. The certificate already bounds
the finest rows and therefore survives every coarsening or refinement.
An extension may also add new matching intervals. Existing intervals
stay valid; a fresh optimal sweep may reroute them while preserving
every original atom and its actual path.

## 5. The unpaid target as an integer coverage deficit

The sharp point-row target is now the following purely positive score:

\[
             \sum_{N=X}^{2X-1}(4N-2C(N))^2\le49X^2.
 \tag{6}
\]

A construction of valid intervals satisfying (6) suffices. It need not
prove that the intervals give the optimal matching, or evaluate eventual
ROOT membership. Formula (6) therefore admits direct arithmetic packing
proofs, with literal source-path receipts as the admissibility condition.
Existence of sufficiently dense such packings remains **OPEN**.

For optimal certificates at two native stages let `U(N)=B_old(N)` and
`D(N)=C_new(N)-C_old(N)>=0`. The later norm is `U-2D` and the exact energy
drop is

\[
                E_{old}-E_{new}=4\sum_N D(N)(U(N)-D(N)).
 \tag{7}
\]

When `E_old>49X^2`, paired half-energy or entry below the floor is therefore
equivalent to

\[
 \sum_N D(N)(U(N)-D(N))
 \ge\tfrac14\min\{E_{old}/2,E_{old}-49X^2\}.
 \tag{8}
\]

This is the earlier weighted-overlap invoice in an explicit interval
representation, not a new proof of its lower bound. Counts of records or
unweighted terminal collisions do not establish (8).

### A weaker useful target recovered from the incoming work

The incoming pre-audit
[THM-4660, conditional Riesz formula](../../01-canon/theorems/THM-4660-collatz-conditional-explicit-formula-for-riesz-means-of-the-parity-sum.md)
uses the weaker hypothesis (W): dyadic discrepancy energy
`E_X=O_epsilon(X^(2+epsilon))` for every positive epsilon. Its hypothesis,
rather than its unaudited contour conclusion, motivates this reroute.

**I3 (PROVED implication; arithmetic premise OPEN).** Suppose that for
every epsilon>0 and every sufficiently large dyadic X one can supply a
finite native controller and valid interval records with

\[
                 \sum_N B(N)^2\le C_\epsilon X^{2+\epsilon}.
 \tag{9}
\]

Then (W) holds for the actual ROOT-supported parity function.
Indeed its terminal values have absolute value at most one, so (3)
bounds `D_X=sum_N(H(4N)-H(N))^2` by (9). Apply
[THM-4658, two-octave energy descent](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md),
statement 3:

\[
 \sqrt{E_{4X}}\le2\sqrt{E_X}+2\sqrt{D_X}+4\sqrt{3X}.
\]

After division by `(4X)^(1+epsilon/2)`, the inherited-energy coefficient
is `2^(-1-epsilon)<1`; the forcing is bounded. Iteration along each of
the two four-adic chains proves `E_X=O_epsilon(X^(2+epsilon))`.

Equivalently, a bound `X^2 B_k` at `X=2^k` with
`log B_k=o(k)` is enough. A fixed constant seven and a linear reset-stage
deadline are stronger optional targets. No quantitative controller bound
has been proved here, and (W) alone neither proves a nonzero residue nor
universal ROOT entry.

## 6. Exact controls and next arithmetic question

[Program](../../04-computation/experiments/collatz_interval_matching_20261016d.py)
and [receipt](collatz_interval_matching_20261016d.json).
The universe is `k=6,8,10,12,14`, every source through `8*2^k`, and exactly
eight first-reset/closure stages, with the initial halving stage retained.
The common controller is inherited unchanged from the compact-reset code.
Every original path is replayed with its full integer clock. A separate
terminal-coefficient event sweep verifies every reported row optimum;
small cases also recompute the full vector norm at every row.

| k | Records after stage 8 | Point energy divided by X squared |
|---:|---:|---:|
| 6 | 334 | 0.99511719 |
| 8 | 1350 | 0.72271729 |
| 10 | 5422 | 2.01906204 |
| 12 | 21690 | 7.88526630 |
| 14 | 87100 | 4.49702376 |

All 35 paired comparisons pass the inherited half-energy-or-floor test.
These are exact tests of a new representation at previously explored
source scales, not an all-scale inference or a held-out-scale forecast.
The receipt retains exact integer energies; decimals in the table are
display values. It stores all intervals, endpoints, clocks and row norms
for k<=8, and deterministic hashes for the larger reconstructed objects.
Malformed source copies, duplicated mass, out-of-range intervals and a
same-parity false clock are rejected. The unrestricted sharp-format
control attains (4a) for every integer X from 1 through 64; its labels
are explicitly excluded from the native-controller tests.

Normal and optimized runs agree byte for byte, with 4,985,024 explicit
checks in this script in addition to the imported native-path audits.
The receipt SHA256 is
`403e10064a88d132cd5405502bc9743d619d218b5445c86a7c9e976061c7e8e7`.

```text
python3 -B 04-computation/experiments/collatz_interval_matching_20261016d.py --k 6 8 10 12 14 --stages 8
python3 -B -O 04-computation/experiments/collatz_interval_matching_20261016d.py --k 6 8 10 12 14 --stages 8
```

The next mathematical obligation is to select explicit families of native
source pairs and prove their *interval coverage deficit*, allowing partner
switches where the three-source hostile demands them. A subpower deficit
bound (9) would already reach a meaningful analytic target. The linear
certificate theorem makes that target concrete without assuming it.
