# Independent row controllers and literal interval coverage

2026-10-16, clock_duality_oct15. **PROVED, scoped:** row transport,
independent-controller gluing, the sharp cone envelope, and sufficient
sampled-energy bounds. **FINITE-EXACT:** the declared finite controls below.
**OPEN:** the uniform sampled bound, (W), independent ROOT positivity,
and universal convergence. Self-audited; no independent referee or priority
claim is made.

## 1. Inheritance and the new degree of freedom

The closest mechanism is I1--I3 of
[optimal interval matching](collatz_interval_matching_20261016d.md).
Keep its sharp `P_X<16X/3` record bound; this note does not replace it with
a weaker count. The hostile is the halving controller, which has linear-size
certificates but cubic energy. The corrected near miss is the second odd
block's loss of a common sign cone in
[odd-ray overlap](collatz_odd_ray_overlap_20261015d.md), not a failure of its
proved first-block estimate. The least-used coordinate here is the location
of a physical row where a particular controller is useful.

The anchor is independent interval-deficit bounds, the niche is the cubic
unit-path ring, and the wildcard is whether local arithmetic witnesses can
be glued without losing their marked source. The board is original atoms /
integer clocks / independently chosen terminals / physical intervals /
prime-ideal branches / marked elliptic torsion. The companion
[cubic adelic note](cubic_adelic_20261016e.md) gives an actual local-to-global
path theorem; it does not transfer that theorem to directed Collatz entry.
META-PATTERNS' controlled-forgetting and shared-label cards apply.

The incoming [THM-4661, odd-degree exteriors](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md)
was read with its nineteen-vertex refinement and pre-audit status. Its
list-criticality argument changes the geometric search floor, not a Collatz
estimate. The older
[adelic price sheet](collatz_pricing_approaches_20260927.md) already proves
an exact real/2-adic/3-adic product formula and documents regenerated
precision. Rephrasing that identity is not a new descent mechanism.
The final pull through `4971cf2b05` adds the finite-exact 22-vertex screen
and THM-4662's forbidden hub. The companion B6 completes its suggested
everywhere-local/global-failure example; B7 changes the carrier to realize
the hub, and B8 exposes a coefficient-field obstruction to a naive cubic
lift. This reinforces the need for A1's literal source transport: local
acceptance by itself does not prove a global compatible certificate.

## 2. A source-preserving transport that costs at most four

Use the shortcut map T, stopped at 1, and

\[
 w(n)=3[n\equiv2\pmod3]-1,\qquad
 A_N=\{(n,c):N<n\le4N,\ 0\le c<|w(n)|\}.
\]

There are exactly `4N` original unit atoms in A_N. A point certificate at N
pairs distinct-source atoms with the same literal endpoint and opposite
weighted clock sign. Each atom's path is a finite actual T path, ends at
an odd integer, and includes its full integer length. An atom is used at
most once. The residual is `B=4N-2(number of pairs)`. No ROOT membership
oracle is part of this definition. If desired, extend every even terminal
by actual halvings to obtain the odd-terminal convention.

**A1 (PROVED).** A point certificate at M transports to every N>=1 with
residual at most `B+4|N-M|`, preserving the actual sources at N, their
weights, and authenticated paths.

**Proof, forward.** In the move N to N+1 put m=N+1. Replace each departing
atom `(m,c)` by `(4m,c)`, retaining its partner and endpoint and prepending
exactly two halvings to its path. The weights agree since `w(4m)=w(m)`,
and the two extra steps preserve the clock sign. Add the atoms of
`4m-3,4m-2,4m-1` unpaired, with their initial halving paths. Their total
absolute weight is four. No surviving atom is duplicated.

**Proof, backward.** In the move N to N-1 put m=N. Delete the atoms of
`4m-3,4m-2,4m-1`, releasing any surviving partners. This increases the
unmatched count by at most their total mass four. Replace each atom of
4m by the corresponding atom of m and strip its first two genuine halving
steps. Its odd endpoint ensures those two steps were present. The weight
and clock sign agree. Iteration proves the assertion. QED.

This is transport of a certificate between different physical source
annuli. It never treats an inverse Collatz choice as a forward orbit step.
The new atom 4m has its own literal source label and two additional native
steps; the old label m is not silently reused.

## 3. Independent controllers glue without identifying their terminals

Let t(n) be the ROOT-supported parity function: `(-1)^d` if n first
reaches 1 after d shortcut steps, and zero if it never does. This definition
does not assume universal convergence. Put

\[
 H(M)=\sum_{n\le M}w(n)t(n),\qquad a_N=H(4N)-H(N).
\]

Every finite ROOT-stopped path obeys `t(n)=(-1)^d t(e)`, so a point
certificate gives `|a_N|<=B`. Independently,

\[
 a_{N+1}-a_N=\sum_{r=1}^{3}w(4N+r)t(4N+r),\qquad
 |a_{N+1}-a_N|\le4.                                      \tag{1}
\]

The omitted fourth entering source cancels the departing one because
`w(4N+4)=w(N+1)` and `t(4N+4)=t(N+1)`.

**A2 (PROVED).** For any nonempty finite set of anchors `S` in `[X,2X)`,
supply a point certificate at each a with residual B_a. The controllers,
terminal sets and depths may differ from anchor to anchor. Then

\[
 |a_N|\le C(N):=\min_{a\in S}\{B_a+4|N-a|\},\qquad
 D_X:=\sum_{N=X}^{2X-1}a_N^2\le\sum_N C(N)^2.             \tag{2}
\]

A1 gives literal point matchings behind every chosen bound, not only a
scalar Lipschitz argument. Choose a minimizing anchor on each physical
row and use its transported matching. Compress consecutive identical
pair/path records into intervals. These form a finite interval certificate:
each row has the exact original atom capacities, and different rows never
consume each other's capacity. Records may have different terminal
controllers on disjoint intervals. Each record retains both actual paths.
Thus every physical-bin refinement or coarsening remains bounded by
`sum_N C(N)^2`, by the triangle inequality and Cauchy--Schwarz.

No gluing theorem for Galois representations or probability measures is
being invoked. Compatibility is proved by the common literal scalar
functional and, more strongly, A1's explicit atom transport. If one also
wants a single terminal-vector Gram matrix, a common native extension must
be supplied and checked. Formula (2) itself does not require it.

**Sharpness of the information used.** Among all real scalar functions f
on the integer octave satisfying `|f(N+1)-f(N)|<=4` and
`|f(a)|<=B_a`, the exact largest possible square sum is `sum C(N)^2`.
Every f is bounded by C. Conversely C is nonnegative, four-Lipschitz,
and satisfies all the anchor upper bounds, so f=C attains it. This is a
sharp abstract envelope, not a claim that C is a native Collatz current.
Further improvement needs new arithmetic information or stronger anchors.

## 4. An explicit sufficient sampled bound

Partition `[X,2X)` into equal integer bins of length b dividing X. Choose
the lower middle row a of each bin, with any independently authenticated
controller there. Define

\[
 Q_S^2=\sum_{a\in S}b B_a^2.
\]

The distances to a within its bin have

\[
 S_1=\lfloor b^2/4\rfloor,\qquad
 S_2=\begin{cases}b(b^2-1)/12&b\text{ odd},\\
                  b(b^2+2)/12&b\text{ even}.\end{cases}
\]

**A3 (PROVED).** The exact assigned-bin score and its simpler upper bound are

\[
 D_X\le\sum_a(bB_a^2+8S_1 B_a+16S_2)
 \le Q_S^2+2b\sqrt X\,Q_S+\frac43X(b^2+2).              \tag{3}
\]

The first bound expands `(B_a+4 distance)^2`; the second uses
`S_1<=b^2/4`, `S_2<=b(b^2+2)/12`, and Cauchy--Schwarz.
The all-anchor minimum in (2) can improve this score.

Consequently `b<=sqrt(X)` and `Q_S<=5X` suffice for `D_X<=49X^2`
(with slack even at X=1). More generally,

\[
 b\le\sqrt X\,X^{o(1)},\qquad Q_S^2\le X^{2+o(1)}       \tag{4}
\]

suffice for `D_X<=X^{2+o(1)}`. Via I3 of
[interval matching](collatz_interval_matching_20261016d.md), using
[THM-4658, two-octave energy descent](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md),
this implies the square-root energy hypothesis (W) for H. The premise (4)
is **OPEN**. This reformulation needs roughly sqrt(X) anchor rows, permits
a different finite native controller at each, and imposes no reset-stage
or path-length deadline. It still demands a genuine all-scale bound on the
weighted residuals. (W) alone does not establish a nonzero residue or ROOT
entry.

The earlier [refinement-debt theorem](collatz_refinement_debt_20261016a.md)
uses means of vector rows in a common controller and already bounds the
cost of refinement. Here the input consists of sampled point certificates
from potentially different controllers. Neither input estimate follows
from the other without supplying its witnesses.

## 5. Why finite endpoint residues cannot pay this bound

For any positive Q divisible by six, set

\[
 a=Q+3,\quad b=2Q+3,\quad n=2a,\quad m=b.
\]

Both sources n and m lie in the physical annulus `(a,4a]`, have weight -1,
and have actual initial-halving receipts

\[
 n\longrightarrow a\quad(d=1),\qquad m=b\quad(d=0).
\]

Their weighted signs are opposite and their endpoints agree modulo Q,
but the endpoints are distinct. They cannot be cancelled by that supplied
controller. Q can contain arbitrarily many prescribed prime powers.
This is an exact native hostile to replacing the common endpoint by a
finite residue signature; it says nothing against a later actual merger.
A height sidecar repairs the elementary test: if endpoints differ by a
multiple of Q and their actual difference has absolute value <Q, they agree.

The geometric companion's local-to-global theorem works because its
unit-path ring is an overring of the full ring of integers, with all
conjugation-fixed prime conditions retained. Neither ring closure nor that
overring premise holds automatically for directed Collatz certificates.
The existing adelic-capacity obstruction in
[AMM/Collatz bridges](procgen_bridges_20260923_lrc_amm_collatz.md), section 2,
also remains in force.

## 6. Reproducible controls and next obligations

[Program](../../04-computation/experiments/collatz_interval_atlas_20261016e.py).
The policy was fixed before the runs: bin width `2^floor(k/2)`, lower-middle
anchors, and reset stage `6+(anchor index mod 3)`. Controllers are snapshots
of the inherited native graph; the theorem allows more general choices,
but these finite controls do not pretend to construct each snapshot from
scratch. Every source through `8X` is replayed at each used stage. The small
controls also replay literal atom matchings and transports in both
directions; the same-parity path issue is inherited from I2's hostile.

For k=6,8,10,12,14,16 the respective cone scores divided by X squared are
`4.109375, 4.29876709, 6.17216873, 16.96984172, 12.39897010, 27.86255435`.
All six are below 49. The larger checks then **refuted the uniform 7X
claim for this fixed schedule**, without refuting A1--A3:

| k | Fixed stages 6--8: cone score / X^2 | Repaired stages | Repaired score / X^2 |
|---:|---:|---|---:|
| 18 | 83.45985295 | 11,12,13 | 16.66345834 |
| 20 | 212.12580486 | 12,13,14 | 20.19465803 |

The repair was chosen after seeing the failures: use base stage
`floor(k/2)+2`, then base plus the anchor index modulo three. It is a finite
repair, not a preregistered prediction or a universal law. Its sampled
energies Q_S^2/X^2 are 9.62529311 and 13.01418127, both below the sufficient
threshold 25 in A3. The first tested failure of the fixed schedule is k=18
in this even-k universe; no least failing scale is claimed. The failed
implication would have been that a depth schedule adequate on six smaller
octaves stays adequate at every height. The missing coordinate is allowed
controller depth as the original source scale grows.

The larger fixed runs replay every original source at all three used
stages, through 8,388,608 sources at k=20. A second
[compositional verifier](../../04-computation/experiments/collatz_interval_atlas_repair_20261016e.py)
checks every original source's endpoint and integer-clock increment at every
extension, then literally replays each distinct old-terminal extension.
It agrees with the full-prefix verifier on all original endpoint/clock
hashes, all anchors, both cone hashes and the final full-row energies at
k=18,20. The repaired runs use this exhaustive compositional method; a
small independent full-prefix replay is also retained. Shared prefixes are
reused, not assumed or discarded.

Normal and optimized executions agree byte for byte for the six small
runs and for the two repaired large runs. The full and compositional large
fixed-schedule checks use different replay organizations. Their agreement
is computational evidence, not an independent mathematical referee.

Receipts: [small fixed atlas](collatz_interval_atlas_20261016e.json),
[large fixed failures](collatz_interval_atlas_20261016e_fixed_large.json),
[compositional comparison](collatz_interval_atlas_20261016e_incremental.json),
[scale-dependent repair](collatz_interval_atlas_20261016e_repair.json),
[SHA256 manifest](interval_cubic_audit_20261016e.json).

```text
python3 -B 04-computation/experiments/collatz_interval_atlas_20261016e.py --k 6 8 10 12 14 16
python3 -B 04-computation/experiments/collatz_interval_atlas_20261016e.py --k 18 20
python3 -B 04-computation/experiments/collatz_interval_atlas_repair_20261016e.py --k 18 20 --base 6
python3 -B 04-computation/experiments/collatz_interval_atlas_repair_20261016e.py --k 18 20
```

Add `-O` after `-B` for the optimized comparisons. The fixed failures remain
part of the saved evidence; no favourable subset replaces them.

The next arithmetic target is (4), or its stronger fixed-constant version:
build bounded-deficit anchors with original weights, then let A1--A3 cover
the gaps. A second useful question is whether a residue-class family of
anchors has a provable native pairing rule with a subpower exceptional
invoice. The hostile above requires exact endpoint or height information;
local densities alone are insufficient. A third is to optimize the choice
of anchor controller under a separate work budget, without charging that
budget as if it were a cancellation proof.
