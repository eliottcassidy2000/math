# Actual ROOT sources: clock correlations and a uniform bounded-clock certificate

**PROVED, elementary reductions and a computer-assisted bounded-clock theorem.**
**FINITE-EXACT:** the separate unrestricted-clock scale census.
**OPEN:** uniformity in the integer-clock cutoff, the full ROOT transport
bound, and universal Collatz entry. No priority claim is made for the
Fourier/triangle inequalities. Proofs are self-audited; the independent
implementations described below are not an independent human/subagent review.

The strongest conclusion here is uniform in the scale: for every integer
`N>=4096`, every native-depth cutoff, every integer-clock cutoff at most 16,
and every frequency, the three-phase energy of the actual first-arrival
source is at most its initial three-phase energy. The restriction to clock
16 is essential. Removing that restriction remains open.

## 1. Inheritance and the representation retained

Closest mechanism: [THM-4690, first-arrival ports and integer clocks](../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md),
and [E1--E4, the Eisenstein source lift](root_eisenstein_first_arrival_20261018j.md).
Canonical hostiles: the latter's expanding phase direction at N=524288 and
small-scale trace expansion at N=1204. Corrected near miss: a unitary model
of a packet does not authenticate the actual source preparation or terminal
readout; see [the cube clock packet](cube_unitary_clock_packets_20261018n.md).
Least-used sidecar: the complete, source-weighted clock autocorrelation.

The earlier [sparse boundary module, R3--R4](root_sparse_boundary_module_20261017i.md)
already gives all-scale dyadic error bounds at fixed native depth by an exact
period in log2(N). We do not re-prove that result. Here the source scale is
any integer N, the cutoff is the literal first-arrival clock, and the object
is the whole stopped-source correlation rather than an ideal-minus-actual
error vector. The [3-adic Wiener theorem W1/W2](collatz_conductor_wiener_20261015d.md)
already gives infinite-depth contraction for regular excess densities, with
an unbounded marked-atom cost. Our use of an absolute Fourier coefficient
sum is on the integer-clock energy polynomial. It neither places the source
in that density space nor evades its atom-cost obstruction.

The board is: actual Eisenstein charges; integer-clock correlation; finite
inverse trees; boundary selectors; marked finite geometry; source-fixed
retractions. Anchor: ROOT energy. Niche: selector logic. Wildcard: the
[5/27/21/9 geometries](marked_geometries_selectors_20261018o.md).
Incoming commits 190ea97bf4 and 269c5e124b were read in full in the relevant
four notes before integration. Their cubic retraction preserves all old
colourings; it supplies no signed-energy contraction. A colour pullback
and a charge pushforward have different consumers and fibre costs.

Use precisely the inherited primitive positive seeds, B(s)=4s+1, the
forward-closed 627-seed frontier V, and its integer seed-clock convention.
For each seed s the literal B-chain source is A_N(s) in Z[omega], with
omega^2+omega+1=0 and

    c_r(s)=2 Re(omega^r A_N(s)),
    Norm(A_N(s)) in {0,1,3},
    S_N=sum_s Norm(A_N(s))=T_source(N)/6.

The row r=0 remains the original ROOT source. An arriving seed retains
entry terminal u, last exterior entry port p, native depth dep, and integer
clock t. Depth-zero sources have no exterior port. Sources not meeting a
declared cutoff are excluded from that stopped polynomial, not declared to
have converged or to carry zero source weight. The denominator S_N includes
the entire original source family.

For a finite native cutoff d and clock cutoff h set

    Q_u(z)=sum_{s: u(s)=u, dep(s)<=d, t(s)<=h} A_N(s) z^t(s).

The unrestricted-clock finite-horizon version has the same definitions.
For uncompressed inputs B^j(s), the physical-clock endpoint correction of
THM-4690 must still be restored. This note never replaces that clock by its
parity or silently identifies the two conventions.

## 2. The exact correlation ledger and sufficient bound

Write q_(u,t)=[z^t]Q_u and c_(u,t)^(r)=2 Re(omega^r q_(u,t)).
Define the real integral polarization

    B((a,b),(c,e))=2ac-ae-bc+2be=2 Re((a+b omega) conjugate(c+e omega)),
    D=sum_(u,t) Norm(q_(u,t)),
    Gamma_l=sum_(u,t) B(q_(u,t),q_(u,t+l)),  l>=1,
    W=D+sum_(l>=1)|Gamma_l|.

All terminal labels remain separate before the indicated sums. Cancellation
between contributions to one lag is retained before taking its absolute
value; replacing Gamma_l by a sum of absolute source-pair products is a
strictly weaker estimate.

**R1 (PROVED).** For every real theta,

    Ecal(theta):=(1/6) sum_(u,r) |sum_t c_(u,t)^(r) exp(i t theta)|^2
              =D+sum_(l>=1) Gamma_l cos(l theta) <= W.             (1)

In particular Ecal(pi)=sum_u Norm(Q_u(-1)). This follows by expanding the
three real source rows and using sum_r c_r(A)c_r(B)=6 Re(A conjugate(B)).
At a nonreal frequency it is important to retain both real phase
coordinates: Ecal(theta) is not generally just sum Norm(Q_u(exp(i theta)))
after choosing one complex embedding of omega.

Thus the following source-specific sufficient target is meaningful:

    (WB) W_(N,d) <= (3/2) S_N for N>=4096 and every finite d. OPEN. (2)

It implies the inherited ET trace bound, simultaneously at every frequency.
It is stronger than a bound just at pi; no equivalence is claimed. It does
not assume independent source phases, arbitrary-input contractivity, or
universal entry. A terminal-only version of (2) is weaker than all prefixes.

**R2 (PROVED exact obstruction).** The stronger all-frequency contraction
Ecal(theta)<=S_N fails for the actual unrestricted-clock source at N=4096.
At z=(-35+12i)/37, an exact rational evaluation gives

    Ecal(arg z) > (1043/1000) S_N,   S_N=9949.

The complete numerator/denominator is in the receipt. Meanwhile the actual
Nyquist value is 9629<S_N. Thus failure is not an arbitrary-input hostile:
it uses the actual weights and first-arrival clocks. Any proposed unitary
realization with norm-preserving preparation and a contractive readout at
every frequency would contradict this witness. A paid readout, frequency
restriction, or larger budget is necessary. The witness does not refute WB.

The divided polynomial Q=z^tau[W_pi+(1+z)R], tau the first occupied clock,
is also measured. At N=4096 its total quotient ell-mass is about 3.005 S_N.
Parity cancellation is not a free removal of the intervening clock interval.
This agrees with the inherited exact odd-lag quotient-cost lemma.

## 3. A finite extraction theorem at each bounded integer clock

Let F_h be all primitive positive seeds whose first arrival in V has clock
at most h, without any upper bound on the seed's numerical size. Let
K_(u,h) count them at terminal u and K2_h=sum_u K_(u,h)^2.

**R3 (PROVED).** F_h is finite, computable by a complete inverse tree, and
for all N, all native-depth cutoffs, and every clock cutoff at most h,

    W <= 3 K2_h.                                                (3)

Every source charge has modulus at most sqrt(3). Expanding the absolute
pair bound separately at each terminal gives
W<=sum_u (sum_{s in F_h(u)} |A_N(s)|)^2<=3 K2_h. Removing any sources by
a native-depth or smaller-clock cutoff preserves this upper bound. This
bound discards cancellation deliberately; its purpose is finite extraction,
not a uniform estimate as h grows.

For completeness of the inverse tree, let y be a target seed. Put v=B^j(y).
If v=1 mod3, its possible predecessor is (4v-1)/3 with clock 2j+2;
if v=2 mod3, it is (2v-1)/3 with clock 2j+1; otherwise there is none.
Exclude predecessors already in V, which stop at depth zero. Keep only
remaining clock budgets. Positive edge clocks make this a finite tree.
These cases are exhaustive by reversing the literal forward equation.
Every generated source is separately replayed forward to its first V hit.

One can also bound its numerical support explicitly. For one native edge,

    3(3s+1)=2^t(3y+1)+(3-2^a),  a in {1,2}.

Hence 3s+1<=2^t(3y+1), and an entire path gives
s<2^h(max(V)+1). No convergence statement is used to enumerate F_h.

**R4 (PROVED).** For every positive integer N,

    S_N >= 15N/8-12.                                           (4)

Primitive seeds in [N,2N) have exactly one lower-band point and norm 3;
those in [2N,4N) have exactly one upper-band point and norm 1. Primitive
seeds occupy three residues modulo 8. The respective counts are at least
3N/8-3 and 3N/4-3. All other norms are nonnegative, giving (4).

Consequently, if

    N >= B_h:=ceil((8 K2_h+32)/5),                              (5)

then W<=S_N for every smaller clock cap and every native-depth cap. Therefore
an all-N contraction statement at bounded clock reduces to the finite range
4096<=N<B_h. This is an independently proved cutoff, not a stopping point
chosen because a scan looked stable. Its dependence on h is the unpaid part
when trying to pass to unbounded clocks.

## 4. Certified all-scale result through clock 16

**R5 (PROVED with the finite certificate below).** For this fixed frontier,

    W_(N,d,h) <= S_N
    for every integer N>=4096, every d>=0, and 0<=h<=16.         (6)

The complete inverse tree has 27,440 sources, largest source 463,596,771,
and K2_16=1,454,166. Its routes retain terminal, clock, depth and exterior
entry port. Formula (5) gives B_16=2,326,672. The finite consumer checks
every scale from 4096 to B_16-1 and all 17*17 pairs of clock/depth caps.
Depth caps above 16 coincide with cap 16 because every native edge costs at
least one clock unit. There are no violations. Beyond B_16, (3)--(5) prove
the result without further orbit computation.

The largest finite-part W/S_N is 6854/9949, at N=4096, clock cap 16 and
depth cap 10. This finite-part maximum is not asserted to be the maximum
over all N; the tail proof gives the upper bound 1 there. Equation (1)
makes (6) a simultaneous all-frequency bound on the actual source.

This is stronger in frequency and scale than a finite scan, and weaker in
clock scope than the desired theorem. It supplies a controlled initial
clock region. It does not pay the interaction between this region and the
unbounded remaining clock tail.

## 5. Separate unrestricted-clock evidence and next target

The C scan uses all first arrivals of seeds s<4N, with overflow guards and a
10,000-native-step failure cap. No route reaches the cap. It updates the
actual source at the three new threshold events 4N-3,4N-1,2N-1; the complex
increments are respectively omega^(m+1), omega^(m+1), -omega^(m+2).

**FINITE-EXACT:** for all 4,190,208 scales 4096<=N<2^22, final-arrival
W<=(3/2)S_N. The largest ratio is 38167/26722 at N=11002, approximately
1.42830. This bound was selected before this exhaustive scan. The observed
proximity to 10/7 does not establish that constant or explain it by the
different inherited 10N/7 sibling cancellation.

A separate complete-B-chain Python implementation agrees on the shared
dyadic and hostile scales, without using the three-event update. It retains
every arrival prefix at its declared scales. The stronger W prefix check
was performed at N=4096,5997,11002,16384; its respective maxima are
13193,19693,38167,54899. These are not every-N prefix proofs.

The remaining useful question is now quantitative: control the correlations
of clock increments beyond 16, including their cross terms with the certified
initial region. Appending separately bounded pieces need not preserve (6).
A sufficient route is a summable majorant for the added Gamma_l contributions
on the actual source, uniform in N and native depth. Generic independent-phase
or unitary-channel arguments cannot replace that source-specific estimate.

## 6. Reproduction and audit scope

Scripts in `04-computation/experiments/`:

```sh
python3 root_clock_channels_20261018o.py --powers 18 --out /tmp/clock.json
cc -O3 -Wall -Wextra root_wiener_scan_20261018o.c -o /tmp/root-wiener
/tmp/root-wiener first_arrival_20261018/V_closed_Y1000.txt 22 > /tmp/wiener.json
python3 root_clock_finite_certificate_20261018o.py --clock 16 --routes /tmp/clock-routes.txt --out /tmp/tree.json
cc -O3 -Wall -Wextra root_clock_finite_certificate_20261018o.c -o /tmp/clock-cert
/tmp/clock-cert /tmp/clock-routes.txt > /tmp/clock-certificate.json
python3 root_clock_certificate_audit_20261018o.py --out /tmp/audit.json
```

All committed receipts have matching names under `05-knowledge/results/`.
The audit reconstructs the clock tree forward, checks every saved peak using
literal B-chain weights and a separate lag aggregation, and compares both
unrestricted source engines on their shared samples. Normal/optimized checks
remain active. The finite loop exhausts its declared universe; selected-peak
replay is an independent check of that implementation, not a second exhaustive
scan of every scale. No independent subagent audit is claimed.

Research moves used from [META-PATTERNS](../../00-navigation/META-PATTERNS.md):
“Type a shared sum as a weighted fiber before transferring it” and
“Re-evaluate a certificate after every fibre-changing operation.” The exact
frequency hostile is the counterindication to transferring the cube norm
without its source/readout map. No new universal meta-pattern is promoted
from this one bounded-clock certificate.
