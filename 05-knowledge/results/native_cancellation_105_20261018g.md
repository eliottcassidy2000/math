# Native zero packets, compatible registers, and fixed-source scale transport

**PROVED:** exact removal of source-authenticated zero packets and transport of
each fixed packet from integer scale N to every 64^k N, and the exact
new-annulus source norm153N. **FINITE-EXACT:** the
four complete source universes below. **OPEN:** uniform native cancellation,
the new outer-source interaction, a global first-arrival bound, and universal
Collatz entry. No optimality or historical-priority claim.

Artifacts: [checker](../../04-computation/experiments/native_cancellation_105_20261018g.py)
and [complete compact receipt](native_cancellation_105_20261018g.json).

## 1. Inheritance and retained source

The charge is exactly E1 of [Eisenstein source charges](root_eisenstein_first_arrival_20261018j.md),
as implemented in the inherited
[source-clock producer](../../04-computation/experiments/root_clock_channels_20261018o.py).
The first-arrival frontier is the same 627-point set from
[THM-4690, first-arrival ladder ports and integer clocks](../../01-canon/theorems/THM-4690-first-arrival-ladder-port-and-integer-clock-decompositions-and-the-nineteen-adic-copy-obstruction.md) and the
inherited [boundary transport](root_sparse_boundary_module_20261017i.md).
The source-specific pairing principle is already present in
[interval cancellation](collatz_interval_matching_20261016d.md). We specialize
it to the actual two-coordinate charges, retaining an exact triple option.

The immediate corrected near miss is [weighted compatible clocks](weighted_crt_20261018g.md):
compatible mod15/mod35 registers, even with small or zero marginal coupling,
do not prove native cancellation. Here the native routes and signed charges
are checked directly. The least-used sidecar is the unabridged integer clock.
The live objects are the original source, its two compatible registers, the
Eisenstein charge, first arrival, packet capacity, and scale boundary.

For integer N>=1 the original source set is

    S_N = {s: 1<=s<4N, s odd, s mod8 !=5}.

Put B(m)=4m+1. Sum the three-real-row weights (0,-3,3), rotated by m mod3,
at B-chain points below 2N, and (-1,-1,2) at points in [2N,4N).
Writing the result as (2a-b,-a-b,-a+2b) defines A_N(s)=(a,b).
Its norm is

    Q(a,b)=a^2-ab+b^2 in {0,1,3}.

After the clock sign, the possible nonzero vectors lie in the twelve-vector
norm-1/norm-3 alphabet. The unsigned charges have nine nonzero possibilities:
the six norm-3 vectors and (0,1),(-1,-1),(1,0). The independent checker
instead sums coordinate weights (1,2),(-2,-1),(1,-1) in the lower band and
(0,1),(-1,-1),(1,0) in the upper band.

The canonical primitive map P uses a=2 when s=1 mod8 and a=1 otherwise, takes
(3s+1)/2^a, and repeatedly replaces a value 5 mod8 by (value-1)/4, adding two
to the clock each time. First arrival at V gives terminal u_s, integer clock
t_s and depth d_s. These are P-edge depths and the inherited integer braid
clocks, not the length or valuation cost of a literal fully odd U-word.
No modular-clock substitution is made. Sources already
in V have depth and clock zero. The file's supplied frontier rank is also
checked: every non-1 frontier point maps to a strictly smaller rank, with
rank(1)=0. The energy below stops at V, rather than aggregating all V to 1.

The rows retain s, s mod15, s mod35, u_s, t_s, d_s, a, b and Q(a,b),
including zero-charge sources. The two registers agree modulo5 and encode
s mod105 jointly; they are not treated as independent choices.

## 2. A zero-packet compiler at the actual Nyquist readout

Define v_s=(-1)^t_s A_N(s) and W_u=sum_{u_s=u}v_s. A packet is a pair or
triple of distinct named sources with a common actual terminal and sum v_s=0.
Its source norm is mu(P)=sum_{s in P}Q(A_N(s)). Polarization gives

    sum_{s<t in P} [2a_s a_t-a_s b_t-b_s a_t+2b_s b_t]
       = -mu(P),                                      (1)

where the coordinates in the bracket are those of the signed vectors v_s.
Thus the packet's complete cross-term debt is exactly paid at z=-1.
Removing any disjoint collection of such packets preserves **every** W_u,
and therefore preserves sum_u Q(W_u), and each of its three real readouts.
Its cross terms with the remaining sources also vanish, since its total
vector in that same terminal fibre is zero. Distinct terminals cannot be
silently merged to manufacture a packet. Disjointness prevents spending a
source's norm twice.

At a finite P-depth horizon d, deleting a whole packet requires
d >= max_{s in P} d_s. Before that horizon, its already-arrived members and
its still-moving members remain in their respective ledgers; a partial
arrival is not a zero packet. The table below uses complete first arrival.

The deterministic algorithm processes terminals in increasing order. It
extracts opposite-vector pairs first, then zero-sum triples in lexicographic
charge order, choosing the least remaining source in each charge bucket.
The twelve-vector alphabet has ten unordered zero-sum triple types.
No maximum-packing claim is needed or made.

The [primitive-circuit companion](eisenstein_circuit_completion_20261018g.md)
extends these same authenticated rows to packets of four to six sources.
It retains source multiplicities and complete-arrival cuts, and distinguishes
primitive-circuit exhaustion from exhaustion of all zero-sum submultisets.

At N=4096, the first pair is

| s | mod15 | mod35 | terminal | clock | A_N(s) | source norm |
|---:|---:|---:|---:|---:|---|---:|
|4257|12|22|3|23|(1,2)|3|
|4791|6|31|3|20|(1,2)|3|

Their signed vectors are opposite; (1) is -6. At z=1, however, their sum
is (2,4), with norm12. Their clock polynomial is (1,2)z^20(1+z^3), which
is not identically zero. This is a concrete barrier to all-frequency use.

The first triple has sources 6977,9303,13955, common terminal15 and clocks
12,14,13. Its signed vectors are (1,-1),(0,1),(-1,0). They sum to zero,
and its source norm5 is paid by cross debt -5.

## 3. The same packet survives every 64-fold scale change

**Fixed-source transport lemma.** For integer N>=1 and every s in S_N,

    A_(64N)(s)=A_N(s).                                 (2)

Consequently, any authenticated zero packet at N is the same zero packet
at 64^k N for every integer k>=0, with the same sources, routes, terminals,
integer clocks, individual source norms and debt (1).

This is a direct corollary of the inherited E1 weight definition; it does
not require first-arrival discovery at any larger scale. Indeed

    B^3(m)=64m+21,   B^3(m)=m mod3.

For integer m, B^3(m)<256N iff m<4N, and B^3(m)<128N iff m<2N.
Thus every old term reappears three B-steps later with the same weight and
band classification. The first three new terms are all below128N, since
B^2(s)=16s+5<=64N-11<128N. Their residues cover all three classes modulo3,
so their lower-band weight vectors sum to zero. This proves (2); induction
proves the claim for 64^k. The paths depend on s, not on N.

The factor cannot simply be changed to4. At N=1, s=1 has charge (-2,-1),
whereas at 4N its charge is (-1,-2). More directly, the actual N=4096 packet
{1137,6065}, terminal3, is zero at N but has signed sum (2,1), norm3, at4N.

The old sources remain fixed. All newly admitted sources in
[4N,4*64^k N) and their interactions with the retained vector are **unpaid**.
Equation (2) does not replicate one packet's mass into new source labels,
bound the energy of the enlarged universe, or prove a uniform arrival rate.

### Exact new-annulus norm and the remaining renewal bill

Let S(N)=sum_{s in S_N}Q(A_N(s)). The leading coefficient17/7 was already
proved, with a logarithmic error, in Proposition E3 of
[the inherited Eisenstein note](root_eisenstein_first_arrival_20261018j.md).
The same charge definition gives the sharper exact scale identity

    S(64N)-S(N)=153N,                                  (3)
    S(64N)-17(64N)/7 = S(N)-17N/7.

Here is a source-only proof, independent of arrival or cancellation. A new
source has4N<=s<256N and at most three B-chain terms below256N. Write U_r
for the upper-shell coordinate weight at residue r, and
L_r=U_r-U_(r+1) for the lower-shell weight. The consecutive residues are
r,r+1,r+2. Integer endpoints give exactly these six cases:

| new-source interval | shell pattern | charge norm | primitive source count | source norm sum |
|---|---|---:|---:|---:|
|[4N,8N)|LLL|0|floor(3N/2)|0|
|[8N,16N)|LLU|1|3N|3N|
|[16N,32N)|LL|3|6N|18N|
|[32N,64N)|LU|1|12N|12N|
|[64N,128N)|L|3|24N|72N|
|[128N,256N)|U|1|48N|48N|

For example B^2(s)=16s+5 is below128N precisely when s<8N,
and below256N precisely when s<16N; the analogous B(s)=4s+1 cutoffs
are s<32N and s<64N. There is no fractional-endpoint relaxation.
The LLL sum is zero; LLU and LU telescope to U_r; LL is
U_r-U_(r+2), of norm3. The remaining norms follow directly.
Every positive-norm interval starts and ends at multiples of8, so the
primitive residues1,3,7 mod8 give the displayed exact counts. Thus the
new source norm is63N+90N=153N. The old sources have their original
norms by(2), proving(3).

In particular, along each fixed integer scale ray N,64N,64^2N,...,
the centered source norm is constant and S(64^kN)/(64^kN) tends to17/7.
This source identity is not a terminal-energy estimate.

The remaining bill can be stated exactly. At a fixed finite P-depth d,
let W_u(N,d) contain only original sources already arrived at u by d,
and let H_u(N,d) contain only new sources4N<=s<256N arrived by d.
The unchanged old routes and(2) give

    W_u(64N,d)=W_u(N,d)+H_u(N,d),
    E(64N,d)-E(N,d)
      = sum_u [Q(H_u)+cross(W_u,H_u)].                 (4)

The bilinear cross expression is the bracket in(1). All the new sources
have total available source norm153N; their still-moving members have
not thereby arrived or cancelled. For a fixed c>=0, an increment bound
by153cN in(4) would transport E(N,d)<=cS(N) to the next scale. Iterating
such a bound only treats the stated64-scale ray; finitely many base scales
do not cover all integers. Neither the increment bound nor new-source
entry is proved here.

## 4. Exact finite results

The declared universe is every source in S_N at exactly the following four
scales. Discovery stops with an error on a repeated state or an unresolved
route of more than10,000 native edges. There were no unresolved routes.

| N | sources | pairs | triples | removed source norm | residual source norm | unchanged terminal energy |
|---:|---:|---:|---:|---:|---:|---:|
|4096|6144|1857|130|6404|3545|9629|
|5997|8996|2923|185|10167|4395|14420|
|11002|16503|6007|289|20611|6111|25611|
|16384|24576|9554|352|32348|7442|35299|

These are Eisenstein norm energies; the sum of the three squared real
readouts is six times the terminal energy. The terminal energies agree
with the inherited producer's Nyquist values. Residual source norm can be
much smaller than terminal energy: cancellation of selected packets does
not make the remaining transport contractive.

The receipt retains all56,219 source rows, every whole selected packet,
all terminal vectors, and38,975 shared route-DAG nodes. Each DAG node keeps
its actual edge and remaining arrival depth/clock. Remaining depth decreases
strictly along every edge until V, so the DAG authenticates first arrival
without a convergence assumption. Source/packet/DAG SHA-256 digests and
the inherited producer hash accompany the full data, rather than replacing
it. No external information is fetched during reproduction.

## 5. Reproduction and controls

```text
python -B -X utf8 04-computation/experiments/native_cancellation_105_20261018g.py
python -B -O -X utf8 04-computation/experiments/native_cancellation_105_20261018g.py --verify
```

Checks include independent shell charges for all four complete universes;
every actual packet and terminal vector; every retained route edge, suffix
clock and decreasing depth; every old source under N->64N; all primitive
sources at integer scales1..128; and selected fixed sources at
N=10^40+7 under64,64^2 and64^5. Wrong residues, clocks, terminals, weights,
Boolean/float identities, nonprimitive/out-of-scale sources and packet
double spending are rejected. The packet-at4N and non-Nyquist examples are
retained as exact hostiles. Normal and optimized runs use explicit checks,
not Python assertions.

The annulus delta additionally checks every new primitive source for
N=1,...,32 against independent shell charges, with exact per-band counts,
and directly sums the entire source universe at4096->262144:
S(4096)=9949 and S(262144)=636637. Analytic integer band counts also run
at the four inherited scales and N=10^40+7, with direct endpoint controls.
No new first-arrival searches are performed for these outer sources.

The finite certificate and scale lemma leave a specific next obligation:
control the original weighted terminal vector of the newly exposed sources,
including its cross terms with the old vector. A small mod15/mod35
projection statistic or repeated use of (2) does not pay that obligation.

Normal and optimized verification both pass **281,135** explicit checks.
The LF-normalized saved JSON SHA-256 is
`d2d4d152593a03eb653d1fd6fecfbc13f573a6eb74e1310f59b7a6cb78b01721`.

Independent proof/API review and optimized saved-receipt replay pass. The
review of the original packet also distinguished the nine unsigned charge
values from the twelve clock-signed alphabet values; the prose above retains
that distinction. The exact annulus delta received a separate independent
proof/API review and optimized saved-receipt replay at281135 checks; both
passed, including the integer endpoints,153N count and remaining energy bill.
