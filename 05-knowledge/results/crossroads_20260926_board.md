# Collatz crossroads: proved transfers, synthesis, and live hypotheses

**Status: PROVED scoped results / CITED inputs / FINITE-EXACT controls /
OPEN Collatz.** Four research rounds on 2026-09-26, starting at
`884d207c9b12` in an isolated worktree and integrating `777e4e137`.
The main result is [THM-4478, critical growth bands and affine integer
capacity](../../01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md),
which proves HYP-9137: the optimal density of edits forcing descent within
L steps is `2^(-(1-H_2(log_3 2))L+o(L))`. This holds for arbitrary edits and
the restricted pairing family. It concerns modified maps with a uniform
finite horizon, not convergence of the unmodified map. No literature-priority
or Lean-formalization claim. The first two rounds below record the path to
the synthesis; their pending obligations are resolved in rounds three/four.

## Inheritance and portfolio

Anchor: Collatz first descent, retaining integer height and affine carry.
Niche: LRC Hall/arrival and AMM fair-extractor capacity mechanisms.
Wildcard: Rule 30 sections, analytic natural boundaries, Padé approximation,
and sparse orbit geometry.

Closest mechanisms: the exact Bernoulli divisibility jumps and repeated-block
budget in `bernoulli_boundary_20260925.md` and `kuratowski_reframe_20260925.md`;
THM-4476 thin injective orbits (under fresh audit); THM-4477 hub moments.
Canonical hostiles: the minus cycles, 5n+1 cycles, all-ones finite prefixes,
and a single surviving boundary atom of a cylinder of vanishing density.
Corrected near misses: finite-center linear ranks fail; projected marginals
lose source-hub incidence; eventual periodicity does not exclude extra cycles.
Least-used sidecars for this session: exact arrival incidence, a prefix
growth-band guard, the finite native-input tail, and the minimum of an entire
future orbit rather than a current-time minimum.

## Six live concepts

| Object / representation | Question and operation | Lost data / next decisive test |
|---|---|---|
| Bernoulli cylinder boundary | Cancel complete guarded packets | A surviving integer can live wholly in the boundary; retain affine threshold |
| Hall/flow source-hub incidence | Condition capacity on the actual undecided source | Second moments tilt toward rare large multipliers; trim growth and remeasure |
| Ballot words and fair extraction | Cycle-rotate near-critical blocks, then concatenate | Word density is not integer termination; use exact positive realizations as hostile |
| Rule 30 / Cartier sections | Track Collatz's online binary-to-parity sections | Finite step carries do not imply finite section memory; retain native terminal zero tail |
| Thin orbit / two-place height | Bound reciprocal tail above a future minimum | Local entropy count loses short-interval residue geometry; probe shell occupancy |
| Power series / Padé | Distinguish mixed scaling from genuine functional equations | Natural-boundary reformulations may be tautological; exact approximation certificate needed |

## Round one pulls (historical discovery state)

- Recovered THM-2163 (radix relation carries), THM-2200 (convex/semigroup/
  finite-place holes), THM-2580 (Hasse/Bockstein integral filling), THM-3152
  (multi-place Newton capacities), THM-4057 (Stern--Brocot depth), and the
  AMM rational-handoff obstruction. Shared words such as carry, rank or
  closure do not transfer their theorems to Collatz. Each needs its native
  operation and integral target.
- Rule 30's lossless current tree does transfer a section construction, but
  does not transfer finite memory to the Collatz itinerary.
- Incoming thin-orbit proof survives the main counting audit. Its harmonic
  corollary overstates full Collatz equivalence; correction is mandatory.
- Exact forward-orbit series has a mixed 2/3 scaling equation. The existing
  `procgen_continuous_20260925_natural_boundary_mahler_harmonic.md` already
  handles the nearby inverse-basin natural-boundary route; avoid calling the
  common Pólya--Carlson mechanism new or inferring two pure Mahler equations.

## Round two: incidence, tilt, and a repair (historical discovery state)

Conditioning hub capacity on actual bad sources improves finite constants.
The attractive hypothesis `E Q_L^2 <= poly(L) rho_L` is nevertheless REFUTED:
its diagonal weights tilt parity to Bernoulli(3/4), above the survival
threshold log_3(2), and leave energy bounded below while rho_L decays.

Repair under examination: keep only bad words whose every prefix multiplier
lies in `[1,K_L]`, with `log K_L=o(L)`. A cyclic-minimum block construction
PROVES there are `2^(h(log_3 2)L-o(L))` such words for
`K_L=3^(O(sqrt(L)))`. Thus this removes exponential self-weight without
changing the population's exponential rate. It does not yet bound collisions
of different sources at the same hub; that is the next exact test.

Concept comparison: the growth band is the missing scale sidecar in the Hall
lane; its exact word count uses the ballot lane; its realizability hostile
uses the Cartier/native-input lane; thinness forbids one infinite integer
from remaining in a fixed band but not these widening finite bands;
Bernoulli counts retain the source residue; Padé would be a separate way to
exclude particular infinite tapes, not a replacement for these guards.

## Round three: actual integer spacing closes HYP-9137

The missing coordinate is the normalized affine offset:

    T^k(n)=w_k(n+h_k),
    h_k=sum_(i<k, odd step i) 1/(3w_i),
    0<=h_k<=k/3.

At fixed hub v, time k and odd count e, all actual sources lie in
`[v/w_k-k/3,v/w_k]`. This interval has at most `floor(k/3)+1` integers.
A prefix-growth cap K permits at most `floor(log_3 K)+1` slopes per time.
Thus one edited hub serves at most polynomially many selected sources,
at a height at most `K(n+L/3)`. Every selected source must encounter an
edit before time L, since its original path never descends.

The resulting first-hit count proves a **lower natural density** bound
`rho_L(K)/(K M_L(K))`, without any density-existence or independence
assumption. Critical blocks retain the entropy rate; choosing block length
of order `sqrt(L log L)` gives

    edit density >= 2^(-eta L-O(sqrt(L log L))),
    eta=1-H_2(log_3 2)=0.050044... .

The inherited pairing upper bound is `2 rho_L`; an arbitrary edit map
has upper bound `rho_L` by sending the actual finite-horizon bad set to 1.
Both optima therefore have exponent eta. The proof was independently
audited by the geometry and automata lanes and promoted as
[THM-4478](../../01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md).

This bypasses the moment barrier by using the incidence data that moments
discard. It does not prove descent for the original map. Constant-factor
or polynomial-factor comparability with rho_L is also not proved.

## Round four: independent arithmetic and geometry gains

**PROVED special-tape certificates.** Thue--Morse and Rudin--Shapiro
selectors between equal-weight blocks `1^(L-1)0` and `1^(L-2)01` yield
explicit Bernstein/Mahler bridges. Exact 2-adic Padé errors beat ordinary
rational height in certified ranges. The displayed Thue--Morse certificate
works for q=3, every L>=3; the Rudin--Shapiro certificate works for q=3,
3<=L<=117. L=118 is a method boundary, not a rationality claim or a claim
about the literature. General Collatz tapes are not classified by these
families. [Automata proofs and controls](crossroads_20260926_automata.md).

The final cross-lane test sharpened this boundary: the two admissible blocks
`11111110000` and `11111101000` both have multiplier 2187/2048.
Periodic selection of the first has a rational source; Thue--Morse
selection has an irrational source by the certificate. Their endpoint
growth clocks are identical. Ordered carry, which the main capacity proof
retains only as an interval, is indispensable for source rationality.

**PROVED thin-tail geometry.** Audited THM-4476 implies any injective tail
above M has reciprocal sum `O_gamma(M^(gamma-1))`, h*<gamma<1. At a future
minimum its normalized real center lies in `[M,M+O_gamma(M^gamma)]`.
This is sublinear width, not a rounding certificate. Bounded products also
repair and strengthen the one-sided discrepancy obstruction in Corollary 9
of the [thin-orbit note](collatz_thin_20260925_thin_divergent_orbits.md).
Incoming `777e4e137` independently repaired the harmonic/cycle wording;
this session preserves that lineage and supplies the explicit refinement.

**CITED plus elementary consequence.** Beukers--Schlickewei's unit-equation
bound forces the multiplicative rank of N distinct paired transitions
`(-3m_j,3m_j+1)` to be at least `(log_2 N)/8-1`. Keep both coordinate slots.
Full rank is a stronger open question. [Geometry and primary source](crossroads_20260926_geometry.md).

**PROVED general capacity lemma.** The affine interval argument extends to
a fixed finite alphabet of positive affine slopes. A lower prefix-slope
bound controls offset width, and a cap controls endpoint height. For a band
that allows contraction this is a hitting-set theorem: a reason that every
path must hit an edit is an additional hypothesis. Sources 4,5,6 reaching
hub 2 at time 5 demonstrate why the lower slope guard matters. [Flow](crossroads_20260926_flow.md).

## Recovered-route ledger

These are mechanism transfers, not assertions that old theorems solve
Collatz. The detailed lane notes give source, target, map, preserved
predicate, destroyed information, required sidecar, and decisive test.

| Recovered route | Transfer or stopping reason | Detailed record |
|---|---|---|
| THM-2545 Hall arrival; THM-2549 future pullback | Retain source-to-hub incidence; separate margins lose the target | [Flow](crossroads_20260926_flow.md) |
| THM-3044 correspondence; THM-2177 unsplittable flow | Count actual deterministic integer ancestors, not formal fractional routing | [Flow](crossroads_20260926_flow.md) |
| AMM fair extraction and rational handoff | Squaring tilts parity to 3/4; trimming succeeds only with a population proof | [Flow](crossroads_20260926_flow.md) |
| Cyclic words / ballot counting | Rotation controls every prefix; fixed block boundaries preserve distinct words | [Ballot](crossroads_20260926_ballot.md) |
| Affine cocycles | Divide carry by slope to recover a short interval of integer sources | [THM-4478](../../01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md) |
| THM-2163 radix carries; Rule 30 sections | Keep both scale and carry; finite step carries do not imply finite full-itinerary memory | [Automata](crossroads_20260926_automata.md) |
| Cartier / coupled Mahler systems | Equal-weight substitutions close a special finite system; retain the transverse channel | [Automata](crossroads_20260926_automata.md) |
| Padé / product-formula reasoning | Exact valuation gain must exceed ordinary height; symbolic contact alone is insufficient | [Automata](crossroads_20260926_automata.md) |
| THM-3152 multi-place Newton capacity | Keep both prime-coordinate slots; joint transition rank is testable | [Geometry](crossroads_20260926_geometry.md) |
| S-unit equations | Finite-rank bounds force rank growth; full chronological rank is open | [Geometry](crossroads_20260926_geometry.md) |
| THM-3743 lattice flatness | Every finite word has arbitrarily large integer lifts; no lattice-free body exists without a height restriction | [Geometry](crossroads_20260926_geometry.md) |
| THM-2991 delayed Newton return; THM-3023 ratio transforms | Local positivity is not global return; killing geometric factors also kills relevant drift | [Geometry](crossroads_20260926_geometry.md) |
| THM-4057 Stern--Brocot / rational approximation | Near-critical rational block frequencies improve finite certificates, not cycle exclusion | [Ballot section 5](crossroads_20260926_ballot.md) |
| THM-2200 support holes; THM-2580 Hasse--Bockstein | Shared words like carry/closure give no map by themselves; no theorem imported | This ledger |
| Pólya--Carlson / natural boundaries | Mixed 2/3 scaling is not two separate pure Mahler equations; analytic continuation is missing | [Existing analytic route](procgen_continuous_20260925_natural_boundary_mahler_harmonic.md) |
| Zsigmondy primitive primes | Fixed-base hypotheses fail for moving Collatz coefficients; fresh-prime-per-rise is false | [Geometry](crossroads_20260926_geometry.md) |

Comparison across the live board: ballot counting provides enough sources;
flow provides incidence; affine carry provides integer spacing; Bernoulli
boundaries prevent loss of endpoints; section constructions check finite
realizability but supply no infinite integer lift; thinness constrains one
infinite orbit without ruling out growing finite families. Padé independently
excludes special tapes. No universal inference is borrowed across these gaps.

## Live hypotheses and next decisive tests

The following session-local labels are research targets, not reserved IDs.

| Label / status | Exact proposed strengthening | Evidence and next test |
|---|---|---|
| C1: REFUTED fixed-slope injectivity | Proposed uniqueness at fixed k,e,v with every prefix slope>=1 | It passed all 654,279 symbolic words through k=24 and 2,388,567 bounded-source prefixes. Exact carry width proves it through k=31, but coordinated odd-position changes give a collision at k=233,e=153 between sources differing by 4. This is not a minimality claim. The polynomial interval bound survives. [Flow](crossroads_20260926_flow.md) |
| G1: OPEN logarithmic shell occupancy | Every injective positive orbit segment has <=C log(2M) points in [M,2M), with uniform C | Constant occupancy is structurally false by balanced-word lifts. Finite complete paths from 50,000 odd sources support the logarithmic target. Would improve center width to O(log M), still not force descent. [Geometry](crossroads_20260926_geometry.md) |
| G2: OPEN chronological valuation rank | N paired transition vectors from a finite injective positive plus-orbit segment have rational rank N | Full rank certified modulo a verified prime on 50,000 complete tested paths. Eight unrelated injective edges have an exact dependence: chronology is essential. [Geometry](crossroads_20260926_geometry.md) |
| P1: OPEN polynomial price gap | Is delta_L >= rho_L/poly(L)? | The proved exponential rate leaves exp(O(sqrt(L log L))) loss. Exact block-gain optimization improves finite constants but does not close this. Need a sharper population/capacity balance. [Ballot](crossroads_20260926_ballot.md) |

A separate bounded method problem is to beat Rudin--Shapiro's Padé ratio
11/7. The full-rank even-denominator ansatz through degree parameter 32
gave no improvement; singular systems and other shapes remain untested.
A success would enlarge the special-tape catalogue, not prove the universal
Periodicity Conjecture.

The final C1 witness reconnects directly to the earlier negative-center
work. Set `u=3^(-153) mod 2^80=186937257649965781719819` and
`N=2^153u-1`. Then N and N-4 merge after 233 steps with 153 odd steps.
For their first 153 parity symbols they shadow, respectively, the negative
fixed point -1 and the cycle -5,-7,-10. Their normalized ordered carries
are `1-(2/3)^153` and `5-(2/3)^153`. Thus a four-unit offset, hidden by
the identical endpoint slope, matters at a huge integer height. This is
an exact version of the micro/macro warning, not merely an analogy.

## Failure ledger, reproduction, and handoff

- Harmonic divergence distinguishes periodicity from divergence, not the
  allowed cycle from other cycles. The upstream repair is retained.
- An endpoint-period average is not automatically a source-period average:
  the pulled-back harmonic observable needed an extra `2^(L-1)` factor.
  Draft values were withdrawn before checkpoint; repaired output checks
  the exact change-of-variables identity.
- Fresh primes need not appear on each rise: `55->83->125` introduces none
  at 125. Generic paired rank also fails on an exact eight-source product
  identity; the geometry note retains its minimal displayed witness.
- Infinite native-input sections do not forbid a weaker first-descent
  quotient or a finite-description proof. The scope is the full converter.
- Natural-boundary facts, density decay, and finite integer lifts do not
  show that a surviving nested cylinder is an ordinary integer.
- Changing the rule for every horizon does not prove convergence of the
  original rule. The settled item is HYP-9137, not Collatz.

Detailed correction lineage: [MISTAKES](../../01-canon/MISTAKES.md).
All six `04-computation/experiments/crossroads_20260926_*.py` scripts use
the standard library; matching `.out` files are here. Universes, filters,
controls, and interpretation limits are in the lane notes. The automata
and independent price audit also pass under `python3 -O`; ballot controls
use ordinary assertions and must run without `-O`. A SHA256 manifest
records exact script/output bytes.

Best next work: classify multiplicities in C1's admissible ternary carry
code, then a height-sensitive proof of G1 or G2. Keep actual chronology and the integer source. The
successful move was to retain an exponentially large manageable population
and count its actual integer incidence; an aggregate marginal or analogy
does not retain that information automatically.
