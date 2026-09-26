# 223 crossroads: results and remaining proof obligations

**Status: PROVED scoped results / CITED inputs / INDEPENDENTLY AUDITED /
FINITE-EXACT controls / OPEN Collatz.** Research session, 2026-09-26,
starting from 0e0301490cc8 and integrating concurrent main-line work.
This replaces the initial research-in-progress board.

The strongest progress is a sharp finite-horizon price theorem, a growing
orbit-prefix rank theorem, and an exact series for two-step pairing cost.
223 supplied useful objects and hostile examples; recurrence of the number
was never treated as a proof of a connection.

## Inheritance and portfolio

Anchor: descent, ordered affine carries, and private versus globally
consistent pairing certificates. Niche: prime223, Fermat sextics and
S-unit specialization. Wildcard: positive transfer matrices, Bernoulli
sine coefficients and the time scale of endpoint connections.

Closest inherited mechanisms: THM-4478 integer capacity,
[THM-4480 peak discount](../../01-canon/theorems/THM-4480-peak-discounted-provability-price.md),
the pairpeak barrier inequality and the LRC (2,2,3) labelled guard.
Hostiles: negative cycles, the previous carry collision at depth233
(a DIFFERENT number), and incompatible private demands from7 and11.
Corrected near miss: P1's polynomial comparison to unweighted bad density
is false for arbitrary edits/private pairing, and open for global pairing.
Recovered sidecars: endpoint phase, halving cost, actual pair ownership,
and prime support incident to one node.

## What the recovered occurrences of 223 meant

| Source | Exact role | Useful transfer or stopping reason |
|---|---|---|
| [THM-2204](../../01-canon/theorems/THM-2204-scalar-depth-223-thirteen-lift-capacity-law.md) | Depth profile (2,2,3), not prime223 | Keep labelled guard state; scalar capacity loses compatibility. |
| [THM-3174](../../01-canon/theorems/THM-3174-projected-k3-z223-terminal-descent-and-cap222.md) | Terminal LRC layer z1=223, cap222 | A finite threshold is not a Collatz drift invariant. |
| [THM-3334](../../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md) | Parameter223 has norm99905=5*13*29*53, eight Gaussian allocations | Factorization explains a record address; prime223 cannot divide these consecutive-square norms. |
| [HYP-3013](../hypotheses/HYP-3013-lrc14-perfect-number-packet-merge.md) | Prime223=14*16-1, deficiency12/223 | Divisor structure matters; composite controls behave differently. |
| [Harmonic shell note](prime_shells_20260921_primes.md) | Sum of the first six harmonic numbers is223/20 | Exact numerator, but no descent-preserving map found. |
| [THM-3394](../../01-canon/theorems/THM-3394-twelve-formerly-missing-hadamard-orders-through-2000.md) | Hadamard order892=4*223 | Row sums alone give no Collatz certificate. |
| [THM-3131](../../01-canon/theorems/THM-3131-prime-resonance-newton-slope-separation.md) | Endpoint of a prime census | Local noncancellation transfers; the endpoint does not. |
| Odd Collatz word(2,3,1,1) | Carry223=27+36+96+64; map(81n+223)/128 | Keep numerator, denominator, exponent order and positive integer realization. |

The incoming [recurrent-numbers atlas](constants_atlas_20260926_recurrent_numbers.md)
independently distinguishes addresses, invariants and coincidences.
Theorem IDs, bibliography pages and scan endpoints were not elevated
to arithmetic identities.

## Six concepts after the experiments

| Concept | Productive operation | Result or boundary |
|---|---|---|
| Finite residue fibres | Repeat a base period before judging closure | Complete lifts preserve cycle-average extrema and inherited weighted growth. |
| Affine carry and real scale | Enumerate valuation words, then realize by CRT | Every223-return of halving cost<=8 descends; a first return expands at9. |
| S-unit primes | Use two incident banks, then the factorial cost of many primes | Growing-prefix rank holds for almost all starts, uniformly over intervals. |
| Positive matrices | Round the mean to return the phase; use TP2 trace bounds | Sharp confinement lower bounds with explicit multiplicity. |
| Bernoulli coefficients and endpoints | Preserve coefficient signs and allow longer connectors | Sharp peak asymptotics for every fixed odd multiplier. |
| Shared pairing bits | Charge outside a paid matching; telescope tree costs | Global29/108 lower bound and an exact minimum-lower-density series. |

## Four canonical results

**[THM-4488: peak prices and private repair](../../01-canon/theorems/THM-4488-private-pairing-peak-price-rational-bridge.md).**
For q=3,

    2 rho_peak(L) <= pi_L <= O(L^3)rho_peak(L),
    N_m(L) <=4e(m+3)L A_(m+5)(L),  L>=3,m>=2.

For every fixed odd q>=3,

    log rho_peak,q(L)=L log D_q-kappa_q L^(1/3)+O_q(log L),
    c=log_q2, D_q=2^(-(1-H_2(c))),
    z_q=q/max(1,(1-c)/c),
    kappa_q=(3/2)(pi^2*c*(1-c))^(1/3)(log z_q)^(2/3).

The arbitrary-edit price has the same asymptotic. The proof combines
[sine supersolutions](crossroads223_20260926_robin.md),
[rational transfer bridges](crossroads223_20260926_bridge.md), and
[polynomial-cost endpoint connectors](crossroads223_20260926_general_peak.md).
Connectors of length M^2/log M cost only a polynomial in M and retain
the endpoint reward when q>=5. A forced O(M) suffix would change the
sharp constant. Mogul'skii is no longer needed for this particular walk.

Concurrent work supplied an independent [q=3 Robin proof](procgen_robin_20260926_robin_inequality.md)
with shift2 and factor(m+2)^17. Both estimates are retained. Neither
proves the constant comparison HYP-9142 or global-price HYP-9140.

**[THM-4489: Fermat sextics and223 returns](../../01-canon/theorems/THM-4489-fermat-sextic-223-collatz-return-budget.md).**
The projective curve x^6+y^6+z^6=0 has no F_p point exactly for
p=7,31,67,79,139,223. For nonzero coordinates the largest exceptional
prime is277. CITED Hasse-Weil/genus bounds reduce the proof to a complete
finite census, independently checked at all84 primes through433.
The missing sextic point gives one forbidden nonzero Collatz coset edge
modulo223, but no sign-selecting drift.

The word(2,3,1,1) gives the rational cycle223/47,179/47,73/47,133/47.
Every223-return with at most8 halvings descends;14495->...->185759
is an expanding first return at9. The [geometry note](crossroads223_20260926_geometry.md)
proves each claim at its actual finite or all-height scope.

**[THM-4490: uniform growing-prefix rank](../../01-canon/theorems/THM-4490-collatz-affine-word-sunit-specialization.md).**
For every interval of H positive integers, almost all odd starts have
full first-slot rank, hence paired rank, through

    N(H)=floor(sqrt(log_2 H * log_2 log_2 H)/8)

odd transitions, uniformly in interval location. Beukers--Schlickewei
gives fixed-word finiteness; incident carry banks and their factorial
prime cost make the bound uniform over growing word families. An exact
Terras tail handles large valuation totals. No practical small threshold
is claimed.

Among starts up to X with no descent through log_2 n half-steps, rank
failures number at most X^(0.46936094...+o(1)). Their entire population
has exponent0.949956..., so most difficult no-descent candidates ALREADY
have full growing-prefix rank. This limits a rank-only approach; it
does not prove descent. See the [growth proof](crossroads223_20260926_geometry_growth.md)
and [independent audit](crossroads223_20260926_geometry_growth_audit.md).

**[THM-4491: exact two-step lower-density price](../../01-canon/theorems/THM-4491-pairing-two-step-tree-extra-density.md).**
Global two-step constraints form a rooted tree. Their prefix optimum is

    C(X)/X -> alpha=sum_(h>=1) A_h/3^(h+1),

with explicit nonnegative integer periodic cost sums A_h. Alpha is the
attained minimum LOWER density among global two-step members;
natural-density attainment remains open. Twenty exact coefficients give

    0.2907539227382626... <= alpha <=0.2953650955221956... .

An independent simple cut gives29/108 by charging extra disjoint clauses
outside the old matching. Its clause fails in a global four-step pairing,
with exceptional path30->45->68->34->17. The earlier9/8 private-cost
penalty counted shared flips once per source, a different functional.
See [global ownership](crossroads223_20260926_flow_global.md) and
[the exact series](crossroads223_20260926_flow_series.md).

## Bernoulli and the original geometric reframe

The useful microcosm/macroscopic move was retaining boundary data while
changing scale. Prior B1 work showed why a surviving endpoint atom cannot
be discarded with vanishing bulk density. Here the sine coefficients are

    a_j=2^(2j-1)|B_(2j)|/[j(2j)!]>0.

Their signs control the remainder and preserve the sharp quadratic term.
Endpoint phase and terminal tilt then determine the usable bridge.
B1's exception alone is not a Collatz proof: the new deduction uses an
explicit operator and its boundary inequalities.

The pasted Kuratowski/Tutte/triangle reframe remains inspiration for
representations. No graph-minor, Pythagorean, Hadamard or finite-residue
quotient was allowed to discard height or repeated returns and still
claim convergence.

## Remaining obligations and failed extensions

- **OPEN HYP-9140:** make near-optimal private certificates globally
  consistent at polynomial peak cost. Shared-bit interference is missing.
- **OPEN HYP-9142:** a uniform constant Robin comparison, especially shift1.
- **OPEN G2:** independence on every injective chronological orbit, plus
  a height-sensitive consequence strong enough to force descent.
- **OPEN two-step refinement:** natural-density attainment of alpha;
  only minimum lower density is presently attained.
- **REFUTED extensions:** every223-return descends; a complete finite
  lift lowers inherited cycle growth; the extra two-step clause holds
  at four steps; fixed-word finiteness implies eventual rank at every
  start across all words. Exact witnesses and surviving statements are saved.

The incoming THM-4487 gamma=1 proof had a false n>=3n intermediate
inequality. An independently audited one-odd-step buffer repairs both
sheets and restores the polynomial lower prefactor. The concurrent block
repair is retained as lineage. See the theorem and the dated
[mistakes entry](../../01-canon/MISTAKES.md).

## Reproduction

Run python 04-computation/experiments/crossroads223_20260926_verify.py.
It runs all12 control programs normally and with -O, requires identical
outputs and writes [the artifact manifest](crossroads223_20260926_manifest.json).
Exact universes include all84 primes through433, 60,648 transfer minors,
617-bit prescribed cylinders, all3,597 short valuation words,65,534
Boolean assignments and twenty tree coefficients. Numerical sine and
probability probes remain labelled VERIFIED. General results rest on
written proofs and explicitly named cited inputs.
