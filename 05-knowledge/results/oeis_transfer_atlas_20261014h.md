# OEIS provenance atlas and an exact run-weight compiler

Status: **PROVED** fixed-weight recurrence below; **FINITE-EXACT** independent
controls and a fixed-snapshot lexical reference census; **HISTORICAL PROPOSALS**
where indicated. This is not a certificate that every local novelty claim is
true, or that any draft was accepted by OEIS.

## 1. What the inventory establishes

The companion [JSON](oeis_transfer_atlas_20261014h.json) inventories every
`A` followed by six digits in tracked `.md`, `.txt`, `.py`, `.out`, and `.tex`
files under `01-canon`, `03-artifacts`, `04-computation`, `05-knowledge`, and
`07-reflections`, at the fixed Git snapshot `de52220938d0`. It excludes files
named `INDEX.md` and `MISTAKES.md`, and this atlas itself. The universe has
**53,338 files, 340 distinct A-identifiers, and 1,843 reference-bearing files**.
Each identifier retains every source path. This includes historical or mistaken
references: the lexical index is a discovery aid, not a mathematical authority.
Agent messages, logs outside these directories, binary artifacts, and future
commits are outside this census. The compiler uses the fixed object so concurrent
edits cannot silently change the reported universe.

Three contribution types must remain separate:

1. A published OEIS entry or extension needs an accepted entry/history identifying
   the contributor. A local `A??????` draft or a statement "should submit" is not
   that evidence.
2. A new interpretation, proof, or faster algorithm for an established sequence
   need not create a new sequence.
3. A finite candidate table needs a stable definition, indexing, enough checked
   terms, and a fresh duplication check before a novelty or submission claim.

The targeted provenance pass found **no independently verified accepted entry or
extension attributable to this repository's owner**. This is an evidence limit,
not a claim that no such contribution exists. On 2026-10-14 the direct OEIS
author/sequence-search endpoints were inaccessible through the browser tool;
ordinary entry pages were readable. Historical "not in OEIS" assertions have
therefore not been upgraded to current absence claims. No submission was made.

### Explicit drafts and proposed families

| Object and indexing | Recovered artifact | Honest status / retained boundary |
|---|---|---|
| $W_n=P_n(2,0)$, $n\ge1$: $1,2,8,32,158,928,6350,49752,\ldots$ | [W draft](../../03-artifacts/oeis/W_n_submission.txt) | Explicit `A??????` submission draft. Its exact definition and first twenty terms are independently checked here; publication unverified. |
| $N(n,j)$: no-unit-descent permutations with $j$ unit ascents | [THM-219](../../01-canon/theorems/THM-219-nud-poisson-euler.md), W draft | Proposed triangle underlying $W_n=\sum_j2^jN(n,j)$. This packet proves the two-variable evaluation mechanism, not all assertions in the old synthesis. |
| Number of distinct skew-Seidel characteristic polynomials `cpK`: $1,1,1,2,2,6,11,50$ | [cpK/cpA draft](../../03-artifacts/oeis-submission-distinct-cpK-klein-S317.md) | Draft's novelty language is superseded by [THM-2010](../../01-canon/theorems/THM-2010-new-tournament-invariant-sequences.md), which identifies existing Seidel-matrix literature. An absent A-number would not restore mathematical novelty. |
| Number of distinct adjacency characteristic polynomials `cpA`: $1,1,2,3,9,28,168,1523$ | Same draft | Candidate; acceptance and present duplication status unverified. |
| Minimum crossings of rotation-invariant two-page drawings of $K_n$, odd $n=5,\ldots,17$: $5,21,72,176,377,705,1224$ | [two-page draft](../../03-artifacts/oeis-min-sym-2page-submission-boxeph.md) | Restricted drawing problem, not the unrestricted crossing number. Reported finite values are inherited, not recomputed here; the draft's asymptotic prose is not adopted. |
| Ten tournament invariant families: adjacency spectra, achievable $H$, cycle data, signed sums, discriminants/spectral variance, arborescence data, maximum signed sum, metagraph edges, and 1-WL classes | THM-2010, above | Reference catalogue from a finite isomorphism census, not ten accepted OEIS entries or ten novelty theorems. |
| Rotational-tournament $H$, minimum $H$ among circulants, strongly-connected $H$ spectra | [circulant catalogue](maxH-circulant-technique-and-new-sequences-opus-20260629.md) | Family counts/lower bounds. Searching only circulants does not certify the unrestricted maximum A038375. A failed bounded recurrence search is not a proof of non-holonomicity. |
| Strong/non-strong fixed-path tilings and exact numbers of bad cuts | [extended tables](new_sequences_extended.out) and [producer](../../04-computation/new_sequences_extended.py) | Defined restricted ensembles. Strong connectivity, self-converse symmetry, and strong-component atoms must not be conflated under the abbreviation `SC`. |
| Arc-reversal metagraph edges, degrees, loops, triangles, score fibres, $H$ support | [metagraph tables](Gn_new_sequences_s169.out) | Finite sequence candidates, retaining the quotient graph and loop/weight convention. |
| Complement-odd first Betti numbers $0,1,7,119,1772$ and signed level multiplicities | [HYP-3565](../hypotheses/HYP-3565-b1minus-closed-form-and-new-proof-sequences.md), [HYP-3540](../hypotheses/HYP-3540-metagraph-level-multiplicity-sequence.md) | Historical candidates. HYP-3540's current header retracts the unsigned edge-subset interpretation: the action has signs. No claim of impossibility of an elementary formula is imported. |
| Four Faulhaber shear rows, shear/product variants, with Fibonacci/Narayana controls | [shear catalogue](shear_catalog_kps_S128c103.out), [producer](../../04-computation/shear_catalog_kps_S128c103.py) | Historical unsubmitted-family claims; each slope and sum/product operation is a distinct sequence. Short prefixes alone do not establish novelty. |
| Variance numerators/denominators and OCR ratios | [historical variance table](prime_transcendence_synthesis.md) | Exact rationals are computable from $W$; the historical universal small-prime cancellation and "transcendence" rhetoric is not a dependency. |
| Paley $H(T_p)$ and affine-orbit quotients $H(T_p)/[p(p-1)/2]$ | [sharing proposal](../../03-artifacts/drafts/research-synthesis-2026-03-10.md) | Proposed sequences indexed by primes $p\equiv3\pmod4$. An orbit quotient requires the free action and the exact source tournament; no accepted A-number was recovered. |

The JSON also lists all filename routes containing `oeis` or `new_sequences`.
This includes historical lookup outputs and stored OEIS b-files. Those copies
are evidence of use or comparison, not evidence of contribution. The table is a
curated proposal inventory; the full 340-ID reference index is the broad recovery
layer. It does not assert that every mathematically possible sequence implicit
in the repository has been named.

### Established entries and computation work

Live primary pages confirm that [A000568](https://oeis.org/A000568) counts
unlabelled tournaments and [A038375](https://oeis.org/A038375) records the maximal
number of Hamiltonian paths. The latter currently lists terms through $n=11$;
this does not certify a circulant-only extrapolation. The recurrence and
succession interpretation of [A000255](https://oeis.org/A000255) agree with
$P_n(1,0)=A000255(n-1)$ for $n\ge1$. These are established sequences reused here,
not owner-created entries.

Recovered algorithm routes include
[sequence_speedups_final.py](../../04-computation/sequence_speedups_final.py),
[f_poly_oeis_fast.py](../../04-computation/f_poly_oeis_fast.py), and
[fast_graph_sequences_s169.py](../../04-computation/fast_graph_sequences_s169.py).
They describe fixed-path/H-weighting reductions, subset dynamic programming,
and direct expectation formulas. This inventory does not certify every old
implementation or timing. The present recurrence is proved and checked
independently; it needs neither enumeration of tournaments nor a fitted
recurrence search.

## 2. Inherited exact object and the missing counting measure

The closest proved mechanism is the oriented-run expansion in
[THM-589, even-run slug](../../01-canon/theorems/THM-589-cvH-even-run-closed-form-and-2-over-n-rate.md).
That slug matters because a different theorem has the same legacy ID. For a
permutation $\pi$ of $\{1,\ldots,n\}$, let $a(\pi)$ and $d(\pi)$ count adjacent
entries whose differences are respectively $+1$ and $-1$. Define

$$P_n(u,v)=\sum_{\pi\in S_n}u^{a(\pi)}v^{d(\pi)},\qquad P_0=1. \tag{1}$$

This depends on the reference ordering of the labels. It is not the ordinary
descent polynomial, which would count all decreasing adjacent entries.
The old canonical hostile is MISTAKE-025: $W_8=49752$, not $49760$.
The least-used useful sidecar is the pair of fixed exact run weights $(u,v)$.

For uniformly labelled tournaments with independent fair arcs,

$$\mathbb E H=\frac{n!}{2^{n-1}},\qquad
\mathbb E H^2=\frac{n!W_n}{4^{n-1}},\qquad
\mathrm{CV}(H)^2=\frac{W_n}{n!}-1,\quad W_n=P_n(2,0). \tag{2}$$

To see the second identity, fix the first Hamiltonian ordering. A second
permutation with a reversed shared arc contributes zero; one with $j$ agreeing
shared arcs has joint probability $2^{-2(n-1)+j}$. Summing over the $n!$ first
orderings proves (2). This is not the variance under uniform isomorphism classes
or uniform fixed-path tilings.

## 3. A rational block kernel gives a finite recurrence

Put $a=u-1$, $b=v-1$. Expand the factors at all consecutive-label adjacencies.
A selected connected run of $\ell-1$ edges must be uniformly ascending or
descending: otherwise an interior vertex would have two selected predecessors
or two selected successors in the permutation. Its weight is
$a^{\ell-1}+b^{\ell-1}$ for $\ell\ge2$; an isolated vertex has weight one.
Contracting the runs leaves $k!$ orders of the $k$ blocks. Hence, in formal
power series,

$$B(z)=z+\sum_{\ell\ge2}(a^{\ell-1}+b^{\ell-1})z^\ell
=\frac{z(1-abz^2)}{(1-az)(1-bz)}=\frac{N(z)}{D(z)},$$
$$Y(z):=\sum_{n\ge0}P_n(u,v)z^n=F(B(z)),\qquad
F(t)=\sum_{k\ge0}k!t^k. \tag{3}$$

These are formal ordinary series. No positive radius of convergence for the
factorial series is asserted. Composition is valid because $B(0)=0$.

Direct coefficient comparison gives

$$t^2F'(t)+(t-1)F(t)+1=0.$$

Let $C=N'D-ND'$. Applying the formal chain rule and clearing denominators gives

$$\boxed{DN^2Y'+C(N-D)Y+CD=0.} \tag{4}$$

Here $N=z-abz^3$, $D=1-(a+b)z+abz^2$, and $C(0)=1$. The first coefficient
of $C(N-D)$ is $-1$, while $DN^2$ has zero constant and linear coefficients.
Thus coefficient $z^n$ determines $P_n$ from earlier coefficients. Degrees are
at most $8,7,6$, respectively, so at most seven previous values occur.

**PROVED algorithmic consequence.** For every fixed pair of rational weights
$u,v$, all values through $n=N$ take $O(N)$ rational arithmetic operations.
Only seven previous values are needed by a streaming implementation; the
provided API returns the whole list and stores its output. This is not an
$O(N)$ bit-complexity claim, nor an $O(N)$ computation of every coefficient
of the full bivariate polynomial. Numerators and denominators grow.

For $W_n=P_n(2,0)$, (4) has coefficients

$$\begin{aligned}
L&=z^2+z^4-z^6-z^8,\\
M&=-1+z-3z^2+5z^3+5z^4+3z^5-z^6-z^7,\\
R&=1+3z^2-5z^4+z^6.
\end{aligned}$$

With $W_j=0$ for $j<0$ and $c_n=[z^n]R$, the explicit recurrence is

$$\boxed{\begin{aligned}
W_n={}&nW_{n-1}-3W_{n-2}+(n+2)W_{n-3}+5W_{n-4}\\
     &+(8-n)W_{n-5}-W_{n-6}+(6-n)W_{n-7}+c_n.
\end{aligned}} \tag{5}$$

It supplies $W_0=1$ and is homogeneous for $n\ge7$. No minimal-order assertion
is made. The run identity is inherited; deriving this explicit differential
equation and fixed-width compiler is the present extension. No literature
priority or OEIS acceptance claim is attached to it.

The separate unweighted sequence $P_n(1,0)$ obeys the familiar two-term
succession recurrence. That recurrence does **not** transfer to $W$ merely
because both come from (1): at $n=3$ it predicts $2W_2+W_1=5$, but $W_3=8$.
The authorized narrow repair is recorded in
[THM-589, variance slug](../../01-canon/theorems/THM-589-metagraph-H-variance-closed-form-nud-succession.md)
and the dated 2026-10-14 succession-weight entry in
[MISTAKES](../../01-canon/MISTAKES.md). Other historical claims were not broadly
audited by this packet.

## 4. What transfers to Paley, and what does not

The concept board is: exact sequence definition, source ensemble, run block,
collision partition, recurrence state, arithmetic terminal, and publication
provenance. Equation (3) compiles repeated line blocks in a permutation-overlap
sum. The new [four-wheel packet](wheel_apery_extension_20261014h.md) instead
retains an oriented graph over a marked finite field and reduces its signed
injective character sum to a Legendre moment and, at primes, an Apéry modular
coefficient. Its collision and field-degree guards are indispensable.

The notation $W_n$ in (1) is not the graph $W_4$ or its finite-field wheel
weight. Their transferable operation is **derive a finite kernel with its
source and boundary data, then eliminate that kernel exactly**. This does not
identify their outputs. Already at $n=3$, transitive and cyclic tournaments
have $H=1$ and $H=3$, while the single ensemble number is $W_3=8$. An aggregate
second moment cannot recover an individual Paley count.

For the actual Hamiltonian run expansion, the parent packet derives the
coefficient of the wheel terminal with its ordered placements; the wheel
packet supplies the arithmetic value. These are compatible interfaces. The
historical [H31 design/HYP-2371 prediction](paley_H31_compute_design_kpc1.md)
requires an independently authenticated count of the specific tournament.
Incoming commit `294e754ecd` reports
$H(T_{31})=19876729650849384218015865$ and confirmation of that finite prediction.
That incoming computation is not re-audited here, and is separate from the
fixed-snapshot inventory and the proof of (5). Neither finite prediction
agreement nor a modular residue proves an all-prime asymptotic expansion.

### The annotation $\tau(E)\ge22$ has a different source

In [the Mersenne lattice guard theorem](mersenne_lattice_child_guards_20261014g.md),
$\tau(E)$ is the number of positive divisors of an **odd exponent** $E$.
Grouping prime divisors of $2^E-1$ by $d=\operatorname{ord}_p(2)\mid E$ gives
$S(2^E-1)<(\tau(E)-1)/5$. The additional full joint Giuga/Korselt hypotheses
force $S>4$, hence $\tau(E)\ge22$. These are necessary conditions for that
marked arithmetic class, not an unconditional statement about Mersenne
exponents or Collatz trajectories. The actual residual exponent with four
divisors fails this eligibility test; that does not supply its ROOT word.

The new [order-budget packet](mersenne_order_budget_20261014h.md) strengthens
the same conditional threshold to `tau(E)>=1920` and `omega(E)>=8` by retaining
distinct-prime spacing inside each order group. The earlier 22 remains the
provenance of the annotation, not the best current bound.

Consequently the number 22 does not by itself map to level 11, a twofold
Langlands representation, or a square-packing capacity. A valid transfer would
need a specified source, map, and preserved constraint. The exact successes
here are the run-block compiler and the independently normalized wheel-to-modular
interface, not a numerical identification of their parameters.

## 5. Reproduction and precise scope

```text
python -B -X utf8 04-computation/experiments/oeis_transfer_atlas_20261014h.py
python -B -O -X utf8 04-computation/experiments/oeis_transfer_atlas_20261014h.py
```

The script rejects Boolean/float aliases for exact counts and weights. Its
**583 always-active checks** compare all permutations through $n=8$ at 36
rational weight pairs; an independent composition-power algorithm through
$n=30$ at five weight pairs; all labelled tournaments through $n=5$; factorial
and unweighted-succession boundaries; the historical $W_8$ typo; and malformed
inputs. It also computes $W_{1000}$, an 8,530-bit integer, storing its byte hash
instead of an unnecessary decimal wall. The output includes the complete
reference census and source routes.

Independent peer review confirms the ODE, recurrence, exact API, and arithmetic
complexity scope. The inventory is not represented as an authorship audit. The
current next OEIS obligation is an identified submission/history record or a
fresh per-candidate duplication review; the next Paley obligation is the exact
source-specific run reduction and arithmetic terminal, not another prefix match.

Normal and optimized output agree exactly. The complete saved JSON has LF
SHA256 `0203e84401fff52ac5026b0b9889e7468df762745091abd04284d0d8dbd1baf0`.
The final proof/scope note received an independent peer PASS after the core replay.
