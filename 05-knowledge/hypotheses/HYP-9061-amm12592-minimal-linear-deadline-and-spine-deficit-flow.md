---
id: HYP-9061
title: "AMM 12592 minimal linear deadline: spine deficit flow and the two-bias rate gate"
status: >
  OPEN (frontier question named; certificate decoded and verified exact;
  context CONFIRMED 2026-07-30: fragment re-supplied with the AMM 12592
  problem statement and minimal-C framing adjacent, margin F byte-identical
  to the repo reconstruction — see
  05-knowledge/results/amm12592-snippet-context-confirmation-deathstar-coinC2.md;
  rival homes HYP-9023/klein-S402-S404/mac-mini-S168 demoted on provenance)
source: death-star-2026-07-30-coinC
related:
  - THM-2160
  - THM-2225
  - THM-2976  # binary-clock parity: checkpoint vanishing, o* = 2^(v2(M+1)), odd-unit-fraction corner ladder
external:
  - "Elliot Glazer, American Mathematical Monthly Problem 12592 (2026)."
script: 04-computation/amm12592_artanh_certificate_decode_deathstar.py
output: 05-knowledge/results/amm12592_artanh_certificate_decode_deathstar.out
---

# HYP-9061 -- the minimal C in the critical-run deadline T(n) <= Cn + D

> **MAJOR UPDATE 2026-08-03 (boxeph) — the general-class golden floor is
> DEMOTED; the bracket is re-opened from below.** A hostile audit
> (MISTAKE-361; `05-knowledge/results/amm12592-golden-floor-audit-boxeph.md`)
> refutes THM-3024's promotion `C*_general = log_5(5 phi^2)` *within its own
> transportation model*: with forward routing at preserved absolute degree
> and an unbounded window, every tail cut with a deeper shell available is
> satisfied for ANY `gamma > 0` (exact, independently confirmed twice), so
> the model yields no general floor at all and the reported binding cuts
> were truncation-edge artifacts. Current honest state:
> `C = log_5(5 phi^2)` attained for ALL n <= 255 exactly and for all
> n <= 2047 with additive slack <= 16 (THM-3029/3302/3329/3330; exact slack
> thresholds D0*(R) = 0,1,5,15,38,89,192,[401..416] for R = 128..16384 —
> LINEAR growth for the plain rule, eps_inf ~ 0.025, so the golden route
> needs the bulk alternation rule; Y-box/Lucas-diagonal structure proved);
> balanced-block-class
> floor `C*_block > 1.5970` exact (THM-3009, audit-hardened, ladder to
> m = 4096, certified rational bracket
> `115939/193882 < gamma* < 105183/175895`); GENERAL-class floor OPEN —
> even `C* = 3/2` is not currently excluded for non-block rules. The repair
> needs a deadline-bounded routing window derived from the extractor axioms:
> the pathwise deadline is exactly what the transportation relaxation
> forgot. Note the symmetry of failure: the construction side's obstruction
> is combinatorial redistribution (THM-3026), the floor side's gap is
> unbounded redistribution — the SAME missing ingredient seen from both
> sides, cf. the integrality-gap framing in
> `05-knowledge/results/amm12592-epoch-closure-nonnegative-transportation-form-boxeph.md`.

> **MAJOR UPDATE 2026-07-31 (opus, THM-3007/3006) — "every known rule has
> `C=2`" is OBSOLETE.** THM-3007: a composition-balanced block `[N,N+l)`
> exists iff `N` and `N+l` are both powers of two, so the dyadic shell
> structure is FORCED and no balanced block has ratio below 2 — this closes
> section 2b's "ratio `rho<2` shell" route outright. But THM-3008 shows that
> block ratio and DEADLINE SLOPE are different: solving the within-block
> normal form gives exactly fair rules with within-shell ratio
> `rho(4)=3/2`, `rho(8)=14/9`, `rho(16)=25/16`, `rho(32)<=11/7`, hence
> `T(n) <= (11/7)n` for `4<=n<64`. The slope-2 envelope of every classical
> construction is an artifact of never optimizing the within-shell profile.
> The live question is now UNIFORMITY: is `sup_r rho(2^r) < 2`? Data
> `1.5000, 1.5556, 1.5625, 1.5714` increases slowly. Sections 2c–2e below
> (the two-bias rate gate, the capacity desert, the lower-bound reading of
> (27)) were premised on `C*=2`-adjacent behaviour and must be re-read in
> this light; the `C* >= 9049/6592 ~ 1.3727` lower-bound reading is NOT
> contradicted by `rho ~ 1.57`, and the two now bracket a plausible answer.

THM-2225 proves the AMM 12592 envelopes `2n` and `max(2,2n-1)`; THM-2160
section 5 reaches `T(n) <= 2n-2` for `n>=3`, and THM-3032 sharpens that to a
rule dominating the checksum pointwise. Every *classical* rule has `C=2`
(superseded for the general question, see the update above).
THM-2160 explicitly does not assert its shell lower bound for cross-shell
rules. The frontier question:

```text
C* = inf{ C : exists D and an exactly fair extractor with
            pathwise deadline T <= C n + D on critical value n }.   (Q)
```

## 1. Spine reformulation (this session; elementary, checkable)

Assign to each undecided finite word `x` its future heads share
`G_x(p) in [0,1]`, with `G_x = p G_{x0} + q G_{x1}`, `G_root = 1/2`, and
`G_x` identically `0` or `1` on decided words. Every word extending
`0^m 1` has critical value exactly `m`, so the whole subtree there is
decided within relative depth

```text
d_m = (C-1) m + D - 1 =: gamma*m + D - 1.
```

Telescoping both constant spines gives the exact equivalence: a deadline-
`(Cn+D)` fair extractor exists iff there are decided-tree polynomials
(integer Bernstein coefficients `0 <= w_{m,k} <= binom(d_m,k)`, values in
`[0,1]` on `(0,1)`) `W_m, V_m` of degree `<= d_m` with

```text
sum_{m>=1} p^m q W_m(p) + sum_{m>=1} q^m p V_m(p) = 1/2   for all p.  (S)
```

`C* = 1 + gamma*` where `gamma*` is the minimal degree-growth rate in (S).
Writing `W_m = 1/2 + Delta_m`, `V_m = 1/2 + Delta'_m`, (S) is the deficit
flow `sum p^m q Delta_m + sum q^m p Delta'_m = 0`: the half-integer parity
deficits of odd Hamming classes (Lucas: forced whenever a shell length is
not dyadic) must be carried to later spine positions and repaid within
budget `binom(d_m, k)`. The carry chain converges only under exponential
rate conditions comparing likelihood ratios `log(p/q)` across the biases
where deficits concentrate: a two-bias gate.

## 2. The decoded certificate

An externally supplied fragment (proposer-side notes, equation (27) in its
own numbering) was decoded and refereed exactly in the companion script:
with rapidity bounds `L(t) <= log((1+t)/(1-t)) <= U(t)`,

```text
t_A = 389/2181,          p_A = 1285/2181,   p_A/q_A = 1285/896,
t_B = 5872957/11821757,  p_B = 8847357/11821757,
                         p_B/q_B = 8847357/2974400,

(27):  (2457/6592) log(8847357/2974400) - log(1285/896) > 1/25,
```

certified by `(2457/6592) L(t_B) - U(t_A) = 1/25 + F` with an explicit
positive rational `F` (margin ~ 0.00474; true slack ~ 0.00573). Note
`2457 = 3^3*7*13`, `6592 = 2^6*103`, `2457/6592 ~ 0.37272`. The gate
shape matches the carry-chain balance of section 1; whether it certifies
an upper construction (`C = 1 + 2457/6592 = 9049/6592 ~ 1.3727`?) or a
dual step in a `C >= c_0` lower bound is exactly what (Q) must settle.

## 2b. The corner-deficit carry mechanism (death-star, same session)

Sharper structure, elementary to check:

1. **Corner reduction.** A shell `[n_j, n_j + l_j)` (0-side, all words decided
   at the shell end) with dyadic tail `l_j = 2^(a_j)` has every interior tail
   Hamming class of even size (Lucas); the only forced parity failure is the
   corner word `0^(n_j) 1^(l_j)`, which exports one deficit `+-(1/2)
   p^(n_j) q^(l_j)`. Mirror for 1-side shells.
2. **Why dyadic ratio 2 is magic.** With `l_j = n_j` (ratio 2) the 0-corner
   and 1-corner land on the same monomial `(h,h)` and annihilate directly:
   this is exactly THM-2160's middle-pair `0^h 1^h / 1^h 0^h` trick. For any
   ratio `rho < 2` the two corners are distinct monomials and can never meet
   at their birth cells.
3. **Carry dynamics.** Since `(1+u)^N = (1+u)^(N-Delta) (1+u^Delta)` over
   `F_2` for dyadic `Delta = 2^a`, pushing a half-deficit at `(z,o)` down by
   `Delta` (multiplying by `(p+q)^Delta = 1`) splits it into half-deficits at
   `(z+Delta, o)` and `(z, o+Delta)` only, with all interior binomial pieces
   `binom(Delta,i)/2` absorbed integrally -- provided the interior pieces fit
   in the integer budgets `binom(d_m, k)` of the receiving spine cells, which
   all live in the cone `o <= gamma z + O(D)`. Coordinates only increase.
4. **Reduction of (Q).** A `C = 1+gamma` scheme exists iff the countable set
   of forced corner deficits admits a pairwise annihilation routing by dyadic
   split-jumps whose interior absorption respects the cone budgets. The
   binding constraint is an exponential-rate comparison (entropy of the jump
   binomial vs budget entropy along the two critical lattice rays), i.e.
   precisely a two-bias log-likelihood gate of the shape (27). Whether the
   gate opens (construction, `C* = 9049/6592`-ish) or closes (lower bound)
   at `gamma = 2457/6592` is the live question.

## 2b0. THE CONSTRUCTION SIDE OPENS (2026-07-31, THM-3002 + lane G5)

`C* = 2` is no longer the favoured answer. Lane G5's checkpoint-closure
reduction (closing the books exactly at every dyadic checkpoint
`M = 2^r - 1` is *sufficient* for `C* <= 1 + gamma`) plus THM-3002 give:

1. **Every dyadic epoch through `[16,31]` closes at `gamma = 1/2`**
   (VERIFIED-EXACT, witnesses re-derived independently): an exactly fair
   extractor exists for all critical values `n <= 31` with
   `T(n) = n + 1 + floor(n/2)`.
2. **Normal form:** block closure forces `F = q^{m_lo-1}H`,
   `G = -p^{m_lo-1}H`; `H = 1` says the epoch's entire imbalance is the
   single middle pair `0^R 1^R` vs `1^R 0^R` — THM-2160's trick promoted
   from one row to a whole epoch — and decouples the two sides.
3. **Sharp capacity criterion:** `max |[p^t]Delta| = binom(d,t)2^t` over the
   Lucas box, so `sum_{i<=t} binom(d_i,t-i)2^{t-i} >= binom(R-1,t)` is
   necessary. Exact ledger: exponentially deficient for `gamma < 1/2`,
   marginal-then-deficient at `gamma = 1/2` (dead by `R = 64`), uniformly
   ample for `gamma >= 3/5`. Asymptotic threshold `gamma* ~ 0.598` from a
   **two-ray entropy comparison** — structurally the shape of (27).
4. So the `gamma = 1/2` successes are finite-size, and the live target is
   `gamma ~ 3/5` (`C = 8/5`): (*) is already solved there for `R = 8, 16`.

Falsifier for `C* = 2`: closure of all epochs at any `gamma < 1`.
Falsifier for the program: an `R` at which (*) is exhaustively infeasible
at every `gamma < 1` and every `H`.

## 2c0. Direction status after the coinC2 session (2026-07-30, evening)

The **evaluation reading of (27) is CLOSED as a class**: THM-2977 (the
evaluation wall) + lane G2 (`laneG2-padic-phi-findings.md`) prove that
every denominator-clearing functional at engineered rational biases has
choice-invariance modulus bounded independently of `M` (K = 6 at the
certificate pair, exactly), that its forced content collapses to one bit
of boundary bookkeeping (blind to the band), and that the proved envelope
covers every residue class from `M = 1` on. Hence (27) can only be the
**numeric gate of a rate/entropy (tropical) dual**: a two-ray
capacity-vs-forced-mass comparison in which the rapidities
`r_A, r_B = log(p/q)` are Legendre-dual ray slopes, `alpha = 2457/6592` a
ray weight, and `1/25` a per-`M` rate margin absorbing `O(D0)`. The
2-adic bias engineering survives as **checksum alignment** (the minimal
coupled one-bit invariant exists iff `s_A != s_B` aligns; SPECULATION as
to intent). Within that reading the LOWER-bound direction
(`C* >= 9049/6592`) remains favored over the construction gate: the
construction side lost its two best mechanisms this session (THM-2976's
corner clock does not rescue transport — greedy at the corner-clocked
rates `1/3, 1/5` freezes harder than at `1/2`, laneD sec. 5b — and the
deep-corner clock fires infinitely often ONLY at `gamma = 1/(2^j - 1)`,
excluding `alpha`, lane G2 T4). Sections 2c/2e below are kept as the
historical route to this conclusion.

## 2c. Working reading of the certificate (dual side; UNVERIFIED direction — evaluation variants now REFUTED, see 2c0)

Truncate the deficit identity `sum p^m q Delta_m + sum q^m p Delta'_m = 0`
at total degree `A`: the shallow part re-expanded at level `A` has integer
(after doubling) coefficients `C^(A)_O = sum_(z+o<=A) c_(z,o)
binom(A-z-o, O-o)`, and the identity forces `C^(A)_O = -(level-A data of
the deep part)`, whose magnitude is bounded by budgeted cone masses beyond
level `A`. The mod-2 cone pairing (`lambda_(A,B) = binom(A-z-o, B-o)`,
Pascal-harmonic, conserved by dyadic split-jumps) can force `|C^(A)_O| >=
1/2` at cells where the parity count is odd. Infeasibility at slope
`gamma` follows wherever forced-half beats the deep-mass exponential rate:
a two-ray tropical comparison, certified at two rational biases -- the
exact shape of (27). Under this reading (27) is a **lower-bound dual step
and C* >= 1 + 2457/6592 = 9049/6592 ~ 1.3727** would be the certified
consequence; the matching upper construction is then the open half.
This direction assignment is a working hypothesis, not verified: the
opposite (construction-gate) reading is not yet excluded.

## 2d. The capacity desert and the Legendre reading of the two biases

At depth rate `gamma < 1`, integer absorption capacity exists only in the
two cones `o <= gamma z + O(D)` (0-side cells `(m + e - k, 1 + k)`) and
its mirror `o >= z/gamma - O(D)`. The band between them is a **desert**:
no W/V cell lives there, so deficit mass crossing it can only split
(Pascal), never settle. Pushing a deficit from level `z+o` to level `L`
spreads it binomially with bulk slope tending to the diagonal `o = z` --
inside the desert -- so only exponentially small tails reach the cones,
and the per-cell requirement `binom(Delta, j)/2 <= binom(d_m, k)` becomes
an entropy race `Delta H(j/Delta)` vs `gamma m H(k/(gamma m))`, tightest
along one or two critical rays. Legendre duality identifies the decoded
certificate's biases with those rays: a binding ray with ones-fraction
`x` has dual bias `q = x` via `H'(x) = log((1-x)/x)`, and indeed
`q_A/p_A = 896/1285 ~ 0.697`, `q_B/p_B = 2974400/8847357 ~ 0.336` --
one ray inside the 0-cone edge (`0.336` vs edge `2457/6592 ~ 0.3727`),
one in the desert. Under this reading (27) certifies the entropy race at
`gamma = 2457/6592` with margin `1/25`, i.e. plausibly the lower bound
`C* >= 9049/6592`, en route to the transcendental threshold where the
desert max-flow balances.
**RETRACTED as evidence (klein-S428, 2026-07-31).** The "capacity
straddle" does not support the specific weight: (27) admits an open
half-line `alpha > (r_A + 1/25)/r_B = 0.3674729...` (e.g. `3/8`
certifies with the *larger* margin `7.21e-3`), and the straddle window has
width `0.33`, containing `3/8, 2/5, 37/100, 7/19, 41/110` alike. So
`2457/6592` is an **output of the source's construction**, not a
consequence of the inequality; and since `r_A/r_B` is irrational (klein's
isolated-prime argument, `257` vs `2949119`), a nonzero floor is free and
only its *size* is open — `1/25` is a chosen safety margin, not an
extracted one. Operational corollary: **do not invert (27) to recover the
construction**; derive the rate from the ledger (THM-3002) and use (27)
only as a verification shape. The single-jump lemma
(`amm12592_single_jump_routing_slack_deathstar.py`: slack exactly
`D = l`, envelope degrades to `2n`) is the finite shadow of the desert:
naive routing pays back the whole shell gain, and `C = 2` is its fixed
point.

## 2e. The 2-adic engineering of the certificate biases (decode advance)

Both certificate biases have odd denominator, odd `p`-numerator, and
`q`-numerator divisible by a high power of two:

```text
q_A = 896/2181,      896 = 2^7 * 7          (s_A = 7),
q_B = 2974400/11821757,  2974400 = 2^6*5^2*11*13^2   (s_B = 6).
```

At such a bias `p = a/b`, `N_M := 2 D_M(p) b^{A_M} = sum 2delta a^z
(b-a)^o` is an integer with per-term 2-adic valuation `v_2(2delta) + s o`:
the minimal forced-odd band position is a non-archimedean leading term, so
a single evaluation can extract one forced half-odd band coefficient,
bypassing the archimedean sign-cancellation that defeats parity-free
tests. The envelope bounds `|N_M| <= (p^{M+1}+q^{M+1}) b^{A_M}`; the
nonzero-integer-below-one contradiction fails at any single bias (rate
check: `gamma log b + log a > 0` always), so the genuine dual must couple
the two biases -- cancel the dominant archimedean mass between
`N^{(A)}, N^{(B)}` while keeping the 2-adic leading terms misaligned
(`s_A = 7` vs `s_B = 6`). The resulting rate inequality has exactly the
shape of (27); reconstruction assigned to wave-2 lane G.

## 3. Cheapest decisive tests

1. Bounded depth (`gamma = 0`, i.e. `T(n) <= n + D`): (S) forces
   `sum_m p^m q W_m` with `W_m` in a finite set; Carlson--Szego rigidity
   (power series with finitely many distinct coefficients analytic past
   the unit circle is rational) should force eventual periodicity and a
   finite refutation. Expected: `C = 1` impossible for every `D`.
2. Extract `Delta_m` for THM-2225's checksum rule: locate where its
   deficit flow uses degree growth `m`, and which part is slack.
3. Exact truncated feasibility (length <= ~24) for envelopes
   `T(4) = 5`, `T(5) = 6`, `ceil(3n/2)+D`: infeasibility certificates via
   rational LP duality are rigorous; feasible prefixes guide construction.

Falsifier for `C* = 2`: any feasible (S) family with `d_m = gamma*m + O(1)`,
`gamma < 1`. Falsifier for `C* < 2`: a dual functional bounding every
truncation away from `1/2` under `gamma < 1` growth.
