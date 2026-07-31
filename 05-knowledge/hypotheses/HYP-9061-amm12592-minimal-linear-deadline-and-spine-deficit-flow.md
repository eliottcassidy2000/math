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

THM-2225 proves the AMM 12592 envelopes `2n` and `max(2,2n-1)`; THM-2160
section 5 reaches `T(n) <= 2n-2` for `n>=3`. Every known rule has `C=2`.
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

## 2c. Working reading of the certificate (dual side; UNVERIFIED direction)

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
desert max-flow balances. The single-jump lemma
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
