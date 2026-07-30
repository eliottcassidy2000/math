---
id: HYP-9061
title: "AMM 12592 minimal linear deadline: spine deficit flow and the two-bias rate gate"
status: OPEN (frontier question named; certificate decoded and verified exact)
source: death-star-2026-07-30-coinC
related:
  - THM-2160
  - THM-2225
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
