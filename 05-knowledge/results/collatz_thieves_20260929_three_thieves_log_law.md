# Three thieves and the planar walk: the probability that a random two-colour necklace admits a fair three-way consecutive split is Θ(1/log m)

**Status: PROVED (Theorem 1, upper and lower bounds of order `1/log m`; the local-limit lemma for hypergeometric prefix counts is derived from the binomial local limit theorem with a Stirling-type remainder, CITED as standard) + FINITE-EXACT / EMPIRICAL for the constants (`P(N_3 >= 1) log m = 2.13 .. 2.25` at `m = 10^3 .. 10^5`, `E[N_3^2] ~ 0.93 log m`). Single session (collatz-necklace-20260929, part 2), not independently audited.**

This closes the open item of [THM-4515](../../01-canon/theorems/THM-4515-fair-consecutive-splits-of-cycle-necklaces-are-circulant.md) section 5: the two-sided order of the three-thieves probability. Scripts: `04-computation/experiments/collatz_necklace_20260929_fairsplit_asymptotics.py`, `..._fairsplit_loglaw.py` (the exact second moment and the Monte Carlo); no new computation is needed for the proof.

## 0. Setting and statement

Fix `rho in (0,1)`. For `m >= 1` let `x = x(m)` with `x/m -> rho` (any such sequence; `x = round(rho m)` say), `K = 3m`, `X = 3x`, and let `w` be a uniformly random binary word of length `K` with `X` ones (the shape of a cycle necklace with `K` steps and `X` odd steps, `gcd(K,X)` divisible by three). Cut positions are `r in Z/m` (a cut at `r` and at `r+m` give the same three arcs); the cut at `r` is *fair* if each of the three windows `[r, r+m)`, `[r+m, r+2m)`, `[r+2m, r+3m)` (indices mod `K`) contains exactly `x` ones. Let `N = N_3(w)` be the number of fair cuts in `[0, m)`.

**Theorem 1.** There are constants `0 < c_1 <= c_2 < infinity` depending only on `rho` (and `m_0(rho)`) such that for all `m >= m_0`

`c_1 / log m <= P(N >= 1) <= c_2 / log m`.

Moreover `E[N] = m binom(m,x)^3 / binom(3m,3x) -> sqrt(3)/(2 pi rho(1-rho))` (THM-4515 section 5), and the exact second moment `E[N^2] = (m/binom(3m,3x)) sum_{d=0}^{m-1} sum_t binom(d,t)^3 binom(m-d,x-t)^3` satisfies `E[N^2] = Theta(log m)`.

The lower bound is the second-moment bound of THM-4515 made asymptotic; the new content is the upper bound, whose mechanism is the return structure of a planar lattice walk: once a fair cut exists, the walk of the window imbalances is at the origin, and a planar walk that has reached the origin returns to it about `c log m` further times, so the expected number of fair cuts (which is bounded) is spread over events of probability `~ 1/log m` carrying `~ log m` cuts each.

## 1. The three prefix processes and the imbalance walk

Number the positions `0, ..., K-1` and call `W_i = [(i-1)m, im)`, `i = 1, 2, 3`, the three windows of the cut at `0`; let `x_i` be the number of ones in `W_i` (`x_1 + x_2 + x_3 = X`) and, for `0 <= u <= m`, let `a_i(u)` be the number of ones among the first `u` positions of `W_i`. For `0 <= r < m` the window `[r, r+m)` contains `(x_1 - a_1(r)) + a_2(r)` ones, `[r+m, r+2m)` contains `(x_2 - a_2(r)) + a_3(r)`, and `[r+2m, r+3m)` contains `(x_3 - a_3(r)) + a_1(r)`. Hence the cut at `r` is fair iff

`a_2(r) - a_1(r) = x - x_1` and `a_3(r) - a_2(r) = x - x_2`; (1.1)

the third window is then automatic. Put `z = (x - x_1, x - x_2) in Z^2` and `W(r) = (a_2(r) - a_1(r), a_3(r) - a_2(r))`. Then `N = #{r in [0,m) : W(r) = z}`.

**Lemma 1.1 (conditional structure).** Conditionally on `(x_1, x_2, x_3)`, the three processes `a_1, a_2, a_3` are independent, and `a_i` is the prefix-count process of a uniformly random arrangement of `x_i` ones in `m` slots; in particular `a_i(u)` is hypergeometric `H(m, x_i, u)` with mean `u x_i/m` and variance `v_i(u) = u(m-u) x_i (m-x_i) / (m^2 (m-1))`. Conditionally on `(x_1,x_2,x_3)` and on the prefixes `(a_i(u'))_{u' <= u}`, the suffix processes `u' -> a_i(u+u') - a_i(u)` are independent prefix-count processes of uniformly random arrangements of `x_i - a_i(u)` ones in `m - u` slots.

*Proof.* A uniformly random word of shape `(K,X)` restricted to the three windows, given the window counts, is a product of three uniform arrangements; and a uniform arrangement of `y` ones in `L` slots, given its first `u` entries, is a uniform arrangement of the remaining ones in the remaining slots. ∎

So `W` is the difference walk of three independent hypergeometric prefix processes; its increments are `(e_2 - e_1, e_3 - e_2)` with `e_i in {0,1}` the next letters. Given a fair cut at `u` (so `W(u) = z`), the event that `u + u'` is also fair reads `b_1(u') = b_2(u') = b_3(u')` where `b_i` are the three suffix prefix-counts; the target `z` has been replaced by the origin, and the suffix counts `y_i = x_i - a_i(u)` satisfy, by (1.1), `y_2 - y_1 = x_2 - x` and `y_3 - y_2 = x_3 - x`.

## 2. Two elementary tools

**Lemma 2.1 (Hoeffding for sampling without replacement).** For a hypergeometric variable `H(L, y, u)` with mean `mu`, `P(|H - mu| >= t) <= 2 exp(-2 t^2/u)`; the same bound holds for the maximum over `u <= L` of `|a(u) - u y/L|` with an extra factor `L` (union bound). (Hoeffding 1963, section 6; sampling without replacement is dominated by sampling with replacement in the convex order.)

Define the *good event* `G = G_B` by: `|x_i - x| <= B sqrt(m log m)` for `i = 1, 2, 3`, and `|a_i(u) - u x_i/m| <= B sqrt(m log m)` for all `i` and all `0 <= u <= m`. By Lemma 2.1 with `t = B sqrt(m log m)`, `P(G^c) <= 6 m^{-2B^2} + 6 m cdot 2 m^{-2B^2} <= 18 m^{1 - 2B^2}`; with `B = 1` this is `18/m`.

**Lemma 2.2 (local lower bound for hypergeometric prefix counts).** For every `rho_0 in (0, 1/2)` there are `u_0` and `c_0 > 0` such that: if `L >= 1`, `y/L in [rho_0, 1 - rho_0]`, `u_0 <= u <= L/2`, and `mu = u y/L`, `v = u (L-u) y (L-y)/(L^2 (L-1))`, then `P(H(L,y,u) = t) >= c_0 / sqrt(u)` for every integer `t` with `|t - mu| <= 2 sqrt(v)`.

*Proof.* Write `p = y/L`. For independent `B_1 ~ Bin(u, p)`, `B_2 ~ Bin(L-u, p)` the conditional law of `B_1` given `B_1 + B_2 = y` is `H(L, y, u)` (whatever `p`), so `P(H = t) = P(B_1 = t) P(B_2 = y - t) / P(B_1 + B_2 = y)`. The binomial local limit theorem with a uniform remainder (e.g. the Stirling/Robbins form: for `k = np + d` with `|d| <= 2 sqrt(np(1-p))` and `np(1-p) >= n_0`, `P(Bin(n,p) = k) = (2 pi n p (1-p))^{-1/2} exp(-d^2/(2np(1-p))) (1 + eta)` with `|eta| <= C/sqrt(np(1-p))`) applied three times — to `B_1` at `t = mu + d` (so `d = t - mu`), to `B_2` at `y - t = (L-u)p - d`, and to `B_1 + B_2 ~ Bin(L, p)` at its mean `y = Lp` — gives, for `u_0` large enough that all three remainders are at most `1/10`,

`P(H = t) >= (8/10)(8/10)/(11/10) cdot sqrt(L / (2 pi u (L-u) p(1-p))) exp(-d^2 (1/(2up(1-p)) + 1/(2(L-u)p(1-p))))`
`= c' (2 pi v')^{-1/2} exp(-d^2/(2 v'))`, `v' = u(L-u)p(1-p)/L`,

and `v' = v (L-1)/L in [v/2, v]`. For `|d| <= 2 sqrt(v)` the exponent satisfies `d^2/(2v') <= 4v/(2 cdot v/2) = 4`, and `v' <= u p(1-p) <= u/4` gives `(2 pi v')^{-1/2} >= (pi u/2)^{-1/2}`. Hence `P(H = t) >= c_0/sqrt(u)` with `c_0 = c' e^{-4} sqrt(2/pi)`. ∎

(The constant is not optimised; only the order `u^{-1/2}` uniformly on the window `|t - mu| <= 2 sqrt(v)` is used.)

**Lemma 2.3 (three-fold coincidence).** Let `b_1, b_2, b_3` be independent prefix-count processes of uniform arrangements of `y_1, y_2, y_3` ones in `L` slots, with `L >= m/2`, `y_i/L in [rho_0, 1-rho_0]`, and `|y_i - y_j| <= Delta := 3 sqrt(m log m)`. Put `C_1 = 576/(rho_0(1-rho_0))`. Then for `u_0 <= u <= m/(C_1 log m)`,

`P(b_1(u) = b_2(u) = b_3(u)) >= c_3 / u`, `c_3 = c_3(rho_0) > 0`.

*Proof.* The three variables `b_i(u)` are hypergeometric `H(L, y_i, u)` with means `mu_i = u y_i/L` and variances `v_i in [u rho_0(1-rho_0)/2, u/4]` (using `u <= L/2`); put `v_min = u rho_0(1-rho_0)/2`. The means differ by `|mu_i - mu_j| = u|y_i - y_j|/L <= 2 u Delta/m = 6 sqrt(u) sqrt(u log m/m) <= 6 sqrt(u)/sqrt(C_1) = sqrt(u rho_0(1-rho_0))/4 <= sqrt(v_min)/2`. So every integer `t` with `|t - mu_1| <= sqrt(v_min)` satisfies `|t - mu_i| <= (3/2) sqrt(v_min) <= 2 sqrt(v_i)` for all `i`, and Lemma 2.2 gives `P(b_i(u) = t) >= c_0/sqrt(u)` for each `i`. Summing over the at least `2 sqrt(v_min) - 1` such `t`: `sum_t prod_i P(b_i(u) = t) >= (2 sqrt(v_min) - 1) c_0^3 u^{-3/2} >= c_3/u` for `u >= u_0`. ∎

## 3. Proof of the upper bound

Let `tau = min{r in [0,m) : cut at r fair}` (`tau = infinity` if `N = 0`), and let `F_u` be the sigma-field generated by `(x_1,x_2,x_3)` and the prefixes `(a_i(u'))_{u' <= u, i = 1,2,3}`. The event `{tau = u}` is `F_u`-measurable (fairness at `r` is determined by the prefixes up to `r`, by (1.1)).

**Step 1 (first cut in the first half).** Since `E[N] <= E_max := sup_m m binom(m,x)^3/binom(3m,3x) < infinity` (THM-4515 section 5; `E_max <= 1.2` for `rho = 1/2`),

`E_max >= E[N] >= sum_{u < m/2} P(tau = u) E[N | tau = u] >= sum_{u < m/2} E[ 1_{tau = u} (1 + N'_u) ]`,

where `N'_u = #{r in (u, u + m/2) : cut at r fair}`; note `u + m/2 < m`, so these cuts lie in `[0, m)`. By Lemma 1.1 and (1.1), conditionally on `F_u` and on `{tau = u}`, the cuts `r = u + u'` with `1 <= u' < m/2` are fair exactly when `b_1(u') = b_2(u') = b_3(u')` for the three independent suffix processes with parameters `(L, y_i) = (m - u, x_i - a_i(u))`. On the good event `G` (Lemma 2.1, `B = 1`), and for `u < m/2`, these parameters satisfy `L >= m/2`, `|y_i - rho L| <= 3 sqrt(m log m)`, hence `y_i/L in [rho_0, 1 - rho_0]` with `rho_0 = min(rho,1-rho)/2` for `m >= m_0`, and `|y_i - y_j| = |x_i - x_j| <= 2 sqrt(m log m)`. Lemma 2.3 (with `rho_0 = min(rho, 1-rho)/2`) then gives

`E[N'_u | F_u] >= sum_{u' = u_0}^{m/(C_1 log m)} c_3/u' >= c_3 (log m - log log m - log C_1 - log u_0 - 1) >= (c_3/2) log m` on `G cap {tau = u}`, `m >= m_0`.

Therefore `E_max >= sum_{u < m/2} P(tau = u, G) (c_3/2) log m`, i.e.

`P(tau < m/2, G) <= 2 E_max / (c_3 log m)`. (3.1)

**Step 2 (rotation).** The map that rotates the word by `m/2` positions (for `m` even; by `floor(m/2)` in general) is a bijection of the words of shape `(K,X)` and sends the cut at `r` to the cut at `r - m/2 mod m`. Hence `P(there is a fair cut in [m/2, m)) = P(there is a fair cut in [0, m/2)) = P(tau < m/2)`, and

`P(N >= 1) <= P(tau < m/2) + P(fair cut in [m/2, m)) = 2 P(tau < m/2) <= 2 P(tau < m/2, G) + 2 P(G^c) <= 4 E_max/(c_3 log m) + 36/m`.

This is the upper bound with `c_2 = 4 E_max/c_3 + 1`. ∎

## 4. The lower bound and the second moment

`P(N >= 1) >= E[N]^2/E[N^2]` (Cauchy–Schwarz). The exact second moment: two fair cuts at `0` and at `d in [1, m)` force the six window counts to be `t, x - t, t, x - t, t, x - t` (the six intervals `[0,d), [d,m), [m, m+d), [m+d, 2m), [2m, 2m+d), [2m+d, 3m)` have counts `a_1..a_6` with `a_1 + a_2 = a_3 + a_4 = a_5 + a_6 = x` from the cut at `0` and `a_2 + a_3 = a_4 + a_5 = a_6 + a_1 = x` from the cut at `d`), so `#{w : both fair} = sum_t binom(d,t)^3 binom(m-d,x-t)^3`, and by rotation symmetry

`E[N^2] = m sum_{d=0}^{m-1} P(0 and d fair) = (m/binom(3m,3x)) sum_{d=0}^{m-1} sum_t binom(d,t)^3 binom(m-d, x-t)^3`.

The term `d = 0` is `E[N]`. For `1 <= d <= m-1` the summand is, by Lemma 1.1 with `u = d` and the three independent hypergeometrics `H(m, x, d)` (here `L = m`, `y_i = x`), `P(a_1(d) = a_2(d) = a_3(d)) = sum_t p_d(t)^3` with `p_d = H(m,x,d)`; the local limit theorem in the form of Lemma 2.2 (two-sided, with the matching upper bound `p_d(t) <= C/sqrt(v_d)` from the same three binomial estimates) gives `sum_t p_d(t)^3 = (1 + o(1)) / (2 pi sqrt(3) v_d)` uniformly for `d` and `m - d` large, `v_d = d(m-d) x(m-x)/(m^2(m-1))`, hence

`E[N^2] = E[N] + (1 + o(1)) (m / (2 pi sqrt 3 rho(1-rho))) sum_{d} 1/(d(m-d)) = E[N] + (1 + o(1)) (log m)/(pi sqrt 3 rho(1-rho)) + O(1)`,

i.e. `E[N^2] ~ (log m)/(pi sqrt 3 rho(1-rho))` (`= 0.735 log m` at `rho = 1/2`; the computed exact values `E[N^2]/log m = 1.05, 1.01, 0.98, 0.95, 0.93` at `m = 10^2 .. 10^4` approach it from above). Then `P(N >= 1) >= E[N]^2/E[N^2] >= c_1/log m`. ∎

**Remark (size-biasing).** `E[N^2]/E[N] = E[N | fair at 0]` is the *size-biased* conditional mean (`~ 0.735 log m`), while `E[N | N >= 1] = E[N]/P(N >= 1)` is the ordinary one (`~ 0.5 log m` by the Monte Carlo); the upper bound of section 3 is exactly a lower bound on the ordinary conditional mean, which the size-biased quantity does not give.

## 5. The constants (EMPIRICAL) and the dimension count

Monte Carlo at `rho = 1/2` (`..._fairsplit_loglaw.out`): `P(N >= 1) log m = 2.13, 2.17, 2.22, 2.19, 2.25` at `m = 10^3, 3 10^3, 10^4, 3 10^4, 10^5`, and `E[N | N >= 1] = 3.5, 4.2, 4.4, 5.2, 5.4 ~ 0.47 log m`; so `P(N >= 1) ~ 2.2/log m` and `E[N] = P(N >= 1) E[N | N >= 1] ~ 2.2 cdot 0.5 = 1.1 = sqrt(3)/(2 pi rho(1-rho))`, consistent. The constant `2.2` is not identified; the proof gives `c_1 = 1.2^2/0.735 = 1.5` asymptotically for the lower bound and an unoptimised `c_2`.

The dimension count behind the trichotomy of THM-4515 section 5: for `j` thieves the imbalance walk lives in `Z^{j-1}`; `j = 2` is a one-dimensional bridge with a forced zero (the discrete intermediate value theorem), `j = 3` is the planar, recurrent-but-marginal case treated here (`E[N] -> const`, `P(N >= 1) ~ c/log m`), and `j >= 4` is transient (`E[N_j] -> 0` like `m^{-(j-3)/2}`, so `P(N_j >= 1) -> 0` by Markov). The Collatz reading (THM-4515): a hypothetical `3x+1` cycle whose shape has `gcd(K,X) = 3` admits the Eisenstein-type factorisation of its clock through a fair 3-split with probability about `2.2/log(K/3)` under the uniform-necklace model; with `gcd >= 4` essentially never.

## 6. Boundary

Not claimed: the exact constant `c` in `P(N_3 >= 1) ~ c/log m`; anything about non-uniform (orbit-generated) words; the corresponding statement for necklaces rather than words (the same order holds since a primitive necklace has `K` rotations, each a word, and fairness is a rotation-invariant property of the cyclic word up to the cut position).
