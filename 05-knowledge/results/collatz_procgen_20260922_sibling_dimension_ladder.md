# Sibling dimension ladder: how big is the exceptional set, and what can be proved

**Status: PROVED: Lemma 1 and Theorems 1, 2, 3 (hand proofs below; Theorem 1(b)'s numerical
enclosures are computer-assisted with exact integer arithmetic and an analytic tail bound).
FINITE-EXACT: every count and census (two code paths agree where marked). CITED: Lagarias 2009,
Alon--Behajaina--Paran 2024, Inselmann 2024/25, Kontorovich--Lagarias 2009/10 (primary sources
read); HMYZ 2008 (via ABP's statement of it) and Behajaina--Paran 2023 (abstract only); THM-2228,
THM-3848 (repo canon). NUMERICAL (not proved): the `C(theta) 2^{hm} m^{-3/2}` law for the Collatz
counts, and the collapse of the `5x+1` choice game (Monte Carlo over random classes, each class
decided exactly). NOT CONFIRMED: the
proposed "drift barrier" (that choice cannot collapse the `q=5` exceptional set); the data point
the other way. OPEN: Collatz, `5x+1` divergence, Mahler's Z-numbers, Erdős's ternary problem,
E-SCC. No independent-agent audit yet.**

Session `collatz-procgen-20260922`, sibling-ladder lane. Inherits the definitions of lane one,
[choice_ladder](collatz_procgen_20260922_choice_ladder.md) (the descent game, `Bad_inf`, the graph
`E`, the partial-choice games `E_S`). Scripts:
`04-computation/experiments/collatz_procgen_20260922_ladder_*.{py,c}`; output
[collatz_procgen_20260922_ladder.out](collatz_procgen_20260922_ladder.out).

**Reproduce:** `python3 04-computation/experiments/collatz_procgen_20260922_ladder_run.py`
(about 9 minutes, peak memory about 0.7 GB, one process at a time; `--quick --out PATH` takes under
a minute with smaller samples). It writes `05-knowledge/results/collatz_procgen_20260922_ladder.out`.

## 1. The ladder

"Exceptional set" means the set of points of the completed space that never descend
multiplicatively (lane one's `Bad_inf`), or the problem's own analogue where that is not defined
(Mahler, Erdős). "All-orbits statement" means the problem itself.

| sibling | exceptional set | dimension or measure | all-orbits statement | source |
|---|---|---|---|---|
| Collatz `3x+1` (shortcut `T`, no choice) | `Bad(3)`: 2-adic `x` with more than `s log_3 2` odd steps in every prefix of length `s` | `dim_H = h(log_3 2) = 0.9499555` (PROVED, Thm 1a); Haar 0; count `N_m ~ C(theta_m) 2^{0.9499555 m} m^{-3/2}` (NUMERICAL to `m=10^5`) | OPEN. Collatz implies `Bad(3)` contains no positive integer; a positive integer in `Bad(3)` would have an unbounded orbit (PROVED, Thm 1c) | this note; measure 0 is classical (Terras 1976, Everett 1977) |
| `3x+r`, any odd `r` (`3x-1`, `3x+5`, ...) | the parity-word preimage of the same word set | the same `h(log_3 2)` (PROVED, Lemma 1) | OPEN; the extra cycles of `3x-1`, `3x+5` are invisible to the dimension | this note |
| `5x+1` | `Bad(5)` | Haar measure `mu_5 = 0.176025784562122695129` (error `< 10^-22`), so dimension 1 (PROVED, Thm 1b; Spitzer formula agrees to `5e-29`) | OPEN, believed false: a density-one set of integers is conjectured to diverge (Kontorovich--Lagarias). A positive integer in `Bad(5)` would diverge (PROVED, Thm 1c); none is known | this note; K--L |
| `7x+1` | `Bad(7)` | `mu_7 = 0.300751353095002575716` (PROVED) | OPEN, as for `5x+1` | this note |
| lane one's `E` (`3x+1` plus the arrow `n -> 3n+1` at evens) | `Bad(E)`, contains `-1`, `-13/9` (PROVED, lane one) | `777` classes mod `2^26`, `908` mod `2^36` (FINITE-EXACT); dimension 0 conjectured (HYP-9120) | relaxed Q1: OPEN (implied by Collatz) | lane one |
| `E_5` (`5x+1` with the same choice) | `Bad(E_5)`, contains `-1` (PROVED, §4) | level-`m` mass `0.0664` at `m=26` (FINITE-EXACT; a rigorous upper bound for the Haar measure, `< mu_5/2.6`), about `4e-5` at `m=64` (NUMERICAL); collapse conjectured | not studied | this note |
| `E_q`, odd `q>=41` | `Bad(E_q)` | Haar measure `>= (1-K_q)/2 > 0` (PROVED, Thm 3) | not studied | this note |
| Applegate--Lagarias multiplier semigroup | the single class `-1 mod 2^j` | one point | weak `3x+1` conjecture PROVED | CITED (A--L 2006, via lane one) |
| `F_2[x]` (HMYZ map) | `Bad_F = {1}` in `F_2[[x]]` | one point; `h(1)=0`; exactly one exceptional class mod `x^k` for every `k` (PROVED, Thm 2; FINITE-EXACT for `k<=20`) | PROVED: every nonzero polynomial reaches `1` (HMYZ), within `d^2+2d` HMYZ steps (re-proved in Thm 2), `O(d^1.5)` (ABP); census `deg<=24` FINITE-EXACT | HMYZ 2008; ABP 2024; Inselmann |
| Mahler `3/2` (Z-numbers) | strict safe carry words `S`; `Phi(S)` in `Z_2` | `log_2(3/2) = 0.5849625`, binary-ultrametric = 2-adic (PROVED, THM-3848 (S10); `Phi` is an isometry by THM-2228); Haar 0; count `a_m ~ 1.5510 (3/2)^m`, no polynomial factor (FINITE-EXACT recount) | OPEN (Mahler 1968): no Z-number iff `Phi(S)` contains no positive integer (THM-2228); a candidate exceeds `2^57` (Dubickas--Mossinghoff) | THM-3848, THM-2228 |
| Erdős, ternary digits of `2^n` | Lagarias's `E(Z_3)`: 3-adic `lambda` with infinitely many `lambda 2^n` omitting the digit 2 | `dim E^(1) = log_3 2 = 0.6309`; `(1/2)log_3 2 <= dim E^(2) <= 1/2`; `(1/6)log_3 2 <= dim E^(3) <= dim E^(2)` (CITED, Lagarias Thm 1.5). `dim E(Z_3) = 0` is his Conjecture B (OPEN) | OPEN (Erdős): equivalent to `1 not in E(Z_3)` | Lagarias, JLMS 79 (2009) |

## 2. The multiplier family `qx+1`: proofs

Fix odd `q>=3` and odd `r`, and let `T(x)=x/2` for even `x` and `T(x)=(qx+r)/2` for odd `x`, on
`Z_2`. Write `v_i(x)=T^i(x) mod 2` and `a_s(x)=v_0+...+v_(s-1)`. Put

```text
Bad(q) = {x in Z_2 : q^(a_s(x)) > 2^s for every s >= 1},     p = log_q 2 = 1/log_2 q.
```

This is the set of infinite coefficient stopping time. `q^a=2^s` is impossible for `s>=1`, so
strict and weak inequalities agree. `Bad_m(q)` imposes the condition for `s<=m`. It is a union of
`N_m(q)` classes mod `2^m`, and `f_m=N_m/2^m`.

**Lemma 1 (parity-vector isometry; classical: Terras, Lagarias 1985 Thm B, Bernstein--Lagarias 1996).**
(i) `T^s(x)=(q^(a_s) x + B_w)/2^s` with `B_w` an integer depending only on `w=(v_0..v_(s-1))`.
(ii) For every word `w` of length `s`, `{x : v(x) starts with w}` is exactly one class mod `2^s`.
(iii) Hence `x -> v(x)` is a bijective isometry from `(Z_2,|.|_2)` onto `{0,1}^N` with
`d(u,u')=2^(-min{i: u_i != u'_i})`, carrying Haar measure to the fair Bernoulli measure.

*Proof.* (i) by induction: an even step keeps `(q^a x+B)/2^(s+1)`, and an odd step gives
`(q^(a+1)x+qB+r2^s)/2^(s+1)`. (ii) by induction on `s`. If the fibre of `w` is `rho+2^s Z_2`,
write `x=rho+2^s t`. By (i), `T^s(x)=T^s(rho)+q^a t == T^s(rho)+t (mod 2)` since `q` is odd. So
`v_s(x)=e` exactly on one class of `t` mod 2, i.e. on one class of `x` mod `2^(s+1)`. (iii)
follows. ∎

Consequently `v(Bad(q)) = Sigma_p = {u : u_0+...+u_(s-1) > ps for every s>=1}` for every odd `r`,
and `N_m(q)` is a ballot number: the number of 0/1 words of length `m` every prefix of which has
more than `p` times its length in ones.

**Theorem 1.**
(a) If `1/2 <= p < 1` then `dim_H Sigma_p = h(p) = -p log_2 p-(1-p)log_2(1-p)`. Among odd `q` this
is only `q=3`: `dim_H Bad(3) = h(log_3 2) = 0.949955527...` in the 2-adic metric, for every odd
`r`, and `Bad(3)` is Haar-null.
(b) If `q>=5` then `Haar(Bad(q)) = mu_q := P(S_s>0 for all s>=1) > 0`, where
`S_s = a_s log_2 q - s` for fair coin flips. Moreover

```text
mu_q = exp( - sum_{n>=1} (1/n) 2^(-n) sum_{k <= floor(n log_q 2)} C(n,k) ),
0 <= f_m - mu_q <= rho^(m+1)/(1-rho),     rho = 2^(t-1)(1+q^(-t)) for any t>0.
```

(c) (`r=1`) A positive integer in `Bad(q)` has an unbounded orbit. In particular Collatz implies
`Bad(3)` contains no positive integer, and a positive integer in `Bad(5)` would give a divergent
`5x+1` orbit.
(d) The positive integers in `Bad(q)` have upper density at most `f_m` for every `m`, hence at most
`mu_q`. For `q=3` the density is 0.

*Proof of (a).* Upper bound. For each `s`, `Sigma_p` is covered by the cylinders of the words of
length `s` with more than `ps` ones. Since `p>=1/2`, `p^k(1-p)^(s-k)` is nondecreasing in `k`, so
`1 >= sum_{k>ps} C(s,k)p^k(1-p)^(s-k) >= #words * 2^(-s h(p))`. There are therefore at most
`2^(s h(p))` cylinders, each of diameter `2^-s`. So `H^t_(2^-s)(Sigma_p) <= 2^(s(h(p)-t)) -> 0`
for `t>h(p)`.

Lower bound. Take `p<p'<1`, `nu` the Bernoulli(`p'`) product measure, and
`Z_s = sum_{i<s}(u_i-p)`. By the strong law, `Z_s/s -> p'-p>0` `nu`-a.s., so `inf_s Z_s > -infinity`
a.s. Choose `K` with `nu(inf Z > -K)>0` and an integer `N>=K/(1-p)`. Every word `1^N w` with
`inf_s Z_s(w) > -K` lies in `Sigma_p`: `Z_s=s(1-p)>0` for `s<=N`, and `Z_s>N(1-p)-K>=0` after.
This set `E_N` has `nu(E_N) = p'^N nu(inf Z>-K) > 0`. Balls in `d` are cylinders, and
`-(1/s)log_2 nu([u_<s]) -> h(p')` for `nu`-a.e. `u` (strong law again). The mass distribution
principle (Falconer, *Fractal Geometry*, Prop. 4.9) gives `dim_H E_N >= h(p')`. Now let
`p' -> p`. Haar-nullity: `h(p)<1`. ∎

*Proof of (b).* By Lemma 1, `Haar(Bad(q))` is the fair-coin probability that `S_s>0` for all
`s>=1`, and `f_m` is the same probability for `s<=m`, exactly. For `q>=5`, `p<1/2`, so the drift
`(log_2 q)/2-1` is positive. The prefix argument of (a) with `p'=1/2` gives `mu_q>0`. The series is
the Sparre Andersen--Spitzer identity `P(tau=infinity)=exp(-sum_n P(S_n<=0)/n)` for the first weak
descending ladder epoch (Feller II, XII.7). Here `P(S_n<=0)=P(S_n<0)=P(a_n <= floor(n log_q 2))`.
For the tail, `f_m-mu_q = P(S_s>0 for s<=m, S_n<0 for some n>m) <= sum_{n>m} P(S_n<0)`, and
Chernoff gives `P(S_n<0) <= E 2^(-t S_n) = rho^n`. ∎

*Proof of (c).* A bounded orbit of a positive integer is eventually periodic. Take a cycle of
length `L` with `a` odd steps. A cycle element `y>0` satisfies `y(2^L-q^a)=B>0` by Lemma 1(i)
(`a>=1`, since `L` halvings cannot fix `y>0`). So `q^a<2^L`. The parity word is eventually periodic
with density of ones `a/L < log_q 2`, so `S_s -> -infinity` and the start is not in `Bad(q)`. For
Collatz: if `n>=2` reaches 1, then at its stopping time `T^s(n)=(3^a n+B)/2^s<n` with `B>=0`, so
`3^a<2^s`; and `1` descends at `s=2`. ∎

*Proof of (d).* `n in Bad(q)` forces `n mod 2^m` into `Bad_m(q)`, a set of density `f_m`; and
`f_m -> 0` for `q=3`. ∎

(c) uses `r>0`. It fails on the minus sheet: for `3x-1` the cycle `5 -> 7 -> 10 -> 5` gains `9/8`
per period, so the positive integer `5` lies in `Bad(3)` of `3x-1` although its orbit is bounded.
The dimension cannot see the sheet; the integer question can.

### 2.1 FINITE-EXACT and NUMERICAL support

* **Lemma 1 checked** as a bijection mod `2^k`, `k<=14`, for `q=3,5,7`. The exact ballot DP
  reproduces lane one's Collatz counts `N_m(3)` for all `m<=26` (`1,037,374` at `m=26`). It also
  reproduces the `mode 0` counts of lane one's C program for `q=3,5,7` at `m=16,20,24` (`q=5` also
  at `m=18,22,26`).
* **The dimension from counts (item 1, `q=3`).** The naive exponent `log_2 N_m/m` is `0.769` at
  `m=26`. It is `0.938` at `m=10^3` and `0.94974` at `m=10^5`. The exact DP (float DP beyond
  `m=2000`; it agrees with the big-integer DP to `7e-13` in `log_2`) fits

  ```text
  log_2 N_m = h m - (3/2) log_2 m + C(theta_m) - 60/m + ...,    theta_m = frac(m log_3 2).
  ```

  On thin slices of fixed `theta` (`m` in `[2000,10^5]`) the four-parameter fit returns
  `alpha = 0.94995552` (`h = 0.94995553`) and `beta = -1.4997, -1.4999, -1.4995`. `C(theta)` takes
  values in `[3.38, 3.49]` (`log_2`). It has downward jumps at `theta = frac(j log_3 2)`: `-0.055`
  at `j=1`, `-0.029` at `j=2,3`, `-0.022` at `j=4`, `-0.015` at `j=5,6`. Within 400 `theta`-bins
  the median spread falls from `0.036` (`m` in `[10^3, 5000]`) to `0.0017` (`[5*10^4, 10^5]`).
  **With `m<=26` alone the exponent cannot be pinned.** A fit with `beta=0` gives `0.871` and one
  with `beta=-3/2` gives `1.016`. The finite-size term is `-1.1` bits at `m=26`, decaying like
  `-60/m`. The exact decomposition at `m=26` is
  `0.7686 = h (0.9500) - 1.5 log_2(26)/26 (0.2712) + R_26/26 (0.0899)`.
* **`mu_q` (item 1, `q>=5`).** The exact DP to `m = 5736` (`q=5`) or fewer (larger `q`), plus the
  Chernoff tail (`< 10^-22`), gives:

  | `q` | `mu_q` (error `< 1e-22`) | `f_20` | `f_24` |
  |---|---|---|---|
  | 5 | `0.1760257845621226951290` | 0.2221 | 0.2130 |
  | 7 | `0.3007513530950025757155` | 0.3098 | 0.3079 |
  | 9 | `0.3861300244112210329781` | 0.3901 | 0.3890 |
  | 11 | `0.4002754722099983914490` | 0.4027 | 0.4016 |
  | 13 | `0.4172213116883590864120` | 0.4183 | 0.4179 |
  | 31 | `0.4637809877408409742713` | 0.4639 | 0.4639 |

  The Spitzer series, evaluated with exact binomial partial sums, agrees for `q=5..13` to
  `<2e-27`. The finite-level fractions approach `mu_5` slowly: `f_26 = 0.2094`, still `0.033` high.
  The tail constant `rho_5=0.99040` explains why counts to `m=26` alone cannot give the limit.

## 3. `F_2[x]`: the degenerate rung

`T(f)=f/x` if `f(0)=0`, and `T(f)=((x+1)f+1)/x` if `f(0)=1`, on `F_2[[x]]` ("odd" means `f(0)=1`).
This is the shortcut of the Hicks--Mullen--Yucas--Zavislak map `C(f)=(1+x)f+1` (odd), `f/x`
(even), since `C(C(f))=T(f)` for odd `f`. Define `Bad_F = {f : T^n(f) is odd for every n>=0}`.
For polynomials, "never an even step" is the same as "the degree never decreases".

**Theorem 2.**
(i) The parity-vector map of `T` is a bijective isometry `F_2[[x]] -> {0,1}^N` for the `x`-adic
metric. The proof is Lemma 1 with `q -> x+1`, a unit `== 1 mod x`, and `2 -> x`.
(ii) (Odd-run law.) If `f, Tf, ..., T^(k-1)f` are odd, then `T^k(f)+1 = (f+1)(x+1)^k/x^k`. The
initial run of odd iterates of `f` has length exactly `v_x(f+1)`.
(iii) `Bad_F = {1}`, and `T(1)=1`.
(iv) (HMYZ 2008.) Every nonzero polynomial reaches `1`, within `d(d+3)/2` steps of `T`, i.e.
`d^2+2d` steps of `C`, where `d=deg f`.

*Proof.* (ii) For odd `f`, `(x+1)f+1 = xf+(f+1)`, so `T(f)=f+(f+1)/x` and
`T(f)+1=(f+1)(1+1/x)=(f+1)(x+1)/x`. Iterating, `T^j(f)+1` has `x`-adic valuation `v_x(f+1)-j`,
because `v_x(x+1)=0`. `T^j(f)` is odd iff that valuation is `>=1`. (iii) By (ii), `f` is in
`Bad_F` iff `v_x(f+1)=infinity` iff `f=1`. Equivalently, `1` is the unique preimage of the word
`1^infinity` under (i). (iv) For odd `f` of degree `d>=0`, `(x+1)f+1` has degree `d+1`, so
`deg T(f)=d`. For even `f != 0`, `deg T(f)=d-1`, and `T(0)=0`. The orbit of a nonzero polynomial
stays among polynomials of degree `<= deg f`, so it is eventually periodic. On a cycle the degree
is constant, so every step is odd and the cycle lies in `Bad_F={1}`. For the count: at degree
`e>=1`, an odd `g != 1` runs `v_x(g+1) <= deg(g+1) = e` odd steps, then one even step lowers the
degree. That is at most `e+1` steps of `T`, or `2e+1` steps of `C`, and summing over
`e=1..d` gives `d(d+3)/2` and `d^2+2d`. ∎

`d^2+2d` is exactly the bound Alon--Behajaina--Paran attribute to HMYZ. We did not access the
*Monthly* article, so no claim is made that this proof differs from theirs. The power-series
statement (iii) is elementary. It is distinct from Behajaina--Paran's theorem that all but
countably many power series have non-eventually-periodic orbits (CITED). Those orbits descend
infinitely often, and only `1` never descends.

**FINITE-EXACT census** (`..._ladder_f2x.c`). All identities pass: (ii) for `2^16` odd `f`, and
`T^k(x^k u+1)=(x+1)^k u+1` for `k<=12`, `deg u<=12`. On all `2^25-1` nonzero polynomials of degree
`<=24`, the degree law, the odd-run law and `t_HMYZ = 2 t_T - deg` (Inselmann's identity) hold, and
every orbit reaches `1`. For `k<=20`, `f mod x^k -> first k parities` is a bijection, and the
exceptional count is exactly 1 (the class of `1`) at every level. Maximal stopping times:

| `d` | 8 | 12 | 16 | 20 | 24 |
|---|---|---|---|---|---|
| max `T`-steps | 21 | 35 | 51 | 71 | 88 |
| max HMYZ steps `sigma(d)` | 34 | 58 | 86 | 122 | 152 |
| `d^2+2d` | 80 | 168 | 288 | 440 | 624 |
| `sigma(d)/(d ln d)` | 2.04 | 1.95 | 1.94 | 2.04 | 1.99 |
| mean `T`-steps `/d` | 1.747 | 1.823 | 1.866 | 1.892 | 1.909 |

The mean ratio tends to 2 (Inselmann's theorem `rho_1(n)/n -> 2`), and `sigma(d)/(d ln d)` stays
near 2, in line with ABP's remark that experiments suggest `sigma(d) >= d log d`. Over all degrees
`<=20`, the maximum is 71 `T`-steps (122 HMYZ steps).

**Why this is the "`log_2 3 -> 1`" case.** For integers the size functional is `log_2|n|`. An odd
step adds `log_2 3 - 1` and an even step subtracts 1, so descent fails only while the density of
odd steps stays above `1/log_2 3 = log_3 2 = 0.63`. The words that manage this carry entropy
`h(0.63)=0.95`. In `F_2[x]` the size functional is the degree and `|x+1| = |x|`. The ratio
`log|x+1|/log|x|` that replaces `log_2 3` is exactly 1, so an odd step changes the degree by
`deg(x+1)-deg(x) = 0`. The threshold density becomes `p=1`: only the all-odd word survives,
`h(1)=0`, and the Cantor set shrinks to one point (Thm 2(iii)). Carries are what make the integer
window positive: `T^k(2^k u-1)=3^k u-1` grows by `(3/2)^k`, whereas `T^k(x^k u+1)=(x+1)^k u+1` has
the same degree. Moreover the one surviving point is the analogue of Applegate--Lagarias's
hostile `-1`, the all-odd word, and in characteristic 2 `-1=1` **is the terminal cycle**. No
transversality question is left, which is why this sibling is a theorem.

## 4. The choice ladder, and the drift barrier that is not there

Counts of exceptional classes, from lane one's program `exceptional_general.c` (FINITE-EXACT). The
no-choice rows equal the exact ballot numbers. The `q=3` rows equal lane one's independent programs
(`e_forward_dp_full` for full `E` at even `m` from 14 to 26, `e_forward_dp` for the two `E_S` rows at
`m=18,22`, and lane one's Collatz counts at every `m<=26`). The exhaustive per-class search
`..._ladder_eq_dfs.c` reproduces all six games checked at `m=16,20`. A scratch run, not in the
`.out`, also gives `1,427,116` for `E_5` at `m=24`.

| game | `m=20` | `m=24` | fraction at `m=24` |
|---|---|---|---|
| `q=3`, no choice (Collatz) | 27,328 | 286,581 | 0.01708 |
| `q=3`, `E_S`, `S={6 mod 8}` | 2,165 | 2,109 | 0.00013 |
| `q=3`, `E_S`, `S={2 mod 4}` | 1,254 | 611 | 0.00004 |
| `q=3`, full choice `E` | 664 | 369 | 0.00002 |
| `q=5`, no choice (`5x+1`) | 232,912 | 3,573,290 | 0.21298 |
| `q=5`, `E_S`, `S={6 mod 8}` | 219,070 | 3,225,154 | 0.19223 |
| `q=5`, `E_S`, `S={2 mod 4}` | 137,762 | 1,551,497 | 0.09248 |
| `q=5`, full choice `E_5` | 133,600 | 1,427,116 | 0.08506 |
| `q=7`, no choice | 324,884 | 5,166,172 | 0.30793 |
| `q=7`, full choice `E_7` | 312,884 | 4,916,894 | 0.29307 |

`E_q` adds the arrow `n -> qn+1` at even `n`, i.e. the excursion `x -> q^2x+q+1` costing
`2 log_2 q` bits. For `q=5` the rising-run entries are the evens `== 2 mod 4` (`(5x+1)/2` is odd
iff `x == 1 mod 4`). Accordingly `S={2 mod 4}` recovers almost all of the full-choice gain, while
`S={6 mod 8}`, lane one's winner for `q=3`, recovers little.

**The proposed drift barrier is not confirmed.** The hypothesis was that the positive drift of
`5x+1` stops choice from collapsing its exceptional set. The fraction does not stay put:

| `m` | 18 | 22 | 26 | 32 | 40 | 48 | 56 | 64 |
|---|---|---|---|---|---|---|---|---|
| no choice `f_m(5)` (exact) | 0.2283 | 0.2159 | 0.2094 | 0.2039 | 0.1963 | 0.1925 | 0.1886 | 0.1867 |
| `E_5` (exact to 26, then MC) | 0.1516 | 0.0995 | 0.0664 | 0.0376 | 0.0088 | 0.0019 | 0.00031 | 0.00004 |

The rows from `m=32` are Monte Carlo over uniform random classes mod `2^m`, 100,000 per level,
fixed seeds. **Each sampled class is decided exactly** by the branch-and-bound search. Controls:
at `m=26` the MC gives `0.0658 +- 0.0008` against the exact 0.06642; the no-choice MC at `m=64`
gives `0.1868 +- 0.0012` against the exact ballot value 0.18668. At `m=64`, 4 of 100,000 classes
survive (approximate 95% upper bound `1e-4`). The growth exponent of the surviving count between
successive levels falls: 0.87, 0.74, 0.73, 0.67, 0.63. Two statements follow.

* Rigorous: `Haar(Bad(E_5)) <= f_26(E_5) = 0.06642 < mu_5/2.6`. Choice provably removes most of
  the no-choice exceptional set, even for `5x+1`.
* Numerical: the level-64 mass is about `4e-5` and still falling. The data favour **collapse**,
  i.e. `Haar(Bad(E_5)) = 0`. That is CONJECTURED, not proved. `Bad(E_5)` is nonempty (it contains
  `-1`, below).

For larger multipliers the decline slows and then stops within the resolution (MC, 20,000 classes
per point, standard error about 0.003; no-choice `mu_q` for reference):

| `q` | `mu_q` | `E_q`, `m=32` | `m=48` | `m=64` | change 32 to 64 |
|---|---|---|---|---|---|
| 7 | 0.3008 | 0.2781 | 0.2451 | 0.2183 | `-0.060` (clear decline) |
| 9 | 0.3861 | 0.3676 | 0.3472 | 0.3374 | `-0.030` (clear decline) |
| 11 | 0.4003 | 0.3927 | 0.3916 | 0.3864 | `-0.006` (1.3 s.e.) |
| 13 | 0.4172 | 0.4127 | 0.4104 | 0.4062 | `-0.007` (1.3 s.e.) |

A slow decline at `q=11,13` cannot be excluded at this resolution.

**Theorem 3 (choice cannot collapse large multipliers).** Let `E_q` be the full-choice game, with
a certificate being a legal path after one of whose halvings `q^a<2^b`. Put
`r_q(t)=2^(t-1)/(1-q^(-t))` and `K_q(t)=2^t q^(-t)/((1-q^(-2t))(1-r_q(t)))`. If `r_q(t)<1`, then
`Haar(Bad(E_q)) >= (1-K_q(t))/2`. At `t=0.713`, `K_41 = 0.98912` and `K_q` decreases in `q`. So
`Haar(Bad(E_q))>0` for every odd `q>=41` (`>=0.0054` at `q=41`, `>=0.29` at `q=101`). The same
bound holds for every partial-choice game `E_q^S` and for the plain map. Also `-1 in Bad(E_q)` for
every odd `q>=3`.

*Proof.* An even `x` descends at once, which accounts for mass `1/2`. For odd `x` the first move is
forced (`a=1`), and every later move sits at an even point except right after a halving. Describe a
certificate by its skeleton `(j_1,e_1,j_2,...,e_(n-1),j_n)`: `j_i>=0` excursions before the `i`-th
halving, and `e_i` the parity of the point just after it. The skeleton is realized exactly on one
class mod `2^n`, because the path is affine with odd multiplier and each `e_i` reads one new bit
(as in Lemma 1). So it has probability `2^-n`, and `a = 1+sum e_i+2 sum j_i`. By the union bound
and `1[q^a<2^n] <= 2^(tn) q^(-ta)`,

```text
P(x odd, x not in Bad) <= q^(-t) 2^(t-1)/(1-q^(-2t)) * sum_{n>=1} r_q(t)^(n-1) = K_q(t)/2,
```

using `sum_j q^(-2tj) = 1/(1-q^(-2t))` and `(1+q^(-t))/(1-q^(-2t)) = 1/(1-q^(-t))`. For `-1`:
every `E_q`-path from `-1` stays at negative integers `(-q^a+B)/2^b <= -1` with `B>=0`, so
`q^a >= 2^b`, and equality is impossible (Applegate--Lagarias's argument, as in lane one). ∎

**The certificate index.** `r_q=min_t r_q(t)` is the growth rate, per bit, of the expected number of
descending certificates. It replaces the no-choice Chernoff index `rho_q=min_t 2^(t-1)(1+q^-t)`,
which is `<1` iff `q>=5`. Values: `r_3=1.485`, `r_5=1.200`, `r_7=1.092`, `r_9=1.032`,
`r_11=0.993`, `r_13=0.965`. For `q<=9` certificates multiply faster than the bits they consume, so
no first-moment argument can protect a positive-measure set. Heuristically this is a supercritical
branching random walk whose minimum runs to `-infinity`. The observed behaviour matches: `q=5`
collapses, `q=7,9` decline slowly, `q=11,13` are flat within resolution. **The threshold for choice
games is heuristically `r_q=1`,
between `q=9` and `q=11`, not the deterministic drift threshold `q=4`.** Positivity for
`11<=q<=39` and collapse for `q=5,7,9` are OPEN.

## 5. Mahler `3/2` (extracted from THM-3848, re-counted)

THM-3848 §6 (S10) and THM-2228 give the following.

* **The set.** `K={c in {0,1}^N : Y_i(c)=sum_j c_(i+j)(2/3)^(j+1) <= 1 for all i}` is the closed
  formal safe-tail shift, and `S` (strict `<1`) removes the countable orbit of the greedy word
  `d`.
* **The dimension.** `h_top(K)=log(3/2)` and `dim_H K = dim_H S = log_2(3/2) = 0.5849625`, by the
  prefix-cylinder formula in the binary ultrametric. `K` is Haar-null (mass `a_m/2^m ~ C(3/4)^m`)
  and nonsofic.
* **The 2-adic reading.** THM-2228's `Phi` sends length-`m` carry words bijectively onto `Z/2^m`
  and conjugates the shift to `a -> (3a+(a mod 2))/2 = ceil(3a/2)`. It is the parity-vector isometry
  of that map (Lemma 1's proof applies verbatim, since both branches multiply by the odd number 3), so
  `Phi(S)` is a Haar-null subset of `Z_2` of 2-adic dimension `log_2(3/2)`.
* **Its meaning.** A positive Z-number exists iff `Phi(S)` contains a positive integer (THM-2228
  (12)-(13)). Mahler's problem thus has exactly the shape of the Collatz question: do the positive
  integers avoid a Haar-null 2-adic Cantor set? Here the dimension is 0.585 instead of 0.950.

Re-count (FINITE-EXACT, `..._ladder_mahler.py`). Three methods give identical counts for `m<=30`:
a DFS over the strict suffix inequalities (S2), the renewal law (S8), and Parry's lexicographic
condition (every suffix `<=_lex` the prefix of `d`). The counts are `1,2,3,5,8,12,18,27,40,60,90,134,...`.
So `K` is the closed `beta`-shift for `beta=3/2`, whose entropy `log beta` is classical
(Rényi--Parry); the identification is standard and is checked here only finitely. Also
`a_m/(3/2)^m -> 1.5510451884` with no polynomial factor. Contrast Collatz's `m^(-3/2)`: `K` is
shift-invariant (a constraint on every suffix), while `Bad(3)` is a ballot set (a constraint
anchored at the start). The number `log_2(3/2)=log_2 3-1` is also the per-odd-step growth of the
Collatz map, which opens the positive-entropy window of §3.

## 6. Erdős's ternary problem (Lagarias 2009), the transversality rung

CITED from J. C. Lagarias, "Ternary expansions of powers of 2", *J. London Math. Soc.* (2) 79
(2009) 562--588, arXiv math/0512006v4 (read):

* **Truncated real system.** `N_lambda(X) <= 25 X^0.9725` for every `lambda>0` (Thm 1.1). The
  truncated real exceptional set has `dim_H = log_3 2 = 0.63092`, with positive
  `log_3 2`-dimensional measure (Thm 1.3).
* **3-adic system.** For nonzero `lambda in Z_3`, `#{n<=X : lambda 2^n omits 2} <= 2X^(log_3 2)`
  (Thm 1.4). For the nested sets `E^(k)(Z_3)` (at least `k` of the `lambda 2^n` omit 2):
  `dim E^(1) = log_3 2`, `(1/2)log_3 2 <= dim E^(2) <= 1/2`, and
  `(1/6)log_3 2 <= dim E^(3) <= dim E^(2)` (Thm 1.5).
* **Intersections of translates.** `dim(Sigma_(3,2bar) cap M^(-1) Sigma_(3,2bar)) <= 1/2` for `M`
  not a power of 3, with `= log_3((1+sqrt5)/2) = 0.438` for `M=7` (Thm 1.6).
* **Conjectures.** Conjectures A and B: the real and 3-adic exceptional sets `E(R+)`, `E(Z_3)` have
  dimension 0. Erdős's conjecture is equivalent to `1 not in E(Z_3)`. Conjecture E: every digit
  pattern eventually appears in `(p^n)_q`.

Lagarias frames this with Furstenberg's transversality (§5): `dim(A cap B) <= max(dim A+dim B-dim X,0)`
for closed sets invariant under transverse semigroups. The ladder's rungs share that shape: a
dynamically defined Cantor set in a `p`-adic space, of known or bounded dimension, and an
arithmetic question about whether a countable orbit-set (the positive integers, or the single orbit
of `1`) meets it. A November 2025 preprint (Roettger--Ren, arXiv 2511.03861) states that Erdős's
conjecture is still open; a quick search on 2026-09-22 found nothing later.

## 7. What this says

* **Dimension is a property of the multiplier, not of the arithmetic.** Lemma 1 makes every `qx+r`
  (`r` odd) isometric to the same word set, so `3x+1`, `3x-1` and `3x+5` have the same
  0.95-dimensional exceptional set. At the 2-adic level even `3x+1` and `5x+1` are conjugate to the
  same shift (Kontorovich--Lagarias §10). What separates them is the archimedean weight `q^a/2^s`
  that defines the exceptional set: dimension 0.95 with measure 0, against measure 0.176. Everything else,
  including which integers are cyclic or divergent, is a **transversality** question: does the
  countable set `Z_{>=1}` meet the Cantor set? Collatz (0.95), Mahler (0.585) and Erdős
  (`<=1/2`, conjecturally 0) are three instances.
* **The solved sibling is the one with no transversality left.** In `F_2[x]` the exceptional set is
  one point, and that point is the terminal cycle (`-1=1`). A smaller positive dimension does not
  help: Mahler's 0.585 is as open as Collatz's 0.95.
* **Choice is too strong to be diagnostic.** The `E` relaxation collapses the exceptional set of
  `5x+1` (numerically), where almost all orbits are expected to diverge, just as it does for
  `3x+1`. Its own threshold is the certificate index `r_q=1` (`9<q*<=11`), not the drift threshold
  `q=4`. So exceptional-set collapse under choice cannot by itself be the mechanism that
  distinguishes `3x+1`. Lane one's reading, that the obstruction is "a no-choice phenomenon",
  stands but is not specific to `3x+1`. The E-SCC program (HYP-9120) should budget for a
  relaxation that is blind to drift.
* **Counting lesson.** Level-`m` counts at `m<=26` understate the Collatz dimension by 0.18
  (`m^(-3/2)` ballot factor plus a `-60/m` transient), and overstate `mu_5` by 0.03. The exact
  ballot DP settles both. Fits of `log(count)/m` from brute force at `m<=26` should not be quoted
  as dimensions.

**Discrepancy to integrate.** The session synthesis
([synthesis](collatz_procgen_20260922_synthesis.md) §1) lists "`5n+1`, with or without choice |
positive measure (DRIFT)". With choice, this is contradicted: it is rigorously at most 0.0664,
against the no-choice value 0.1760, and numerically about `4e-5` at level 64. It should read
"without choice: positive measure `mu_5=0.17603` (PROVED); with full choice: collapses
numerically".

## 8. Files

* `04-computation/experiments/collatz_procgen_20260922_ladder_ballot.py`: exact ballot DP, `q=3`
  dimension fit to `m=10^5`, `mu_q` enclosures, Spitzer check.
* `04-computation/experiments/collatz_procgen_20260922_ladder_f2x.c`: `F_2[x]` identities,
  census, and parity bijection.
* `04-computation/experiments/collatz_procgen_20260922_ladder_mahler.py`: the THM-3848 language
  counted three ways.
* `04-computation/experiments/collatz_procgen_20260922_ladder_choice.py`: the choice table via lane
  one's `exceptional_general.c` (unchanged), with cross-checks.
* `04-computation/experiments/collatz_procgen_20260922_ladder_eq_dfs.c`: an independent exact
  per-class search, exhaustive or Monte Carlo up to `m=64`.
* `04-computation/experiments/collatz_procgen_20260922_ladder_choice_mc.py`: validation, the MC
  tables, and the certificate indices.
* `04-computation/experiments/collatz_procgen_20260922_ladder_run.py`: the runner.

## Sources

* K. Hicks, G. L. Mullen, J. L. Yucas, R. Zavislak, "A polynomial analogue of the 3N+1 problem",
  *Amer. Math. Monthly* 115 (2008) 615--622 ([JSTOR](https://www.jstor.org/stable/27642557));
  bound `d^2+2d` as reported by ABP.
* G. Alon, A. Behajaina, E. Paran, "On the stopping time of the Collatz map in F2[x]",
  [arXiv:2401.03210](https://arxiv.org/abs/2401.03210) (read): `O(deg^1.5)`, Question 4.1, and the
  `d log d` remark.
* M. Inselmann, "On the average stopping time of the Collatz map in F2[x]",
  [arXiv:2401.12781](https://arxiv.org/abs/2401.12781) (read): `rho_1(n)/n -> 2`,
  `tau_0 = 2 tau_1 - deg`.
* A. Behajaina, E. Paran, "The Collatz problem in F_p[x] and F_p[[x]]", *Finite Fields Appl.* 91
  (2023) 102265 (abstract read); also [arXiv:2312.00390](https://arxiv.org/abs/2312.00390).
* A. V. Kontorovich, J. C. Lagarias, "Stochastic models for the 3x+1 and 5x+1 problems",
  [arXiv:0910.1944](https://arxiv.org/abs/0910.1944) (read §1, §10 summary).
* J. C. Lagarias, "Ternary expansions of powers of 2", *J. London Math. Soc.* 79 (2009) 562--588,
  [arXiv:math/0512006](https://arxiv.org/abs/math/0512006) (v4 read).
* C. Roettger, X. Ren, "Ternary digits of powers of two",
  [arXiv:2511.03861](https://arxiv.org/abs/2511.03861) (abstract read: Erdős's conjecture still open).
* Classical, not re-read this session: R. Terras (1976), C. J. Everett (1977), J. C. Lagarias
  (1985) Thm B, D. J. Bernstein and J. C. Lagarias (1996), W. Feller vol. II §XII.7
  (Sparre Andersen--Spitzer), K. Falconer, *Fractal Geometry* Prop. 4.9, H. G. Eggleston (1949),
  A. Rényi (1957) and W. Parry (1960) on `beta`-shifts. Mahler 1968 and Dubickas--Mossinghoff
  2009 via the repo ledger [CORE-PAPERS-MAHLER-THREE-HALVES](../reference/CORE-PAPERS-MAHLER-THREE-HALVES.md).
