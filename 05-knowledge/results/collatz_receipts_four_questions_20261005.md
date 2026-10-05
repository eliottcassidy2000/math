# The four receipt questions (2026-10-05): surplus receipts contain the actual orbit, so "deficits below the source" is the stopping time and the product identity plays no role; the Applegate--Lagarias coverage is a prefix-code descent with multipliers from the `3x-1` trunk on `20.3%` of odd residues, and those insertion points are exactly the composite defects of the universal receipts; product joins are hub coincidences, synchronous fusion is impossible in one step (F1), empty in two steps to `3001`, sporadic deeper; the first deficit of a multiple of three is below it iff `m = 3 mod 4`; a substitution chain carries one obligation; exact entry is terminal-family membership

**Session:** opus, `pascal-lyapunov-20261005` (second task; worktree `codex/session-pascal-lyapunov-20261005`), 2026-10-05.
Owner's directive: "pursue Collatz coverage creatively and prioritize the four helper questions" of the
[universal weak receipts note](collatz_universal_weak_receipts_20261005.md) (section 8 there; the owner's
mail restates them): (Q1) parameterized edge replacements eliminating whole families of composite defects;
(Q2) deficits of source-anchored multiples of 3 forced into certified families or below an induction bound;
(Q3) compiled controller substitutions preserving a small family of defect coordinates; (Q4) several rooted
terminal families turning local coverage into exact entry.
**Inherits (cited):** that note (W1, D1, D2, A1-A3, G1, G2, F1, F2, S1); Applegate--Lagarias,
[The 3x+1 semigroup](https://arxiv.org/abs/math/0411140) (Theorem 1.1, Lemmas 2.1-2.3, 3.1, Tables 1-3);
[THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md); the
[partition-cover obstruction](collatz_partition_cover_20261004.md); the
[dependency kernel](collatz_recursive_dependency_kernel_20261004.md); the
[Pascal-tower note](collatz_pascal_tower_lyapunov_20261005.md) section 5 (the no-descent rate).

**Status: PROVED (elementary; Theorems 1-7 below, each a few lines); CITED (the semigroup theorem and its
proof structure, read from the paper); FINITE-EXACT (two scripts, section 8: `11,325` pairs for joins, `250,000`
for synchronous fusion, `2,250,000` for two-step fusion, `500,000` multiples of three, the `26` semigroup
classes replayed); CONJECTURED (no synchronous two-step fusion, section 4). Author-audited only; audit OWED.
Nothing here is a Collatz step; every question is shown to be a known equivalent of the problem or a
bookkeeping fact, and the creative content is the identification of where the universal receipts' defects
actually sit.**

---

## 0. What is new, in one screen

1. **A surplus receipt contains the actual orbit (Theorem 1).** If `d c(n) > 0`, the actual trajectory of `n`
   follows edges of `c` until it meets a deficit vertex. Hence "a receipt with surplus at `n` whose deficits are
   all below `n`" exists iff `n` has finite stopping time; the product identity (W) never enters. Q2 is
   therefore the stopping-time conjecture on multiples of three (plus closure), and no receipt-side
   manipulation changes that.
2. **Defect elimination is rootedness (Theorem 2).** A family of composites whose defects can be eliminated by
   nonnegative edge replacements is a rooted family (G2). What replacements can do is transport a root
   certificate along a join; the depth-one joins `U(ab) = U(a)` are exactly F2 (complete classification,
   Theorem 3), the depth-two joins form explicit congruence families (Proposition 4, examples for every
   `a <= 51` but `19, 35`), and all of them are predecessor-tree membership.
3. **Multiplicative joins are hub coincidences, and synchronous fusion is empty where it matters.** The orbit of
   `ab` contains a product of orbit points of `a` and `b` for `60.8%` of the `11,325` pairs `3 <= a <= b <= 301`,
   but the join points are the hubs `175, 65, 35, 377, 319, 325, ...`; synchronous fusion `U^i(ab) = U^i(a)
   U^i(b)` is impossible for `i = 1` (F1), has **no solution for `i = 2`** with `a <= b <= 3001` (Conjecture 5:
   none at all; the size ratio is confined to `{8, 16}`), and is sporadic deeper (`103` pairs below `1001`, e.g.
   `U^3(37^2) = U^3(37)^2 = 17^2`).
4. **The semigroup theorem, dissected (section 5).** Its universal coverage is a prefix code of residue classes
   mod `2^12`: on `51/64` of the odd residues the actual orbit descends below `(76/79) x` within at most five
   odd steps; on `415/2048 = 20.26%` a multiplier `m in {5, 7, 11, 13, 23, 29, 43}` (or `25, 35`) is inserted at an
   odd point `y` of the orbit and the descent continues from `m y`; the last cell `-1 mod 4096` is handled by the
   multiplier `m_j = (2^j+1)/3`, for which `U(m_j x) = oddpart(x + (x+1)/2^j)` whenever `x = -1 mod 2^j`
   (Theorem 6) -- `m_j` is the trunk of the `3x-1` map (`3 m_j - 1 = 2^j`). **The composite defects of the
   universal receipts are precisely the labels `[m y]` at these insertion points**, and the `20.26%` is the
   no-descent set at depth twelve (the renewal value is `P(Bin(19,1/2) >= 12) = 17.96%`): the universal weak
   coverage is plain descent where plain descent exists and a multiplier where it does not.
5. **Multiples of three (Theorem 7).** The receipt of `n = 3m` must start with `3m -> U(3m)`, and that first
   deficit is below `n` iff `m = 3 mod 4` (`500,000/500,000`); beyond that the deficits-below-`n` condition is
   the Terras stopping time (half the multiples of three need more than one step, `6.5%` more than ten, the
   deepest below `3 10^6` needs `141`). An induction that stays inside the multiples of three exists (the
   3-parent rank, section 6) and needs a descent about `2^6` deeper; its hypothesis is again Collatz.
6. **Substitutions carry one obligation (Theorem 8).** A compiled common-future substitution `n -> d` is the
   signed packet with boundary `[n] - [d]`, defect `[1] - [d]`: one coordinate, preserved by composition;
   fusion appears only under multiplication of receipts. Q3 is answered, and the whole difficulty is the rank
   of the one obligation, which is the repo's credit-potential problem.
7. **Exact entry is membership (Theorem 9).** Replacing the terminal ray by the family `T_j` of integers
   reaching `1` within `j` odd steps, the S1 shadow is exact (`r = n`) iff `U^m(n) in T_j`. The least-exponent
   shadow is sometimes smaller than `n` (`n = 27, 255, 999` give `r = 3`) and sometimes `2^12 n`
   (`n = 97`); it is never a route from `n`.

---

## 1. Objects (from the receipts note)

`U(u) = (3u+1)/2^(v_2(3u+1))` on positive odd `u`. An edge multiset `c` (nonnegative integer multiplicities
on edges `u -> U(u)`) has value `V(c) = prod (u/U(u))^c(u)`, boundary `d c = sum c(u)([u] - [U(u)])` in the
free abelian group on labels, and, for a target `n`, defect `D_n(c) = d c - [n] + [1]`. A weak receipt for `n`
has `V(c) = n` (W1: every odd `n` has one). Surplus/deficit vertices are those with positive/negative
boundary. G1: if `d c(n) > 0` and every deficit vertex is rooted, `n` is rooted. G2: `n` reaches `1` iff some
nonnegative `c` has `d c = [n] - [1]`. `K_u`, `R(a,b) = [ab] + [1] - [a] - [b]` are the composite coordinates
and the fusion relation (D1, D2).

## 2. Theorem 1: a surplus receipt contains the actual orbit; the product identity is idle

**Theorem 1.** Let `c` be any finite nonnegative edge multiset of a deterministic map on a set `X`, and let
`n` have `d c(n) > 0`. Then the actual trajectory `n, U(n), U^2(n), ...` uses only edges with positive stock
in `c` until it first visits a deficit vertex, and it does visit one. Consequently, for the Collatz map:
`n` admits a receipt (with or without (W)) whose surplus is at `n` and whose deficit vertices other than `1`
are all `< n` **iff** `n` has finite stopping time (some `U^k(n) < n`, or `n = 1`).

*Proof.* This is G1's argument with `R` = the set of deficit vertices: if the supported trajectory from `n`
never meets a deficit vertex, let `S` be the set of vertices it visits (through its first repetition or
through the vertex where it stops for lack of stock); no positive-stock edge leaves `S`, so
`sum_(v in S) d c(v) <= 0`, while every vertex of `S` has boundary `>= 0` and `n` has `> 0`: contradiction. For
the equivalence: a trajectory prefix cut at the first `U^k(n) < n` is a receipt with one deficit (at that
point, or at `1`); conversely the trajectory reaches a deficit vertex `v`, which is `1` or `< n`. The
identity (W) was not used. □

**Reading.** The ladder "nonempty weak fibre -> marked prefix -> surplus -> grounded deficits" has its whole
content in the last rung, and that rung is the classical one: with the integer as rank, Q2's "every deficit
endpoint below a valid induction bound" is Terras's finite stopping time, for which the receipt calculus
offers no new leverage; with "an already certified family" it is the join problem of section 3.

## 3. Theorem 2, Theorem 3, Proposition 4: what edge replacements can and cannot do (Q1)

**Theorem 2.** A nonnegative receipt for `ab` with zero composite defect is a root path for `ab` (G2). Hence a
"family of composite defects eliminated by nonnegative edge replacements" is a rooted family, and the
replacement is a root-certificate transport along a join `U^i(ab) = U^j(a)` (or with `b`, or with any rooted
number). □

**Theorem 3 (depth-one joins; F2 is complete).** For odd `a, b > 1`: `U(ab) = U(a)` iff `b = 2^k + (2^k - 1)/(3a)`
with `ord_(3a)(2) | k`; `k` is then even (so this is F2 with `4^k`), and `b` is automatically odd.
*Proof.* Equality of the odd parts of `3ab+1 > 3a+1` means `3ab + 1 = 2^k (3a+1)` with `k >= 1`, i.e.
`b = (2^k(3a+1) - 1)/(3a) = 2^k + (2^k - 1)/(3a)`, integral iff `3a | 2^k - 1`; since `3 | 3a` and
`ord_3(2) = 2`, `k` is even; `(2^k - 1)/(3a)` is odd, so `b` is odd. □

**Proposition 4 (depth-two joins).** `U^2(ab) = U(a)` with first valuation `alpha` at `ab` iff
`b = [2^A (3a+1) - 3 - 2^alpha]/(9a)` is an odd integer `> 1` for some `A >= alpha + 1` and the valuation of
`3ab + 1` is exactly `alpha` (then `A = alpha + alpha_2 - v_2(3a+1)`). For each `a` the admissible `A` form
residue classes modulo the order of `2` modulo `9a/gcd`, so each solvable `(a, alpha)` gives an infinite family.
Script D (section 8) lists, for every odd `3 <= a <= 51` except `19, 35` (none with `alpha <= 13`, `A <= alpha
+ 60`), the smallest member: e.g. `a = 3`: `b = 23` (`69 -> 13 -> 5 = U(3)`), `a = 7`: `b = 11` (`77 -> 29 -> 11 =
U(7)`), `a = 9`: `b = 11` (`99 -> 149 -> 7 = U(9)`), `a = 31`: `b = 43`, `a = 41`: `b = 43`. These are depth-two
predecessors of `U(a)` divisible by `a`; the general depth-`i` family is the predecessor tree of `U^j(a)`
intersected with `a Z`, i.e. the inherited inverse rays in multiplicative coordinates. □

**What this settles for Q1.** Yes, infinitely many parametrized families of composite defects are
eliminable, one for every join; and no, this never reaches a composite whose orbit does not meet a rooted
orbit, because elimination is rootedness (Theorem 2). The multiplicative coordinate adds the possibility of a
join at a product point, tested next.

## 4. Product joins and synchronous fusion (Q1, the multiplicative option)

A product join `U^i(ab) = U^j(a) U^l(b)` with all factors `> 1` would transport the obligation `R(a,b)` to
`R(U^j(a), U^l(b))` at a composite that is smaller when both factors have descended. Census
(`collatz_receipt_joins_census_20261005.py`, Part A; `3 <= a <= b <= 301`, `11,325` pairs):

| statistic | pairs | share |
|---|---|---|
| ordinary join (an orbit point `> 1` of `a` or `b` lies on the orbit of `ab` at depth `>= 1`) | `10,623` | `93.8%` |
| ... above the trunk ray `(4^e - 1)/3` | `7,525` | `66.4%` |
| product join, all three factors `> 1` | `6,881` | `60.8%` |
| ... with the join point `z <= max(a,b)` (hub) / `max(a,b) < z < ab` (mid) / `z >= ab` (upward) | `4,471 / 4,773 / 258` | `39.5 / 42.1 / 2.3%` |
| synchronous (`i = j = l`) | `22` | `0.19%` |

The join points are overwhelmingly the hubs through which most orbits pass: `z = 175` (`2,938` hits), `65`
(`1,397`), `35` (`1,256`), `377, 319, 325, 445, 55, 155, 91, 121, 425`; `x = 5 = U(3)` and `y in {7, 13, 17, 35,
65}` make up almost every listed example. A product join is a coincidence of two popular small numbers,
not a multiplicative structure of the map: it transports `R(a,b)` to `R(5, 7)`, which is "35 is rooted".

**Synchronous fusion.** F1 (cited) forbids `U(ab) = U(a) U(b)` because the size ratio lies strictly between
`3` and `4`. Two steps: `U^2(ab) = U^2(a) U^2(b)` forces the ratio `rho = (9a + 3 + 2^alpha)(9b + 3 + 2^beta) /
(9ab + 3 + 2^gamma)` (`alpha, beta, gamma` the first valuations of `a, b, ab`) to be a power of two, and the
bounds `2^alpha <= 3a + 1`, `2^gamma <= 3ab + 1` give `6.75 < rho < 16 (1 + ((a+b)/3 - 4/9)/(ab + 5/9))`, so
`rho in {8, 16}`: the one-step size argument no longer excludes everything. Search (addendum (b)): **no
solution with `3 <= a <= b <= 3001` and all three factors `> 1`** (`2.25 10^6` pairs).

**Conjecture 5 (no synchronous two-step fusion; no HYP number assigned).** `U^2(ab) != U^2(a) U^2(b)` for all
odd `a, b > 1` with `U^2(a), U^2(b) > 1`. Evidence: the search above; the two candidate ratios `8, 16` each
require a 2-adic coincidence between `3 + 2^gamma` and the product of the two first-step constants that the
census never produced. Deeper synchronous fusions do exist and are sporadic: `103` pairs with `a <= b <= 1001`,
at depths `i in {3, 5, 8, 10, 34, 39, ...}`, e.g. `U^3(37^2) = U^3(37)^2 = 289` (`1369 -> 1027 -> 1541 -> 289`,
`37 -> 7 -> 11 -> 17`), `U^5(31 * 115) = U^5(31) U^5(115) = 121 * 7`. They carry no induction (the depths are
large and the pairs isolated).

## 5. The semigroup theorem dissected: where the universal receipts keep their defects (Q1-Q3)

Applegate--Lagarias prove Theorem 1.1 (the semigroup generated by `2` and `(2k+1)/(3k+2)` is `{a/b : 3 does
not divide b}`) by a see-saw induction on `k >= 12` with three hypotheses: (1) every `x > 1` with `x != -1 mod
2^k` has `s` in the inverse ("wild") semigroup `W` with `s x` an integer `<= (1235/1264) x`; (2) the weak
conjecture to `2^k - 2`; (3) the wild-numbers conjecture to `(2^k - 1)/189`. Hypothesis (1) is their Lemma 2.1:
Tables 1-2 give a prefix code of residue classes mod `2^k`, `k <= 12`, covering every class except `-1 mod
4096`, each with a path of `T`-steps and multiplications by elements of `H = {5, 7, 11, 13, 23, 29, 43}` (or
`25, 35`), ending below `(76/79) x`. In receipt terms a multiplier class is the nonnegative packet

```
{ x -> y } + { m y -> z } + w_m,       V = (x/y)(m y/z)(1/m) = x/z,
d = [x] - [y] + [m y] - [z] + d w_m,   defect relative to [x] - [z]:   [m y] - [y] + d w_m,
```

where `w_m` is the wild certificate of `1/m` (Table 3 of the paper; for `m = 5` the packet `{7:2, 11, 17, 55,
65, 83}` of the receipts note). **The composite coordinate of the universal receipt is `[m y]`, the label at
which the paper multiplies.** Replayed with the accelerated map (addendum (c)): all `26` listed classes
descend for `t < 64` with the stated worst ratios (`0.5039` to `0.9543`), the plain classes within `<= 5` odd
steps, the multiplier classes within `<= 3` odd steps after the insertion; the composite defect labels for
`t = 0, 1` are e.g. `533, 3029` (class `27 mod 128`, `27 -> 41`, `13 * 41 = 533 -> 25`), `35, 355` (class `7 mod
64`, `5 * 7 = 35 -> 53`), `22517, 23221` (class `2047 mod 4096`, two multiplications by `11`).

| odd residues mod `4096` | density | content |
|---|---|---|
| plain descent (nine classes, `1 mod 4` to `79 mod 256`) | `51/64 = 79.69%` | actual orbit below `(76/79) x` within `<= 5` odd steps; zero defect |
| multiplier classes (seventeen) | `415/2048 = 20.26%` | one or two composite defects `[m y]`, `m in H` |
| `-1 mod 4096` | `1/2048 = 0.05%` | Lemma 2.3: multiplier `m_j = (2^j+1)/3`, then the table |

For comparison the renewal computation of the no-descent set at depth twelve (Pascal-tower note, section 5)
is `P(Bin(19, 1/2) >= 12) = 17.96%`; the paper's cruder prefix code pays `20.26%` with multipliers. **So the
universal weak coverage is: plain Collatz descent exactly where bounded-depth descent exists, and a
multiplicative shortcut exactly where it does not** -- the same cells the partition-cover theorem shows no
bounded-depth portfolio can pay.

**Theorem 6 (the `-1` cell and the `3x-1` trunk).** For odd `j` let `m_j = (2^j + 1)/3`. Then `3 m_j - 1 = 2^j`
(`m_j` is on the trunk of the `3x-1` map, `U_-(m_j) = 1`), and for every `x = -1 mod 2^j`,
`3 m_j x + 1 = 2^j x + (x + 1)` is divisible by `2^j`, so `U(m_j x) = oddpart(x + (x+1)/2^j)`: one odd step
from `m_j x` lands at `x (1 + 2^-j) + 2^-j`, a number of the size of `x` whose binary tail is `-1` only modulo
`2^(k-j)` when `x = -1 mod 2^k`. *Proof.* Direct. (Checked for `k <= 15`, odd `j <= k`, five `x` each.) □

This is what the paper's Lemma 2.2 uses (with `T`): the deepest growth cell of `3x+1` is neutralized by
multiplying with a trunk number of `3x-1`. In the receipt calculus the obligation of that step is the single
composite coordinate `[m_j x]`, and eliminating it means the actual orbit of `x` itself reaching a rooted
vertex, i.e. the plain descent of the `-1 mod 2^k` cell, which climbs for `k` steps and is exactly the cell
the partition-cover obstruction exhibits. Q1 for the universal receipts is therefore equivalent to the
descent of those cells; the multiplier families are "parameterized edge replacements" that pay them weakly
and cannot be made exact without that descent.

## 6. Multiples of three (Q2)

**Theorem 7.** Let `n = 3m` with `m` odd. Every receipt with surplus at `n` contains the edge `3m -> U(3m)`
(A2), and the deficit `U(3m)` is below `n` iff `m = 3 mod 4`. *Proof.* `U(3m) = (9m+1)/2^a` with `a = v_2(9m+1)`;
`9m + 1 = 0 mod 4` iff `m = 3 mod 4` (`9 = 1 mod 4`); then `U(3m) <= (9m+1)/4 < 3m`; if `a = 1`, `U(3m) =
(9m+1)/2 > 3m`. □ (`500,000/500,000` checks.)

So exactly half of the multiples of three are paid by their first edge; for the other half Theorem 1 applies
and the condition "all deficits below `n`" is the stopping time in odd steps, whose distribution over
`n = 3m <= 3 10^6` is `1: 50%, 2: 12.5%, 3: 12.5%, 4: 4.69%, 5: 5.47%, ...`, `6.45%` beyond ten steps, `0.12%`
beyond fifty, maximum `141` at `n = 2,252,031`: the Terras law, nothing specific to the factor three.

**An induction confined to multiples of three.** For a 3-unit `v` let `a(v)` be the least `a >= 1` with
`2^a v = 1 mod 9` (it exists and is at most `6`, since `2` generates the units mod `9`) and
`pi(v) = (2^a(v) v - 1)/3`: the least multiple of three with `U(pi(v)) = v`. If for every multiple of three
`n > 21` some `v_k = U^k(n)` has `pi(v_k) < n`, then all multiples of three are rooted (induction on `n` within
the multiples of three: `Root(pi(v_k))` gives `Root(v_k)` gives `Root(n)`; base cases `3, 9, 15, 21` by
direct routes), hence all odd integers (A2). The hypothesis is `2^a(v_k) v_k < 9m + 1`, a descent up to
`2^6` deeper than the plain one; its rank `k'` (least such `k`) has distribution `1: 1.56%, 2: 16.4%, 3: 11.7%,
4: 12.0%, ...`, `22.6%` beyond ten, maximum `144` at `n = 1,501,353`. It is a valid scheme with a harder
hypothesis, and that hypothesis is implied by Collatz (the orbit reaches `1`, and `pi(1) = 21`). Q2's "already
certified family" can be taken to be the multiples of three themselves; the obligation is unchanged.

## 7. Substitutions (Q3) and exact entry (Q4)

**Theorem 8.** A checked common-future receipt `U^u(n) = U^v(d)` (words `u`, `v`, endpoint `e`) is the signed
packet `c_u - c_v` with boundary `([n] - [e]) - ([d] - [e]) = [n] - [d]`; relative to the ideal `[n] - [1]` its
defect is `[1] - [d]`, i.e. one composite coordinate (`K_d`) when `d` is composite and none when `d` is prime
or `1`. Composition `n -> d -> d'` adds boundaries and keeps one open coordinate; a product of two
substitution packets for `n, n'` acquires `-R(n, n')` (D2) and is not a substitution. The compiled controllers
of the repo (`H`, `G`, `L`, `A`, `B`, the five-letter alphabet) are substitutions, so their symbolic defect is
always `[1] - [d(t)]` with `d(t)` affine in the parameter; the semigroup multipliers are substitutions glued
by a product, with the extra coordinate `[m y]`. *Proof.* Cancellation. □ So Q3 holds trivially, and the
difficulty it was meant to isolate is entirely the rank of the single obligation (the
[credit potential](adaptive_credit_potential_20261004.md) `E(x,k) = (x+5)(9/8)^k`), not a proliferation
of coordinates.

**Theorem 9.** Let `T_j` be the set of odd integers reaching `1` within `j` odd steps. Replacing the terminal
ray `(4^e - 1)/3 = T_1` in S1 by `T_j` (exponent vectors of length `j`), a shadow of `n` with prefix length `m`
is exact (`r = n`) iff `U^m(n) in T_j`, i.e. iff `n in T_(m+j)`. *Proof.* The construction prescribes the
prefix and lands the endpoint in the family; `r = n` says the endpoint is `U^m(n)`. □ Several rooted terminal
families therefore reduce the mismatch exactly to the extent that `n` already belongs to their union, and
exact entry for all `n` is `union_j T_j = all odd`, the conjecture. The least-exponent S1 shadow (`m = h = 1`,
Part F) is `r = 3` for `n = 27, 255, 999` (the shadow is a smaller rooted number of the right class, not a
route), `r ~ 2^12 n` for `n = 97`, `r ~ 2^11 n` for `n = 7`.

## 8. Reproduction

```bash
python 04-computation/experiments/collatz_receipt_joins_census_20261005.py 301 1001 3000000    # 3 s
python 04-computation/experiments/collatz_receipt_joins_addendum_20261005.py                   # 1 s
```
Outputs `.out` beside the scripts. Parts: A/B product and synchronous joins; C multiples of three; D depth-two
joins; E the trunk-multiplier identity; F shadow mismatch; addendum (a)-(d) join sizes, two-step fusion to
`3001`, the semigroup classes replayed (Tables 1-2 transcribed as `(residue, modulus, [(odd index, m)])`), and
the densities.

## 9. Verdicts

| claim | status |
|---|---|
| a surplus receipt contains the actual orbit; "deficits below `n`" is the stopping time; (W) is idle (Theorem 1) | PROVED |
| defect elimination is rootedness; F2 is the complete depth-one join; depth-two joins are congruence families (Theorems 2-3, Proposition 4) | PROVED + FINITE-EXACT |
| product joins are hub coincidences (`60.8%` of pairs, join points `175, 65, 35, ...`) | FINITE-EXACT |
| no synchronous two-step fusion to `3001`; conjectured for all; deeper fusions sporadic | FINITE-EXACT + CONJECTURED |
| the universal receipts' composite defects are the multiplier insertion points of the semigroup proof; `20.26%` of odd residues; the `-1` cell uses the `3x-1` trunk (Theorem 6) | CITED + PROVED + FINITE-EXACT |
| first deficit of `3m` below `3m` iff `m = 3 mod 4` (Theorem 7); the 3-parent induction and its rank | PROVED + FINITE-EXACT |
| substitutions carry one obligation (Theorem 8); exact entry is terminal-family membership (Theorem 9) | PROVED |
| Collatz, universal grounded coverage | OPEN |
