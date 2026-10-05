# The refuel bill is a codelength (2026-10-05): every solution of the base criterion is dominated by `2^-bill(b)`, the probability of `b` under the base-tree code, so the bill `log_2(32/15) + 4k(b)` is the exact information a certificate must carry per base change (not a tunable budget); the size rank pays none of it at the atomic prior's own exponent and at most half of the edges at any exponent (the valuation-one base edges with `k = 0` are pure climbs); the chain code costs `5.85 log_2 b` bits on average (`6.04` measured, `364` bits at `837799`), the Kraft sum converges to `1.031`, and the Codex kernel's `88` bases carry essentially all of it

**Session:** opus, `pascal-lyapunov-20261005` (sixth task), 2026-10-05.
Owner's directive (pasted): "The remaining global target is now explicit. Decode `U(b) = S^k(b)(G(b))`. We need
positive summable base weights satisfying `g(b) <= rho(1-r) r^k(b) g(G(b))`. For `r = 1/16, rho = 1/2` the rank must
drop by at least `log_2(32/15) + 4k(b)` at each change of base. This is the exact refuel bill. The finite kernel pays
it on its certified families; establishing suitable weights everywhere remains open."
**Inherits (cited):** the Codex [three-bit sibling-flow note](collatz_three_bit_sibling_flow_20261005.md) (sections 5-7:
bases `B = {b odd : v_2(3b+1) in {1,2}}`, the decomposition `n = S^k(b)`, the base preimage `b0(y)`, the fibre
elimination C3, the exact equivalence C4, the finite kernel of `88` bases), the Codex
[prefix-mass note](collatz_effective_prefix_mass_20261005.md) (P4 mass-one coverage, P6 discounted flow, P7-P8 the
climb register and its refuel bill), the [information-dimension bridge](collatz_information_dimension_bridge_20261005.md)
(as corrected), the [two-sheet note](collatz_two_sheet_receipts_20261005.md) (sibling ray `S(u) = 4u + 1`).

**Status: PROVED (Theorem 1, the forcing/maximality statement and its Kraft bound, elementary); FINITE-EXACT (integer
checks of the decomposition, the child map, the `(a, k)` census, the chain lengths, the kernel; all bases below
`2^20`); VERIFIED (finite numerical: the bills and Kraft partial sums are float64 sums of exact rational terms);
typed (section 5). Author-audited only; audit OWED. Collatz OPEN; this note changes the reading of the target,
not its status.**

---

## 0. What is new, in one screen

1. **The bill is a codelength (Theorem 1).** For a rooted base `b` with base chain `b = b_0 -> G(b_0) = b_1 -> ... -> b_T = 1`
   and sibling depths `k_i = k(b_i)`, put `bill(b) = sum_i [log_2(1/(rho(1-r))) + k_i log_2(1/r)]`. Then
   `g*(b) := 2^-bill(b)` solves the base inequalities with equality, is the pointwise-largest solution with
   `g(1) <= 1`, and `sum_B g* <= 1/(1-rho)` (Kraft). Hence **every** solution of the criterion satisfies
   `g(b) <= 2^-bill(b)`: the weights cannot be chosen more generously than the chain's own code probability, and
   "establishing suitable weights everywhere" is literally "every base has a finite chain", with the weights forced
   up to downward slack. The bill is the exact information a certificate must carry per base change; it is not a
   budget that a cleverer rank could lower.
2. **The code and its machine.** The map `b -> (k_0, ..., k_(T-1))` is injective (each base has one child per
   admissible depth, `child_k(c) = b0(S^k c)`, missing exactly for the one residue of `k` modulo `3` with
   `3 | S^k c`); decoding runs upward from `1`. So the bill is the length of a prefix-free description of `b` by
   the self-delimiting machine that reads depths under the geometric law `(1-r) r^k` and a stop symbol of
   probability `1 - rho`; Collatz is the statement that this machine's range is all of `B`. This is the prefix-free
   form of P4's mass-one coverage, with an explicit, non-atomic code.
3. **What the bill costs against size (sections 3-4).** On a base edge `a = v_2(3b+1) in {1, 2}` and the size drop is
   exactly `log_2(b/G(b)) = a + 2k - log_2 3 + carry`. Under the Haar law of base edges (`a = 1` with probability
   `2/3`, `k` geometric with `P(k) = (3/4)(1/4)^k`; measured to four digits) the expected drop is `0.415` bits per
   base change while the codex bill costs `2.43`: `bill(b)/log_2 b -> 5.85` (measured mean `6.04`, median `5.84`,
   `364` bits at `837799`). With `g = b^-s`: at the atomic prior's own exponent `s = 2` **no** base edge is paid
   (the sibling term is exactly neutral and the base-change charge `1.093` exceeds the largest possible base drop
   `0.415`); for `s >= 4` exactly half are paid; the unpaid half are the edges `(a, k) = (1, 0)`, pure climbs, for
   every `s` and every `(r, rho)`. The critical parameters `r = 1/4`, `rho -> 1` give the bill `2 - log_2 3 + 2k`,
   whose per-edge balance against `s = 1` size is exactly `a - 2`: zero on `a = 2`, minus one bit on every
   valuation-one base edge.
4. **The kernel.** The Codex kernel's `88` bases (chains to length `42`, bills to `81.9` bits at `b = 97`) have
   `sum 2^-bill = 1.0311`, which is also the Kraft partial sum over all `393,216` bases below `2^20` to four digits:
   the kernel carries essentially the whole code mass, because `2^-bill` decays like `b^-5.85`. Mass is not coverage:
   the remaining mass is tiny while the remaining bases are all of them.

---

## 1. The objects (from the sibling-flow note, verified)

`B = {b odd >= 1 : v_2(3b+1) in {1, 2}}`; every odd `n = S^k(b)` with `b in B` unique and `k = floor((v_2(3n+1) - 1)/2)`
(`0` errors on `131,072` values); for an odd target `y` with `3` not dividing `y` the unique base preimage is
`b0(y) = (2y-1)/3` (`y = 5 mod 6`) or `(4y-1)/3` (`y = 1 mod 6`), all predecessors being `S^j(b0(y))`. For `b in B
\ {1}`: `U(b) = S^k(b)(G(b))`, `G(b) = base of U(b)`. The child-by-depth map `c -> child_k(c) = b0(S^k c)` is defined
exactly when `3` does not divide `S^k c` (one residue of `k` modulo `3` is excluded) and `G(child_k(c)) = c`,
`k(child_k(c)) = k` (`53,247` children checked, `0` errors).

## 2. Theorem 1: the maximal solution is the code, and the code is Kraft-bounded (PROVED)

Fix `0 < r, rho < 1` and write `d(b) = rho (1-r) r^k(b)`. Call `g > 0` on `B` a solution if `g(b) <= d(b) g(G(b))`
for all `b in B \ {1}` (C3's base inequalities).

**Theorem 1.** (a) If every base is rooted, `g*(b) := prod_(i < T(b)) d(b_i)` along the chain of `b` (so `g*(1) = 1`)
is a solution with equality, and every solution `g` with `g(1) <= 1` satisfies `g(b) <= g*(b)` for all `b`.
(b) `sum_(b in B) g*(b) <= 1/(1 - rho)`. (c) `bill(b) := -log_2 g*(b) = sum_i [log_2(1/(rho(1-r))) + k_i log_2(1/r)]`
is the length of a self-delimiting description of `b`: the sequence `(k_0, ..., k_(T-1))` determines `b` by
`b_(i) = child_(k_i)(b_(i+1))` upward from `1`, and the map is injective.
*Proof.* (a) Equality along the chain is the definition; for any solution, induction on the chain gives `g(b) <=
d(b_0) d(b_1) ... g(1) <= g*(b)`. (b) Consider the random walk up the base tree from `1`: at each node stop with
probability `1 - rho`, otherwise draw `k` with probability `(1-r) r^k` and move to `child_k` if it exists (else
fail). The probability of stopping at `b` is `(1 - rho) g*(b)`, and these events are disjoint, so
`(1-rho) sum g* <= 1`. (c) Uniqueness of `child_k(c)` (section 1) gives injectivity; the walk reads a
self-delimiting input. □ If some base is not rooted, no solution exists (C4), and the chain of a non-rooted base
has no finite code.

**Reading.** The criterion has no freedom upward: the only choices are how much to lower `g` below `g*`. A proof of
"suitable weights everywhere" must supply, for every base, information at least equal to its chain code -- the
orbit itself, as C4 already says ("an exact equivalence, not a solution"). What Theorem 1 adds is that no
re-parametrization (`r`, `rho`, or a different rank) changes this: the bill per base change is the code's cost
of one more step, and any rank `R = -log_2 g` satisfies `R(b) >= bill(b) - log_2 g(1)`.

## 3. The census of base edges (FINITE-EXACT, bases below `2^20`)

| `(a, k)` | measured fraction | Haar `(2/3 or 1/3) x (3/4)(1/4)^k` |
|---|---|---|
| `(1, 0)` | `0.5000` | `0.5000` |
| `(1, 1)`, `(1, 2)`, `(1, 3)` | `0.1250, 0.0313, 0.0078` | `0.1250, 0.0312, 0.0078` |
| `(2, 0)` | `0.2500` | `0.2500` |
| `(2, 1)`, `(2, 2)`, `(2, 3)` | `0.0625, 0.0156, 0.0039` | `0.0625, 0.0156, 0.0039` |

So on base edges `E[a] = 4/3`, `E[k] = 1/3`, and since `U^i(b) = S^(k_(i-1))(G^i b)`, the actual orbit valuation at
odd step `i` is `a_i + 2k_(i-1)` with mean `2`: the base chain is the odd orbit with each sibling depth moved one step
earlier. The size drop per base change is `a + 2k - log_2 3` plus the carry `log_2((1 + 1/(3b))/(1 + 2^-(a+2k)...))`,
which is `O(1/b)`; its mean is `0.415` bits, the Terras drift, so `T(b)` (the chain length = odd root time) is about
`2.41 log_2 b` (measured mean `46.75` for bases below `2^20`, maximum `194` at `837799`).

## 4. The bill against size (VERIFIED on the exact chains)

| parameters `(r, rho)` | bill per base change | `E[bill]` | `bill/log_2 b`: mean, median, 90%, max | Kraft partial sum `sum_(b <= 2^20) 2^-bill` |
|---|---|---|---|---|
| codex `(1/16, 1/2)` | `1.093 + 4k` | `2.43` | `6.04, 5.84, 8.74, 18.8` | `1.0311` (bound `2`) |
| `(1/4, 1/2)` | `1.415 + 2k` | `2.08` | `5.21, 4.99, 7.93, 17.9` | `1.1233` (bound `2`) |
| critical `(1/4, 1)` | `0.415 + 2k` | `1.08` | `2.69, 2.61, 3.83, 8.10` | `1.3848` |

Size weights `g = b^-s` pay the edge `b` iff `s log_2(b/G(b)) >= bill`. Fractions of base edges paid (first edge of
every base below `2^20`):

| parameters | `s = 2` | `s = 3` | `s = 4` | `s = 6` | `s = 8` | chains with no overdraft at any prefix (`s = 6`) |
|---|---|---|---|---|---|---|
| codex | `0.0000` | `0.3750` | `0.5000` | `0.5000` | `0.5000` | `6.9%` (median worst overdraft `-34.7` bits) |
| `(1/4, 1/2)` | `0.1250` | `0.2500` | `0.5000` | `0.5000` | `0.5000` | `8.5%` |
| critical | `0.5000` | `0.5000` | `0.5000` | `0.5000` | `0.5000` | `19.7%` (median `-7.8` bits) |

**Exact statements behind the table.** With the codex parameters and `s = 2`: the sibling term is neutral
(`2 x 2k = 4k`) and the base-change charge `log_2(32/15) = 1.093` exceeds the largest base drop `2 - log_2 3 =
0.415` (`a = 2`, `k = 0`), so no base edge is size-paid. For any parameters and any `s`, the edges `(a, k) = (1, 0)`
have drop `-0.585 + carry < 0` and are never paid: half of all base edges are climbs. The critical bill
`2 - log_2 3 + 2k` makes the `a = 2` edges break even at `s = 1` and charges exactly one bit on each `a = 1` edge
(`balance = a - 2`), so along any chain the size-rank deficit under the critical bill is the number of
valuation-one base steps. The 90th-percentile bill is `8.7 log_2 b`; the record `837799` (`194` base changes) costs
`364` bits against `19.7` bits of size.

## 5. The kernel and what "pays it on its certified families" means (typed)

The Codex kernel certifies `88` bases (all bases `<= 127` and their parents up to `2051`; chains to length `42`)
by literal checked equations `U(b) = S^k G(b)` ending at `1`, and builds `g` on the kernel from finite path
contributions. Its bills under the codex parameters run from `0` (`b = 1`) to `81.9` bits (`b = 97`, chain length
`42`), and `sum_(kernel) 2^-bill = 1.0311`, equal to the Kraft partial sum over all `393,216` bases below `2^20` to
four digits: the code mass is concentrated on the smallest bases because `2^-bill(b) ~ b^-5.85`. "Pays the bill on
its certified families" means exactly that these `88` chains are finite and each sibling family inherits its base's
chain through `S`; the infinite families are the sibling rays of finitely many bases, and the criterion's
summability is automatic on them (the Kraft bound), not an achievement.

* **Structural (PROVED):** the forcing/maximality statement; the base-tree code; the per-edge drop identity; the
  critical balance `a - 2`.
* **Observed (VERIFIED):** the Haar law of base edges to four digits below `2^20`; the bill statistics; the Kraft
  partial sums; the size-payability fractions.
* **Unchanged (OPEN):** the existence of positive summable `g` on every base, which Theorem 1 shows is the same as
  every base having a finite chain, with `g` forced to be at most the chain's code probability. No rank that is a
  function of size, residue, or any finite-depth datum can be a solution, because the `(1, 0)` edges are unpaid for
  every such rank: a solution encodes the orbit.

## 6. Reproduction

```bash
python 04-computation/experiments/collatz_refuel_bill_code_20261005.py 1048576     # 19 s
```
Output `.out` beside the script; the kernel is read from `collatz_three_bit_sibling_flow_20261005.json`.

## 7. Verdicts

| claim | status |
|---|---|
| decomposition `n = S^k(b)`, the base map, the child-by-depth map (one child per admissible depth) | FINITE-EXACT (reproduced) |
| the chain product `2^-bill(b)` is the pointwise-largest solution of the base criterion and is Kraft-summable (Theorem 1) | PROVED |
| the bill is a prefix-free codelength of `b` (base-tree machine); Collatz = totality of its range | PROVED (reading of C4/P4) |
| base-edge law `(a, k)`: Haar to four digits; `E[drop] = 0.415`, `E[codex bill] = 2.43` per base change; `bill/log_2 b` mean `6.04` | VERIFIED |
| size weights pay no base edge at `s = 2` (codex), at most half at any `s`, never the `(1,0)` edges; critical balance `a - 2` | PROVED + VERIFIED |
| kernel: `88` bases, bills to `81.9` bits, Kraft mass `1.0311` = the mass of all bases below `2^20` | VERIFIED |
| suitable weights everywhere | OPEN (= Collatz; weights forced to be at most the chain code) |
