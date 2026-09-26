# The exact price of periodic Collatz edits is a feedback-vertex problem on the de Bruijn graph; Golomb–Mykkeltveit fails for density thresholds

**Status.**
- **PROVED** (full hand proofs below; every lemma is also checked by code):
  - **Theorem 1 (Golomb's conjecture; CITED: J. Mykkeltveit, *J. Combin. Theory Ser. B* 13 (1972) 40–45; proof re-derived in §1).** The minimum number of nodes of `B(2,k)` meeting every cycle is `Z(k)`, the number of binary necklaces of length `k`. For `k >= 3` an explicit minimum set is `M_k = S* ∪ {one node of each zero-weight necklace}`, where `S*` holds the nodes at which a necklace's sine weight crosses the positive real axis. **Corollary 1.6:** `FVS_c(k) <= Z(k) - 1` for every threshold `c >= 0`.
  - **Theorem 3 (Q3, dynamics).** Let `G` be a periodic edit of level `k`: `G = v_0` (a constant) on the residues `R`, `G = T_q` elsewhere.
    - `G` has bounded-lookahead descent iff `R` meets every expanding cycle of `G_0 = B(2,k)`. The "only if" direction holds for arbitrary edit values.
    - When it holds, every positive orbit is eventually periodic, in a cycle through `v_0` or in a positive (contracting) `T_q`-cycle.
    - **Hence the minimal Haar density of a provable periodic edit of level `k` is exactly `FVS_c(k)/2^k`.**
  - **Theorem 4 (Q4, DRIFT).**
    - `q = 5`: `1/(2k) <= price_5(k) <= (Z(k)-1)/2^k <= 1/k + 2^(-k/2)` for every `k >= 1`, and `k * price_5(k) -> 1`.
    - `q = 3`: `2^(-(1-h)k)/(3k^2) <= price_3(k) <= |Bad_k|/2^k <= 2^(-(1-h)k)`.
    - So positive drift turns the exponential rate `2^(-0.05k)` into the polynomial rate `1/k`, and both prices tend to 0.
  - **Proposition 5 (Q5).** `delta_k >= FVS^odd_c(k) >= FVS_c(k) >= nu_c(k) >= N_c(k)`, in words: sign-flip distance >= odd deletion price >= deletion price >= packing number >= necklace count. This holds in the 3-cube and in the 5-cube.
  - **Theorem 2 (Q2: the Golomb analogue is false, and the gap does not vanish).**
    - (a) **Explicit counterexamples.** At `k = 3`, for every `c in [1/3, 1/2)` (including `log_5 2`): the cycles `(1), (01), (0011)` are disjoint and expanding, but `N = 2`. At `k = 5`, `c = log_3 2`: the cycles `(1), (110), (111100), (111010)` are disjoint and expanding, but `N_5 = 2`.
    - (b) **Fractional `(k+1)`-packing.** `FVS_c(k) >= ceil(Nprim_c(k+1)/2)` for every `c` and `k`.
    - (c) **For `1/2 < c < 1`**, at every `k` with `ceil(c(k+1)) = ceil(ck)`, bound (b) exceeds `N_c(k)` by the factor `1 + (2c-1)/(2(1-c)) - o(1)`, with explicit error terms.
      - For `c = log_3 2` this factor is **1.3548**, on a set of `k` of density `1 - c = 0.369`.
      - In fact `FVS_(log_3 2)(k) > N_k` at **every** such `k >= 8`, with ratio `>= 1.35` from `k = 10` on. The proof: exact counts for `k <= 60`, the explicit bound (`>= 1.358`) computed for `61 <= k <= 1000`, and an elementary estimate beyond.
      - So `FVS_(log_3 2)(k)/N_k` does **not** tend to 1.
    - (d) **For `c < 1/2`** (e.g. `q = 5`): `FVS_c(k)/N_c(k) -> 1`.
  - **Theorem 6 (closed thresholds where Golomb–Mykkeltveit does hold).** The minimum number of nodes meeting every cycle of density `>= j/k` equals the number of necklaces with at least `j` ones, for `j in {0, 1, k-2, k-1, k}` and every `k`. The sets are explicit; the case `j = k-2` uses a gap-sequence argument.
- **FINITE-EXACT.** Every value is re-verified by the runner from frozen certificates.
  - *How the certificates are checked.* The optimal set by Bellman–Ford plus an edge-checked integer potential. The lower bound by a family of pairwise disjoint expanding cycles, or by the MaxSAT solver RC2 (independent of the HiGHS MIP used in the search), or, where RC2 exceeds its time limit, by a HiGHS re-solve.
  - **Collatz, `c = log_3 2`:** `FVS = 1, 2, 2, 4, 5, 8, 12, 20, 37` for `k = 2..10`, with `FVS^odd` equal. `FVS^odd_11 = 58`.
    - `56 <= FVS_11 <= 58`, re-verified; the search's longer HiGHS run gives `>= 57`.
    - `95 <= FVS^odd_12 <= 102`, re-verified (the lower bound solver-free); the search's HiGHS bound is 96.
    - `95 <= FVS_12 <= 102`.
  - **5n+1, `c = log_5 2`:** `FVS = 2, 3, 4, 6, 9, 15, 27, 45, 84` for `k = 2..10`, with `FVS^odd` equal, and `146 <= FVS_11 <= 148`.
  - **Step functions.** `c -> FVS_c(k)` is computed on all of `[0,1)` for `k <= 9`. The equality set `{c : FVS_c = N_c}` has length `0.700, 0.638, 0.493, 0.477, 0.328` for `k = 5..9`; at `k = 9` it misses all of `[2/9, 3/4)`.
  - **Closed thresholds.** `FVS(>= j/k) = N(>= j/k)` for all `j` when `k <= 7`. It FAILS at `(k,j) = (8,4)`, `(9, 3..6)` and `(10, 4..6)`. For `j = 2` it holds for `3 <= k <= 16` (negative-axis crossing set).
  - **Min-max.** `nu < FVS` at `(q,k) = (3,7)`: `7 < 8`, via a certified fractional cover of weight `22/3`. Also at `(3,9)`: `18 <= nu <= 19 < 20`; and at `(5,9)`: `44 < 45`. Everywhere else with `k <= 10`, `nu = FVS`, certified by a packing.
  - Mykkeltveit's set has `|M_k| = Z(k)` and leaves `B(2,k)` acyclic for `3 <= k <= 16`.
  - **Dynamics.** Each optimal odd set comes with its exact lookahead `L+1` and threshold `n_0`, and a census of all `n <= n_0`. Hence **every positive integer** reaches `{1,2}` under the minimal provable edit, for `q = 3` and `k <= 11`; the same holds at `k = 12` for the best edit found (102 residues). For `q = 5` and `k <= 10`, every positive integer reaches a cycle through 1 or one of the 5n+1 cycles of 13 and 17.
  - **New lower bounds for the sign-flip distance:** `delta_11 >= FVS^odd_11 = 58` (was `52 <= delta_11 <= 72`). Likewise `delta_12 >= FVS^odd_12 >= 95` (re-proved solver-free; 96 by the search's MIP bound), against `N_12 = 70`.
- **EMPIRICAL.** `FVS^odd_c(k) = FVS_c(k)` in every case where both are exact (all `k <= 10`, `q = 3, 5`); there is no proof.
- **OPEN.**
  - `FVS_11` (57 or 58 by the search's bound); `FVS_12` and `FVS^odd_12` (in `[95, 102]`; 96 by the search); `FVS_(log_5 2)(11)` (in `[146, 148]`);
  - the exact asymptotic constant of `FVS_c/N_c` for `c > 1/2`;
  - a proof of the `j = 2` closed-threshold equality;
  - whether `FVS^odd = FVS` always;
  - `nu` at `(3,9)` (18 or 19);
  - the conjecture, supported by the shrinking equality sets, that for every fixed `c in (0,1)`, `FVS_c(k) > N_c(k)` for all large `k`.
- **REFUTED.**
  - "`FVS_c(k) = N_c(k)` for every threshold `c` and every `k`" (the density-threshold analogue of Golomb–Mykkeltveit).
  - Its closed-threshold version "`FVS(>= j/k) = N(>= j/k)` for all `j, k`".
  - The min-max version "`nu_c(k) = FVS_c(k)`".
- No HYP or THM file was created. Nothing here bears on Collatz itself (§7).

Session `collatz-procgen-20260922`, lane "mykk", 2026-09-26.
- Scripts: `04-computation/experiments/procgen_mykk_20260926_{lib,search,run,certs}.py`.
- Output: [procgen_mykk_20260926.out](procgen_mykk_20260926.out).
- Parents:
  - [THM-4474](../../01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md) (strategy cube, Theorem A);
  - the [cube-distance note](procgen_cubedist_20260925_distance_to_provability.md) (HYP-9138, `delta_k`, `Bad_k`, the necklace bound);
  - [THM-4478](../../01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md) (non-periodic edits).

## 0. The answer in brief

* **The problem.** Collatz's level-`k` parity graph is the de Bruijn graph `B(2,k)` (Terras). A periodic edit of level `k` replaces `T_q` on a set `R` of residues mod `2^k`. It is provable by bounded lookahead exactly when `R` meets every *expanding* cycle, i.e. every cycle of odd density `> c = log_q 2` (Theorem 3). So the price of provability is `FVS_c(k)/2^k`, a feedback-vertex number of `B(2,k)` restricted to dense cycles.
* **What is classical.** For *all* cycles this is Golomb's conjecture, proved by Mykkeltveit: the minimum is `Z(k)`, one node per necklace, because the necklaces are disjoint cycles and a clever "sine weight" picks one node on each so that every cycle is caught (Theorem 1). Consequently `FVS_c(k) <= Z(k) - 1 ~ 2^k/k` for every threshold.
* **The Golomb analogue fails.** The necklaces of density `> c` still give `FVS_c(k) >= N_c(k)`, but equality is false in general, already for `k = 3`, and for Collatz's own threshold at `k = 5, 7, 8, 9, 10, ...` A cycle of density `> c` may pass through only one or two nodes of dense necklaces and borrow the rest of its length from sparse ones; several such cycles can share the nodes of one dense necklace. The fractional packing by the `(k+1)`-cycles proves the ratio `FVS/N` exceeds 1.35 infinitely often for `c = log_3 2` (Theorem 2).
* **Drift.** For `q = 5` (`c < 1/2`) the dense necklaces are almost all necklaces, and Mykkeltveit's bound meets the necklace bound: the price is `(1 + o(1))/k`. For `q = 3` both bounds are `2^(-(1-h)k)` up to polynomial factors. The drift changes the rate (polynomial versus exponential), not the limit 0.
* **Flips versus deletions.** A sign flip at an odd residue keeps the node and rewires it; any flip set that proves descent must still meet every expanding cycle of `G_0`, so `delta_k >= FVS^odd_c(k)`. The exact values improve the known lower bounds for `delta_k` where it is not known (§6).

## 1. Setting, and Mykkeltveit's theorem

**Notation.**
* `B(2,k)`: nodes are words `x = x_0 x_1 ... x_(k-1)`, edges `x -> x_1 ... x_(k-1) b`. A node is *odd* iff `x_0 = 1`; `|x|` is its number of ones. `L(x) = x_1 ... x_(k-1) x_0` is the rotation (the edge with `b = x_0`).
* A closed walk of length `p` is the same as a `p`-periodic letter sequence `(z_t)`, where `z_t` is the first letter of the `t`-th node; its nodes are the `k`-windows `z_t ... z_(t+k-1)`. Its **density** is `a/p`, `a` = number of odd nodes = number of ones per period. The average window weight is `k * density`.
* **Terras.** For odd `q`, the first `k` parities of `T_q^j(n)` depend only on `n mod 2^k` and give a bijection `W : Z/2^k -> {0,1}^k`; the parity graph `G_0` (edges from `s` to the two lifts of `T_q(s) mod 2^(k-1)`) becomes `B(2,k)`, and odd residues become odd words (checked in §A of the output for `q = 3, 5`, `k <= 12`).
* A closed walk is **`c`-expanding** iff its density exceeds `c`; for `c = log_q 2` this is `q^a > 2^p`, and equality never occurs.
* `FVS_c(k)` = min `|R|` over node sets meeting every expanding cycle; `FVS^odd_c(k)` = the same with `R` inside the odd nodes; `nu_c(k)` = max number of pairwise node-disjoint expanding cycles; `N_c(k)` = number of necklaces (rotation classes of `k`-words) of density `> c`; `Nprim_c(n)` = number of *primitive* necklaces of length `n` of density `> c`; `Z(k) = (1/k) sum_(d|k) phi(d) 2^(k/d)`.
* A closed walk of density `> c` contains a simple cycle of density `> c` (a closed walk is a union of simple cycles, and its density is a mediant of theirs). So "meet every expanding cycle" and "meet every expanding closed walk" are the same requirement.
* **Trivial chain.** The necklaces are node-disjoint cycles (a necklace of period `d` is a simple cycle of length `d`), so `N_c(k) <= nu_c(k) <= FVS_c(k) <= FVS^odd_c(k)`.

**Closed thresholds.** For `0 <= j <= k`, `FVS(>= j/k)` denotes the minimum number of nodes meeting every cycle of density `>= j/k`. For `c` slightly below `j/k` (above every simple-cycle density that is `< j/k`), the cycles of density `> c` are exactly the cycles of density `>= j/k`, so `FVS_c = FVS(>= j/k)` there; and `N_c` = number of necklaces with at least `j` ones.

### Theorem 1 (Golomb's conjecture; Mykkeltveit 1972)

Let `k >= 3`, `omega = exp(2 pi i/k)`, and for a node `x` put

    w(x) = sum_(j=0)^(k-1) x_j omega^j,     s(x) = Im w(x) = sum_j x_j sin(2 pi j/k).

Let `S* = { y : s(y) <= 0 < s(L^(-1) y) }` (equivalently: `w(y) != 0` and `arg w(y) in (-2 pi/k, 0]`). Call a necklace *zero-weight* if `w = 0` on it. Then `M_k = S* ∪ {one arbitrary node of each zero-weight necklace}` meets every cycle of `B(2,k)`, and `|M_k| = Z(k)`. Hence the minimum feedback vertex set of `B(2,k)` has exactly `Z(k)` nodes (also for `k = 1, 2`, by inspection).

*Proof.*
* **(1.1) One node per necklace.** Along a necklace `w` is multiplied by `omega^(-1)`: `w(Lx) = omega^(-1) w(x)`. If `w != 0` on it, the necklace has period `k` (a word of period `d < k` has `w = sum_(j<d) x_j omega^j * sum_(i<k/d) omega^(di) = 0`), and its `k` values are `k` points equally spaced on a circle, exactly one of which lies in the half-open sector `(-2 pi/k, 0]`. A zero-weight necklace contains no node of `S*`. So `|M_k| = Z(k)` and `M_k` has one node on every necklace.
* **(1.2) Edge rule.** For an edge `x -> y = x_1 ... x_(k-1) b`: `w(y) = omega^(-1) (w(x) + e)` with `e = b - x_0 in {-1, 0, 1}`. (Drop `x_0`, shift, append `b` at position `k-1`, where `omega^(k-1) = omega^(-1)`.) So an edge is a horizontal translation by `e` followed by a clockwise rotation by `2 pi/k`. Translations do not change `Im`.
* **(1.3) Zero mean.** For every closed walk with letter sequence `z` of period `p`: `sum_t w(x^(t)) = sum_j omega^j sum_(t<p) z_(t+j) = a * sum_j omega^j = 0`. In particular `s` has mean 0 on every closed walk.
* **(1.4) Crossing.** If `x -> y` is an edge with `s(x) > 0 >= s(y)`, then `y in S*`. Indeed `v = w(x) + e` has `Im v = s(x) > 0`, so `arg v in (0, pi)` and `arg w(y) = arg v - 2 pi/k in (-2 pi/k, pi - 2 pi/k)`; as `s(y) <= 0`, `arg w(y) in (-2 pi/k, 0]`.
* **(1.5) Flat walks.** If a closed walk avoids `S*` and `s = 0` on all its nodes, then `w` is real, and real `w > 0` would put the node in `S*` (argument 0). So `w <= 0` is real. For an edge `x -> y` of the walk, `w(y) = omega^(-1) (w(x) + e)` is real while `omega^(-1)` is not (`k >= 3`), so `w(x) + e = 0` and `w(y) = 0`. Every node is the head of an edge of the walk, so `w = 0` everywhere, `e = 0` on every edge, i.e. `b = x_0`: the walk goes round a zero-weight necklace.
* **Conclusion.** Let `C` be a cycle avoiding `M_k`. If `s` is not identically 0 on `C`, it takes both signs by (1.3); walking along `C` from a node with `s > 0`, the first edge into `s <= 0` lands in `S*` by (1.4). So `s = 0` on `C` and, by (1.5), `C` is a zero-weight necklace; but `M_k` contains a node of it. The lower bound `Z(k)` is the necklace packing. ∎

**Remarks.**
* The proof uses only the edge rule (1.2) (position 0 carries the real weight 1) and the zero sum (1.3). Crossing the negative real axis instead gives a second minimum set `S** ∪ {...}`, `S** = { y : s(y) >= 0 > s(L^(-1) y) }` (arg in `(pi - 2 pi/k, pi]`).
* Mykkeltveit's embedding is `sum x_i omega^i`, as recalled by Zheng–Kingsford–Marçais (arXiv:2001.06550, Definition 3). Their "Mykkeltveit set" (Definition 4 and Lemma 9) is the negative-axis version in a convention shifted by one factor of `omega`, with the same monotonicity argument. They also prove that the longest path in `B(2,k) - M_k` has `O(k^3)` nodes (their Theorem 3), which bounds the lookahead of the full Mykkeltveit edit (§2).
* Champarnaud–Hansel–Perrin (*Internat. J. Algebra Comput.* 14 (2004) 241–251) give a second construction of a minimum set. Golomb's conjecture is from S. W. Golomb, *Shift Register Sequences* (Holden-Day, 1967).
* Verified by the runner for `3 <= k <= 16`: `|M_k| = Z(k)`, one node per necklace, `B(2,k) - M_k` acyclic (Kahn), the crossing lemma on every edge. For prime `k` the only zero-weight necklaces are `0^k` and `1^k`. `S*` is closed under setting the first or the last letter to 1; it contains even nodes (28% at `k = 16`).

**Corollary 1.6.** For every `c >= 0` and `k >= 1`: `FVS_c(k) <= Z(k) - z_c(k)`, where `z_c(k) >= 1` is the number of zero-weight necklaces of density `<= c` (always including `0^k`). *Proof:* by (1.5), a cycle avoiding `S*` and the chosen nodes of the zero-weight necklaces of density `> c` is one of the other zero-weight necklaces, of density `<= c`. ∎

## 2. Theorem 3 (Q3): the exact price of a periodic edit

**Periodic edits.** Fix odd `q >= 3`, `c = log_q 2`, `k >= 1`, a set `R` of residues mod `2^k`, and an integer `v_0 >= 1`. Let

    G(n) = v_0          if n mod 2^k lies in R,
    G(n) = T_q(n)       otherwise            (n >= 1),

with `T_q(n) = n/2` (`n` even) and `(qn+1)/2` (`n` odd). The Haar density of the edited set is `|R|/2^k`. "Provable" means **bounded-lookahead descent**: there are `L, n_0` such that every `n > n_0` has `G^j(n) < n` for some `1 <= j <= L`.

**Theorem 3.** The following are equivalent:
* (a) `W(R)` meets every expanding cycle of `B(2,k)`;
* (b) `G` has bounded-lookahead descent.

If (a) fails, (b) fails for **every** edit with the same `R` and arbitrary values on the classes of `R`. If (a) holds, every positive `G`-orbit is eventually periodic, `G` has finitely many cycles, and each cycle either contains `v_0` or is a positive `T_q`-cycle, whose parity word satisfies `q^a < 2^p`. Consequently

    min { |R|/2^k : some periodic edit of level k with edited residues R is provable } = FVS_c(k) / 2^k,

attained with `v_0 = 1` on a minimum set.

*Proof.* Write `x^(j)(n) = W(T_q^j(n) mod 2^k)`: it is the window `j .. j+k-1` of the parity sequence of `n`, so consecutive windows are joined by edges of `B(2,k)` (Terras). Let `H = B(2,k) - W(R)`.

* **(a) ⇒ (b).** All cycles of `H` have density `< c` (never `= c`). Let `l(x) = log(q/2)` for odd `x` and `-log 2` for even `x`, and `mu = max over simple cycles γ of H of l(γ)/|γ|` (`mu < 0`; if `H` is acyclic take `mu = -1`).
  * *Potential.* `phi(x) = max` over walks of `H` starting at `x` (the empty walk included) of `sum (l - mu)` over the walk's nodes. Every closed walk has `sum (l - mu) <= 0`, so the maximum is attained on a simple path: `0 <= phi <= D := 2^k max|l - mu|`, and `phi(x) >= l(x) - mu + phi(y)` on every edge `x -> y` of `H`. Summing along a walk `x_0 -> ... -> x_j` of `H`: `sum_(i<j) l(x_i) <= j mu + D`.
  * *Affine form.* `T_q(y) = e^(l) y + 1/2` for odd `y` and `e^(l) y` for even `y`. Hence, while the orbit follows `T_q`, `T_q^j(n) = M_j n + E_j` with `M_j = exp(sum_(i<j) l(x_i)) <= e^(j mu + D)` and `E_j = sum_(i<j, x_i odd) (1/2) exp(sum_(i<s<j) l(x_s)) <= (j/2) e^D`.
  * *Lookahead.* Let `L_0` be least with `e^(L_0 mu + D) <= 1/2`, and `n_0 = max(v_0, L_0 e^D)`. Take `n > n_0` and let `tau` be the first `i` with `T_q^i(n) mod 2^k in R`. If `tau >= L_0`, the windows `x_0 .. x_(L_0 - 1)` form a walk of `H`, so `G^(L_0)(n) = T_q^(L_0)(n) <= n/2 + (L_0/2) e^D < n`. If `tau < L_0`, `G^(tau+1)(n) = v_0 < n`. ∎(a ⇒ b)
* **(b) ⇒ (a)**, for arbitrary values on `R`. Let `γ` be an expanding cycle of `H`, with letter sequence `z` of period `p` and `a` ones. The partial sums `P(t) = sum_(s<t) (z_s log q - log 2)` satisfy `P(t+p) = P(t) + (a log q - p log 2)`, with a positive increment, and they are pairwise distinct (`log q/log 2` is irrational). Start at `t_0` = the minimiser of `P` on `[0,p)` (cycle lemma): then `P(t_0 + j) > P(t_0)` for all `j >= 1`, i.e. `q^(a_j) > 2^j` for every prefix of the shifted sequence. For any `L`, let `n` be any positive integer whose first `k+L` parities are `z_(t_0) z_(t_0+1) ...` (a full residue class mod `2^(k+L)`, by Terras). For `j <= L` the window of `T_q^j(n)` is a node of `γ`, so it is not in `W(R)` and `G^j(n) = T_q^j(n) >= (q^(a_j)/2^j) n > n` (the additive term `E_j` is `>= 0`). So no `L, n_0` work. ∎(b ⇒ a)
* **Orbits under (a).** Every `n > n_0` descends within `L_0` steps, and `G(y) <= ((q+1)/2) y + v_0`. Starting from `n`, the successive "record descents" decrease until the orbit enters `[1, n_0]`, and all values in between are at most `K n` with `K = ((q+1)/2 + v_0)^(L_0)`; after the orbit enters `[1, n_0]`, all later values are at most `K B` with `B = max(n_0, max_(m <= n_0) G(m))`. So the orbit is bounded, hence eventually periodic. The minimum of a cycle cannot descend below itself, so every cycle has a point `<= n_0`: finitely many cycles. A cycle meeting a class of `R` contains `v_0`. A cycle avoiding them is a `T_q`-cycle of positive integers; it contains an odd point `x`, and `T_q^p(x) = x` reads `x (2^p - q^a) = c_w > 0`, so `2^p > q^a`. ∎
* **The price.** By (b ⇒ a) every provable periodic edit has `|R| >= FVS_c(k)`, and a minimum set with `v_0 = 1` is provable by (a ⇒ b). ∎

**The computed lookahead and threshold.** The proof's `L_0`, `n_0` are crude. The runner computes the exact ones: a *ballot walk* of length `j` is a walk `x_0 .. x_(j-1)` of `H` whose prefix multipliers `q^(a_i)/2^i` all exceed 1. Let `L` be the maximal ballot-walk length (finite under (a)); then every `n > n_0` descends within `L + 1` steps, where `n_0` is the maximum over walks that first drop below multiplier 1 at step `j+1` of `floor(c_w/(2^(j+1) - q^(a_(j+1))))`, with `T_q^(j+1)(n) = (q^a n + c_w)/2^(j+1)` on that class. Both are computed by an exact dynamic program over (node, number of odd steps), carrying the maximal `c_w` (`c_(j+1) = q c_j + 2^j` on an odd step, unchanged on an even one).

**FINITE-EXACT (output §H).** For each optimal odd set `R` (Terras-mapped to residues) and `v_0 = 1`, the runner computes the exact `L` and `n_0`. The DP proves descent within `L+1` steps for every `n > max(n_0, 1)`; the runner also re-checks this for `n <= n_0 + 20000`. It then follows every orbit from `n <= n_0 + 20000`. By strong induction this covers **all** positive integers:

| `q` | `k` | `|R|` | `L+1` | `n_0` | every positive orbit ends in |
|---|---|---|---|---|---|
| 3 | 2–6 | 1, 2, 2, 4, 5 | 2–16 | `<= 22` | `{1,2}` |
| 3 | 7 | 8 | 92 | 639 | `{1,2}` |
| 3 | 8 | 12 | 65 | 762 | `{1,2}` |
| 3 | 9 | 20 | 181 | 1911 | `{1,2}` |
| 3 | 10 | 37 | 206 | 1964 | `{1,2}` |
| 3 | 11 | 58 | 549 | 45780 | `{1,2}` |
| 3 | 12 | 102 (best found, not known optimal) | 1083 | 46585 | `{1,2}` |
| 5 | 2–7 | 2, 3, 4, 6, 9, 15 | 1–15 | `<= 17` | a cycle through 1 (`{1}`, `{1,3}` or `{1,3,8,4,2}`) |
| 5 | 8 | 27 | 157 | 1772 | `{1,3,8,4,2}` |
| 5 | 9 | 45 | 23 | 20 | `{1,3,8,4,2}`, the 5n+1 cycles of 13 and of 17 |
| 5 | 10 | 84 | 252 | 1814 | `{1,3,8,4,2}`, the 5n+1 cycles of 13 and of 17 |

The converse is also checked. For three nodes `r` of each optimal set (`k = 4..9`), `B(2,k) - (R - r)` has an expanding cycle, and the least positive integer on its ballot periodic point mod `2^(k+40)` rises for 40 consecutive steps under the weakened edit.

**Remarks.**
* With the full Mykkeltveit set (`R = W^(-1)(M_k)`, density `~ 1/k`), `H` is acyclic. Every orbit meets `R` within (longest path of `B - M_k`) `+ 1` steps. That longest path has `3, 7, 11, 21, 27, 43, 55, 83, 89, 141, 143, 215` nodes for `k = 3..14`, about `k^2`; Zheng–Kingsford–Marçais prove `Omega(k^2)` and `O(k^3)`. So this single edit is provable **for every odd `q` at once**.
* THM-4478 prices *non-periodic* edits with a *fixed* horizon `L` (density `2^(-(1-h)L + o(L))` for `q = 3`). Theorem 3 prices *periodic* edits with an *unbounded but finite* lookahead. For `q = 3` the exponents agree (Theorem 4).

## 3. Q1: exact values

**Method** (`procgen_mykk_20260926_search.py`, library `..._lib.py`). Implicit hitting set with lazy cycle constraints:
* **Master.** Minimum hitting set of the expanding cycles found so far (HiGHS MIP, one thread). Its optimum is a lower bound.
* **Separation.** Given the master's set `R`, find expanding cycles of `B(2,k) - R`, in order:
  1. the pool of *all* simple expanding cycles of length `<= 2k+2` (at most 20–22), built from primitive necklaces of that length whose `k`-windows are distinct;
  2. minimal-length expanding cycles from an all-pairs dynamic program over walk lengths;
  3. Bellman–Ford with exact integer weights, extracting node-disjoint positive cycles.

  The weights are exact: `v * [odd] - u`, where `u/v` is the best lower approximation of `c` with denominator `<= 2^k`, so a simple cycle is expanding iff its weight is positive.
* **Speed-ups.** Every new cycle is also added reversed (word reversal is an anti-automorphism of `B(2,k)` preserving density). After each master solve a greedy repair (MaxHS-style) extends `R` until feasible, collecting cycles on the way; the repaired set is made inclusion-minimal and becomes the incumbent (upper bound).
* **Termination.** The first master optimum with no expanding cycle left is optimal.
* **Certificates** (`..._certs.py`, re-verified by the runner).
  * The optimal set `R` is checked by Bellman–Ford and an edge-checked integer potential.
  * The lower bound is checked by a packing of `FVS` pairwise disjoint expanding cycles when one was found. Otherwise the runner re-proves the minimum hitting set of a (shrunk) list of expanding cycles: with RC2 (pysat MaxSAT, independent of HiGHS) in all cases with `k <= 10`, and with a fresh HiGHS solve for `k = 11` (RC2 times out).
  * For runs stopped by their time limit, the runner checks the feasible set (upper bound). For the lower bound it checks an exact LP-dual fractional packing, which is solver-free once the weights are rationalised, plus whatever HiGHS re-proves in 300 s. The longer search's MIP bound is reported separately.

**Collatz, `c = log_3 2`.** `nu` is the size of the best packing found; `FVS` and `FVS^odd` are certified as described above.

| `k` | `N_k` | `nu` | `FVS` | `FVS^odd` | `Z(k)-1` | `FVS/N` | price `FVS/2^k` | `ceil(Nprim(k+1)/2)` |
|---|---|---|---|---|---|---|---|---|
| 2 | 1 | 1 | 1 | 1 | 2 | 1 | 0.250 | 1 |
| 3 | 2 | 2 | 2 | 2 | 3 | 1 | 0.250 | 1 |
| 4 | 2 | 2 | 2 | 2 | 5 | 1 | 0.125 | 1 |
| 5 | 2 | 4 | **4** | 4 | 7 | 2 | 0.125 | 2 |
| 6 | 5 | 5 | 5 | 5 | 13 | 1 | 0.0781 | 2 |
| 7 | 5 | **7** | **8** | 8 | 19 | 1.6 | 0.0625 | 2 |
| 8 | 6 | 12 | **12** | 12 | 35 | 2 | 0.0469 | **7** |
| 9 | 16 | **18–19** | **20** | 20 | 59 | 1.25 | 0.0391 | 9 |
| 10 | 19 | 37 | **37** | 37 | 107 | 1.95 | 0.0361 | **26** |
| 11 | 52 | `>= 53` | **56–58** (57 by the search) | **58** | 187 | 1.08–1.12 | 0.0283 | 32 |
| 12 | 70 | | **95–102** | **95–102** (96 by the search) | 351 | 1.36–1.46 | 0.0232–0.0249 | 42 |

* **`k = 11`.** `FVS^odd = 58` is exact: 7 IHS rounds, 20 minutes. Its lower bound is re-proved by HiGHS on an 878-cycle certificate, because RC2 does not finish in 5 minutes. For the all-node `FVS`, the MIP bound is 57 after two 30-minute master solves, and the upper bound is the odd set (58). The runner itself re-proves only `>= 56` in 300 s (HiGHS; the exact LP-dual packing gives 55). So `57 <= FVS_11 <= 58` by the search, `>= 56` as re-verified. `FVS_11 = 57` would be the first case with `FVS < FVS^odd`.
* **`k = 12`.** After 90 minutes (8 IHS rounds, MIP capped at 15 minutes per round), `96 <= FVS^odd_12 <= 102` by the search: the upper bound is an edge-checked feasible odd set, the lower bound the HiGHS dual bound. The runner re-proves `FVS^odd_12 >= 95` without any MIP solver, by an exact LP-dual fractional packing. The same cycles give `FVS_12 >= 95` for all nodes as well, so `95 <= FVS_12 <= 102`. `N_12 = 70`, and the pruned sign-flip set of the cube-distance note has 131 flips.

**5n+1, `c = log_5 2`.**

| `k` | `N_c` | `nu` | `FVS` | `FVS^odd` | `Z(k)-1` | `k*FVS/2^k` | `k*N/2^k` | `k(Z-1)/2^k` |
|---|---|---|---|---|---|---|---|---|
| 2 | 2 | 2 | 2 | 2 | 2 | 1.000 | 1.000 | 1.000 |
| 3 | 2 | 3 | **3** | 3 | 3 | 1.125 | 0.750 | 1.125 |
| 4 | 4 | 4 | 4 | 4 | 5 | 1.000 | 1.000 | 1.250 |
| 5 | 4 | 6 | **6** | 6 | 7 | 0.938 | 0.625 | 1.094 |
| 6 | 9 | 9 | 9 | 9 | 13 | 0.844 | 0.844 | 1.219 |
| 7 | 10 | 15 | **15** | 15 | 19 | 0.820 | 0.547 | 1.039 |
| 8 | 23 | 27 | **27** | 27 | 35 | 0.844 | 0.719 | 1.094 |
| 9 | 44 | **44** | **45** | 45 | 59 | 0.791 | 0.773 | 1.037 |
| 10 | 67 | 84 | **84** | 84 | 107 | 0.820 | 0.654 | 1.045 |
| 11 | 136 | | **146–148** | | 187 | 0.784–0.795 | 0.730 | 1.004 |

* **`k = 11` (5n+1).** A 70-minute run (6 IHS rounds, MIP capped at 15 minutes) gives `146 <= FVS <= 148`, against `N = 136` and the Mykkeltveit bound `187`. The runner re-proves 146 without any MIP solver, by an exact LP-dual fractional packing.

Bold `FVS` marks `FVS > N`. Bold `nu` marks `nu < FVS`, proved by a certified fractional cover (§4.6). A bold `(k+1)`-bound marks where that bound alone proves `FVS > N`.

## 4. Q2: is `FVS_c(k) = N_c(k)`? No — and the gap does not vanish

### 4.1 Explicit counterexamples (Theorem 2a)

* **`k = 3`, any `c in [1/3, 1/2)`** (this contains `log_5 2 = 0.4307`). Expanding necklaces: `[011]` and `[111]`, so `N_c(3) = 2`. The cycles

      (1) = {111},   (01) = {010, 101},   (0011) = {001, 011, 110, 100}

  are pairwise disjoint and have densities `1, 1/2, 1/2 > c`. So `FVS_c(3) >= 3`, and `{111, 101, 110}` works: `FVS_c(3) = 3 > 2`. By hand: a 2-set must contain `111` and one node of `[011]`; `(01)` forces `101`; then `(0011)` is missed.
* **`k = 5`, `c = log_3 2` (Collatz).** `N_5 = 2` (necklaces `[11111]`, `[01111]`). The periodic words (first letter = parity of the first step)

      (1),   (110),   (111100),   (111010)

  have densities `1, 2/3, 2/3, 2/3 > 0.6309`, and their 5-windows are pairwise disjoint:

      {11111};  {11011, 10110, 01101};  {11110, 11100, 11001, 10011, 00111, 01111};  {11101, 11010, 10101, 01011, 10111, 01110}.

  So `FVS_{log_3 2}(5) >= 4 > 2`; the value is exactly 4 (§3). The necklace `[01111]` (5 nodes) is shared by three of these cycles; each of them borrows the rest of its length from the non-expanding necklaces of weight 3.
* **Why Mykkeltveit's argument breaks.** Every cycle that is not a zero-weight necklace crosses the real axis at a node of `S*`. But the crossing node may lie on a *non-expanding* necklace while the cycle is expanding, so the restriction of `S*` to expanding necklaces is not enough; the counterexamples show that *no* one-node-per-necklace selection is.

### 4.2 The fractional packing by `(k+1)`-cycles (Theorem 2b–d)

**Theorem 2(b).** For every threshold `c` and `k >= 1`: `FVS_c(k) >= ceil(Nprim_c(k+1)/2)`.

*Proof.*
* A primitive necklace `u` of length `k+1` gives a simple cycle of `B(2,k)`: a `k`-window of the `(k+1)`-periodic sequence `u^inf` is `u` read cyclically with one letter deleted; two windows at different positions agree in `k` letters, hence (same weight) in all `k+1`, which contradicts primitivity. It is expanding iff `u` has density `> c`.
* A node `x` lies on such a cycle only if the period read from `x` is `x b`, `b in {0,1}`: at most two of these cycles pass through `x`.
* If `R` meets all of them, double counting gives `|R| >= #cycles/2`. ∎

**Theorem 2(c).** Let `1/2 < c < 1`, `k >= 2`, `a = ceil(ck)`, `T(k,a) = sum_(j>=a) C(k,j)`, `rho = (k-a)/(a+1)`. If `ceil(c(k+1)) = a`, then

    FVS_c(k)/N_c(k)  >=  k (2T(k,a) + C(k,a-1) - 2^(floor((k+1)/2)+1)) / ( 2(k+1) (T(k,a) + (k-1) 2^(floor(k/2))) ),

and `C(k,a-1)/T(k,a) >= (a/(k-a+1)) (1 - rho)`. As `k -> infinity` along such `k`, the right side tends to `1 + (2c-1)/(2(1-c))`.

*Proof.*
* `T(k+1,a) = sum_(j>=a) (C(k,j) + C(k,j-1)) = 2T(k,a) + C(k,a-1)` (Pascal).
* Non-primitive words of length `n` have a period `d <= n/2`, so there are fewer than `2^(floor(n/2)+1)`; a primitive necklace has `n` words. So `Nprim_c(k+1) >= (T(k+1,a) - 2^(floor((k+1)/2)+1))/(k+1)`.
* Burnside for rotations of the words with `>= a` ones: `N_c(k) = (1/k) sum_(r<k) Fix(r)`, `Fix(0) = T(k,a)`, and `Fix(r) <= 2^(gcd(r,k)) <= 2^(floor(k/2))` for `r != 0`.
* Ratio of consecutive binomials: `C(k,j+1)/C(k,j) = (k-j)/(j+1) <= rho` for `j >= a`, so `T(k,a) <= C(k,a)/(1-rho)`, and `C(k,a-1)/C(k,a) = a/(k-a+1)`.
* Divide, using Theorem 2(b). For the limit: `a/k -> c`, so `(a/(k-a+1))(1-rho) -> (c/(1-c))(1 - (1-c)/c) = (2c-1)/(1-c)`, while `T(k,a) >= C(k,a) >= 2^(k h(a/k))/(k+1)` with `h(c) > 1/2` swamps the error terms `k 2^(k/2)`. ∎

**For `c = log_3 2`:** the limit is `1 + 0.26186/0.73814 = 1.35475`. The condition `ceil(c(k+1)) = ceil(ck)` means `frac(ck) < 1 - c`. Since `c` is irrational, `frac(ck)` is equidistributed, so the condition holds on a set of `k` of natural density `1 - c = 0.3691`. Hence

    limsup_(k -> inf) FVS_(log_3 2)(k) / N_k  >=  1.3547,

and `FVS > N` for infinitely many `k`.

**In fact `FVS_(log_3 2)(k) > N_k` at every admissible `k >= 8`.**
* `8 <= k <= 60`: exact counts. The bound (b) beats `N_k` exactly at the admissible `k >= 8`: `k = 8` (`7 > 6`), `10` (`26 > 19`), `13` (`123 > 85`), `16` (`641 > 434`), `18, 21, 24, 27, 29, ...`. The ratio is `>= 1.368` for admissible `k >= 10` and lies in `[1.399, 1.477]` for `13 <= k <= 60`.
* `61 <= k <= 1000`: the explicit inequality of Theorem 2(c) gives ratio `>= 1.358` (output §G).
* `k > 1000`: here `a/(k-a+1) >= (c/(1-c))(1 - 1/((1-c)k)) >= 1.7048` and `1 - rho > 1 - (1-c)/c = 0.4150`, so the main term is `>= 1.3524`. The error terms are below `2^(-0.4k)`, so the ratio stays `>= 1.35`.

**Theorem 2(d).** For `c < 1/2`: `(2^k - 2^(h(c)k))/k <= N_c(k) <= FVS_c(k) <= Z(k) - 1 <= (2^k + (k-1)2^(k/2))/k`, so `FVS_c(k)/N_c(k) -> 1`, exponentially fast but at the slow rate `2^(-(1-h(c))k)` (`1 - h(log_5 2) = 0.0139`). *Proof:* each necklace has at most `k` words and at most `2^(h(c)k)` words have `<= ck` ones; the upper bound is Corollary 1.6 plus Burnside. ∎

So for `c < 1/2` the Golomb analogue holds asymptotically in ratio (not exactly: `FVS > N` at `k = 3, 5, 7, 8, 9, 10` for `c = log_5 2`, e.g. `84 > 67` at `k = 10`), while for `c > 1/2` the ratio stays away from 1 along a positive-density set of `k`.

### 4.3 Where equality holds: the complete step functions (FINITE-EXACT, `k <= 9`)

`FVS_c(k)` is a non-increasing step function of `c`. Its breakpoints are densities of simple cycles of `B(2,k)`, while `N_c(k)` only jumps at the multiples of `1/k`. The runner computes `c -> FVS_c(k)` on all of `[0, 1)` by descent: from the current value `V`, the next breakpoint is `min over |R| = V of rho_max(B(2,k) - R)`, found by alternating Karp's maximum cycle mean with exact solves of `FVS(>= rho)`. Every segment is certified:
* an upper set `R_i` with `|R_i| = V_i` and exact `rho_max(B - R_i)` equal to the segment's bottom;
* a list of cycles of density `>=` the segment's top whose minimum hitting set is `V_i` (RC2, or HiGHS where RC2 exceeds 20 s).

| `k` | `FVS_c(k)` on `[0,1)`: `[lo, hi)`: value | equality set `{c : FVS_c = N_c}` | length |
|---|---|---|---|
| 5 | `[0,1/3)`:7 `[1/3,1/2)`:6 `[1/2,2/3)`:4 `[2/3,4/5)`:2 `[4/5,1)`:1 | `[0,1/5) ∪ [1/3,2/5) ∪ [1/2,3/5) ∪ [2/3,1)` | 0.700 |
| 6 | `[0,1/4)`:13 `[1/4,1/3)`:12 `[1/3,2/5)`:11 `[2/5,3/7)`:10 `[3/7,1/2)`:9 `[1/2,4/7)`:7 `[4/7,3/5)`:6 `[3/5,2/3)`:5 `[2/3,5/7)`:4 `[5/7,3/4)`:3 `[3/4,5/6)`:2 `[5/6,1)`:1 | `[0,1/6) ∪ [1/4,1/3) ∪ [3/7,1/2) ∪ [3/5,2/3) ∪ [3/4,1)` | 0.638 |
| 7 | `[0,1/4)`:19 `[1/4,1/3)`:18 `[1/3,2/5)`:17 `[2/5,1/2)`:15 `[1/2,5/9)`:11 `[5/9,3/5)`:10 `[3/5,5/8)`:9 `[5/8,7/11)`:8 `[7/11,2/3)`:7 `[2/3,3/4)`:5 `[3/4,7/9)`:3 `[7/9,6/7)`:2 `[6/7,1)`:1 | `[0,1/7) ∪ [1/4,2/7) ∪ [2/5,3/7) ∪ [5/9,4/7) ∪ [2/3,5/7) ∪ [7/9,1)` | 0.493 |
| 8 | `[0,1/5)`:35 `[1/5,1/4)`:34 `[1/4,1/3)`:33 `[1/3,5/13)`:30 `[5/13,2/5)`:29 `[2/5,3/7)`:28 `[3/7,4/9)`:27 `[4/9,5/11)`:26 `[5/11,1/2)`:25 `[1/2,5/9)`:20 `[5/9,4/7)`:18 `[4/7,3/5)`:16 `[3/5,5/8)`:13 `[5/8,2/3)`:12 `[2/3,7/10)`:8 `[7/10,5/7)`:7 `[5/7,3/4)`:6 `[3/4,7/9)`:5 `[7/9,4/5)`:3 `[4/5,7/8)`:2 `[7/8,1)`:1 | `[0,1/8) ∪ [1/5,1/4) ∪ [1/3,3/8) ∪ [3/5,5/8) ∪ [5/7,3/4) ∪ [4/5,1)` | 0.477 |
| 9 | `[0,1/5)`:59 `[1/5,1/4)`:58 `[1/4,2/7)`:57 `[2/7,3/10)`:56 `[3/10,1/3)`:55 `[1/3,2/5)`:52 `[2/5,5/12)`:49 `[5/12,3/7)`:48 `[3/7,5/11)`:45 `[5/11,10/21)`:43 `[10/21,1/2)`:42 `[1/2,7/13)`:34 `[7/13,6/11)`:33 `[6/11,5/9)`:32 `[5/9,9/16)`:31 `[9/16,4/7)`:30 `[4/7,3/5)`:28 `[3/5,8/13)`:23 `[8/13,5/8)`:21 `[5/8,7/11)`:20 `[7/11,2/3)`:17 `[2/3,7/10)`:14 `[7/10,8/11)`:10 `[8/11,11/15)`:9 `[11/15,3/4)`:8 `[3/4,4/5)`:6 `[4/5,5/6)`:3 `[5/6,8/9)`:2 `[8/9,1)`:1 | `[0,1/9) ∪ [1/5,2/9) ∪ [3/4,7/9) ∪ [5/6,1)` | 0.328 |

(`N_c(k)` is the number of necklaces with more than `ck` ones, constant on each `[j/k, (j+1)/k)`. For `c = log_3 2 = 0.6309` read the segments containing it: `4, 5, 8, 12, 20` at `k = 5..9`. For `c = log_5 2 = 0.4307`: `6, 9, 15, 27, 45`.)

**Structure.**
* On each necklace interval `[j/k, (j+1)/k)`, `N_c` is constant. `FVS_c` starts high just above `j/k` and falls through breakpoints at the densities of short non-necklace cycles (`1/4, 1/3, 2/5, 1/2, 5/9, 3/5, 5/8, 7/11, 2/3, 3/4, 7/9, ...`).
* Equality `FVS_c = N_c` holds exactly on a final segment `[beta_j, (j+1)/k)` when that segment exists.
* It exists for every `j` when `k <= 7`, but not for `(k,j) = (8,4)`, `(9, 3..6)` and `(10, 4..6)` (§4.4).
* For prime `k = 5, 7`, the value just above `j/k` equals the number of necklaces with at least `j` ones (`1 <= j <= k-2`): crossing `j/k` from below does not lower `FVS` at all, although `N` drops by the number of weight-`j` necklaces.

### 4.4 Closed thresholds: Theorem 6 and its failure in the middle

**Theorem 6.** Let `F_j(k)` be the minimum number of nodes meeting every cycle of density `>= j/k`, and `N_(>=j)(k)` the number of necklaces with at least `j` ones. Then `F_j(k) = N_(>=j)(k)` for
* `j = 0` and `j = 1` (every `k`): `Z(k)` and `Z(k) - 1`, by Theorem 1 and Corollary 1.6;
* `j = k` (every `k`): only the loop `1^k` has density 1;
* `j = k-1` (`k >= 2`): `F = 2`, with the set `{1^k, 1^(k-1)0}`;
* `j = k-2` (`k >= 3`): `F = 2 + floor(k/2)`, with the set `R` below.

*Proof for `j = k-1`.* Avoiding `1^k` means every window has at most `k-1` ones. Density `>= (k-1)/k` means average window weight `>= k-1`. So every window has exactly `k-1` ones, hence `z_(t+k) = z_t` for all `t`: the cycle is the necklace `[0 1^(k-1)]`, which contains `1^(k-1)0`. ∎

*Proof for `j = k-2`.* Let

    R = { 1^k,  1^(ceil(k/2)-1) 0 1^(floor(k/2)) }  ∪  { 1^(k-1-d) 0 1^(d-1) 0 : 1 <= d <= floor(k/2) }

(one rotation of each necklace with at most two zeros; `|R| = 2 + floor(k/2)`). Let `C` avoid `R`.
* **Gaps.** `C` contains no window `1^k`, so its zeros exist. Label them cyclically with gaps `1 <= g_i <= k`: `g_i` is the distance from zero `Z_i` to `Z_(i+1)`. Density `>= (k-2)/k` means average gap `>= k/2`.
* **(A)** The word `1^(ceil(k/2)-1) 0 1^(floor(k/2))` occurs iff some zero has `g_(i-1) >= ceil(k/2)` and `g_i >= floor(k/2) + 1`.
* **(B)** The word `1^(k-1-d) 0 1^(d-1) 0` occurs iff `g_i = d` and `g_(i-1) >= k - d`.
* **Big and small gaps.** Call a gap *big* if `> floor(k/2)`, *small* otherwise. If `g_i` is big, its successor `g_(i+1)` is small: otherwise avoiding (A) at `i+1` gives `g_i <= ceil(k/2) - 1 <= floor(k/2)`. Then avoiding (B) at `i+1` gives `g_i + g_(i+1) <= k - 1`.
* **Counting.** Pair each big gap with its (small) successor. With `P` pairs among `r` gaps, `sum g <= (k-1)P + floor(k/2)(r - 2P) <= (k/2) r - P`.
* **Conclusion.** An average `>= k/2` forces `P = 0` and all gaps equal to `k/2`. But then (B) with `g_(i-1) = g_i = k/2` gives `k <= k-1`, a contradiction. (With a single zero, `g = p` is its own predecessor; the same two rules give `2g <= k-1`.) So every cycle avoiding `R` has density `< (k-2)/k`.

The lower bounds are the disjoint necklaces. ∎

**FINITE-EXACT.**
* `j = 2`: `F_2(k) = Z(k) - 2` for `3 <= k <= 16`. The set is the negative-axis crossing set `S**` restricted to the necklaces with at least two ones (plus their zero-weight necklaces). I have no proof. The positive-axis set `S*`, restricted the same way, fails for every `4 <= k <= 16`.
* All `j` for `k <= 7`: `F_j(k) = N_(>=j)(k)`.
* **The equality fails in the middle.**
  * `F_4(8) = 25 > 23`;
  * `F_j(9) = 55, 45, 32, 17 > 54, 44, 30, 16` for `j = 3, 4, 5, 6` (from the step functions of §4.3);
  * `F_j(10) = 91, 72, 45 > 89, 67, 41` for `j = 4, 5, 6` (separate runs). At `k = 10` it holds for `j = 3` (101) and `j = 7, 8, 9` (19, 7, 2).

  So even the closed-threshold version of the Golomb–Mykkeltveit statement is **REFUTED**. It fails on a growing middle range of `j`.

### 4.5 Mykkeltveit-type selections (what fails)

One node per expanding necklace can only work where `FVS = N`. Among the natural rules (the tests are in the runner, §G3):
* `S*` or `S**` restricted to expanding necklaces;
* the lexicographically least or greatest *ballot* rotation (all prefixes expanding: the height-walk / cycle-lemma choice);
* the lexicographically greatest rotation (a Christoffel-type choice).

At `c = log_3 2` some rule works exactly at `k = 3, 4, 6`: all five rules at `k = 3`; `S*` and ballot-lexmin at `k = 4`; `S**` and ballot-lexmin at `k = 6`. No one-per-necklace rule can work at `k = 5, 7, 8, 9, 10`. At `c = log_5 2` a rule works only at `k = 4`; at `k = 6` (where `FVS = N = 9`) none of these rules does. The ballot-lexmin rule is the most robust: it is the set of Theorem 6 at `(k-2)/k`, and it works at every closed threshold `(k-2)/k` and `(k-1)/k`.

### 4.6 Min-max: even `nu = FVS` fails

Mykkeltveit's theorem is a min-max statement: the packing number (`Z(k)` disjoint necklaces) equals the covering number. For density thresholds even this fails:
* `q = 3, k = 7`: a certified fractional cover of total weight `22/3` shows `nu_c(7) <= nu* = tau* <= 22/3 < 8`, and a packing of 7 cycles exists. So `nu = 7 < 8 = FVS`.
* `q = 3, k = 9`: `18 <= nu <= 19 < 20 = FVS` (cover of weight `6989/360 = 19.41`). A scratch CP-SAT search (EMPIRICAL; not part of the runner) found no 19-packing among the 94,418 simple expanding cycles of length `<= 24` plus the certificate cycles (18 is optimal within that pool), so probably `nu = 18`.
* `q = 5, k = 9`: `nu = 44 < 45 = FVS`.
* In every other computed case with `k <= 10` (`q = 3, 5`) the packing found equals `FVS`, which is then certified by the packing alone. At `k = 11` (`q = 3`) the packing found has 53 cycles against `FVS^odd = 58`, and `nu` is not determined.

*Certification of a cover.* A fractional cover `x` (integer units over a common denominator `D`) is valid iff no expanding cycle has fewer than `D` units. Let `P` be the support of `x`.
* A positive cycle inside `B - P` is found by Bellman–Ford.
* Otherwise, longest walks through `B - P` between nodes of `P` give a max-plus matrix, and a dynamic program over (node of `P`, units used `<= D-1`) finds any positive closed walk. Its densest simple sub-cycle is then an expanding cycle of `x`-weight `< 1`.

The runner re-runs this exact separation for every stored cover (§G2).

## 5. Q4 (DRIFT): polynomial price for `q = 5`, exponential for `q = 3`

**Theorem 4.** Let `price_q(k) = FVS_(c_q)(k)/2^k`, the minimal Haar density of a provable periodic edit of level `k` (Theorem 3).
* **(a) `q = 5`** (`c = log_5 2 = 0.430677 < 1/2`, positive drift). For every `k >= 1`,

      1/(2k)  <=  N_c(k)/2^k  <=  price_5(k)  <=  (Z(k)-1)/2^k  <=  1/k + 2^(-k/2),

  and `N_c(k) >= (2^k - 2^(h(c)k))/k` with `h(c) = H_2(log_5 2) = 0.986083`. Hence `k * price_5(k) -> 1`: **`price_5(k) = (1 + o(1))/k`**.
* **(b) `q = 3`** (`c = log_3 2 > 1/2`, negative drift). For every `k >= 2`,

      2^(-(1-h)k) / (3k^2)  <=  N_k/2^k  <=  price_3(k)  <=  |Bad_k|/2^k  <=  2^(-(1-h)k),     h = H_2(log_3 2) = 0.9499555,

  so `log_2 price_3(k) = -(1-h)k + O(log k)`, with `1 - h = 0.0500445`.

*Proof.*
* (a), lower bounds. A necklace has at most `k` words, so `N_c(k) >= (1/k) #{words with more than ck ones}`. Since `ck < k/2 <= ceil(k/2)`, all words with at least `ceil(k/2)` ones count, and by the symmetry `j <-> k-j` they are at least half of all words: `N_c(k) >= 2^(k-1)/k`. Also at most `2^(h(c)k)` words have `<= ck` ones (`c <= 1/2`, entropy bound).
* (a), upper bounds. Corollary 1.6 (`k >= 3`; `k = 1, 2` directly: `FVS = 1, 2 = Z - 1`). Burnside: `Z(k) = (1/k) sum_(d|k) phi(d) 2^(k/d) <= (2^k + (k-1) 2^(k/2))/k`, because `2^(k/d) <= 2^(k/2)` for `d >= 2` and `sum_(d|k, d>=2) phi(d) = k - 1`.
* (b). Lower bound: the necklace bound and `N_k >= 2^(hk)/(3k^2)` (cube-distance note, Theorem 2; re-checked for `k <= 12`). Upper bound: let `Bad_k` be the words all of whose prefixes `j = 1..k` satisfy `3^(a_j) > 2^j`. Every expanding cycle contains a node of `Bad_k`: start it at its ballot rotation (cycle lemma, as in Theorem 3). And `Bad_k` lies inside the words with more than `ck` ones, at most `2^(hk)` of them (`c > 1/2`). All words of `Bad_k` are odd. ∎

**Reading.**
* **The drift changes the rate, not the limit.** `price_3(k) = 2^(-0.050 k + O(log k))` and `price_5(k) ~ 1/k`. Both tend to 0: in both cubes the provable periodic edits accumulate, in Haar measure, on the unedited map.
* **Where each bound comes from.** For `q = 3` the sparse `Bad_k` (the undecided classes) already suffices, and the dense necklaces are exponentially rare. For `q = 5` almost every necklace is expanding (`N_c/Z -> 1`), so one node per necklace is both necessary (packing) and, up to `1 + o(1)`, sufficient (Mykkeltveit). `Bad_k(5)` is useless there: its density tends to about 0.176.
* **Exact values** (§3). `k * price_5(k) = 0.94, 0.84, 0.82, 0.84, 0.79, 0.82` for `k = 5..10` and `0.78–0.80` at `k = 11`, still below the limit 1. The upper bound `k(Z(k)-1)/2^k` is already `1.045` at `k = 10`, but the necklace bound `k N_c(k)/2^k = 0.654` converges slowly: the words with at most `ck` ones make up `P(Bin(k,1/2) <= ck)`, which decays only like `2^(-0.0139 k)`.

**Comparison with the sign-flip data for 5n+1** (cube-distance note §6: pruned class-(i) flip sets with `k * delta/2^(k-1) ~ 3.2` for `7 <= k <= 14`, exact `delta_7 = 29`).
* Proposition 5 gives the rigorous floor `delta^(5)_k >= FVS^odd_(c_5)(k) >= N_(c_5)(k) >= 2^(k-1)/k`, i.e. `k * delta/2^(k-1) >= 1` for every `k`, and `>= 2 - o(1)` asymptotically. The deletion price explains the order `1/k`.
* It does **not** explain the constant. In the odd-residue normalisation, `k * FVS^odd/2^(k-1) = 1.64, 1.69, 1.58, 1.64` for `k = 7..10` (tending to 2 if `FVS^odd ~ FVS`), against `3.17` for the exact flip distance at `k = 7` (`29` flips versus `15` deletions) and about `3.1–3.3` for the pruned flip sets.
* The extra cost has a clear source. A flip does not delete a node; it re-routes the node to the `5n-1` branch, which can close new expanding cycles. That is why no class-(i) strategy exists at all for `k <= 6`, while deletion sets always exist.
* Whether the flip distance of 5n+1 tends to 0 remains open (another lane); deletions tend to 0 like `1/k` (Theorem 4a).

## 6. Q5: flips versus deletions

**Proposition 5.** In the `q`-cube (`T_sigma(n) = (qn + sigma(n mod 2^k))/2` on odd `n`), let `delta_k` be the least number of sign flips (relative to `qn+1`) giving a class-(i) strategy. Then

    delta_k  >=  FVS^odd_c(k)  >=  FVS_c(k)  >=  nu_c(k)  >=  N_c(k).

*Proof.* Let `F` be the set of flipped (odd) residues of a class-(i) strategy `sigma`. The edges of the parity graph out of a residue `s` depend only on `T_sigma(s) mod 2^(k-1)`, and `T_sigma = T_q` at every even residue and every odd residue outside `F`. So every cycle of `G_0` avoiding `F` is a cycle of `G_sigma` with the same parities. Class (i) forbids expanding cycles in `G_sigma` (THM-4474 Theorem A, whose proof uses only that `q` is odd). Hence `W(F)` is an odd node set meeting every expanding cycle of `B(2,k)`: `|F| >= FVS^odd_c(k)`. The other inequalities are §1. ∎

The cube-distance note proved `delta_k >= N_k` by the same argument restricted to necklaces. The link to `FVS^odd` is strictly stronger wherever `FVS > N`.

**Collatz (`q = 3`).** `delta_k` is taken from the cube-distance note; the other columns are from §3.

| `k` | `delta_k` | `FVS^odd` | `FVS` | `nu` | `N_k` | `delta > FVS^odd` | `FVS > N` |
|---|---|---|---|---|---|---|---|
| 2 | 1 | 1 | 1 | 1 | 1 | no | no |
| 3 | 2 | 2 | 2 | 2 | 2 | no | no |
| 4 | 2 | 2 | 2 | 2 | 2 | no | no |
| 5 | 4 | 4 | 4 | 4 | 2 | no | **yes** |
| 6 | 5 | 5 | 5 | 5 | 5 | no | no |
| 7 | 9 | 8 | 8 | 7 | 5 | **yes** | **yes** |
| 8 | 14 | 12 | 12 | 12 | 6 | **yes** | **yes** |
| 9 | 23 | 20 | 20 | 18–19 | 16 | **yes** | **yes** |
| 10 | 40–44 | 37 | 37 | 37 | 19 | **yes** | **yes** |
| 11 | 52–72 → **58–72** | 58 | 56–58 | `>= 53` | 52 | ? | **yes** |
| 12 | `<= 131` (pruned) → **`>= 95`** | 95–102 | 95–102 | | 70 | ? | **yes** |

* **Where each inequality is strict.**
  * `delta_k > FVS^odd_k` at `k = 7, 8, 9, 10`; equality for `k <= 6`.
  * `FVS^odd_k = FVS_k` wherever both are known.
  * `FVS_k > N_k` at `k = 5` and at every `k >= 7` computed; equality at `k = 2, 3, 4, 6`.
  * At `k = 5` the whole excess `delta_5 - N_5 = 2` is already forced by deletions. At `k = 7, 8, 9, 10` the deletions force `3/4, 6/8, 4/7` and `>= 18/25` of the excess `delta_k - N_k`.
* **New lower bounds.** `delta_11 >= FVS^odd_11 = 58` improves the cube-distance note's `52`. `delta_12 >= FVS^odd_12 >= 95` (solver-free; 96 by the search's MIP bound) improves `N_12 = 70`.
* **5n+1.** `delta^(5)_7 = 29 >= FVS^odd_7 = 15` (§5). The level-`<= 6` emptiness of the 5n±1 provable class has no deletion analogue: deletion sets always exist.

## 7. What this does and does not say

* **What it does.**
  * It identifies the price of provability by periodic edits with a classical object: for every odd `q` it is a feedback-vertex number of the de Bruijn graph, restricted to cycles denser than `log_q 2`. The price is exact (Theorem 3), and exactly computable for small `k`.
  * It explains the cube-distance lane's lower bound `delta_k >= N_k` as the necklace half of a Golomb-type problem, and it shows that the other half fails. `FVS_c > N_c` is typical, and for `q = 3` the ratio exceeds 1.35 on a positive-density set of levels. The mechanism is explicit: dense cycles that pass through only one or two nodes of dense necklaces.
  * It separates the two drifts cleanly. The deletion price is `(1 + o(1))/k` for `5n+1` and `2^(-(1-h)k + O(log k))` for `3n+1`; both tend to 0.
* **What it does not do.**
  * Nothing here is about the unedited map. Collatz stays in class (iv) of THM-4474 at every level: the loop at `-1` (word `1^k`) is an expanding cycle of `G_0` that no unedited map can remove.
  * The edits are *periodic*, and the lookahead depends on `k` and `R` (`L + 1 = 181` at `q = 3, k = 9`, and 1083 at `k = 12`). THM-4478 handles non-periodic edits with a fixed horizon; the two results have the same exponent for `q = 3`, but neither implies the other.
  * `FVS_11` (all nodes), `FVS_12`, `FVS^odd_12` and the `q = 5, k = 11` value are bounds, not exact values (§3). `FVS^odd = FVS` and the `j = 2` closed-threshold equality are empirical.
* **Failure modes and caveats.**
  * The IHS master (a set-cover MIP) becomes the bottleneck at `k = 11` (minutes per solve, LP gap about 2).
  * RC2, independent of HiGHS, re-proves most lower bounds for `k <= 10` in seconds. Where it exceeds its time limit, the runner re-solves with HiGHS: 8 step-function segments, `k = 11`, and possibly the largest `k = 10` closed-threshold certificates (the output does not record which solver closed those). For the time-limited runs it also checks a solver-free LP-dual fractional packing, which is weaker.
  * The step functions were computed only for `k <= 9`.
  * The cycle-pool enumeration keeps only simple cycles, and a closed walk is used only through its simple sub-cycles, which is justified in §1.

## 8. Reproduction

```bash
python3 04-computation/experiments/procgen_mykk_20260926_run.py > 05-knowledge/results/procgen_mykk_20260926.out
```

* **Runner.** It re-verifies everything from the frozen certificates in `procgen_mykk_20260926_certs.py` (a data module: base64 of gzipped JSON) and prints only `check(...)`-guarded claims. It needs `numpy`, `highspy` (HiGHS), `pysat` (RC2) and `ortools`. `ortools` must be imported first: it ships its own HiGHS, and loading `highspy` first breaks `ortools`' symbol resolution.
* **Searches.** These are not repeated by the runner. Results go to `scratch/procgen_mykk/results/`, then are frozen:

  ```bash
  S=04-computation/experiments/procgen_mykk_20260926_search.py
  python3 $S stepfun:2 stepfun:3 stepfun:4 stepfun:5 stepfun:6 stepfun:7 stepfun:8 stepfun:9
  python3 $S log3:2 log3odd:2 ... log3:10 log3odd:10 log5:2 log5odd:2 ... log5:10 log5odd:10
  python3 $S "log3odd:11@mip=900@time=10800" "log3odd:12@mip=900@time=5400" "log5:11@mip=900@time=3600"
  python3 $S "log3:11@mip=1800@time=7200"      # stopped after 2 rounds; then: python3 $S bounds:log3:11
  python3 $S "ge:3/10:10@time=900" ... "ge:9/10:10@time=900"
  python3 $S cover:3:2 ... cover:3:9 cover:5:2 ... cover:5:9
  python3 $S freeze
  ```

  Instances are checkpointed (`scratch/procgen_mykk/ckpt_*.json`) and resumable; `@mip=` caps each master MIP, `@time=` the whole run (bounds are then recorded).
* **Cost.**
  * Runner: wall 1859 s on the shared 8-core machine. Peak RSS 513 MB (`ru_maxrss`; `/usr/bin/time -l`: 538 MB max resident, 384 MB peak footprint). One process; HiGHS on one thread; RC2 in forked, time-limited subprocesses.
  * Searches, run one or two at a time, each under 540 MB RSS:
    * small instances `k <= 10`: minutes;
    * step functions `k <= 9`: 22 min;
    * `log3odd:11`: 36 min;
    * `log3:11`: 60 min, stopped, bounds taken from its checkpoint;
    * `log3odd:12`: 93 min, time-limited;
    * `log5:10`: 13 min;
    * `log5:11`: 69 min, time-limited;
    * closed thresholds at `k = 10`: 65 min;
    * covers and packings: under a minute each.
* **SHA-256** (raw bytes):
  * `procgen_mykk_20260926_lib.py` `2b5ada6a2dd7dcfbaaca6bb62ce16be3119f3f34934e6fedbc8d2eb735771eb6`
  * `procgen_mykk_20260926_search.py` `6b6ba966c0f81a830560b1a5f08dc5dee8304f5cad1e7dcb814c4add8cec8772`
  * `procgen_mykk_20260926_run.py` `c3de1c5579cc352c417e94af825644b8d599254677a12a112859600780502c75`
  * `procgen_mykk_20260926_certs.py` `60895e01b664d8ae6929100b52a42a16057358270ea39b40566c2f9efb7e8eb0` (frozen certificates, 166 KB)
  * `procgen_mykk_20260926.out` `e5d870f1c942b36241437051b905cc4a98f1fc26d5e53d61b43e9adcf3aa9ea0`
