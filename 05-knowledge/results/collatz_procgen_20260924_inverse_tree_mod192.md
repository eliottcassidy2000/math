# The owner's inverse tree mod 192: an exact type automaton, a proof that residue types cannot prove Collatz, and what the types do prove

**Status.**
* **PROVED** (hand proofs below; every finite ingredient is re-checked by the scripts):
  1. **Theorem 1 (equivalence).** The closure `R` of the root cycle `{1,2}` under `D(x) = 2x` and `E(x) = (2x-1)/3` (legal iff `x = 2 mod 3`) is exactly the `T`-basin of `{1,2}`. So Collatz holds iff `R = Z_{>0}`, and `R` is a unicyclic functional graph.
  2. **Proposition 2 (types).** `D`, `E` and `T` act on types `n mod 3*2^k` with an exact precision flow. A residue mod 192 decides `E(n) mod 128`, `D(n) mod 384`, the parity window `x_0..x_5`, the mod-3 types of `T^0(n), ..., T^6(n)` and `T^6(n) mod 3^(1+a)`. It does **not** decide `E(n) mod 3`, which needs `n mod 9`; the full type `E(n) mod 192` needs `n mod 288`.
  3. **Proposition 3 (the automaton, every `k`).** `A_k` has one strongly connected component, the `2^(k+1)` non-multiples of 3, plus `2^k` singletons (the multiples of 3; only `{0}` carries a loop). Its spectral radius is `4/3`, with Perron vector `(1, 4/3)` on the classes `1, 2 mod 3`.
  4. **Proposition 5 (order in the window).** For `n > 979` the whole order pattern of `T^0(n), ..., T^6(n)` is a function of `n mod 64`, and it is the same on both sheets.
  5. **Theorem 6 (transport) and Corollary 7 (SHEET no-go).** `x -> -x` is an isomorphism of the entire type structure of `3n+1` onto that of `3n-1`. This covers all residues at once, Haar measure, the integer predicate, every `l`-adic absolute value and the archimedean **size** `|x|`. Hence every proof of Collatz uses a hypothesis about `T` that changes truth under `x -> -x`, i.e. one that refers to the fixed positive cone. In exact form this is **side-awareness**: residues, densities, integrality and `l`-adic data are side-blind, so they prove for the positive half exactly what they prove for the negative half, where `T_+` has three more cycles.
  6. **Proposition 8 (DEFECT).** Two explicit planted maps agree with `T` off a density-zero set of starting points and have divergent orbits:
     * `T1(3*2^m) = 3*2^(m+1)`;
     * `T2`, which agrees with `T` mod `6^i` at its `i`-th planted point.
  7. **Propositions 9 and 10 (exact Haar identities).**
     * `E[N_(d,e)] = C(d,e) 3^-e`.
     * The averaged size series is `1/(1 - 2^-s - (3/2)^s/3)`. It has poles at `s = 1, 2`, and its residue at `s = 1` is `1/log(2/sqrt3) = 6.9521`, the reciprocal of the drift. The residues are `0, R, 2R` on the classes `0, 1, 2 mod 3`.
     * So the Haar-averaged, carry-free model tree of a root `a` of class `c` has density `c R/a` (Wiener–Ikehara). This is a statement about the model, not about any single integer tree (§3.2).
  8. **Proposition 11 (self-similarity).**
     * `E D^2 = S E` with `S(p) = 4p+1`, and `3 S^j(p) + 1 = 4^j (3p+1)`.
     * The predecessor tree splits exactly into a doubling chain and the trees of one sibling ladder.
     * Every ladder converges 2-adically to `-1/3 = E(0)`, and mod 192 ends in the 3-cycle `{21, 85, 149}`.
  9. **Proposition 12 (incongruent classes).** The owner's mod-3 refinement removes about half of the undecided classes, but not the exponent: the refined exceptional set still has dimension `h(log_3 2) = 0.9500`.
  10. **Proposition 13 (the trap is SHEET).** The 16-point trap `F` of `E_{6 mod 8}` is `-(orbit of 1 under the 3n-1 relaxation)`. It contains the 3n-1 cycles `{1,2}` and `{5,14,7,20,10}`. The full E-orbit of `-1` also meets the `-17` cycle.
  11. **Proposition 14 (renewal identity).** `N_R(X) = N_R(X/2) + #{n in R cap G : E(n) <= X}` holds exactly for every basin. Its main term is neutral, so every density solves it.
  12. **Theorem 15 (reformulation).** Given the (proved) sign law, `Collatz AND (no divergent 3n-1 orbit)` is equivalent to the sheet-symmetric statement `(i) AND (ii)` about all integers. Collatz itself is the side-relative statement `H(T)` of §2.2. Counting by size `|x|` can never separate the sheets.
* **FINITE-EXACT:**
  * SCC tables (`k <= 10`); the mod-192 decision tables; the `nu`-isomorphism of the automata (`k <= 12`).
  * The 3n-1 basins at `10^7`: `0.32738, 0.32450, 0.34812`. Each is residue-equidistributed mod 3, 8, 9 and 192 (chi-square consistent with uniform).
  * `Bad_+ cap [-10^7, 10^7] = {-17, -5, -1}`, and its mirror `Bad_- cap [-10^7, 10^7] = {1, 5, 17}`.
  * The integer periodic points of every period `<= 18`.
  * Incongruent-class counts: all three notions for `k <= 18`; notions (I) and (II) for `k <= 384`.
  * Tree densities of 119 integer roots (4 of them multiples of 3) on both sheets at `X = 10^7`; Z-Collatz on `[-10^6, 10^6]`.
* **CITED:**
  * Krasikov–Lagarias 2003: arXiv abstract and journal data read.
  * Wirsching 1993, 1996, 1997, 1998: bibliographic data via Crossref; content `[R]` via Lagarias's annotated bibliography I, entries 182–187, read. The Springer pages sit behind a bot-check, which was not bypassed.
  * Applegate–Lagarias, Terras, Tao and Gonçalves–Greenfeld–Madrid, as recorded in the [barrier atlas](collatz_procgen_20260922_barrier_atlas.md).
* **OPEN:**
  * Collatz.
  * No divergence for `3n-1`.
  * The growth exponent of the chain-refined counts (III).
  * Whether `(i) AND (ii)` yields to any sheet-symmetric method.
* **REFUTED / CORRECTED:**
  * The reading "a residue mod `3*2^k` determines the predecessor's residue mod `2^k`" is CORRECTED to "mod `2^(k+1)` for `E`, mod `3*2^(k+1)` for `D`, and **not** mod 3 for `E`".
  * A sheet separation by tree growth counted by size `|x|` is REFUTED: `nu` preserves `|x|`. Counted with the sign, the separation is the sign law.

Session `collatz-procgen-20260922` (mac-mini), lane "inverse tree mod 192", 2026-09-24. No HYP or THM file was created.

* **Scripts:** `04-computation/experiments/collatz_procgen_20260924_tree_{automaton,controls,counts}.py`, run by `..._tree_run.sh`.
* **Output:** [`collatz_procgen_20260924_tree.out`](collatz_procgen_20260924_tree.out). It has sections A1–A8, C1–C4 and K1–K7, and every check raises on failure. The run takes about 20 s, peaks at about 400 MB, and uses one process.

## 0. The owner's proposal, read exactly

The owner proposes to build the tree backwards: double every odd number, double the evens, and send each number `2 mod 3` somewhere "mod 192". Three remarks guide the reading:
* this is "logically equivalent to `(3N+1)/2`";
* `A` points to `2A` and `4A` (a "4x fractal recursion") but also to `(4A-1)/3 = A + (A-1)/3`;
* `A/2` is "the branch responsible for the hidden choice between both odds `(N-1)/3` and `(4N-1)/3`".

The owner also asks to study "incongruency itself, singletons or small sets", and "what features represent". Dictionary (all PROVED in §1 and checked in A6):

| owner's words | exact object |
|---|---|
| doubling | `D(x) = 2x`, the even predecessor; always legal |
| "each number 2 mod 3" | `E(x) = (2x-1)/3`, the odd predecessor. It is legal iff `x = 2 mod 3`, and then it is odd. So `T^-1(x) = {2x} U {E(x)}` |
| "equivalent to `(3N+1)/2`" | `T(D(x)) = T(E(x)) = x`: the moves are exactly the two inverse branches of the shortcut map |
| "4x fractal recursion" | `E D^2 = S E` with `S(p) = 4p+1`: doubling twice moves the odd branch one rung up the sibling ladder `p, 4p+1, 16p+5, ...` (Proposition 11) |
| "`(N-1)/3` versus `(4N-1)/3`" | consecutive rungs, `(4N-1)/3 = S((N-1)/3)`, for `N = 1 mod 3` |
| "`A/2` is the branch" | if `N` is even, `(N-1)/3 = E(N/2)` is a genuine odd predecessor of `N/2`. If `N` is odd, `(N-1)/3` is **even**: it is the phantom rung `j = -1`, and the arrow `(N-1)/3 -> N` is exactly the extra arrow of the E-graph relaxation ([choice ladder](collatz_procgen_20260922_choice_ladder.md)). So the owner's "hidden choice" is the sibling index, and E-SCC is the ladder extended one rung down |
| "deciding where it goes mod 192" | the type automaton `A_6` on `Z/192` (§1.3) |

**The coordinator's reading, verified.**
* The branching at `2 mod 3` is TRUE.
* The shared Syracuse image of the siblings `p, 4p+1, 16p+5, ...` is TRUE.
* `192 = 3*64` is TRUE.
* The claim "`n mod 3*2^k` decides the predecessor mod `2^k`" is CORRECTED: it decides `E(n)` mod `2^(k+1)`, one bit more, but not `E(n) mod 3`.
* "Incongruency = the exceptional set" is TRUE, and its dimension `h(log_3 2)` was PROVED in the choice ladder. The mod-3 refinement is in §3.3.
* The foundry's SHEET and DEFECT typing of residue arguments becomes Corollary 7 and Proposition 8 below.

## 1. The exact programme

### 1.1 The equivalence theorem

**Theorem 1.** Let `R` be the smallest subset of `Z_{>0}` that contains `{1,2}` and is closed under `D` and under `E` at legal points. Then
`R = B := {n >= 1 : T^j(n) in {1,2} for some j >= 0}`,
so **Collatz (every positive integer reaches 1 under `T`) holds iff `R = Z_{>0}`**. Moreover `T(R)` is contained in `R`, and the functional graph of `T` on `R` has exactly one cycle, `{1,2}`: its edges are `D(1) = 2` and `E(2) = 1`. Every other node has a finite height `h(n) = min{j : T^j(n) in {1,2}}`, and the owner's generation by moves is breadth-first search by height: `G_(h+1) = T^-1(G_h) \ {1,2}`.

*Proof.*
1. `D` and `E` map positive integers to positive integers: `(2x-1)/3 >= 1/3` for `x >= 1`, and an integer.
2. `T(D(x)) = x` for every `x`. For legal `x`, `E(x)` is odd (since `2x-1` is odd), so `T(E(x)) = (3E(x)+1)/2 = x`.
3. `R` is contained in `B`: `B` contains `{1,2}` and is closed under `D` and `E`, because `T(D(n)) = T(E(n)) = n` lies in `B`.
4. `B` is contained in `R`, by induction on `j`. If `T^j(n)` lies in `{1,2}` with `j >= 1`, then `T(n)` lies in `B` with index `j-1`, hence in `R`.
5. Now `n` lies in `T^-1(T(n)) = {D(T(n))} U {E(T(n))}`, so closure puts `n` in `R`.
6. Reaching 1 is the same as reaching `{1,2}`, since `T(2) = 1`.
7. Every `n` has exactly one `T`-image, so the graph is functional. A second cycle inside `R` would be a cycle disjoint from `{1,2}` all of whose points reach `{1,2}`, which is impossible. ∎

The Syracuse (odd-only) form of this statement is Proposition 4.1 of the [reverse-tree note](collatz_mod6_20260917_reverse_tree_pieces.md) (inherited). Checks (A7): every `n <= 10^6` reaches `{1,2}`. The `D/E`-closure of `{1,2}` inside `[1, 10^5]` contains exactly the `n <= 3000` whose orbit stays `<= 10^5`: 2964 of them. The cap binds at `703, 937, 1055, ...`, because inverse moves do not grow monotonically in size.

### 1.2 The sibling ladder and the exact "fractal recursion"

**Proposition 11 (PROVED; A3, A6).**
1. `E(D^2 x) = (8x-1)/3 = 4(2x-1)/3 + 1 = S(E(x))`, and `x` is legal iff `4x` is.
2. Siblings share the Syracuse image, `U(4p+1) = U(p)`, and `3S^j(p) + 1 = 4^j(3p+1)`. So `v_2(3S^j(p)+1) = v_2(3p+1) + 2j`: the ladder index is read off `v_2(3p+1)`.
3. **Predecessor decomposition.** Let `x` be a positive integer, not a multiple of 3 and not on the root cycle. Let `h0` in `{0,1}` satisfy `2^h0 x = 2 mod 3`, and let `p_0 = E(2^h0 x)` be the minimal odd predecessor. Then
   `Pred*(x) = {2^i x : i >= 0}  disjoint-union  (disjoint union over j >= 0 of Pred*(S^j(p_0)))`,
   where `Pred*` is the full backward tree, node included. For `3 | x` the tree is the doubling chain `{2^i x}` alone. So the owner's "4x recursion" says exactly this: **every backward tree is a doubling chain carrying the backward trees of one sibling ladder**, and the ladder is generated by `S`.
4. **Limits.**
   * `S^j(p) -> -1/3` in `Z_2`, with `-1/3 = E(0)` and `T(-1/3) = 0 = D(0)`.
   * Mod 192, `S` has the single periodic orbit `{21, 85, 149}`: all three residues mod 3 of `-1/3 mod 64`. Every ladder enters it within 3 rungs.
   * The trunk `(4^i-1)/3 = E D^(2i-1)(1) = S^(i-1)(1)` reads `1, 5, 21, 85, 149, 21, ...` mod 192.
   * Mod 9 the rows of a ladder rotate with period `3 = ord_9(4)` (inherited, reverse-tree note Theorem 1.1).

*Proof of 3.* Unfold `Pred*(x) = {x} U Pred*(2x) U Pred*(E(x))` along the doubling chain. `E(2^i x)` is legal iff `i = h0 mod 2`, and `E(2^(h0+2j) x) = S^j(p_0)` by item 1. The pieces are disjoint because `T` is a function and `x` is not on a cycle. ∎

### 1.3 The type automaton at modulus `3*2^k`

**Proposition 2 (PROVED; A1, A4).** For types `n mod 3*2^k`:
* **Forward.** `T(n) mod 3*2^(k-1)` is a function of `n mod 3*2^k`. More precisely, one `T`-step consumes one bit and each odd step *produces* one trit: `n mod 3^m 2^k` decides `T^j(n) mod 3^(m+a_j) 2^(k-j)`, where `a_j` counts the odd steps.
* **Backward.**
  * `D` sends `n mod 3*2^k` to `2n mod 3*2^(k+1)`: it gains a bit and keeps the trit.
  * For `n = 2 mod 3`, `E` sends `n mod 3*2^k` to `E(n) mod 2^(k+1)`, a bijection onto the odd classes. It gains a bit but consumes the trit: `E(n) mod 3*2^j` needs `n mod 9*2^(j-1)`, and the three lifts of a class mod 192 give `E(n) = 0, 1, 2 mod 3` once each.

*Proof.* Write `T^j(n) = (3^(a_j) n + c_j)/2^j`. The first `j` parity bits, hence `c_j`, depend on `n mod 2^j` (Terras). Then `3^(a_j) n mod 3^(m+a_j)` needs only `n mod 3^m`, and `2^j` is a unit mod `3^(m+a_j)`. Backward, `2n - 1 mod 3*2^(k+1)` is determined by `n mod 3*2^k`, and division by 3 leaves `2^(k+1)` bits. The lifts `n + 3*2^k t`, `t = 0, 1, 2`, change `E(n)` by `2^(k+1) t`, which runs through the residues mod 3. ∎

**The automaton `A_k`.** Its states are the elements of `Z/3*2^k`. Every state has the edge `r -D-> 2r`. A state `r = 2 mod 3` also has three edges `r -E-> s`, one for each `s` with `s = E(r) mod 2^k`, that is one per residue of `s mod 3`, each of Haar weight 1/3. The automaton is a *sound over-approximation*: if `m` lies in `T^-1(n)`, then `type(m)` lies in `delta(type(n))`. The nondeterminism sits exactly in the trit that `E` consumes.

**Proposition 3 (PROVED for every `k`; FINITE-EXACT check `k <= 10`, A2).**
1. The strongly connected components of `A_k` are:
   * the core `C_k`, made of the `2^(k+1)` non-multiples of 3;
   * `2^k` singletons, the multiples of 3, with a loop only at `0`.
2. The multiples of 3 form a closed region with `D`-moves only, and `D` is nilpotent on its 2-part: it drains every multiple of 3 to `0` within `k` steps. Every type is reachable from the root types `{1, 2}`.
3. The weighted matrix (weight 1 on a `D`-edge, `1/3` on each `E`-edge) has spectral radius `4/3`. On the core its Perron vector is `1` on `r = 1 mod 3` and `4/3` on `r = 2 mod 3`, from `3L^2 - L - 4 = 0` (roots `4/3` and `-1`).
4. The forward automaton, `r -> T(r)` with the top bit free, is the mirror image. The same core is one component; the multiples of 3 are singletons, and forward the core never re-enters `3Z`. So `3Z` is a **sink backward and a source forward**: these are the leaf cones.

*Proof.*
* *Reachability of every type from the core.* From a core state, apply `D` at most once to make `E` legal. After an `E`, choose the target trit so that the next `E` in the desired word is legal; `D` swaps `1 <-> 2` and fixes `0`.
* *Hitting a prescribed `s`.* The last `k` moves are the Terras word of `s mod 2^k`, and the trit chosen at the last `E` gives `s mod 3`.
* *Structure of `3Z`.* From `3Z` only `D` applies, and `2^k r = 0` on the 2-part.
* *Perron vector.* Put `v = (1, 4/3)`. For `r = 1`, `(Wv)(r) = v(2r) = 4/3`. For `r = 2`, `(Wv)(r) = 1 + (1 + 4/3)/3 = 16/9`. A positive eigenvector of an irreducible block gives its Perron root, and the `3Z` block has spectral radius 1. ∎

**At modulus 192** (A2, A3):
* 64 branch states (`r = 2 mod 3`) have out-degree 4, and the other 128 states have out-degree 1.
* The mean weighted branching is `4/3`: `1` from `r = 1 mod 3`, `1 + 2/3` non-leaf from `r = 2 mod 3`. An `E`-child is a leaf with probability 1/3.
* **Finite shadows of the generators:**
  * `D` has the fixed type `0` and the 2-cycle `{64, 128}`;
  * `E` has a self-loop at `191 = -1`, because `E(-1) = -1`;
  * `S` has the 3-cycle `{21, 85, 149}`;
  * the root is the 2-cycle `1 -D-> 2 -E-> 1`.

### 1.4 What a residue mod 192 decides

**Proposition 4 (PROVED, Proposition 2 with `k = 6`; checked on 59 lifts of every class on both sheets, A4).** A residue `n mod 192` decides:
* the parity window `(x_0, ..., x_5)`, since `Z/64 -> {0,1}^6` is a bijection;
* the mod-3 types of `T^0(n), ..., T^6(n)`: a type is `2` right after an odd step and alternates `1, 2` along halvings;
* `T^6(n) mod 3^(1+a)` (and not mod `3^(2+a)`);
* backward, `D(n) mod 384` and `E(n) mod 128`.

The classes by the number `a` of odd steps in the window number `3 C(6,a)`.

On a tree node `m = w(root)`, `m mod 64` is the **last six moves** of its path (`E` = odd, `D` = even) and `m mod 3` is its own branching type. So **the owner's type is one branching trit plus a six-move window.** The two halves concern different primes: the window is 2-adic and decides the forward parities, while the trit is 3-adic and decides one backward branching. Counting predecessors is a purely 3-adic question, because every word is 2-adically realizable.

**Proposition 5 (PROVED by a finite enumeration of all words of length `<= 6`; A5).** In a window `x_i -> x_j` with `m = j - i <= 6` and `a` odd steps, `2^m x_j = 3^a x_i + b c`, where `c >= 0`. A comparison can contradict the word prediction ("up iff `3^a > 2^m`") only if `x_i` is below the window's gate. The largest gate is `76/5` on the plus sheet (word `00111`) and `260/17` on the minus sheet (word `001111`). Every value in the window is `>= n/64`, so **for `n > 979` the full order pattern of `T^0(n), ..., T^6(n)` is the word's generic pattern**: a function of `n mod 64`, identical on the two sheets after `r -> -r`. The census up to `10^6` shows that the true exceptions are `n in {1,2,3,4,5,8,10,16,32}` on the plus sheet and 21 starts `<= 80` on the minus sheet. This is the owner's window seen through the [order-laws note](collatz_procgen_20260922_order_laws.md)'s short-lag genericity. There every window of at most 12 odd steps is sheet-blind, and the first sheet-separating window, `165 -> 163`, has 12 odd steps and 19 halvings, far beyond the six steps a residue mod 192 sees.

## 2. The no-go

### 2.1 Transport

**Theorem 6 (PROVED).** Let `X` be `Z`, or the rationals with denominator prime to 6, or `Z_2 x Z_3`. The maps `T_b`, `D` and `E_b(x) = (2x-b)/3` are defined on all of these. Let `nu(x) = -x`. Then:
* (a) `nu T_+ = T_- nu`, `nu D = D nu`, and `x` is `E_+`-legal iff `-x` is `E_-`-legal (`2 mod 3` corresponds to `1 mod 3`), with `nu E_+ = E_- nu`;
* (b) `nu` is an automorphism of the additive group. It sends each class `r + MZ` to `-r + MZ`, preserves Haar measure on `Z_2 x Z_3`, maps `Z` onto `Z`, and preserves `|x|_l` for every prime `l` and the real size `|x|_inf`;
* (c) `nu` reverses the order, and `nu(Z_{>0}) = Z_{<0}`.

Consequently the structures `M_b = (X; +, residue predicates, T_b, D, E_b, Z, Haar, (|.|_v)_v)` for `b = +1` and `b = -1` are isomorphic via `nu`. **Every statement about `M_b`, first- or higher-order, has the same truth value on the two sheets.** In particular `r -> -r` is an isomorphism `A_k^+ = A_k^-` of the owner's automata for every `k` (FINITE-EXACT check `k <= 12`, A8).

*Proof.*
* *Even `x`:* `T_-(-x) = -x/2 = -T_+(x)`.
* *Odd `x`:* `T_-(-x) = (-3x-1)/2 = -T_+(x)`, since `x` and `-x` have the same parity.
* *The odd branch:* `E_-(-x) = (-2x+1)/3 = -E_+(x)`, and `3 | 2x-1` iff `3 | -2x+1`.
* *Items (b) and (c):* `|-1|_v = 1` at every place, and `nu` is the Haar-preserving group automorphism `x -> -x`. ∎

### 2.2 The SHEET no-go

**Corollary 7 (PROVED).** Write `Coll(T)` for "every positive integer's `T`-orbit reaches the cycle through 1". Suppose a proof derives `Coll(T)` from hypotheses `H_1(T), ..., H_m(T)` about the map `T`, used as a parameter, together with mathematics that does not mention `T`. Suppose further that each `H_i` is a statement about `M_T`. It may use `|x|` and points definable from `T`, but not the fixed cone `Z_{>0}` or the fixed point `1`. Then all `H_i` hold for `T_-`, by Theorem 6 and their truth for `T_+`. The same proof would then give `Coll(T_-)`, but `T_-` has the cycle `5 -> 7 -> 10 -> 5`. **Hence every proof of Collatz uses a hypothesis `H` with `H(T_+)` true and `H(T_-)` false. By Theorem 6(c), such an `H` refers to the fixed positive cone, that is, to the sign of points; the remark below makes this exact as side-awareness.**

Every statement of the type language transports, in any combination: residues mod `3*2^k` for one `k` or for all `k` at once (the whole of `Z_2 x Z_3`), the automata `A_k` and their spectra, Haar densities, the integer predicate, the `l`-adic structure at any set of primes, and the size `|x|` (counting by size, logarithmic density, drift). The conclusion `Coll(T)` names the fixed positive cone `Z_{>0}` and the point `1`, which is exactly why no `nu`-invariant premises can imply it.

**What exactly is missing: the side.** Corollary 7 does not say that Collatz is inexpressible in the type-plus-size language. A `T`-definable half-line exists: the side of the odd point `p` of the unique cycle with word `(1,0)`, namely `{x : |x + p| = |x| + |p|}`. Relative to it, `H(T)` := "every nonzero integer on the side of `p` reaches `p`'s cycle" is a statement about `M_T`. `H(T_+)` is Collatz, and `H(T_-)` is its transport (the negative half-line of 3n-1). So the corollary's real content is this:
* A proof must be **side-aware**: it must use a property that holds on the side of the contracting 2-cycle and fails on the opposite side. On the opposite side `T_+` has the cycles `-1`, `-5`, `-17`, the mirror images of the positive 3n-1 cycles.
* Residue types, Haar densities, integrality and `l`-adic data are **side-blind**: every residue class contains positive and negative integers in the same proportion. Whatever they prove about the positive half they prove about the negative half, where it is false. Formally this is Corollary 7 again: `nu` carries the negative half of 3n+1 onto the positive half of 3n-1. This is the SHEET barrier in exact form, and the owner's automaton is side-blind.
* The size `|x|` is `nu`-invariant too. It becomes side-aware only together with a reference point, and then the side-aware input is the sign law of §2.4 or an equivalent statement.
* The last link, "the side of `p` is `Z_{>0}`", is the trivial fact `1 > 0`.

The control data (C1–C2, FINITE-EXACT):

| sheet | positive cycles | multiplier `3^L/2^K` | basin densities at `10^7` | `Bad cap [-10^7, 10^7]` |
|---|---|---|---|---|
| `3n+1` | `{1,2}` | `3/4` (contracting) | `1` | `{-17, -5, -1}` |
| `3n-1` | `{1}`, `{5,7,10}`, `{17, ..., 136}` | `3/2, 9/8, 2187/2048` (expanding) | `0.32738, 0.32450, 0.34812` | `{1, 5, 17}` |

**The owner's type cannot see that there are three basins.** Each 3n-1 basin is residue-equidistributed to sampling accuracy. The chi-square statistics against a uniform profile are:
* mod 3 (2 d.o.f.): `0.3, 0.5, 0.5`;
* mod 9 (8 d.o.f.): `2.2, 4.0, 5.0`;
* mod 8 (7 d.o.f.): `6.8, 16.0, 4.7`;
* mod 192 (191 d.o.f.): `197, 218, 195`.

Every class mod 192 (and mod 8 and mod 9) meets each basin in that basin's global proportion, to sampling accuracy. Seen by types, a basin of density about 1/3 looks like a uniformly random third.

The sign asymmetry already sits in the owner's own move. `E_b(x) + b = (2/3)(x + b)`, so `E_b` contracts toward its fixed point `-b`:
* on the plus sheet the fixed point `-1` lies outside `Z_{>0}`, so for every positive `n` in a legal class the multiplier `2/3` is actual descent (`E(n) < n` iff `n > -1`);
* on the minus sheet the fixed point `+1` **is** the positive root cycle, and the certificate `E_-(n) < n` fails exactly at `n = 1`.

The type `191 = -1 mod 192` is decided by the owner's refinement; the hostile point `-1` is not (§3.3, §3.4).

### 2.3 The DEFECT controls

**Proposition 8 (PROVED; C3).**
1. **`T1(n) = 2n` if `n = 3*2^m` (`m >= 1`), and `T1(n) = T(n)` otherwise.**
   * A multiple of 3 is never `E`-legal, so `T^-1(3*2^m) = {3*2^(m+1)}` and `T1^-1(S1)` is contained in `S1 = {3*2^m : m >= 1}`.
   * Hence every orbit that starts outside `S1` is its `T`-orbit, while `S1` itself is a divergent `T1`-orbit.
   * `S1` has `floor(log2(X/3))` elements up to `X`.
   * `T1 - T = 9*2^(m-1)` at `3*2^m`, so at level `k` the type map of `T1` equals that of `T` except at `k-1` integers.
2. **`T2`, with agreement at every 2- and 3-adic level.**
   * Construction: `s_1 = 12`, then `s_(i+1) = s_i/2 + 6^i k_i` with the least `k_i` such that `2^(i+2) | s_(i+1)` and `s_(i+1) > s_i^2`. This gives `12, 168, 28272, 799306368, ...`.
   * Set `T2(s_i) = s_(i+1)` and `T2 = T` elsewhere.
   * Then `T2(s_i) = T(s_i) mod 6^i`, so the discrepancy tends to 0 in `Z_2` **and** in `Z_3`.
   * Each `s_i` is a multiple of 3, and the affected starting points are the disjoint doubling chains `{s_i 2^j}`: `O(log X log log X)` of them up to `X`.

Every statement about `T` that tolerates exceptions on a density-zero set of starting points holds verbatim for `T1` and `T2`. This covers densities, logarithmic densities (Terras, Korec, Tao, GGM), "almost all" statements, lower bounds `pi_a(x) >= x^gamma` for `gamma < 1`, and residue-class statements with density-zero exceptions. Yet `Coll(T1)` and `Coll(T2)` are false. **So a proof must also use pointwise (every-orbit) arithmetic.** *Proof:* the listed invariances, checked in C3. ∎

### 2.4 What any proof must use, precisely

1. **The sign (order).** By Corollary 7 this is needed. Size `|x|` is not enough.
   * The canonical sign-sensitive fact is the **sign law**: for every `x` and every prefix word `w` of its orbit (length `p`, `a` odd steps), `2^p T^p(x) - 3^a x = b c_w`, where `c_w >= 0` depends only on `w` and `c_w > 0` when `a >= 1` (PROVED, one line). Read against `x > 0` it says `T_+^p(x) >= 3^a x/2^p` and `T_-^p(x) <= 3^a x/2^p`: a hypothesis true for `T_+` and false for `T_-`, as Corollary 7 demands.
   * Equivalently: a positive `T_b`-cycle is contracting (`2^K > 3^L`) on the plus sheet and expanding on the minus sheet. The five integer cycles of `T_+` confirm this (C4): `0` and `{1,2}` are contracting and nonnegative, while `{-1}`, `{-5,-7,-10}` and `{-17,...}` are expanding and negative.
2. **The cycle gates.** A primitive word `w` has the unique rational periodic point `x_w = c_w/(2^p - 3^a)`.
   * Its **sign** is `sign(2^p - 3^a)`: archimedean.
   * Its **integrality** is `(2^p - 3^a) | c_w`: a congruence **modulo a number prime to 6**, invisible to every type mod `3*2^k`.
   * Each class-level certificate becomes actual descent only above the gate `c_w/(2^p - 3^a)`. For all periods `p <= 18` the integer periodic points are exactly the five known cycles (FINITE-EXACT, K5; classical for small periods).
3. **The density-zero exceptional set.**
   * `Bad_+`, the parity sequences with every prefix above the line, has dimension `h(log_3 2) = 0.9500` (PROVED, choice ladder), and a positive integer in `Bad_+` has a divergent orbit (choice ladder §5b).
   * So no-divergence **requires** `Bad_+ cap Z_{>0}` to be empty. That is a statement about integers in a null set: DEFECT forbids deriving it from densities, and SHEET forbids deriving it from order-free data.
   * On `[-10^7, 10^7]`, `Bad_+ cap Z = {-17, -5, -1}`: exactly the points closest to 0 on the three expanding integer cycles. Their mirror is `{1, 5, 17}`, the minima of the 3n-1 cycles. **The order-free content is the same on both sheets; only the signs of three integers differ.**

## 3. What the automaton can prove

### 3.1 Exact Haar identities for predecessor counts

**Proposition 9 (PROVED; K1).** Let `N_d(a) = |T^-d(a)|`, where legality is 3-adic only. It is a function on `Z/3^d`, given by `N_d(a) = N_(d-1)(2a) + [a = 2 mod 3] N_(d-1)((2a-1)/3)`. For Haar-random `a` in `Z_3`:
* `E[N_(d,e)] = C(d,e) 3^-e` (with `e` the number of `E`-steps), `E[N_d] = (4/3)^d`, and `E[sum N_(d,e) z^d u^e] = 1/(1 - z(1 + u/3))`;
* per class, `m0 = 1`, `m1_d = m2_(d-1)` and `m2_d = m1_(d-1) + (4/3)^(d-1)`.

*Proof.* If `x` is Haar-uniform on `{x = 2 mod 3}`, then `E(x)` is Haar-uniform on `Z_3`, since `E` scales measure by 3 onto `Z_3`. Also `D` permutes classes. So every `E`-step is legal with conditional probability `1/3`, independently. ∎

FINITE-EXACT, `d <= 14`:

| `d` | `E[N_d]` | min over units | at class | `E[N^2]/E[N]^2` |
|---|---|---|---|---|
| 4 | 3.16 | 3 | 5 | 1.070 |
| 8 | 9.99 | 8 | 25 | 1.078 |
| 12 | 31.57 | 25 | 25 | 1.088 |
| 14 | 56.12 | 43 | 25 | 1.089 |

Depth counts are tightly concentrated: at `d = 14` the minimum is half of the unit mean `83.7`. The thinnest trees belong to the classes of small integers (5, 7, 25, 34).

**Wirsching (CITED `[R]`, Lagarias bibliography I #182–187; bibliographic data via Crossref).**
* Wirsching bounds the number `|P_T^n(a)|` of predecessors of `a` in `[2^n a, 2^(n+1) a)` from below by `s_n(a) = sum_l e_l(n + floor(l log2(3/2)), a)`, built from counting functions `e_l(k, a)` (LNM 1681, Theorem II.4.9). These functions extend to 3-adic `a`.
* `s_n = +inf` on a dense set of 3-adic `a`; the 3-adic mean satisfies `liminf s̄_n/2^n > 0` (Theorem III.5.2).
* His *Heuristic Principle* is `s_n(a) > c_1 s̄_n`.
* His 1997 3-adic equidistribution hypothesis would give `pi_a(x) >= x^(1-eps)`; his 1993 paper proves `x^0.48`.

**Krasikov–Lagarias (CITED `[P]`, arXiv math/0205002 abstract; Acta Arith. 109 (2003) 237–258).** Their difference inequalities on residue classes mod `3^k`, solved by computer-aided nonlinear programming, give `pi_1(x) >= x^0.84` (and, per the atlas, `pi_a(x) >= x^0.84` for every `a` prime to 3 and `x >= x_0(a)`). Both lines are exactly "type automaton at modulus `3^k` + size counting". By Corollary 7 both are sheet-blind, and the atlas records K–L as sheet-blind.

### 3.2 The averaged size series: the drift is a residue

**Proposition 10 (PROVED; K2).** Ignore carries, so that a depth-`d` predecessor with `e` `E`-steps has relative size `2^d/3^e`. The Haar-averaged series is

`F(s) = sum C(d,e) 3^-e (2^d/3^e)^-s = 1/(1 - g(s)),   g(s) = 2^-s + (3/2)^s/3`.

* `g(1) = 1/2 + 1/2 = 1`: the arithmetic mean of the step factors `1/2, 3/2` is 1.
* `g(2) = 1/4 + 3/4 = 1`, and `g < 1` exactly on `(1, 2)`.
* `-g'(1) = (1/2) log(4/3) = log(2/sqrt 3) = 0.143841`, **the Collatz drift**, so the residue at `s = 1` is `R = 1/log(2/sqrt3) = 6.95212`.
* The pole at `s = 1` is the only one on `Re s = 1`: `|g(1+it)| = 1` would force `2^-it = (3/2)^it = 1`.
* By Wiener–Ikehara (non-lattice, since `log_2 3` is irrational) the averaged number of predecessors of relative size `<= X` is `~ R X`.
* Per class mod 3 the residues are `0, R, 2R`: a root `= 2 mod 3` has twice the tree of a root `= 1 mod 3`, because orbits visit `2 mod 3` twice as often (stationary law `2/3, 1/3`).

Numerically, `A(2^n)/2^n = 6.9724, 6.9688, 6.9556, 6.9431` at `n = 20, 40, 60, 80`. These oscillate slowly around `R`, the near-lattice oscillation coming from the convergents `19/12`, `84/53`, ... of `log2 3`.

**Reading (the owner's AM–GM principle in the tree).**
* AM = 1 is the pole at `s = 1`: the averaged tree has positive density.
* The AM–GM gap `log(2/sqrt3)` is the reciprocal of the residue. It is the renewal constant: `1/(y log(2/sqrt3))` visits per unit length at size `y`.
* For `qx+1`, `g_q(s) = 2^-s (1 + q^(s-1))` has roots `{1, 2}` for `q = 3`, but `{0.6509, 1}`, `{0.3735, 1}` and `{0.2581, 1}` for `q = 5, 7, 9`. So **the automaton sees DRIFT exactly** (positive versus zero density), matching the K–L prediction `eta_5 ~ 0.650` quoted in the atlas.
* `g` does not contain the offset `b`, so the automaton is SHEET-blind.

**Against actual trees (FINITE-EXACT, `X = 10^7`, K3).** Let `rho = a * density / (c R)`.

| roots | mean `rho`, plus | mean `rho`, minus | spread |
|---|---|---|---|
| `4..20` | 0.71 | 0.23 | archimedean boundary |
| 64 consecutive units near `10^4` | `1.02 +- 0.17` | `0.91 +- 0.15` | `rho` in `0.06..6.8` |

* For large roots the model is right on average **on both sheets**. Individual trees scatter by two orders of magnitude with their 3-adic expansion; Wirsching's `s_n(a)` is unbounded on a dense set of 3-adic `a`.
* The density is governed by predecessors with about half `E`-steps (the `s = 1` tilt), a large-deviation regime of the Haar statistics (HEURISTIC reading).
* Near the root the model fails by order-one factors. On the plus sheet `tree(4)` is everything `>= 3`, so `rho = 4/R = 0.575`.
* Small 3n-1 roots are capped by their basin (`0.327`, `0.3245`, `0.348`). **SHEET is visible only at the small end, where carries are as large as the values.**

### 3.3 Incongruent classes mod `3*2^k`, with the mod-3 refinement

A class `(n mod 3, parity word of length k)` is *decided* if the class determines a certificate with multiplier `< 1`.
* **(I) Forward only (Terras):** some prefix has `3^(a_j) < 2^j`. The mod-3 digit plays no role, so the undecided count is exactly `3 |Bad_k|` (PROVED: the word depends on `n mod 2^k` alone).
* **(II) Owner-refined:** (I); or `n = 2 mod 3`, with the certificate `E(n) = (2n-1)/3 < n`; or an even step `j-1 -> j` with `T^j(n) = 2 mod 3` and `(2/3) 3^(a_j)/2^j < 1`. The last is the "A/2 branch" `E(T^j n) = (T^(j-1) n - 1)/3`.
* **(III) Chain-refined:** (II), plus at every window point the non-retracing predecessor branch, followed along its minimal-child chain for as many `E`-steps as the class decides. Recall that `T^j(n)` is known mod `3^(1+a_j)`: odd steps manufacture 3-adic precision.

| `k` | `3*2^k` | (I) | (II) | (III) | (II)/(I) | (III)/(I) | `log2(I)/k` | `log2(II)/k` |
|---|---|---|---|---|---|---|---|---|
| 6 | 192 | 24 | 14 | 14 | 0.583 | 0.583 | 0.764 | 0.635 |
| 12 | 12288 | 678 | 384 | 288 | 0.566 | 0.425 | 0.784 | 0.715 |
| 18 | 786432 | 22485 | 12684 | 9200 | 0.564 | 0.409 | 0.803 | 0.757 |
| 64 | | `8.26e16` | `4.44e16` | | 0.537 | | 0.878 | 0.864 |
| 384 | | `2.48e107` | `1.31e107` | | 0.527 | | 0.929 | 0.927 |

Checks: the DP agrees with brute force for `k <= 12`, and (III) is a class invariant on three lifts per class for `k <= 10`. At the owner's modulus the undecided classes are:
* under (I), 24 classes: the three lifts of `{7, 15, 27, 31, 39, 47, 59, 63}` mod 64;
* under (II), 14 classes: `{7, 27, 31, 39, 63, 91, 103, 111, 123, 127, 135, 159, 175, 187}`.

All of them begin with two odd steps and are `0` or `1 mod 3`.

**Proposition 12 (PROVED).** The (II)-exceptional set has 2-adic dimension `h(log_3 2) = 0.9500`, and `log2(#(II)-undecided mod 3*2^k)/k -> h(log_3 2)`.

*Proof.*
* The (II)-exceptional set is contained in `Bad`, which gives the upper bound.
* Lower bound: let `W` be the set of parity sequences with `M_j = 3^(a_j)/2^j >= 3/2` for every `j >= 1`. On `W` no (II) certificate fires, since branch certificates need `M_j < 3/2`.
* Take the prefix `1^N`, then Bernoulli(`p'`) with `p' > log_3 2`. The walk `log M_j` has positive drift and stays above `log(3/2)` with positive probability. So `W` has positive `mu_(p')`-measure, and `dim W >= h(p')` (Besicovitch–Eggleston). Let `p' -> log_3 2`.
* The count statement follows from the same Chernoff and ballot bounds as for `Bad_k`. ∎

**So "incongruency counted with the mod-3 refinement" is still `2^(0.95k)`:** the owner's trit removes about half of the undecided classes and never the exponent. With unbounded 3-adic lookahead the certificates become the undirected game, which is equivalent to Collatz and does not collapse the exceptional set ([choice ladder](collatz_procgen_20260922_choice_ladder.md) §6). The exponent of (III) is FINITE-EXACT only (OPEN).

### 3.4 Singletons and small sets

| point or set | where it lives | in the owner's automaton | what it represents |
|---|---|---|---|
| `0` | fixed point of `D`; word `(0)`, contracting, gate `0` | type `0`, the `D`-loop; sink of the leaf cones | the doubling chain of a multiple of 3 |
| `-1` | fixed point of `E_+` and of `T_+`; word `(1)` | type `191`, the `E`-loop; the class is decided under (II), the point never descends | Applegate–Lagarias's single class; the extreme point of `Bad_+`; the mirror of the 3n-1 root `1` |
| `-1/3` | `E(0)`; 2-adic limit of every sibling ladder; `T(-1/3) = 0` | the `S`-cycle `{21, 85, 149}` | the "4x recursion" pushed to its limit |
| `1/2` | `E(1/2) = 0` in `Z_3`; hostile point of Q2 | reached by the trunk only 3-adically: the trunk isometry `i -> (4^i-1)/3` takes the value `1/2` at an irrational 3-adic `i*` ([trunk note](collatz_procgen_20260923_trunk_rh.md)) | the exits of Q2's hard neighbourhood are the trunk |
| `{1, 2}` | the root cycle | the 2-cycle `1 -D-> 2 -E-> 1` | contracting (`3/4`): the sign law's positive side |
| `{-1, -5, -17}` | `Bad_+ cap Z` on `[-10^7, 10^7]` | classes of `Bad_k` for every `k` | expanding integer cycles; the mirror `{1, 5, 17}` is the set of 3n-1 cycle minima |
| rational points of `Bad` | `1, 0, 1, 2, 3, 6, 12, 16, 36, 60, 127, 216, 366, 721, 1290, 2095, 4227, 7451` for periods `p = 1..18` (K5) | the rational skeleton of the incongruent set, about `2^(0.95p)/p` per period | all negative (sign law); mirror: positive rational 3n-1 cycles |
| the 16-point trap `F` of `E_{6 mod 8}` | `-F = {1, 2, 4, 5, 7, 8, 10, 14, 16, 20, 29, 32, 43, 64, 86, 128}` | (Proposition 13) | the 3n-1 cycles seen from inside the plus sheet |

**Proposition 13 (PROVED, finite check; K6).**
* `-F` is closed under the minus-sheet relaxation: `x` odd goes to `3x-1`; `x` even goes to `x/2`, and also to `3x-1` when `x = 2 mod 8`, the mirror of `6 mod 8`.
* `-F` is the orbit of `1` under this relaxation, and it contains the 3n-1 cycles `1 -> 2 -> 1` and `5 -> 14 -> 7 -> 20 -> 10 -> 5` (C-form), glued by the excursions at `2` and `10`.
* The trap's rate `2/3` is the rate of the 3n-1 cycle `{5, 7, 10}` (2 odd steps, 3 halvings, multiplier `9/8`).
* In the full relaxation `E`, `-1 -> -2 -> -5 -> -14 -> -41` uses `3x+1` at the even `-14`, and `-41` lies on the C-form `-17` cycle.

**So the trap that refuted HYP-9121 ([hypothesis sweep](collatz_procgen_20260923_hypothesis_sweep.md) §2) is SHEET made visible inside the plus sheet: the hostile neighbourhood of `-1` is the mirror image of all three 3n-1 cycles.**

### 3.5 The renewal identity, and why types cannot fix a density

**Proposition 14 (PROVED; K7 checks it exactly for the plus basin and each 3n-1 basin at four values of `X`).** Let `R` satisfy `T^-1(R) = R` (a basin or a union of basins) and let `G` be the guard class. Then `R = D(R)` disjoint-union `E(R cap G)`, and for every `X`
`N_R(X) = N_R(floor(X/2)) + #{n in R cap G : E(n) <= X}`.

If `R` has density `delta` and is equidistributed mod 3, as all four basins are (§2.2), the main terms read `delta = delta/2 + (3/2)(delta/3)`, which holds **for every** `delta`. The identity is neutral because the arithmetic mean of the step factors is 1. So the automaton's renewal structure cannot tell density 1 (plus) from `0.327, 0.3245, 0.348` (minus). Those values are boundary data from the small numbers, i.e. non-residue information.

### 3.6 What each feature represents

| feature | exact meaning | sees | blind to |
|---|---|---|---|
| branching at `2 mod 3` | 3-adic legality of the odd inverse branch; the image class of odd steps. Probability 1/3, mean branching 4/3 | DRIFT (via `4/3` and `g(s)`) | SHEET (becomes `1 mod 3`), 2-adic data |
| the `4p+1` ladder | `E D^2 = S E`; siblings share `U`; the ladder index is `(v_2(3p+1) - v_2(3p_0+1))/2`; 2-adic limit `-1/3 = E(0)`; exact decomposition of `Pred*` | the 3-adic row rotation (period 3) | order |
| the exceptional (incongruent) set | parity sequences above the line `a/j = log_3 2`; dimension `h(log_3 2)`; skeleton of negative rational cycles | where certificates fail | whether a positive integer lies in it |
| the drift | AM–GM gap `log(2/sqrt3)` `= -g'(1)`, the reciprocal residue of the tree series | DRIFT exactly | SHEET, DEFECT |
| `sign(2^K - 3^L)` | sign of the rational periodic point of a word (contracting iff positive on the plus sheet) | the **only** sheet-separating datum (Corollary 7) | integrality |
| the gate `c_w/(2^p-3^a)` | fixed point of the word's affine map; threshold above which a certificate descends; integrality is a congruence mod `2^p - 3^a` | cycles, crossings | residues mod `3*2^k` |
| the trunk `(4^i-1)/3` | the sibling ladder of the root, `S^(i-1)(1) = E D^(2i-1)(1)`; a 3-adic isometry (trunk note) | Q2 exits | SHEET, DRIFT |

## 4. A reformulation in which the missing input is clean

**Three layers.**
1. *Types* (residues, Haar measure). These prove measure and dimension statements and the exact averaged identities of §3. They are blind to SHEET, INTEGRAL and DEFECT.
2. *Types + integrality* (the predicate `Z`). This layer can *state* statements about every integer orbit, but it is still `nu`-invariant (Theorem 6).
3. *+ sign.* Only this layer contains Collatz.

**Theorem 15 (PROVED).** Consider the two statements about a map `T`:
* **(i)** no integer has an unbounded `T`-orbit;
* **(ii)** every `T`-periodic integer whose period word is contracting (`2^p > 3^a`) is `0` or lies on the cycle with word `(1,0)`.

Both are statements of layer 2, so `(i) AND (ii)` has the same truth value for `T_+` and `T_-`. With the proved sign law,
`(i) AND (ii)  <=>  Collatz  AND  (no 3n-1 orbit of a positive integer diverges)`.

*Proof.*
* (`<=`) Collatz bounds the positive orbits. No-divergence for 3n-1, transported by `nu`, bounds the negative orbits, and `0` is fixed; this gives (i). For (ii), a contracting periodic integer `x` satisfies `x(2^p - 3^a) = c_w`, which is `>= 0` and `> 0` unless the word is `(0)` (then `x = 0`). So `x > 0`, and Collatz forces `x` into `{1,2}`, whose word is `(1,0)`.
* (`=>`) By (i) every positive orbit is eventually periodic. It stays positive, so it enters a positive cycle. By the sign law that cycle is contracting, and by (ii) it is `{1,2}`. Negative orbits are bounded by (i), which is the 3n-1 statement after `nu`. ∎

FINITE-EXACT (C4): `(i) AND (ii)`, and even the stronger five-cycle Z-Collatz, hold on `[-10^6, 10^6]`. The words are `0`, `10` (contracting) and `1`, `110`, `11110111000` (expanding).

**What this says.**
* Side-blind methods cannot distinguish the two half-lines of `T_+`, so they can only aim at **two-sided** targets. The natural one is `(i) AND (ii)`. It **contains the 3n-1 no-divergence conjecture** and must allow the three expanding cycles `{-1}`, `{-5,-7,-10}` and `{-17,...}` as exceptions.
* Once such a target is reached, the whole order input is the proved sign law; for the stronger five-cycle Z-Collatz it is a finite sign check of the named cycles. This is the "clean, checkable form" asked for.
* Side-aware methods can aim at `H(T)` of §2.2 directly, with no 3n-1 overhead. The side-aware input they need is the sign law or an equivalent statement. Theorem 15 does not claim that side-blind methods must prove the 3n-1 statement; it gives the natural two-sided target and its exact cost.
* What remains in `(i) AND (ii)`, and in `H`, is not a sheet problem. It is INTEGRAL (which rational periodic points are integers: congruences mod `2^p - 3^a`) and DEFECT (pointwise avoidance of the null set `Bad`).

**Tree growth counted by size.**
* **Counted by `|x|` it cannot separate the sheets (PROVED, Theorem 6(b)):** the `T_+`-tree of `{1,2}` and the `T_-`-tree of `{-1,-2}` have identical size-counted growth.
* The sheets differ only in *which side* carries which tree. On the plus sheet the positive side is one tree (density 1). On the negative side there are three trees, of densities `0.3274, 0.3245, 0.3481` by `nu` from C1, each residue-equidistributed.
* The sheet-symmetric statement "the tree of the unique non-zero contracting integer cycle has density 1 on its own side" is true for both sheets (conjecturally). With the sign law it gives the density-one form of Collatz, which is still DEFECT-blind (Proposition 8).
* No size-counted statement separates the sheets without the sign, and with the sign it is the sign law.

## 5. Assessment

* **What the reformulation does.**
  * It makes the owner's programme exact: an equivalence theorem, a finite automaton at every modulus `3*2^k` with an explicit precision flow, and a complete description of `A_6`.
  * It turns the foundry's SHEET and DEFECT typing into theorems. The sign, and only the sign, separates the sheets; size, integrality and every `l`-adic datum are sheet-blind.
  * It proves exact identities: Haar predecessor counts; the averaged tree density `cR/a`, whose constant is the reciprocal drift; and the renewal identity.
  * It identifies the refuted HYP-9121 trap and the integer points of `Bad` as the mirror of the three 3n-1 cycles.
* **What it cannot do.**
  * It proves nothing new about Collatz.
  * Every exact identity above is an average over 3-adic classes. Individual trees scatter by two orders of magnitude, so "average density `cR/a`" is not a density theorem for any given root.
  * The owner's refinement does not change the `0.95` exponent.
  * The reformulation shows the remaining work is pointwise and integral (`(i) AND (ii)`), not residue-theoretic. It does not supply that work.
* **Cross-checks.**
  * The Haar and averaged identities were verified three ways: by the `N_d` recursion, by explicit word enumeration, and by lattice sums.
  * The incongruent-class DP agrees with brute force.
  * The basin labels come from a vectorized census, and the renewal identity was checked from independently computed pure-Python labels at `10^6`.
  * A scratch recomputation (`scratch/procgen_tree/indep_check.py`, not a deliverable) agrees. It iterates `T_+` directly on `[-10^5, -1]`, giving basin fractions `0.3303, 0.3210, 0.3487` for `-1`, `-5`, `-17`, and recomputes the 3n-1 basin counts at `10^6` in pure Python.
  * The automata isomorphism was checked edge by edge.
  * No independent-agent audit yet.

## 6. Reproduction

```bash
cd <worktree>
sh 04-computation/experiments/collatz_procgen_20260924_tree_run.sh     # writes the .out (about 20 s)
# or individually (stdout = sections, stderr = timing and memory):
python3 04-computation/experiments/collatz_procgen_20260924_tree_automaton.py   # A1-A8
python3 04-computation/experiments/collatz_procgen_20260924_tree_controls.py    # C1-C4 (3n-1 census to 10^7)
python3 04-computation/experiments/collatz_procgen_20260924_tree_counts.py      # K1-K7 (tree densities at 10^7)
```

Peak memory is about 400 MB (the `counts` script). Each script ends with `ALL CHECKS PASSED`, and every check is an explicit `raise`.
