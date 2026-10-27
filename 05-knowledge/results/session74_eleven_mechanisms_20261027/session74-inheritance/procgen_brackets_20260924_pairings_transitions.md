# The pairing ladder: Collatz is a perfectly balanced up/down pairing of lengths, the graceful/Skolem structure cannot see tree-ness, and the square-sum problem shares its microcosm with the odd-square brackets, not with Collatz

**Status.**
- **PROVED** (elementary; each proof is in the text and each statement is also checked by code):
  - **(A) pair-sum invariant.** `T(2i-1) + T(2i) = 4i-1 = (2i-1) + 2i`. The product shrinks by exactly `3/4 + 1/(8i-4)`. Only `{1,2}` maps onto itself, and `sum_(n<=2M) T(n) = M(2M+1)`. Among the maps `(qn+-1)/2`, only `q = 3` admits a sum-preserving perfect matching, or any zero-sum partition into blocks of bounded diameter.
  - **(B) graceful form.** The down-edges `{i,2i}` and the up-edges `{2i-1,3i-1}` each realise every difference `1, 2, 3, ...` exactly once. The halving part on `[1,2n]` is a Skolem graceful forest. No Collatz subtree with two or more edges is graceful under its own labels, and a `k`-near-graceful one has at most `3k+3` edges.
  - **Escape sets.** The exact integer threshold for every multiplier. The owner-requested formula `floor((sqrt r+1)/(2(sqrt r-1)))` is the real-interval threshold; for integers it overshoots by one exactly on the windows `r in (g(m), h(m)]`. Every move of either shortcut map from `n >= 54` changes bracket. The complete list of `(k,p)` with `kp` in `p`'s bracket. The 13 skip-free primes.
  - **The pairing family.**
    - `3n-1` is the all-swapped Collatz pairing, conjugated by `n -> n+1`.
    - No periodic pairing admits a bounded-lookahead descent certificate: the length-0 pair creates a 2-adic repeller, the analogue of Applegate–Lagarias's `-1`.
    - Every "landing ⇒ down" member is a tree.
    - A density-zero set of flips produces a divergent orbit.
- **FINITE-EXACT.**
  - Pairing family: periodic census (`K <= 3` to `10^6`, `K = 4` to `10^4`), random census, and single flips to `10^7`. For `T`, 24 fragile pairs, all `<= 2308`; for `3n-1`, 12, all `<= 12029`.
  - Landing family: exact window optimum `0.2908`; the greedy member has flip density `1/3`.
  - Square-sum: existence and counts agree with OEIS A090461, A071983 and A071984. Threshold scaling across target sets.
  - Cycles of `(qn+-1)/2`; bracket prime counts to `3.6*10^9`.
- **CITED**: see section 7. The square-sum proof (Gerbicz 2018) is forum-only and was **not read** (login-gated); its mechanism is quoted from a README and from OEIS.
- **VERDICTS on the owner's analogy.**
  - "graceful tree : 3N+1" is an **ANALOGY**. The difference structure is real and exact, but it transfers nothing, and it is blind to tree-ness.
  - "square-sum : odd-square brackets" is **REAL**: the same square-density mechanism, and both microcosms end near 24–25.
  - "square-sum : Collatz microcosm" is **NUMEROLOGY**. Collatz-type exceptions are Diophantine, pinned to approximations of `log_2 3`, and they recur at every scale.
- **OPEN**: Collatz; a prime in every bracket (implied by Legendre, but out of reach of BHP and of RH); the graceful tree conjecture; the minimal flip density of `L`-step-provable pairings for `L >= 4`.

Session `collatz-procgen-20260923`, brackets lane `procgen_brackets_20260924` (mac-mini). Scripts: `04-computation/experiments/procgen_brackets_20260924_{escape_sets,pairings,provable_rung,square_sum,legendre,run}.py` with C helpers `..._pairings.c` and `..._squaresum.c`. Output: [procgen_brackets_20260924.out](procgen_brackets_20260924.out), about 7 minutes via `python3 04-computation/experiments/procgen_brackets_20260924_run.py`. Parents: [odd-square brackets](procgen_numerology_20260923_odd_square_brackets.md), [snippet dispatch](procgen_numerology_20260923_snippet_dispatch.md), [choice ladder](collatz_procgen_20260922_choice_ladder.md), [Q1 mirror](collatz_procgen_20260922_q1_mirror.md), [foundry](collatz_procgen_20260922_foundry.md). Independent cross-check: the census code's Collatz exceptional-class count at `2^20` is `27,328`, matching the foundry probe row.

## 0. Headline

`T(n) = n/2` (n even), `(3n+1)/2` (n odd); `|T(n) - n| = ceil(n/2)`. The pair `{2i-1, 2i}` shares the length `i`: `2i` steps down by `i` to `i`, and `2i-1` steps up by `i` to `3i-1`. For `3n-1` (`U`) the pairs are `{2i, 2i+1}`: `2i -> i` and `2i+1 -> 3i+1`, with `{0,1}` as the length-0 pair of two fixed points.

**(A) Pair-sum invariant (PROVED).**
- **Sums.** `T(2i-1) + T(2i) = (3i-1) + i = 4i-1 = (2i-1) + 2i`, and `U(2i) + U(2i+1) = 4i+1`. So each sheet preserves the sum of every consecutive pair of its own pairing: the two sheets are exactly the two ways of pairing consecutive integers.
- **Products.** `T(2i-1)T(2i)/((2i-1)2i) = (3i-1)/(4i-2) = 3/4 + 1/(8i-4)`. The arithmetic mean of each pair is preserved and the geometric mean shrinks by `sqrt(3/4 + 1/(8i-4)) -> sqrt3/2`. This is the exact combinatorial form of the AM–GM drift.
- **Fixed pairs.** `{i, 3i-1} = {2i-1, 2i}` forces `i = 1`, so the only pair mapped onto itself is `{1,2}`, the trivial cycle; for `U` it is `{0,1}`.
- **Corollary.** `sum_(n<=2M) T(n) = sum_(n<=2M) n = M(2M+1)`.
- **The coordinator's facts, all verified.** `|T(n) - n| = ceil(n/2)` (checked to `2*10^6`). For `3n-1` the pairing shifts to `{2i, 2i+1}`. With the consecutive pairing, the up-length of `2i-1` under `(qn+1)/2` is `((q-2)(2i-1)+1)/2`, which equals the down-length `i` of `2i` for all `i` iff `q = 3`. Equivalently, the arithmetic mean `(q+1)/4` of the step factors `q/2` and `1/2` is 1. The pair-sum invariant is the exact, pairwise form of "arithmetic-mean factor 1".
- **Only `q = 3` (PROVED).** For `(qn+-1)/2`, a sum-preserving pair can be neither two evens (displacements `-m/2 < 0`) nor two odds. For `q >= 3` the odd displacements `((q-2)n+-1)/2` are `>= 0`, zero only for `U(1) = 1`; for `q = 1` every displacement is `<= 0`. So the pair is `{odd n, even m}` with `m = (q-2)n +- 1`. This map is a bijection from the positive odds onto the needed evens only for `q = 3` (checked for all odd `q <= 99`). Moreover the displacement over `[1,2M]` equals `(q-3)M^2/2` on the `+` sheet (and `((q-3)M^2 - 2M)/2` on the `-` sheet). A partition into zero-sum blocks of diameter `<= D` would make it `O(DM)`, so none exists for `q != 3`. For `q >= 5`, a greedy "Hilbert hotel" partition into finite zero-sum blocks of unbounded diameter does exist (checked for `q = 5, 7` on `[1, 20000]`), so the bound on the diameter is essential.

**(B) Graceful form (PROVED).**
- **Two perfect difference systems.** The down-edges `D = {{i,2i}}` and the up-edges `C = {{2i-1,3i-1}}` each use every difference `1, 2, 3, ...` exactly once. They share only the edge `{1,2}`. So the Collatz graph is the union of two perfect difference systems, glued along the trivial cycle. The conjecture is equivalent to `D ∪ C` being connected, and then it is a spanning tree of the positive integers.
- **Truncations.**
  - On `[1,M]` every difference `d <= (M+1)/3` occurs twice (once up, once down) and every `d in ((M+1)/3, M/2]` once. Checked: `M = 10^2, 10^4, 10^6`.
  - The backward tree `R_M` of 1 inside `[1,M]` (the `n <= M` whose orbit stays `<= M`) has `|R_M|/M = 0.690, 0.611, 0.596` for `M = 10^2, 10^4, 10^6`. Of the differences `d <= M/2`, 35% occur twice, 50% once and 15% not at all (`M = 10^6`).
  - So a truncation is "graceful up to multiplicity 2": each difference occurs at most once per class. But its labels spread over about twice the range a graceful labelling allows.
- **Literal gracefulness is impossible.**
  - **Graceful.** Say a subtree with `m` edges is labelled by its own integers inside a window of length `m`. Difference 1 forces the edge `{1,2}`, the only edge of length 1. Difference `m` forces an edge from `2m-1` or `2m`. So `2m-1 <= 1+m`, i.e. `m <= 2`. Brute force over all intervals `[a,b]`, `a < 60`, `b < 200`, finds only `[1,2]`.
  - **`k`-near-graceful** (labels in a window of length `m+k`, differences distinct). Some difference `<= k+1` forces a vertex `<= 2k+2`, and some difference `>= m` forces a vertex `>= 2m-1`. Hence `m <= 3k+3`.
- **Skolem and Rosa.** The halving forest `H_n = {{i,2i} : i <= n}` on `[1,2n]` has `p = 2n` vertices labelled `1..2n` and `q = n` edges labelled `1..n`. So it is **Skolem graceful** in the sense of Lee–Shee (1991). A Skolem sequence of order `n` is the same object with `H_n` replaced by the perfect matching `nK_2`. By Rosa's translation argument, the `2n+1` translates of `H_n` mod `2n+1` decompose `K_(2n+1)` (checked for `n < 60`).

**The pairing ladder (the new axis).**
- **The family.** For every `i`, a bit `eps_i` chooses which member of the pair goes up: `n` goes up iff (`n` odd) XOR `eps_i`. Collatz is `eps = 0`.
- **The invariants are shared.** Every member `F_eps` preserves every pair sum and has the two perfect difference systems of (B).
- **Symmetries.** Shifting by one and complementing every bit exchanges the two offsets (`F^0_eps(n) + 1 = F^1_(1-eps)(n+1)`, checked). So `3n-1` is `eps = 1` shifted by one: **Collatz and `3n-1` are antipodal corners of the cube `{0,1}^N`**. The all-swapped Collatz pairing has the cycles `{4,6,9}` and `{16,...,135}`, the `3n-1` cycles shifted by one.
- **So (A) and (B) cannot decide tree-ness.** The family contains:
  - trees;
  - maps with extra cycles: the antipode, 24 single flips of `T`, and 64% of random members;
  - maps with divergent orbits: an explicit density-zero modification, and up-heavy periodic members.

| rung (all members satisfy (A) and (B)) | tree-ness | provability of tree-ness |
|---|---|---|
| uniform random `eps` (i.i.d. 1/2) | 36.3% trees (`n <= 10^5`, 300 maps); at most 5 cycles, minima `<= 416` | none |
| periodic `eps_i = c(i mod 2^K)` | 36.7% trees (`K = 3`); 29 of 256 have orbits beyond `4*10^18` | **never** by bounded lookahead (PROVED) |
| **Collatz** `eps = 0` | tree (conjecture; verified to `2^71`, cited in the choice ladder) | OPEN |
| `3n-1` (= `eps = 1`, shifted) | 3 cycles | — |
| Collatz + one flipped pair | a new cycle for exactly 24 pairs, all `<= 2308` (to `10^7`) | — |
| Collatz + a density-zero set of flips | divergent orbit | FALSE, PROVED |
| landing ⇒ down (every up-move lands on a down-mover) | tree | PROVED; needs `>= 29.08%` of pairs flipped (exact window), `1/3` for the explicit member |
| `L`-step descent, `L = 4, 5` | tree on the window only | window designs (not infinite members) with 14.5% and 5.9% flips |

So Collatz sits at density distance **zero** from a divergent map (PROVED). The price of an `L`-step proof falls quickly with `L`: a 29% floor at `L = 2` (exact), 14.5% and 5.9% window designs at `L = 4, 5`. This is only evidence that the price tends to 0.

**No property of the pairing itself (pair sums, the two difference systems, the density of the flip set) can decide the conjecture**: the divergent member has all of them. This is the foundry's DEFECT control ("planted density-zero modifications"), made concrete inside the one family where the owner's balance and gracefulness hold exactly.

Orbit statistics are different. A chain planted at scale `c_0` captures a fraction of orbits that falls roughly like `1/c_0`: 99.97%, 11.6%, 0.93% and 0.001% of `n <= 10^5` for `c_0 = 3, 27, 1001, 100003`. The chain from 3 runs through `8`, which carries almost every Collatz orbit. So almost-all theorems see a defect at a small scale but cannot exclude one planted at a large scale.

## 1. Escape sets for every multiplier (PROVED; tables FINITE-EXACT)

Brackets are `B_m = ((2m-1)^2, (2m+1)^2]`, with `B_0 = {1}`. `S_r = {n : n` and `rn` lie in the same bracket`}` (`rn` real).

**1.1 Thresholds.**
- **Real interval.** `((2m-1)^2, (2m+1)^2/r]` is nonempty iff `r < h(m) = ((2m+1)/(2m-1))^2`, i.e. `m < X(r) = (sqrt r+1)/(2(sqrt r-1))`.
- **Integers.** An integer `n in B_m` with `rn in B_m` exists iff `r <= g(m) = (2m+1)^2/((2m-1)^2+1)`, since the smallest candidate is `(2m-1)^2+1`.
- Both are strictly decreasing, and `g(m+1) < h(m+1) < g(m)`; the last inequality is `(2m+1)^4 - (2m+3)^2((2m-1)^2+1) = 28m^2+20m-17 > 0`. So `S_r` meets exactly the brackets `1..M_int(r)`, where `M_int(r) = max{m : r <= g(m)}`.
- **The requested formula.** `m_max(r) = floor((sqrt r+1)/(2(sqrt r-1)))` is the real-interval threshold, except that it overshoots by one when `X(r)` is an integer (`r = 9, 25/9, 49/25, ...`, where the interval degenerates). For integers it overshoots by exactly one iff `r in (g(m), h(m)]` with `m = m_max(r)`. Examples: every `r in (4.5, 9]` (formula 1, truth 0) and `r in (2.5, 25/9]` (formula 2, truth 1). On 4,410 rationals `p/q <= 10`, `q <= 40`, it is right for 1,989 and one too large for 2,421, all in those windows. A brute-force check confirms `M_int` for every tested `r > 1.1`.

| `m` | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|
| `h(m)` | 9 | 2.778 | 1.960 | 1.653 | 1.494 | 1.397 | 1.284 | 1.222 |
| `g(m)` | 4.5 | 2.5 | 1.885 | 1.620 | 1.476 | 1.385 | 1.279 | 1.218 |

**1.2 Tables (complete).**

| `r` | brackets | `S_r` | primes in `S_r` |
|---|---|---|---|
| 2 | 1–2 | `{2,3,4} ∪ {10,11,12}` | `{2,3,11}` |
| 3 | 1 | `{2,3}` | `{2,3}` |
| 4 (also 9/2) | 1 | `{2}` | `{2}` |
| 5 | none | ∅ | ∅ |
| 3/2 | 1–4 | `2–6, 10–16, 26–32, 50–54` | `2,3,5,11,13,29,31,53` |
| 5/3 | 1–3 | `2–5, 10–15, 26–29` | `2,3,5,11,13,29` |
| 4/3 | 1–6 | 50 numbers, max 126 | 12 primes, max 89 |
| 5/4 | 1–8 | 92 numbers, max 231 | 22 primes, max 229 |

**1.3 The actual Collatz steps, and the microcosm.**

| step | `n` with the image in `n`'s own bracket (complete) | proof |
|---|---|---|
| `(3n+1)/2`, n odd | `3,5,11,13,15,27,29,31,51,53` | `3((2m-1)^2+2)+1 <= 2(2m+1)^2 <=> m^2-5m+2 <= 0 <=> m <= 4` |
| `n/2`, n even | `4,6,8,20,22,24` (`= 2 S_2`) | `2(2m-1)^2+2 <= (2m+1)^2 <=> m <= 2` |
| `(3n-1)/2`, n odd | `1,3,5,11,13,15,17,27,29,31,33,51,53` | `2m^2-10m+3 <= 0 <=> m <= 4`; `1` is fixed |
| `3n+1`, `5n+1`, `(7n+1)/2` | ∅ | `m = 1` fails |
| `3n-1`, `(5n+-1)/2` | `{3}` | |

- **The microcosm (PROVED).** Every move of either shortcut map from `n >= 54` changes bracket; in particular every move inside a bracket `m >= 5` does, which covers all `n > 81 = 9^2`. For large `n` the bracket index is multiplied by `sqrt(3/2) = 1.22474` after an odd step and by `1/sqrt 2` after halving (measured `1.22475`, `0.70711`). The owner's `81` is thus the exact bracket boundary of the region where Collatz moves can stay inside a bracket.
- **Fully internal pairs.** The pairs whose up- and down-move both stay inside a bracket are exactly `i = 2, 3`: the pairs `{3,4}` and `{5,6}`, inside `B_1`. For `3n-1` it is only `{4,5}`.

**1.4 One function encodes every escape set, and the singleton structure.** Put `rho(n) = (2m(n)+1)^2/n`, the ratio of the next odd square to `n`. Then `S_r = {n : rho(n) >= r}`.
- **Escape order.** Listing `n` by decreasing `rho` gives the order in which numbers join as `r` decreases: `2 (4.5), 3 (3), 10 (2.5), 11 (2.273), 4 (2.25), 12 (2.083), 13 (1.923), 26, 27, 5 (1.8), 14, ...`.
- **Primes.** The prime escape sets are nested: `{2}` for `r in (3, 4.5]`, `{2,3}` on `(2.273, 3]`, and **`{2,3,11}` on `(25/13, 25/11] = (1.923, 2.273]`**. Then `13` (1.923), `5` (1.8), `29`, `31`, `53`, `17`, `83`, ... join. The owner's doubling `r = 2` sits inside the `{2,3,11}` window.
- **All `(k,p)` with `kp` in `p`'s bracket** (`k >= 2`, `p` prime, complete because `k <= rho(p) <= 4.5`): `(2,2), (3,2), (4,2), (2,3), (3,3), (2,11)`. For all integers `n` add `(2,4), (2,10), (2,12)`. In singleton language: **every prime outside `{2,3,11}` is the only multiple of itself in its own bracket.**
- **Skip-free primes.** These are the primes whose multiples meet every bracket from their own onwards: `2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 53` (13 primes). *Proof of completeness.* If `p in B_m`, the multiples of `p` up to `(4m+1)^2` number at most `(4m+1)^2/(2m-1)^2`. This is `< m+1`, the number of brackets `B_m..B_2m`, for every `m >= 5`. So every prime above 81 skips some bracket (e.g. `41 in B_3`, `82 in B_5`), and the finite check below 81 is complete.

**1.5 The sieve-window view (REAL, classical).** For `m >= 2` the primes of `B_m` are exactly the elements coprime to all primes `<= 2m+1`. The new sieving prime `q = 2m+1` first strikes `q^2`, the top of `B_m`, so each bracket is one exact round of Eratosthenes. The naive sieve count `8m prod_(p<=2m+1)(1-1/p)` overshoots by the classical factor `2e^(-gamma)`. Measured `#primes(B_m)/(naive)`: `0.8958, 0.8933, 0.8918` on `m in [100,500), [500,1500), [1500,4000)`, against `e^gamma/2 = 0.8905`. This is the one genuinely structural way "multiples of other numbers fall into the brackets"; its constant comes from Mertens' theorem, not from the brackets.

**1.6 Is `{2,3,11}` special beyond `m <= 2`? No (NUMEROLOGY), with one real but unrelated fact.**
- **Wieferich.** `11` is a base-3 Wieferich prime (`3^5 = 2*11^2 + 1`; below `2*10^6` the list is `11, 1006003`). The solutions of `b^(p-1) = 1 (mod p^2)` are the `p-1` roots of unity mod `p^2`, one per nonzero class mod `p`. So a fixed base is Wieferich for `p` with chance about `1/p`, a 1-in-11 event. Smallest Wieferich bases: `2 -> 5`, `3 -> 8`, `5 -> 7`, `7 -> 18`, `11 -> 3`, `13 -> 19`.
- **Trunk.** `3, 11` are the minus-trunk (Wagstaff) numbers `(2^(2n+1)+1)/3` inside `(1, 12.5]`. That window holds only the trunk values `1, 3, 11`, and a random pair of the odd primes `{3,5,7,11}` equals `{3,11}` with chance `1/6`.
- **Ljunggren (REAL, unrelated).** `(3^k-1)/2` is a square only for `k = 1, 2, 5` (Ljunggren 1943; checked `k <= 400`). So `B_5 = (81,121] = (3^4, (3^5-1)/2]` is the unique bracket bounded by a power of 3 and a base-3 repunit. The `E`-graph climb `1 -> 4 -> 13 -> 40 -> 121` meets an odd square `> 1` only at `11^2`. This concerns `11^2`, not `11 in B_2`, and no mechanism ties it to the escape property.
- **Figurate identities at the bracket tops (exact).** `T^2((2m+1)^2) = 3m^2+3m+1 = (m+1)^3 - m^3`: odd squares reach the centred hexagonal numbers in two steps, so they have stopping time 2. Likewise `U((2m+1)^2) = 6m(m+1)+1`, a star number.

## 2. The Skolem/graceful view and the pairing ladder

**2.1 Precise relation to the literature (CITED definitions; the relations are PROVED above).**
- **Skolem sequences** (Skolem 1957). The two copies of `k` sit `k` apart, for `k = 1..n`. They exist iff `n = 0, 1 (mod 4)`; O'Keefe (1961) did the hooked case. They are exactly the Skolem graceful labellings of the matching `nK_2`: the positions of equal entries are the edges.
- **Langford pairings** (Langford 1958) put the copies `k+1` apart. They exist iff `n = 0, 3 (mod 4)` and are the `(2,1)`-Skolem graceful matchings.
- **Graceful labellings** (Rosa's `beta`-valuations, 1967; the name is Golomb's). Labels are distinct in `[0,m]` and the differences are exactly `1..m`. **Skolem graceful** (Lee–Shee 1991) uses labels exactly `1..p`; a tree is Skolem graceful iff it is graceful. Skolem sequences are the standard tool for gracefully labelling 2-regular graphs and windmills (Abrham 1991; Alkasasbeh–Dyer–Howell, arXiv:2112.04265).
- **The Collatz case.** In the Collatz graph each class `D`, `C` is a perfect difference system, but not a matching: `D` is a forest of rays `u, 2u, 4u, ...`, `C` a forest of paths. `D ∩ [1,2n]` is Skolem graceful, while `C`'s labels spread over `[1,3n-1]`.
- **Harmonious labellings** (Graham–Sloane 1980) use sums mod `q`. The class sums are `3i` (down) and `5i-2` (up), distinct within each class, but they collide across classes (`{6,12}` and `{7,11}` both sum to 18). With unbounded labels no modular condition applies, so nothing more is available.
- **The two conjectures.** The **graceful tree (Ringel–Kotzig) conjecture** asks for *some* labelling of *every* tree. It is verified to 35 vertices (Fang 2010, as quoted in arXiv:2605.02303) and remains open. Its decomposition version (Ringel's conjecture) was proved for large `n` (Montgomery–Pokrovskiy–Sudakov 2020, per the Wikipedia summary; not read). The Collatz question asks whether *one fixed* labelled graph with the difference structure is a tree. No theorem in either direction is known or suggested, and sections 0 and 2.3–2.4 show that the structure itself cannot decide it.

**2.2 The family.**
- **Offsets.** Offset 0 has pairs `{2i-1,2i}` with length `ceil(n/2)`; offset 1 has pairs `{2i,2i+1}` with length `floor(n/2)` and the length-0 pair `{0,1}`.
- **Moves.** `n` goes up by its length iff (`n` odd) XOR `eps_i`. Each map is affine on each residue class: `(3n+(n mod 2))/2` up, `(n-(n mod 2))/2` down, with the signs flipped on the minus sheet.
- **Uniform cycle gate (PROVED).** Each step is `n -> (c n + e)/2` with `c in {3, 1}` and `|e| <= 1`. So a cycle with `a` ups and `K` steps satisfies `n(2^K - 3^a) = B = sum_t e_t 3^(u_t) 2^t`, where `u_t` counts the ups after step `t`. Since `u_t <= K-1-t`, `|B| <= K 2^(K-1) (3/2)^a`. So large cycles need `|2^K - 3^a|` small in every member, as for Collatz.
- **Offsets are conjugate.** `F^0_eps(n) + 1 = F^1_(1-eps)(n+1)`, so the offsets are conjugate after complementing every bit. Negation gives `T(-n) = -U(n)`.

**2.3 Census (FINITE-EXACT; `..._pairings.py` sections D–F).**
- **Periodic members, `eps_i = c(i mod 2^K)`, `2^(2^K)` maps per offset.**
  - *Measure preservation.* A member is Haar measure-preserving iff every class mod `2^K` has exactly two preimage classes: 8 of 16 at `K = 2`, 32 of 256 at `K = 3`, 128 of 65,536 at `K = 4`.
  - *Trees.* Offset 0: 1 of 2, 2 of 4, 8 of 16, **94 of 256 (36.7%)** and 23,492 of 65,536 (35.8%, `N = 10^4`) are trees. Offset 1 has fewer (12.5%, 13.3% and 11.1% for `K = 2, 3, 4`). There `1` is always a fixed point, and whenever pair 1 is swapped `2 <-> 3` is a second cycle, so half of those maps cannot be trees.
  - *Cycles.* At most 5 attracting cycles at `K = 3` (15 at `K = 4`). The largest cycle minimum found is 11,651 (92,326 at `K = 4`).
  - *Divergence.* **29 of 256** `K = 3` members have orbits beyond `4*10^18`. For the heavily escaping ones the up-frequency of Haar-random 2-adic orbits exceeds `log2/log3 = 0.631`: for mask 83 at least 52,385 of the `n <= 10^6` escape (up-frequency `0.720`); for masks 98 and 99 at least about 24,500 (`0.665`). These are `5n+1`-like divergent members. Non-escaping maps lie in `[0.335, 0.561]`. So **the pair-sum invariant controls the drift only for measure-preserving members.** A non-measure-preserving pairing can drift upward although it preserves every pair sum.
  - *The pair-0 obstruction (PROVED).* Near the 2-adic points `-1, 0` (offset 0), or `0, 1` (offset 1), the pair index is `= 0 (mod 2^K)`, so the bit is `c(0)`. Exactly one of the two neighbourhoods then rises: `x -> x + i` stays near the same point with one digit fewer, for as many steps as the digits allow. So every periodic member has an exceptional class at every level (measured minima `1,942` at `2^20`, `K = 3`; `14` at `2^12`, `K = 4`), and none satisfies landing ⇒ down (0 of all 131,628 periodic maps checked). **No periodic pairing can be proved a tree by bounded lookahead.** This is Applegate–Lagarias's "`-1` resists elimination", now forced by the length-0 pair itself.
  - *Thinner exceptional sets.* Some members have much thinner sets than Collatz: masks 39 and 54 count `41,326` classes at `2^27`, against `2,292,648` for Collatz. They are still exponential.
- **Random members (`eps_i` i.i.d.).**

  | flip prob. `p` vs Collatz | trees (300 maps, `n <= 10^5`) | mean attracting cycles | cycle minima: 90% / max |
  |---|---|---|---|
  | 0.5 (uniform) | 36.3% | 1.81 | 15 / 416 |
  | 0.1 | 47.7% | 1.70 | 19 / 447 |
  | 0.01 | 84.3% | 1.16 | 15 / 155 |
  | 0.001 | 97.0% | 1.03 | 1 / 91 |

  - *Cycles are small.* No member had an escaping orbit. `97.0%` at `p = 0.001` is about `(0.999)^24`: tree-ness is decided by the 24 fragile small pairs.
  - *Long cycles sit at approximations of `log_2 3`.* All 13 cycle types of length `>= 20` found (lengths 22 to 149; `p = 0.5, 0.1, 0.01`) have `|K - a log_2 3| <= 0.23`, i.e. `2^K` and `3^a` agree within a factor 1.17. Examples: `65/41`, `84/53`, `149/94`.
- **Single flips (the microcosm of the pairing).**
  - *`T` (flip one pair of Collatz, `i <= 10^7`).* A new cycle appears for **exactly 24 pairs**: `i = 1, 4, 5, 10, 11, 13, 20, 22, 40, 61, 84, 122, 126, 167, 189, 217, 244, 325, 334, 433, 445, 577, 1154, 2308`. None lies in `(2308, 10^7]`.
    - Examples: `i = 1` gives `1 -> 0` and the cycle `(2 3 5 8 4)`; `i = 4` gives `(3 5 8 12 6)`.
    - The new cycles have `(a, K) = (3,5), (5,8), (10,16), (17,27), (29,46), (34,54), (46,73)`. All are **upper** approximations of `log_2 3` (`2^K > 3^a`).
  - *`U = 3n-1`.* 12 fragile pairs, `i = 1, 2, 20, 30, 41, 43, 45, 86, 205, 365, 410, 12029`. The cycles have `(a, K) = (1,2), (2,3), (7,11), (12,19), (53,84)`: after the trivial `(1,2)`, all are **lower** approximations. The last, at `i = 12029`, is the convergent `84/53`: a cycle of length 84 through `1253..158284`.
  - *Consistency with Q1/Q2.* This matches the Q1-mirror dichotomy (Q2's upper best approximations, Q1's lower ones).
  - *Mechanism.* After flipping pair `i`, a new cycle is a path `3i ~> 2i` (or `i-1 ~> 2i-1`) of the old map, i.e. `i (2^(K+1) - 3^(a+1)) = B_w`. This is the cycle gate again, so fragility dies out like cycles do. It can recur at later convergents, which is why no finiteness is claimed.

**2.4 The provable rung (`..._pairings.py` G, `..._provable_rung.py`).**
- **Landing ⇒ down (PROVED to give trees).** Suppose every up-move lands on a down-mover. Then `n -> ceil(3n/2) -> floor(ceil(3n/2)/2) < n` for `n >= 2`, so every orbit reaches `{1,2}` (or `0` if `eps_1 = 1`).
- **The constraints.** The rule is equivalent to: `eps_3k = 1 - eps_2k` for `k >= 1`; `eps_(2k+1) = 0 ⇒ eps_(3k+1) = 0`; `eps_(2k+1) = 1 ⇒ eps_(3k+2) = 1`.
- **Collatz violates it at every even `i`.** There `3i-1` is odd, which gives the rising runs.
- **The explicit member.** Free bits set to 0 give a recursion. Its flip density is **`1/3`** (`0.333325` at `2^24`; tree verified to `10^6`, and two-step descent for every `n <= 10^6`). The flips are `1/6` on indices with even `v_3` and `5/6` on odd `v_3`.
- **Lower bounds.** The `{2k,3k}` constraints alone need density `>= 1/4`: `{i : v_3(i) odd}` meets each edge `{2k,3k}` exactly once and is a minimum vertex cover (sizes agree to `10^7`). The exact minimum over the whole landing family on a window, with sources `i <= X` and all constraints (CP-SAT, OPTIMAL), is **`0.29133, 0.29097, 0.29080`** for `X = 3000, 3*10^4, 10^5`. So every provable-by-landing pairing differs from Collatz on at least about 29% of the pairs in `[1, 10^5]`.
- **`L`-step descent.** The class `P_L` asks that every `n >= 3` falls below itself within `L` steps.
  - First descents occur at the Terras lengths `floor(1 + a log_2 3) = 1, 2, 4, 5, 7, 8, 10, 12, 13, ...`. The lengths 6, 9 and 11 occur only for `n <= 19`, so `P_6 = P_5` away from small `n`.
  - Collatz itself fails `P_4` exactly on `n = 7, 11, 15 (mod 16)`, i.e. on `3/16` of all `n`; the re-check code reproduces this count.
  - CP-SAT window designs on `[1,3000]` (FEASIBLE, not proven optimal) flip **14.5%** at `L = 4` and **5.9%** at `L = 5`.
  - A sequential greedy that certifies `n = 3, 4, ...` in turn and freezes what it uses flips 38.8%, 24.2%, 15.7%, 12.4% and 9.0% at `L = 2, 4, 5, 7, 8` (`n <= 2*10^4`; stable to `10^5`). An independent re-check confirms that every `n <= N` descends within `L` steps under the final choice, except the one `n` the greedy could not certify at `L = 8`.
  - The window designs are not infinite members, so "the price tends to 0" is **evidence, not a theorem**.
- **The falsifying rung (PROVED).** Take `c_0 = 3`, `c_(t+1) = c_t + ceil(c_t/2)`, and flip the pair of every even chain member. The chain pairs are distinct, so the modified map follows `3, 5, 8, 12, 18, 27, 41, ...` upward forever (verified to `10^30`: 82 flips among 169 chain terms). The flip set has at most `log_(3/2) X + 1` elements below `X`. Every pair sum and every difference system is untouched.
  - The same construction works from any `c_0`. The captured fraction of `n <= 10^5` falls roughly like `1/c_0` (section 0).
  - This is the positive-density-orbit, density-zero-pairing contrast behind the verdict that pairing statistics are DEFECT-blind.

**2.5 Where Collatz sits, and what would change it.**
- **Falsified by:** one flip at one of 24 small pairs, or a density-zero infinite flip set.
- **Made provable by:** a positive-density aperiodic change. It is never periodic (pair-0 obstruction). Its price is about 29% at `L = 2`, with window designs of 14.5% and 5.9% at `L = 4, 5`.
- Collatz's own violations of `P_L` are the `n` with stopping time `> L` (for `L = 4`: `n = 7, 11, 15 (mod 16)`). They are dominated by the 2-adic neighbourhood of `-1`, the long rising runs, which is the neighbourhood the choice ladder singled out. There the `6 mod 8` excursion alone recovers 96.5% of the collapse; here the landing rule must intervene at every even `i`.
- **The two ladders are complementary axes.** Choice adds nondeterminism and collapses the exceptional set. Pairing keeps determinism and every balance statistic, and changes the truth value at density zero. Neither carries the sign or the drift.
- **Candidate hypotheses** (for numbering by the coordinator; no files created):
  - (P1) The minimal flip density `delta_L` of `P_L` near Collatz tends to 0 as `L -> infinity`.
  - (P2) `delta_2 = 0.2907...` exists as a limit.
  - (P3) Collatz has finitely many fragile pairs. Evidence: none in `(2308, 10^7]`. Against: the `84/53` example for `3n-1` shows late ones can occur.

## 3. The square-sum problem

**3.1 Literature status.**
- **Chains.** Chains (Hamiltonian paths of `G_n`: vertices `1..n`, `a ~ b` iff `a+b` is a square) exist exactly for `n = 1, 15, 16, 17, 23` and all `n >= 25`.
- **Loops.** Loops exist exactly for `n >= 32`. For `n <= 30` a vertex has degree `<= 1` (18's only partner is 7), and `n = 31` was done by hand (Dobbelaere, in A071984).
- **Proof status.** The all-`n` statement is **claimed PROVED** by R. Gerbicz, Mersenneforum, January 2018 (OEIS A090461: "The conjecture has been proved: every k >= 25 is in the sequence, moreover for k >= 32 there is a Hamiltonian cycle"). The argument is elementary plus computer-checked tables. It builds "nice pairs" (a chain of length `n` and one of length `n+1` whose common entries occupy positions of equal parity) for `49n + r`, `r = 24..72`, from those for `n`, and covers `n <= 2032` by lookup tables. The proof was popularised by a 2023 video (HexagonVideos, per OEIS); the thread is login-gated and was **not read**. No refereed publication was found (arXiv search negative). Status: **PROVED per an unrefereed forum proof with published code**, consistent with every computation here.
- **Counts.** Counts are in A071983 (chains) and A071984 (loops). Cubic loops first exist at `n = 473` (Rivera's Puzzle 311, via OEIS).

**3.2 Independent data (FINITE-EXACT).**
- **Existence.** Exhaustive search reproduces the chain set for `n <= 60` and loops exactly for `32 <= n <= 60`. The chain counts for `n = 15..38` equal A071983 term by term, and the loop counts for `32..38` equal A071984.
- **Obstructions.** Isolated vertices for `n = 2..6`. At least three leaves for `n = 7..14, 18`. Exhaustive search only for `n = 19..22, 24` (paths) and `31` (loop).
- **Growth.** `log(#chains)/n` rises through `0.09, 0.14, 0.19, 0.24, 0.26` at `n = 25, 28, 31, 34, 37`.

**3.3 The transition quantified: a degree threshold, powered above by exact self-similarity.**
- **Degrees.** A vertex near `n` has about `c n^alpha` partners: `c = sqrt2 - 1`, `alpha = 1/2` for squares. Measured minimum/mean degree is `1/2.56`, `2/3.72`, `4/5.42`, `8/11.0` at `n = 25, 50, 100, 400`.
- **Mechanism test 1 (target sets).** The last loop exception sits where the typical degree reaches about 2.3:

  | target | `c` | last `n` without chain (last degree obstruction) | last `n` without loop (last degree obstruction) | `c n^alpha` at the loop threshold |
  |---|---|---|---|---|
  | squares | 0.414 | 24 (18) | 31 (30) | 2.34 |
  | triangular `k >= 2` | 0.586 | 8 (8) | 14 (11) | 2.27 |
  | squares + doubled squares | 0.707 | 12 (9) | 20 (15) | 3.24 |
  | cubes (literature) | 0.260, `alpha = 1/3` | — | 472 | 2.03 |

  Thinning by congruence classes disconnects the graph instead. For example, odd-square sums are `1 (mod 8)` and split `G_n` into four residue-pair components, so no chain exists for any `n >= 8`. This is why the test uses unions and triangular numbers.
- **The macrocosm engine (PROVED identity).** If `a + b = s^2` then `(49a + c) + (49b - c) = (7s)^2`. From a chain of `[1,25]`, the 49 alternating chains `(49a_1 + c, 49a_2 - c, ...)`, `c = -24..24`, are square-sum chains that partition `[25, 1249]` exactly (checked). Gerbicz's gluing of these 49 chains with `1..24` is the part not re-verified. So the square-sum problem has an **exact** microcosm-to-macrocosm map, scaling by the square 49 with alternating offsets. Collatz has no known analogue.

**3.4 Comparison with Collatz-type small-number phenomena.**

| phenomenon | last exception | governing quantity | mechanism |
|---|---|---|---|
| square-sum chains / loops | 24 / 31 | degree `~ 0.414 sqrt n` | square density (degree threshold) |
| bracket halving escape `n/2` | 24 | odd squares in `(n/2, n]` | square density (exact inequality) |
| bracket-internal Collatz moves | 53 (brackets `<= 4`, `n <= 81`) | `((2m+1)/(2m-1))^2` vs `3/2` | archimedean ratio |
| `3n-1` cycles (= negative `3n+1` cycles) | minima 1, 5, 17; max 136 | `(a,K) = (2,3), (7,11)` | `|2^K - 3^a|` (Diophantine) |
| E-game rigid precisions (Q1 mirror) | `K = 5..9` | no `a0`-loop through `-1` for `K <= 9` | Diophantine + carries |
| single-flip fragility, `T` / `U` | 2308 / 12029 | `(a,K)` at approximations of `log_2 3`; `84/53` | Diophantine |
| cycles of `(qn+-1)/2`, `q <= 39` and 181 | e.g. `q = 181`: minima 27, 35, elements to 55,296 | `2^15 - 181^2 = 7` | Diophantine |

- **Mechanism test 2.** Every nontrivial cycle of `(qn+-1)/2` with minimum `<= 20000` sits at a near-coincidence `2^K ~ q^a`. Examples: relative gap `2.3*10^(-2)` for the `5n+1` cycles, `2.1*10^(-4)` for `q = 181`. The size of the "microcosm" jumps with continued-fraction accidents rather than scaling smoothly (`q = 1093`, a base-2 Wieferich prime, has no cycle with minimum `<= 10^6`).
- **The two classes behave differently.** Square-type thresholds scale smoothly with the target density (`c n^alpha ~ 2-3`) and end for good. Collatz-type ones are pinned to `log_2 3` approximations and can reappear at later convergents: the `84/53` fragility at `i = 12029` came after a gap from 410.
- **The one Collatz microcosm of square type is trivial.** It is the bracket-internal region (an archimedean inequality) and carries no orbit information.
- **The shared number 24.** The bracket halving escape set ends at 24 and the last square-sum chain exception is also 24. Both belong to the "`(sqrt2-1) sqrt x ~ 2`" regime: `(sqrt2-1) sqrt x > 2` iff `x > 23.3`, and the last `x` with at most one square in `(x,2x]` is 17. The exact coincidence comes from different sub-mechanisms (an interval inequality against an exhaustive-search exception), so the **order of magnitude is REAL; the equality 24 = 24 is NUMEROLOGY**.

## 4. Legendre-type statements for the brackets

**4.1 Status.**
- **Legendre.** Each bracket contains two Legendre intervals `((2m-1)^2, (2m)^2]`, `((2m)^2, (2m+1)^2]`, so Legendre's conjecture gives at least 2 primes per bracket. Oppermann's conjecture splits it into four intervals `((2m-1)^2, (2m-1)2m], ..., (2m(2m+1), (2m+1)^2]` and gives at least 4.
- **Brocard.** For odd primes `p < p'`, the interval `(p^2, p'^2]` is a union of whole brackets. So Brocard's conjecture (at least 4 primes between consecutive prime squares) is a bracket statement. It follows from Oppermann, but not from Legendre at twin primes, where the union is a single bracket.
- **Unconditional results do not reach.** Brackets have length `8m = 4 sqrt x`. Baker–Harman–Pintz (2001) give primes in `[x - x^0.525, x]` for large `x`, with `x_0` not explicit. `x^0.525 > 4 sqrt x` once `x > 4^40 = 1.2*10^24`, so BHP never implies a prime in every bracket.
- **Neither does RH.** Under RH, Cramér's gap bound `O(sqrt x log x)` is still a logarithm too long.
- **Almost-all results.** Bazzanella (2000) bounds the exceptional set of Legendre unconditionally and under RH (abstract read; exponents not checked). Under Lindelöf, all but `O(N^eps)` of the intervals `[n^2,(n+1)^2] ⊂ [1,N]` have the expected number of primes (Bazzanella 2013, abstract read). So brackets inherit "all but `O(N^eps)`" under Lindelöf.
- **Verification.** Legendre is verified to `n^2 ~ 2*10^19` via maximal prime gaps (Oliveira e Silva–Herzog–Pardi 2014, per Wikipedia; not read), so every bracket below `2*10^19` holds at least 2 primes.
- **The cube analogue is PROVED:** Ingham (1937) for large `n`, and Dudek (2016) explicitly.

**4.2 Data (FINITE-EXACT, `m <= 30000`, odd squares to `3.6*10^9`).**
- Brackets 1–3 hold 4, 5 and 6 primes. From `m = 4` on, every bracket holds at least 7; the minimum over `m >= 2` is 5 (at `m = 2`).
- Each Legendre half holds `>= 2`, and each Oppermann quarter `>= 1`; for `m >= 100` the quarter minima are `16, 15, 16, 13`.
- Brocard unions hold `>= 5`.
- `count/(4m/log(2m+1))` has mean `1.0000` and minimum `0.8835` (`m > 100`).

**4.3 Links to Collatz orbits.**
- **Provable, and trivial:**
  - the microcosm theorem of section 1.3;
  - the figurate identities of section 1.6;
  - the bracket index is multiplied by `sqrt(3/2)` or `1/sqrt 2` per step, so Collatz is equivalent to "the bracket index of every orbit reaches 1".
- **No substantive link, as a PROVED meta-statement.** Every pairing map of section 2 has the same brackets and the same primes, yet the family contains trees, cycles and divergent orbits. Hence no statement about primes in brackets can distinguish Collatz from a density-zero mutant with a divergent orbit, and none can imply the conjecture by an argument that holds across the family.
- **Negative control.** The detrended mean total stopping time of `B_m` against the detrended normalised prime count of `B_m`, for `11 <= m <= 1000`, correlates at `r = -0.0056` (noise level `0.064`).

## 5. Verdicts on the owner's texts

| claim | verdict | test |
|---|---|---|
| `{2,3,11}` are the only primes whose double lies below the next odd square | PROVED (inherited), sharpened: `{2,3,11}` is the prime escape set for every `r in (25/13, 25/11]`; the full `(k,p)` list; singleton reading | section 1 |
| "11 and 2,3" are special beyond `m <= 2` | NUMEROLOGY | Wieferich base 3 is a 1-in-11 event; the trunk match has chance 1/6; Ljunggren's `11^2 = (3^5-1)/2` is REAL but concerns `11^2`, not the escape |
| mapping multiples into brackets | REAL but classical | brackets are Eratosthenes windows; the Mertens factor `e^gamma/2 = 0.8905` measured `0.892`; skip-free primes finite (13) |
| Collatz is a balanced, graceful-type pairing | REAL, exact: (A) pair sums, (B) two perfect difference systems, `q = 3` unique | section 0 |
| "graceful tree conjecture : 3N+1" | **ANALOGY** | the structure is shared by the whole pairing family, which contains trees, extra cycles and divergence, so it carries no information on tree-ness; literal gracefulness impossible beyond `m <= 3k+3` |
| "square-sum : [its partner]": start one way, transition to another | **REAL for the odd-square brackets** (both are square-density thresholds near 24–25); **NUMEROLOGY for Collatz's own microcosm** | threshold scaling `c n^alpha ~ 2.3` across targets vs Diophantine pinning of every Collatz-type exception (cycles, fragility, rigid precisions, `qn+1`) |
| "complex microcosm–macrocosm transition" | REAL for square-sum: an exact 49-scaling engine drives the macrocosm; for Collatz the macrocosm is statistical (AM–GM drift) and the microcosm recurs at each convergent | sections 2.3, 3.3 |

## 6. Reproduction

`python3 04-computation/experiments/procgen_brackets_20260924_run.py` (about 7 minutes, peak memory `< 700 MB`, CP-SAT with 2 workers) writes `procgen_brackets_20260924.out`.
- Part scripts: `..._escape_sets.py` (section 1), `..._pairings.py` + `..._pairings.c` (section 2), `..._provable_rung.py` (section 2.4), `..._square_sum.py` + `..._squaresum.c` (section 3), `..._legendre.py` (section 4).
- Random maps use `splitmix64` hashing with fixed seeds, so they are deterministic. The time-limited CP-SAT designs at `L = 4, 5` are FEASIBLE solutions and can vary slightly between runs.

## 7. Sources

**Read this session:**
- OEIS A090461, A071983, A071984 and A090460 entries (JSON).
- The README of github.com/AlasdairWilkins/reversed-square-sum, quoting Gerbicz's forum post.
- Wikipedia: Langford pairing / Skolem sequence, Graceful labeling, Graph labeling, Legendre's conjecture (and its reference list).
- Abstracts or bibliographic records:
  - Lee–Shee, "On Skolem graceful graphs", *Discrete Math.* 93 (1991) 195–200;
  - Abrham, *Discrete Math.* 93 (1991) 115–121;
  - Skolem, *Math. Scand.* 5 (1957) 57–68;
  - O'Keefe, *Math. Scand.* 9 (1961) 80–82;
  - Bazzanella, *Arch. Math.* 75 (2000) 29–34, and the Lindelöf paper (*Period. Math. Hungar.*, 2013);
  - arXiv:2112.04265 (Alkasasbeh–Dyer–Howell);
  - arXiv:2605.02303 (Niu; quotes Fang 2010).

**Not read (cited from standard knowledge or secondary summaries):**
- Langford, *Math. Gazette* 42 (1958) 228.
- Rosa (1967), `beta`-valuations.
- Graham–Sloane, *SIAM J. Algebraic Discrete Methods* 1 (1980) 382–404.
- Montgomery–Pokrovskiy–Sudakov (2020).
- Baker–Harman–Pintz, *Proc. LMS* 83 (2001) 532–562.
- Ingham, *Quart. J. Math.* 8 (1937) 255–266.
- Dudek, *Funct. Approx.* 55 (2016) 177–197.
- Oliveira e Silva–Herzog–Pardi, *Math. Comp.* 83 (2014) 2033–2060.
- Ljunggren (1943).
- Gerbicz's Mersenneforum thread (login-gated).
- Applegate–Lagarias, "The 3x+1 semigroup" (via the choice-ladder note).
