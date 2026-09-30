# The mediant tree and the set `K`, unified through the prefix fixed points: the periodic point of a prefix is the 2-adic convergent of `n` and the threshold of THM-4512, the orbit is the clock times the gap to it (the observer identity), the two-place product identity, the precision residual is empty to `10^6`, the 2-adic distance to the never-descending set is `2^{-⌊τ log_2 3⌋}`, the record function and a Terras-minimum conjecture (a Liouville inequality against a 2-adic Cantor set), the mediant tree as the tree of all `3x + d` cycles, and the `5x+1` control

**Session:** opus, `collatz-posets-zeta5-20260927` (S15, fourteenth note), 2026-09-30.
**Owner's directive:** "pursue the mediant tree and the set K and any other
creative and useful syntheses you can come up with after exploring seemingly
unrelated past work in the repo extensively, especially honing in on novel
math we have done and our own unique ways of looking at problems."

**Status: PROVED (elementary) for Theorems M1, M2, F and Lemma K1;
FINITE-EXACT for every census (the residual to `10^6`, the records to
`2·10^6`, the denominators to `100`, the `5x+1` control); CONJECTURE G stated
with data and typed against the barrier atlas; CITED for THM-4512, THM-4476,
THM-4478, THM-4495, the first, fourth, sixth, twelfth and thirteenth notes,
the precision-residual note, META-PATTERNS, THM-972/998 (analogy only) and
Lagarias 1990 (from memory); Collatz OPEN. Independent audit OWED (weekly
subagent limit); self-checks: 300 random `(n, j)` for the identities, the
residual and the record census by two loops, the `3x+d` cycles by direct
orbit search against the tree.** Script
`04-computation/experiments/collatz_mediant_tree_K_20260930.py` (outputs `.out`,
`_p5.out`).

**What the repo already had, and what is new here.** The threshold
`N(w) = S/(2^A - 3^j)` with "descent on the exact cylinder iff `n > N(w)`" is
in the [precision-residual note](collatz_precision_residual_20260926.md) §2
and THM-4512. The two-place corollary (real Bernstein value against 2-adic
value, `R(d) < n` for divergent orbits) is THM-4476. The universal size price
`n ≥ ρ_K(x)` is Theorem 1 of the [sixth note](collatz_generic_price_topology_20260927.md).
The weighted-mediant law is Theorem C1 of the
[thirteenth note](collatz_pentagon_fixed_points_argument_styles_20260930.md).
New: (i) `N(w)` is the *fixed point* `x_w` of the prefix word, hence a
weighted mediant of the letter points and a periodic point of the Collatz map
on `Z_(2)`; (ii) `x_{w_j} ≡ n (mod 2^{A_j})`, the prefix fixed points are the
2-adic convergents of `n`, and the orbit displacement is the clock times the
gap (Theorem M1); (iii) the residual is empty to `10^6` and the reason is
Diophantine; (iv) `dist_2(n, K) = 2^{-⌊τ log_2 3⌋}` (Lemma K1), which turns the
record stopping indices into a Liouville-type statement and gives Conjecture G;
(v) the tree of fixed points is the tree of all `3x + d` cycles (Theorem F);
(vi) the `5x+1` control.

## 1. The observer identity: prefix fixed points are the 2-adic convergents of `n` (PROVED)

For an odd integer `n` write `w_j = (v_1, ..., v_j)` for the first `j`
valuations of its Syracuse orbit `m_0 = n, m_1, ...` (`m_i = U^i(n)`),
`A_j = v_1 + ... + v_j`, `S_j` the carry, `D_j = 2^{A_j} - 3^j` the clock of the
prefix, and `x_{w_j} = S_j/D_j` its fixed point (the periodic point of the
Collatz map on `Z_(2)` with word `w_j` repeated, the weighted mediant of the
letter points along `w_j`).

**Theorem M1.** For every `j ≥ 1`:
1. `x_{w_j} ≡ n (mod 2^{A_j})` as 2-adic integers (`D_j` is odd);
2. `U^j(n) - n = D_j (x_{w_j} - n) / 2^{A_j}`;
3. `|x_{w_j} - n|_2 · |x_{w_j} - n|_∞ = |U^j(n) - n|_2 · |U^j(n) - n|_∞ / |D_j|`.

*Proof.* The Terras identity `2^{A_j} U^j(n) = 3^j n + S_j` gives
`S_j - n D_j = 2^{A_j}(U^j(n) - n)`; divide by `D_j` for (2); reduce modulo
`2^{A_j}` for (1); (3) is (2) with both absolute values of a rational number
whose denominator is odd. ∎ (Checked on 300 random `(n, j)`; for `n = 27` the
first eight prefix fixed points are `-1, -5, -23/11, -85/49, -287/179,
-925/601, -2903/1675, -9221/4513`, all negative because the first 36 prefixes
of 27 have negative clocks; the first positive-clock prefix is `j = 37`, where
27 descends.)

**Reading.** This is the repo's observer lens ("keep observer, recurrence
class, and finite head", META-PATTERNS): `n` is the marked observer, the
finite head is the prefix, and the recurrence class is represented by the
periodic point that shares the head. Every displacement of the orbit is the
clock of the head times the gap between the observer and that periodic
point. Along `j` the periodic points converge to `n` 2-adically at rate
`2^{-A_j}` and, when the orbit descends (clock `≈ 2^{A_j}`), they track the
orbit archimedeanly: `x_{w_j} ≈ U^j(n)`. Item (3) is the same shape as the
first note's product identity `|ξ - r_L|_∞ |r_L|_2 = η_L/3^L` and THM-4476's
two-place corollary, now for the convergents themselves: the product of the
two sizes of the approximation error is the product of the two sizes of the
displacement, divided by the clock.

## 2. THM-4512's threshold is the prefix fixed point; the precision residual is empty (PROVED + FINITE-EXACT)

**Theorem M2.** Let `w` be a prefix of `n`'s word with clock `D_w`. If
`D_w < 0` then `U^j(n) > n`. If `D_w > 0` then `U^j(n) < n` iff `n > x_w`,
`U^j(n) = n` iff `n = x_w` (then `n` is a cycle point with word `w`), and
`U^j(n) > n` iff `n < x_w`. *Proof.* Theorem M1(2), `D_w`'s sign. ∎ So the
threshold `N(w)` of the precision-residual note and THM-4512 is the fixed
point of `w`, and the three sets of that note read: word-level descent at
`j` is `D_{w_j} > 0`; value-level descent is `D_{w_j} > 0` and `n > x_{w_j}`;
the *precision residual* is the set of `n` that lie **at or below the fixed
point of their own prefix** while its clock is positive: integers hiding
under a rational cycle point of `3x + d` (§4) that shares their head.

**Census (FINITE-EXACT).** For every odd `3 ≤ n ≤ 10^6`, at the first prefix
with positive clock, `n > x_{w}`: the residual is empty, and the first
word-level descent index `τ(n)` equals the first value-level descent index
`σ(n)` (THM-4512 reports the same coincidence to `10^7` in its notation).
The only integer at its own threshold is `n = 1` (`w = (2)`, `x = 1`).

**Why the residual is a Diophantine near-criticality effect (PLAUSIBLE).**
At the first positive-clock prefix the clock is `D = 2^A - 3^j` with
`2^A ∈ (3^j, 2^v 3^j]`; the carry is `S ≤ Σ_{i<j} 3^{j-1-i} 2^{d_i} ≤ c 3^{j-1}`
with `c = Σ_i 2^{d_i}/3^i = O(1)` typically (the prefix is no-descent, so each
ratio is `≤ 1`). Hence `x_w = S/D ≈ c/(3θ)` with `θ = 2^A/3^j - 1 ∈ (0, 2^v - 1]`:
`x_w` is `O(1)` unless the crossing is near-critical. Rhin's bound
`|A log 2 - j log 3| ≥ A^{-13.3}` gives `θ ≥ A^{-13.3}`, so `x_w ≤ c A^{13.3}`,
while a residual member must be a positive integer `≤ x_w` in an exact class
modulo `2^{A+1}`. With `2^{hA}` no-descent classes and least representatives
spread over `[0, 2^{A+1})`, the expected number of residual members at depth
`A` is `≲ A^{13.3} 2^{-(1-h)A}`, summable in `A`: the residual is expected to
be finite, and empirically it is `{1}`. This is the size price of the sixth
note doing the work: a small integer cannot afford a deep near-critical head.

## 3. The set `K`: distance, records, and a Terras-minimum conjecture

Let `K ⊂ Z_2` be the odd 2-adic integers whose word never descends
(`2^{d_j} ≤ 3^j` for all `j`; the sibling ladder's `Bad(3)`, Hausdorff
dimension `h = h(log_3 2) = 0.9499555`). For odd `n ∉ K` let `τ(n)` be the
first `j` with `2^{d_j} > 3^j`.

**Lemma K1.** `dist_2(n, K) = 2^{-⌊τ(n) log_2 3⌋}`.
*Proof.* A point of `K` cannot share the first `τ` valuations of `n` (they
descend). It can share the first `τ - 1` and have `τ`-th valuation `v' ≤ v^* :=
⌊τ log_2 3⌋ - A_{τ-1}` (append ones afterwards: `2^{A+m} ≤ 3^{τ+m}` stays true).
The set of odd `x` with the first `τ - 1` valuations of `n` and `τ`-th
valuation `≥ v^*` is one class modulo `2^{A_{τ-1} + v^*}` (the exact word needs
one extra bit, and `v_2(3x+1) ≥ v^*` is a class modulo `2^{v^*}`; `v^* ≥ 1`
because the head is no-descent), and `n` lies in it since `v_τ(n) > v^*`. So the
closest points of `K` are at distance exactly `2^{-(A_{τ-1} + v^*)} =
2^{-⌊τ log_2 3⌋}`. ∎ (`n = 1, 3, 7, 27` give exponents `1, 3, 6, 58`.)

**Records (FINITE-EXACT, odd `n < 2·10^6`).** The record holders of `τ`:

| `n` | 3 | 7 | 27 | 703 | 10087 | 35655 | 270271 | 362343 | 381727 | 626331 | 1027431 | 1126015 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `τ` | 2 | 4 | 37 | 51 | 66 | 85 | 103 | 104 | 109 | 111 | 115 | 141 |
| `⌊τ log_2 3⌋` | 3 | 6 | 58 | 80 | 104 | 134 | 163 | 164 | 172 | 175 | 182 | 223 |

So `dist_2(n, K) ≥ n^{-κ}` holds for all `n < 2·10^6` with `κ = 12.2`, the
worst case being `n = 27`; the ratio `τ/log_2 n` at the last records is
`5.6–7.0`. Define the **record function** `g(A) = min{n ≥ 1 odd : dist_2(n, K)
≤ 2^{-A}}`, the least positive representative among the no-descent residue
classes of depth `A` (Terras's classes); it is non-decreasing, `g(A) = 27` for
`7 ≤ A ≤ 58`, `g(223) = 1126015 = 2^{20.1}`, and `log_2 g(A)/A` runs
`0.08, 0.12, 0.13, 0.11, 0.11, 0.09` at the record depths `58, 80, 104, 134,
163, 223`.

**Two exact equivalences.** `K ∩ N = ∅` (no positive integer diverges or lies on
a non-trivial cycle's minimum, i.e. the Collatz conjecture over `N` apart from
the cycle `{1}`) iff `g(A) → ∞`; and `τ(n) ≤ C log_2 n` for all `n` iff
`g(A) ≥ 2^{A/(C log_2 3)}` for all `A`. (Both are rewritings of the
definitions through Lemma K1.)

**Conjecture G (Terras minimum).** `liminf_A log_2 g(A)/A ≥ 1 - h = 0.0500445`:
the least positive representative of a no-descent class of depth `A` is at
least `2^{(1-h-o(1))A}`. *Heuristic.* There are `2^{hA + o(A)}` no-descent
classes modulo `2^{A+1}` (Terras; THM-4495 for the exact order), and if their
least representatives were spread over `[0, 2^{A+1})` the minimum would be
`2^{(1-h)A}`; the records lie above this line (`0.08–0.13` against `0.05`), so
the classes look *less* crowded near zero than at random. *Consequences.*
`g(A) → ∞`, hence `K ∩ N = ∅`; and `τ(n) ≤ (1 + o(1)) log_2 n / ((1-h) log_2 3)
= 12.6 log_2 n`, against the observed `5.6–7.0`. *Typing against the atlas.*
G is a statement about least representatives of residue classes (a size
statement), not about descent of a class (it uses Terras's memorylessness
only through the count `2^{hA}`), so it passes UNIFORM and STICKY; it is
sheet-sensitive through the sign (the classes of `-1, -5, -17` have least
positive representatives `2^{A+1} - 1, 2^{A+1} - 5, 2^{A+1} - 17`, the price
`A` of the sixth note, while the *minimum* over classes is what G bounds);
and it fails for `5x+1` (§5), as DRIFT demands. It is a quantified form of
the folklore "the glide is `O(log n)`", which the repo had not written down
as a statement about residue classes.

**A field to look at (pointer).** Lemma K1 makes the divergence half of
Collatz an *intrinsic Diophantine approximation* problem on a 2-adic Cantor
set: how well can the positive integers approach `K`? Approximation of and by
points of Cantor sets (Mahler's question on the middle-third set; results of
Bugeaud, Levesley–Salp–Velani and others, CITED from memory, unverified here)
is a developed subject in the real case; the 2-adic, one-sided (sign-aware)
version with the arithmetic constraint "`n` is a positive integer of size
`n`" is what Conjecture G asks for.

## 4. The mediant tree is the tree of all `3x + d` cycles (PROVED + FINITE-EXACT)

**Theorem F.** Let `x = m/d` in lowest terms with `d` odd and `3 ∤ d` (every
fixed point `x_w` has this form: denominators divide clocks, which are
`≡ ±1 mod 6`). The Syracuse map on `Z_(2)` sends `m/d` to `m'/d` with
`m' = oddpart(3m + d)`, and `gcd(m', d) = 1`. Hence the fixed points with
denominator exactly `d` are `1/d` times the cycle points of the map
`x ↦ oddpart(3x + d)` on positive integers coprime to `d`, negative `d`
giving `3x - |d|` (the minus sheet: `d = -1` gives `-1, -5, -17`). *Proof.*
`3(m/d) + 1 = (3m + d)/d`; `gcd(3m + d, d) = gcd(3m, d) = 1`. ∎ These are
Lagarias's rational cycles (1990, CITED from memory), and the weighted-mediant
law of the thirteenth note says the cycle point of a concatenated word is a
weighted mediant of cycle points of *different* `3x + d` problems:
`-17 = (α(-1) + β · 103/175)/(α + β)` with `α = -1539`, `β = 1400` joins the
`3x - 1` cycle `-1` (repeated thrice) to a cycle of `3x + 175`.

**Census (FINITE-EXACT).** Over all words with `A ≤ 22`, every `d ≤ 100`
coprime to 6 occurs as a denominator except `53, 65, 67, 79`, and a direct
orbit search finds cycles of `3x + 53`, `3x + 65`, `3x + 67`, `3x + 79`
(least elements `103`; `19`; `17`; `1, 7, 233, 265`), so every `d ≤ 100`
coprime to 6 is the denominator of a rational cycle. For `d = 5` the tree
points `m/5` with `A ≤ 22` have `m ∈ {1, 19, 23, 29, 31, 37, 49}`, exactly the
elements of the `3x+5` cycles through `1`, `19`, `23` (the cycles through
`187` and `347` need longer words); for `d = 13` the tree reaches the cycles
through `1, 211, 227, 251, 259, 283, 287, 319` and not yet `131`. Conjecture
(Lagarias's framework): every `d` coprime to 6 occurs, i.e. every `3x + d`
has a cycle on integers coprime to `d`; the eleventh note's prime-occurrence
lattice says which clocks `d` divides, and equidistribution of carries modulo
`D/d` (thirteenth note §3.4) is the heuristic behind existence.

## 5. The `5x+1` control (FINITE-EXACT)

For `5x + 1` (positive drift; the never-descending 2-adic set has Haar measure
`μ_5 = 0.176`, fourth note), among odd `n ≤ 3000`, `511` (a third) show no
word-level descent within 400 odd steps or before exceeding 600 bits,
starting with `7, 9, 21, 23, 25, 29`; the record function `g_5` stalls at `7`.
Conjecture G is therefore a `3x+1` statement: it encodes `log_2 3 < 2` through
`h < 1`, exactly as THM-4476's thin-divergence proof does.

## 6. Syntheses with the repo's own lenses

| repo object | what it is here | type |
|---|---|---|
| observer lens (marked basepoint; META-PATTERNS "keep observer, recurrence class, finite head"; "correct the object: shadow/marked observer") | Theorem M1: the orbit seen from `n` is the clock times the gap to the periodic point sharing `n`'s head; `K` is the observer's own set | EXACT |
| first note's product identity; THM-4476 two-place corollary | Theorem M1(3): the two-place product for the convergents `x_{w_j}` | EXACT |
| THM-4512 / precision-residual note (Gilbreath thread S11) | `N(w) = x_w`; residual = integers under their prefix's rational cycle point; empty to `10^6`; Diophantine finiteness heuristic | EXACT / PLAUSIBLE |
| sixth note's universal size price; THM-4478 provability price | `g(A)` is the price of depth `A`; Conjecture G is the size price stated for the *cheapest* class | EXACT reformulation |
| hedgehog identity `T^{2j}(1 + 4^j t) = 1 + 3^j t` (fourth note) | the prefix `(2)^j` has fixed point `1` and clock `4^j - 3^j`: M1(2) verbatim | EXACT |
| SHEET control | the sign of `d` in Theorem F; `-1, -5, -17` as `3x - 1` cycles; cross-sheet mediants | EXACT |
| eleventh note's prime-occurrence lattice; thirteenth note's carry equidistribution | which `d` divide clocks; why every `d` should occur | EXACT / heuristic |
| THM-972 (LRC relation lock: witnesses inherit light integer relations of the speeds) | the carry congruence `D_u S_{uv} ≡ 2^{A_u} β (mod D_{uv})`: carries inherit the clock relation `D_{uv} = 2^{A_u} D_v + 3^{p_v} D_u` modulo the clock | ANALOGY |
| THM-998 (LRC Farey-circle deep law: the deep set is a Farey dissection) | the tree's points organised by denominator `d` are the `3x+d` cycles; small `d` = "major arcs" of the periodic points | ANALOGY |
| barrier atlas (UNIFORM, STICKY, SHEET, DRIFT) | Conjecture G passes UNIFORM and STICKY, is sheet-aware, fails for `5x+1` | typed |

## 7. Verdicts

| claim | status |
|---|---|
| Theorem M1 (convergents, displacement, two-place product) | PROVED |
| Theorem M2 (threshold = fixed point; three sets) | PROVED (the threshold itself is THM-4512's) |
| residual empty for odd `3 ≤ n ≤ 10^6`; `σ = τ` there | FINITE-EXACT (THM-4512 to `10^7`) |
| Diophantine finiteness of the residual | PLAUSIBLE (Rhin's bound + Terras count) |
| Lemma K1 | PROVED |
| records to `2·10^6`; `κ = 12.2`; `g` values | FINITE-EXACT |
| Conjecture G and its two consequences | CONJECTURE (consequences PROVED from it) |
| Theorem F; denominators `≤ 100` all occur | PROVED + FINITE-EXACT |
| `5x+1` control | FINITE-EXACT |
| Collatz | OPEN |

## 8. Directions

* **D44.** Compute `g(A)` from the integer side to `10^9` (C, the record
  holders of `τ`) and from the class side to depth `A ≈ 40` (enumerate
  no-descent words, least representatives), and test Conjecture G's exponent;
  compare with Roosendaal's glide records (value-level), which coincide with
  `τ` records wherever the residual is empty.
* **D45.** Make the residual's finiteness effective: an explicit `A_0` beyond
  which `x_w < n_w` for every first positive-clock prefix `w`, from Rhin's
  exponent and the least-representative distribution; combined with the
  census this would prove `σ = τ` for all `n`, i.e. word-level and value-level
  stopping indices coincide.
* **D46.** Prove that every `d` coprime to 6 is a denominator (every `3x+d`
  has a cycle): choose a shape with `d | D` by the prime lattice and a word
  with `S ≡ 0 (mod D/d)` by equidistribution, then control the gcd.
* **D47.** Intrinsic approximation on `K`: the exponent of approximation of
  `n` by the periodic points `x_{w_j}` is 1 (M1(1) with height `2^{A_j}`);
  study the exponent of approximation of points of `K` by rationals with odd
  denominator — the analogue of Mahler's question — and whether `h` enters.
* **D48.** The observer identity for the parallel session's THM-4521 (base
  cycles cohere to every level): the periodic points `x_w` are its cycle
  closers; the level-`k` residues of `x_{w_j}` are those of `n` for `k ≤ A_j`.
