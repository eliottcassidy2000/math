# Posets and DAGs of a Collatz orbit: the value-time poset, the excursion forest and its spine, tight spine blocks, the two-place descent tree, and why every DAG finishing move is the transversality statement

**Session:** opus, `collatz-poset-dag-20260927`, 2026-09-27.
**Owner's directive:** "consider deeply posets and their directed acyclic
graphs as you push these final steps towards a creative Collatz proof
finishing move."
**Inherits (read, cited, not re-derived):**
[THM-4503](../../01-canon/theorems/THM-4503-collatz-poset-bridge-height-selection.md)
(no-descent words of a cell are the linear extensions of a width-2 poset;
height selection is not a poset operation),
[THM-4495](../../01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md)
(the Spitzer identity; its Step 1 is the ladder decomposition
`W = 1/(1 - P)` into positive min-ending words, its Step 2 the reversal),
[THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md)
(thin divergence; `sum 1/m_l < infinity`; no bounded strip, no narrow
log-band),
[THM-4506](../../01-canon/theorems/THM-4506-landing-multiplicity-exact-worst-case-and-recursion-saturation.md)
(dippers of a landing point lie in one dyadic shell),
[THM-4512](../../01-canon/theorems/THM-4512-collatz-coefficient-descent-classes-certified-thresholds.md)
(coefficient-descent classes certified above the threshold `N(w)`; Terras's
equality `kappa = sigma` for odd `n <= 10^7`),
the S13 note [`collatz_directions_20260926.md`](collatz_directions_20260926.md)
(the no-descent fractal `E_inf`, Proposition 4: no prefix rank),
the S8 note [`collatz_crossings_20260926_potential_and_seeds.md`](collatz_crossings_20260926_potential_and_seeds.md)
(section 2: the 1/3–2/3 construction; section 2.4: the orbit poset
speculation this note makes exact; `m_l = 2^(-d_l) S_(l-1) mod 3^l`),
the codex [excursion forest](forest_20260926_excursions.md) (dyadic-height
crossings with ordered carries; a different but cognate forest),
the [crossroads poset barrier audit](crossroads_poset_20260926_barrier.md)
(local-event floor; the Fibonacci boundary escape), the
[crossroads geometry note](crossroads_20260926_geometry.md) section 3 (a
real interval at each future minimum), the
[board note](entry_20260927_board.md) section 7 (obligation 2: "a
well-founded rank across completed excursions"; obligation 3: coverage),
the [barrier atlas](collatz_procgen_20260922_barrier_atlas.md)
(Yolcu–Aaronson–Heule: Collatz iff termination of a string rewriting
system; Kurtz–Simon), Terras 1976 and Lagarias 1985 (from memory).

**Status: PROVED elementary propositions (sections 2–5; the one genuinely new
exact object is Theorem 3, the tightness of spine blocks) + FINITE-EXACT
(descent tree to `2 * 10^6`, two-place in-degree formula to `10^5`, cells to
`a + b <= 12`, block enumeration to length 16 and generating functions to
length 60, controls) + DIRECTION (section 6, the finishing-move shapes, each
typed with its obstruction) + AUDIT pending (section 9). Collatz OPEN. Nothing
here proves anything about all `n`; the honest content is that every DAG or
poset finishing move reduces to the same transversality statement, and the
note says exactly where.** Script
`04-computation/experiments/collatz_posets_dags_20260927.py`, output
`collatz_posets_dags_20260927.out`. Reserved ID: THM-4514 (RESERVED stub until
audit).

## 0. The answer in one paragraph

A Collatz orbit carries three natural order structures, and all three say the
same thing. (1) The **value-time poset** (`i` below `j` iff `i` is earlier and
larger) has dimension at most two; its minimal elements are the leaders, its
maximal elements the strict future minima, and a positive-integer orbit has a
strict future minimum **iff it diverges** (Proposition 1). (2) The
**excursion forest** (`i` an ancestor of `j` iff the height walk stays above
`h_i` on `(i, j]`) has the strict lower records (the descent chain `n > D(n)
> D^2(n) > ...`) as its roots and, on a divergent orbit, exactly one infinite
tree whose spine is the sequence of strict coefficient future minima; the
spine is **the orbit's intersection with the no-descent fractal `E_inf`**, so
a positive integer lies in `E_inf` iff it is the minimum of a divergent orbit,
and then its orbit meets `E_inf` infinitely often along an increasing chain
(Proposition 2). (3) The spine cuts the word into **spine blocks** (positive
min-ending words, THM-4495's ladder blocks), and these are **tight**: a block
of length `l >= 2` starts with an odd step, has exactly `ceil(l log_3 2)` odd
steps, and has height in `(0, log_2 3 - 1)`; blocks exist only at the lengths
`l = 1 + floor(b log_2 3/(log_2 3 - 1))`, a shifted Beatty sequence of density
`1 - log_3 2 = 0.369` (Theorem 3). Hence every spine step multiplies the
value by less than `3/2` (plus a carry of at most `3l/4`): a divergent orbit
would climb its infinite ladder in steps of ratio `2^({a log_2 3})`, the
rotation by `log_2 3` seen at its smallest positive phases. (4) Reading the
descent tree `m -> D(m)` (the DAG of stopping times) with both places: `D(m)
in [m/2, m)`, so the lower records hit every dyadic shell (at least `floor(log_2
n) + 1` records, equality exactly at powers of two), and each first-descent
word is an affine bijection from a **2-adic** source class onto a **3-adic**
landing class `m_0(w) + t 3^(o(w))`; the in-degree of `m` is the number of
first-descent words whose 3-adic landing class contains `m` (exact to `10^5`,
modulo Terras's equality which is verified far beyond), with limiting mean
`c_D = sum_w 3^(-o(w)) in [1.6696, 1.7005]` (Proposition 4). What none of
this gives: a rank. The board's "well-founded rank across completed
excursions" is a function on the spine that must decrease, S13's Proposition 4
says it cannot be a function of the word prefix, and Theorem 3 sharpens that:
the prefix at a spine point is a concatenation of tight blocks, and the same
block sequence is realized by every integer in a residue class, so the rank
must read the integer, which only `sigma` itself does. The well-quasi-order
route (Kruskal/Higman, the only termination principle that needs no rank)
requires the spine of every divergent orbit to be a *bad* sequence, and on the
`5x+1` control (believed divergent) the spine has Higman-good pairs within at
most 141 blocks. Compactness gives nothing because `E_inf` is closed, nowhere
dense and of measure zero while `Z^+` is dense. The one place where the DAG
language is more than bookkeeping is HYP-9161's whole-orbit multiplicity:
the dippers of a landing point at depth `D` form a chain in the
`D`-coarsened excursion order, so the multiplicity is a chain length and the
whole-orbit average is a local-time statement about one orbit (section 6.6).

## 1. Inheritance pass and concept board

Closest proved mechanism: THM-4495's ladder decomposition (positive
min-ending words generate the no-descent words freely). Canonical hostile:
the residue class realizing any given block sequence (THM-4503's sidecar:
words are classes, not integers). Corrected near miss: S8's "orbit poset"
speculation, which took the relations from the descent data and found the
missing oscillation lemma; here the relations are the order itself, and the
poset is a permutation poset of dimension two with no dynamical content
beyond its extremal elements. Least-used sidecar: the 3-adic landing class of
a first-descent word (S8's `m_l = 2^(-d_l) S mod 3^l`, applied at the landing
time).

| Object | Representation | Predicate / invariant | Operation | What the quotient loses |
|---|---|---|---|---|
| Value-time poset `P_val(n)` | `i < j` and `x_i > x_j` (dimension `<= 2`) | maximal = strict future minima; minimal = leaders | intersection of two linear orders | the arithmetic step (any permutation is some poset) |
| Excursion forest `F(n)` | `i ⊑ j` iff `h_k > h_i` on `(i, j]` | roots = lower records = descent chain; infinite tree iff divergent; spine = strict coefficient future minima | König: infinite locally finite tree has a branch | carries (the forest is a function of the word) |
| Spine blocks | positive min-ending words (THM-4495 Step 1) | tight: `ceil(l log_3 2)` ones, height `< log_2 3 - 1`, lengths `1 + floor(b beta)` | free monoid; reversal to first-passage words | which integer realizes the block sequence |
| Descent tree `D` | `m -> D(m) = T^(sigma(m)) m in [m/2, m)` | spanning tree rooted at 1 iff Collatz; two-place edges | 2-adic source class -> 3-adic landing class (affine bijection per word) | nothing on the class; the threshold `N(w)` at the least residue |
| Cell lattice `L(a, b)` (THM-4503) | interval `[empty, g]` of Young's lattice | carry `C_w` strictly increasing; least residue not monotone | add a cell = move an odd step one place later | height selection (THM-4503) |
| No-descent tree / `E_inf` | prefix tree of no-descent cylinders | infinite paths = `E_inf`; closed, nowhere dense, dimension `h*` | König gives `E_inf != empty`; transversality with `Z^+` | everything: the conjecture is the transversality |

## 2. The value-time poset and the excursion forest (PROVED)

Notation. `T(x) = x/2` (`x` even), `(3x+1)/2` (`x` odd); parity word
`w_1 w_2 ...`; `o_j` = number of ones among the first `j` letters; height
`h_j = o_j log_2 3 - j`, so `T^j(n) = n 2^(h_j) C_j` with the carry factor
`C_j = prod_(l < j, x_l odd) (1 + 1/(3 x_l)) >= 1` non-decreasing in `j`
(equivalently `T^j n = (3^(o_j) n + C_(w, j))/2^j`). Heights at distinct
times are distinct and `h_j != 0` for `j >= 1` (`log_2 3` irrational). `E_inf
= {x in Z_2 : h_j(x) >= 0 for all j >= 1}` (S13). All comparisons of heights
below are exact (`3^a` against `2^b`).

**Definition.** For an injective orbit segment `x_0, x_1, ..., x_K` (or the
whole orbit when injective) the *value-time poset* is `i ≺ j` iff `i < j` and
`x_i > x_j`. It is the intersection of the time order with the reversed value
order, hence a poset of dimension at most two (a permutation poset). The
*excursion order* on the same times is `i ⊑ j` iff `i <= j` and `h_k > h_i`
for every `k in (i, j]`.

**Proposition 1.** (a) The minimal elements of `P_val` are the strict upper
records (leaders, `x_i > x_l` for all `l < i`); the maximal elements are the
strict future minima within the segment (`x_j < x_k` for all `k > j`).
(b) `⊑` is a forest order: transitive, and the ancestors of any `j` form a
chain. Its roots are the strict lower records of the height walk. A node `j`
has infinitely many descendants iff `j` is a strict coefficient future minimum
(`h_k > h_j` for all `k > j`). If `h_j -> infinity` the forest is locally
finite (every node has finitely many children), so by König's lemma it has an
infinite tree iff it has an infinite branch, and every infinite branch is a
tail of the sequence of strict coefficient future minima (the *spine*).
(c) For the orbit of a positive integer `n`: a strict value future minimum
exists iff the orbit diverges, and then there are infinitely many; strict
coefficient future minima are strict value future minima, and strict value
lower records are strict height lower records (positive carry); the step out
of any value future minimum is an odd step (the value rises).

*Proof.* (a) is the definition. (b) Transitivity: if `h > h_i` on `(i, j]`
and `h > h_j` on `(j, k]` with `h_j > h_i`, then `h > h_i` on `(i, k]`. Chain:
if `i ⊑ j` and `i' ⊑ j` with `i < i'` then `i' in (i, j]` so `h_(i') > h_i`,
and `h > h_(i') > h_i` on `(i', j]`, whence `i ⊑ i'`. Roots: `j` is a root iff
no `i < j` has `h > h_i` on `(i, j]`; if `h_j < h_i` for all `i < j` this
holds (take `k = j`); otherwise the minimum of `h` over `[0, j]` is attained
at some `i* < j` and `h > h_(i*)` on `(i*, j]`, so `i*` is an ancestor.
Descendants of `j` are the `k >= j` with `h > h_j` on `(j, k]`, an infinite set
iff `h_k > h_j` for all `k > j`. Children of `j` are the running minima of the
walk on `(j, infinity)` restricted to heights above `h_j`; if `h -> infinity`
these are finitely many (the running minimum stops changing once the walk
exceeds its own minimum on the tail, attained at the next spine point). The
branch through a spine point `f` continues through the next spine point (the
last running minimum of the tail from `f` is the tail's minimum, a spine
point), and any infinite branch consists of nodes with infinitely many
descendants, i.e. spine points. (c) If the orbit is eventually periodic its
global minimum is attained infinitely often, so no time is a strict value
future minimum. If it diverges, every value is undercut only finitely often;
the minimum of each tail is attained, at a strict future minimum, and the
tails' minima tend to infinity, so there are infinitely many. Inclusions:
`x_k/x_j = 2^(h_k - h_j) C_k/C_j` with `C_k >= C_j` for `k > j`, so
`h_k > h_j` implies `x_k > x_j`, and `x_j < x_i` with `i < j` implies
`2^(h_j) < 2^(h_i) C_i/C_j <= 2^(h_i)`. The step out of a value future
minimum must increase the value, and only odd steps do. ∎

*Remark (Terras along the chain).* The strict value lower records of the
orbit of `n` are the descent chain `n > D(n) > D^2(n) > ...` and the strict
height lower records are the coefficient-descent chain; by induction they
coincide exactly when Terras's equality `kappa = sigma` holds at every element
of the chain. Since every chain element is at most `n`, the two chains
coincide for all `n <= 10^7` (THM-4512); the script confirms set equality for
all `n <= 10^4` and on the eight record orbits of section 8.

**Proposition 2 (the spine is the orbit's intersection with `E_inf`).** For a
positive integer `n` with `T`-orbit `(x_j)`, `x_j in E_inf` iff `j` is a strict
coefficient future minimum. Consequently: (i) `n in E_inf` implies that the
orbit of `n` is not eventually periodic, hence `x_j -> infinity`, hence (by
THM-4476) `sum_j 1/x_j < infinity`, `C_infinity < infinity`, `h_j ->
infinity` and `sum_j 2^(-h_j) < infinity`, and the spine of `n` is infinite;
(ii) `Z^+ ∩ E_inf = {m : kappa(m) = infinity}` is contained in the set of
positive integers that are the minimum of a divergent orbit (`sigma(m) =
infinity` and `m` not a cycle minimum), with equality iff Terras's equality
`kappa(m) = sigma(m)` holds at every minimum `m` of a divergent orbit; it is
closed under the *next-spine-point* map `s(m) = x_(f')` (the next strict
coefficient future minimum of the orbit of `m`), and `s(m) > m`; (iii) an
orbit that reaches `1`, or any eventually periodic positive orbit, meets
`E_inf` nowhere.

*Proof.* `x_j in E_inf` iff the shifted word has all heights `>= 0`, i.e.
`h_(j+k) - h_j >= 0` for all `k >= 1`, i.e. (distinct heights) `h_(j+k) > h_j`
for all `k >= 1`. (i) If the orbit of `n in E_inf` were eventually periodic
with period word `u`, then `H(u) != 0`; `H(u) < 0` sends the heights to
`-infinity`, contradicting `n in E_inf`; `H(u) > 0` is impossible for a
positive cycle because `T^p(m) = m` on the cycle gives `(3^(o(u)) - 2^(|u|))
m + C_u = 0` with `C_u > 0`, i.e. `3^(o(u)) < 2^(|u|)`. So the orbit is
injective; an injective sequence of positive integers tends to infinity;
THM-4476 gives `sum 1/x_j < infinity`, so `C_j` increases to a finite
`C_infinity`, and `h_j = log_2(x_j/(n C_j)) >= log_2(x_j/(n C_infinity)) ->
infinity`, with `sum 2^(-h_j) = sum n C_j/x_j <= n C_infinity sum 1/x_j <
infinity`. A walk tending to `+infinity` has infinitely many strict future
minima. (ii) `m in E_inf` iff `kappa(m) = infinity` by definition of `E_inf`.
Since `kappa <= sigma` (an actual descent `T^k m < m` forces `3^(o_k) m +
C_(w,k) < 2^k m`, hence `3^(o_k) < 2^k`), `kappa(m) = infinity` gives
`sigma(m) = infinity`: `m` is the minimum of its orbit, which by (i) is not a
cycle, so `m` is the minimum of a divergent orbit. Conversely a minimum `m` of
a divergent orbit has `sigma(m) = infinity`, and it lies in `E_inf` iff
`kappa(m) = infinity`, i.e. iff `kappa(m) = sigma(m)`; a coefficient descent
`3^(o_k) < 2^k` at a time `k` with `T^k m > m` is possible in principle (it
needs `(2^k - 3^(o_k)) m < C_(w,k)`, i.e. `m` below THM-4512's threshold), and
Terras's equality is exactly the statement that this does not happen. Closure
under `s`: `s(m)` is a spine point of the same orbit, hence in `E_inf`, and
`s(m) > m` because spine points are value future minima (Proposition 1(c))
and `s(m)` comes later. (iii) Eventually periodic positive orbits have no
strict coefficient future minima: the period height is negative by the cycle
identity, so the heights tend to `-infinity` and every time is undercut
later. ∎

*Scope note.* The unconditional content of (ii) is the inclusion `Z^+ ∩ E_inf
⊆ {minima of divergent orbits}` together with closure under `s`; equality
needs Terras's equality only at orbit minima of divergent orbits. This is
the same quantifier that separates THM-4512's certified classes from the
classical residual (S11 correction), and the note keeps it visible.

## 3. Spine blocks are tight (PROVED; the new exact object)

**Definition.** A *spine block* is a finite word `w` of length `l >= 1` with
`H(w) := h_l(w) > 0` and `h_j(w) > H(w)` for `1 <= j < l`: a positive word
ending at its strict minimum over positive times. These are THM-4495's
positive min-ending words `P_l` (its Step 1, `W(t) = 1/(1 - P(t))`).

**Theorem 3.** (i) *Decomposition.* Every finite no-descent word is uniquely a
concatenation of spine blocks (cut at the strict future minima within the
word), and every concatenation of spine blocks is a no-descent word with
that decomposition; every infinite word with `h_j -> +infinity` and all `h_j
> 0` is uniquely an infinite concatenation of spine blocks, cut at its spine.
Hence `sum_l b_l t^l = 1 - 1/W(t) = 1 - exp(-sum_(n >= 1) B_n t^n/n)` with
`B_n = sum_(3^j > 2^n) C(n, j)` (THM-4495 (A)).
(ii) *Tightness.* A spine block of length `l >= 2` begins with the letter `1`
and satisfies `0 < H(w) < log_2 3 - 1`. Consequently its number of ones is
forced, `a(l) = ceil(l log_3 2)` (the least `a` with `3^a > 2^l`), its height
is `H_l = a(l) log_2 3 - l = (1 - {l log_3 2}) log_2 3`, and spine blocks of
length `l >= 2` exist iff `{l log_3 2} > log_3 2`, iff `l - 1` is not of the
form `floor(a log_2 3)` (`a >= 1`), iff `l = 1 + floor(b beta)` for some `b >=
1`, where `beta = log_2 3/(log_2 3 - 1) = 2.70951...`. The block lengths have
density `1 - log_3 2 = 0.36907`. The word `1^(a(l)) 0^(l - a(l))` is a block
at every such length.
(iii) *Spine growth.* Along the spine `s_1 < s_2 < ...` of a divergent orbit,
with blocks `w_i` of lengths `l_i`,

```text
2^(H(w_i)) s_i < s_(i+1) < 2^(H(w_i)) s_i (1 + l_i/(2 s_i)),   0 < H(w_i) < log_2 3 - 1  (l_i >= 2),
```

so `s_(i+1) < (3/2) s_i + (3/4) l_i`; a one-letter block gives exactly
`s_(i+1) = (3 s_i + 1)/2`. The spine heights `h_(f_i)` increase by
`H(w_i) in (0, log_2 3 - 1]` per block, i.e. by the fractional parts `{a
log_2 3}` that fall below `log_2 3 - 1`.

*Proof.* (i) Let `w` be finite with all `h_j > 0` for `1 <= j <= k`. Its strict
future minima within `[0, k]` are `0 = f_0 < f_1 < ... < f_r = k` (`0` because
all later heights are positive, `k` vacuously). For `f_i < j < f_(i+1)`, `j` is
not a future minimum, so some later `k' <= k` has `h_(k') < h_j`; the minimum
of `h` over `[j, k]` is attained at a future minimum `>= f_(i+1)`, whose height
is `>= h_(f_(i+1))`; so `h_j > h_(f_(i+1)) > h_(f_i)`. Read the segment
`(f_i, f_(i+1)]` as a word: its prefix heights are `h_j - h_(f_i) > h_(f_(i+1))
- h_(f_i) > 0`, a spine block. Conversely a concatenation of blocks has, at the
end of each block, a strict new minimum over positive times (each block's
interior stays above its end and starts above the previous end), so all its
heights are positive and its future minima are exactly the block ends. The
infinite case is the same with the spine in place of the window minima (a walk
tending to `+infinity` has infinitely many strict future minima, and each
block between consecutive ones is finite). The generating function is
THM-4495's Step 1 with its identity (A). (ii) If the first letter were `0`,
`h_1 = -1 < 0 < H`, contradicting `h_1 > H`. So `h_1 = log_2 3 - 1 > H(w)`.
The numbers `a log_2 3 - l` for consecutive `a` differ by `log_2 3 > log_2 3 -
1`, so at most one `a` has `a log_2 3 - l in (0, log_2 3 - 1)`, and it is the
least `a` with `3^a > 2^l`, i.e. `a(l) = ceil(l log_3 2)`; the condition `a(l)
log_2 3 - l < log_2 3 - 1` reads `(a(l) - 1) log_2 3 < l - 1`. Existence: `1^a
0^(l-a)` with `a = a(l)` has heights `j(log_2 3 - 1) >= log_2 3 - 1 > H_l` for
`j <= a` and `a log_2 3 - j`, decreasing, `> H_l` for `a <= j < l`, ending at
`H_l > 0`. Beatty form: `l - 1 = floor(a log_2 3)` for some `a >= 1` iff the
interval `[(l-1) log_3 2, l log_3 2)` (length `log_3 2 < 1`) contains an
integer iff `floor(l log_3 2) >= (l - 1) log_3 2` iff `{l log_3 2} < log_3 2`
(equality impossible). Rayleigh's theorem gives the complement of `{floor(a
log_2 3)}` in `Z^+` as `{floor(b beta)}` with `1/log_2 3 + 1/beta = 1`. The
density of `{l : {l log_3 2} > log_3 2}` is `1 - log_3 2` by equidistribution.
(iii) `s_(i+1) = (3^a s_i + C_w)/2^l = 2^(H) s_i (1 + C_w/(3^a s_i))` with
`C_w = sum_i 2^(j_i) 3^(a - i)` over the positions `j_i` of the ones; the
height after the `i`-th one is `i log_2 3 - (j_i + 1) > H > 0`, so `2^(j_i)/3^i
< 1/2` and `C_w/3^a < a/2 <= l/2`. The factor `2^H` is below `3/2` by (ii). ∎

*Numerical face (script, part 2).* `b_l` for `l <= 24`: `1, 0, 1, 0, 0, 2, 0,
0, 7, 0, 30, 0, 0, 113, 0, 0, 525, 0, 2652, 0, 0, 11433, 0, 0`; the two
generating-function forms agree exactly to `l = 60`; the block lengths up to
`60` are `1, 3, 6, 9, 11, 14, 17, 19, 22, 25, 28, 30, 33, 36, 38, 41, 44, 47,
49, 52, 55, 57, 60 = {1} ∪ {1 + floor(b beta)}`; on block lengths `b_l/W_l in
[0.0705, 0.1124]` for `40 <= l <= 60`; every block of length `<= 16` has
exactly `a(l)` ones and starts with `1`.

*What (iii) says about a hypothetical divergent orbit.* Its spine is an
increasing ladder of elements of `E_inf` climbing by factors `2^({a log_2 3})`
restricted to the phases `{a log_2 3} < log_2 3 - 1`: the rotation by `log_2 3`
(S8's phase rotation, the exact Benford law) is visible on the spine at its
smallest positive phases only. This is a structural constraint on the *word*
of a divergent orbit, not on the integer; the same block sequence is the word
of every integer in the corresponding residue class (THM-4503's sidecar), so
it cannot by itself exclude integers.

## 4. The descent tree with both places (PROVED + FINITE-EXACT)

For `m >= 2` with `sigma(m) < infinity` let `D(m) = T^(sigma(m))(m)` be the
first orbit value below `m`. The *descent tree* has the edges `m -> D(m)`;
Collatz holds iff it is a spanning tree of `Z^+` rooted at `1` (every `m >= 2`
has `sigma(m) < infinity`; the chain `m > D(m) > ...` then reaches `1`).

**Proposition 4.** (a) `D(m) in [ceil(m/2), m - 1]`: the step into the landing
is a halving of a value `>= m`.
(b) Every dyadic shell `[2^j, 2^(j+1))`, `0 <= j <= floor(log_2 m)`, contains
a strict lower record of the orbit of `m`; hence the number of strict lower
records (the length of the descent chain to `1`) is at least `floor(log_2 m) +
1`, with equality iff `m` is a power of two.
(c) *Two-place structure.* Let `w` be a word of length `k` with `o` ones whose
first coefficient descent is at `k` (`3^(o_j) > 2^j` for `j < k`, `3^o < 2^k`),
carry `C = C_w`, source residue `r = -C 3^(-o) mod 2^k`, threshold `N(w) =
C/(2^k - 3^o)`, and `m_0 = (3^o r + C)/2^k` (an integer). The affine map `m'
-> (3^o m' + C)/2^k` restricts to a bijection

```text
{ r + t 2^k : t >= 0, r + t 2^k > N(w) }  ->  { m_0 + t 3^o : same t },
```

and these are exactly the descent-tree edges out of sources whose actual
first descent coincides with their coefficient first descent (`sigma =
kappa`). Hence, for every landing `m`,

```text
indeg_D(m) = #{ w : m ≡ m_0(w) (mod 3^(o(w))), the source (2^k m - C_w)/3^(o(w)) exceeds N(w) }
           + #{ m' in (m, 2m] : sigma(m') > kappa(m'), D(m') = m },
```

where the second term is the Terras-defect count, empty for all landings `m
<= 10^6` (Terras's equality is verified here to `2 * 10^6` and in THM-4512 to
`10^7`).
(d) *Sources per landing.* `#{m' : D(m') <= X}/X -> c_D := sum_w 3^(-o(w)) =
sum_w 2^(-|w|) 2^(|H(w)|)` as `X -> infinity` (sum over all first-descent
words; `|H(w)| = k - o log_2 3 in (0, 1)` is the overshoot at the first
descent), and `1 < c_D < 2`; the partial sum over `|w| <= 26` is `1.669582`,
the tail is at most `2(1 - 0.984542) = 0.0309`, so `c_D in [1.6696, 1.7005]`.

*Proof.* (a) `T^(sigma) m < m <= T^(sigma - 1) m`; an odd step increases the
value, so the last step is a halving, `D(m) = T^(sigma - 1)(m)/2 >= m/2`.
(b) The chain decreases with consecutive ratios in `[1/2, 1)`: when it passes
from a value `>= 2^(j+1)` to a value `< 2^(j+1)` the landing is `>= 2^j`.
Equality forces exactly one record per shell. Induction upward: the record
in shell `0` is `1`; the record in shell `1` is `2` or `3`, and `D(3) = 2`
would put a second record in shell `1`, so it is `2`; if the record in shell
`j` is `2^j` and `y` is the record in shell `j + 1`, then `D(y)` is the next
record, which lies in shell `j` (shell `j + 1` has no other record), so `D(y)
= 2^j`, and `D(y) >= y/2` gives `y <= 2^(j+1)`, hence `y = 2^(j+1)`. So `m =
2^J`. Conversely powers of two halve straight down, one record per shell.
(c) Sources in the class `r mod 2^k` have `w` as their length-`k`
parity prefix (Terras's bijection), the coefficient descent at `k` is an actual
descent iff `(2^k - 3^o) m' > C`, i.e. `m' > N(w)`, and then `sigma(m') = k =
kappa(m')` (no earlier coefficient descent, and an earlier actual descent would
be an earlier coefficient descent since the carry is positive); the image of
`r + t 2^k` is `m_0 + t 3^o`. Every descent-tree edge from a source with
`sigma = kappa` arises this way from its own first-descent word. (d) Summing
(c) over words of length `<= K`: `#{m' : D(m') <= X} = sum_(|w| <= K) (X
3^(-o(w)) + O(1)) + R_K(X)`, where `R_K(X)` counts sources `m' <= 2X` with
`sigma(m') > K`; these lie in no-descent classes mod `2^K` (`kappa > K`) or
are Terras-defect sources of words of length `<= K` (finitely many per word,
below `N(w)`), so `R_K(X) <= 2X W_K 2^(-K) + O_K(1)`, and `W_K 2^(-K) -> 0`.
The series converges because `3^(-o) < 2^(1 - k)` for a first-descent word
(`3^o > 2^(k-1)` from the no-descent prefix), and the same inequality with
`3^o < 2^k` gives `2^(-k) < 3^(-o) < 2^(1-k)`, whence `1 = sum_w 2^(-|w|) <
c_D < 2` (the first equality is the measure-zero of `E_inf`, `W_k 2^(-k) ->
0`). ∎

*Numerical face (script, part 3, `T`-map, `2 <= m <= 2 * 10^6`).* `D(m) in
[m/2, m)` for all `m`; maximal `sigma = 224` at `m = 1126015`; no Terras
defect; the record-count bound holds for all `m` with equality exactly at the
twenty powers of two; in-degrees over landings `1..10^6`: `1: 531332, 2:
337608, 3: 78502, 4: 28890, 5: 15341, 6: 4769, 7: 2102, 8: 874, 9: 348, 10:
129, 11: 61, 12: 26, 13: 13, 15: 3, 16: 1, 17: 1`, mean `1.6903`, maximum `17`
at `m = 293501`; `190069` first-descent words of length `<= 26` (`F_k = 2
W_(k-1) - W_k` checked); the 3-adic prediction matches the direct in-degree for
every landing `m <= 10^5` (sources with `sigma > 26` landing there: `2071`,
reconciled separately); example `w = 11100`: `o = 3`, `C = 19`, `r = 23 mod
32`, `N = 3.8`, landings `20, 47, 74, 101, 128, 155 = 20 + 27t`.

*Reading.* The descent tree is a union of affine bijections between 2-adic
source classes and 3-adic landing classes, one per first-descent word: the
Terras structure (sources are classes mod `2^k`) has an exact dual on the
landing side (landings are classes mod `3^o`). This is S8's "the past is the
3-adic address" (`m_l = 2^(-d_l) S mod 3^l`) evaluated at the landing time.
It makes obligation 3 of the board ("coverage") literal: Collatz iff every `m
>= 2` lies in some source class above its threshold, i.e. `kappa(m) <
infinity` (`m not in E_inf`, the transversality) and `m > N(w_m)` (the
uncertified-representative question of THM-4512, S11).

## 5. THM-4503's cells are intervals of Young's lattice (PROVED, small)

Write a word of the cell `(a, b)` (a ones, b zeros, every nonempty prefix
with `i` ones and `j` zeros having `3^i > 2^(i+j)`) by the partition `e_1 <=
... <= e_a`, `e_i` = number of zeros before the `i`-th one. THM-4503's
condition "the `j`-th zero is preceded by at least `f(j) = m_j` ones" is
`e_(f(j)) <= j - 1`; with `g_i = min({j - 1 : f(j) >= i} ∪ {b})` the cell is
empty if `f(b) > a` and otherwise equals the set of partitions `e <= g`
componentwise: the interval `[empty, g]` of Young's lattice inside the `a x
b` box, whose covering relation "add one cell" moves one odd step one place
later past a halving (THM-4503's swap `10 -> 01`). The carry `C_w = sum_i
3^(a-i) 2^(i-1) 2^(e_i)` increases by `3^(a-i) 2^(i-1) 2^(e_i) > 0` under
that cover: **`C_w` is a strict order embedding of the cell lattice into
`(Z, <)`.** The least residue `r_w = -C_w 3^(-a) mod 2^(a+b)` is not
monotone: for `(5, 1)` the chain `111110 ⊂ 111101 ⊂ 111011 ⊂ 110111` has
residues `31, 47, 39, 27` (THM-4503's numbers), and THM-4503's height
selection `{27, 31}` is the top and bottom of the chain, which is neither an
ideal nor a filter. Verified for all cells with `a + b <= 12` (870 covering
pairs). The pointwise dominance order on all words of length `k` (`o_j <=
o'_j` for all `j`) is Young's lattice in the staircase, the no-descent words
are the principal filter above the upper mechanical word of slope `log_3 2`
(section 7), and `C_w` is strictly *decreasing* along dominance (odd steps
earlier, carry smaller). Nothing here is height selection; it is recorded so
that the lattice is not re-derived.

## 6. The finishing-move shapes, typed (DIRECTION)

Every DAG or poset termination argument is one of the following. For each:
source, target, map, preserved predicate, destroyed information, needed
sidecar, cheapest decisive test.

**6.1 Well-foundedness = rank (the excursion forest).** A divergent orbit is
an infinite tree of the excursion forest with spine `s_1 < s_2 < ...`; a
proof by well-foundedness is a function `Phi` on spine points, valued in a
well-order, with `Phi(s_(i+1)) < Phi(s_i)` (the board's "rank across
completed excursions"; the completed excursions are the blocks). Source:
the word; target: the integer. Preserved: the block sequence. Destroyed: the
integer (every integer in the class realizes the same blocks). S13's
Proposition 4 forbids `Phi` to depend on the finite prefix; Theorem 3 makes
the witness concrete: the prefix at `s_i` is a concatenation of tight blocks
(the S13 word `(1, 1, 2)^N` in valuation coding is the block sequence
`(1)(110)` repeated), and the residue class of that concatenation contains
integers with every possible future. Sidecar: the whole integer. The only
known such `Phi` is `sigma`. Cheapest test: none needed; this is a theorem.

**6.2 Well-quasi-orders (Kruskal, Higman; Dershowitz's simplification
orders; Yolcu–Aaronson–Heule's rewriting system).** The one termination
principle that needs no rank: if a relation is contained in a well-founded
*simplification order* (containing the embedding of a wqo), it terminates,
because an infinite run would be a bad sequence. Source: the spine (or the
block sequence, or the excursion trees); target: a wqo `<=` on `Z^+` with
`s_i <= s_j` for no `i < j` on any divergent spine. Preserved: nothing yet;
this is exactly what must be proved. Destroyed: the sheet and the drift (any
wqo statement about words alone is sign- and drift-blind). Obstruction, as
a computation: on the `5x+1` orbit of `7` (believed divergent, the DRIFT
control) the first 300 spine points contain Higman-good pairs (binary
expansion of `s_i` a subsequence of that of `s_j`) for 165 of them within the
window, with waits `j - i` up to `141`, mean `61`, never `1`. A wqo proof for
`3x+1` would therefore have to show that *its* hypothetical spines avoid, for
ever, what the `5x+1` spine does within 141 blocks; that is a sheet- and
drift-sensitive statement about integers, which is the conjecture again.
Sidecar: an embedding order that sees `3` against `2` (none known). Cheapest
test done.

**6.3 Balanced pairs (1/3–2/3; THM-4503; the crossroads barrier).** Source:
linear extensions of the cell poset (section 5: an interval of Young's
lattice); target: a descent statement for one integer. Preserved: the word
law. Destroyed: the height selection (THM-4503), the source identity, the
conditioning (barrier note section 4). The value-time poset of an orbit is a
permutation poset; its balanced pairs are pairs of times whose *order* is
undetermined by the descent data, which is S8's oscillation lemma (HYP-9161
in crossing form), not a new handle. Cheapest test: THM-4503's `(5, 1)`
example already refutes a transfer.

**6.4 König / compactness on the no-descent tree.** The prefix tree of
no-descent cylinders is infinite and binary, so `E_inf != empty` (König);
`E_inf` is closed, nowhere dense (every cylinder contains a descending word),
of Haar measure zero (`W_k 2^(-k) -> 0`) and box dimension `h*`; `Z^+` is
dense in `Z_2`. A dense set and a closed nowhere-dense null set are disjoint
all the time; no compactness argument separates them. Preserved: the
cylinder structure. Destroyed: finiteness of the binary expansion, which is
the only property distinguishing integers (S13 section 2). Cheapest test:
the Sturmian element (section 7) is in `E_inf` and is not an integer for a
*different* reason (divergent real series), so `E_inf` has at least two kinds
of non-integer points, and the integer points, if any, are the spines.

**6.5 Two-place self-consistency along the spine.** At a spine time `f` with
`d_f >= log_2(n + 1) + 1`, the value is the least positive residue `x_f =
(2^(-d_f) S_(f-1) mod 3^f)` (S8), a function of the first `f` letters alone;
and the tail word from `f` is the parity word of `x_f`. So the word of a
divergent orbit is a fixed point: from some time on, each tail is the word
of the least residue of the 3-adic address of its past. "No infinite block
sequence is the word of its own 3-adic least residue" is an exact
reformulation and not a method; Proposition 4(c) is its finite shadow (the
landing of a class is a 3-adic class). Cheapest test: none; it restates the
conjecture.

**6.6 Where the DAG language is more than bookkeeping: HYP-9161.** Define the
`D`-coarsened excursion order `i ⊑_D j` iff `i <= j` and `h_k > h_i - D` on
`(i, j]`. A dipper `i` of a landing point `j` (THM-4506) is exactly a node
whose `⊑_D`-subtree ends at `j - 1`; two dippers `i < i'` of the same `j`
satisfy `i ⊑_D i'`. So **the dippers of a landing point form a chain of the
`D`-coarsened forest, the multiplicity `m(j)` is a chain length**, THM-4506's
shell lemma says the chain lies in a height window of width one, and its
odd-letter separation is the statement that consecutive chain elements are
separated by a rise. HYP-9161's whole-orbit average multiplicity is the mean
length of these chains along one orbit, a local-time statement (how long the
walk of one orbit hovers in unit windows before dropping `D`). The forest
gives the exact bookkeeping and the correct object; it does not give the
bound, because the chain lengths are realized by residue classes at every
scale (THM-4506 (H)). Cheapest next probe: compute the chain-length
distribution of the `D`-coarsened forest on long real orbits and on the
`5x+1` control and compare with the ballot prediction `sqrt(L)` of S9.

**Ranked next steps.** (1) 6.6: the `D`-coarsened chain-length statistics on
one orbit (the only orbit-coupled object here). (2) Theorem 3 (iii) on the
`5x+1` control: the spine ratios `2^({a log_2 5})` restricted to phases below
`log_2 5 - 1` (verified on the window: forced ones count and height bound
hold); ask whether the phase sequence along a real divergent spine is
equidistributed on `(0, log_2 q - 1)` in the rotation measure, which would be
the spine form of the exact Benford law of S8. (3) The Terras-defect term of
Proposition 4(c): a proof that `N(w) < 2^(|w|)` for every first-descent word
(THM-4512 has it for `|w| <= 5000` and S11's correction records that the
all-`j` claim needs an effective cutoff) would make the in-degree formula
unconditional for every `m`. None of these is a finishing move.

## 7. Controls (FINITE-EXACT)

**The Sturmian element of `E_inf`.** The upper mechanical word of slope `log_3
2` (`o_j = ceil(j log_3 2)`) has heights `h_j = (1 - {j log_3 2}) log_2 3 in
(0, log_2 3)` (checked to `j = 20000`: minimum `0.000063`, maximum
`1.584621`). It is in `E_inf` and has **no** strict future minimum (the
fractional parts are dense), so its spine is empty and Proposition 1(b)'s
infinite tree does not exist for it; within a window its future minima are
the one-sided best approximations of `log_3 2`, with block lengths `1, 19,
84, 569, 1054` up to `20000` (denominators of the convergents `12/19, 53/84,
665/1054` and the semiconvergent `359/569`); its real series diverges
linearly (every term `2^(-h_j) >= 1/3`), so by THM-4476 it is not a positive
integer, and by Proposition 2 an integer point of `E_inf` has an infinite
spine and a convergent series: the bounded-height part and the integer part
of `E_inf` are disjoint for two independent reasons. Its cylinders' least
residues `r_k` sit at `0.98, 0.36, 0.74, 0.79, 0.67, 0.25, 0.005, 0.68` times
`2^k` for `k = 8, 16, ..., 64`: nothing special about them.

**The `5x+1` orbit of `7` (DRIFT control).** In `3000` steps: `512` strict
coefficient future minima, equal to the strict value future minima within
the window (inclusion is a theorem, equality here an observation); block
lengths mean `5.87`, maximum `174`; every block of length `>= 2` starts with
an odd step, has exactly `ceil(l log_5 2)` odd steps and height below `log_2
5 - 1` (Theorem 3 (ii) holds verbatim for `qx+1`: the window `(0, log_2 q -
1)` has length below the spacing `log_2 q`); the value at step `3000` has
`451` bits. Higman-good pairs on the spine: section 6.2.

**Record orbits (`T`-map).** For `n = 27, 703, 6171, 77031, 837799, 8400511,
63728127, 670617279`: strict lower records `8, 14, 20, 27, 33, 35, 38, 44`
against the bound `floor(log_2 n) + 1 = 5, 10, 13, 17, 20, 24, 26, 30`;
leaders `16, 17, 16, 15, 23, 27, 25, 19`; forest roots equal the lower
records in all eight cases (Terras along the chain); excursion-forest depth
`17, 15, 15, 14, 23, 25, 26, 25`; within the window the only strict value
future minimum is the final `1` (Proposition 1(c)).

## 8. Reproduction

```text
python 04-computation/experiments/collatz_posets_dags_20260927.py
```

Parts: P1 value-time poset and forest on the eight record orbits and set
equality of value/height lower records for `n <= 10^4`; P2 block enumeration
to length 16, both generating-function forms to length 60, the block-length
law and the forced ones count; P3 the descent tree to `2 * 10^6`, the
record-count bound and its equality set, in-degrees, the two-place formula to
`10^5`; P4 the cells to `a + b <= 12`; P5 the Sturmian point; P6 the `5x+1`
control with the Higman probe. Every assertion is exact integer arithmetic;
decimals are display only. Output: `collatz_posets_dags_20260927.out`.

## 9. Independent audit

Pending (a blind re-derivation by an auditor subagent with its own script is
run before the letter; corrections are applied above and listed here).
