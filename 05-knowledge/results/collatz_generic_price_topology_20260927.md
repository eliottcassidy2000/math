# D14 and the topological reframe: the empirical price of generic growth, the universal size price K − z_K(x) and the least-representative tautology, what the necklace-splitting mechanism would need and why the sheet symmetry forbids it, and the fixed points of the conjugacy map

**Session:** opus, `collatz-posets-zeta5-20260927` (sixth note), 2026-09-27.
**Owner's directive:** "keep going, pursue D14 pricing the generic words,
consider the stolen necklace problem and the techniques that are used to
solve it, and think in terms of similar clever topological analogous
insights, and how Collatz parts can be rethought in terms of Borsuk–Ulam and
fixed points."
**Inherits (cited):** the fifth note (price sheet; D14 as stated there, now
corrected), THM-4512 (exact classes), Terras 1976, Lagarias 1985 and
Bernstein–Lagarias 1996 (the conjugacy map `Q`; from memory), the Kuratowski
note (sheet symmetry `T_+(-m) = -T_-(m)`, Kempe chains, the SHEET control),
the parallel spine note (Sturmian element of `E_inf`), the first note
(Proposition 9: the residue of a word's 2-adic point), Goldberg–West 1985 and
Alon 1987 (necklace splitting via Borsuk–Ulam; from memory), Lovász 1978
(Kneser graphs; from memory), Filos-Ratsikas–Goldberg 2018 (necklace
splitting is PPA-complete; from memory).

**Status: PROVED elementary (Theorem 1, Propositions 2–4) + FINITE-EXACT
(prices for odd `n < 10^6`; approach classes of the Sturmian point and of
2000 random no-descent fixed points; fixed points of `Q` to depth 44) +
CITED + DIRECTION. Collatz OPEN. D14 is resolved as a tautology (the cheap
growth of an integer is its own word) with a correction to the fifth note
logged in MISTAKES; the necklace mechanism is typed as sheet-blind, so no
Borsuk–Ulam-type argument can decide Collatz without breaking the sheet
symmetry by the sign of the intercept; the exact fixed-point structures
that exist are recorded.** Scripts
`04-computation/experiments/collatz_generic_price_topology_20260927.py` and
`..._fixq.py`, outputs beside them. Independent audit OWED (as for the
session's other notes).

## 0. The answer in one paragraph

Generic growth is cheaper than any cycle shadow: along the record
excursions below `10^6` a bit of net growth costs `1.14`–`1.49` bits of
source (`log_2 n / log_2(peak/n)`; peak exponents `1.67`–`1.88`), against
`1.71` for the `-1` shadow, and the smallest integers are cheaper still
(`27` buys `6.8` bits with `4.75`). The exact statement behind this is a
tautology, not a lack of price: for any 2-adic point `x`, every positive
integer in the depth-`K` class `x + 2^K Z` is at least `x mod 2^K`, so the
price of a depth-`K` approach is `K - z_K(x) - O(1)` bits where `z_K(x)` is
the run of zero bits of `x` just below position `K` (Theorem 1). For the
negative integer cycle points `z_K = 0`; for the trivial cycle the price is
zero; for the Sturmian point of `E_inf` and for the fixed points of random
no-descent words `z_K` is a geometric fluctuation. Hence a small integer `m`
approaches deeply only the points that look like it (`x = m mod 2^K`), the
cheap growth of `m` is the growth of its own word, and the residual
obstruction of the fifth note was misstated: it is not that generic points
are approachable by small integers, it is that an integer's own growth
class has it as least member, and counting such classes is Terras's density
theorem. The necklace-splitting proof has four moves (relax the cuts to a
sphere with a free `Z_2` action, build an odd imbalance map, force a zero by
Borsuk–Ulam, round to bead boundaries); Collatz has one `Z_2` symmetry, the
sheet involution `m -> -m`, and it conjugates `3n+1` to `3n-1`, which has
three cycles, so every argument invariant under it is false on one sheet
(Proposition 2): Borsuk–Ulam-type arguments are sheet-blind unless the sign
of the intercept enters. The rounding move is the cycle question, and its
failure for long words is Baker: Collatz is anti-rounding. What fixed-point
structure exists is exact and small: `Fix(T) = {0, -1}` on `Z_2`, and the
conjugacy map `Q` (parity vector as a 2-adic number) has `Fix(Q) ⊇ {0} ∪
{-2^j} ∪ {2^j/3}`, the doubling chains of `-1` and of the point `1/3` that
feeds the trivial cycle (Proposition 4; `Q(1/3) = 1/3` because `1/3 -> 1 ->
2 -> 1`), with `Q(1) = -1/3` and `Q(-1/3) = 1` a 2-cycle; to depth 44 a few
other long-lived solutions remain unidentified.

## 1. D14 resolved (PROVED + FINITE-EXACT)

**Empirical prices (script part 1).** For odd `n < 10^6`, with `peak(n)` the
maximum of the Syracuse orbit and `h(n) = log_2(peak/n)` the net growth:

| record `n` | peak | `log peak / log n` | bits of source per bit of growth |
|---|---|---|---|
| 4255 | 2270045 | 1.752 | 1.331 |
| 77671 | 523608245 | 1.783 | 1.277 |
| 159487 | 5734125917 | 1.876 | 1.142 |
| 270271 | 8216025965 | 1.825 | 1.212 |
| 704511 | 18997161173 | 1.758 | 1.320 |

The cheapest growth of at least 4 bits: `27` (`0.70` bits per bit, `6.8` bits
of growth), `31` (`0.75`), `41` (`0.86`), `47` (`0.92`), `55` (`1.00`). The
Lagarias–Weiss heuristic `peak <= n^(2+o(1))` says the price tends to `1`
along records; the `-1` shadow's price is `1/0.585 = 1.71`; generic words
are cheaper than any cycle shadow, and no pointwise bound is known (a
theorem `peak <= n^C` would be a strong form of no-divergence).

**Theorem 1 (universal size price).** Let `x in Z_2` and `K >= 1`, and let
`rho_K(x) = x mod 2^K in [0, 2^K)` be the least non-negative residue. Every
positive integer `m` with `m = x mod 2^K` satisfies `m >= rho_K(x)` (and `m >=
2^K` if `rho_K(x) = 0`). Writing `rho_K(x) = 2^(K - 1 - z_K(x)) (1 + o(1))`
with `z_K(x)` the number of zero bits of `x` at positions `K-1, K-2, ...`
down to the first one, the price of a depth-`K` approach to `x` is
`log_2 rho_K(x) = K - 1 - z_K(x) + O(1)` bits. In particular: for `x = -1,
-5, -17`, `rho_K = 2^K - |x|` and the price is `K - O(2^(-K))`; for `x = 1`
(the trivial cycle) `rho_K = 1` and the price is `0`; for `x` a positive
integer the price is `log_2 x` for all `K > log_2 x`.

*Proof.* The class `x + 2^K Z` meets the positive integers first at
`rho_K(x)` or `2^K`. ∎

FINITE-EXACT (script part 2): for the Sturmian point of `E_inf` (`d_k =
ceil(k log_2 3)`) the zero runs `z_K` at `K = 10, 20, 30, 40, 60, 80, 120, 160`
are `0, 1, 0, 6, 1, 0, 0, 0`; for the fixed points of 2000 random no-descent
words of length 30 at `K = 40`, `P(z_K >= z) = 1, 0.50, 0.26, 0.135, 0.070,
0.037, 0.017` for `z = 0..6`, i.e. the uniform law `2^(-z)`. So generic
points of `E_inf` are priced by size up to a geometric deficit, exactly like
the cycle points.

**Corollary 2 (the least-representative tautology).** An integer `m`
approaches `x` to depth `K > log_2 m + 1` iff `x = m mod 2^K`, i.e. iff the
bits of `x` from position `floor(log_2 m) + 1` to `K - 1` vanish. Hence the
only 2-adic points a small integer approaches deeply are the points that
look like it, the cheap growth of `m` is the growth of `m`'s own word, and
"an integer below `2^K` follows a no-descent word of precision `K`" means
"`m` is the least positive member of that class". The number of such
classes at precision `K` is the no-descent count `|W|` of THM-4495, and
their density `|W_k| 2^(-k) -> 0` is Terras's theorem; the pointwise
statement is the conjecture.

**Correction.** The fifth note's D14 said the points of `E_inf` that are
not negative integers are approachable by small integers; Theorem 1 shows
they are priced like everything else up to `z_K`. The direction is corrected
in place and the misstatement is logged in MISTAKES.

## 2. The necklace mechanism and its Collatz instances (PROVED typing)

The necklace-splitting theorem (two thieves, `k` bead types with even
counts, `k` cuts) is proved by: (i) relaxing cut configurations to points
of the sphere `S^k` (fractional cuts with signs), a space with a free
antipodal action; (ii) building the odd map `S^k -> R^k` of imbalances
between the two shares; (iii) Borsuk–Ulam forces a zero; (iv) rounding the
zero to bead boundaries because the measures are atomic and the dimension
count leaves room. Balance is forced by symmetry and dimension.

**Proposition 2 (odd-map arguments are sheet-blind).** The sheet involution
`iota(m) = -m` on `Z \ {0}` satisfies `T_+(-m) = -T_-(m)` (Kuratowski note,
Theorem F), where `T_+` is `3n+1` and `T_-` is `3n-1`. Any statement about
orbits that is invariant under `iota` holds for `3n+1` iff it holds for
`3n-1`; since `3n-1` has the three positive cycles through `1`, `5`, `17`,
no `iota`-invariant argument can conclude that every positive orbit of
`3n+1` reaches `1`. In particular an argument whose only topological input
is an odd map (Borsuk–Ulam, Tucker, ham sandwich, necklace splitting on the
sheets) cannot decide Collatz unless the sign of the intercept enters.

*Proof.* Conjugation by `iota`. ∎

This is the repo's SHEET control in topological clothing. The three other
moves of the necklace proof have Collatz counterparts, none useful: (i) the
relaxation of a discrete configuration to a continuous parameter is the
passage from an integer to a 2-adic point, and `Z_2` is a Cantor set on
which odd maps to `{+-1}` exist (the bit above the lowest set bit is flipped
by `x -> -x`), so no Borsuk–Ulam obstruction lives there; (ii) the natural
"imbalance" of a word is its height `L log_2 3 - d`, whose sign is
drift-driven, not symmetry-driven; (iv) the rounding move is the cycle
question: every word has a continuous zero, the rational cycle point `x_w`,
and it rounds to an integer only for the six known cycles, the failure for
long words being Baker's theorem (first note, section 3.5). Collatz is
anti-rounding. Necklace splitting itself does apply to cycle words: the
doubled `-17` word (14 beads, types `1, 2, 4`) splits with the cuts `(1, 2, 7)`
into two shares of type counts `{1: 5, 2: 1, 4: 1}` (script part 4, Alon's
bound `k(q-1) = 3` attained); the shares are unions of intervals, not orbit
segments, and carry no dynamics.

**Fixed-point theorems.** Brouwer, Lefschetz and Sperner produce a fixed
point of a self-map of a compact connected space; Collatz needs global
attraction to a known cycle. On `N` the only invariant compact sets are
finite unions of cycles (fourth note, Proposition 1), and on `Z_2` the
map is a shift (Lagarias): there is no stage on which a fixed-point
theorem could act. The PPA view (the parity argument on graphs of maximum
degree two, of which necklace splitting is a complete instance) has as
Collatz counterpart the Kempe chains of the Tait-coloured Collatz graph
(doubling orbits and rising runs, the Kuratowski note's Theorems K2–K4);
their parity bookkeeping is done there, and Collatz is a reachability
question in a functional graph, not a "find the other endpoint" question.

## 3. The fixed points that exist (PROVED + FINITE-EXACT)

Let `T(x) = x/2` (x even), `(3x+1)/2` (x odd) on `Z_2`, and `Q(x) = sum_i
w_i(x) 2^i` with `w_i` the parity of `T^i(x)` (the 3x+1 conjugacy map: a
measure-preserving homeomorphism conjugating `T` to the shift).

**Proposition 3.** `Fix(T) = {0, -1}` in `Z_2`.

*Proof.* `x/2 = x` gives `0`; `(3x+1)/2 = x` gives `-1`. ∎

**Proposition 4 (fixed points of the conjugacy).** (a) `Q(2y) = 2 Q(y)`, so
`Fix(Q)` is closed under doubling and under halving of even members. (b)
`Q(-1) = -1` and `Q(1/3) = 1/3`; hence `{0} ∪ {-2^j : j >= 0} ∪ {2^j/3 : j >=
0} ⊆ Fix(Q)`. (c) `Q(1) = -1/3` and `Q(-1/3) = 1`: a 2-cycle. (d) An odd
fixed point `y = 2z + 1` satisfies `Q(3z + 2) = z`; a node `(a, b)` of the
equation `Q(az + b) = z` is dead if `b` is odd and otherwise has the two
children `(a, b/2)` (for `z` even) and `(3a, (3a + 3b + 1)/2)` (for `z`
odd); the chain of `-1` is `(3^i, 3^i - 1)`, live at every depth.

*Proof.* (a) the word of `2y` is `0` followed by the word of `y`. (b) the
orbit of `-1` is constant and odd; `1/3` is odd and `T(1/3) = 1`, whose
orbit `1, 2, 1, 2, ...` gives `Q(1) = 1 + 4Q(1)`, i.e. `Q(1) = -1/3`, and
`Q(1/3) = 1 + 2Q(1) = 1/3`. (c) is (b). (d) `Q(y) = 1 + 2Q(T(y))` for odd `y`
and `T(2z+1) = 3z + 2`; the two cases of the parity of `z` give the two
children. ∎

FINITE-EXACT (`..._fixq.py`, depth 44). The number `N_k` of solutions of
`Q(x) = x mod 2^k` is `2, 4, 6, 10, 16, 22, 30, 44, 58, 68, 80, 96, 122,
144, ...` (linear growth, as for a critical branching tree conditioned to
survive). The solutions mod `2^10` that persist to depth 40 are exactly `0`,
`-2^j` (`j = 0..9`) and `2^j/3` (`j = 0..7`); at levels 20, 24, 28 the
persistent odd solutions are `-1`, `1/3` and three to six residues not
identified as rationals of height below 300 — long-lived branches of the
equation tree or genuine further fixed points, undecided here. Reading: the
self-describing 2-adic integers are the doubling chains (the Kempe chains
`<b,c>` of the Kuratowski note) of the two fixed points of `T` and of the
one preimage `1/3` of the trivial cycle; the conjugacy's own fixed points
say nothing about the integers.

## 4. Verdict

D14 is closed: the price of generic growth is universal up to a
geometric deficit, and what looked like a missing price is the tautology
that an integer's cheap growth is its own word; the pointwise content is
the conjecture and the average content is Terras. The topological reframe is
closed at the level of mechanisms: Borsuk–Ulam-type arguments are
sheet-blind (Proposition 2), the rounding move is the cycle question, and
fixed-point theorems have no stage. What remains exact and new: Theorem 1's
price formula, the fixed-point set of the conjugacy map (Proposition 4 and
the computation), and the necklace split of a cycle word as a typed
non-instance.

## 5. Directions (DIRECTION; none pursued)

* **D16.** Decide `Fix(Q)`: prove that the equation tree of Proposition
  4(d) has only the `-1` chain (and the `1/3` chain, which enters through
  the even branch `(a, b/2)`) as infinite paths, or exhibit a fourth
  family. A 2-adic self-consistency argument on `(a, b)` mod `2^s`, or a
  deeper persistence computation (depth 80), would settle the three
  unidentified residues at level 20.
* **D17.** The only symmetry-breaking parameter is the sign of the
  intercept. An odd function on `Z_2 \ {0}` exists (`f(x) = h_+(x) - h_-(x)`,
  the plus-sheet height minus the minus-sheet height, is odd by Theorem F),
  but the Cantor set forces no zero of an odd map; a Borsuk–Ulam argument
  would need a connected relaxation of `Z_2` on which `f` stays odd and
  continuous, and the real extensions of the map (Chamberland's
  interpolation) are not odd. Nothing to try without a new idea.

## 6. Reproduction

    cd 04-computation/experiments
    python3 collatz_generic_price_topology_20260927.py > collatz_generic_price_topology_20260927.out
    python3 collatz_generic_price_topology_20260927_fixq.py 40 > collatz_generic_price_topology_20260927_fixq.out

Standard library only; about six minutes for the first script (peaks of
`5·10^5` orbits).
