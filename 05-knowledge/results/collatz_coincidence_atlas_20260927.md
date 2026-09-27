# A coincidence atlas: the polygonal triangle's hidden region is the minus sheet, Euler's pentagonal exponents are indexed by the Syracuse core, the golden solids fold to Petersen and K_6, and the number 139 across cycles, tournaments and modular curves

**Session:** opus, `collatz-posets-zeta5-20260927` (third note), 2026-09-27.
**Owner's directive:** relate the `139` numerology of the second note to the
repo's Langlands reading, the tournament "hidden values" and the discrete
graph-theory triple `{Petersen, K_5, K_{3,3}}`; merge the golden-ratio facts
(Lucas rounding, `phi^n = (L_n + F_n sqrt5)/2`, the icosahedron's three
golden rectangles, the cube in the dodecahedron with volume ratio
`2/(2+phi)`, `[Q(zeta_5) : Q(sqrt5)] = 2`) and the triangular arrangement
`{1},{2,1},{3,3,1},{4,6,4,1},{5,10,9,5,1},{6,15,16,12,6,1},{7,21,25,22,15,7,1}`
with its edges, its "region of interest" `{6},{10,9},{15,16,12},{21,25,22,15}`
and the "hidden region behind the triangle" with the single value `2`; look
for microcosm/macrocosm bridges wherever coincidences occur, and say which
are exact.
**Inherits (read, cited):** the S8 crossings note
[`collatz_crossings_20260926_potential_and_seeds.md`](collatz_crossings_20260926_potential_and_seeds.md)
section 1.2 (the `{2,3,11}` theorem, `139 = 3^7 - 2^11` pairing the
exponents `7, 11`, and the Langlands reading through `X_0(11)` and the
theta series, with its verdict "nothing here connects it to the Collatz
map"); the Kuratowski note
[`procgen_kuratowski_20260925_tait_kempe_triples.md`](procgen_kuratowski_20260925_tait_kempe_triples.md)
(Kohl's group as a Tait-coloured Collatz graph; "dodecahedron to Petersen
is an ANALOGY"; the triple dictionary "dual pair + self-dual" and "twins +
container"; no finite excluded-minor list); the flip-rank thread
(HYP-3798, HYP-3805: `k(7) = 12`, the best 11-arc configuration's `2^11`
completions reach 454 of 456 classes and miss the Paley heptagon of score
sequence `3^7`); HYP-9162 (the Sierpinski tournament tower, `T_3` = the
Paley heptagon, `Aut = F_21`); THM-261 (Petersen = the `A_4` root
orthogonality graph), THM-262 (tournaments as `so(n)` and `A_(n-1)`
carriers); HYP-3728 (Ihara discriminants `9 - 4n` of `K_n`, Heegner at `n =
3, 4, 5, 7, 13`); HYP-3771/3772 (the five Platonic solids from the `(2,3,n)`
angle-defect spine and the genus of `X_0(2p)`); the wave-10 note's sheet
symmetry `T_+(-m) = -T_-(m)`; the first two notes of this session.

**Status: EXACT identifications (section 1: the triangle, its hidden region
and the two sheets; section 2: the eta indexing; section 3: the golden
facts and the antipodal folds, all verified by the script) + CITED
placements against the repo (sections 4–5) + NUMEROLOGY typed (section 5:
`139`, and the tournament `3^7` is a score sequence, not a power) +
DIRECTION (section 7). Nothing here bears on the truth of the Collatz
conjecture; the exact statements are dictionaries. Script
`04-computation/experiments/collatz_coincidence_atlas_20260927.py`, output
beside it. Independent audit OWED (as for the session's other notes).**

## 0. The answer in one paragraph

The owner's triangle is the table of polygonal numbers read by antidiagonals,
`T(n, j) = P_(j+2)(n - j)`: the left edge is the "2-gonal" numbers (the
naturals), the right edge `P_k(1) = 1`, the next diagonal `P_k(2) = k` (the
shifted naturals), and the region of interest `1 <= j <= n - 3` is the
triangular, square, pentagonal, hexagonal, ... numbers from their second
term on. The hidden region behind the triangle is exact: the generalized
polygonal numbers at negative index, `P_k(-m) = P_k(m) + (k - 4) m`, a mirror
triangle whose pentagonal column `2, 7, 15, 26, ...` begins with the owner's
single value `2 = P_5(-1)`. This is where the owner's "analogous to `3n-1`"
becomes literal: `P_5(n) = sum_(j<n) (3j + 1)` and `P_5(-n) = sum_(j<=n)
(3j - 1)` are the partial sums of the plus and minus sheets' affine forms,
Euler's two families of pentagonal exponents are the two sheets, `24 P + 1 =
(6n -+ 1)^2`, and `eta(tau) = sum chi_12(m) q^(m^2/24)` is indexed by the odd
numbers coprime to 3, which are exactly the Syracuse nodes with odd
predecessors, split into `6k - 1` (predecessor valuations odd) and `6k + 1`
(even). That is the one microcosm/macrocosm bridge in this note that is
exact on both ends: the partial sums of `3j -+ 1` on the small side, a weight
`1/2` modular form on the large side, the same mod-6 split in between. The
golden facts are all exact (verified in `Q(sqrt5)`), and their graph shadows
are the antipodal folds dodecahedron to Petersen and icosahedron to `K_6`,
which remove exactly the golden eigenvalues `+-sqrt5`; the Kuratowski note
already records that the Collatz sign fold is not such a covering, so the
Petersen side of the owner's triple stays an analogy there. The number `139
= 3^7 - 2^11` is the `-17` cycle's denominator (7 odd steps, 11 halvings, `S
= 17 * 139`), the S8 Langlands reading of the pair `(7, 11)` stands as
written, and in the tournament thread `3^7` is the *score sequence* of the
Paley heptagon while `2^11` counts subcube completions: the subtraction is a
pun of notation, so the tournament face of `139` is numerology, typed as
such.

## 1. The triangle, its edges and its hidden region (EXACT)

Let `P_k(n) = ((k-2) n^2 - (k-4) n)/2` be the `n`-th `k`-gonal number,
defined for all integers `n` (the generalized polygonal numbers).

**Proposition 1.** (a) The owner's rows are `T(n, j) = P_(j+2)(n - j)` for
`0 <= j < n` (checked rows 1–7). (b) Edges: `T(n, 0) = P_2(n) = n`; `T(n,
n-1) = P_(n+1)(1) = 1`; `T(n, n-2) = P_n(2) = n`, the shifted naturals the
owner saw beside the ones. (c) The region of interest `1 <= j <= n-3`
consists of the `k`-gonal numbers for `k >= 3` from their second term on:
`{6}, {10, 9}, {15, 16, 12}, {21, 25, 22, 15}` are triangular, square,
pentagonal and hexagonal; its rows have `1, 2, 3, ...` entries, so its size is
a triangular number. (d) The hidden region: `P_k(-m) = P_k(m) + (k - 4) m`;
the squares are their own mirror, the triangular numbers shift down by one
index (`P_3(-m) = T_(m-1)`), the pentagonal column becomes `2, 7, 15, 26`
(the second pentagonal numbers), the hexagonal column `3, 10, 21, 36 =
T_(2m)`. The apex is `P_5(-1) = 2`, the only value of the hidden region
that the visible triangle also shows at its own apex row (`T(2,0) = 2`).

*Proof.* Substitution in the closed form; the script checks (a)–(d)
exactly for `k, m < 30`. ∎

**Proposition 2 (the two sheets).** `P_5(n) = sum_(j=0)^(n-1) (3j + 1)` and
`P_5(-n) = sum_(j=1)^(n) (3j - 1)`: the visible pentagonal numbers are the
partial sums of the plus sheet's affine form at `j = 0, ..., n-1`, and the
hidden ones the partial sums of the minus sheet's form at `j = 1, ..., n`.
Moreover `24 P_5(n) + 1 = (6n - 1)^2` and `24 P_5(-n) + 1 = (6n + 1)^2`.

*Proof.* Arithmetic progressions with difference 3; the square identities
expand. Checked to `n = 60`. ∎

So the owner's "hidden region behind the triangle" is the minus sheet in
the only sense the numbers allow: the same quadratic form read at negative
index, with `3n - 1` in place of `3n + 1`. The value `2` at the apex is `3 *
1 - 1`, the first minus-sheet value. This is a dictionary, not a
mechanism: the polygonal numbers are values of quadratic forms, the Collatz
sheets are affine maps; what they share is the pair of affine forms `3j +-
1` and the sign flip `n -> -n` that exchanges them (on the Collatz side
this is the wave-10 note's exact sheet symmetry `T_+(-m) = -T_-(m)`).

## 2. Euler's pentagonal exponents and the Syracuse core (EXACT)

**Proposition 3.** (a) `prod_(n>=1) (1 - q^n) = sum_(k in Z) (-1)^k
q^(k(3k-1)/2)` (Euler; checked to order 400). (b) With `m = 6k - 1` for `k >
0` and `m = 6|k| + 1` for `k < 0`, the exponent is `(m^2 - 1)/24` and the sign
`(-1)^k` equals the character `chi_12(m)` (`+1` for `m = +-1 mod 12`, `-1`
for `m = +-5 mod 12`); this is the classical `eta(tau) = q^(1/24) prod (1 -
q^n) = sum_(m>=1) chi_12(m) q^(m^2/24)`. (c) The index set `{6k +- 1}` is the
set of odd integers coprime to 3, which is exactly the set of Syracuse nodes
possessing odd predecessors: `m = 6k - 1` has the predecessors `(2^v m -
1)/3` with `v` odd, `m = 6k + 1` with `v` even, and multiples of 3 have none
(checked for `m < 400`).

*Proof.* (a) is Euler's theorem (the script verifies the coefficients); (b)
is Proposition 2 plus the residues of `6k -+ 1` mod 12; (c) is `2^v m = 1
mod 3`. ∎

The microcosm/macrocosm reading. On the small side, the pentagonal numbers
are partial sums of `3j -+ 1`; on the large side, `eta` is a weight-`1/2`
modular form whose Fourier support is the Syracuse core `6k +- 1` and whose
signs are a quadratic character; the modularity `eta(-1/tau) = sqrt(tau/i)
eta(tau)` is the same micro/macro duality the S8 note recorded for the
owner's theta series `theta_3^2`, now for the theta series of the core
squares. What is exact: the index set and its mod-6 split coincide with the
predecessor-parity split of the Syracuse tree. What is not: no
Collatz-dynamical statement is encoded in `eta` (the map is affine, the
exponents are quadratic, and the dynamics on `6k +- 1` mixes the classes
freely). The bridge is the residue class, and it stops there. The
character level is also where the tournament thread touches this one: the
Paley tournament's arcs are the quadratic character mod 7 and `eta`'s signs
the character `chi_12`; both are Dirichlet characters, and that is the whole
of the shared mechanism.

## 3. The golden facts and their graph shadows (EXACT)

Verified exactly in `Q(sqrt5)` or numerically to twelve digits (script part
3): `phi^n = (L_n + F_n sqrt5)/2`; `round(phi^n) = L_n` for `n >= 2` (with
`phi^0, phi^1` against `L_0 = 2, L_1 = 1` reversed, as the owner says);
`zeta_5 + conj(zeta_5) = 1/phi = phi - 1`, a golden integer, so `[Q(zeta_5) :
Q(sqrt5)] = 2` inside the degree-4 cyclotomic field; the cube inscribed in
the dodecahedron has edge `phi` times the dodecahedron's (a pentagon
diagonal) and volume ratio `phi^3 / ((15 + 7 sqrt5)/4) = 1 - sqrt5/5 = 2/(2 +
phi)` exactly; three mutually perpendicular golden rectangles have the twelve
icosahedron vertices `(0, +-1, +-phi)` and cyclic permutations as corners.

**Proposition 4 (antipodal folds).** The dodecahedron graph (20 vertices,
spectrum `3, sqrt5^3, 1^5, 0^4, -2^4, -sqrt5^3`) has antipodal quotient the
Petersen graph (spectrum `3, 1^5, -2^4`; Petersen is determined by its
spectrum), and the icosahedron graph (spectrum `5, sqrt5^3, -1^5, -sqrt5^3`)
has antipodal quotient `K_6` (spectrum `5, -1^5`). The golden eigenvalues
`+-sqrt5`, of multiplicity 3 (the two three-dimensional icosahedral
representations), are exactly the eigenvalues killed by the fold.

*Proof.* Built from the coordinates in the script, quotient by `v ~ -v`,
spectra compared. Standard (the dodecahedron is the canonical double cover of
Petersen, the icosahedron of `K_6`). ∎

Placement. Petersen and `K_6` are the two best-known members of the
Petersen family (the linkless-embedding obstructions); Petersen is THM-261's
`A_4` root-orthogonality graph, and THM-262 carries 5-vertex tournaments in
the same `A_4` root lattice, so the owner's triple sits in the `A_4`/`S_5`
world whose golden double cover is the icosahedral `H_3`. The Kuratowski note
has already tested the one Collatz-side candidate for such a fold: the
sign quotient of the two sheets is a fold, not a covering, so "dodecahedron
to Petersen" is an ANALOGY for Collatz (its Theorem F), and Petersen
obstructs nothing in a single sheet (every component of a functional graph
has at most one cycle, so the sheet is planar). Nothing in this note changes
that verdict.

## 4. The `139` thread, placed (CITED + NUMEROLOGY)

| where `139` or `(7, 11)` appears | what it is there | status |
|---|---|---|
| Collatz, the `-17` cycle | `p = 7` odd steps, `A = 11` halvings, `-17 = S/(2^11 - 3^7)`, `S = 2363 = 17 * 139`; `11/7` is the mediant of `3/2` and `8/5`, a semiconvergent of `log_2 3` | EXACT (second note, section 6; Steiner's route) |
| S8 crossings note, `{2, 3, 11}` theorem | `139 = 3^7 - 2^11` pairs the exponents `7, 11`; Langlands reading: `11` is the conductor of `X_0(11)` (genus 1; `X_0(7)` has genus 0, `X_0(139)` genus 11), theta series `theta_3^2` as micro/macro duality; verdict "nothing connects it to the Collatz map" | CITED; verdict unchanged |
| flip-rank thread (HYP-3798, HYP-3805) | `k(7) = 12`; the best 11-arc configuration's `2^11` completions reach 454 of 456 classes and miss the Paley heptagon, written `(3^7)` for its score sequence (seven vertices of score 3) | CITED; `3^7 - 2^11` is a pun: a score sequence minus a count |
| HYP-9162 tower | `T_3` = the Paley heptagon, `Aut(T_k) = F_21` for `k = 3..9`; `139 = 3 mod 4`, so `P_139` exists with `|Aut| = 139 * 69 = 9591` | CITED / EXACT arithmetic; no relation to the cycle |
| Zenodo zeta(5) preprint | decay constant `139/5` | NUMEROLOGY (second note) |
| HYP-3728 | Ihara discriminants `9 - 4n` of `K_n`: `K_5` gives `-11` (Heegner) | CITED; `11` as an integer |

Verdict. The pair `(7, 11)` has three independent causes: in Collatz the
near-miss `2^11 < 3^7` of `log_2 3` (so the `-17` cycle closes with
denominator 139); in tournaments the arc count `C(7,2) = 21` and the free-arc
count 11 of a configuration that happens to miss the most symmetric
tournament; in the modular world the smallest prime levels of genus 0 and
1. No mechanism links them, and the note records the coincidence as the
owner asked, typed NUMEROLOGY. The one exact cross-thread fact is
arithmetic: the same semiconvergent `11/7` that makes the `-17` cycle close
is the reason `2^11` and `3^7` are near each other at all.

## 5. Microcosm/macrocosm, in one table

| microcosm | macrocosm | the bridge | exact? |
|---|---|---|---|
| partial sums of `3j -+ 1` (pentagonal numbers, both signs) | `eta(tau)`, weight `1/2`, support `(6k +- 1)^2/24`, signs `chi_12` | the index set `6k +- 1` = the Syracuse core, split by predecessor parity | EXACT dictionary, no dynamics |
| the `-17` cycle `(7, 11)` | the continued fraction of `log_2 3`; Baker's bounds (shape (B)) | `2^A - 3^p` divides the carry `S` | EXACT (classical) |
| Petersen, `K_6` (Petersen family), `A_4`, `S_5` | dodecahedron, icosahedron, `H_3`, `Q(sqrt5)`, `Q(zeta_5)` | antipodal double covers; the golden eigenvalues are the covering's kernel | EXACT (Proposition 4); for Collatz an ANALOGY (Kuratowski note, Theorem F) |
| the Paley heptagon, `F_21` | quadratic characters, theta series with character | both `eta` and Paley arcs are governed by a Dirichlet character | ANALOGY (shared object type, no shared theorem) |
| the five Platonic solids | the `(2,3,n)` angle-defect spine; `X_0(2p)` genus (HYP-3771/3772) | the sign of `1/2 + 1/3 + 1/n - 1` | CITED synthesis, not new here |
| the polygonal triangle's hidden region | the generalized polygonal numbers `P_k(-m)` | `P_k(-m) = P_k(m) + (k-4) m` | EXACT (Proposition 1) |

## 6. What this does for Collatz

Nothing on the conjecture, and the note says so. The two exact bridges (the
sheets in the pentagonal numbers; the Syracuse core as the support of
`eta`) are residue-class facts: they see `m mod 6` and the sign of `m`, which
is exactly the information the repo's SHEET control shows to be
insufficient (the `3n - 1` sheet has the same mod-6 structure and three
cycles). Any use of the modular side would have to feed the *dynamics*, not
the residue classes, into a modular object; the first note's Proposition 10
(the per-orbit series `F_n`) and the wave-10 basin series are the two
candidates on record, and neither is modular.

## 7. Directions (DIRECTION; none pursued)

* **D9.** The hidden-region identity `P_k(-m) = P_k(m) + (k-4) m` says the
  visible and hidden `k`-gonal numbers differ by a multiple of `m` with
  slope `k - 4`; for the pentagonal case the slope is 1, the sheet shift. A
  quadratic analogue of the sheet symmetry (`T_+(-m) = -T_-(m)`) for the
  Collatz *carries* `S_(L-1)` (which are sums of `3^a 2^b`, not quadratic
  forms) does not exist; the pentagonal dictionary is a dictionary of
  values, not of maps.
* **D10.** The theta series of the Syracuse core, `sum_(m coprime to 6)
  chi(m) q^(m^2)`, could be twisted by orbit data (for instance by the sign
  of the first descent) to test whether any orbit statistic is modular; the
  expected answer is no, and the test is cheap (compute the first few
  hundred coefficients and check for a functional equation numerically).
* **D11.** The Paley heptagon's `F_21` and the `-17` cycle's `(7, 11)` share
  the number 7 for unrelated reasons; a real test of "hidden values" in the
  owner's sense would be whether the tournament tower `T_k` of HYP-9162 has
  any Collatz reading beyond the Gilbreath sea (the S10 note found none),
  and this note adds none.

## 8. Reproduction

    cd 04-computation/experiments
    python3 collatz_coincidence_atlas_20260927.py > collatz_coincidence_atlas_20260927.out

Needs numpy for the spectra; everything else is exact integer or
`Q(sqrt5)` arithmetic.
