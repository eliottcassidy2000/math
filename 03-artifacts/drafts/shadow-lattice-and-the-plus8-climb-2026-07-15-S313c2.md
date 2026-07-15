# The (r, g) Shadow Lattice, the Missing 32nd Region at Every Gap, and the +8 Climb as a 3-Cycle Counter

**Instance:** klein-2026-07-15-S313 (cont.2) — Part II of
`vandermonde-tail-polygonal-polyhedral-2026-07-15-S313.md`
**Status:** PROVED (master theorem has a one-line proof; all claims verified in
`04-computation/shadow_lattice_gap_rank_klein_S313.py`, 14/14 checks, output in
`05-knowledge/results/shadow_lattice_gap_rank_klein_S313.out`)
**Claims:** HYP-6946 (master theorem + missing-region law), HYP-6947 (new-sequence bank),
HYP-6948 (k(T) = c3 and the solvability inversion), T1533.

---

## 1. The two knobs: rank r (how much dimension) and gap g (how steep the diagonal)

Part I had two objects (polygonal = rank-1, simplex = rank-∞) and two functionals (row sums
g = 0, shallow diagonals g = 1). Both generalize:

- **Rank-r figurate table** `T_r(k, j) = Σ_{i=0}^{r} C(k−1, i)·C(j, i+1)` (column 0 = ones).
  r = 1: polygonal numbers. r = 2: a genuine interpolating family ("polyhedral-gonal";
  agrees with tetrahedral numbers at k = 3 exactly, then bends flatter). r ≥ k−1: simplex.
- **Gap-g sums** `D_{r,g}(d) = Σ_k T_r(k, d − (g+1)k + 1)` — g = 0 rows, g = 1 skipping one
  row, g ≥ 2 skipping more (the new axis requested).

**Full-rank reference (r = ∞):** the gap-g sums of Pascal are the g-bonacci family
`a(d) = a(d−1) + a(d−g−1)`, GF `1/(1 − x − x^{g+1})`:
2^n (g=0), Fibonacci (g=1), Narayana's cows A000930 (g=2), A003269 (g=3), A003520 (g=4).
Verified.

## 2. Master theorem (uniform over the whole plane, one-line proof)

With layer GFs `L_i^{(g)} = x^{(g+2)(i+1)−1}/[(1−x^{g+1})^{i+1}(1−x)^{i+2}]`,

    D_{r,g}(x) = [ (1−x)^{r+2}(1−x^{g+1})^{r+1} − x^{(g+2)(r+2)−1} ]
                 / [ (1−x)^{r+2}(1−x^{g+1})^{r+1}(1 − x − x^{g+1}) ],

and the numerator is **divisible by the g-bonacci kernel** `1 − x − x^{g+1}`:

*Proof.* At any kernel root, `x^{g+1} = 1 − x`, hence `1 − x^{g+1} = x` and
`x^{(g+2)(r+2)−1} = (x^{g+1})^{r+2}·x^{r+1} = (1−x)^{r+2}·x^{r+1}` — equal to the other
term. ∎

So the exponential pole cancels for EVERY rank and EVERY gap, and

> **D_{r,g} is a quasi-polynomial of degree 2r + 2 with period g + 1.**

(The Part-I certificate `x⁸ ≡ 13 − 21x (mod 1−x−x²)` is the (1,1) shadow of this; the
Fibonacci numbers in it were the value of `(1−x)^{r+2}x^{r+1}` in disguise.) The geometric
series over layers uses only `(1−x^{g+1})(1−x) − x^{g+2} = 1 − x − x^{g+1}` — the same
identity — so the whole structure is one algebraic move repeated.

## 3. The missing-region law: Moser's 32nd region at every (r, g)

Define the deficit `t_{r,g}(d) = full_g(d) − D_{r,g}(d)` (what the shadow misses).

> **Missing-region law (verified r ≤ 3, g ≤ 4, and forced by the GF):** the deficit is 0
> for `d < d₀` and **exactly 1** at `d₀ = (g+2)(r+2) − 1`.

Every shadow has its own "missing 32nd region" moment, always of deficit exactly one —
the first firing of Vandermonde mode r+1, a single unit cell, before the tail floods in.
The r = 1 row of first misses:

| g | d₀ | full | shadow | the "missing region" |
|---|----|------|--------|----------------------|
| 0 | 5 | 32 | 31 | **Moser's classical missing 32nd region** |
| 1 | 8 | 34 | 33 | the missing 34th (Fibonacci shadow, Part I) |
| 2 | 11 | 41 | 40 | missing 41st (Narayana shadow) |
| 3 | 14 | 50 | 49 | missing 50th |
| 4 | 17 | 60 | 59 | missing 60th |

In the circle picture, Moser's missing region is mode i = 2 firing once — the unit cell
where the tetrahedron first outgrows the square (10 vs 9 at (j,k) = (3,3), Part I §1a).
The law says this is not an accident of rows: every gap and every rank truncation has its
single phantom cell, at the predictable address `(g+2)(r+2) − 1`.

## 4. The new-sequence bank (all OEIS-checked 2026-07-15)

**Not in OEIS (new; submission candidates), joining D = D_{1,1} from Part I:**

| name | first terms | law |
|------|-------------|-----|
| `D_{1,2}` | 1,1,1,2,3,4,6,9,13,19,28,**40**,56,78,106,141,186,241,307,388,484,596,729,883 | quartic qp, period 3; leaves Narayana at d = 11 |
| `D_{1,3}` | 1,1,1,1,2,3,4,5,7,10,14,19,26,36,**49**,65,85,111,143,181,226,281,346,421 | quartic qp, period 4; leaves A003269 at d = 14 |
| `D_{1,4}` | 1,1,1,1,1,2,3,4,5,6,8,11,15,20,26,34,**45**,59,76,96,120,150,186,228 | quartic qp, period 5; leaves A003520 at d = 17 |
| `D_{2,1}` | 1,1,2,3,5,8,13,21,34,55,89,**143**,228,358,554,841,1255,1837,2644,3739,5207,7140,9659,12893 | SEXTIC qp, period 2; tracks Fibonacci one gap deeper (through F(12)) |
| `D_{2,2}` | 1,1,1,2,3,4,6,9,13,19,28,41,60,88,**129**,188,272,391,556,781,1086,1492,2024,2717 | sextic qp, period 3 |
| `t_{1,1}` (deficit) | 1,4,13,33,76,159,315,594,1084,1923,3343,5714,9645,16115,26720,44035 (d ≥ 8) | GF `x⁸/((1−x²)²(1−x)³(1−x−x²))` — "Fibonacci minus quasi-quartic" |
| hyper-pyramidal Moser (s=2) | 1,2,6,17,43,99,211,421,793,1420,2432,4005 | `1 + C(n+3,4) + C(n+3,6)` |
| CD-level c3-max | 1, 5, 30, 204, 1496 | `(n³−n)/24` at `n = 2^k+1` = `C(2^k+2, 3)/4` |

**Found in OEIS (cross-connection):** the s = 1 integration ladder ("pyramidal Moser")
`1,2,5,12,27,57,113,211,373,628,1013,1574,2367,3459,…` = **A362193** (Grassmannian
permutations avoiding a pattern), shifted by one — verified against the b-file through
3459. Pattern-avoiding Grassmannian permutations = row sums of the pyramidal-number
triangle: a bijection is presumably readable off their formula; backlogged.

**Integration ladder (proved):** s-fold partial-summing the polygonal columns before
building the triangle gives row sums `1 + C(n+s+1, s+2) + C(n+s+1, s+4)` — Moser (s=0),
A362193 (s=1), new at s ≥ 2. Truncation and integration commute: the ladder is the
rank-1 shadow of Pascal's own integration ladder.

**Brown completeness:** every `D_{1,g}` (g ≤ 4) and `D_{2,1}`, `D_{2,2}` passes Brown's
criterion — the whole lattice consists of no-holes complete sequences. Completeness is
robust to BOTH knobs.

## 5. The +8 climb is a 3-cycle counter: k(T) = c3(T)

**Theorem (proved; exhaustive n = 5, random n ≤ 9).** For any tournament,
`x(T) = Σ_v d_v² = (n³−n)/3 − 8·c3(T)` where c3 = number of cyclic triangles.
Consequently THM-866's canonical distance to transitivity is

> **k(T) = ((n³−n)/3 − x(T))/8 = c3(T)** — the +8 quantum is literally ONE cyclic
> triangle, and every tie-splitting flip destroys exactly one net 3-cycle.

*Proof.* `Σd² = 4Σs² − 4(n−1)Σs + n(n−1)² = 4Σs² − n(n−1)²` (scores sum to C(n,2)).
Kendall–Babington Smith: `c3 = C(n,3) − Σ C(s_v,2)`, i.e. `Σs² = 2C(n,3) + C(n,2) − 2c3`.
Substitute: `x = 8C(n,3) + 4C(n,2) − n(n−1)² − 8c3 = (n³−n)/3 − 8c3`. ∎
(Ingredients classical; the reading of THM-866's lattice coordinate as an exact 3-cycle
count, and of tie-splits as single-3-cycle annihilations, is the new content. Walk
verified step-by-step: every tie-split has Δx = +8 AND Δc3 = −1 simultaneously.)

So the axis level lattice of THM-866 **is** the 3-cycle-count lattice `{0, 1, …, c3_max}`
read backwards, and "no holes in the levels" = "every 3-cycle count between 0 and
(n³−n)/24 is realized" — a Rédei-flavored intermediate-value theorem for c3.

## 6. The solvability inversion (extends 07-reflections/the-solvability-tower s703)

The S703 reflection identified: A₅ is perfect because **two cyclic triangles sharing one
vertex** regenerate each other; the smallest carrier is the round tournament C₅; hence the
quintic threshold IS a tournament. The new results sharpen this into a clean dichotomy:

- **Tournament dynamics always sorts.** The tie-splitting walk kills one 3-cycle per move
  and reaches the transitive order in exactly c3(T) moves (§5). In the tournament world,
  the 3-cycle atoms can be annihilated one at a time, always.
- **Group dynamics cannot.** In A₅ the same atoms (3-cycles) interlock — `[(123),(345)]`
  is again a 3-cycle — so the derived series never descends. Abel–Ruffini.
- **The precise divider is parity of symmetry.** |Aut(T)| is odd for every tournament (an
  involution would have to reverse the arc inside a swapped pair; exhaustive n ≤ 5
  reverified; classical). By **Feit–Thompson, every tournament automorphism group is
  solvable**: the unsolvability obstruction is structurally un-hostable as tournament
  symmetry. The quintic's A₅ (order 60, even) can act on no 5-vertex tournament; the
  solvable-quintic criterion Gal ⊆ F₂₀ = AGL(1,5) has largest odd subgroup ℤ₅ — which is
  exactly Aut(C₅), realized by the n = 5 floor. **Solvable quintics have tournament-sized
  symmetry; unsolvable ones have symmetry that no tournament can carry.**
- **At n = 5 the floor is the quintic tournament itself**: all 24 labelled floor
  tournaments (x = 0) are copies of C₅, each with c3 = 5 = k_max(5) and |Aut| = 5 =
  max over all n = 5 tournaments. The farthest-from-sorted point of the n = 5 landscape,
  the LRC cyclotomic witness, and the quintic threshold coincide.

**Three towers, one shape** (with opus-S316): the radical/derived tower (Galois), the
2-adic OCF digit tower above Rédei (S316: digit 0 constant, digit 1 dies at n = 7 —
"the OCF is itself the Toda tower"), and the +8 level tower (THM-866, = c3 descent).
The first two have perfect/irreducible cores where descent halts; the third always
descends — because its atoms are counted with multiplicity (a lattice ℤ-grading), not
conjugated against each other (a group commutator). The x-axis is what the derived series
would be if 3-cycles commuted.

## 7. Cayley–Dickson resonances (numerology — flagged, NOT claims)

At the repo's CD levels n = 2^k + 1 (R→C→H→O→S at n = 2,3,5,9,17):

| n | 3 | 5 | 9 | 17 | 33 |
|---|---|---|---|----|----|
| ceiling (n³−n)/3 | 8 | 40 | **240** | 1632 | 11968 |
| c3_max = k_max | 1 | 5 | **30** | 204 | 1496 |

Observed: at the octonion level n = 9 the transitive ceiling is **240 = #E8 roots** and
c3_max = **30 = #icosahedron edges**; at the quaternion level n = 5, c3_max = 5 = |Aut(C₅)|
and the binary icosahedral group 2I (order 120 = |S₅|) sits inside the unit quaternions,
with the icosian ring generating E8 (Conway–Sloane) — the classical bridge that would have
to carry any real mechanism. The quintic's topological obstruction (π₁ of the Poincaré
sphere = 2I, perfect) lives at exactly this H-level. No mechanism is claimed; recorded
because the icosian/E8 chain is the one place where "8 per 3-cycle" and "rank-8 lattice"
share an address. The honest content stands in §§5–6 without it.
(Also flagged: `x = Σd²` is a quadratic form and its mod-8 level arithmetic smells of
Gauss–Milgram signatures; OPEN, no construction.)

## 8. Where this touches the live problems

The master theorem's mechanism — *kernel divides numerator ⟺ exponential mode cancels,
leaving a quasi-polynomial* — is the exact algebraic shape hunted in the LRC(14) signed
Fourier walls (the memory "the closing estimate must be signed" is about forcing such a
cancellation); here we possess the cancellation with a one-line certificate because the
truncation is a Vandermonde prefix. S703's inversion (solvable = cyclotomic = TIGHT in
LRC) composes with §6: the LRC-hard configurations are the C₅-like ones — maximally
symmetric, floor-level, odd-Aut — i.e. **LRC hardness concentrates exactly where
tournament symmetry is largest and where the +8 climb starts farthest from sorted.**
