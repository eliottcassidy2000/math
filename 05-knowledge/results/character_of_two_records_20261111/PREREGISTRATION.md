# Pre-registration (mac-mini-2026-11-11-atoms)

Written and committed before the runs named below were started. Notation: `H(T)` = number of Hamiltonian paths;
`c_k` = number of directed `k`-cycles; `P_q` = Paley tournament on `GF(q)`, `q = 3 (mod 4)` (`a -> b` iff `b - a` is
a non-zero square); `psi(x) = -x`, an anti-automorphism of `P_q` fixing only `0`; `F(T, psi)` = number of directed
odd cycles `C` with `psi(C) = C` (they pass through the fixed vertex); `F_k` = the number of those of length
`2k + 1`. KNOWN (Session29 converse-congruences, section 3): `H = 1 + 2F (mod 4)` for a tournament with an
anti-involution with one fixed vertex. `CD_k` = the sign tournament of the Cayley-Dickson algebra of dimension `2^k`
(THM-5133 (G), kind-pasteur S23; it is the doubling tower of THM-448). `G_35(R) = sum_X (-1)^((c_3 + c_5)(R + v_X))`
(THM-5132 (3)).

## What is known when this is written (in-sample)

- `H(P_q)` for `q = 3, 7, 11, 19, 23`: `3, 189, 95095, 1172695746915, 15760206976379349` (repository); residues
  mod 4: `3, 1, 3, 3, 1`, i.e. `1` for `q = 7 (mod 8)` and `3` for `q = 3 (mod 8)`. `H(CD_4) = 198335025 = 1 (mod 4)`.
- Proved on paper before any run: if `a -> psi(a)` iff `a -> v_0` for every `a != v_0` (no invariant cyclic triple),
  then `F` is even (fixed-point-free involution `P -> psi(reverse(P))` on the half-paths); this holds for `P_q` with
  `q = 7 (mod 8)`. For `q = 3 (mod 8)` every invariant triple `{0, a, -a}` is cyclic, `F_1 = (q - 1)/2` is odd, and
  `F_2`, `F_3` are even by regularity and double regularity (hand computation for prime `q`).
- The kernel formula of THM-5132 (3) and the classical facts on binary quadratic-residue codes suggest R4.
- Nothing has been computed for `q = 27, 31, 43, 47`, for the self-converse census of Run B, for the kernels of
  `A(CD_k)` over `F_2`, or for anti-automorphisms of `CD_4`, `CD_5`.

## Run A: invariant odd cycles of Paley tournaments (`paley_sym.c`)

Exact dynamic programme over (set of antipodal pairs used, last vertex) for `q = 3, 7, 11, 19, 23, 27, 31, 43, 47`.

- **R1.** `F(P_27)` and `F(P_43)` are odd (so `H(P_27) = H(P_43) = 3 (mod 4)`).
- **R2.** For `q = 11, 19, 27, 43`, `F_k(P_q)` is even for every `k >= 2`.
- **Control (proved).** `F(P_q)` is even for `q = 7, 23, 31, 47`; `F(P_q)` is odd for `q = 3, 11, 19` (known `H`).

Conjecture to be recorded if R1 holds: `H(P_q) = q(q-1)/2 (mod 4)` for every prime power `q = 3 (mod 4)`, i.e.
`(-1)^((H(P_q) - 1)/2) = (2/q)`; equivalently the number of orbits of Hamiltonian paths under `Aut(P_q)` is
`1 (mod 4)` for prime `q`.

## Run B: one-fixed-point anti-involutions in general (`selfconv_census.py`)

All self-converse classes of orders 5, 7, 9 (one anti-involution each, they are conjugate), with the numbers of
invariant odd cycles by length and the pattern of cyclic invariant triples.

- **R3.** If every invariant triple `{v_0, a, psi(a)}` is cyclic, then `F = (n - 1)/2 (mod 2)`, i.e. `H = n (mod 4)`.
  (A guess from Paley; it may fail.)
- Exploratory, not registered: the law of `F mod 2` for mixed patterns.

## Run C: pentagon cube sums of Paley and Cayley-Dickson tournaments (`family_cubes.py`, kernel formula)

- **R4.** `G_35(P_q) = 2^q` for `q = 7, 23, 31, 47` (every one-vertex extension of `P_q` has `c_3 + c_5` even), and
  `G_35(P_q) = 2^((q+1)/2)` for `q = 3, 11, 19, 43`.
- **R5.** `G_35(CD_k) = 2^(2^k - 1)` for `k = 3, 4, 5`: the kernel of the adjacency matrix of `CD_k` over `F_2` has
  dimension `2^(k-1) - 1` and all its weights are `0 (mod 4)`.

## Run D: the Cayley-Dickson tower mod 4

- **R6.** `CD_4` and `CD_5` are self-converse, and `H(CD_5) = 1 (mod 4)` (through `F` even).

## Refutation

Any of R1-R6 outside its statement is reported as FAILED. A failure of a control refutes the program or the proof.

---

# Second block (same session): R7

Written after Runs A-D were read (R1, R4, R5 CONFIRMED; R2, R3 FAILED as stated and replaced by a proved law; R6
FAILED: `CD_4` has no anti-automorphism) and before Run E was started.

In-sample: `H(P_q) = 1 (mod 4)` for `q = 7 (mod 8)` and `3 (mod 4)` for `q = 3 (mod 8)` is now proved (quotient
formula); `H(CD_4) = 198335025 = 1 (mod 4)` although `CD_4` is not self-converse.

## Run E: doubly regular tournaments of order 15 from one skew Hadamard matrix (`drt_descendants.py`)

For each of the 16 indices `c` of the Cayley-Dickson sign matrix of order 16, normalise row and column `c` and take
the tournament on the other 15 indices (the "descendants"; each is doubly regular).

- **R7.** Every descendant has `H = 1 (mod 4)`.

Conjecture to be recorded if R7 holds: every doubly regular tournament of order `n` has `H = 1 (mod 4)` if
`n = 7 (mod 8)` and `H = 3 (mod 4)` if `n = 3 (mod 8)`.

---

# Third block (same session): R8

Written after Run E was read (R7 holds but is uninformative: all 16 descendants are isomorphic to `CD_4`) and before
Run F was started. `H(P_27)` was not found in the repository (grep of canon and hypotheses for the value).

In-sample: `H(P_q) mod 8 = 3, 5, 7, 3, 5` for `q = 3, 7, 11, 19, 23`, which equals `q(q-1)/2 mod 8` in all five
cases (so the number of `Aut`-orbits of Hamiltonian paths is `1 (mod 8)` for these primes). Mod 4 this is now a
theorem; mod 8 it is an observation.

## Run F: the Paley tournament on `GF(27)` (`paley27.c`)

Number of Hamiltonian paths starting `0 -> 1 -> ...` by a layered subset dynamic programme modulo 256;
`H(P_27) = 27 * 13 * (that number)`.

- **R8.** `H(P_27) = 27 * 13 = 351 = 7 (mod 8)`.

Conjecture to be recorded if R8 holds: `H(P_q) = q(q-1)/2 (mod 8)` for every prime power `q = 3 (mod 4)`.

---

# Fourth block (same session): R9, R10

Written after Runs A-F were read and the note was drafted (the folding formula, Theorem Q, is proved; so is
`H(P_q) = q(q-1)/2 (mod 4)` for every prime power `q = 3 (mod 4)`), and before Run G was started. Checked first:
canon has no statement on `H mod 4` of doubly regular tournaments and no path count of a non-Paley doubly regular
tournament of order 19, 23 or 27 (grep of THM-447, 448, 451, 4557 and of "doubly regular" with "mod 4" / "mod 8").

In-sample: Paley tournaments (`q = 3, 7, 11, 19, 23, 27`: `H = q(q-1)/2 (mod 8)`; mod 4 proved for all `q`) and one
non-Paley class, `CD_4` of order 15 with `H = 198335025 = 1 (mod 8) = 105 (mod 8)`. The two doubly regular
tournaments of order 15 are `CD_4` and its converse, so order 15 is exhausted. Nothing is known here about a
non-Paley doubly regular tournament of order 19, 23 or 27.

## Run G: doubly regular tournaments that are not Paley (`drt_zoo.py`, `hampaths_layered.c`)

- **G1 (order 19).** All skew Hadamard matrices of order 20 given by the Goethals-Seidel array with four circulants
  of order 5, the first of skew type; all 20 descendants of each (normalise a row and its column, delete them); exact
  `H` by subset dynamic programming. If every tournament obtained is isomorphic to `P_19`, G1 is uninformative and a
  search for another doubly regular tournament of order 19 (SAT, random seeds) replaces it.
- **G2 (order 23).** The double `[[M, M], [-M^T, M^T]]` of the skew Hadamard matrix `M` of order 12 of `P_11`; its
  24 descendants; `H mod 256` by the layered subset programme.
- **G3 (order 27).** Goethals-Seidel with circulants of order 7; for the descendants that have a converse symmetry,
  `H mod 4` by the folding formula. (Without a converse symmetry: not computed, unless time allows the direct
  programme.)
- Controls: the programmes must return `H(P_7) = 189`, `H(P_11) = 95095`, `H(P_19) = 1172695746915`,
  `H(CD_4) = 198335025`, and `H(P_23) = 5 (mod 8)`.

- **R9.** Every tournament of G1 has `H = 3 (mod 4)`; every tournament of G2 has `H = 1 (mod 4)`; every tournament
  of G3 whose residue is computed has `H = 3 (mod 4)`. (The doubly regular form of the Paley law: `H = n(n-1)/2
  (mod 4)`.)
- **R10.** Every tournament of G1 has `H = 3 (mod 8)` and every tournament of G2 has `H = 5 (mod 8)`
  (`H = n(n-1)/2 (mod 8)`).

A single doubly regular tournament outside these residues refutes the corresponding statement; R9 failing would show
that the Paley law of this note is a statement about the quadratic character and the converse symmetry, not about
double regularity.

---

# Fifth block (same session): R11

Written after the sampled Run G was read (G1 on every 40th Goethals-Seidel solution: one class, not Paley,
self-converse, `H = 1172732961771 = 3 (mod 8)`; G2: one class, not Paley, not self-converse, `H = 5 (mod 8)`; G3: 39
classes from 40 matrices, three of them self-converse with `H = 3 (mod 4)` by folding; R9, R10 hold so far; the full
G1 is running) and after an exploratory look at switching classes, which is in-sample:

- In-sample (exploratory, not registered): let `T^+` be `T` with a new source. For `T = P_3, P_7` (all members) and
  `P_11, CD_4, P_19` (60, 60, 30 random members) every member of the switching class of `T^+` has the same `H mod 4`
  (3, 1, 3, 1, 3). `H mod 8` is not constant on the class of `CD_4^+` (1 and 5 occur). For general tournaments
  `H mod 4` is not a switching invariant (random tests, orders 5-10).
- Proved between the blocks (two lines, by counting the vertices between `a` and `psi(a)`): in a doubly regular
  tournament of order `4t + 3` with a converse symmetry, every invariant triple is cyclic if `t` is even and none is
  if `t` is odd; with the folding formula, `H = n(n-1)/2 (mod 4)` for every self-converse doubly regular tournament.

Members of the switching class of `T^+` are exactly the tournaments `U` on `N = n + 1` vertices whose `+-1` skew
matrix `S` satisfies `S^2 = -(N-1) I` and which are switching-equivalent to `T^+`; the class is the "skew Hadamard
class" of `T`.

## Run H: switching classes at new orders and direct counts at order 27 (`switch_class.py`, `drt27_direct.py`)

- **H1 (`N = 24`).** 40 random members (seeded) of each of the classes of `P_23^+` and of `D(P_11)^+` (the class of
  G2); `H mod 256` by the layered programme.
- **H2 (`N = 20`).** 40 random members of the class of the non-Paley doubly regular tournament of order 19 of G1;
  exact `H`.
- **H3 (order 27, direct).** For the three self-converse classes of G3 and the first six classes of G3 without a
  converse symmetry: `H mod 256` as the sum over the start vertex of 27 layered programmes on 26 vertices. Control:
  the same method on `P_27` must return 95.

- **R11.** Every member in H1 has `H = 1 (mod 4)`; every member in H2 has `H = 3 (mod 4)`.
- **R9 (continued).** Every tournament of H3 has `H = 3 (mod 4)`; the three self-converse ones agree with the
  folding formula. **R10 (continued).** Every tournament of H3 has `H = 7 (mod 8)`.

Conjecture to be recorded if R11 and R9 hold: for every tournament `U` on `N` vertices with `S^2 = -(N-1) I`,
`H(U) = 1 + N/2 (mod 4)`.

---

# Sixth block (same session): R12

Written after Runs G, H1, H2 and part of H3 were read and the note was drafted (pre-audit checkpoint `1ca9c5f83b`;
the independent audit is running), and before Run I was started. Checked first: the Sierpinski lane
(`04-computation/experiments/sierpinski_ocf_20261022/drt*`) and canon contain no cycle counts by length for doubly
regular tournaments of order 19, 23 or 27 (grep for `c_7`, "7-cycles", "cycle counts", "odd cycles" in the lane
files and for "doubly regular" near cycle counts in canon).

In-sample (exploratory, `cycle_parities` run of this hour): for every member examined of the switching classes of
`P_3^+` (all 8), `P_7^+` (all 128), `P_11^+` (40 random) and `CD_4^+` (`CD_4` itself and 12 random members) the number
`c_(2k+1)` of directed cycles of each odd length has the parity of `[N = 4 mod 8] C(N/2 - 1, k)`: the parity vectors
`(c_3, c_5, ...)` are `1`, `000`, `10011`, `0000000`. Proved: the same for every doubly regular tournament with a
converse symmetry (Theorem E), and for `k = 1, 2` for every tournament whose scores have one parity (THM-5132).

## Run I: cycle counts by length (`cycle_lengths.c`, `cycle_parities.py`)

Numbers of directed cycles of every length modulo 256, by a subset programme rooted at the smallest vertex of the
cycle.

- **I1 (order 23).** The thirteen classes of HYP-9371 (`drt23_classes.json`: `P_23` and twelve classes without a
  converse symmetry) and the double of `P_11`.
- **I2 (`N = 24`).** Ten random members (seeded) of each of the switching classes of `P_23^+` and `D(P_11)^+`.
- **I3 (`N = 20`).** `D_19`, `P_19` with a source, and twenty random members of each of their switching classes.
- **I4 (order 27).** `P_27`, the three self-converse classes and the six classes without converse symmetry of H3.
- Controls: `c_3 = n(n^2 - 1)/24` for regular tournaments; the Hamiltonian cycle counts `24, 5505, 34358763933,
  374127973957716` of `P_7, P_11, P_19, P_23` (THM-4747) modulo 256; `1 + 2 (sum of the odd cycle counts) = H (mod 4)`
  where `H` is known.

- **R12.** In I1 and I2 every odd cycle count is even (lengths 3 to 23). In I3 the parities of
  `(c_3, c_5, ..., c_19)` are `C(9, k) mod 2 = 1, 0, 0, 0, 0, 0, 0, 1, 1` for every member. In I4 the parities of
  `(c_3, ..., c_27)` are `C(13, k) mod 2 = 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1` for every tournament.

Conjecture to be recorded if R12 holds: if `S^2 = -(N-1) I` then `c_(2k+1) = [N = 4 mod 8] C(N/2 - 1, k) (mod 2)` for
every `k`. One tournament with one odd length outside this refutes it (the total parity may then still agree).

---

# Outcomes (appended after the runs; the blocks above are unchanged)

Outputs: `out/`. Program names changed after registration: `paley_sym.c` is `halfpaths.c`, `selfconv_census.py` is
`antiinv_census.c`, `paley27.c` is `hampaths_layered.c`.

| | verdict | what happened |
|---|---|---|
| R1 | CONFIRMED | `F(P_27)`, `F(P_43)` odd. Now a corollary of the folding formula. |
| R2 | FAILED | `F_k(P_q)` is odd for `k = 4, 5` already at `q = 11`; the parities are `C((q-1)/2, k) mod 2`. |
| R3 | FAILED | With all invariant triples cyclic, `F` is odd at every order tried (`n = 3, 5, 7, 9, 11`), not `(n-1)/2 mod 2`. Replaced by `F_k = C(m, k) (mod 2)`, proved. |
| R4 | CONFIRMED | Kernel formula, elimination and (orders <= 11) direct count agree. Now a theorem for all doubly regular tournaments. |
| R5 | CONFIRMED | `k = 3, 4, 5` (also `6, 7`). Same theorem. |
| R6 | FAILED | `CD_4`, `CD_5` have no anti-automorphism. **This part was decidable from canon** (THM-448, THM-451: the tower is not self-converse beyond order 7) and should not have been registered. `H(CD_5) mod 4` was not computed. |
| R7 | CONFIRMED, uninformative | All 16 descendants are isomorphic to `CD_4`. |
| R8 | CONFIRMED | Paths starting `0 -> 1` in `P_27`: `1 (mod 256)`; `H(P_27) = 95 (mod 256)`. |
| R9 | CONFIRMED | G1: all 12,000 Goethals-Seidel matrices give one class, `D_19`, `H = 1172732961771 = 3 (mod 4)`. G2: one class, `H = 189 (mod 256)`. G3: 39 classes from 40 matrices, three with a converse symmetry, `H = 3 (mod 4)` by folding. H3: nine classes of order 27 (three with, six without converse symmetry), `H = 127, 191, 111, 151, 55, 23, 247, 127, 23 (mod 256)`, all `3 (mod 4)`; the folding formula agrees on the first three. |
| R10 | CONFIRMED | G1: `3 (mod 8)`; G2: `5 (mod 8)`; H3: `7 (mod 8)` for all nine. |
| R11 | CONFIRMED | 40 random members of each class: `H = 1 (mod 4)` for both classes at `N = 24`, `3 (mod 4)` at `N = 20`. |
| R12 | CONFIRMED | I1: 14 of 14 classes of order 23 have every odd cycle count even. I2: 20 of 20. I3: 42 of 42 with the vector `100000011`. I4: 10 of 10 with the vector `1001100110011` (six classes without converse symmetry). Extension after R12 was read, same criterion, not separately registered: all 39 classes met at order 27 (`out/cycle_parities_i4_all.out`), 39 of 39. The vertex-local form (cycles through one vertex) fails, so the law is global. |

**Registration errors.** (1) R6, as above. (2) R9 and R10 at orders 19 and 23: the repository already had the exact
path counts of `D_19` (HYP-9365) and of the doubled Paley class of order 23 (HYP-9371,
`04-computation/experiments/sierpinski_ocf_20261022/drt23_doubled_check.out`), written by this machine on
2026-10-22; the check before registering looked at four theorem files and for the words of the predicted
property, not for the objects. So G1 and G2 were not out-of-sample. The same lane holds the exact counts of twelve
further classes of order 23 (all `5 mod 8`), used in the note as in-sample evidence. Out-of-sample in the strict
sense: G3, H1, H2, H3.

**Deviations.** G1 was first run on every 40th solution (`out/drt_zoo_g1_stride40.out`), then on all
(`out/drt_zoo_g1.out`). G3 used every 2000th of the 246,960 solutions, 40 matrices, all 28 descendants of each. An
exploratory look at switching classes (`out/switch_class_insample.out`) preceded the fifth block and is declared
there as in-sample.

**Controls.** All held: `F(P_q)` even at `q = 7, 23, 31, 47` and odd at `3, 11, 19`; exact `H` of `P_7, P_11, P_19, CD_4`;
`H(P_23) = 213 (mod 256)`; `H(P_27) = 95 (mod 256)` by the 27-start method.
