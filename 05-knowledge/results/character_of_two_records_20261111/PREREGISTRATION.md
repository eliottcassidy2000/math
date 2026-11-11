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
