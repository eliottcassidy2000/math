# Pre-registration - flip and spin for path counts (kind-pasteur-2026-11-11-S23)

Owner directive: another session on Chan-Pak, with pasted texts on "flippy" and "spinny" algebras (split-complex
versus complex, Cayley-Dickson doubling, octonions and sedenions, the scale 5, Fermat primes, A_5) as sources of
topics, to be combined with results of the repo. This file is written and pushed BEFORE the runs E1-E3 below.

## Working claims already derived or checked (not predictions)

Notation of THM-5130: pair `(B, O)`, state `(r, M, h, E)`, closing `C(B, O; X) = (B => X => v, v -> O)`;
`C_T(k)` = cover polynomial (ordered `k`-tuples of disjoint possibly empty paths covering `T`; `C_T(1) = H`,
`C_T(2) = Pi_2`; KNOWN: Session197 bin polynomial, THM-3166, THM-4756).

1. **Closing operator (checked by brute force on 360 cases, script `opcalc.py`).** For every pair `(B, O)` and every
   tournament `X`,

       C_{C(B,O;X)}(k) = k [ u_k C_X(k+1) + r_k C_X(k) + s_k C_X(k-1) ],

   where, over `(k+1)`-tuples `(L_1, ..., L_(k+1))` of disjoint possibly empty paths covering `B`: `u_k` counts those
   with `L_2` nonempty starting in `O`; `-s_k` counts those with `L_2` empty and `L_1` nonempty ending in `O`; and
   `r_k = a_k - 2 u_k + b'_k` (`a_k`: `L_2` empty or starting in `O`; `b'_k`: `L_2` nonempty starting in `O` and `L_1`
   empty or ending outside `O`). So `u_k >= 0 >= s_k`, `u_1 = M`, `r_1 = r`. A closer without back arc (`O` empty)
   acts by multiplication (`u = s = 0`); any back arc gives a second-order difference operator whose two shift
   coefficients have opposite signs. Examples: the 3-cycle (`B` a point, `O = B`) is
   `Phi = k (E + k - E^(-1))`; `(TT_2, both)`: `(2k+1, k^2 - 2, -(2k-1))`; `(C_3, one vertex)`:
   `(k^2+k+1, k^3+2k, -(k^2-k+1))`; `(C_3, two vertices)`: `(2(k^2+k+1), k^3, -2(k^2-k+1))`.
   The transitive tournament on `2j+1` vertices with the `j` nested longest arcs reversed has
   `(Phi^j k)(1) = 1, 3, 15, 123, 1479, 23955, 496479, 12725547` paths (`j <= 4` recounted literally).
2. **Flip and spin as bilinear forms (derived by hand).** In the coordinates `v_X = (H(X), Pi_2(X))` the ordered join
   is the componentwise product (the split algebra `Z x Z`), and the 3-cycle with two modules is
   `H(C_3[A, B, 1]) = x_A x_B + 2 h_A h_B` with `x = Pi_2 - H`, i.e. the real part of `z_A conj(z_B)` for
   `z = x + sqrt(-2) h` in `Z[sqrt(-2)]` (the quadratic form `x^2 + 2y^2` on the diagonal is Session260's). The pair
   (all of `A` selected) has state `(r, M) = (3 H(A) - Pi_2(A), Pi_2(A) - H(A))`.
3. **Module-side sources (derived by hand).** For a seed of slope `M` and any tournament `X`,
   `H(C(B, O; X)) = H(X) + M Pi_2(X)`. So a tournament `X` with `H(X) = 2^s + 1` and `Pi_2(X) = pi 2^e`, `pi` odd,
   gives the grid points `1 + 2^s + d 2^(r+s)` for every odd `d = pi c` with `c` the odd part of a seed slope, for all
   sufficiently large `r`. Known by hand: `C_3` (`s = 1`, `Pi_2 = 12`), the strong 4-tournament (`s = 2`,
   `Pi_2 = 32`), `C_3 => C_3` (`s = 3`, `Pi_2 = 144`).
4. **Cayley-Dickson dictionary (derived by hand; the algebra is classical).** A twisted group algebra on `Z_2^n` with
   `e_a^2 = -1` and `e_a e_b = -e_b e_a = +-e_(a+b)` is a tournament on the nonzero vectors. A line `{a, b, a+b}`
   spans a quaternion algebra iff it is a cyclic triangle; on a transitive line two of the three left
   multiplications are reflections. A zero divisor `(e_a + e_b)(e_c +- e_d) = 0` exists iff `a + b = c + d` and the
   four arcs between `{a, b}` and `{c, d}` contain an even number in each direction. The octonions are the Paley
   tournament on 7 points; counting gives 16 labelled copies of it with a given Fano plane as cyclic lines.

## What is on disk and what has not been computed

- Class lists of orders at most 9 hold `H` for every class but no value of `Pi_2`. No pair `(H, Pi_2)` beyond the
  three in claim 3 and the transitive ones has been computed.
- `mf_3.txt ... mf_8.txt` (scratch files of S21/S22) hold, for every kernel of order at most 8 and every ordered pair
  of ports, the coefficients of the family `al + be 2^a + ga 2^b + de 2^(a+b)`, from which the matrix
  `G = [[h11, h12], [h21, h22]]` of E2 can be read. Neither `G`, its determinant nor any discriminant has been
  computed or looked at. A reader may count E2 as not blind in the strict sense; I record my guesses anyway.
- Nothing about Cayley-Dickson sign tournaments has been computed.

## E1 - Fermat modules

For every class `X` of order at most 9: `(H(X), Pi_2(X))`. One-apex words of chain length at most 22 with
`H = 2^s + 1`.

- **F1.** For at least 7 of the 9 levels `s = 3, ..., 11` there is a tournament of order at most 9 with
  `H = 2^s + 1` and odd part of `Pi_2` dividing 9 (a module-side source for the 9-grid).
- **F2.** For at least half of the levels `3 <= s <= 16` there is a one-apex word of chain length at most 22 with
  `H = 2^s + 1` and odd part of `Pi_2` dividing 9.

## E2 - two-port kernels (not blind in the strict sense, see above)

`G` for kernels of order at most 6 whose two ports are visited at most twice.

- **K1.** Some kernel has `det G < 0`.
- **K2.** Among kernels with `h12 = h21` the form of discriminant `-4` (`b^2 - ac = -1`) occurs.

## E3 - Cayley-Dickson sign tournaments

- **C1.** The sign tournament of the sedenions (15 vertices, standard doubling) is doubly regular.
- **C2.** Of the 128 tournaments on the Fano plane all of whose lines are cyclic, exactly 16 have 189 paths (a count,
  not a guess), and exactly these have no even zero-sum quadrangle.

## Reading

F1 or F2 confirmed: module-side sources are a second supply for the grids, and the search for a uniform family
moves to words. K1 confirmed: the kernels come in both signatures and the split between definite and indefinite
forms is a structural invariant worth classifying; K1 refuted: every two-port kernel is definite or degenerate, a
positivity statement to prove. C1 refuted: the doubling leaves the doubly regular class at dimension 16, in line
with the appearance of zero divisors.

---

## Outcomes (added after the runs; nothing above this line was changed)

Runs: `runs/e1_pair_spectrum.txt`, `runs/e1_words_f2.txt`, `runs/e2_kernel_signatures.txt`,
`runs/e2_kernel_identities.txt`, `runs/e3_cd_signs.txt` and the follow-up files listed in `README.md`.

| prediction | outcome |
|---|---|
| F1 (at least 7 of the levels 3..11 have a class of order at most 9 with `H = 2^s + 1` and odd part of `Pi_2` dividing 9) | **FAILED**: 2 of 9 (levels 3 and 7) |
| F2 (at least half of the levels 3..16 have such a one-apex word of chain length at most 22) | **FAILED**: none of the levels 3..21 |
| K1 (some two-port kernel of order at most 6 has `det G < 0`) | **CONFIRMED**: none at orders 3, 4; 4 at order 5; 62 at order 6 |
| K2 (among symmetric kernels the discriminant -4 occurs) | **CONFIRMED**: the strong 4-tournament, `G = [[5, 2], [2, 1]]` |
| C1 (the sign tournament of the sedenions is doubly regular) | **CONFIRMED**, and proved for every `n` (README section 3.5) |
| C2 (of the 128 cyclic orientations of the Fano plane exactly 16 have 189 paths, and exactly these have no even zero-sum quadrangle) | **CONFIRMED**: 16 with (0 even 4-sets, `H = 189`), 112 with (6, `H = 171`) |

Caveats, in the direction of less credit:

1. **The prior-art sweep returned after this file was pushed.** Three of the four "working claims" were already in
   canon: claim 1 (closing operator on cover polynomials) is Theorem 1 of Session193 `profile_substitution/PROOF.md`,
   with the 3-cycle operator and the sequence 1, 3, 15, 123, ... in Sessions 159 and 197 and THM-4900; the bilinear
   law of claim 2 is eq. (1) of Session260 `QUASIPERFECT_NORMS.md` (I had credited only its diagonal); the count 16
   in claim 4 and in C2, with the values 189 and 171, is S14-S16 of
   `collatz_mod6_20260922_paley_fano_octonion_design.md`. So the first half of C2 was not a prediction, and the
   sentence "Nothing about Cayley-Dickson sign tournaments has been computed" was true of this session only.
   What was open in C2 is the second half (the 16 are exactly the orientations without an even quadrangle).
2. C1 is a consequence of a short induction (the sign table is a Hadamard matrix); I do not claim it is new.
3. E2 was declared not blind in the strict sense; K1 was moreover derivable by hand from THM-5130 (the anti-seed
   `(TT_2, both)` with its two ports joined is the hyperbolic kernel of order 5), which I saw only after the run.
4. The first run of E1 read the files `cls9_1.txt ... cls9_9.txt` (orders 1..9 of one generator run) as nine pieces
   of order 9: 7,412 duplicate rows, a wrong class count (206,360 for 198,948). No table changed; the stored run is
   the corrected one.
