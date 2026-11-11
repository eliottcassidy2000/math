# S91 (2026-11-10): the local densities of F1, the numbers of an Euler-Sombor abstract, and qualifying path counts

Session `opus-2026-11-10-S91`. Canon: THM-5128. Hypotheses: HYP-9660. Scripts and outputs:
`04-computation/experiments/s91_f1_local_densities_20261110/`. The Chan-Pak conjecture (every odd number except 7
and 21 is a Hamiltonian path count) is NOT proved here and remains OPEN; this session settles the local (residue
class) side of the best-understood join family and says precisely why that side cannot decide the conjecture.

Owner directive: "work a session on the local densities of F1 and see how prior work on the numbers 2, 3, 5, 7, 13,
19, 27 and 1039 relate deeply with symmetry structure we have been studying regarding tournaments and primes and
constructions like collatz which relate deeply to prime geometry. think of qualifying path counts in k-vertex
extensions of tournaments and how chan-pak can be proven by modifying that to Hamiltonian path counts in
tournaments themselves", with four attached papers (section 3).

## 0. Inheritance and portfolio

*Anchor:* the by-value local densities of F1 (obligation 2 of the S90 note). *Niche:* k-vertex extensions and rooted
("qualifying") counts. *Wildcard:* the four attached papers and the eight numbers.

Inherited, and credited throughout (the startup sweep found the first three, which S90 had missed: MISTAKE-758):

| What | Where | Status there |
| --- | --- | --- |
| Two favoured classes `0` and `-1`, confinement at `d = U = 0 (mod p)`, limits `1 + (3p-4)/(4p^2)`, `1 - 1/p^2 + 1/(4p)` | THM-4916 (kind-pasteur S15, 2026-11-07), letter-count model, one root change | PROVED; by-value limit "NOT CLAIMED" |
| Census factors, "Fermat boost", conjecture that value weighting closes the gap | HYP-9506 (A), (B), (C) | OPEN; pre-registered R2 FAILED |
| 3-adic law; "a joint real/residue transfer theorem ... is still required" | Session223 | PROVED at 3; the transfer theorem missing |
| Digit chain on `(v mod 2p, Q(v) mod p)`, uniform stationary law; by-value equidistribution of the one-apex count | Session206 | PROVED for fixed `p` |
| Mixing of `H mod M`; joint limit `(UV, R_M)` | Session193 modular_sync Theorems 3, 4; THM-4886 (C) | PROVED |
| Twisted transfer matrices of the count, rates (0.9320 at 1/7, `sqrt(3)/2` at 1/3) | S82 note section 7 | VERIFIED, `q <= 120` |
| Readout mass `(X/a) ln(1 + a/c)`, mean `(pi^2/12) ln N` for unordered pairs | THM-4885 | PROVED |
| By-value class means at 13 primes; presence / absence of the excess | HYP-9657 (A) and its audit (S90) | NUMERICAL |
| Mod-4 quadratic form of a vertex cube | THM-5127 (C) (mac-mini, 2026-11-10) | PROVED |

## 1. The local theory of F1 in one page (THM-5128)

Rooted state of the one-apex tournament with word `h`: `x = 2 psi(h) + 1` (its path count), `y = 2h + 1`. F1 is
`N = x1 x2 + 2^s y1 y2 - [s = 0]`: a bilinear form in two rooted states, on two sheets (`s = 0`, `s >= 1`).

1. **One biased class, carried by divisible states.** For any two independent states modulo `p` and a nondegenerate
   bilinear gluing of rank `r`, the class law is exactly
   `1/p + (pi_1 + pi_2 - p^(-r)) ([c = 0] - 1/p) + (second order)`, `pi_i` the mass of states that vanish modulo `p`.
   For F1 (`r = 2`) a state vanishes modulo `p` when `p | gcd(x, y)`; then `N` is forced into the class `0`
   (`s >= 1`) or `-1` (`s = 0`), whatever the other factor is. Everything else is flat.
2. **By value the law is first order up to a power saving.** Counting triples with `N <= X`:
   `share(c) = 1/p + ([c = c_s] - 1/p) F(X) + O(X^(-delta))`, `F(X)` the share of triples whose smaller state is
   p-divisible. The proof cuts the larger word into cells (top digits fix the size, bottom digits fix `x/y`, the
   middle digits mix the residues at the rate of the twisted transfer operator).
3. **The limit is generic, the approach logarithmic.** `F(X) -> 1/p^2`, so the local factors are
   `1 - 1/p^2 + [c = c_s]/p`, an Euler product over the odd primes (singular series
   `(8/pi^2) prod_{p | N} (1 + p/(p^2 - 1))` on the sheet `s >= 1`, with `N + 1` on the sheet `s = 0`). But a single
   p-divisible state `h` keeps a share of order `1/ln X` of the mass. So at every feasible scale the profile is
   decided by *how early a p-divisible state occurs*:
   - `x = y` exactly at `h = 2^j`, where `x = y = 2^(j+1) + 1`. The prime factors of `2^m + 1` therefore have a
     divisible state at `h = 2^(m-1)`: `3` (h = 1), `5` (2), `17` (8), `11` (16), `13` (32), `43` (64), `257` (128),
     `19` (at 256, and already at 218), `41`, `65537` (2^15). These primes show an excess far above the generic value.
   - Primes whose first divisible state is late (`29`: 1435, `37`: 1350, `31`: 8664, `61`, `73`, `127`: 109283) show
     none, or a deficit: without divisible states the distinguished class is *under*-represented by `1/p`
     (two non-zero vectors of `(Z/p)^2` are orthogonal with probability `1/(p + 1)`, not `1/p`).
4. **The prime 2 is exact.** Modulo `2^j` the twisted operator is `sqrt 2` times a unitary matrix: square-root
   cancellation with no loss, no 2-adic factor. The mod-4 law of the one-apex family is a band quadratic form whose
   symbol over `F_2` is `(1 + z)^2 (1 + z + z^2)`; its rank has period 6 in the word length.
5. **Tested before looking.** Predictions for the unseen window `[2^24, 2^27)` were computed from the list of
   divisible states and the uniform law of `x/y` alone, and pre-registered (commit 1f97e27805): 8 tests, 279
   intervals, all CONFIRMED. The Fermat prime 65537, never used before, shows `+0.14 / +0.17` (predicted
   `+0.18 / +0.21`, standard error `0.022`); its neighbours 65521 and 65539 show nothing.

| p | first divisible state | rate per digit | excess, sheet 0: measured / predicted | sheet 1: measured / predicted | generic limit |
| --- | --- | --- | --- | --- | --- |
| 3 | 1 | 0.8660 | +0.539 / +0.515 | +0.544 / +0.532 | 0.375 |
| 5 | 2 | 0.7418 | +0.3130 / +0.3120 | +0.3280 / +0.3284 | 0.208 |
| 7 | 38 | 0.9320 | +0.233 / +0.174 | +0.230 / +0.175 | 0.146 |
| 11 | 16 | 0.8587 | +0.1690 / +0.1705 | +0.1707 / +0.1749 | 0.092 |
| 13 | 32 | 0.8272 | +0.0441 / +0.0459 | +0.0446 / +0.0448 | 0.077 |
| 17 | 8 | 0.8289 | +0.2107 / +0.2081 | +0.2296 / +0.2288 | 0.059 |
| 19 | 218 | 0.8053 | +0.0510 / +0.0511 | +0.0495 / +0.0497 | 0.053 |
| 29 | 1435 | 0.7894 | +0.0051 / +0.0032 | +0.0020 / +0.0007 | 0.034 |
| 31 | 8664 | 0.9002 | +0.0041 / -0.0061 | -0.0038 / -0.0089 | 0.032 |
| 43 | 64 | 0.8342 | +0.0885 / +0.0877 | +0.0930 / +0.0941 | 0.023 |
| 257 | 128 | - | +0.1859 / +0.1815 | +0.2023 / +0.2044 | 0.004 |
| 65537 | 32768 | - | +0.144 / +0.183 | +0.168 / +0.206 | 0.00002 |

(Window `[2^24, 2^27)`; excess = distinguished class over the mean of the other classes; "predicted" used no window
data. At 3, 7, 31 the difference is the second-order term, which the theorem bounds by a power of `X`.)

## 2. What this settles, and what it corrects

- *HYP-9657 (A)* ("the classes 0 and -1 are favoured at some primes and not at others; the size is not `1/p`") is
  explained: presence, absence and size are the mass of p-divisible states. The S90 file overstated its novelty
  (the two classes and the confinement mechanism are THM-4916's): corrected there, MISTAKE-758.
- *The census question of THM-4916 / HYP-9506 for one root change.* By value the limit has the same shape as the
  letter-count law but a different split of the boost `1/p` between the classes `0` and `-1`:
  `(0.80359, 0.19641)` for word pairs against `(3/4, 1/4)`. The constants are dilogarithm values
  (`w_0 = (pi^2/12 + 2 Li_2(-1/2) - Li_2(-1/4)) / (pi^2/12)`), because `x/y` is uniform on `[1/2, 1]` and the sum over
  the padding telescopes. So value weighting does not "close the gap" to the letter-count limits: it replaces them.
- *The Kronecker question of S90* ("is the mass a class-number-type quantity?"): no. In the limit the local factor
  at `p` depends only on whether `p` divides `N` (or `N + 1`); at finite size it depends on the 2-adic order
  structure of `p` through the divisible states, not on a quadratic character of `N`.
- *The "Fermat boost"* (HYP-9506 (B)) is the first-order term of (C): quantified, and decaying like `1/ln X`.

## 3. The numbers 2, 3, 5, 7, 13, 19, 27, 1039

**Where they come from (EXACT).** They are the integers in the abstract of J. Allikvere, "Minimal Trees for the
Euler Sombor Index Revisited: A Correction and the First Thirteen Minima", MATCH 97 (2026) 1281-1293, one of the
four attached papers. The Euler-Sombor index of a graph is the sum over its edges of
`sqrt(i^2 + ij + j^2) = |i - j omega|`, `i, j` the end degrees: the length of an Eisenstein integer. `3, 7, 13, 19,
27` are the norms at the degree pairs `(1,1), (1,2), (1,3), (2,3), (3,3)`; `2` is the degree of the bulk (thread)
vertices, with baseline weight `2 sqrt 3`; `5` is the coefficient in `eps = 2 sqrt 19 - 5 sqrt 3 = sqrt 76 - sqrt 75`,
one of the two gaps that order the first thirteen minima; `1,039 = 2 + 5 + 13 + 33 + 89 + 240 + 657` is the number
of unicyclic graphs with 4 to 10 vertices used in the verification.

**The same numbers on the 2-clock (EXACT arithmetic; MECHANISM for F1).** `a^2 + ab + b^2 = (a^3 - b^3)/(a - b)`, so
`7 = 2^3 - 1`, `13 = (3^3 - 1)/2`, `19 = 3^3 - 2^3`, `27 = 3^3`. In base 2:

    3 = 2 + 1,      5 * 13 = 2^6 + 1,      27 * 19 = 2^9 + 1,      7 = 2^3 - 1,
    ord_p(2) = 2, 4, 3, 12, 18, 18, 519   for   p = 3, 5, 7, 13, 19, 27, 1039.

The list splits by the position of `-1` on the 2-clock. *Fermat type* (`-1` is a power of 2): `3, 5, 13, 19, 27`.
*Mersenne type* (odd order): `7` and `1039`. This is exactly the dichotomy of THM-5128:
- `2^m + 1` is the Hamiltonian path count of the transitive tournament with one reversed arc of span `m + 1` (the
  one-apex word `h = 2^(m-1)`, the only words with `x = y`). So a prime factor of `2^m + 1` has a balanced
  p-divisible state: `3` at `h = 1`; `5` and `13` together at `h = 32` (`x = y = 65`); `19` and `27` together at
  `h = 256` (`x = y = 513`). These states carry the first-order excess of F1 at those primes.
- A Mersenne-type prime never divides `2^m + 1`. Its divisible states are late and unbalanced (`7`: first at
  `h = 38`), and the short cycles of the doubling map make it a slow mixer: `7` is the slowest modulus found (rate
  `0.932` per digit), followed by divisors of `2^j - 1` with small `j` (`31`, `15`, `63`, `127`, `73`). They carry the
  second-order term.
So on the list the Fermat-type members feed the first-order term and `7` the second-order one. The doubling map
`eta -> 2 eta + b` on `Z/m` under the digit recursion is the same base dynamics as in the Collatz threads
(THM-4553, THM-4629: "characters hear near-returns"); `ord_27(2) = ord_19(2) = 18` with `2^9 = -1` modulo both is
the repository's mod-18 / mod-19 pairing, here realised inside one tournament (`H = 513`). The two missing values
are `7 = 2^3 - 1` and `21 = (2^6 - 1)/3`, the repunits `111` in bases 2 and 4 (recorded:
seven_twentyone_mersenne_openai_math_20261006.md); `(2^m - 1)/3` for even `m` are the numbers sent to a power of
two by `3x + 1`, `(2^m + 1)/3` for odd `m` by `3x - 1`, and the latter are a third of the balanced states.

**Quadratic forms in the join calculus (EXACT).** A one-arc join is a bilinear form in rooted states,
`H = x_h y_h - x_e y_e + x_R y_R` (signature `(2,1)`); on the diagonal (a one-apex tournament joined to itself with
padding `s`) it is `x^2 + 2^s y^2 - [s = 0]`, the norm form of `Z[sqrt(-2^s)]` at the rooted state: Gaussian for
`s = 0`, `Z[sqrt(-2)]` for `s = 1`. With `x, y` odd these are `1, 3, 5, 1, 1, ... (mod 8)`: a self-join is never
`7 (mod 8)`, and `7` and `1039` are `7 (mod 8)`. `27 = 3^2 + 2 * 3^2` is a self-join (`h = 1`, `s = 1`). The orders
that occur are those of discriminant `-4` and `-8`, the 2-power family. For gluings through a core `Z` of
interface vertices (`bridges/b7_selfjoin_eisenstein.py`: every core with `|Z| <= 4`, each form confirmed on
explicit tournaments; FINITE-EXACT) the literal Euler-Sombor form `x^2 + xy + y^2` never occurs. With one interface
vertex the diagonal is `x^2 + xy + 2y^2 - x - y - 1`, whose quadratic part is the norm form of `Q(sqrt(-7))`;
`x^2 + 3y^2` first occurs with two; a form properly equivalent to the Eisenstein form first occurs with three,
`3x^2 + 9xy + 7y^2 = X^2 + XY + Y^2` at `(X, Y) = (x + y, x + 2y)`, exactly when `Z` is a cyclic triangle, and
always with lower-order terms. So the triangle is where `Z[omega]` enters, and `7` enters before `3`. Where the
cubic `x^2 + x + 1` also enters is over `F_2`: in THM-5125 (4b) (nullity of `A^2 + A + I`), and, new, as a factor of the symbol
`(1 + z)^2 (1 + z + z^2)` of the mod-4 form of the one-apex family (THM-5128 (E)), which gives the period 6. So `3`
acts on these tournaments through `F_4 = F_2[omega]`, not through the lattice `Z[omega]`.

**As values and as moduli of F1.**

| n | H-value | one-apex words | `r_F1` | as a modulus: first divisible state, rate |
| --- | --- | --- | --- | --- |
| 3 | yes | h = 1 | 1 | 1, `0.866 = sqrt(3)/2` |
| 5 | yes | h = 2, 3 | 3 | 2, `0.742` |
| 7 | no (THM-1370) | none (`psi != 3`) | 0 | 38, `0.932` (slowest) |
| 13 | yes (order 5) | none (`psi != 6`) | 0 | 32, `0.827` |
| 19 | yes | h = 10, 13 | 4 | 218, `0.805` |
| 27 | yes | h = 17, 23 | 5 | (prime power; rate `0.866`) |
| 1039 | yes (least order 9) | none (`psi != 519`) | 0 | 405,729, `0.764` (58 divisible states below `2^26`, `62.2` expected) |

`1039 + 1 = 2^4 * 5 * 13` lies in the boosted class of the sheet `s = 0` at both `5` and `13`, and F1 still misses
`1039`: local factors are averages (65 of the 512 odd numbers in `[1025, 2049)` are F1 misses).

**The other three papers** (typed; details in the digests kept with the session).
- *Multiset dimension* (Allikvere, arXiv:2607.28813): the finite values are all positive integers except `2`: a
  theorem of Chan-Pak shape whose exclusion is a one-line symmetry (two landmarks can be swapped), where 7 and 21
  need a strong floor and deletion. LOOSE. Its strip proof (three-column blindness, a transfer automaton on column
  words, tight transitions as a regular language) has the architecture of the cut-state automata of THM-4878 and
  of the cells of THM-5128. MECHANISM; lost: min-plus optimisation instead of counting. One exact arithmetic
  fact was added by this session (`bridges/b4_strip_automaton.py`, FINITE-EXACT, three independent confirmations
  of the spectrum): the recurrent core of that automaton (19 states) has characteristic polynomial
  `x^14 (x - 3)(x - 1)(x + 1)(x^2 + x + 2)`, so the paper's counts satisfy `a(n) = 2a(n-1) + a(n-2) + 6a(n-3)`
  and `42 a(n) = 32 * 3^n - 9 V_n - 21 U_n` with `V_n^2 + 7 U_n^2 = 2^(n+2)`: the second factor is the minimal
  polynomial of `(-1 + sqrt(-7))/2`, of norm 2, and `|U_n| = 1` exactly at `n = 1, 2, 3, 5, 13` (checked to
  `n = 20,000`; the Ramanujan-Nagell equation `x^2 + 7 = 2^m`). So `Q(sqrt(-7))`, the field of the one-vertex
  core above and of the repository's `{1, 2, 4} mod 7` thread, is also the field of this strip count.
- *Edge multiset dimension of hypercubes* (arXiv:2608.09983; already THM-4525, THM-4534). The one-apex count
  depends on its word only through the multiset of gaps between a 1 and a later 0 (THM-4876 (6)): a
  distance-multiset invariant on a path. EXACT restatement, nothing new.
- *Malik's leaf-constrained minimisers*: a `2 x 2` exchange inequality with a classification of equality, as in
  THM-5124, for an additive index. MECHANISM; lost: `H` is not additive over local pieces.
- *Allikvere's correction* itself: "bulk weight times size plus a finite sum of defect penalties", and an error
  that came from giving an adjacent-defect link the bulk weight. F1 has the same anatomy (`2^s` for the padding,
  the correction `-[s = 0]` when the two defects touch), and yesterday's over-general sentence about it was
  repaired the same way. MECHANISM; lost: additive versus bilinear.

## 5. Where this leaves the proof

1. *The local side is closed for F1 and is not where the difficulty lies.* Local densities exist, form an Euler
   product, are positive for every odd `N`, and are a thin-interface effect: for a gluing of seam rank `r` the bias
   is `p^(-r)` (THM-5128 (A)). The two missing values are not local.
2. *What F1 does at feasible sizes.* Mean multiplicity `0.748 ln N`; relative variance about `0.75` on a sheet; miss
   rates `1.38 %`, `0.82 %`, `0.62 %`, `0.35 %` on `[2^18, 2^21)`, `[2^21, 2^23)`, `[2^23, 2^24)`, `[2^24, 2^27)` (and
   `0.087 %` on `[2^31, 2^32)`, THM-4875): a decay like a small power of `N`. The misses sit where the singular series
   is smallest: `N = 1 (mod 3)`, the only class without a boost at 3, is missed three times as often (S78), and
   integers with `N (N + 1)` free of small primes are the lightest (HYP-9657 (B)); the Euler product has minimum
   `8/pi^2 = 0.811` of the mean, and the S82 audit measured `0.809` on the `N` coprime to `255255`.
3. *What completeness still needs is unchanged:* a collision (second-moment) bound, the minor-arc side. What this
   session adds to that program is the major-arc side for any bilinear join family (a count of divisible states)
   and, for one state, a level of distribution with a uniform power saving (THM-5128 (B4)).
4. *The extension lane* is section 4.

## 6. Principles, named (triggers and counter-indications)

- **Confinement is the only local bias of a bilinear gluing.** Trigger: a residue-class profile of a join family
  looks irregular in `p`. Action: do not fit a formula in `p`; list the states that vanish modulo `p` and weigh them
  (THM-5128 (A)): exactly one class is biased and its excess is that mass. Counter-indication: at slowly mixing
  moduli (`7`, `31`, `3`, `127`, `73`) the second-order term is visible at every feasible size.
- **Top digits fix the size, bottom digits fix the ratio, middle digits mix.** For the one-apex state, `y = 2h + 1`
  is read off the top digits, `x/y = 1 - phi(h)/2` off the bottom digits (binary radical inverse), and the residues
  modulo `p` are mixed by the digits in between at the rate of the twisted transfer operator. Trigger: a statement is
  known in a letter-count model and wanted by value. Counter-indication: the by-value limit need not equal the
  letter-count limit; here the shape is the same and the constants are not.
- **On a hyperbola, small states never leave: finite-size effects decay like `1 / ln X`.** A fixed state keeps a
  share of the mass of order `1/ln X`, so the profile at any computable size is dominated by the first few
  divisible states. Trigger: an "excess that decays slowly with N". Action: compute the first-order term from the
  small states before extrapolating.
- **Thin interfaces carry a singular series, thick ones do not.** The bias of a rank-`r` gluing with equidistributed
  states is `p^(-r)`: visible for `r = 2` (F1), negligible for rich grammars. Local densities therefore cannot
  explain the two missing values 7 and 21, which are global.
- **The base prime is exact.** At 2 the digit recursion is unitary up to `sqrt 2`: use Rudin-Shapiro-type identities
  there, not mixing estimates.

## 7. Reproduction

    cd 04-computation/experiments/s91_f1_local_densities_20261110
    python3 twisted_spectrum.py 16        > twisted_spectrum.out         # rates, checked against brute force
    python3 twisted_spectrum_2adic.py 18  > twisted_spectrum_2adic.out   # the prime 2: unitary, quadratic form
    python3 mod4_symbol.py                > mod4_symbol.out              # rank of the band form, period 6
    python3 box_identity.py               > box_identity.out             # exact class law (rationals), boxes
    python3 divisible_states.py 24 300    > divisible_states.out         # p-divisible states, limit law of x/y
    python3 constants.py                  > constants.out                # dilogarithm constants, sheet shares
    python3 first_order_exact.py 18 20 22 24 26 > first_order_exact.out  # the two forms of the first-order law
    python3 window_split.py 21 18         > window_split_insample_21_18.out
    python3 window_split.py 23 21         > window_split_insample_23_21.out
    python3 window_split.py 27 24 write-prereg prereg_window_27_24.json   # predictions only
    python3 window_split.py 27 24 prereg_window_27_24.json > window_split_prereg_27_24.out   # about 90 s, 1.5 GB

Every script exits nonzero on a failed check. numba is required for the last three.
