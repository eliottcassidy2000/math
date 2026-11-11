# Literature and OEIS check for the audit of `mod4_race_square_term_20261111` (part D, outside sources only)

Scope: generic web searches and OEIS queries only (no project detail in any request). Statements of other authors
are paraphrased, with location; one short verbatim quote is used (section 5). "Not found" means: not found by the
searches listed, which is weaker than "not in the literature". Items marked [memory] are bibliographic data I did not
see on a page today.

## 0. Method

- OEIS: `curl -s -A "Mozilla/5.0 (research; math-repo)" "https://oeis.org/search?q=...&fmt=json"`; a `null` body is a
  non-match. The mechanism was validated three times during the session with known sequences (A002785, A038375,
  A054946 all returned).
- Web: about 25 generic searches, plus reading of: Moon, *Topics on Tournaments* (Project Gutenberg edition), sections
  5-6 (pages 9-19); the arXiv abstract pages or HTML of the papers cited below; the slides of the IPAC talk on the
  Redei-Berge symmetric function.

## 1. Parity of the number of cyclic triangles

**Classical (found).**
- `c_3 = C(n,3) - sum_i C(s_i, 2)`: Moon, *Topics on Tournaments* (1968), section 5, corollary to Theorem 4 (p. 11 of
  the Gutenberg edition), attributed there to Kendall and Babington Smith (1940), Szele (1943), Clark; exercise 5 of
  that section gives the form `n(n^2-1)/24 - (1/2) sum (s_i - (n-1)/2)^2`. Also Proposition 1.1 of M. Burnham,
  "Structure of singular and nonsingular tournament matrices", arXiv:2210.12371 (2022), called well known there.
- Distribution of `c_3` over labelled tournaments: Moon section 6 (tables of Kendall-Babington Smith for `3 <= n <= 7`,
  Alway for `8 <= n <= 10`; asymptotic normality, Moran 1947; third and fourth moments). OEIS: A000142 (no triad),
  A090672 (one triad, index shifted by one), A357242, A357257, A357248, A357266 (exactly 2, 3, 4, 5 circular triads), from I. R. Harris and
  R. P. A. McShane, "Counting Tournaments with a Specified Number of Circular Triads", J. Integer Sequences 27 (2024),
  Article 24.8.7, and J. B. Kadane, Ann. Math. Statist. 37 (1966) 488-494. None of these OEIS entries, nor the landing
  page/abstract of the JIS paper, mentions parity. (Full text of the JIS paper not read.)
- Score parities as invariants: C. Waldrop Jr., "An arc-reversal theorem for tournaments", Congr. Numer. XVII (1976)
  501-507, as stated in W. H. W. Wong and E. G. Tay, arXiv:2108.12605 (Discrete Math. 2022), Theorem 4.2: two
  tournaments of the same order have the same number of even-score vertices if and only if one is obtained from the
  other by reversing transitive triples; Lemma 4.3 there: reversing a directed `k`-cycle preserves all score parities.
  J. Bang-Jensen, F. Havet, F. Hoersch, C. Rambaud, A. Reinald, C. Silva, arXiv:2511.22562 (2025), Theorem 1.3: for
  odd `p` a tournament on `n >= p + 2` vertices can be made acyclic by inversions of `p`-sets iff it has `ceil(n/2)`
  vertices of even out-degree (out-degree parities are invariant under such inversions).
- Switching classes of tournaments: L. Babai and P. J. Cameron, Electron. J. Combin. 7 (2000) #R38; OEIS A049313.
  The abstract does not mention "even tournaments", scores or parity (full text not read).

**Not found.** No source states "the parity of `c_3` depends only on the number of odd scores mod 4", the formula
`c_3 = C(q,2) + floor(n/4) (mod 2)`, the value `E (-1)^(c_3) = 2^(-floor((n-1)/2))`, or the count
`2^(m-1) + 2^(ceil(m/2)-1)` of one-vertex extensions with an even number of triangles. Moon sections 5-6 (text and
exercises) contain no parity statement. OEIS has no entry for the labelled counts with `c_3` even
(6, 48, 640, 20480, 1179648) or odd (2, 16, 384, 12288, 917504).

**How close to classical (my derivation, checked by a short script, not a quotation).** The parity law is a two-line
corollary of the classical formula: `sum_i C(s_i,2) = (sum s_i^2 - C(n,2))/2` and `sum s_i^2 = q (mod 4)`, so
`c_3 = C(n,3) + (q - C(n,2))/2 (mod 2)`, which equals `C(q,2) + floor(n/4)` for every `n <= 64` and admissible `q`.
With the (standard) uniformity of the score-parity vector on its coset and the binomial sum
`sum_{q = C(n,2) (2)} C(n,q) (-1)^C(q,2) = (-1)^floor(n/4) 2^floor(n/2)`, this gives `E (-1)^(c_3) = 2^(-floor((n-1)/2))`
(reproduced for `n <= 20`). The classical tables agree: even/odd counts 6/2, 48/16, 640/384, 20480/12288 at
`n = 3, ..., 6`.

## 2. Parity of the number of 5-cycles

**Found (related only).**
- N. Komarov and J. Mackey, "On the number of 5-cycles in a tournament", J. Graph Theory 86 (2017) 341-356
  (arXiv:1410.6828): exact integer formula for `c_5` in terms of edge scores; bounds. No parity statement in the
  abstract.
- S. V. Savchenko, "On 5-cycles and strong 5-subtournaments in a tournament of odd order n", arXiv:2403.19555 (2024):
  restates the Komarov-Mackey formula (`8 c_5 = 6 C(n,5) + sum over arcs of a function of the four intersection
  numbers`); states that the trace of the `m`-th power of the adjacency matrix counts closed `m`-walks and that for
  `m = 3, 4, 5` every such walk is a shift of an `m`-cycle, so `tr A^m = m c_m`. So **`tr A^5 = 5 c_5` is in the
  literature**. No parity or congruence statement for `c_5` or `c_3`; no use of `A^2 1`.
- `A^2 1` as an integer vector is the classical second-level score of Wei (1952) and Kendall (1955) (sum of the scores
  of the beaten players; Sonneborn-Berger), see e.g. arXiv:2210.17300. Nothing found on its reduction mod 2.

**Not found.** Any statement on `c_5 mod 2`, on `c_3 + c_5 mod 2`, or on the vector `A^2 1 mod 2` in connection with
cycle counts. OEIS has no entry for the labelled counts with `c_5` even (624, 17408, 1093632) or for the signed sums
(224, 2048, 90112, 11534336; 416, 8192, 548864, 70254592).

## 3. Rank and kernel of tournament matrices over GF(2)

**Found.**
- D. de Caen, "The ranks of tournament matrices", Amer. Math. Monthly 98 (1991) 829-831 [memory for the pages]: rank
  at least `(n-1)/2` over every field; restated for `Z_2` in R. Fulek and J. Kyncl, arXiv:1903.08637, section 4.1,
  equation (3), with the proof from `A + A^T = J + I` and subadditivity of the rank; also cited in arXiv:2108.10871.
  T. S. Michael, same title, Amer. Math. Monthly 102 (1995) 637-639. Doering, Michael, Shader, "Even and odd tournament
  matrices with minimum rank over finite fields", Electron. J. Linear Algebra 22 (2011) (search snippet only: for even
  `n` the rank over characteristic `p` is at least `n/2`, with a characterisation of equality through `A A^T`).
- Burnham, arXiv:2210.12371: number of 3-cycles against singularity over the reals (not mod 2).
- Coding theory parallel (standard; only the first half confirmed by a search today): for `p = -1 (mod 8)` the extended
  binary quadratic-residue code is self-dual and doubly even, hence the code of length `p` has weights `0` or
  `3 (mod 4)`.

**Not found.** The statement that kernel vectors of a tournament matrix over `F_2` have weight `0` or `1 (mod 4)` (it is
one line from `u^T A u = C(|u|,2)`), any link between the `F_2` kernel and counts of short cycles mod 2, or a Gauss-sum
evaluation of the kind in Theorem P (3).

## 4. OEIS

| query | result |
|---|---|
| 1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784 | **A002785**, "Number of self-complementary oriented graphs with n nodes"; a comment (McKay, 2020) identifies it with self-converse tournaments. Confirmed. |
| labelled self-converse: 8,48,704,7040; 704,7040,380416; 380416,6737920; 1,2,8,48,704(,7040); 2,8,48,704,7040; 48,704,7040,380416; 7040,...,958578688; single terms 6737920, 958578688, 31596429312 | no match (11 queries) |
| signed sums 1,2,4,32,416,8192; 4,32,416,8192(,218368); 416,8192,218368; 218368,11149312; 32,416,8192 | no match |
| 2, 4, 32, 416 | A225170 (fanout-free Boolean functions of AND rank 1): coincidence of four terms (next term 7552, not 8192) |
| odd-order self-converse signed sums 4, 256, 37824, 20782592 | no match |
| labelled counts with `H = 1 (mod 4)`: 6, 48, 720, 20480, 1157760; with `H = 3`: 2, 16, 304, 12288, 939392 | no match |
| signed strong counts 2,24,176,4080,123952,8135680 (also with the sign, and shorter windows) | no match |
| numerators 13,67,263,1007; 13,67,263; 853,1361,63103(,650471); 650471; 260546929; 1706,2722; 7,11,23,49 | no match |
| class race 8,36,264,1268,22600; 36,264,1268,22600 | no match |
| text: "self-converse tournaments" | only A002785 |
| 1, 0, 2, 24, 544, 22320 | A054946 (strongly connected labelled tournaments), used as control |

## 5. `H mod 4`, `H mod 8`, and the conjecture

- **`H = 1 + 2 c_odd (mod 4)` has a literature source**: D. Grinberg and R. P. Stanley, "The Redei-Berge symmetric
  function of a directed graph", arXiv:2307.05569 (2023); the abstract announces a modulo-4 refinement of Redei's
  theorem (section 7 of the paper). Talk slides of the first author (IPAC seminar, 2023-08-15), section 6: for a
  tournament, the number of Hamiltonian paths is `1 + 2 (number of odd cycles of length > 1) (mod 4)`, with the remark
  "We have not seen this anywhere in the literature" (Grinberg, slides).
- **The mod-8 formula used in Theorem S (2)** (`H = 1 + 2 alpha_1 + 4 alpha_2`) follows from their second main theorem
  (Theorem 1.38 of the paper; "Main Theorem II" of the slides): for a tournament the Redei-Berge function is the sum,
  over permutations all of whose cycles are odd and whose non-trivial cycles are cycles of the tournament, of
  `2^(number of non-trivial cycles)` times a power-sum monomial; specialising gives `H = sum_k 2^k alpha_k`.
- **Conjecture.** S. H. Chan and I. Pak, "Computational complexity of counting coincidences", arXiv:2308.10214,
  section 6.5, **Conjecture 6.2**: the function `h` takes all odd values except 7 and 21. The paper attributes it to
  the MathOverflow user bof and G. Royle (March 2016; checked up to 80557); the slides above refer to MathOverflow
  question 232751 (range given there as 80555).
- **Not found:** any work on the distribution of `H mod 4` or `mod 8` over all tournaments; any statement on `H mod 4`
  of self-converse tournaments. On self-converse tournaments: W. J. R. Eplett, "Self-converse tournaments", Canad.
  Math. Bull. 22 (1979) 23-27 (structure and score vectors); M. R. Sridharan, Indag. Math. 32 (1970) 441-447
  (unlabelled enumeration; linked from A002785); E. Thornblad, arXiv:1611.09579 (limits). No labelled count found.

## 6. Spot-checks of the external facts in section 8 of the note

| fact in the note | check |
|---|---|
| primes mod 4: first lead of class 1 at 26861 | confirmed (Wikipedia, "Chebyshev's bias") |
| logarithmic density 0.9959 | confirmed as the Rubinstein-Sarnak value, **conditional on GRH and the simplicity hypothesis**; the note states it without the condition |
| Polya: first failure 906,150,257 | confirmed (Tanaka 1980) |
| `pi(x) < li(x)` up to `10^19` | confirmed through secondary citations of Buthe |
| Mertens checked to `10^16` | confirmed (Hurst, arXiv:1610.08551) |
| Fermat numbers composite for `5 <= k <= 32` | confirmed; `F_33` open |
| square-sum paths: `n = 15, 16, 17, 23` and all `n >= 25` | confirmed (A090461; proof by Gerbicz, 2018) |
| Busch's minimum for strong tournaments | confirmed: A005517 (3, 5, 9, 15, 25, 45, 75 at `n = 3..9`; formula `3*5^((n-3)/3)`, `5^((n-1)/3)`, `9*5^((n-5)/3)`), comment citing A. H. Busch, Electron. J. Combin. 13 (2006) N3 |
| maximum `H`: 1, 3, 5, 15, 45, 189 at `n = 2..7` | confirmed: A038375 (1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095) |
| ball volume peak | confirmed: real dimension 5.2569464..., integer maximum at 5 (value 5.2638) |

## 7. Judgment for the auditor

1. **Theorem T and the first formula of Theorem P (1) are elementary corollaries of the Kendall-Babington Smith
   formula** (section 1 above). I found no source that states them, so "not found" is accurate, but the typing should
   say so: e.g. "elementary; immediate from `c_3 = C(n,3) - sum C(s_i,2)` since `sum s_i^2 = q_1 (mod 4)`; not found
   stated in the literature searched". The outcome line "(PROVED) The triangle race has a permanent, exact bias" and
   "NEW here: (T)" in THM-5131 are acceptable only with that qualifier ("new to the repository").
2. **`tr A^5 = 5 c_5`** is in the literature (Savchenko, arXiv:2403.19555), not only in the repository.
3. **The pentagon law** (`c_3 + c_5 = C(q_2,2) + floor(n/4)`, the cube sum, the closed forms): nothing found. The
   note's "no claim about the literature" is the right wording; the searches give no reason to change the typing.
4. **Weights mod 4 of kernel vectors**: one line from `A + A^T = J + I`; the rank bound of de Caen is the classical
   neighbour and could be cited next to Theorem P (3).
5. **`H = 1 + 2 c_odd (mod 4)` and `H = sum 2^k alpha_k`** have an outside source (Grinberg-Stanley 2023); the note
   credits repository files only. Worth checking that those files, or the note, name it.
6. **OEIS**: A002785 confirmed for the class counts; the labelled self-converse counts, the signed sums and the
   strong signed counts are absent.
7. **Section 8**: every external fact checked is right; the density 0.9959 should carry "under GRH and the simplicity
   hypothesis", and the conjecture's origin (bof and Royle, MathOverflow 2016) could be named beside Chan-Pak.

Not done: full-text reading of the JIS paper, of Babai-Cameron, of de Caen's and Michael's notes, of Komarov-Mackey;
MathSciNet/zbMATH searches (no access).
