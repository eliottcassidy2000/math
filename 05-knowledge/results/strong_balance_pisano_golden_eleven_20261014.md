# Strong tournaments, Pisano periods and the golden eleven: a majority-cycle route to 1/3–2/3, golden Wieferich primes, and an icosahedral analogy for the eleven-square torus

2026-10-14, session opus-2026-10-14-S36 (worktree `/tmp/math-wt-strong-balance-pisano-36`). Audited (section 9; MISTAKE-614).

**Owner's prompt (verbatim).** Attached: "The Fibonacci geometry of eleven squares" (PINGYOU LTD, 2026; UNREFEREED).

"try to prove the tournament 1/3-2/3 conjecture for strong tournaments. consider the idea that the square sum problem (arranging a list of integers such that the sums of each consecutive number is a square) and the graceful tree conjecture are each analogous roughly with the 3n-1 and normal collatz conjectures, respectively. try to figure out how to prove all square packing bounds simultaneously. consider the fibonacci series under all possible bases and the mapping of a modulus to a pisano period, how a modulus 11 maps to a length 10 and 19 to 18 and 12 to 24 and 20 to 60, and it feels like primes go one below them while composites go to increasingly composite numbers, and how the period of the modulus divided by the modulus maxes out at 6, where 10 maps to 10 and where 50 maps to 300 ... etc. and see how this relates deeply to path of least resistance recursion we have studied at the heart of number theory. (arXiv 2607.23926) reports 128 poset classes at order 14 with balance constant 1/3, matching Narayana's cows number c(14)−1 and linking to the 1/3–2/3 problem. consider the identities involving 11 coordinates and 55 pairs. consider the golden rectangle and the dodecahedron and icosahedron. think of how https://en.wikipedia.org/wiki/Wall%E2%80%93Sun%E2%80%93Sun_prime realte with our work regarding wieferich primes, think rigitdity circuits"

**Concurrent work.** mac-mini-2026-10-14-tourn worked on the same brief in parallel and reached several of the same results first.
- **THM-4655**: 1/3–2/3 for every strong transitive-plus-one tournament, plus a rotation lemma.
- **THM-4656**: Wall–Sun–Sun primes as golden Wieferich primes. This includes the orbit law reserved here as THM-4657, which is now ceded.
- **HYP-9279**: vertex-wise balance, and `δ >= 4/11`.
- Exhaustive balance checks for all strong `n <= 9`.

This note credits those results and keeps only what is this session's own.

**Status.**
- **PROVED:**
  - the reductions in HYP-9280: majority cycle ⇒ 1/3–2/3, and the ladder (SS) ⇒ `D_in > 0` ⇒ vertex lemma ⇒ no majority source ⇒ `M(T)` not transitive ⇒ a balanced pair;
  - `Pr[w first] < 1/2` for every non-source;
  - the vertex lemma for in-degree or out-degree `<= 2`;
  - the weak form when `T` has at most one inverted arc;
  - the exact identity `C(u, w) − C(w, u) = H(T − w − o)` at an out-degree-1 vertex `w → o` when `o` beats all other in-neighbours of `w` (section 11);
  - the Lucas–Fibonacci identities of section 3, several of them the eleven-squares paper's own (credited);
  - the exact verification of the Erdős Problem 106 counterexample.
- **FINITE-EXACT:**
  - HYP-9280 (`F(T)` strongly connected), the vertex lemma and the sharp `T_n` minima for every strong tournament with `n <= 9`. This session covered labelled `n <= 7` and the `n = 8` majority sweep; the independent audit covered `n = 8, 9` by isomorphism classes.
  - The same, together with "no majority source", for every strong tournament of feedback arc number `<= 2` (`n <= 16`), `<= 3` (`n <= 14`), `<= 4` (`n <= 12`) and `<= 5` (`n = 10`) (section 11).
  - (SS) for `n <= 8`, and its failure at `n = 9`;
  - the failure of the averaged form `D_in > 0` at `n = 11` (section 11);
  - the PSL(2, 11) structure;
  - Pisano data;
  - no Wall–Sun–Sun prime below `2·10⁶`;
  - the square-sum data to 40.
- **NUMERICAL:** annealing searches, `n = 8..12`.
- **OPEN:** HYP-9280 (which would prove the tournament 1/3–2/3 conjecture for strong tournaments); HYP-9281 (an interpretation).
- **CITED:**
  - Rédei (1934); Moon (1968);
  - Wall (1960); Sun–Sun (1992); Bragman–Rowland (2022); McIntosh–Roettger (2007); Dorais–Klyve (2011); PrimeGrid (2020–2022); Freyd–Brown (1992);
  - Galois (1832 letter); Klein (icosahedron); Liebeck–Praeger–Saxl (1990); Luo (1989);
  - Erdős–Graham (1975); Roth–Vaughan (1978); Chung–Graham (2009, 2020); Wang–Dong–Li (2016); Bui (2025); McClenagan (2026);
  - Gerbicz (2018; OEIS A090461);
  - the linear-extension-majority literature (section 1.6).
- **ANALOGY:** the icosahedral picture of the torus colouring; the packing and square-sum lenses; rigidity circuits.
- **UNREFEREED inputs:** the eleven-squares paper; Gupta, arXiv:2607.23926; Sprite143's GitHub counterexample.

## 0. Answers in brief

1. **The tournament 1/3–2/3 for strong tournaments is not proved here.** It reduces, by proved steps, to the majority-cycle conjecture HYP-9280, and further to purely local statements.
   - **The conjecture.** In a uniform random Hamiltonian path, the arcs whose tail usually comes first form a strongly connected spanning subgraph. The majority tournament is then strong, with at least `n − 2` cyclic triples. Each cyclic triple has its three precedence probabilities summing to at most 2, so it contains a pair in `(1/2, 2/3]`.
   - **The weakest local condition found.** It would suffice that no vertex of a strong tournament precedes every other vertex with probability `> 1/2` ("no majority source"). A transitive majority tournament has such a source. That condition is implied by the vertex lemma (some in-neighbour of `w` precedes `w` with probability `> 1/2`).
   - **Proved.**
     - The vertex lemma for in-degree `<= 2`, and dually for out-degree `<= 2`.
     - The first vertex of a random path is a given non-source with probability `< 1/2`.
     - With at most one inverted arc the forward digraph has a cycle.
   - **Evidence.** The conjecture holds for every strong tournament with `n <= 9`: this session to `n = 7` (and the majority statements at `n = 8`), the independent audit at `n = 8, 9`. Annealing to `n = 12` found no counterexample, and inverted arcs are rare: at most 3 at `n = 9`.
   - **Sharp constants.** The transitive tournament with its extreme arc reversed is the unique extremiser for `n <= 9` of the 3-cycle certificate, the vertex-lemma value and `D_in`, with value `(2^(n−3) + 1)/(2^(n−2) + 1)` or margin `1/(2^(n−2) + 1)`. It is not the extremiser of the balance constant itself (`4/11`, `62/157`).
   - **What died.** Two tempting routes. The stochastic dominance (SS) is false at `n = 9`. The averaged form "`E[N] > d/2` at every vertex" is false at `n = 11`, at an out-degree-1 vertex of a "double chain"; there the vertex lemma survives.
   - **What is left.** The weakest open rung is (SS) at `k = 1`. No failure is known, and it would give the vertex lemma for in-degree and out-degree `<= 4` (section 11).
2. **Pisano periods.**
   - `π(n) <= 6n`, with equality exactly at `n = 2·5^k`, `k >= 1` (10, 50, 250, 1250, …). Your "10 maps to 10" should read `π(10) = 60`.
   - Primes split by their behaviour in `Q(√5)`:
     - `p ≡ ±1 (mod 5)` have `π(p) | p − 1`, the "one below" (`11 → 10`, `19 → 18`, but `29 → 14`);
     - `p ≡ ±2 (mod 5)`, `p >= 3`, have `π(p) | 2(p + 1)` but not `p + 1` (`7 → 16`, `47 → 32`);
     - `π(2) = 3` and `π(5) = 20`.
   - Composites take least common multiples of their prime-power periods. Here `π(p^k) = p^(k−1)π(p)` whenever `p` is not Wall–Sun–Sun, and no Wall–Sun–Sun prime is known.
3. **Wall–Sun–Sun versus Wieferich.**
   - Wall–Sun–Sun primes are exactly the Wieferich primes to the base `φ`. This is known (Bragman–Rowland 2022; "Fibonacci–Wieferich primes").
   - The tournament form is mac-mini's THM-4656, which this session found independently (THM-4657, ceded). The circulant tournaments on `Z_(p²)` invariant under the odd part of `φ` explode exactly at Wall–Sun–Sun primes, as the doubling-invariant ones explode at Wieferich primes (THM-4644).
   - None exists below `2·10⁶` (this run), below `9.7·10^14` (Dorais–Klyve), or below `2^64` (PrimeGrid).
   - At the golden Paley primes `11, 19, 31, 59, …` the invariant circulants are exactly the Paley tournaments.
4. **11 coordinates and 55 pairs.**
   - `11 = L_5 = 2F_5 + 1`, and `55 = F_10 = F_5 L_5 = C(11, 2) = |Aut P_11|`.
   - `C(L_n, 2) = F_(2n)` exactly for `n = 2, 4, 5`, which gives Luo's triangular Fibonacci numbers 3, 21, 55 at `L_n = 3, 7, 11`.
   - The eleven-squares paper's "exactly when `n = 5`" (its Theorem 6.4) is `|R_n| = 2n + 1`, which needs `L_n = 2n + 1`.
   - That all of this is one "golden fixed point `F_5 = 5`" structure is an interpretation (HYP-9281). mac-mini types `55 = F_10` as a coincidence.
5. **The golden rectangle and the icosahedron.**
   - Inside `PSL(2, 11)` (order 660), an icosahedral `A_5` and `Aut(P_11)` (order 55) multiply to the whole group, `PSL(2, 11) = A_5·Aut(P_11)`, with intersection `C_5`. This is classical; it is Galois' 11-point action.
   - On the 12 points of `P¹(F_11)` the `A_5`-invariant valency-5 graph is the icosahedron, with `∞` and `0` antipodal and pentagons the squares `{1, 3, 4, 5, 9}` and the non-squares.
   - The eleven-square torus colours its cells by exactly this pattern `{0}, H, −H`. Beyond that orbit pattern nothing matches: the torus adjacency is 6-regular, and `A_5` does not act on the cells. So "the torus is an icosahedron with a marked vertex" is an ANALOGY.
6. **Square packing.**
   - The paper's (and Sprite143's) counterexample to Erdős Problem 106 (`f(17) > 4`) checks out in exact rational arithmetic. erdosproblems.com already lists the problem as disproved.
   - "All square-packing bounds simultaneously" is open. Exact `s(n)` is known only sporadically. The wasted-space exponent lies between `1/2` (Roth–Vaughan) and `3/5` (McClenagan 2026; Bui).
   - A unifying lens is proposed in section 4, typed ANALOGY.
7. **Square-sum ~ 3n−1, graceful ~ 3n+1.**
   - Square-sum orderings exist exactly for `n ∈ {15, 16, 17, 23}` and all `n >= 25` (Gerbicz 2018; checked here to 40).
   - The analogy pairs a finite exceptional set, like the extra `3n−1` cycles, with no known exception. It is typed ANALOGY.

## 1. Strong tournaments: the majority-cycle route (HYP-9280)

**Notation.**
- `p(x, y) = Pr[x before y]` in a uniform random Hamiltonian path. Since `H` is odd (Rédei), `p ≠ 1/2`.
- `M(T)` is the majority tournament.
- An arc of `T` is forward or inverted as `p(tail, head)` is `>` or `<` `1/2`.
- `F(T)` is the forward digraph.

**1.1 Reductions (PROVED).**
- Every linear order satisfies one or two of the three cyclic precedences of a triple. So a directed 3-cycle of `M(T)` has `Σ p <= 2` with each `p > 1/2`, and one of them lies in `(1/2, 2/3]`: a balanced pair.
- `F(T) ⊆ M(T)`. If `F(T)` is strongly connected and spanning, `M(T)` is strong and has at least `n − 2` directed 3-cycles (Moon).
- **The ladder.** At every vertex `w`, let `d = indeg(w)` and let `N` be the number of in-neighbours before `w`. Then

  (SS) `Pr[N >= d − k] >= Pr[N <= k]` for all `k` ⇒ `D_in(w) = Σ_(u→w)(2p(u,w) − 1) > 0` ⇒ the vertex lemma (some in-neighbour precedes `w` with probability `> 1/2`) ⇒ no vertex is a majority source ⇒ `M(T)` is not transitive ⇒ a balanced pair.
  - The step from `D_in` to the vertex lemma uses that `E[N] >= d/2` forces some `p(u, w) >= 1/2`, hence `> 1/2` by parity.
  - The first two rungs are false in general: (SS) at `n = 9`, `D_in > 0` at `n = 11` (section 11). The implications still hold, and the live rungs start at the vertex lemma.
  - Parity does not make `E[N] > d/2` automatic for even `d`.
  - The fourth implication holds because a strong `T` has no source.
  - "No majority source" has to hold only at one vertex, the source of `M(T)` if there is one. It allows any precedent, not just an in-neighbour. At `T_n` it is exactly as tight as the vertex lemma.
- Without a balanced pair, the strong-majority relation `p > 2/3` is a linear order (HYP-9277). Then every 3-cycle of `T` has an arc with `p < 1/3`.

**1.2 Evidence (FINITE-EXACT, `n <= 9`; NUMERICAL to `n = 12`).**

| `n` | strong labelled tournaments | `M(T)` transitive | no forward 3-cycle in `T` | inverted-arc distribution (labelled) | source |
|---|---|---|---|---|---|
| 3–5 | 2, 24, 544 | 0 | 0 | always 0 | session |
| 6 | 22,320 | 0 | 0 | 0: 20,880; 1: 1,440 | session, audit |
| 7 | 1,677,488 | 0 | 0 | 0: 1,633,808; 1: 43,680 | session, audit |
| 8 | 236,522,496 | 0 | 0 | 0: 231,227,136; 1: 4,932,480; 2: 362,880 | session, audit |
| 9 | 64,026,088,576 | 0 | 0 | 0: 63,298,691,584; 1: 693,286,272; 2: 33,626,880; 3: 483,840 | audit |

- `F(T)` is strongly connected in all of them. The vertex lemma and `D_in > 0` hold at every vertex, and the fewest 3-cycles of `M(T)` is `n − 2`.
  - For `n <= 7` this is this session's `tournament_inverted_arcs_structure_20261014.py` and `tournament_vertex_lemma_20261014.py`.
  - For `n = 8, 9` it is the audit's independent isomorphism-class sweep. This session's own `n = 8` sweep tested only the majority statements.
- Inverted arcs `x → y` at `n <= 7` have degree type `(outdeg y, indeg x) = (1, 1)`, a "squeezed triangle" `y → z → x`. At `n = 7` the types `(2, 4)` and `(4, 2)` also occur, 10,080 times each.
- **Annealing (`n = 8..12`)** found at most 3 inverted arcs. The minimum number of 3-cycles of `M(T)` was `n − 2` (Moon's bound) for every `n`, and `F` was strongly connected at every optimum.
- The HYP-9277 Claim A counterexamples have one inverted arc each and 10 or 12 cyclic majority triples.

**1.3 Sharp extremal family (FINITE-EXACT, `n <= 9`).** `T_n` is transitive with the extreme arc reversed, with `H(T_n) = 2^(n−2) + 1`. For every `n <= 9` the following extremes over strong `n`-tournaments are attained only at `T_n`:
- the best 3-cycle certificate (max over 3-cycles of the smallest arc probability): minimum `(2^(n−3) + 1)/(2^(n−2) + 1)`;
- the vertex-lemma value `min_w max_(u→w) p(u, w)`: the same value;
- `min_w D_in(w)` and `min_w D_out(w)`: `1/H(T_n)`;
- `max_w Pr[w first]`: `2^(n−3)/(2^(n−2) + 1) < 1/2`.

Annealing to `n = 12`, and the exhaustive near-transitive sweeps of section 11, find the same minima, at `T_n`.

The balance constant `δ(T) = max over pairs of min(p, 1 − p)` is different. Its minimum over strong tournaments is `4/11` (`n = 5, 6, 7`) and `62/157` (`n = 8, 9`), attained by tournaments other than `T_n`; mac-mini's HYP-9279(A) conjectures `δ >= 4/11`.

**1.4 Proved pieces** (the injections are spelled out in HYP-9280).
- **(SS) at `k = 0`.** Move `w` from the front of a path to just after its last in-neighbour. This is a valid injection into paths with `N = d`.
- **First-vertex bound.** `Pr[w first] < 1/2` for every non-source, by the same injection and parity.
- **Vertex lemma for in-degree `<= 2`** (and out-degree `<= 2`). For `d = 2`, `E[N] = Pr[N >= 1] + Pr[N = 2] >= 1`. So in an acyclic `F(T)`, and at a majority source, the vertex has in-degree `>= 3`.
- **Weak form.** If `T` has at most one inverted arc, `F(T)` contains a directed cycle. This yields a balanced pair, but not the strong connectivity of `F`.

**1.5 What fails.**
- "Every arc has `p(tail first) >= 1/3`" is false at `n = 8` (HYP-9277).
- **The stochastic symmetry (SS)** holds for all strong `n <= 8`. It fails at `n = 9` in exactly 5 isomorphism classes (1,572,480 labelled tournaments).
  - Every failure sits at a vertex of in-degree 7 (out-degree 1), with `k = 3`.
  - In each, `E[N] > d/2` and the vertex lemma still holds. This session's sampled examples are two of the five.
- "Best in-neighbour precedence `+ Pr[w first] >= 1` per vertex" is false: 10,420,928 violating (tournament, vertex) pairs at `n = 7`.
- The pointwise complement pairing behind (SS) is false: 2,427 violations in 16,426 random tests.
- The averaged form `D_in(w) > 0` is false at `n = 11` (section 11).
- The full conjecture remains OPEN.

**1.6 Literature.**
- For posets, the linear-extension-majority (LEM) relation and its cycles are classical: Kislitsyn 1968; Fishburn 1974, 1976; Ewacha–Fishburn–Gehrlein 1990; De Loof–De Baets–De Meyer.
- Gupta (arXiv:2607.23926) finds LEM cycles of length at most 8 at order 14. The owner's "128 classes at order 14 = `c(14) − 1`" is S35's Narayana equality (HYP-9277).
- No prior tournament version of the majority relation was found.
- Posets can have a transitive LEM relation; HYP-9280 says strong tournaments cannot.

## 2. Pisano periods, Wall–Sun–Sun and Wieferich (THM-4656)

**2.1 Data (FINITE-EXACT; facts CITED).**
- `π(n)/n <= 6` up to 2000, with equality at `n = 10, 50, 250, 1250` (`= 2·5^k`, `k >= 1`). The bound is Freyd–Brown (1992) / K. S. Brown (OEIS A001175).
- Sample values:

  | `n` | 10 | 11 | 12 | 19 | 20 | 50 | 55 | 89 | 121 | 125 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | `π(n)` | 60 | 10 | 24 | 18 | 60 | 300 | 20 | 44 | 110 | 500 |

- Split primes have `π(p) | p − 1`. Inert primes `p >= 3` have `π(p) | 2(p + 1)` but not `p + 1` (checked for every prime from 7 to 2000; `p = 3` gives 8). `π(2) = 3` does divide `p + 1`, and `π(5) = 20`.
- The ratio `6 = (3/2)·4` comes from `π(2) = 3` and `π(5^k) = 4·5^k`.

**2.2 The orbit law (CITED and THM-4656; independently re-proved here as THM-4657, ceded).**
- For split `p` and a `p`-adic golden ratio `φ ∈ Z/p²`: Wall–Sun–Sun ⟺ `φ^(p−1) ≡ 1 (mod p²)`. This is known: Bragman–Rowland 2022, Theorem 10.
- Let `g` be the odd part of `φ`. The `g`-invariant circulant tournaments on `Z_(p²)` number `2^((p−1)/d')`, or `2^((p+1)(p−1)/(2d'))` exactly at Wall–Sun–Sun primes (THM-4656 part 1).
- The formula is checked by enumeration at all 40 (prime, root) pairs with split `p <= 199`. The inert case is THM-4656 part 2.
- No Wall–Sun–Sun prime exists below `2·10⁶` (this run), below `9.7·10^14` (Dorais–Klyve), or below `2^64` (PrimeGrid, 2020–2022).

**2.3 Golden Paley primes (THM-4656 part 4).** For `p ≡ 3 (mod 4)` and `p ≡ ±1 (mod 5)` where the square golden root generates the squares, the invariant circulants on `Z_p` are exactly `P_p` and its reverse.
- They are `11, 19, 31, 59, 71, 79, 131, …` (52 primes below 2000).
- `p = 11` is the eleven-square torus group: `G_5 = ⟨x + 1, 4x⟩ = Aut(P_11)`.
- There `ψ = −1/φ` is a primitive root.

**2.4 Least resistance (ANALOGY with exact anchors).**
- The Fibonacci matrix is the smallest hyperbolic element of `GL_2(Z)`. Its square, with trace 3, gives the shortest closed geodesic of the modular surface, of length `4 log φ`. This is the eleven-squares paper's "shortest closed modular geodesic".
- The Pisano period is the order of this least-resistance dynamics modulo `n`. It records the splitting of `p` in the golden field: the Frobenius sits in `Gal(Q(√5)/Q)`.
- The Wieferich and Wall–Sun–Sun conditions are the same failure of lifting from level `p` to level `p²`, in the bases 2 and `φ`. That failure is the level-2 jump of THM-4631 and THM-4644.
- THM-4656 part 3 adds the mirror image: golden rigidity happens exactly when `π(p)` is a power of 2, which includes the Mersenne primes `2^q − 1` with `q ≡ 3 (mod 4)`.

## 3. Eleven coordinates, fifty-five pairs, the golden rectangle and the icosahedron

**3.1 The golden numerology made exact** (PROVED, elementary; items marked [paper] are the eleven-squares paper's Lemma 5.1, Table 2 and Theorem 6.4, re-derived here).
- [paper] `N(φ^n − 1) = (−1)^n + 1 − L_n`, so `R_n = Z[φ]/(φ^n − 1)` has order `|L_n − 1 − (−1)^n|`: `1, 1, 4, 5, 11, 16, 29, 45, 76, 121, …`.
  - That is `L_n` for odd `n` and `L_n − 2` for even `n`.
  - `R_5 ≅ F_11` because `L_5 = 11` is prime. Also `R_4 ≅ F_5` and `R_2 = 0`.
- [paper] For `n >= 3`, `|G_n| = n|R_n|` for `G_n = ⟨x + 1, φx⟩`. Regularity on unordered pairs needs `|R_n| = 2n + 1`.
  - For odd `n` that is `L_n = 2n + 1`, true only at `n = 5`.
  - For even `n` it is `L_n = 2n + 3`, never true.
  - The pre-audit formula `|G_n| = nL_n` was wrong for even `n`.
- `C(L_n, 2) = F_(2n)` iff `L_n = 2F_n + 1` iff `F_(n−3) = 1` iff `n ∈ {2, 4, 5}`.
  - The cases are `C(3,2) = 3 = F_4`, `C(7,2) = 21 = F_8` and `C(11,2) = 55 = F_10`.
  - These are three of Luo's four triangular Fibonacci numbers (the fourth is 1).
  - `P_3, P_7, P_11` maximise the number of Hamiltonian paths at `n = 3, 7, 11` (A038375, known for `n <= 11`). That these are exactly the `L_n` with `C(L_n, 2) = F_(2n)` is NUMEROLOGY within the known range. Paley tournaments need not stay maximal: THM-135 shows `P_19` is not even the best circulant on `Z_19`.
- `F_n = n` only for `n ∈ {0, 1, 5}`. So at `n = 5` the two conditions `L_n = 2n + 1` and `L_n = 2F_n + 1` coincide.
- **Interpretation (HYP-9281, not proved).** The eleven-square coincidences `11 = L_5`, `55 = F_10 = F_5 L_5 = C(11, 2)` and `π(11) = 10` may be one structure at the Fibonacci fixed point `F_5 = 5`. mac-mini's THM-4656 types `55 = F_10` as a coincidence (by Luo's theorem, 55 is the largest triangular Fibonacci number). No mechanism is proved.

**3.2 Inside `PSL(2, 11)` (FINITE-EXACT; classical).**
- `|PSL(2, 11)| = 660`.
- The Borel subgroup `B = {x ↦ sx + t : s a nonzero square}` is `Aut(P_11)`, of order 55, and acts regularly on the 55 arcs.
- There are two non-conjugate subgroups `A ≅ A_5` containing the 5-fold rotation `x ↦ 4x`. Each is transitive on `P¹(F_11)` with stabilisers `C_5`, and `A ∩ B = C_5`.
- `PSL(2, 11) = A·B` as a product of subsets (`60·55/5 = 660`). It is not an exact factorisation, since `A ∩ B ≠ 1`; the exact one is `A_5·C_11`. This is Galois' 11-point action, and factorisations with an `A_5` factor are tabulated by Liebeck–Praeger–Saxl (1990).
- `A`'s orbitals on ordered pairs have sizes 12, 60 and 60. The two valency-5 symmetric orbitals are the icosahedron and its distance-2 graph.
  - In the icosahedron, `∞`'s neighbours are `{2, 6, 7, 8, 10}` and `0`'s are `{1, 3, 4, 5, 9}`, or the reverse in the other orbital.
  - `∞` and `0` are antipodal, and the two pentagons are the non-squares and the squares.
- The 12 vertices carry 15 golden rectangles, one for each two of the 6 antipodal pairs, forming 5 orthogonal triples. The dodecahedron is the dual, with the 20 triangles as its vertices.
- **The torus colouring (ANALOGY).** The paper colours `H = {1, 3, 4, 5, 9}` gold, `−H` teal and `0` grey. That is the `C_5`-orbit pattern of `Z_11` under `x ↦ 4x`, the pattern of the icosahedron seen from the vertex `∞`. Nothing more matches:
  - the paper's "mark" is a boundary point of the flat torus (the canonical lift encoding the side), not a twelfth cell;
  - the cell adjacency is the 6-regular circulant `C_11(1, 3, 4)`, not the icosahedron minus a vertex (degrees 4 and 5);
  - `A_5` does not act on the cells.
- **Galois' exceptional actions** (CITED). `PSL(2, p)` has a subgroup of index `p` only for `p <= 11`: `p = 2, 3` (trivially), `5` (tetrahedral `A_4`), `7` (octahedral `S_4`) and `11` (icosahedral `A_5`).
  - The Paley primes among them, 3, 7 and 11, carry the 3-cycle, the Fano plane and the 11-point biplane.
  - Associating the golden rectangle with 11 and the Klein quartic with 7 is an association, not an inference. It is NUMEROLOGY beyond the shared group theory.

## 4. Square packing

**4.1 Erdős Problem 106 (exact verification).** Credit: erdosproblems.com lists the disproof as due to Claude Opus 5 (prompted by Silverstein), linking Sprite143's repository. The rational rotation `(143/145, 24/145)` is the eleven-squares paper's refinement, and the paper cites Baek–Koizumi–Ueoro for the axis-parallel case. The 16 axis-parallel squares and the square rotated by `(143/145, 24/145)` (a 24-143-145 Pythagorean rotation) were checked in exact rational arithmetic.
- They lie in the unit square.
- Their interiors are pairwise disjoint, checked by the separating-axis test.
- Their side lengths sum to `4 + 101/1153304`.

So `f(17) > 4`. The grid extension (the 17-square packing scaled into a `4/k` corner of the `k × k` grid) gives `f(k² + 1) > k` for all `k >= 4`. The audit confirmed the check by an independent clipping-area method.

**4.2 "All bounds simultaneously" (status + a lens; ANALOGY).**
- **Exact values.** `s(n)` is known for small `n` and sporadic families.
  - For `s(11) <= 3.8770835900…` (Trump), the paper shows this is the unique minimum of a two-parameter contact family (its Theorem A).
  - An AI-assisted proof of global optimality, dated 29 Sep 2026, has been reported and is under review.
- **Asymptotics.** For the wasted space `W(x)` in packing unit squares into a square of side `x`:
  - upper bounds: `x^(7/11)` (Erdős–Graham 1975), `x^((3+√2)/7) log x` (Chung–Graham 2009) and `x^(5/8)` (Shuang Wang, T. Dong and J. Li 2016);
  - `O(x^(3/5))` is proved by McClenagan (arXiv:2602.01484, 2026) and, by a different method, by Bui (arXiv:2508.04603). Chung–Graham's 2020 argument for `3/5` contains an error, per McClenagan;
  - the lower bound is `W(x) > 10^(−100)·√(x‖x‖)` whenever `x(x − ⌊x⌋) > 1/6` (Roth–Vaughan 1978).

  So the exponent lies in `[1/2, 3/5]`; the true value is open.
- **The proposed lens.**
  - An optimal packing is a contact framework with inequality constraints, i.e. a tensegrity. Local optimality is certified by a self-stress, an LP dual.
  - In S35's language the certified configurations are rigidity circuits of the contact graph: minimally over-braced, one stress.
  - The eleven-square case realises the circuit as a Fibonacci torus, with contact matrix `I + wQ`.
- **A "simultaneous" programme** would classify circuit tori by their integral contact matrices `M ∈ GL_2(Z)` (hyperbolic, trace ordered) and compute each tilt from `M`'s characteristic equation. This is the analogue of reading every packing side off a Markov/geodesic spectrum, with golden first.
- This is a research programme, not a result.

## 5. Square-sum and graceful versus `3n ∓ 1` (ANALOGY)

**Square-sum (KNOWN; checked here to 40).**
- Orderings of `1..n` in which every two consecutive numbers sum to a square exist exactly for `n = 15, 16, 17, 23` and all `n >= 25`; cycles exist for all `n >= 32` (Gerbicz 2018; OEIS A090461).
- Our exhaustive check to 40 agrees: the exceptions are `2–14`, `18–22` and `24`.

**The analogy, read structurally.**
- A finite exceptional set, then universal success, matches `3x − 1` with its extra small cycles.
- The graceful tree conjecture (no known exception; verified for small trees) matches `3x + 1` with one cycle.
- Both are Hamiltonicity or labelling problems for arithmetic graphs. Collatz is a functional-graph problem, as in S30's carriers dictionary.

No transfer of proofs is claimed.

## 6. Rigidity circuits

- S35's HYP-9278 put the arithmetic of Paley pattern sums at 2D rigidity circuits. There `|E| = 2|V| − 2`, the double cover is Calabi–Yau, and the motive weight is `|V| − 2`.
- In square packing the circuits are the certified contact frameworks of section 4.2.
- In tournaments the analogue of "least resistance versus circuit" is the majority-cycle structure of section 1. The forward arcs are the low-resistance channels, and a single cyclic triple carries the balance.

This is ANALOGY, recorded as a possible unification.

## 7. Hostile controls

1. HYP-9280 is not proved. Its evidence is exhaustive to `n = 9` (with `n = 8, 9` from the audit's independent code), and annealing to `n = 12` is a heuristic search.
2. The sharp constants are FINITE-EXACT to `n = 9`. The formula `(2^(n−3)+1)/(2^(n−2)+1)` is conjectural beyond that, and the balance constant `δ` is not extremised by `T_n`.
3. The orbit law belongs to mac-mini's THM-4656, and Part 1 of it is in the literature. This session's THM-4657 was a duplicate and is ceded. No Wall–Sun–Sun prime is known, so the explosion branch is vacuous so far.
4. The icosahedral model of `P¹(F_11)` and the factorisation `PSL(2,11) = A_5·Aut(P_11)` are classical. Their link to the eleven-square torus is only the `C_5`-orbit pattern (ANALOGY).
5. The Erdős 106 counterexample is not ours. We give an independent exact check.
6. Sections 4.2, 5 and 6 are ANALOGY.

## 8. Open threads

1. Prove "no majority source", or the vertex lemma, for in-degree `>= 3`. That would prove the tournament 1/3–2/3 conjecture for strong tournaments. (SS) and `D_in > 0` are false. The next provable-looking rung is (SS) at `k = 1`, which would give in-degree `<= 4`; (SS) at `k = 2` already fails at `n = 10`.
2. Prove the out-degree-1 average: the in-neighbours that beat `w`'s out-neighbour average `p > 1/2`. Its `U− = ∅` case is proved in section 11.
3. Push the HYP-9280 check to all of `n = 10`; this also decides whether `D_in > 0` first fails at `n = 10` or at `n = 11`. mac-mini has an exact `n <= 10` DP (tbal2.c) and offered to run the majority statements on it.
4. Classify the contact-matrix tori behind the other known optimal square packings (`n = 5, 10, 19, …`).

## 9. Audit

An independent adversarial audit ran with its own code (`04-computation/experiments/s36_audit_20261014/`; report `REPORT.md` there).
- **Confirmed:**
  - the HYP-9280 reductions, the in-degree-1 lemma, (SS) at `k = 0`, and the vertex lemma for in-degree and out-degree `<= 2`;
  - all exhaustive tournament data for `n <= 7`;
  - the `T_n` family, which it extended to `n = 8, 9` together with HYP-9280 itself;
  - the two `n = 9` (SS) counterexamples, by brute force over `9!` orders;
  - the Claim A data;
  - the orbit-law proof and numerics, and the Pisano facts;
  - the identities of section 3.1 (except `|G_n| = nL_n`) and the PSL(2,11) facts;
  - the Erdős 106 check, by an independent clipping-area method;
  - the square-sum data.
- **Errors fixed in this version:**
  - THM-4657 duplicated mac-mini's earlier THM-4656. It is now ceded and cited, and Part 1 is typed CITED.
  - The PrimeGrid bound was stale; it is `2^64`.
  - `|G_n| = nL_n` was false for even `n`.
  - The torus-icosahedron statement was an overclaim; it is retyped as ANALOGY.
  - "Every sharp constant is attained by `T_n`" was false for `δ`.
  - The one-inverted-arc case gives only the weak form.
  - "`E[N] > d/2` by parity" was wrong for even `d`.
  - Several claims had no shipped script; they are now shipped (`tournament_inverted_arcs_structure_20261014.py`, `tournament_pointwise_pairing_test_20261014.py`, the annealing output).
  - Packing citations: the Chung–Graham 2020 error, McClenagan and Bui, Roth–Vaughan's actual bound, and Shuang Wang's name.
  - The inert Pisano law fails at `p = 2`, and the `6n` equality needs `k >= 1`.
  - The A038375 range is known only to `n = 11`; the Lucas match is NUMEROLOGY.
- **Imprecisions fixed:**
  - the first-vertex proof wording;
  - "exact factorisation" misused;
  - 15 golden rectangles, not three;
  - Galois' index-`p` statement;
  - typing in section 3.1;
  - the HYP-9281 predictions confused `L_n` with `|R_n|`;
  - credits to Gerbicz, Freyd–Brown and erdosproblems.com;
  - the meaning of "golden-invariant";
  - the `s(11)` status;
  - the LEM literature;
  - the cross-references to mac-mini's THM-4655 and HYP-9279;
  - the ranges, updated to `n <= 9`.
- **Ledger:** MISTAKE-614.

## 10. Reproduction and references

**Scripts** (`04-computation/experiments/`, outputs in `05-knowledge/results/`):
- `tournament_strong_balance_rules_20261014.py`;
- `tournament_majority_cycle_20261014.py`;
- `tournament_inverted_arcs_structure_20261014.py`;
- `tournament_majority_search_20261014.py`;
- `tournament_vertex_lemma_20261014.py`;
- `tournament_inneighbour_symmetry_20261014.py`;
- `tournament_ss_counterexamples_n9_20261014.py`;
- `tournament_vertex_lemma_search_20261014.py`;
- `tournament_pointwise_pairing_test_20261014.py`;
- post-audit: `tournament_near_transitive_sweep_20261014.py`, `tournament_double_chain_dfs_check_20261014.py`, `tournament_outdegree_one_vertex_lemma_20261014.py`, `tournament_witness_rules_20261014.py`, `tournament_ss_k1_refinement_test_20261014.py`;
- `fibonacci_pisano_wall_sun_sun_20261014.py`;
- `eleven_icosahedron_paley_psl211_20261014.py`;
- `erdos106_seventeen_squares_exact_check_20261014.py`;
- `square_sum_paths_20261014.py`;
- audit: `s36_audit_20261014/` (independent code, `.out` files, `REPORT.md`).

**Literature and sources:**
- Tournaments and posets:
  - L. Rédei (1934);
  - J. W. Moon, *Topics on Tournaments* (1968);
  - S. S. Kislitsyn (1968);
  - P. C. Fishburn, J. Combin. Theory B (1974, 1976);
  - K. Ewacha, P. C. Fishburn and W. V. Gehrlein, *Order* (1990);
  - K. De Loof, B. De Baets and H. De Meyer (mutual rank probabilities);
  - Gupta, arXiv:2607.23926 (UNREFEREED).
- Fibonacci, Pisano and Wall–Sun–Sun:
  - D. D. Wall, Fibonacci series modulo m, *Amer. Math. Monthly* 67 (1960);
  - Z.-H. Sun and Z.-W. Sun, *Acta Arith.* 60 (1992);
  - R. J. McIntosh and E. L. Roettger, *Math. Comp.* 76 (2007);
  - F. G. Dorais and D. Klyve, *J. Integer Seq.* 14 (2011);
  - B. Bragman and E. Rowland, arXiv:2202.00704 (2022);
  - PrimeGrid Wieferich/Wall–Sun–Sun search (2020–2022);
  - P. Freyd and K. S. Brown, *Amer. Math. Monthly* 99 (1992);
  - OEIS A001175, A244801, A000930, A038375, A090461.
- Groups and geometry:
  - E. Galois, letter to Chevalier (1832);
  - F. Klein, *Lectures on the Icosahedron* (1884);
  - M. W. Liebeck, C. E. Praeger and J. Saxl, *Mem. AMS* 432 (1990);
  - M. Luo, *Fibonacci Quart.* 27 (1989).
- Square packing:
  - P. Erdős and R. L. Graham, J. Combin. Theory A 19 (1975);
  - K. F. Roth and R. C. Vaughan, J. Combin. Theory A 24 (1978);
  - F. Chung and R. Graham (2009, 2020);
  - Shuang Wang, T. Dong and J. Li (arXiv:1603.02368);
  - Bui, arXiv:2508.04603;
  - R. McClenagan, arXiv:2602.01484 (2026);
  - PINGYOU LTD, "The Fibonacci geometry of eleven squares" (2026, UNREFEREED);
  - erdosproblems.com/106; github.com/Sprite143/erdos-106-counterexample;
  - Wikipedia, "Square packing".
- Wikipedia, "Wall–Sun–Sun prime".

**This repository:**
- mac-mini's THM-4655, THM-4656 and HYP-9279;
- HYP-9277, HYP-9278, THM-4644, THM-4631, THM-135, THM-4592;
- the S16 eleven-square octic note and the S35 notes.

## 11. After the audit: the averaged form fails, wider evidence, out-degree-1 vertices

This continues the proof attempt after the audit fixes. Details are in HYP-9280, "Post-audit continuation". A second independent audit checked this continuation (11.5; MISTAKE-617).

**11.1 `D_in > 0` is false (FINITE-EXACT, two independent methods).**
- **The family.** The double chain `D(n, j)` is the transitive order on `0..n−1` with the arcs `0 → j` and `j → n−1` reversed.
- **The failing vertex.** At its out-degree-1 vertex `w = 0`, the top vertex is an in-neighbour that precedes `w` only rarely: `p = 33/289` in `D(11, 5)`. That drags `E[N] − d/2` down to `−23/578`.
- **The vertex lemma survives.** The in-neighbours `6..9` have `p = 161/289`.
- **Verification.** Explicit enumeration of all 289 Hamiltonian paths confirms the values, in agreement with the subset DP.
- **Growth.** Along the family, `min D_in` reaches `−0.864` at `n = 18`.
- **Where.** Every failure found is at a vertex of in-degree `>= 9`. For `d <= 6`, `E[N] − d/2 > 0` in every sweep: at least `0.0471` for `3 <= d <= 6`, tending to 0 for the proved degrees `d = 1, 2`.
- **The annealing caveat.** The earlier annealing had reported `min D_in = 1/(2^(n−2) + 1)` at `n = 8..12`. At `n = 8, 9` that is the exact minimum and at `n = 10` nothing refutes it. At `n = 11, 12` it was a search artefact.

**11.2 Exhaustive near-transitive sweep (FINITE-EXACT).** Every extremal and failing example so far is nearly transitive. The sweep therefore covered every strong tournament of feedback arc number:
- `<= 2` for `n <= 16`;
- `<= 3` for `n <= 14`;
- `<= 4` for `n <= 12`;
- `<= 5` for `n = 10`.

That is 22 runs and 401,442 strong arc sets, with every vertex checked.
- `M(T)` is never transitive and `F(T)` is always strong.
- The vertex lemma and "no majority source" hold. The minimum vertex-lemma value, 3-cycle certificate and `Pr[w first]` never reach `T_n`'s values (`T_n`, the only strong feedback-1 tournament, is checked separately). The `D_in` minimum fails from `n = 11`.
- (SS) never fails at `k <= 1`, and first fails at `k = 2` at `n = 10`.
- The vertex-lemma margin decays like `2^(−n)` even at vertices of in- and out-degree `>= 3`. So any proof must be exact, like the degree-`<= 2` injections, rather than a bound with slack.

**11.3 Out-degree-1 vertices.** Let `w → o` be the only out-arc of `w`.
- **PROVED.** If `o` beats every vertex of `U = V − {w, o}` (equivalently `indeg(o) = 1`), then `C(u, w) − C(w, u) = H(T[U])` for every in-neighbour `u`.
  - The paths ending at `w` are `o Q' w`.
  - The others are `A w o B'` with `(A, B')` free, and swapping `A` and `B'` is an involution.
  - This gives `T_n`'s extremal value as a special case.
- **FINITE-EXACT, no failure** (session runs plus the second audit's exhaustive `n = 8`). A witness always lies among the in-neighbours that `o` beats. The in-neighbours that beat `o` have positive average `2p − 1`, the conjectured next case.

**11.4 Negative results that constrain a proof.** No natural local choice of witness works:
- the predecessor-weighted average (first failure found at `n = 10`);
- minimum out-degree, minimum or maximum score in `T[U]`, and the most frequent predecessor (all by `n = 7`);
- "beats fewest out-neighbours of `w`" (`n = 8`, feedback 3; second audit).

The per-in-neighbour refinement of (SS) at `k = 1` fails in 41.3% of the triples `(T, w, u)` with `d >= 3` at `n = 7`. The vertex lemma looks true, with margin exactly `1/H` at `T_n`. These negative results rule out the simplest local witnesses and refinements, not every local argument.

**11.5 Second audit (MISTAKE-617).** An independent referee reproduced everything above with its own code and found three errors, now fixed:
- a claimed slack bound `0.047` that holds only for `3 <= d <= 6`;
- `T_n` minima claimed over the feedback families, including the `D_in` minimum, which fails there;
- an annealing value called an artefact at `n = 8..10`, where it is exact or unrefuted.

It also corrected eleven imprecisions (wording of the identity's hypothesis, script coverage, the `G` rule failing already at `n = 8`, the local-move typing, denominators). Report: `04-computation/experiments/s36_postaudit_audit_20261014/REPORT.md`.
