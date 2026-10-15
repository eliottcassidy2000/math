# Strong-tournament balance, golden Wieferich primes, Pisano periods and the icosahedron

2026-10-14, mac-mini-2026-10-14-tourn.

**Status.**
- **PROVED:** [THM-4655](../../01-canon/theorems/THM-4655-strong-tournament-balance-lemmas.md) (balance of transitive-plus-one strong tournaments; rotation lemma) and [THM-4656](../../01-canon/theorems/THM-4656-wall-sun-sun-golden-orbit-explosion.md) (Wall–Sun–Sun primes as golden Wieferich primes).
- **FINITE-EXACT:** all strong tournaments on `n <= 10` vertices; Pisano facts to `n = 20000`.
- **HYP:** [HYP-9279](../hypotheses/HYP-9279-strong-tournament-balance.md).
- **OPEN:** the tournament 1/3–2/3 conjecture HYP-9277 in general.

**The owner's brief (verbatim):** "try to prove the tournament 1/3-2/3 conjecture for strong tournaments. consider the idea that the square sum problem (arranging a list of integers such that the sums of each consecutive number is a square) and the graceful tree conjecture are each analogous roughly with the 3n-1 and normal collatz conjectures, respectively. try to figure out how to prove all square packing bounds simultaneously. consider the fibonacci series under all possible bases and the mapping of a modulus to a pisano period, how a modulus 11 maps to a length 10 and 19 to 18 and 12 to 24 and 20 to 60, and it feels like primes go one below them while composites go to increasingly composite numbers, and how the period of the modulus divided by the modulus maxes out at 6, where 10 maps to 10 and where 50 maps to 300 ... etc. and see how this relates deeply to path of least resistance recursion we have studied at the heart of number theory. (arXiv 2607.23926) reports 128 poset classes at order 14 with balance constant 1/3, matching Narayana's cows number c(14)−1 and linking to the 1/3–2/3 problem . consider the identities involving 11 coordinates and 55 pairs. consider the golden rectangle and the dodecahedron and icosahedron. think of how https://en.wikipedia.org/wiki/Wall%E2%80%93Sun%E2%80%93Sun_prime realte with our work regarding wieferich primes, think rigitdity circuits"

**Concurrent work.** opus-2026-10-14-S35 received an overlapping prompt earlier today. Its results are in [one_twenty_eight_atlas_balance_gold_20261014.md](one_twenty_eight_atlas_balance_gold_20261014.md) and [paley_least_resistance_runs_one_modular_form_20261014.md](paley_least_resistance_runs_one_modular_form_20261014.md):
- HYP-9277;
- the reductions to strong tournaments and to a majority order;
- the Narayana reading of the 128;
- HYP-9278 on rigidity circuits.

This note builds on those and does not repeat them.

## 1. The tournament 1/3–2/3 conjecture for strong tournaments

**Setup (S35).** Take a uniformly random Hamiltonian path of a tournament `T`, and let `δ(T)` be the best balance `min(p, 1 − p)` over pairs.
- HYP-9277 says that `δ >= 1/3` for every non-transitive `T`, with equality exactly at transitive compositions of singletons and `C_3`. That gives `A000930(n) − 1` classes, 128 at `n = 14`, which is the tournament mirror of the 128 extremal 14-element posets of arXiv:2607.23926.
- S35 proved that `δ` is the maximum over the strong components, so the whole conjecture lives on strong tournaments.

**New here.**

1. **A proved subclass (THM-4655 (1)).** Take a transitive tournament `u_1 → … → u_m` plus one vertex `z` with any pattern `w`. It is strong iff `w_1 = 1` and `w_m = 0`.
   - Every Hamiltonian path is two increasing runs around `z`.
   - The pair `(z, u_2)` is always balanced: `1/3 < P(u_2 before z) <= 2/3`, with equality only at `w = 1 0 1^(m−3) 0`. The lower bound `1/3` is approached by `w = 110^(m−2)`.
   - So the conjecture holds, with an explicit witness pair, for every tournament whose strong components are transitive after deleting one vertex.
   - The proof is a two-case count with two lemmas: an injection of 01-boundaries into later zeros, and a geometric bound on crossed weights. Every step is checked by computer for `m <= 14`, and the statement by brute force for `m <= 20` (524,286 words).
2. **A rotation lemma (THM-4655 (2)).** The `n` paths obtained by cutting one Hamiltonian cycle cannot sustain a `> 2/3` majority order: their average concordance with any linear order is at most `(2n−1)/(3n)`.
   - The proof is one line: from any vertex the forward distances are exactly `1, …, n−1`.
   - Corollary: a strong tournament with no balanced pair has fewer than `n/(n+1)` closable paths.
   - This is weak. Non-closable paths can carry a majority order, and in the transitive-plus-one family most pairs among the `u`'s really are unbalanced (`≈ 3/4`). The balance must be found on special pairs.
3. **The Kahn–Saks route is closed (THM-4655 (3)).** Stanley's log-concavity, which gives `3/11` for posets, fails for tournament Hamiltonian paths already at `n = 4`. Every strong tournament at `n = 6, 7` fails it for gaps.
4. **Exhaustive to `n = 10` (THM-4655 (4)).** All 9,355,949 strong 10-tournaments satisfy `δ >= 62/157 = 0.3949`, and HYP-9277's inequality holds for every tournament on `n <= 10` vertices. The strong minima are:

   | `n` | minimum `δ` |
   |---|---|
   | 3 | `1/3` |
   | 4 | `2/5` |
   | 5, 6, 7 | `4/11` |
   | 8, 9, 10 | `62/157` |
5. **Vertex-wise balance (THM-4655 1(e)–(g), HYP-9279(B)).** In every strong tournament checked, *every vertex* has a balanced partner. That holds for all strong tournaments with `n <= 10` (9,355,949 at `n = 10`).
   - On the transitive-plus-one family it is proved. The extra vertex `z` is balanced against every interior `u_i`, by the same wide-crossing mechanism. `u_1` and `u_m` are balanced against `z` or their transitive neighbour.
   - The unique Hamiltonian cycle `z u_1 … u_m z` carries at least two balanced arcs.
   - The vertex-wise constant is exactly `1/3` for every `n = 5..9`. The tight vertices sit in a `C_3` block hung at the end of the order, the local version of the global extremals.
   - The stronger statement (B⁺), that the balanced graph is connected, holds for `n <= 9`.
6. **Cut vertices: the transitive-plus-one mechanism in general (THM-4655 (5), HYP-9279 (Z)).** Let `z` be a vertex whose deletion leaves `T − z` non-strong, with strong components `D_1 ⇒ … ⇒ D_k`.
   - The Hamiltonian paths factor through an exact 4-state transfer chain over the blocks. Each block enters only through its own subset path counts.
   - The interior blocks of the crossing window are fair coins against `z`, because `Σ_X h(X) h(D ∖ X)` is symmetric under complement.
   - Conjecture (Z): `z` is balanced against every vertex of `D_2 ∪ … ∪ D_(k−1)` and always has a balanced partner.
     - Exhaustive to `n = 9`, sampled to `n = 46`, sharp at `1/3`.
     - Proved when all blocks are singletons.
     - It would settle HYP-9277 for every strong tournament with a cut vertex: about 30% of strong tournaments at `n = 10`.
   - A local refinement of (B), (B_loc), also holds to `n = 9`: every vertex lies on a forward arc with `1/2 < p <= 2/3` on a directed 3-cycle. This ties (B) to S36's forward-arc digraph.
7. **A refuted strengthening.** "Every Hamiltonian cycle of a strong tournament carries a balanced arc" holds for all 133,632 Hamiltonian cycles at `n = 8` and on the transitive-plus-one family. It fails at `n = 9`.
   - The unique counterexample among the 178,133 strong 9-tournaments has a Hamiltonian cycle whose nine arcs all have `p(tail first) ∈ (0.671, 0.831)`: a `> 2/3` majority cycle.
   - The tournament itself has `δ = 100/213` and every vertex balanced.
   - So cycles can carry a majority orientation while vertices and pairs cannot avoid balance. This is the right warning for S36's majority-cycle route (HYP-9280): the majority structure lives on cycles, the balance on vertices.
8. **HYP-9279(A).** Every strong tournament on `n >= 4` vertices has `δ >= 4/11`, a strict gap above 1/3. The exhaustive minima are 4/11 at `n = 5, 6, 7` and 62/157 at `n = 8, 9, 10`. Within the transitive-plus-one family, the minimum rises towards 1/2:

   | `n` | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | min `δ` | 2/5 | 4/11 | 8/19 | 2/5 | 4/9 | 4/9 | 113/241 | 8/17 | 116/241 | 16/33 | 400/817 |

**Where a proof of the general strong case would have to go.**
- The minimisers are near-transitive, with feedback arc sets of sizes 1, 2, 3, 5, 4, 6 for `n = 4..9`. The `n = 8` minimiser is transitive plus 2.
- The natural induction is "transitive plus `k`": delete a vertex whose removal keeps the tournament strong (Moon), and track the balanced pair.
- The transitive-plus-one theorem is the base case and shows the mechanism: the extra vertex is fair against the second element of the order, and in fact against every interior element.
- The vertex-wise form (B) is the better induction hypothesis. It is sharp at `1/3`, it is local, and it survives exhaustively to `n = 10`. The Hamiltonian-cycle form is false.

## 2. Wall–Sun–Sun primes and the Wieferich work (THM-4656)

- **The dictionary.** THM-4644 (opus S33) made Wieferich primes the orbit explosion of doubling-invariant circulant tournaments on `Z_(p²)`. Replace the unit `2 ∈ Z` by the golden unit `φ ∈ Z[φ]`, or by any metallic mean `μ_a`, root of `x² − ax − 1`.
- **Split primes.** A root `φ₁ ∈ Z/p²` is a Wieferich base exactly at Wall–Sun–Sun primes, and THM-4644's count transfers verbatim.
- **Inert primes.** On `Z[φ]/(p²) ≅ (Z/p²)²` the exponent jumps from `(p+1)(p²−1)/(2π')` to `(p²+1)(p²−1)/(2π')`.
- **Evidence.** No Wall–Sun–Sun prime exists below `2^64` (PrimeGrid), so the jump was tested at metallic analogues. For the silver ratio (`a = 2`, Pell), `p = 13` and `p = 31` are Wall–Sun–Sun-type primes, and brute-force orbit counts show the predicted jump: 4080 against 336, and 64 against 4.
- **Rigidity (the mirror of Fermat primes).**
  - Golden rigidity holds iff the Pisano period `π(p)` is a power of 2, iff `p | L_(2^j)`. The primes below `10^6` are 3, 7, 47, 127, 1087, 2207, 4481, 21503, 34303, 119809, 524287.
  - This list contains every Mersenne prime `2^q − 1` with `q ≡ 3 (mod 4)`, with the maximal period `2^(q+1)`. That is exactly the Lucas–Lehmer test with seed 3.
  - Fermat primes (`p − 1 = 2^k`) make the split torus a 2-group. Mersenne primes (`p + 1 = 2^q`) make the golden nonsplit torus a 2-group.
- **"Think rigidity circuits."** Here "rigidity" is THM-4644's multiplier rigidity. S35's HYP-9278 develops the other sense, 2D rigidity circuits (`|E| = 2|V| − 2` with every proper subgraph Laman-independent) as the Calabi–Yau cores of Paley pattern sums, with motivic weight `|V| − 2`.
  - The K6 − P6 terminal of HYP-9276 is such a circuit: 10 = 2·6 − 2 edges, and no `K_4` inside, since the complement of `P_6` has independence number 3. Its conjectural motive `Sym²(ρ_g)` has weight 4 = 6 − 2, consistent with HYP-9278.
  - The wheel `W_4` (5 vertices, 8 edges, Apéry weight-4 form, motivic weight 3) is the previous rung.
  - See also §5.

## 3. Pisano periods

**Facts (KNOWN, re-checked to `n = 20000`; `pisano_facts.out`).**
- Split primes (`p ≡ ±1 mod 5`) have `π(p) | p − 1`. Equality, "one below", holds at 11, 19, 31, 41, 59, …
- Inert primes (`p ≡ ±2 mod 5`) have `π(p) | 2(p+1)`. Equality holds at 7, 13, 17, 23, 37, …; for example `7 → 16`.
- `π(5) = 20`.
- For composites, `π(n) = lcm_p π(p^k)` with `π(p^k) = p^(k−1)π(p)` unless `p` is Wall–Sun–Sun. This is the "increasingly composite" growth.
- The ratio `π(n)/n` is at most 6, with equality exactly at `n = 2·5^k` (10, 50, 250, 1250, 6250). The owner's examples:

  | `n` | 10 | 11 | 12 | 19 | 20 | 50 |
  |---|---|---|---|---|---|---|
  | `π(n)` | 60 | 10 | 24 | 18 | 60 | 300 |

  (`π(10)` is 60, not 10; it is the first ratio-6 case.)
- **Why the extremal moduli are `2·5^k` (an explanation of the equality case, not a proof of the bound).**
  - 5 is the ramified prime of `Q(√5)`: `φ ≡ 3 (mod 5)` is a double root, and `π(5^k) = 4·5^k`. Ratio 4 is the largest among prime powers; every other prime power has `π(p^k)/p^k <= 8/3`, attained at `p = 3`.
  - 2 is inert, with `π(2) = 3`.
  - `lcm(3, 4·5^k) = 12·5^k = 6·(2·5^k)`. For `p != 5`, `π(p^k)/p^k <= 2(p+1)/p <= 8/3`, and `π(2^k)/2^k = 3/2`.

**Least resistance (the corpus's sense, S35).** Least resistance means series-parallel reduction. Continued fractions are its arithmetic.
- The `n`-rung ladder of unit resistors has resistance `F_(2n+1)/F_(2n)` (standard termination), which tends to `φ` (KNOWN).
- The Pisano period is the period of this least-resistance recursion modulo `n`.
- In S35's Paley analysis, series-parallel patterns give polynomials. The first resistances, the wheel and `K_6 − P_6`, give modular forms.
- So the golden ratio is the fixed point of least resistance, and its arithmetic (Pisano periods, Wall–Sun–Sun primes) measures how least resistance fails modulo `p`.

## 4. The golden rectangle, dodecahedron and icosahedron

- **STRUCTURAL (KNOWN in substance).** The Fibonacci matrix `Q = [[1,1],[1,0]]` satisfies `Q^5 ≡ 3I (mod 5)`. Since `det Q = −1 = 2²` is a square mod 5, `Q` is projectively an element of order 5 in `PSL(2,5) ≅ A_5`, the rotation group of the icosahedron: a 5-fold vertex rotation.
- **The two exceptional isomorphisms `A_5 ≅ PSL(2,4) ≅ PSL(2,5)`** are `Z[φ]` reduced at its two special primes:
  - the inert prime 2, where `F_4 = F_2[φ]` and `φ` has order `π(2) = 3`;
  - the ramified prime `√5`, where `φ ≡ 3`.
- The icosians, the quaternion order over `Z[φ]` containing the 120 unit elements of the binary icosahedral group, carry this structure globally and give `E_8` (CITED: Conway–Sloane). This links to S35's `E_8` half-spin 128.
- The three mutually perpendicular golden rectangles span the 12 icosahedron vertices `(0, ±1, ±φ)`, the projective image of the golden ratio.
- **COINCIDENCE of orders.** `π(10) = 60 = |A_5|`, but `⟨Q mod 10⟩` is cyclic, not `A_5`. Also `π(5) = 20` is the number of dodecahedron vertices.

## 5. Eleven coordinates and 55 pairs

- **STRUCTURAL.** 11 is a golden Paley prime: the residue root `φ₁ = 4` generates `QR_11`. So the golden-invariant circulant tournament on `Z_11` is the Paley tournament `P_11`. `Aut(P_11) = {x ↦ ax + b : a ∈ QR_11}` has order 55 and acts regularly on its 55 arcs, as `Aut(P_p)` does for every Paley prime.
- **STRUCTURAL (rigidity).** 11 points in `R^10` with all 55 distances form the 10-simplex: `K_11` is minimally rigid in `R^10`, since `d·n − C(d+1,2) = 110 − 55 = 55`. `K_12` is a rigidity circuit there. This is the sense in which "11 coordinates and 55 pairs" are balanced.
- **COINCIDENCE.** `55 = F_10 = F_(π(11))` is also `C(11,2) = T_10`. By Luo's theorem (1989), 55 is the largest triangular Fibonacci number, so this cannot recur.

## 6. The analogies of the brief (typed)

- **Square-sum ↔ 3n−1.**
  - The square-sum path exists for `n = 15, 16, 17, 23` and all `n >= 25`, and is reported proved by Gerbicz (2018; OEIS A090461, forum-level, with cycles for `n >= 32`). Its finitely many exceptional small `n` mirror the finitely many extra cycles of `3n − 1` (`{1,2}`, `{5,14,7,20,10}`, `{17,…}`).
  - ANALOGY only. The square-sum problem is solved by a constructive extension algorithm; `3n−1` is open.
- **Graceful trees ↔ Collatz.** Both are universal statements ("every tree", "every orbit"), verified far out and open. ANALOGY.
  - A concrete bridge to test: is the Collatz tree truncated at `N` gracefully labelable by a rule compatible with its parity structure?
- **"All square packing bounds simultaneously."** Not attempted; no unifying mechanism was found. The nearest coincidence is the cannonball identity `1² + … + 24² = 70²` (Watson's unique solution; the Leech-lattice vector) beside the square-sum threshold 25. COINCIDENCE.

## 7. Reproduction

| Script | What it does | Time |
|---|---|---|
| `gentourng -c n \| ./tbal n` | exhaustive strong-tournament balance | 19 s at `n = 9` |
| `gentourng -c n \| ./tbal2 n` | the same, plus per-vertex partners and cut vertices | about 45 min at `n = 10` |
| `tbal3.c`, `tbal4.c`, `tbal5.c` | (B⁺), cycle-arc and sharpness checks; `c_counterexample.py` | a few minutes at `n = 9` |
| `tpo_cycle_arcs.py`, `tpo_vertexwise_proof_check.py` | THM-4655 1(e)–(g) | |
| `tbal6.c`, `tbal7.c` | (B_loc); (Z) exhaustive | about 1 min at `n = 9` |
| `cutvertex_chain.py`, `cutvertex_search.py` | THM-4655 (5); (Z) sampling and adversarial search | minutes |
| `trans_plus_one.py`, `tpo_proof_steps.py`, `tpo_u2_check` | the transitive-plus-one theorem | |
| `rotation_lemma_check.py`, `logconcave_test.py` | rotation lemma; log-concavity failures | |
| `wss_orbits.py`, `golden_fermat_check.py`, `pisano_facts.py` | THM-4656 and §3–4 | |

The first four rows are in `04-computation/experiments/tournament_balance_20261014/`, the last in `04-computation/experiments/wall_sun_sun_20261014/`.
