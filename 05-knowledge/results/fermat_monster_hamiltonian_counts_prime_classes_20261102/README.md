# Fermat primes, Monster layers and Hamiltonian path counts, merged with the prime square-sum structure: what is forced, what is counted, what is only numerals

klein-2026-11-02 (worktree `codex/session-ramsey-20261030`), continuing [`machines_automata_fractran_tournaments_dictionary_20261101`](../machines_automata_fractran_tournaments_dictionary_20261101/README.md) (THM-4860: the prime square-sum graph, its residue structure mod 4, 12, 36, and the arrangement of the first 50,396 primes with consecutive square sums). **Status: PROVED elementary + FINITE-EXACT (THM-4861: the five Fermat primes as Hamiltonian path counts of tournaments with their exact least orders 3, 4, 6, 8, 11; a vertex-transitive tournament has a prime number of Hamiltonian paths only if it is `C_3`; Fermat primes sit in the free classes of the prime square-sum graph), typed verdicts for the Monster/Baby Monster layers (COINCIDENCE with null probabilities), and HYP-9452 (open items). Every Paley count and every census value through order 9 is inherited (THM-4747, THM-4810, CONJ-002, `tournament_pseudoprimes_20261021`) and merely re-confirmed here; new are the Fermat least orders, the `C_3` corollary and the residue-class facts.**

Directive (2026-11-02, verbatim): "keep going how you were but also merge in ideas like the monster and baby monster, the fermat primes, and tournament hamiltonian path counts to make advances in all these areas simultaneously".

## 0. Inheritance (read newest first)

- **THM-4747** (mac-mini-2026-10-21): `|Aut(T)|` divides `H(T)` (free action, LEM-003); Fermat and Wilson congruences for circulant and Paley tournaments; reflection-reversed paths; subset DP to `n = 23`; the four Paley values `H(P_3), H(P_7), H(P_11), H(P_23)`. [`CONJ-002`](../../../01-canon/theorems/CONJ-002-paley-ham-paths.md) (REFUTED): `H(P_q) = |Aut| · 3^{(q−3)/2}` fails at 11 and 19; its table holds `H(P_19) = 1,172,695,746,915`. My C counter reproduces all five values exactly (§2).
- **THM-1370** (`h-spectrum-omits-7-21-all-n`): the Hamiltonian-path spectrum of tournaments omits exactly 7 and 21 for all orders, is monotone in `n`, and by `n = 8` covers every odd value in `[1, 609]` except 7 and 21; strong floors `f(m) = min{3^a 5^b : 2a + 3b = m − 1} = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375, …`. **THM-4810** (atlas, orders 3–8). **THM-4751** (circulant H-landscape: Paley peak at class-number distance from the rotational peak). **THM-133** (Paley maximises at order 7). OEIS A038375 (CITED): maximal `H` over tournaments of order `n`: `1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095` for `n = 1..11`.
- **THM-4644** (Fermat primes, rigid circulant tournaments: Akin's theorem that the multiplier action on circulant presentations of order `n` is free iff `n` is a squarefree product of Fermat primes) and the concurrent [`clocks_carriers_genus_20261102/FERMAT_NEIGHBOURS.md`](../clocks_carriers_genus_20261102/FERMAT_NEIGHBOURS.md) (orders `F_m − 2 = ∏_{j<m} F_j` and `F_m` have equal multiplier-orbit counts `2^{2^m−1}` … with unequal symmetry; the pairs (3,5), (15,17), (255,257), (65535,65537)). THM-4855: "3 and 5 Fermat primes is NUMEROLOGY" for the pure-blue circulants `K_1, C_3, T_5`.
- **The concurrent square-sum ⇄ tournament bridge**: [`joint_incidence_clocked_frames_20261102/PATH_KERNEL.md`](../joint_incidence_clocked_frames_20261102/PATH_KERNEL.md) defines the guarded response `H_G(T)` (directed Hamiltonian orders of the square-sum graph `G = Q_n` that are paths of the tournament `T`) with an exact spectrum compiler; 25 is the first order where two guard orientations share a path. That is the right place for "square-sum paths versus tournament path counts"; nothing here competes with it.
- Monster lanes: the Ogg split `Fi22 / B∖Fi22 / M∖B` is STRUCTURE (Aricheta; `monster_levels_square_sum_7_21_20261022.md`), every square-sum link to it NUMEROLOGY; THM-4752 (scores in `II_{25,1}`), THM-4758/4759 (frame shapes); the 15-prime square-sum tree (session32), re-proved by side imbalance in THM-4860 (c).

## 1. Fermat primes as Hamiltonian path counts (THM-4861 (a)–(c))

Rédei's theorem makes every `H(T)` odd, so every prime is a candidate value, and THM-1370 says every odd number except 7 and 21 occurs at some order. The exact **least orders** for the five Fermat primes, from the census of all tournaments of orders 3–9 ([`scripts/hampath.c`](scripts/hampath.c), nauty's `gentourng`; [`census_spectrum.py`](scripts/census_spectrum.py), [output](data/census_spectrum.out)) and the hill-climb witnesses ([`climb_65537_11.out`](data/climb_65537_11.out), [`climb_65537_12.out`](data/climb_65537_12.out); all witnesses recounted in pure Python, [`paley_circulant_counts.py`](scripts/paley_circulant_counts.py) (3)):

| Fermat prime | least order | tournaments of that order with this `H` | witness (upper triangle, row by row) |
|---|---|---|---|
| 3 | 3 | 1 (`C_3`) | `101` |
| 5 | 4 | 1 (the strong 4-tournament) | — |
| 17 | 6 | 2, both strong, scores `(1,1,2,3,4,4)` | `111011111111110`, `101011111111101` |
| 257 | 8 | 30 | `1111101111001111111111110110` (scores `2,2,2,3,4,4,5,6`) |
| 65537 | 11 | ≥ 3 found (two score sequences) | `0010010010100011000001001100011011000110101001101110111` |

Minimality: 5 > 3 = max at order 3; 17 > 15 = max at order 5; 257 > 189 = max at order 7; 65537 > 15,745 = max at order 10 (A038375; also the elementary bound `max H(n) ≤ n · max H(n−1) = 10 · 3357`). So `n(F_m) = 3, 4, 6, 8, 11` for `m = 0..4`.

**Vertex-transitive tournaments never have a prime count except `C_3` (THM-4861 (c)).** If `Aut(T)` is transitive on `n` vertices then `n | H(T)` (THM-4747), so `H(T)` prime forces `H(T) = n`; a vertex-transitive tournament is regular, hence strong, hence `H(T) ≥ f(n)` (THM-1370), and `f(n) > n` for `n ≥ 5` (`f(5) = 9`), while `n = 3` gives `C_3`. In particular no Fermat prime `≥ 5` is the count of a circulant tournament; the census agrees (`R_5` has 15, all order-7 circulants have 175 or 189, order 9: 3159, 3267, 3357, order 11: 92,411 … 95,095, each a multiple of the order with odd quotient).

**The order-9 census (re-confirmed, not new).** 191,536 tournaments, 1,520 distinct values, maximum 3,357 (attained by three of the sixteen circulant connection sets and matching A038375); the odd values below the maximum that are not attained are 7, 21 and a top band starting at 2,883. All of this was already computed by mac-mini-2026-10-21 ([`tournament_pseudoprimes_20261021/analyze_census.out`](../../../04-computation/experiments/tournament_pseudoprimes_20261021/analyze_census.out), cited by THM-4747); the present run reproduces it with an independent counter. Orders 6, 7, 8 show the same shape (gaps 35, 39 at order 6; a band from 63 at order 7; from 611 at order 8; THM-4810).

**Paley re-confirmation.** `H(P_q) = 3, 189, 95,095, 1,172,695,746,915, 15,760,206,976,379,349` for `q = 3, 7, 11, 19, 23` ([`paley.out`](data/paley.out)), all equal to the canon values (THM-4747, CONJ-002, THM-336). Quotients by `|Aut| = q(q−1)/2`: `1, 9, 1729, 6,857,869,865, 62,293,308,207,033`, odd as THM-4747 requires. `1729 = H(P_11)/55` is the taxicab number — NUMEROLOGY, recorded so that nobody promotes it.

## 2. Fermat primes and Monster layers inside the prime square-sum graph (THM-4861 (d)–(e))

THM-4860 split the primes into the bipartition `A = {2} ∪ {p ≡ 1 (4)}`, `B = {p ≡ 3 (4)}` and, mod 36, into **restricted** classes `1, 13, 25` (the squares mod 36, each adjacent to a single partner class `35, 23, 11` and to 3) and `7, 19, 31` (`7 ×` squares; partners `29, 17, 5` and 2), and **free** classes `5, 11, 17, 23, 29, 35`. ([`fermat_monster_classes.py`](scripts/fermat_monster_classes.py), [output](data/fermat_monster_classes.out).)

- **Fermat primes are free.** For `m ≥ 1`, `F_m = 2^{2^m} + 1 ≡ 1 (mod 4)` and `2^{2^m} ≡ 4, 7 (mod 9)` alternately, so `F_m ≡ 5, 17 (mod 36)` alternately: `5, 17, 257 ≡ 5, 65537 ≡ 17`. Two Fermat numbers are never adjacent (`F_m + F_k ≡ 2 (mod 4)`). `F_0 = 3` is the shared vertex of the three class-1 pieces. First square-sum prime partners: `3–13 (16)`, `5–11 (16)`, `17–19 (36)`, `257–67 (18²)`, `65537–2063 (260²)`; degrees below `10^5`: 45, 56, 82, 57, 31.
- **Monster layers.** `M∖B = {29, 41, 59, 71}` are all free (`29, 5, 23, 35`); `Fi22` and `B∖Fi22` hold two restricted primes each (`7, 13` and `19, 31`); the ten non-supersingular primes below 100 hold seven restricted ones, among them the three isolated primes `37 (1), 43 (7), 67 (31)`, whose first partners `107, 101, 257` are free (`35, 29, 5`) — the isolation of restricted small primes is the mechanism of THM-4860 (c), and 257 is a Fermat prime by coincidence. Null model for the Monster split: the 13 odd supersingular primes above 3 contain 4 restricted ones (`P ≈ 0.13` under a fair coin per prime), and `M∖B` is restricted-free (`P = 1/16`). **COINCIDENCE**, no map: Ogg's list is defined by `|M|`, the classes by squares mod 36.
- **Heegner aside (NUMEROLOGY).** The class-number-one primes `3, 7, 11, 19, 43, 67, 163` split as supersingular (`≤ 19`, exactly the Paley orders whose `H` the repository knows) and non-supersingular (`43, 67, 163`); `43` and `67` are two of the three isolated primes of the square-sum graph on the primes below 100. Typed and shelved.

## 3. What does not connect (typed, so it is not redone)

| chain | verdict |
|---|---|
| Fermat primes 3, 5 ↔ strong floors `f(m) = 3^a 5^b` | the floors come from `C_3` and the strong 4-tournament as atoms (THM-1370); the Fermat reading is NUMEROLOGY (THM-4855 already typed it for the pure-blue circulants) |
| Fermat primes ↔ `H(T)` | every odd number but 7, 21 is a count; the content is the least-order sequence `3, 4, 6, 8, 11` and the `C_3` corollary (THM-4861) |
| Monster primes ↔ `H(P_q)` | the supersingular Paley orders `3, 7, 11, 19, 23` are exactly those within reach of the subset DP; `31, 47, 59, 71` are not, and nothing distinguishes supersingular from other orders in the computed values |
| Monster layers ↔ residue classes of the square-sum graph | COINCIDENCE (§2) |
| square-sum Hamiltonian paths ↔ tournament `H` | the guarded response `H_G(T)` of PATH_KERNEL is the structural bridge; the prime square-sum path of THM-4860 (f) is a bipartite-flow object with no tournament behind it |

## 4. Reproduction

```bash
cd 05-knowledge/results/fermat_monster_hamiltonian_counts_prime_classes_20261102
cc -O2 -o scripts/hampath scripts/hampath.c
python3 scripts/census_spectrum.py            # orders 3-9 via gentourng (about 1 min); writes data/census_*.out
for q in 3 7 11 19 23; do scripts/hampath paley $q; done   # q = 23 needs about 0.8 GB
scripts/hampath climb 11 65537 1              # a witness with 65537 Hamiltonian paths (seconds)
python3 scripts/paley_circulant_counts.py     # pure-Python recounts (Paley q <= 11, circulants, all witnesses)
python3 scripts/fermat_monster_classes.py     # residue classes of Fermat and Monster primes
shasum -a 256 -c SHA256SUMS
```

Receipts in [`receipts/`](receipts/) are byte-identical under `python3` and `python3 -O`.
