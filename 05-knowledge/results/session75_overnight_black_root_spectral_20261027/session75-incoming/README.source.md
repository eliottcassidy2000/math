# Order-10 prime cycle-index types, the order-11 ledger, spectral parse obstruction, black-core distance of regular tournaments

2026-10-27, klein (worktree `codex/session-overnight-20261027`, seventh directive: "work to see the
order-10 prime cycle index types fall out as you pursue doing a long overnight session merging in
related to ... regular tournament on 2m+1 >= 7 vertices is exactly m black moves from a singleton
plus a strong core ... X = 3s+1 ... H(T) versus H(T)/|Aut(T)| ... {7,21} / {2,5} ... Eisenstein
1+2ω ... dyadic ledger repair ... trace norm O(n log n) versus Paley (n-1)√n"). Continues
[`a000568_intransitivity_families_20261026/`](../a000568_intransitivity_families_20261026/README.md)
(THM-4803, HYP-9413 (d)); canonical record: THM-4814, HYP-9422. Not independently audited.

Concurrency: while this packet was being written, `mac-mini-2026-10-27-potentials` reserved THM-4811–4813 and HYP-9421 (commit `3613d30439`) for "black terminal distance of arbitrary strong sources; tc fibres 7/21 beyond eleven vertices; Collatz numerator coordinate and positive carries; spectral exception budgets", the same directive themes from another angle; no file of this packet overlaps theirs, and the ids here were moved to 4814 / 9422.

| Part | Content | Status |
|---|---|---|
| A. Order-10 prime types (§1) | the nine automorphism / anti-automorphism cycle-index types of the 7,906,161 prime tournaments of order 10 (five with `|Aut| = 3`, three with `|Aut| = 5`, two rigid), by a registry-free nauty census in eight parts; the two cycle lemmas that explain every row; the self-converse decomposition identity `sc(n) = sc_nonstrong + sc_strong_nonprime + sc_prime` validated at every `n <= 10` | FINITE-EXACT; lemmas PROVED |
| B. Order-11 ledger (§2) | what the species identity determines at `n = 11` without the order-11 primes (non-strong 19,288,658; strong non-prime 87,588,161; labelled strong non-prime 3,453,658,999,907,072; self-converse non-prime 20,661); hence with OEIS A000568(11), A002785(11): **A280772(11) = 796,876,429** prime tournaments, 259,307 self-converse primes, **c(11) = 447,114,461** intransitivity classes (strong 442,336,617), **l(11) = 16,901,802,400,214,078** labelled, P_lab(11) = 31,808,433,417,949,440 labelled primes; HYP-9413 (d) resolved in this form | FINITE-EXACT given the CITED OEIS inputs |
| C. Spectral obstruction (§3) | `ν(T) = ‖S_T‖_*` (nuclear norm of the skew adjacency): `ν <= (n-1)√n` (n odd) with equality iff T is a one-vertex deletion of a skew-Hadamard tournament, `ν <= n√(n-1)` (n even) with equality iff T is skew-Hadamard (`S Sᵀ = (n-1)I`); the exact census of maximizers `n <= 9` (n = 8: the four skew-Hadamard tournaments, one of them `source + P_7`; n = 9: 216 classes at `6√7 + 2√15`); `ν(TT_n) = (2/π) n (ln n + γ + ln(4/π)) + o(n)`; the linear-parse bound and the exception lower bound `e >= Δ/4`, `e >= Δ²/(8n)`; comparison with THM-4810's score bound: the score bound wins for linear parses, the spectral bound is the only one of the two that survives regular tree parses (C3 towers); explicit table | PROVED (bounds, equality cases, comparison); FINITE-EXACT (census); the all-tree `O(n log n)` bound is a concurrent packet's, CITED as live |
| D. Black-core distance (§4) | from a regular tournament on `2m+1` vertices, every native path complement changes `|s(v) - m|` by at most one, so a vertex of score `0` or `2m` needs `>= m` moves (the endpoint rule and parity follow); all-black words of length exactly `m` to a singleton plus a strong core exist for every regular class at `n = 7` (3 classes; 30,716 / 41,046 / 32,928 black geodesic words) and `n = 9` (15 classes), and for `P_11` and twelve random regular 11-tournaments; at `n = 5` the regular tournament is the pure-blue atom `T5` (no black move at all; distance 2 with blue moves) | lower bound PROVED; witnesses FINITE-EXACT (n = 7, 9 complete), sampled at 11 |
| E. Lattices, ledgers, ROOT coordinate (§5) | a doubly regular tournament of order `q` makes the root lattice `A_{q-1}` a module over `Z[(1+√-q)/2]` through `θ = (I+S)/2` (the Eisenstein `1+2ω` is `q = 3`); the four-bin averaging identity and why no single pairing attains it; the corrected native edge identity `2^t Y = 3X + 2^a` checked on 75,000 seeds and tied to THM-4581 / THM-4793's merges; `{7, 21}` versus `{2, 5}` restated through the C3 fibre (H is C3-invariant at every order, H/|Aut| only through 11) | PROVED (lattice), CITED + checked (ledger, identity), DICTIONARY (gap pairs) |

What is new: the order-10 cycle-index types (A280772-level data beyond any published census), the
order-11 ledger (first values of A280772(11), of the intransitivity counts at 11, of the labelled
prime count at 11), the self-converse decomposition identity, the spectral equality
characterisation and its census, the black-move distance theorem for regular tournaments with its
witness census, and the `O_{-q}`-module structure. What is cited: OEIS A000568 / A002785 / A051337 /
A054946 (inputs and cross-checks), THM-4803 (fibre, species), THM-4810 and session74 (native moves,
score bound), THM-4730 (pure-blue alphabet), the paid-ROOT packet R1 and A1 and the concurrent
dyadic-module packet (edge identity, Eisenstein rows, four-bin average, all-tree trace-norm bound),
THM-1370 / THM-4790 / THM-4805 (gap pairs), S63's atlas (the gap-pair dictionary).

## 1. The order-10 prime cycle-index types

`scripts/tourn_types.c` reads `gentourng 10 i/8` (nauty 2.9.3), keeps the strong tournaments,
runs the modular decomposition of `tourn_enum.c` (THM-4803) and, for the prime ones, enumerates
`Aut(T)` from nauty's generators and the coset of anti-automorphisms (`ψ[lab[i]] = labc[i]` when the
canonical forms of `T` and `T*` agree). No registry is needed: `gentourng` lists each class once.
Each part took about a minute; the eight summaries (`data/census_n10_part*.txt`) add up to 9,733,056
classes, 9,355,949 strong, 7,906,161 prime, 5,525 self-converse primes, `|Aut|` histogram of the
primes `1: 7,905,266; 3: 884; 5: 11` (the earlier run of THM-4803 had the same totals). The merged
types (`data/prime_types_n10.txt`; full table `n = 3..10` in `data/prime_cycle_types_n3_n10.tsv`):

| `|Aut|` | Aut cycle types (elements) | self-converse | anti-automorphism cycle types | classes |
|---|---|---|---|---|
| 1 | `1^10` | no | – | 7,899,762 |
| 1 | `1^10` | yes | `2^5` (1) | 5,504 |
| 3 | `1^10` (1), `1^4 3^2` (2) | no | – | 684 |
| 3 | `1^10` (1), `1^1 3^3` (2) | no | – | 184 |
| 3 | `1^10` (1), `1^4 3^2` (2) | yes | `2^2 6^1` (2), `2^5` (1) | 8 |
| 3 | `1^10` (1), `1^4 3^2` (2) | yes | `2^5` (3) | 8 |
| 5 | `1^10` (1), `5^2` (4) | no | – | 6 |
| 5 | `1^10` (1), `5^2` (4) | yes | `10^1` (4), `2^5` (1) | 3 |
| 5 | `1^10` (1), `5^2` (4) | yes | `2^5` (5) | 2 |

Two elementary lemmas account for every row here and in the `n <= 9` table (PROVED):

- **Cycle lemma.** An automorphism of prime order `p` of a prime tournament on `n > p` vertices has
  at least two `p`-cycles. (Automorphisms of tournaments have odd order; if `g` had a single
  `p`-cycle `C`, every fixed vertex `v` would satisfy `v -> c iff v -> g(c)`, so `v` relates uniformly
  to `C`, and `C` would be a nontrivial module.) Hence the types `1^4 3^2`, `1^1 3^3`, `5^2` and never
  `1^7 3^1` or `1^5 5^1`.
- **Anti-automorphism lemma.** An anti-automorphism `ψ` of a tournament has at most one fixed
  point and no cycle of odd length `>= 3` (`ψ^l` is an anti-automorphism for odd `l` and would fix
  the cycle pointwise; two fixed points `u, v` give `u -> v iff v -> u`). Hence `2^5`, `2^2 6^1`,
  `10^1` at `n = 10` and `1^1 2^4`, `1^1 2^1 6^1` at `n = 9`.

**Self-converse decomposition (PROVED; `scripts/species11b.py`).** A tournament with prime
quotient `Q` and modules `M_i` is self-converse iff some anti-automorphism `ψ` of `Q` satisfies
`M_{ψ(i)} ≅ M_i*` for all `i`. Counting `Aut(Q)`-orbits of module assignments fixed by the
dualising action of `Aut±(Q)` by Burnside (`fixed = 2·#Aut±-orbits − #Aut-orbits`; a `ψ`-cycle of odd
length forces a self-converse content, counted by `sc(x^l)`, an even one is free, `a(x^l)`), and
counting non-strong self-converse tournaments as dual-palindromic sequences of strong classes
(even length: free half; odd length: a self-converse strong centre), gives
`sc(n) = sc_nonstrong(n) + sc_strong_nonprime(n) + sc_prime(n)`. Evaluated from the type census it
returns the census values `1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784` at every `n <= 10`
(`data/species11b.out`), which is an independent check of the anti-automorphism data of all 32 types.

## 2. The order-11 ledger

The species identity `A = SEQ(S)`, `S = X + Σ_Q Z_{Aut Q}[A − 1]` (session56; THM-4803) contains the
order-`n` primes themselves as the term `x^n` of each `Z_{Aut Q}` with `|Q| = n`, so it cannot predict
`A000568(11)` without `P(11)`; what it determines exactly from the orders `<= 10` is everything else.
`scripts/species11.py` (Fraction arithmetic over the 32 types) and `species11b.py` give:

| quantity at `n = 11` | value | source |
|---|---|---|
| tournament classes `A(11)` | 903,753,248 | OEIS A000568, CITED |
| non-strong classes | 19,288,658 | `SEQ(S)` from `S(<= 10)` |
| strong classes `S(11)` | 884,464,590 | `A − 19,288,658`; equals OEIS A051337(11) |
| strong classes with a nontrivial module (prime quotient of order `<= 10`) | 87,588,161 | Pólya over the 32 types |
| **prime tournaments `P(11)`** | **796,876,429** | `S(11) − 87,588,161` (A280772 stops at `n = 8` on OEIS) |
| self-converse classes `SC(11)` | 279,968 | OEIS A002785, CITED |
| self-converse non-strong / strong non-prime | 2,768 / 17,893 | §1 identity |
| **self-converse primes** | **259,307** | difference; `P(11) − 259,307` is even, as converse pairs require |
| classes up to converse `(A + SC)/2` | 452,016,608 | |
| strong intransitivity classes `c_s(11)` | 442,336,617 | `43,768,749 + (796,876,429 + 259,307)/2` |
| **intransitivity classes `c(11)`** | **447,114,461** | `48,546,593 + 398,567,868`; `c/A = 0.4947` (THM-4803: `-> 1/2`) |
| labelled strong `s_lab(11)` | 35,262,092,417,856,512 | from `2^55` by `T = 1/(1 − S_lab)`; equals OEIS A054946(11) |
| labelled strong non-prime | 3,453,658,999,907,072 | `Σ_k P_lab(k)/k! · 11! [x^11] T₊^k` |
| **labelled primes `P_lab(11)`** | **31,808,433,417,949,440** | difference; `P_lab/2^55 = 0.8829` (0.82 at 10) |
| **labelled intransitivity structures `l(11)`** | **16,901,802,400,214,078** | `997,585,691,239,358 + P_lab(11)/2`; `l/2^55 = 0.4691` |

The same programs reproduce every `n <= 10` value of THM-4803 (`A`, `S`, `c`, `c_s`, `l`, `P_lab`)
from the type files alone. The labelled identity is `l_s(x) = x + ½ Σ_k P_lab(k) l(x)^k / k!`,
`1 + l = exp(l_s)` (a prime quotient contributes `1/|Aut±|`, and `Σ_Q 1/|Aut±(Q)| = ½ Σ_Q 1/|Aut Q|`
whether or not `Q` is self-converse), so `l(11)` needs only the labelled prime count, which is exact.
The Burnside defect `11!·P(11) − P_lab(11) = 323,623,157,760` corresponds to about 12,200 non-rigid
order-11 primes, in line with 895 at order 10. HYP-9413 (d) asked for `A000568(11)`, `c(11)` from the
order-10 types: the honest form is this ledger (the identity gives the non-prime side; the known
`A(11)` then gives the primes). `A(12)` would need the order-11 types.

## 3. The spectral obstruction: nuclear norm of the skew adjacency

Let `S_T` be the skew adjacency matrix (`+1` for `i -> j`), `ν(T) = ‖S_T‖_*` the sum of its singular
values (`= Σ|λ|` over the eigenvalues of the Hermitian `iS_T`). `scripts/nuclear.py`, `nucmax.py`,
`nucmax9sv.py`, `parse_bound.py`; outputs in `data/`.

**Theorem (PROVED).** `‖S_T‖_F² = n(n−1)` for every tournament, so by Cauchy–Schwarz over the
nonzero singular values: for odd `n`, `ν(T) <= (n−1)√n`, with equality iff the `n−1` nonzero singular
values all equal `√n`, i.e. `S Sᵀ = nI − wwᵀ` with `w` a `±1` vector in `ker S` (the diagonal
`n−1 = n − w_i²` forces `w_i = ±1`); bordering by `w` gives `S' S'ᵀ = nI`, so equality holds iff `T`
is a one-vertex deletion of a skew-Hadamard tournament of order `n+1` (then `n ≡ 3 mod 4`), the
doubly regular ones being the deletions with `w = ±1`. For even `n`, `ν(T) <= n√(n−1)` with equality
iff `S Sᵀ = (n−1)I` (skew-Hadamard; `n ≡ 0 mod 4` or `n = 2`). Census (`data/nucmax_n*.out`): `n = 4`:
two maximizers (`C3` with a source or a sink); `n = 7`: six (`P_7` and the five other deletions);
`n = 8`: four, the skew-Hadamard tournaments of order 8, with score sequences `3^7 7` (source + `P_7`),
`0 4^7`, `2^3 4^4 6`, `1 3^4 5^3`; no equality at `n = 5, 6, 9` (ratios 0.9789, 0.9636, 0.9842); the 216
maximizers at `n = 9` all have singular values `√7 (×6), √15 (×2), 0`, `ν = 6√7 + 2√15`. Random
tournaments have `ν ≈ (8/3π)(n−1)√n` (semicircle; ratio 0.851 at `n = 191`), so Paley is extreme but
not isolated. `ν(TT_n) = Σ_k |cot((2k−1)π/2n)|` (cot formula checked against SVD to `n = 400`)
`= (2/π) n (ln n + γ + ln(4/π)) + o(n)` (constant 0.52125 per vertex at `n = 10^3..10^5`). The `C3`
tower has `ν = (2/√3) n log_3 n` (regular-block orthogonality, THM-1960), and the caterpillar
`C3[v, v, C3[v, v, …]]` stays at about `0.9 n ln n` (`n = 640`): it is `TT_n` with a perfect matching
of antipodal pairs reversed, so `ν <= ν(TT_n) + 2n`.

**Linear parses (PROVED).** For `T' = TT_k[B_1, …, B_k]` with `|B_i| = m_i <= b`,
`S_{T'} = (S_{TT_n} − ⊕S_{TT_{m_i}}) + ⊕S_{B_i}`, so `ν(T') <= ν(TT_n) + Σ_i (ν(TT_{m_i}) + ν(B_i))
<= ν(TT_n) + κ(b) n`, `κ(b) = max_{m<=b} (ν_max(m) + ν(TT_m))/m`. For every bounded-arity
substitution tree (any shape) the concurrent packet `dyadic_module_spectrum_20261026` (another
agent's worktree, read 2026-10-27, not on main when this was written) proves `ν <= (2b−3) n ⌈log₂ n⌉`
by ordering each node around a heavy child; the caterpillar numbers above are consistent with it.

**Exception bound (PROVED).** If `T` and `T'` differ on `e` pairs, `X = S_T − S_{T'}` is skew with
`±2` entries on `2e` positions, `‖X‖_F² = 8e`, `rank X <= min(n, 2e)`, hence
`Δ := ν(T) − ν(T') <= ‖X‖_* <= min(4e, √(8en))`: `e >= Δ/4` when `Δ <= 2n`, and `e >= Δ²/(8n)` beyond.
For `P_q`, `Δ = (q−1)√q − ν(T')`, so every fixed bank needs `e >= q²/8 − (q^{3/2}/π) ln q (1 + o(1))`
exceptions. **Comparison with the score bound (THM-4810 Theorem L).** For a linear parse the sorted
scores of `T'` deviate from `0..n−1` by at most `b−1` each, so `L(T') <= n(b−1)/2` and
`e >= L(P_q) − L(T') >= (q²−1)/8 − q(b−1)/2`: the same leading constant with a linear instead of a
`q^{3/2} log q` error; for `b = 1` it is the classical `f(P_q) >= (q²−1)/8`. The score bound is void
for tree parses with regular quotient nodes (`C3` towers are regular, `L(T') = L(P_q)` is possible)
while the spectral bound is unchanged there: the two obstructions have different scopes, and neither
dominates. Explicit values (`data/parse_bound.out`): transitive parse `e >= 3, 9, 28, 66` at
`q = 11, 19, 43, 47` (true `f = 20, 64, …`; score bound 15, 45, 231, 276); the join/`C3` tree grammar
with the `(2b−3) n ⌈log₂ n⌉` constant is nonvacuous from `q ≈ 10³` (`e >= 482` at `q = 1019`,
`37,521` at `2039`, `2,714,845` at `8191`).

## 4. Black moves from regular tournaments to a singleton plus a strong core

Native move (THM-4810, session74): `τ_Q(T) = T^op xor E(Q)` for a directed Hamiltonian path `Q`
(keep the arcs of `Q`, reverse every other pair); `Q` is blue in `T` iff its positional reversal is an
anti-automorphism of `T`, black otherwise (THM-4729 / THM-4732 dictionary). Scores transform by
`s'(v) = n − 1 − s(v) + [v = start Q] − [v = end Q]`.

**Theorem (PROVED).** Let `T` be regular on `n = 2m+1` vertices. Any sequence of native moves
(black or blue) after which some vertex `a` has score `0` or `2m` has length `>= m`: with
`d(v) = s(v) − m`, one move sends `d(v)` to `−d(v) + [start] − [end]`, so `|d(v)|` grows by at most one
per move, and a singleton strong component needs `|d(a)| = m`. Along a word of length exactly `m` the
vertex `a` is an endpoint of every path (the start when `d(a) <= 0`, the end when `d(a) >= 0`), `|d(a)|`
equals the depth, and `|d(a)| ≡ depth (mod 2)`.

**Witnesses (`scripts/blackcore.py`; FINITE-EXACT).** Searching only over black paths with the
endpoint rule and requiring the final `T − a` to be strong: at `n = 7` every regular class has a
black word of length 3 (exhaustive counts of black geodesic words over all choices of `a`: 30,716,
41,046, 32,928; `data/../receipts/blackcore7-*.txt`); at `n = 9` all 15 regular classes have a black
word of length 4 (`data/bc9.out`); at `n = 11`, `P_11` and twelve random regular tournaments (a
3-cycle-reversal walk from `P_11`) have black words of length 5 (`data/bc11p.out`, `bc11r.out`). So
the distance is exactly `m` there, with the paths alternating `a` as start and end (the "two distinct
path starts"). At `n = 5` (`data/bc5.out`) the regular tournament is `T5`, the strong pure-blue atom
of THM-4730: all 60 endpoint-rule two-move words use a blue path (50 of them do land on a singleton
plus a strong core), and the black-move graph from `T5` is empty, so its black distance is infinite
while the unrestricted distance is 2. This is why the claim starts at `2m+1 >= 7`. HYP-9422 (1)
states the all-order claim.

## 5. Lattices, ledgers and the ROOT coordinate

**Paley root lattices (PROVED; `scripts/misc_checks.py`).** For a doubly regular tournament of order
`q` (`S² = J − qI`, `S1 = 0`) and `v` in the sum-zero lattice `Λ = Zⁿ ∩ 1^⊥ ≅ A_{q−1}`,
`((I+S)v)_i = 2v_i + 2Σ_{j ∈ N⁺(i)} v_j`, so `θ = (I+S)/2` is an integral endomorphism of `Λ` with
`θ² − θ + (q+1)/4 = 0` (`= (S² + qI)/4 = J/4 = 0` on `1^⊥`): `Λ` is a module of rank `(q−1)/2` over the
maximal order `Z[(1+√−q)/2]` of `Q(√−q)`, and `S = 2θ − 1 = √−q` on it. Checked for `q = 3, 7, 11, 19,
23`. For `q = 3` this is the paid-ROOT / dyadic-module statement that the directed three-cycle's
skew operator is `±(1+2ω)` on the centred triangular lattice `A_2 = Z[ω]` (sign = orientation
convention; `(1+2ω)² = −3`); the three ROOT rows `(2a−b, −a−b, −a+2b) = ι((2+ω)(a+bω))` are that
packet's L6, CITED. The Paley case is the natural generalisation: the "shared spectrum" is
`√−q` acting on a rank-`(q−1)/2` `O_{−q}`-module.

**Four-bin ledger repair (CITED, elementary).** For any quadratic form `Q` and four vectors, the
three pairings `p` have `Σ_p M_p = E + 2S` (each cross term once, each diagonal three times), so the
average middle energy is `(E + 2S)/3`; it lies between `S` and `E`, and
`(E − M̄)₊ + (M̄ − S)₊ = (E − S)₊` (paid-ROOT R3 / dyadic R1–R2). Added remark: `M_p` is an integer for
integer charges while `(E + 2S)/3` is integral iff `E ≡ S (mod 3)`, so for `E ≢ S (mod 3)` no single
pairing can attain the average (two of the three random Eisenstein examples in `data/misc_checks.out`
are of this kind). The directive's "433 versus 79" (pairwise positive parts against the positive part
after a level sum) is the inequality `Σ_k (Δ_{2^k})₊ <= Σ_pairs (β)₊`; the two numbers are not in
any committed packet and are not reproduced here.

**The ROOT coordinate and Collatz (checked; links).** With `X = 3s + 1`, `a = 2` for `s ≡ 1 (mod 8)`
and `a = 1` otherwise, `(3s+1)/2^a = B^j(y)`, `B(x) = 4x + 1`, `t = a + 2j`, `Y = 3y + 1`, the paid-ROOT
identity `2^t Y = 3X + 2^a` holds for all 75,000 odd seeds `s < 200001`, `s ≢ 5 (mod 8)` (`j <= 8`);
the uncorrected constant `2^a + 3` fails already at `s = 3` (`32 = 30 + 2 ≠ 35`). In this coordinate
every compressed edge is the affine map `X -> (3X + 2^a)/2^t` with positive carry, and a word of `d`
edges has `2^{T_d} X_d = 3^d X_0 + C_d`, `C_d > 0`. Two links to this repository's Collatz canon:
(i) the compressed map `P` strips exactly the `B`-chains, so the "4x+1 merges" of THM-4556 (ii) and
THM-4793 / Theorem M* (`S(a) = 4S(b) + 1`) are coincidences `P(s) = P(s')` of the native dynamics:
merging there is landing on one `B`-chain; (ii) on the Fermat line (THM-4793) `X_k = 3s_k + 1 = 4 s_{k+1}`
for `2^{n−1} − 1` steps, i.e. every edge has `a = 2`, `j = 0`, clock `t = 2` (with the odd-seed
convention for the root `y`), before the valuation-3 landing. The energy identity (paid-ROOT (3))
is the Eisenstein two-port polarisation applied to `b = 3X/(3X + 2^a)`; it is cited, not re-derived.

**`{7, 21}` and `{2, 5}` (DICTIONARY).** Raw `H` omits 7 and 21 at every order (THM-1370, PROVED);
`tc = H/|Aut|` omits them through order 11 (THM-4790, FINITE-EXACT), reproduced here through
order 9 from the nauty census (`H`, `|Aut|` columns of THM-4803's tables: `7, 21 ∉ tc(n)`, `n <= 9`,
and at `n = 6` also `27, 35, 39` are missing below 60). The structural reason the two scopes differ
after THM-4803's correction (session67, MISTAKE-680 family): `H` is a function of the cyclic-triple
hypergraph at every order, `|Aut|` only through order 11 (first failure `81/27` at 12), so `H/|Aut|`
is not a C3-invariant beyond 11. On the aliquot side, 2 and 5 are untouchable (finite), and "no other
prime is untouchable" follows from every even number `>= 8` being a sum of two distinct primes
(`p = 1 + q + r = s(qr)`; `3 = s(4)`, `7 = s(8)`), which THM-4805 records as CONDITIONAL; S63's
dictionary (proof design transfers, no map does) is adopted unchanged.

## Reproduction

```
cc -O2 -I/opt/homebrew/include/nauty scripts/tourn_types.c -L/opt/homebrew/lib -lnauty -o tourn_types
for k in 3 4 5 6 7 8 9; do gentourng -q $k | ./tourn_types $k --types t$k.txt; done      # < 1 s each
for i in 0 1 2 3 4 5 6 7; do gentourng -q 10 $i/8 | ./tourn_types 10 --types t10_$i.txt > s10_$i.txt & done   # ~1 min each
python3 - <<'X'   # merge the parts: sum counts per type line -> t10.txt
X
python3 scripts/species11.py 11      # A, S, c, c_s, l for n <= 11 from t3..t10.txt (47 s)
python3 scripts/species11b.py        # self-converse decomposition, labelled strong and prime counts
python3 scripts/nuclear.py ; for n in 4 5 6 7 8 9; do python3 scripts/nucmax.py $n; done ; python3 scripts/nucmax9sv.py
python3 scripts/parse_bound.py ; python3 scripts/misc_checks.py
gentourng -q 7 | python3 scripts/blackcore.py 7 gentourng 0 all ; gentourng -q 9 | python3 scripts/blackcore.py 9 gentourng 60 first
python3 scripts/blackcore.py 11 paley 40 ; python3 scripts/blackcore.py 11 randreg 40 12
```

`receipts/*-normal.txt` and `*-optimized.txt` (`python3 -O`) are byte-identical for every Python
script (`receipts/receipts_check.txt`). `data/` holds the type files for `n = 3..10` (and the eight
`n = 10` parts with their census summaries), the human-readable type table, and every program output
quoted above. Python 3.10 with numpy 2.2 (nuclear norms only); everything else is standard library.
