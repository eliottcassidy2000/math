# Three exact pieces of the {7, 21} programme: the bounded-growth lemma for literal span-3 words, the seam product decoder in reverse, and two-apex completeness below 2^15

**Session:** klein-2026-11-15-ribbon (branch `codex/session-ribbon-20261115`). The owner's directives were: choose factor pairs from the desired count so that their exact remainders land on the target; build the span-3 rank-potential certificate and prove the bounded-growth lemma (one added vertex multiplies H by at most a constant, observed 39); and ask, on the conflict graphs tournaments can actually have, whether the independence polynomial at 2 takes every odd value except 7 and 21 — one apex gives cliques and half the odd numbers, the other half must come from a second apex or a disjoint cycle at controlled overlap. Sources: [HYP-9802](../../hypotheses/HYP-9802-literal-span-3-certificates-need-a-bounded-growth-lemma-and-rank-adaptive-memory-short-cheap-cycles-are-cut-within-nine-letters.md) and wave 39 of [`chanpak_waves_20261113`](../chanpak_waves_20261113/README.md); [`chanpak_signed_seam_bank_20261116/PRODUCT.md`](../chanpak_signed_seam_bank_20261116/PRODUCT.md); [`chanpak_apex_ocf_literal_cones_20261115/UNIFYING.md`](../chanpak_apex_ocf_literal_cones_20261115/UNIFYING.md) §3.

## 0. Status

- **PROVED (elementary, §1.1).** For a literal span-3 word `W` and a letter `ℓ`, `H(Wℓ) = e_ℓ(W) + s_ℓ(W) + Z_ℓ(W)` with `e, s ≤ H(W)` (paths of `W` to which the new vertex is appended, resp. prepended) and `Z_ℓ(W)` the number of two-path covers of `W` bridged by the new vertex (one fragment's head among the last three vertices and beaten by the new vertex, the other fragment's tail not beaten by it). So the bounded-growth lemma `H(Wℓ) ≤ K·H(W)` is equivalent to `Z_ℓ(W) ≤ (K − 2)·H(W)`, and the observed `K = 39` is `1 + 0 + 38`: the maximisers have `H(W) = 1`, one appendable path, no prependable one and 38 bridgeable two-covers.
- **FINITE-EXACT (§1.2).** The maximum of `H(Wℓ)/H(W)` over all words of length `n` is 1, 3, 5, 11, 23, 39, 39, 39 for `n = 0..7` (independent recomputation of wave 39's values), attained for `n ≥ 5` only by `W = 0^(n−3) 7 6 0`, `ℓ = 6` (letters coded by bits: 7 = {1,2,3}, 6 = {2,3}, 0 = ∅); the two-fragment part `Z/H` is at most 38, attained there; lengths 8–9 are the pre-registered Q1–Q2 (§5). The exact envelope `env[H] = max H(Wℓ)` over **all** words with `H(W) = H` is computed for every attained `H ≤ 40` by closing the set of transfer states with `H ≤ 40` (16,794 states, §1.3): it is exactly `39·H` at every attained `H`, because three sinks seal a word and `W' 0 0 0 7 6 0` followed by {2,3} reproduces the ratio 39 at `H = H(W')`. So the sharp constant is attained at unboundedly large `H`, and the lemma holds with `K = 39` for every word with `H(W) ≤ 40`, of any length.
- **PROVED obstructions (§1.4).** Two natural certificates for the lemma fail, so any proof must use correlations between profile coordinates: (i) the reachable span of the 119-state transfer has dimension 117 (annihilator of dimension 2) and the LP `h A_ℓ − ψ ≤ K h` with `ψ` in the annihilator is infeasible for every `K`; (ii) the coordinate-cap cone `{x ≥ 0 : x_s ≤ c_s H(x)}` with the observed maxima `c_s` is not letter-invariant (478 violations). The lemma is not proved here.
- **FINITE-EXACT (§2).** The seam count formula (★) of `PRODUCT.md`, `H = (2^k − b)(2^r + a)·2^t + 1 + rev_k(b) − rev_r(a) − J(b) − J(a)`, re-verified against brute-force Hamiltonian-path counts of 1,922 literal tournaments (`k, r ≤ 5`, `t ≤ 2`) including both packet controls; the longest run of consecutive odd remainders at a fixed product is 3 for all `k, r ≤ 11` (the packet's bound 22 on `k + r` extended); the family's value set has density 0.42–0.46 among odd numbers up to `2^20` (`k, r ≤ 11`, `t ≤ 12`); the decoder run in reverse lists every realisation of a target (1,623,871 has 11; 262,145 has 75; 99,999 and 1,000,001 have none).
- **FINITE-EXACT (§3), the new result.** Over transitive spines with two apexes, every odd `N < 2^16` other than 7 and 21 is attained: with spine ≤ 11 over all `2^23` tapes (first miss above 21 is 77,533), and with spine ≤ 12 when the second apex beats at most three spine vertices (2,449,408 tapes; first miss 117,345). One apex alone misses half of the odd numbers; a second apex beating one spine vertex brings the coverage below `2^14` to 56 %, two vertices to 98 %, three to 100 %; below `2^16` only 0.24 % of the values need overlap three. The orientation `A → B` suffices for every value below `2^16`. Spine 13 is the pre-registered second half of Q4.
- **OPEN (HYP-9803).** (a) Two-apex completeness: every odd `N ∉ {7, 21}` is the Hamiltonian-path count of a transitive spine with two apexes, the second beating at most three spine vertices. (b) The bounded-growth lemma in its sharp form `Z_ℓ(W) ≤ 38·H(W)`. (c) The rank-potential certificate (§4): designed, not built; the finiteness of the image orbit it would need is the pre-registered Q7.

Notation. Literal span-3 word: a first vertex, then letters `ℓ ⊆ {1, 2, 3}`, the new vertex beating the older vertices at the gaps in `ℓ` and losing to every other older vertex (THM-4880; `ribbon.py`'s `transfer_system(3)`: 119 profile states, 14 complete). `H(W)` = number of Hamiltonian paths. `A_ℓ` = the letter matrices, `x_W = A_W v₀` the state, `h` the complete-state indicator, `H(W) = h·x_W`.

## 1. The bounded-growth lemma

### 1.1 The exact decomposition (PROVED)

Let `T = T_W` on vertices `1..q` and `v = q + 1` the new vertex with letter `ℓ`; `G = {q + 1 − g : g ∈ ℓ}` is the set of older vertices `v` beats (its out-neighbours; every other vertex beats `v`). Every Hamiltonian path of `T + v` has `v` last, first, or between two vertices `w, u` with `w → v → u`; so `u ∈ G`. Deleting `v`:
- `v` last: a Hamiltonian path of `T` ending at an in-neighbour of `v` — `e_ℓ(W) ≤ H(W)` of them;
- `v` first: a Hamiltonian path of `T` starting in `G` — `s_ℓ(W) ≤ H(W)`;
- `v` between `w` and `u`: if `w → u` the remainder is a Hamiltonian path of `T` in which `u` is preceded by a non-member of `G`; if `u → w` it is a two-path cover `(F₁, F₂)` with `tail(F₁) = w ∉ G`, `head(F₂) = u ∈ G`. In the profile transfer both cases are the bridging transitions of two-fragment states, and their total is `Z_ℓ(W) = Σ_s x_s·b_ℓ(s)` over the two-fragment states `s`, `b_ℓ(s)` the number of bridges `ℓ` offers `s`.

Hence `H(Wℓ) = e_ℓ + s_ℓ + Z_ℓ ≤ 2 H(W) + Z_ℓ(W)`, and `H(Wℓ) ≤ K H(W)` for all `W, ℓ` iff `Z_ℓ(W) ≤ (K − 2) H(W)` up to the two end terms. Since `u ∈ G` lies among the last three vertices and `u → w` forces `w` into the last six, the bridged covers have both junction vertices in the last six vertices; the difficulty is that the fragments themselves are unconstrained.

### 1.2 Exact maxima (FINITE-EXACT; `scripts/hratio.c`, `hratio2.c`)

Enumerating all `8^n` words by depth-first search on the 119-dimensional integer states (`scripts/ribbon_setup.py` exports the matrices to `data/mats3.txt`; the counts were checked against brute force on random words):

| length `n` | words | `max_{W,ℓ} H(Wℓ)/H(W)` | maximiser `W`, `ℓ` | `max Z_ℓ/H` |
|---|---|---|---|---|
| 1 | 8 | 3 = 3/1 | `0`, {2} | 1 |
| 2 | 64 | 5 | `00`, {3} | 3 |
| 3 | 512 | 11 | `0 3 6`, {1,2} | 9 |
| 4 | 4,096 | 23 | `0 3 6 0`, {2,3} | 21 |
| 5 | 32,768 | 39 = 39/1 | `0 0 7 6 0`, {2,3} | 37 |
| 6 | 262,144 | 39 | `0 0 0 7 6 0`, {2,3} | 38 |
| 7 | 2,097,152 | 39 | `0 0 0 0 7 6 0`, {2,3} | 38 |

The maximiser is a transitive run, then a vertex beating the last three, a vertex beating the second and third last, a sink; it has a unique Hamiltonian path and 38 two-covers that the next vertex (beating the sink's two predecessors) bridges. Lengths 8 and 9: §5.

### 1.3 The lemma is exact for every word with `H(W) ≤ 40`, and the constant is attained at every such `H` (FINITE-EXACT; `scripts/bounded_states.py`)

Breadth-first search on the exact 119-dimensional states, keeping those with `H ≤ 40`, closes after 13 letters with 16,794 states: every word of any length with `H(W) ≤ 40` has one of these states. Over them the maximum of `H(Wℓ)/H(W)` is 39 and the maximum of `Z_ℓ/H` is 38, and both are attained at **every** attained `H ≤ 40` (1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39): `env[H] = 39H`, `max Z = 38H`. The reason is sealing: after three empty letters every earlier vertex beats every later one, so `H(W' 0 0 0 X) = H(W')·H(0 0 0 X)` for any continuation `X` (the Hamiltonian paths concatenate), and `X = 7 6 0` with `ℓ = {2, 3}` gives the ratio 39 on top of any `W'`. Hence the sharp lemma `Z_ℓ(W) ≤ 38 H(W)` is tight at unboundedly large `H`, and the apparent decrease of `env[H]/H` over words of length ≤ 7 (receipt `hratio2_7`) is a bounded-length artifact. The range `H ≤ 40` is extended in §5.

### 1.4 Why two cheap certificates fail (PROVED obstructions; `scripts/ribbon_setup.py`, `hratio3.c`)

A certificate for `h A_ℓ x ≤ K h x` on the reachable cone must exclude the pure two-fragment generators `e_s`, where `h e_s = 0 < (h A_ℓ)_s` (HYP-9802's readout obstruction).
- *Annihilator.* The reachable span from `v₀` has dimension 117 of 119. The LP "`h A_ℓ − ψ_ℓ ≤ K h` coordinatewise with `ψ_ℓ` in the 2-dimensional annihilator" is infeasible for every `K`. So span and positivity alone do not see the lemma.
- *Coordinate caps.* With `c_s = max x_s/H(x)` over words of length ≤ 7 (largest 14, at a two-fragment state), the cone `P = {x ≥ 0 : x_s ≤ c_s·h x}` contains every reachable state of length ≤ 7 but is not invariant: `max_{x ∈ P, hx = 1} (A_ℓ x)_s − c_s h A_ℓ x` is positive for 478 pairs `(ℓ, s)` (largest excess 56). The reachable cone is thinner than any product of coordinate bounds; a proof must relate two-fragment weights to complete weights through the fragments' ages (a young head is fed within three letters or the state dies, unless it is prepended; a blocked tail frees itself within three letters) — the case analysis is not carried out here.

## 2. The seam product decoder in reverse (`scripts/seam_decoder.py`, `seam_spectrum.py`)

Family of `PRODUCT.md`: tape `TT_k → TT_t → TT_r`, sentinel, head, left marks `b` (backward, `0 ≤ b < 2^(k−1)`), right marks `a` (forward, `0 ≤ a < 2^(r−1)`); `H = P·2^t + C`, `P = (2^k − b)(2^r + a)`, `C = 1 + rev_k(b) − rev_r(a) − J(b) − J(a)`.
- **Verification.** Brute-force Hamiltonian-path counts of the literal tournaments agree with (★) for all 1,922 cases with `k, r ≤ 5`, `t ≤ 2` (at most 14 vertices), and the controls `(k, r, t, a, b) = (5, 4, 1, 5, 13), (5, 4, 1, 3, 11)` give 801, 803 as recorded (receipt `seam_verify`).
- **Windows.** For every `k, r ≤ 11` and every product `P`, the set of remainders `C` over the admissible factor pairs `(2^k − b, 2^r + a)` contains at most three consecutive odd values; the record 3 is the packet's `P = 812,130` at `k = r = 10` (`C = −391, −389, −387`), and the maximum number of factor pairs sharing a product is 20 (`k = r = 11`) (receipt `seam_windows`). So within this family a target window longer than three consecutive odds does not exist up to these widths.
- **Spectrum.** Forward enumeration (`k, r ≤ 11`, `t ≤ 12`) attains 222,796 of the 524,288 odd `N < 2^20`; the density is 40.6 % below 256 and 42.3 % in `[2^18, 2^20)`, stable (receipt `seam_spectrum`). The misses below 400 include 1–25 odd, 29, 35, 37, 41, 43, 47, 49, 55, …: the family is a "half" family like the one-apex family.
- **Decoding a target.** For odd `N`, every `(k, r, t, a, b)` with `P·2^t + C = N` is found by enumerating `t`, the window `P ∈ [(N − C_max)/2^t, (N + C_max)/2^t]` and the divisor pairs: 1,623,871 has 11 realisations with `k, r ≤ 11` (three with `(k, r, t) = (7, 11, 3), (8, 11, 2)`), 262,145 = `2^18 + 1` has 75, while 99,999 and 1,000,001 have none (receipt `seam_decode`). Choosing the factor pairs "from the desired count" is therefore exact but reaches only the family's 43 %.

## 3. Two apexes over a transitive spine attain every odd number below 2^15 except 7 and 21 (`scripts/two_apex.c`, `two_apex_small.c`)

Family: spine `0 < 1 < … < m − 1` transitive; apex `A` beats spine vertex `i` iff bit `i` of `a` is set, apex `B` likewise with `b`; the arc `A → B` (orientation 0) or `B → A` (1). All `2^(2m+1)` tapes are counted by the subset DP; the one-apex family (`B` a sink: `b = 0`, orientation 0) is included, and padding a source in front shows the value sets are nested in `m`. In the OCF frame (UNIFYING.md §3) the count is `1 + 2(a + b + c) + 4d` with `a, b, c` the odd cycles through `A` only, `B` only, both, and `d` the disjoint pairs.

| spine `m` | tapes | distinct values | max `H` | odd `N` below the maximum attained | first misses |
|---|---|---|---|---|---|
| 4 | 512 | 19 | 45 | — | — |
| 6 | 8,192 | 202 | 509 | — | — |
| 8 | 131,072 | 1,967 | 6,069 | — | — |
| 10 | 2,097,152 | 19,876 | 69,485 | every odd `N < 2^14` except 7, 21 | 7, 21, 27141, 27639, … |
| 11, second apex beating ≤ 3 spine vertices | 950,272 | 48,963 | 194,529 | every odd `N < 2^15` except 7, 21 | 7, 21, 50023, 51673, … |
| 11, all tapes | 8,388,608 | 63,113 | 227,295 | every odd `N < 2^16` except 7, 21 | 7, 21, 77533 |
| 12, second apex beating ≤ 3 spine vertices | 2,449,408 | 141,026 | 611,117 | every odd `N < 2^16` except 7, 21 | 7, 21, 117345, 128617 |

- **How much overlap is needed.** Below `2^14` (spine ≤ 10): 414 values are one-apex counts; 4,157 need a second apex beating one spine vertex; 3,462 need two; 157 need three (the smallest 8,725); none needs more. Below `2^15` (spine 11, overlap ≤ 3): 825 / 9,093 / 6,359 / 105; below `2^16` (spine 12, overlap ≤ 3): 1,605 / 19,467 / 11,615 / 79 (0.24 % need three). Over all spine-11 tapes no value below `2^15` needs overlap above three. So the "controlled overlap" is at most three vertices.
- **Orientation.** Every value below `2^16` is attained with `A → B` at spine 12 (at spine 11 only 32,507 below `2^15` needed `B → A`).
- **Least spine.** Every odd `N ≤ 33` other than 7, 21 is attained with spine 4 (`N ≤ 6` vertices); 35, 39, 47–59 need spine 5; 61, 63, 87, 99, 107, 113, 119 spine 6 (receipt `two_apex_analysis`).
- **Reading.** One apex gives the cliques `K_c` (`N = 1 + 2c`) for about half of all `c`; the second apex, beating at most three spine vertices, supplies every missing residue — in fact every missing value — below `2^16`. The lemma to aim for (UNIFYING.md §3) is confirmed in its strongest finite form, and HYP-9803 (a) states it for all `N`.

## 4. The rank-potential certificate: design and the open step

HYP-9802 (iii) / LITERAL_CONES.md §4: automaton states `(short suffix tag, rank of the prefix product)`, cones = reachable cones of products of that rank, an LP charging a bounded loss per rank drop and requiring growth `ρ` between drops; with the bounded-growth lemma replacing the readout cap. What this session adds: (a) the exact decomposition §1.1 says which quantity the lemma must control (`Z_ℓ`); (b) the two obstructions §1.4 say the lemma is not a consequence of span and positivity, so the certificate cannot import it as a linear side constraint on raw generators — it must be proved separately and then used as `|Φ| ≤ C·H` on actual states; (c) the exactness of the rank transition relation needs the image subspaces of prefix products (LITERAL_CONES.md's open step): whether their orbit under the eight letters is finite is the pre-registered Q7 (§5). The certificate is **not built**; wave 39's growth-only LP (6,262 states) did not terminate at `ρ = 1.2`, and no smaller sound automaton is known to cut the slow cycles of period 7–12.

## 5. Pre-registered tests (PREREGISTERED.md, committed 18:52 before any output was read)

| # | prediction | verdict |
|---|---|---|
| Q1 | max `H(Wℓ)/H(W)` = 39 at lengths 8 and 9, maximiser `0^(n−3) 7 6 0`, `ℓ = 6` | **length 8 confirmed** (16,777,216 words; receipt `hratio3_8`); length 9: PENDING_L9 |
| Q2 | max `Z_ℓ/H` = 38 at lengths 8, 9 | length 8: implied by Q1's maximiser (`H = 1`, 39 = 1 + 0 + 38); length 9: PENDING_L9Z |
| Q3 | all spine-11 tapes attain every odd `N < 2^15` but 7, 21; first miss in `[45,000, 70,000]` | **confirmed and exceeded** (every odd `N < 2^16` but 7, 21); the first-miss guess **refuted**: 77,533 |
| Q4 | overlap ≤ 3: spine 12 covers every odd `N < 2^16` but 7, 21; spine 13 every odd `N < 2^17` but 7, 21 | **spine 12 confirmed** (first miss 117,345); spine 13: PENDING_S13 |
| Q5 | values needing overlap exactly 3 form 0.1–2 % below `2^16`; overlap ≤ 3 suffices | **confirmed**: 79 of 32,766 (0.24 %); over all spine-11 tapes no value below `2^15` needs overlap above 3 |
| Q6 | `A → B` suffices for all but < 10 values below `2^16` | **confirmed**: for all of them |
| Q7 | the image orbit is not finite within the cap | **confirmed**: 9, 61, 381, 2,213, 10,337, 50,329 distinct image subspaces at depths 1–6 (119-state representation, modulo a 31-bit prime), about fivefold growth per letter, ranks from 1 to 119 (receipt `image_orbit`) — an exact rank-tagged automaton on image subspaces is infinite, so the rank potential must be coarsened |

The closed state sets of §1.3 were computed after the pre-registration and were not predicted.

## 6. The owner's Collatz framing

Recorded as framing, not as a result: exact constructions here are rigid like Collatz cycle equations (exponential Diophantine in 2 and 3: the seam family's `P·2^t + C`, the C3 register's base-four digits); the statistical layer (density-one coverage by positive semigroups, UNIFYING_FRAME.md) plays the role of "almost all orbits"; both problems lack the passage from almost all to every target, and the carry of Collatz corresponds to the cross term `4d` (the disjoint-pair count) in the two-apex formula. Nothing in this packet proves any of these correspondences.

## 7. Files

- [`scripts/ribbon_setup.py`](scripts/ribbon_setup.py) (transfer export, annihilator LP), [`scripts/hratio.c`](scripts/hratio.c) / [`hratio2.c`](scripts/hratio2.c) / [`hratio3.c`](scripts/hratio3.c) (maxima, two-fragment split and envelope, coordinate maxima), [`scripts/bounded_states.py`](scripts/bounded_states.py) (closed state set for `H ≤ M`), [`scripts/image_orbit.py`](scripts/image_orbit.py), [`scripts/seam_decoder.py`](scripts/seam_decoder.py) (verify / decode / windows), [`scripts/seam_spectrum.py`](scripts/seam_spectrum.py), [`scripts/two_apex.c`](scripts/two_apex.c), [`scripts/two_apex_small.c`](scripts/two_apex_small.c).
- [`data/`](data/): `mats3.txt` (the eight 119×119 letter matrices, `v₀`, `h`), `frag.txt` (fragment counts of the states), `cmax_7.txt`, the two-apex value lists `tapex_<m>.txt` (`V h count minpop orient`), `tapexs_11.txt(.gz)`.
- [`receipts/`](receipts/): run outputs named in the text; the Python receipts are byte-identical under `python3` and `python3 -O`.
