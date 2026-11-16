# The junction inequality: H is multiplicative up to a constant on bounded-span ribbons, and what that does to the growth certificate

Session mac-mini-2026-11-15-junction. Canon: THM-5252. **Chan–Pak remains OPEN. The literal span-3 lower-growth
certificate is not built.**

Owner prompt (2026-11-15, sent to several sessions): choose the factor pairs from the desired count so that their
exact remainders land on the target; build the span-3 rank-potential certificate and prove the bounded-growth lemma;
density zero with a power saving for the odd numbers missed by span-3 ribbons; decide whether the conflict graphs of
tournaments give every odd value except 7 and 21 once a second apex is allowed; the Collatz parallel.

Two other sessions answered the same prompt while this one ran, and what they proved is not repeated here:

- [chanpak_target_ribbon_audit_20261116](../chanpak_target_ribbon_audit_20261116/README.md): the bounded-growth
  lemma with the **sharp** constant, `H(T) <= H(T+x) <= 39 H(T)` for every span-3 prefix (tail deletion), the
  target-first product decoder and its limits, the rank alias, the density gate, the two-apex bridge
  `H = hj + PQ - 1`.
- [chanpak_eight_braid_20261116](../chanpak_eight_braid_20261116/README.md): the eight sharp one-letter factors, and
  a lower growth bound for one return language by cutting it into blocks and multiplying.

This note adds the statement those two results are special cases of, and draws its consequences for the certificate.

## 0. Results

| | statement | status |
|---|---|---|
| A | **Junction inequality.** For words `u`, `v` of span `w`: `H(u) H(v) <= H(uv) <= K_w^2 H(u) H(v)`, with a constant depending only on `w` (`K_3 = 26,425`). The left side is known; the right side is new. Largest ratio found in span 3: 3,965. | PROVED (general tournament form in §1) |
| B | Every raw coordinate of the 119-state transfer vector of a word `W` is at most `K_3 H(W)`. So any potential built from continuation readouts is bounded by a constant times `H`, for continuations of any length. | PROVED |
| C | **Floor for literal words.** A span-3 factor whose operator has rank at least 2 has a strong component missing at most 12 of its vertices, so `H >= 5^((n-13)/3)`. The bottleneck-free growth rate `rho_min` exists as a limit, equals `sup_m g_m^(1/m)`, and `5^(1/3) <= rho_min <= rho* = ((9+sqrt 113)/2)^(1/4) = 1.7700005`. | PROVED (with Busch's bound, CITED) |
| D | HYP-9742 for literal span 3 (`rho_min > 5^(1/3)`) is equivalent to each of: some `g_m > 5^(m/3)`; the tree of slow factors (§3) is finite at some rate above `5^(1/3)`. | PROVED (equivalence only) |
| E | Periodic words, exhaustive: among all primitive periods 2 to 8 the only bottleneck-free cycles with growth below 1.7705 are 30 necklaces of period 4 and 185 of period 8, all with growth exactly `rho*` (nonzero spectrum `x^2 - 9x - 8`, and its square image `x^2 - 97x + 64`). No cycle of period at most 8 beats `rho*`. | FINITE-EXACT |
| F | A certificate made of linear potentials on cones of single state vectors cannot see the condition "rank at least 2", whatever the finite memory. The slow-factor tree is the sound replacement; its first six levels have 8, 64, 504, 3,976, 31,292, 162,702 nodes at rate 1.72. | argument + FINITE-EXACT data; certificate NOT built |

## 1. The junction inequality

`H(T)` is the number of Hamiltonian paths of a tournament `T`; `d+`, `d-` are out- and in-degrees.

**Known (abundancy; Grinberg–Stanley odd-cycle form, THM-002; THM-5130 (H4)).** For every vertex set `S`,
`H(T) >= H(T[S]) H(T - S)`: packings of disjoint odd cycles on the two supports combine.

**Lemma 1 (several targets).** Let `D` be a tournament and `y_1, ..., y_r` distinct vertices. The number of tuples
`(F_0, F_1, ..., F_r)` of pairwise disjoint directed paths covering `D`, with `F_0` possibly empty and `F_i` starting
at `y_i`, is at most

    sum over pairwise disjoint S_i ⊆ N+(y_i) of  H(D - ∪S_i) · ∏ H(D[S_i])   <=   ∏ 2^(d+(y_i)) · H(D).

*Proof.* Put `R_0 = F_0`. For `i = 1, ..., r`: let `σ_i` be the longest terminal segment of `R_(i-1)` all of whose
vertices are beaten by `y_i`, `S_i` its vertex set, and `R_i = (R_(i-1) minus σ_i) · F_i`. This is a path: if
something is left of `R_(i-1)`, its last vertex is not beaten by `y_i`, so it beats `y_i`, the first vertex of
`F_i`. At the end `R_r` is a Hamiltonian path of `D - ∪S_i`, each `σ_i` is a Hamiltonian path of `D[S_i]`, the `S_i`
are pairwise disjoint and `S_i ⊆ N+(y_i)`.
The data `(S_i, σ_i, R_r)` determine the tuple: for `i = r, ..., 1`, cut `R_i` just before `y_i`; the part from
`y_i` on is `F_i`, and the part before it followed by `σ_i` is `R_(i-1)`. So the number of tuples is at most the
first sum, and each term is at most `H(D)` by abundancy. ∎

For `r = 1` this is the tail deletion of the concurrent sharp-growth proof; the point here is that it iterates.

**Theorem A (junction inequality).** Let `V(T) = A ⊔ B`. Call an arc from `B` to `A` a back arc; let `A°` be the set
of heads of back arcs and `B°` the set of their tails. Then

    H(T[A]) H(T[B])  <=  H(T)  <=  κ+(A) κ-(B) · H(T[A]) H(T[B]),
    κ+(A) = sum over R ⊆ A° of |R|! ∏_{y in R} 2^(d+_{T[A]}(y)),     κ-(B) = sum over R ⊆ B° of |R|! ∏_{z in R} 2^(d-_{T[B]}(z)).

*Proof.* Cut a Hamiltonian path `Q` at every arc joining `A` and `B`. Its maximal segments alternate between `A` and
`B`. A step from a `B`-segment to an `A`-segment uses a back arc. So `Q = [F_0] G_1 F_1 G_2 F_2 ... G_r F_r [G_(r+1)]`
with `F_0` and `G_(r+1)` possibly empty, each `F_i` (`i >= 1`) starting at a distinct vertex of `A°` and each `G_j`
(`j <= r`) ending at a distinct vertex of `B°`. The two tuples determine `Q`. By Lemma 1 the number of tuples
`(F_0; F_1, ..., F_r)` with prescribed starts is at most `∏ 2^(d+(y_i)) H(T[A])`; summing over ordered choices of the
starts gives `κ+(A) H(T[A])`. The same lemma in the converse of `T[B]` gives `κ-(B) H(T[B])`. ∎

**Corollary (bounded span).** For span-`w` words `u`, `v` (the letters of `v` may reach back into `u`; `H(v)` is
the count of the tournament induced on the vertices of `v`):

    H(u) H(v) <= H(uv) <= K_w^2 H(u) H(v),     K_w = sum over R ⊆ {0, ..., w-1} of |R|! ∏_{i in R} 2^(w+i).

Indeed `A°` lies in the last `w` vertices of `u`, and the vertex `i` places from the end has at most `w` back arcs
and `i` later out-neighbours inside `u`; symmetrically for the first `w` vertices of `v`. `K_3 = 26,425`.

- The case `|v| = 1` is the bounded-growth lemma of HYP-9802 (ii): there `r <= 1`, and the proof gives
  `1 + 8 + 16 + 32 = 57` in place of the sharp 39 proved in SHARP_GROWTH.md.
- **Corollary B.** A raw coordinate of the transfer vector of `W` counts path systems on `T_W` with one free path
  and further paths starting in the last three vertices. By Lemma 1 it is at most `K_3 H(W)`. A readout
  `H(W s)` is at most `K_3^2 H(s) H(W)` for a continuation `s` of **any** length, where the one-letter bound gives
  `39^|s| H(W)`. Observed: over all words with at most 8 vertices the largest raw coordinate is `14 H(W)`
  (`scripts/check_theorem_a.py`).
- **Observed constant** (`scripts/qm_ratio.py`, `qm_beam.py`). Over all standalone `u` with at most 6 vertices and
  all `v` with at most 4, the largest `H(uv) / (H(u) H(v))` is 39, 179, 603, 2,071 for `|v| = 1, 2, 3, 4`, and
  the smallest is 1. A beam search to `|v| = 12` saturates at **3,965** from `|v| = 5` on (two transitive pieces of
  6 and 5 vertices joined by four back arcs have 3,965 paths). The proved constant `K_3^2 ≈ 7·10^8` is far from
  sharp.

**What it gives at once.**

1. *Limits exist.* For a set of words closed under taking factors, `m -> min H` over its words of `m` vertices is
   supermultiplicative, so its `m`-th root converges to its supremum. `m -> max H` over all words satisfies both
   inequalities, so `max H = λ^(m + O(1))` with `λ = sup_m (max_m)^(1/m)`.
2. *Pressure with bounded defect.* `Z_m(s) = sum over standalone words of m vertices of H^s` satisfies
   `64 Z_m Z_m' <= Z_(m+m')(s) <= 64 K_3^(2s) Z_m Z_m'` for `s >= 0` and `m, m' >= 3` (reversed for `s < 0`; 64
   counts the letters of `v` that reach into `u`). So `log Z_m(s) = m P(s) + O(1)`: `log H` is an almost additive potential in the
   sense of the thermodynamic formalism, and the word counts by size of `H` are governed by `P`. This is the
   hypothesis under which Zaremba-type counting arguments start; it says nothing yet about distinct values (the
   density gate of the concurrent packet).

## 2. The floor for literal words, and `rho_min`

Conventions: a span-3 word is a first vertex followed by letters; a *factor* `u` is a string of consecutive
letters with its operator `A_u` in the minimal representation (dimension 61, after the first vertex); `T[u]` is
the tournament induced on its vertices. A word is bottleneck-free when its operator has rank at least 2, which
holds exactly when every factor has rank at least 2 (the rank of a product is at most the rank of each factor).

**Lemma 2 (an interior cut is a bottleneck).** Let `u` be a factor with `m` vertices. If `T[u] = A ⇒ B` is an
ordinal sum with `7 <= max A <= m - 3`, then `rank A_u <= 1`.

*Proof.* Every vertex of `B` is at least `max A - 3 >= 4`, and every vertex of `A` is at most `m - 3`. So in any
word `p u s` no back arc joins `B` to `p` or `s` to `A`, and `T[p u s] = (p ∪ A) ⇒ (B ∪ s)`. Hence
`H(p u s) = H(T[p ∪ A]) · H(T[B ∪ s])` is a product of a function of `p` and a function of `s`. The prefix states
span the minimal representation and the suffix readouts span its dual, so the operator has rank at most 1. ∎

**Theorem C.** If `rank A_u >= 2` and `u` has `m` vertices, then `T[u]` has a strong component with at least
`m - 12` vertices, and `H(T[u]) >= 5^((m-13)/3)`.

*Proof.* Let `C_1 ⇒ ... ⇒ C_r` be the strong components and `A_i = C_1 ∪ ... ∪ C_i`. Each `A_i` with `i < r` is a
cut, so `max A_i <= 6` or `max A_i >= m - 2`. A vertex more than 3 below `max A_i` beats it, hence lies in `A_i`.
Let `i` be the last index with `max A_i <= 6` (`A_0` empty). Either `C_(i+1)` is the last component, with at least
`m - 6` vertices, or `|A_(i+1)| >= m - 6` and `C_(i+1) = A_(i+1) - A_i` has at least `m - 12`. The count of a
strong tournament on `k` vertices is at least `5^((k-1)/3)` (Busch 2006, as cited in THM-4810). ∎

A check on 1,500 random factors of 11 to 16 vertices (`scripts/cutcheck.py`): all 416 with an interior cut have
rank at most 1; among the 665 of rank at least 2 the largest strong component never misses more than 6 vertices.

**Definition.** `g_m` is the least `H(T[u])` over factors `u` of `m` vertices with `rank A_u >= 2`.
Its first values (`scripts/gmin.py`, ranks modulo a prime): `g_m = 1` for `m <= 8` (there are 4,384 factors of 8
vertices whose induced tournament is transitive and whose operator has rank at least 2, and none of 9),
`g_9 = g_10 = g_11 = 3`, `g_12 = 5`, `g_13 >= 9`. So a bottleneck-free factor can start a factor 190 below
`rho*^m` (at `m = 12`).

**Corollary D.**

1. `g_(m+m') >= g_m g_m'`, so `rho_min := lim g_m^(1/m) = sup_m g_m^(1/m)` exists, and every bottleneck-free word
   with `n` vertices has `H >= g_m^floor((n-1)/m)` for every `m`.
2. `5^(1/3) <= rho_min <= rho*`. The lower bound is Theorem C. The upper bound is the periodic word `(0,2,1,2)^k`
   (letters as bit masks of the back arcs), whose operator has stable rank 2 and spectral radius `rho*^4`.
3. The following are equivalent: (i) `rho_min > 5^(1/3)`, which is HYP-9742 for literal span 3; (ii) `g_m > 5^(m/3)`
   for one `m`; (iii) for some `rho > 5^(1/3)` there are only finitely many factors `u` with `rank A_u >= 2` and
   `H(T[u']) < rho^|u'|` for every prefix `u'` of `u`.

*Proof of 3.* (i) ⟺ (ii) by the supremum formula. (iii) ⟹ (i): cut a bottleneck-free word greedily into shortest
prefixes with `H(T[u']) >= rho^|u'|`; each exists within the depth `L` of the finite tree, the pieces are factors,
and abundancy gives `H >= rho^(n-L-1)`. (i) ⟹ (iii): for `rho < rho_min`, `g_m >= rho^m` from some `m` on. ∎

So the non-strict form of the literal growth bound is a theorem, the strict form is one finite statement, and no
potential or cone is needed to state it.

## 3. Why the certificate is still not built

**Linear potentials cannot see the rank.** A cone certificate attaches to each memory state a cone of state
vectors and a linear potential, and asks for growth on every generator along every allowed edge. A slow cycle `c`
(growth below the target) contracts every single vector to its leading direction after boundedly many periods, and
from then on that vector grows at the rate of `c`. What forbids `c^k` for large `k` is that the *operator* drops
to rank 1, which is a statement about pairs of vectors, not about any one of them. The set of operators of rank at
least 2 is not convex, and the cone generated by the reachable operators contains rank-one generators. So the
inequalities fail on generators unless the memory itself is long enough to have removed the edge, and the memory
needed is the bottleneck length of the slowest cycle below the target rate. The rank alias of the concurrent
packet (same suffix, same rank, different next rank) is the finite form of the same fact.

**The sound replacement is the slow-factor tree** of Corollary D (iii): the block being built is its own memory,
its exact operator rank is tested, and abundancy allows a restart with empty memory after each good block.
`scripts/badtree.py` builds it level by level (ranks modulo a prime; a full certificate would re-check each discard
exactly). At rate 1.72:

| vertices | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| candidates | 8 | 64 | 512 | 4,032 | 31,808 | 250,336 |
| rank at most 1 | 0 | 0 | 8 | 56 | 516 | 4,348 |
| already good | 0 | 0 | 0 | 0 | 0 | 83,286 |
| slow nodes kept | 8 | 64 | 504 | 3,976 | 31,292 | 162,702 |

Two things make the tree large, and both are now quantified.

- *The generic bulk.* A typical word (growth 2.24 per letter) needs about 10 letters to get above `1.72^i`, so the
  tree first grows by nearly 8 per level.
- *Slow starts.* A factor may begin with a stretch of count 1 to 3 over 8 to 12 letters and still have rank at
  least 2. Below such a node the tree contains every continuation whose count stays below `e^D · rho^i`, where
  `e^D` is the deficit of the start. Its depth is about `D / log(rho_min / rho)`: hundreds of letters for any rate
  between `5^(1/3)` and `rho*`.

**A design that should work (not built).** Two nodes with the same column space of the operator have the same
future ranks. If in addition the standalone state vector of one is coordinatewise at most that of the other, every
continuation of the first has the smaller count. So only the minimal state vectors need to be kept for each column
space. This is the dominance pruning that made the span-2 level sets tractable (THM-4880). The slow starts are
exactly the nodes with the smallest vectors, so they would dominate the bulk instead of multiplying it.

## 4. Periodic words to period 8 (exhaustive)

`scripts/cheap_periodic.py`: all primitive necklaces of a given period over the 8 letters, pruned exactly by
abundancy (a periodic word with growth `g` per letter has `g^P >= H(T[v])` for every prefix `v` of its period, and
three consecutive empty letters make a cut). For each survivor the growth is the spectral radius of the
minimal-representation product (floating point), and the stable rank of the operator is computed exactly
(`python-flint`, the rank of `O M^64 Rb`).

| period | necklaces after pruning | growth below 1.7705 | of these with stable rank at least 2 |
|---|---:|---:|---:|
| 2 | 28 | 7 | 0 |
| 3 | 168 | 20 | 0 |
| 4 | 1,001 | 117 | 30 |
| 5 | 6,496 | 576 | 0 |
| 6 | 26,247 | 2,944 | 0 |
| 7 | 105,191 | 16,306 | 0 |
| 8 | 524,596 | 97,391 | 185 |

The 30 necklaces of period 4 all have growth `rho*`; the nonzero part of the characteristic polynomial of their raw
119-state product is `(x - 1)(x^2 - 9x - 8)`, and their stable rank is 2. They include `(0,2,1,2)`, the block
found in wave 28. The 185 necklaces of period 8 also have growth exactly `rho*` and stable rank 2; their nonzero
spectrum is `(x - 1)(x^2 - 97x + 64)`, the squares of the period-4 roots, so they are two period-4 minimising
blocks in a row (primitive as words of period 8). This extends the exhaustive range of HYP-9802 (i) (period 6,
growth below 1.72) and the code-unit census of LITERAL_CONES.md to every primitive word of period at most 8, at
the critical rate itself. It supports `rho_min = rho*`; it does not prove it, because the slow starts of §3
belong to longer periods.

## 5. The other parts of the prompt

Nothing new is proved here on these; the present state, with credit, is as follows.

- **Factor pairs.** The exact remainder of a one-apex tape is the sum of the checksums of its two factors:
  `p(A·2^s + rc_s(t)) = 2·ceil(At/2) + p(A) + p(t)` for all `A >= 0`, `0 <= t < 2^s` (`rc_s` the reversed
  complement on `s` bits). This is the split identity of UNIFYING.md with THM-5232's invariance substituted; it
  makes the one-apex spectrum the set `1 + 2·ceil(At/2) + p(A) + p(t)` and shows the problem is self-similar: the
  remainder that must "land on the target" is again a pair of values of the same function. The target-first
  decoder and its limits are in PRODUCT_INVERSION.md.
- **Second apex.** The clique case is closed and has been since October: tournaments with no two disjoint
  triangles have exactly the one-apex spectrum and two dyadic rays
  ([intercyclic_rays.md](../feedback_word_hamilton_spectrum_20261020/intercyclic_rays.md),
  [INTERCYCLIC_CLASSIFICATION.md](../session336_overlap_inverse_and_intercyclic_20261114/INTERCYCLIC_CLASSIFICATION.md)).
  So complete conflict graphs give the one-apex half and `O(log X)` further values, and the missing half needs
  disjoint cycles: `N = 1 + 2c + 4d` with `d >= 1` (TWO_APEX.md §5; the bridge there is the sub-case `d = ab`).
  Whether the tournaments with at most two disjoint triangles already give every odd number except 7 and 21 is
  open; the finite evidence is the two-apex census of opus S93 (to `2^26`) and the carrier census of THM-5192.
- **Density.** Not established here or elsewhere. §1 gives the almost-additive potential; the collision estimate of
  DENSITY_GATE.md is the missing input.
- **Collatz.** ANALOGY, with one exact point. The checksum is an additive cocycle over the binary shift,
  `(w, p) -> (2w + b, p + (1 - b)(w + (w mod 2)))`, as the Collatz remainder is an additive cocycle over the parity
  shift. In both the cocycle is understood modulo every fixed modulus (every admissible residue of the checksum
  is attained, Session25), and in both the open step is to prescribe its exact value: the cross term `A t` above
  is the carry.

## 6. Reproduction

All scripts read `lit_base.pkl`, written by `../chanpak_waves_20261113/scripts/w39_lit_setup.py` (run it in the
working directory first; it needs `04-computation/experiments/ribbon_spectrum_20261104/ribbon.py`).

    python3 scripts/check_theorem_a.py       # Lemma 1 and Theorem A on random tournaments; Corollary B to 8 vertices
    python3 scripts/qm_ratio.py 6 4          # all short u, v: ratios between 1 and 2071
    python3 scripts/qm_beam.py 6 12 40       # saturation at 3965
    python3 scripts/cutcheck.py              # Lemma 2 and Theorem C on random factors
    python3 scripts/gmin.py 5 30 120000      # g_m for m <= 12 (needs badtree.py beside it)
    python3 scripts/cheap_periodic.py 8 1.7705   # one period per run; period 8 takes about 30 minutes
    python3 scripts/badtree.py 1.72 60000 30     # the slow-factor tree, first levels

## 7. Not claimed

- Chan–Pak; density one or positive density of the ribbon spectrum; `rho_min > 5^(1/3)`; `rho_min = rho*`.
- The sharp one-letter constant 39 and the eight sharp mask factors (other sessions, same day).
- Sharpness of `K_3^2`; the true constant is at least 3,965 and probably of that order.
- Exactness of the growth values in §4 (floating point; the ranks and the pruning are exact).
