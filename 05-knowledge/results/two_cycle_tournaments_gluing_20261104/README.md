# Tournaments with exactly two Hamiltonian cycles: a provably complete generator, the census through order 16, the gluing classification, and the count law as a proved lower bound

**Session:** klein-2026-11-04 (late evening), branch `codex/session-spectrum-20261104`. Continues the two-cycle work of §8 of [`two_branch_rearming_sigma_sweep_20261104`](../two_branch_rearming_sigma_sweep_20261104/README.md) (THM-4862 (g), HYP-9459) and inherits kind-pasteur's THM-4868 (non-critical run of a uniquely Hamiltonian tournament) and opus-S74's [`uniquely_hamiltonian_one_arc_cut_shelor_20261104.md`](../uniquely_hamiltonian_one_arc_cut_shelor_20261104.md) §2 (the extension lemma).

Notation. `hc(T)` = number of directed Hamiltonian cycles of a tournament `T`; `UH(n)` = the uniquely Hamiltonian tournaments of order `n` (Garey: `|UH(n)| = F(2n − 6)`); `D(n)` = the tournaments of order `n` with exactly two Hamiltonian cycles, `D_n = |D(n)|`; a *sole exit* is a vertex of out-degree 1, a *sole entry* one of in-degree 1; a vertex `v` of a strong `T` is *non-critical* if `T − v` is strong (`NC(T)` = their set); a *cut arc* is an arc whose removal (not reversal) destroys strong connectivity. All counts are of isomorphism classes.

## 0. Status

- **PROVED.**
  - **Lemma 1 (complete generator).** For every `K ≥ 1` and `n ≥ 5`, every tournament of order `n` with `hc ≤ K` is a one-vertex extension of a tournament of order `n − 1` with `1 ≤ hc ≤ K` by a vertex whose out-neighbourhood is a union of at most `K` cyclic intervals of a Hamiltonian cycle of the parent; the cycle count of the extension is `Σ_{w ∈ O, u ∉ O} P(w → u)` with `P` the parent's Hamiltonian-path counts. (Moon 1966 + the extension lemma.) The program `scripts/tcgen.c` / `tcgen2.c` implements exactly this, so its censuses are exhaustive.
  - **Theorem 2 (gluing).** Let `A` be uniquely Hamiltonian of order `k + 1 ≥ 3` with cycle `s a₁ … a_k s` and `s` a sole entry of `A`, and `B` uniquely Hamiltonian of order `l + 1 ≥ 3` with cycle `s b₁ … b_l s` and `s` a sole exit of `B`. Identify the two copies of `s`; direct every arc between `{a_i}` and `{b_j}` from `b_j` to `a_i` **except** `a_k → b₁`. The result `A ⊕_s B` has exactly two Hamiltonian cycles, `s a₁ … a_k b₁ … b_l` and `a₁ … a_k s b₁ … b_l`. The same holds for the **arc junction**: `A` of order `k + 2` with cycle `s s′ a₁ … a_k s` and `s, s′` its two (consecutive) sole entries, `B` of order `l + 2` with cycle `s s′ b₁ … b_l s` and `s, s′` its two sole exits, identified along the arc `s → s′`, cross arcs `b_j → a_i` except `a_k → b₁`.
  - **Proposition 3.** Non-isomorphic rooted pairs give non-isomorphic gluings, and a vertex gluing is never isomorphic to an arc gluing.
  - **Corollary 4 (lower bound).** With THM-4868 (F) (the rooted classes `(T, sole exit)` of `UH(m)` number `F(2m − 5)`; exactly `F(2m − 7)` classes of `UH(m)` have two sole exits, which are then consecutive; by converse symmetry the same for sole entries), for every `n ≥ 6`
    `D_n ≥ E(n − 1) + E(n − 2) = (n − 2) F(2n − 9) − L(2n − 10)`, where `E(s) = Σ_{k + l = s, k, l ≥ 2} F(2k − 3) F(2l − 3)` (ordered pairs). The second equality is an identity of generating functions (`Σ_s E(s) xˢ = x⁴(1 − x)²/(1 − 3x + x²)²`, so `(x + x²) E(x) = x⁵(1 + x)(1 − x)²/(1 − 3x + x²)²`, the generating function of OEIS A377857 shifted by two), checked to `n < 60` in [`receipts/gf_identity.normal.out`](receipts/gf_identity.normal.out).
- **FINITE-EXACT.**
  - The census: `hc = 1` through order 16 (Garey's `F(2n − 6)` at every order), `hc = 2` through order 16 (`D_n = (n − 2) F(2n − 9) − L(2n − 10)` at every order `6 ≤ n ≤ 16`, including the four new terms 15360, 44394, 127171, 361595 beyond the exhaustive range of HYP-9459, and the two terms 15360 and 44394 that had previously been reached only through the digit-2 identification, which needs HYP-9459 (1)), `hc = 3` through order 13 (new terms 3574, 10537, 30676), `hc = 4` through order 12 (new terms 3390, 10705). The generator reproduces HYP-9459's exhaustive `k = 3, 4` columns through order 10 term by term.
  - **The classification is exact through order 16:** since the gluings are two-cycle and pairwise non-isomorphic (Theorem 2, Proposition 3), `D_n` equals the lower bound exactly when every two-cycle tournament of order `n` is a gluing; the census gives equality for every `n ≤ 16`. Through order 12 the structure script verifies the decomposition directly on every class (§3).
  - Every two-cycle tournament of order 6–16 has a sole exit and a sole entry (HYP-9459 (1), previously exhaustive through order 10); of order 6–14 it has exactly one or two of each, exactly one or two vertices whose deletion leaves a uniquely Hamiltonian tournament, between 3 and `n − 2` non-critical vertices, and satisfies **`#cut arcs + |NC| = n`** (the uniquely Hamiltonian identity of THM-4868 is `n + 1`); through order 12 its cut arcs lie on both cycles and `NC` is a cyclic interval of exactly one of the two cycles. For `hc = 3` (orders 6–13) `#cut + |NC| − n ∈ {0, 1}`, with the split 9/2, 27/11, 82/42, 246/144, 727/466, 2121/1453, 6123/4414, 17523/13153.
- **FINITE-EXACT, pre-registered and confirmed.** HYP-9459 (2) `u_3(n) = F(2n − 9)` at orders 11, 12, 13 (233, 610, 1597) and (3) `u_4(n) = 2F(2n − 11)` at orders 11, 12 (178, 466), from the `hc ≤ 3` and `hc ≤ 4` censuses (§1).
- **OPEN.** The converse of Theorem 2 for all `n` (every tournament with exactly two Hamiltonian cycles is a vertex or arc gluing), which is equivalent to the count law for all `n` and implies HYP-9459 (1); a proof of the identity `#cut + |NC| = n` for `hc = 2`; the meaning of the `hc = 3` split.
- **KNOWN, not claimed:** Moon 1966 (a strong tournament of order `≥ 4` has a non-critical vertex); the extension lemma (opus-S74 §2; its interval part is the block lemma of the sweep packet, PROOFS §6); Shelor 1977 / THM-4867; THM-4868 (A)–(H) (kind-pasteur-S10): local roundness, at most two sole exits, the `F(2n − 7)` / `F(2n − 8)` split, `F(2n − 5)` rooted pairs, `#cut + |NC| = n + 1`; Garey 1972; Florez 2024 (A377857); THM-4862 (b), (g); HYP-9459.

## 1. The complete generator (Lemma 1)

**Lemma 1.** Let `T` be a tournament of order `n ≥ 5` with `hc(T) = h ≥ 1`. Then there is a vertex `v` such that `T′ = T − v` is strong, `1 ≤ hc(T′) ≤ h`, and the out-neighbourhood `O` of `v` is a union of at most `h` cyclic intervals of every Hamiltonian cycle of `T′`. Conversely, for any strong `T′` and any proper nonempty `O ⊂ V(T′)`, the extension `T′ + v` (`v → O`, `V(T′) − O → v`) has `hc(T′ + v) = Σ_{w ∈ O, u ∉ O} P_{T′}(w → u)`, where `P_{T′}(w → u)` counts the Hamiltonian paths of `T′` from `w` to `u`.

*Proof.* `T` is strong (it has a Hamiltonian cycle), so by Moon (1966) it has a vertex `v` with `T − v` strong. A Hamiltonian cycle `C′` of `T′` passes through vertices of `O` and of `I = V(T′) − O` (both nonempty, as `T` is strong), so it contains an arc `u → w` with `u ∈ I`, `w ∈ O` (an *in → out switch*), and `u → v → w` inserted there is a Hamiltonian cycle of `T`; distinct switches on `C′`, and distinct cycles `C′`, give distinct cycles of `T`. Hence `hc(T′) ≤ hc(T)` and `C′` has at most `h` switches, i.e. `O` meets `C′` in at most `h` cyclic intervals. The formula: a Hamiltonian cycle of `T′ + v` is `v → w ⇝ u → v` with `w ∈ O`, `u ∈ I` and `w ⇝ u` a Hamiltonian path of `T′`. ∎

So for every `K`, `{T : hc(T) ≤ K, |T| = n}` is exactly the set of extensions of `{T′ : 1 ≤ hc(T′) ≤ K, |T′| = n − 1}` by out-neighbourhoods with at most `2K` in/out switches along one chosen Hamiltonian cycle of the parent, filtered by `hc ≤ K`. The program computes, per parent, one Hamiltonian cycle and the full path-count matrix `P` (one subset DP per start vertex), then evaluates `Σ P(w → u)` for every pattern (`C(m, 2) + C(m, 4) + …` switch sets, two phases each), emits the survivors, and the driver canonicalises them with nauty's `labelg` (digraph6) and deduplicates. Seeds are all tournaments of order 5 with `hc ≤ K` (from `gentourng`). A class with `hc = 1` can only arise from a parent with `hc = 1` and a one-interval `O`, which is opus-S74's generator; the program does not assume this and reproduces Garey's counts.

| order `n` | `hc = 1` | `F(2n − 6)` | `hc = 2` | `(n − 2)F(2n − 9) − L(2n − 10)` | `hc = 3` | `hc = 4` | time (K = 2) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 3 | 3 | 2 | (1) | 1 | 0 | seeds |
| 6 | 8 | 8 | 5 | 5 | 11 | 6 | |
| 7 | 21 | 21 | 18 | 18 | 38 | 25 | |
| 8 | 55 | 55 | 60 | 60 | 124 | 91 | |
| 9 | 144 | 144 | 191 | 191 | 390 | 315 | |
| 10 | 377 | 377 | 589 | 589 | 1193 | 1049 | |
| 11 | 987 | 987 | 1775 | 1775 | 3574 | 3390 | |
| 12 | 2584 | 2584 | 5257 | 5257 | 10537 | 10705 | 1 s |
| 13 | 6765 | 6765 | 15360 | 15360 | 30676 | – | 4 s |
| 14 | 17711 | 17711 | 44394 | 44394 | – | – | 20 s |
| 15 | 46368 | 46368 | 127171 | 127171 | – | – | 103 s |
| 16 | 121393 | 121393 | 361595 | 361595 | – | – | 574 s |

Order 5 is the known exception of the law (the regular tournament `RT₅`). The `hc = 3` and `hc = 4` columns agree with the exhaustive table of HYP-9459 through order 10 (11, 38, 124, 390, 1193 and 6, 25, 91, 315, 1049). Logs: [`data/tcgen_6_12.log`](data/tcgen_6_12.log), [`data/tcgen_13_16.log`](data/tcgen_13_16.log), [`data/tcgenK3_12_13.log`](data/tcgenK3_12_13.log), [`data/tcgenK4_6_12.log`](data/tcgenK4_6_12.log). Class lists (`"hc digraph6"`, canonical): [`data/cls_5.txt`](data/cls_5.txt) … [`data/cls_13.txt`](data/cls_13.txt) (`hc ≤ 2`), `clsK3_5..13` (`hc ≤ 3`), `clsK4_5..12` (`hc ≤ 4`); orders 14–16 (3–23 MB) are not committed and regenerate in 12 minutes. Every two-cycle class of orders 6–16 has a vertex of out-degree 1 and one of in-degree 1 (the driver checks all classes; see the logs).

**HYP-9459 (2), (3) at the new orders.** From the class lists, the numbers `u_k(n)` of `k`-cycle tournaments with no vertex of out-degree 1 ([`scripts/uk_counts.py`](scripts/uk_counts.py), [receipt](receipts/uk_counts.normal.out)) are `u_3(11), u_3(12), u_3(13) = 233, 610, 1597 = F(2n − 9)` and `u_4(11), u_4(12) = 178, 466 = 2F(2n − 11)`, exactly as HYP-9459 predicted from orders ≤ 10 (`u_1 = u_2 = 0` throughout, i.e. HYP-9459 (1) through order 13 for `k ≤ 2` from these lists and through 16 from the driver's check).

## 2. Structure of the two-cycle tournaments (FINITE-EXACT)

From [`scripts/moonstat2.c`](scripts/moonstat2.c) (per class: sole exits, sole entries, cut arcs, non-critical vertices, and for each non-critical `v` the cycle count of `Y − v`, `|O|`, the number of intervals of `O`, the parent's sole exits inside `O`) tabulated by [`scripts/moon_tab.py`](scripts/moon_tab.py) ([receipt](receipts/moon_tab.normal.out), orders 6–12; orders 13–14 in [`receipts/moon_tab_13_14.single_run.out`](receipts/moon_tab_13_14.single_run.out)):

| order | classes | (sole exits, sole entries) = (1,1), (1,2), (2,1), (2,2) | unique-cycle deletions 1 / 2 | `#cut + |NC| − n` |
|---:|---:|---|---|---|
| 10 | 589 | 154, 160, 160, 115 | 266 / 323 | 0 for all |
| 11 | 1775 | 437, 474, 474, 390 | 865 / 910 | 0 for all |
| 12 | 5257 | 1233, 1385, 1385, 1254 | 2706 / 2551 | 0 for all |
| 13 | 15360 | 3461, 4003, 4003, 3893 | 8240 / 7120 | 0 for all |
| 14 | 44394 | 9671, 11467, 11467, 11789 | 24598 / 19796 | 0 for all |

So (orders 6–14, exhaustive): a two-cycle tournament has exactly one or two sole exits and one or two sole entries; every non-critical vertex `v` has `hc(Y − v) ∈ {1, 2}` (Lemma 1), and exactly one or two of them have `hc(Y − v) = 1`; `|NC| ∈ [3, n − 2]`; and `#cut arcs + |NC| = n`. In particular the deletion grammar can be restricted to uniquely Hamiltonian parents: `D(n)` is exactly the set of two-cycle extensions of `UH(n − 1)` by at most two intervals (orders ≤ 14). The per-parent counts are far from uniform ([`receipts/perparent_7_11.single_run.out`](receipts/perparent_7_11.single_run.out): a parent of `UH(11)` has between 2 and 16 distinct two-cycle children, and there are 7,808 (parent, child) incidences for 5,257 children), so the natural bijection is not "delete a canonical vertex" but the gluing of §3.

For comparison, THM-4868 gives `#cut + |NC| = n + 1` for `hc = 1`; for `hc = 3` the excess over `n` is 0 or 1 ([`scripts/identity_tab.py`](scripts/identity_tab.py), [receipt](receipts/identity_tab.normal.out), orders 6–11; 12–13 in [`receipts/identity_tab_12_13.single_run.out`](receipts/identity_tab_12_13.single_run.out)): the numbers with excess 0 are 9, 27, 82, 246, 727, 2121, 6123, 17523 and with excess 1 are 2, 11, 42, 144, 466, 1453, 4414, 13153 (orders 6–13; neither in OEIS as far as checked).

## 3. The segment swap and the gluing theorem

[`scripts/two_cycle_structure.py`](scripts/two_cycle_structure.py) ([receipt](receipts/two_cycle_structure_6_12.normal.out), all 7,895 two-cycle classes of orders 6–12) enumerates both Hamiltonian cycles `C₁, C₂` of each class and finds, without exception:

1. **Segment swap.** `C₁` and `C₂` differ in exactly three arcs; the three runs of `C₁` between the changed arcs have sizes `(j, k, l)` with `min = j ∈ {1, 2}`: `C₁ = P J S R`, `C₂ = P S J R` where `J` (the *junction*) has one or two vertices. Counts of sorted size triples at order 12: `(1,2,9)` 1220, `(1,3,8)` 932, `(1,4,7)` 890, `(1,5,6)` 884, `(2,2,8)` 466, `(2,3,7)` 356, `(2,4,6)` 340, `(2,5,5)` 169.
2. **Vertex junction** `J = {s}`, `C₁ = P s S R`, `C₂ = P S s R`, `S = a₁ … a_k`, `R ∪ P = b₁ … b_l` (in cycle order from `R`): `A = Y[S + s]` and `B = Y[P ∪ R + s]` are **uniquely Hamiltonian**, `s` is a **sole entry of `A`** and a **sole exit of `B`**, and the **only arc from `S` to `P ∪ R` is `a_k → b₁`** (every other cross arc goes from `P ∪ R` to `S`). The rooted pair `((A, s), (B, s))` determines the class, and **every** pair occurs: for each `(k, l)` with `k + l = n − 1` the number of classes equals the number of distinct rooted pairs equals `F(2k − 3) F(2l − 3)` (order 12: 610, 466, 445, 442, 442, 445, 466, 610 for `k = 2 … 9`).
3. **Arc junction** `J = (s → s′)`: `A = Y[S + s + s′]` and `B = Y[P ∪ R + s + s′]` are uniquely Hamiltonian, `s, s′` are the two (consecutive) sole entries of `A` and the two sole exits of `B`, the only arc from `S` to `P ∪ R` is `a_k → b₁`, and the classes with `(|S|, |P ∪ R|) = (k, l)`, `k + l = n − 2`, number `F(2k − 3) F(2l − 3)` = the number of distinct rooted pairs (order 12: 233, 178, 170, 169, 170, 178, 233 for `k = 2 … 8`).
4. `NC(Y)` is a cyclic interval of exactly one of `C₁, C₂`; every cut arc lies on both cycles; `#cut + |NC| = n`.

**Theorem 2 (gluing; PROVED).** With `A`, `B`, `s` (or `s, s′`) as in §0, `Y = A ⊕_s B` (resp. `A ⊕_{s → s′} B`) has exactly two Hamiltonian cycles.

*Proof (vertex junction).* Write `X = {s, a₁, …, a_k}`, `Z = {b₁, …, b_l}`; `k, l ≥ 2` since `A, B` have order `≥ 3`. In `Y` the arcs from `X` to `Z` are `s → b₁` (`s` is a sole exit of `B`) and `a_k → b₁` (every other `a_i → b_j` is reversed by construction), so a Hamiltonian cycle `C` of `Y` enters `Z` exactly once, at `b₁`, then follows a Hamiltonian path `Q` of `Y[Z] = B − s` from `b₁` to some `b_t`, then returns to `X`. Since `s → b₁` is the only out-arc of `s` in `B`, every `b_t` with `t ≠ 1` beats `s`, so `Q + (b_t → s → b₁)` is a Hamiltonian cycle of `B`; by uniqueness `Q = b₁ b₂ … b_l`. The rest of `C` is a Hamiltonian path `R` of `Y[X] = A` from a vertex `x₀` with `b_l → x₀` to `x₁ ∈ {s, a_k}`.
If `x₁ = s`: the last arc of `R` enters `s`, hence is `a_k → s` (sole entry), and `R − s` is a Hamiltonian path of `A − s` from some `a_i` (`i < k`) to `a_k`; as `s → a_i` (the only in-neighbour of `s` is `a_k`), `(R − s) + (a_k → s → a_i)` is a Hamiltonian cycle of `A`, so `R − s = a₁ … a_k`, `x₀ = a₁`, and `C = a₁ … a_k s b₁ … b_l` (closing with `b_l → a₁`, a default cross arc).
If `x₁ = a_k`: `s` lies on `R` and its predecessor on `R` could only be `a_k`, the last vertex, so `s = x₀` and `R` is a Hamiltonian path of `A` from `s` to `a_k`, which closes with `a_k → s` to a Hamiltonian cycle of `A`: `R = s a₁ … a_k` and `C = s a₁ … a_k b₁ … b_l` (closing with `b_l → s`, an arc of `B`).
Both candidates are Hamiltonian cycles of `Y`, so `hc(Y) = 2`. For the arc junction replace `s` by the path `s → s′` (`s′ → b₁` and `a_k → b₁` are the only arcs from `X ∪ {s′}` to `Z`; `s′` is a sole entry of `A` with in-neighbour `s`, and `s` a sole entry with in-neighbour `a_k`; the same two cases give `a₁ … a_k s s′ b₁ … b_l` and `s s′ a₁ … a_k b₁ … b_l`). ∎

**Proposition 3 (injectivity; PROVED).** The pair of cycles of `Y` determines the three changed arcs, hence the junction `J`, the runs `S` and `P ∪ R`, and the pieces `A = Y[S ∪ J]`, `B = Y[P ∪ R ∪ J]` with their roots; an isomorphism of gluings maps cycles to cycles, hence junctions to junctions and, since it preserves in- and out-degrees, the entry piece to the entry piece. So distinct rooted pairs give distinct classes, and vertex and arc gluings are never isomorphic (their junctions have different sizes). ∎

**Corollary 4 (PROVED, given THM-4868 (F)).** `D_n ≥ E(n − 1) + E(n − 2) = (n − 2) F(2n − 9) − L(2n − 10)` for all `n ≥ 6`, because the vertex gluings with `|S| = k`, `|P ∪ R| = l`, `k + l = n − 1`, number (rooted `(A, s)` with `A ∈ UH(k + 1)`, `s` a sole entry) × (rooted `(B, s)`, `B ∈ UH(l + 1)`, `s` a sole exit) `= F(2(k + 1) − 5) · F(2(l + 1) − 5) = F(2k − 3) F(2l − 3)`, and the arc gluings with `k + l = n − 2` number (`A ∈ UH(k + 2)` with two sole entries) × (`B ∈ UH(l + 2)` with two sole exits) `= F(2(k + 2) − 7) · F(2(l + 2) − 7) = F(2k − 3) F(2l − 3)` (a class with two sole entries has exactly one ordered pair of them, THM-4868 (C)). ∎

**Consequence.** `D_n = (n − 2) F(2n − 9) − L(2n − 10)` holds at order `n` **if and only if** every two-cycle tournament of order `n` is a gluing. The census (§1) therefore proves the classification for every `n ≤ 16`, while §3's direct verification covers `n ≤ 12`. The classification would also give HYP-9459 (1) (a gluing has a sole exit, namely a sole exit of `A` other than `a_k`; that such a vertex always exists is a statement about `UH` tournaments rooted at a sole entry, exhaustive through order 16 here but not derived from THM-4868).

**Why the generating function has a squared denominator.** `Σ_{k ≥ 2} F(2k − 3) xᵏ = x²(1 − x)/(1 − 3x + x²)` is the generating function of the rooted classes, and `D(x) = (x + x²) · (that)²`: the two-cycle tournaments are, up to the finitely many exceptions that do not exist, pairs of rooted uniquely Hamiltonian tournaments, glued at a vertex or along an arc. The equivalent form `D_n = a_{n−1} + (n − 4) b_{n−1}` with `a_m = F(2m − 8)`, `b_m = F(2m − 7)` (the one-/two-sole-exit classes of `UH(m)`) is verified in the moon_tab receipt but is not the shape of the bijection (§2).

## 4. What remains

1. **The converse of Theorem 2** (every tournament with exactly two Hamiltonian cycles of order `≥ 6` is a vertex or arc gluing). Equivalent to the count law for all `n`; exhaustive through order 16. A proof would presumably run by Lemma 1 on a uniquely Hamiltonian parent with THM-4868's local roundness, showing that a two-interval or one-interval extension with `hc = 2` is a gluing.
2. The identity `#cut arcs + |NC| = n` for `hc = 2` (and `n + 1 − ε`, `ε ∈ {0, 1}`, for `hc = 3`): for gluings it should follow from THM-4868's `n + 1` for the two pieces.
3. The `hc = 3` and `hc = 4` censuses can be pushed (order 14 for `hc ≤ 3` takes minutes; `hc ≤ 2` at order 17 about ten hours); the two `hc = 3` sequences by excess are new.

## 5. Reproduction

```
cc -O2 -o scripts/tcgen scripts/tcgen.c && cc -O2 -o scripts/tcgen2 scripts/tcgen2.c && cc -O2 -o scripts/moonstat2 scripts/moonstat2.c
cd scripts && python3 seed5.py && python3 tcgen_driver.py 6 16            # hc <= 2 (cls_n.txt; order 16 ~ 10 min, 23 MB)
python3 seedK.py 3 && python3 tcgen_driver2.py 6 13 3                      # hc <= 3 (clsK3_n.txt)
python3 seedK.py 4 && python3 tcgen_driver2.py 6 12 4                      # hc <= 4 (clsK4_n.txt)
grep '^2 ' cls_12.txt | ./moonstat2 > ../data/moon2_12.txt                  # deletion statistics (any order)
cd .. && python3 scripts/moon_tab.py 6 12 && python3 scripts/identity_tab.py 6 11
python3 scripts/two_cycle_structure.py 6 12 && python3 scripts/gf_identity.py && python3 scripts/perparent.py 7 11
shasum -a 256 -c SHA256SUMS
```

Receipts for the Python scripts are byte-identical under `python3` and `python3 -O`; the generator logs and the order 13–14 tabulations are single-run records. Tools: nauty `gentourng` (seeds) and `labelg` (canonical digraph6, `-q`).
