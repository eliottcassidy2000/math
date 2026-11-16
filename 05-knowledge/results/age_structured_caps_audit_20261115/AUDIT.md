# Independent audit of the sharp growth lemma for span-3 ribbons

Session klein-2026-11-15-age. Audited objects: [`SHARP_GROWTH.md`](../chanpak_target_ribbon_audit_20261116/SHARP_GROWTH.md)
(the uniform bound `H(T) ≤ H(T+x) ≤ 39·H(T)`) and [`SHARP_EIGHT.md`](../chanpak_eight_braid_20261116/SHARP_EIGHT.md)
(the eight per-letter sharp constants), both written on 2026-11-15/16 by other sessions and not previously audited.
The author of this audit proved the decomposition `H(Wℓ) = e + s + Z_ℓ` and the finite-exact data of THM-5243 in a
parallel session that did not know of these proofs; THM-5243 and HYP-9803 (b) called the lemma open. Verdict first.

**Verdict: SOUND.** Both proofs are correct as written; their checkers replay (receipts `replay_check_sharp_growth.out`,
`replay_check_sharp_eight.out`); and the exhaustive data of this session over all `8^n` words of length `n ≤ 9`
(134,217,728 words at length 9) gives exactly the eight constants as the per-letter maxima of `H(Wℓ)/H(W)`:

| letter `ℓ` (bits = beaten offsets) | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|---|
| sharp constant `C_ℓ` (SHARP_EIGHT, proved) | 1 | 9 | 17 | 19 | 33 | 33 | 39 | 31 |
| max `H(Wℓ)/H(W)`, all words of length ≤ 9 (this session) | 1 | 9 | 17 | 19 | 33 | 33 | 39 | 31 |
| max `Z_ℓ/H`, all words of length ≤ 9 (this session) | 0 | 8 | 16 | 18 | 32 | 32 | 38 | 30 |

(receipts `hratio5_9_b*`, features 130–137 and 122–129 at history length 0). The third row shows that at every per-letter
maximiser the bridged two-covers carry all but one of the paths (`e + s = 1`). Consequences recorded below.

## 1. The argument, step by step

Notation of SHARP_GROWTH: old sites `0, …, n−1` in insertion order, new vertex `x = n`, `A = N⁺(x)` among the last three
old sites. `H` is the Hamiltonian-path count, `H(∅) = 1`.

1. **Product inequality.** For every `S ⊆ V(T)`, `H(T−S)·H(T[S]) ≤ H(T)`. Proof used: the odd-cycle packing identity
   (THM-002 / Grinberg–Stanley), `H = Σ_packings 2^(#cycles)`; a packing of `T−S` and a packing of `T[S]` have disjoint
   supports, their union is a packing of `T` of weight the product, and distinct pairs give distinct unions. Checked.
   The same argument gives monotonicity `H(T−S) ≤ H(T)` and, for pairwise disjoint `S_1, …, S_m`,
   `H(T − ∪S_i)·Π H(T[S_i]) ≤ H(T)`, which CAPS.md uses.
2. **Paths ending at `x`.** Deleting `x` is an injection into the Hamiltonian paths of `T`: at most `H(T)`. Checked.
3. **Tail deletion.** A path not ending at `x` is `P₁, x, P₂` with `P₂` starting at some `b ∈ A` (the arc `x → b`). Let `S`
   be the maximal terminal block of `P₁` all of whose vertices are beaten by `b` (`S = ∅` allowed; `S = V(P₁)` allowed).
   If `P₁ − S` is nonempty, its last vertex is not beaten by `b`, hence beats `b` (tournament), so `(P₁ − S), P₂` is a
   Hamiltonian path of `T − S`; if `P₁ − S` is empty, `P₂` is one. Checked.
4. **Recovery.** Fix `b` and the set `S`. From the Hamiltonian path `Q` of `T − S` and the ordered block (a Hamiltonian
   path of `T[S]`) the original path is recovered by inserting the block and then `x` immediately before `b` in `Q`.
   So the paths with these data number at most `H(T−S)·H(T[S]) ≤ H(T)`. Checked. (Maximality of the block and the
   condition "the vertex before the block beats `b`" only shrink the image; they are not needed for the bound.)
5. **Counting the data.** `S ⊆ N⁺(b)`, and a nonempty `S` ends at a vertex that precedes `x` in the path, hence beats `x`,
   hence lies outside `A`. For `b = n−j` (`j ∈ {1,2,3}`), `b` beats at most `3` older vertices (span) and the `j−1` newer
   window vertices; writing `d` for its out-degree inside the window and noting that at most `j` of its older
   out-neighbours lie left of the window, `|N⁺(b)| ≤ j + d`, so the admissible `S` number at most `1 + 2^(j+d) − 2^c`
   with `c = |N⁺(b) ∩ A|`. Checked. (A cleaner count used in CAPS.md: a window vertex of age `a` has `|N⁺| ≤ 3 + a`.)
6. **The table and the single overcount.** Summing `1 + Σ_b (1 + 2^(j+d) − 2^c)` and maximising over the eight
   orientations of the three window arcs gives 1, 9, 17, 20, 33, 36, 40, 39 for the masks 0–7. The 40 occurs only for
   `A = {n−2, n−3}` with the window transitive and forward; there the subset `S = {n−2, n−6}` for `b = n−3` is counted
   but impossible (its only Hamiltonian order is `n−6 → n−2`, which ends inside `A`). Removing it gives 39. We replayed
   `check_sharp_growth.py` ("64 local cases, unique excluded 40 case, sharp 39 witness … passed"). Checked.
7. **Sharpness.** The fixture (old transitive order `3, 0, 4, 1, 2, 5`, `x = 6` beating `{3, 4}`) has `H(T) = 1`,
   `H(T+x) = 39`. It is the block `B` of THM-5243 (2): `TT_6` plus a vertex beating the vertices of ranks 1 and 3, whose
   count `5·2^3 − 1 = 39` THM-5243 (2′) proves in closed form. The two sessions' maximisers agree.

**SHARP_EIGHT.** The weighted refinement replaces the subset count by `Σ_S e_A(S)/H(T[S])` with `e_A(S)` the number of
Hamiltonian paths of `T[S]` ending outside `A`, which follows from step 4 by multiplying and dividing by `H(T[S])`. Only
the last six old sites enter, since `b` is within three of `x` and beats older sites only within three of itself; the
padding of earlier sites by universal sources preserves `A`, every weight and both counts because every site left of
the six beats all six. The replay enumerates the `2^12` orientations of the twelve possible back arcs among six sites
("32768 exact fractional bounds and eight independent equality witnesses: [1, 9, 17, 19, 33, 33, 39, 31]"). Checked.

## 2. What the exhaustive data add

- **Exact agreement to length 9, letter by letter.** The per-letter maxima over all words of length ≤ 9 equal the proved
  constants (table above), so each constant is attained already at length 5 and never exceeded through length 9. This
  is a second, independent route to the eight constants (the proof is local on six sites; the data are global).
- **The two-fragment part.** `max Z_ℓ/H = C_ℓ − 1` for every letter (0, 8, 16, 18, 32, 32, 38, 30). The lemma of
  HYP-9803 (b), `Z_ℓ ≤ 38·H`, is therefore PROVED by SHARP_GROWTH (paths not ending at `x` are `s_ℓ + Z_ℓ ≤ 38·H`)
  and sharp exactly where `e_ℓ + s_ℓ = 1`.
- **Maximisers.** THM-5243 (2) lists every maximising pair to length 9: the sealed gadgets `T_(W′) ⇒ B`. The proof's
  fixture is the `W′ = ∅` case, and sealing (three empty letters) shows why the bound is attained at every count
  `H(W′)`: the deletion data for the gadget do not see `W′` at all.

## 3. The record, corrected

- HYP-9802 (ii) already recorded the lemma as PROVED via SHARP_GROWTH; HYP-9803 (b) and THM-5243 (status, (2), (3))
  still called it OPEN and proposed an "age-structured potential on the 119 states" as the next step. Both entries are
  updated in this commit: HYP-9803 (b) is PROVED (by SHARP_GROWTH, audited here), the proposed potential is withdrawn
  as a proof device, and THM-5243's finite-exact statements stand as corroboration and sharpness data.
- The obstruction statements of THM-5243 (3) remain correct and now have a conceptual explanation (CAPS.md §4): the
  deletion proof uses the product inequality `H(T−S)H(T[S]) ≤ H(T)`, a statement about sub-tournaments that the
  119-dimensional transfer state does not carry; no linear functional of the state can express it, which is why the
  annihilator LP and the coordinate-cap cone had to fail while a combinatorial injection succeeds.

## 4. What stays open

The upper growth bound is settled. Open are the lower-growth certificate on rank-constant literal stretches
(HYP-9742 for literal span 3, HYP-9802 (iii)), two-apex completeness (HYP-9803 (a)), and the sharp per-profile caps
of CAPS.md where the deletion bound is not attained.
