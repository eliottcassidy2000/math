# LRC difficulty by n is non-monotone — it is the apex prime's arithmetic; and the tournament reconstruction wall at 7 vertices is the combinatorial shadow of LRC(14)

*kind-pasteur-2026-07-01-S17. Applying this session's tournament-reconstruction work to the LRC, and surveying the LRC across n above and below 14. The result: LRC difficulty is NOT monotone in n — it is governed by the arithmetic of the apex prime `p` (with `n=2p`), and the places that are surprisingly easier/harder/different are exactly the places where `p`'s three pillars (Mersenne / Heegner / 3-mod-4) turn on or off. Complements mac-mini-S87's six merged-metagraph↔LRC bridges (which live at n=14) with the by-n variation.*

## The apex-prime table (the backbone)

The LRC(`n=2p`) apex object is the Paley/heptagon tournament on `p` vertices, and the repo's three proof tools are the three arithmetic properties of `p`:
- **Mersenne** (`p=2^k−1`) → the 2-adic descent `2p=2·p` peels one `Z₂` (THM-580).
- **Heegner** (`Q(√−p)` class number 1) → the Fejér–Bochner SOS minorant on the totally-real cyclotomic (S75e/HYP-3535), the floor.
- **3 mod 4** → the Paley *tournament* exists and Borsuk–Ulam odd-degree applies (THM-581).

| `p` | `n=2p` | pillars | atoms `φ(2p)` | note |
|---|---|---|---|---|
| **3** | **6** | M·H·3 | 2 | hard regime (`p≡3 mod4`) but small → **PROVED**; trivial tournament space (2 classes) |
| 5 | 10 | — | 4 | `p≡1 mod 4`: reflection = **automorphism → Brouwer fixed point → EASIER** (SOS) |
| **7** | **14** | M·H·3 | 6 | all three → **first hard case**; tournament space on 7 vtx = the **reconstruction wall** |
| 11 | 22 | ·H·3 | 10 | `p≡3`: loses **Mersenne** (no 2-adic descent) |
| 13 | 26 | — | 12 | `p≡1 mod 4`: Brouwer → **easier** |
| 19 | 38 | ·H·3 | 18 | `p≡3`: loses Mersenne — but n=19 itself is a **cleaner control** than 14 |
| 23 | 46 | ··3 | 22 | `p≡3`: loses Mersenne+Heegner, but **tight-equiangular d=23** |
| 31 | 62 | M·3 | 30 | `p≡3`: loses **Heegner** (no `√−31` SOS) |

**The primary difficulty axis is `p mod 4`, not pillar count** (survey correction). For `p≡1 mod 4` the complement reflection is an *automorphism* — a **Brouwer fixed point** exists, directly witnessing the lonely runner (SOS suffices) — so these are **easier** (conjecturally uniform). For `p≡3 mod 4` the reflection is an *anti*-automorphism — a **free ℤ₂**, no fixed point, needing **Borsuk–Ulam odd degree** — the hard regime. Within the hard regime, the three pillars (Mersenne/Heegner/3-mod-4) are the *tools*, all co-present only at `p∈{3,7}`. `n=6` is proved (small); `n=14` is the **first hard case**.

## Easier / harder / clearer — and the surprises

- **Easiest / clearest: `n=6`** (apex 3). All three pillars *and* the tournament space on 3 vertices is trivial (2 classes, reconstruction unique). Everything aligns; it is proved.
- **Richest but first-open: `n=14`** (apex 7). All three pillars, maximal symmetry (Paley heptagon, `D₇`, tight-equiangular `d=7`) — yet the tournament space on 7 vertices is exactly where this session's **reconstruction wall** sits: 90% OCF-cospectral, `(I(Ω,x),d)` fails, local invariants stop determining the class. **The gentlest place to break the LRC open is also where the combinatorics first becomes irreducibly complex.** That is not a coincidence to me — the same apex-7 richness that supplies all three tools is what makes the class-space first-degenerate.

**Surprise 1 (the big one) — the cases that *lack* the number-theoretic pillars are EASIER, not harder.** `p ≡ 1 mod 4` (`n=10, 26, 34`) has none of the three pillars — yet it is the **easy** regime, because the complement reflection is an automorphism (Brouwer fixed point → SOS witnesses the lonely runner directly). The pillars are not "difficulty"; they are the **tools you only need when the reflection is a free ℤ₂** (`p≡3 mod 4`, Borsuk–Ulam). So "fewer pillars" cuts *both* ways: `p≡1 mod 4` (no pillars, easy) vs `p≡3 mod 4` losing a pillar (`22, 62`, hard *and* under-tooled). Difficulty is governed first by `p mod 4`.

**Surprise 2 — non-monotone within the hard regime; `n=19` can be a *cleaner* case than `n=14`.** Above 14 you *lose tools*: `n=22` (apex 11) loses Mersenne, `n=62` (apex 31) loses Heegner. But it is not simply "bigger = harder": `n=19` (with `C=37` prime, both `n` and `C` prime — no side-channel apex debt) is flagged as a **cleaner laboratory** than 14 or 21 for isolating which parts of an LRC(14) argument are structure-specific vs overfitted to the apex-7 residue. Difficulty is not sorted by runner count.

**Surprise 3 — tight-equiangular `{7, 23}`.** `n=14` and `n=46` are the two apices whose dimension is a *tight* equiangular-line dimension (absolute bound `C(d+1,2)` attained). `n=46` has only the 3-mod-4 pillar but carries the equiangular/tight-frame structure (my equioscillation↔equiangular reflection) — a **different special mechanism**, marking `n=46` as the next "geometrically special" case after 14, orthogonal to Mersenne/Heegner.

**Surprise 4 — the open core is FINITE (and covering sets are *too* lonely).** MSS 2024 (arXiv:2411.06903) bounds all LRC(14) counterexample speeds by `91^12 ≈ 3.2×10²³` — refuting the repo's long-held "unbounded far configs are the difficulty" premise: LRC(14) is a *finite* high-dimensional packing/satisfiability problem, and Rosenfeld-style prime-filtering already closed n=8,9,10. And the "loneliness paradox": *covering* 13-sets are over-lonely (`M≈0.095 ≫ 1/14`) — more structure makes them *harder* to be counterexamples, not easier.

## Where this session's tournament work plugs in

- **The Brouwer/Borsuk–Ulam split IS this session's sgn/parity finding (S11).** The whole `p mod 4` difficulty axis is the sign of the complement permutation. My S11 result: for a *tournament*, `|Aut|` is odd so `Aut ⊆ A_n` (automorphisms are even permutations) and the complement/converse `ι` is an **anti-automorphism = odd permutation** — a fixed-point-free (free) `ℤ₂`. That free `ℤ₂` is *exactly* why `p≡3 mod 4` (Paley tournament) needs **Borsuk–Ulam** (no fixed point). For `p≡1 mod 4` the complement is an *automorphism* (a Paley *graph*, `ι` even/fixed) → **Brouwer**. So the sgn character `even-Aut / odd-anti-Aut` **is** the Brouwer/Borsuk–Ulam dichotomy **is** the LRC easy/hard split. My "even degree vs even parity" separation was the combinatorial statement of the topological obstruction.
- **The reconstruction wall at 7 vertices is the shadow of LRC(14).** Below 7 (apex 3 → 3 vertices) the tournament space is trivial and LRC(6) is proved; at 7 the class-space first becomes locally-undecidable and LRC(14) is first-open. The tournament-space complexity *by number of vertices* tracks the LRC difficulty *by apex prime*.
- **The half-tiling is the LRC witness side.** This session's stance — the SC/blue spine lives on the half-tiling (σ = complement) — is exactly mac-mini-S6's two-index: half-tiling = Rédei-**odd** witness, LRC = ι-**even** existence. My blue-spine connectivity / odd-parity results are structural facts about that witness half.
- **Bridge-2 convergence (spectral separation is universal).** My twin separator (S15): the OCF/H-cospectral tournament twins `{4,6}` are separated only by the *spectrum* (`d=det(I+S)`, skew char-poly). mac-mini-S87's Bridge 2 (grounded): the LRC champions (construction/AP/GW) are **trace-cospectral** (all `1.857 = |S|·2r`, the 1st moment) but **Gram-spectrally distinct** (top eig `0.5006/0.4992/0.4803`). **The same move — the 2nd moment / spectrum separates 1st-moment-cospectral objects — resolves twins in *both* worlds.** The determinant/Gram spectral coordinate is the shared "second axis" (the determinant-lens `d ⊥ H`).

## Mental model

**LRC difficulty is a function of the apex prime's arithmetic, not of `n` — on two axes.** *Axis 1 (the coarse one): `p mod 4`.* `p≡1 mod 4` → complement is an automorphism → Brouwer fixed point → **easy** (SOS); `p≡3 mod 4` → complement is a free `ℤ₂` → **Borsuk–Ulam**, the hard regime. This axis *is* the sign of the complement permutation — my S11 even-Aut/odd-anti-Aut finding. *Axis 2 (within the hard regime): pillar count.* All three (Mersenne/Heegner/3-mod-4) only at `p∈{3,7}`; `n=6` proved, `n=14` first-hard, larger apices lose tools. And the tournament-space reconstruction wall at `p` vertices is the thermometer: trivial at `p=3`, first-irreducible at `p=7`. **Don't read difficulty off `n`; read it off `p mod 4`, then `p`'s pillars.**

## Honest status & next

- **Grounded:** the apex-pillar table (arithmetic); the tournament-space reconstruction facts (this session, p≤7); the Bridge-2 spectral convergence (my twins ↔ mac-mini's LRC Gram, both 1st-moment-cospectral / 2nd-moment-distinct).
- **Survey-corrected:** LRC(n) = n runners = n−1 speeds; classical proved to n≤7, recent frontier n=8,9,10 (Rosenfeld 2025 +), so **14 is the first case the repo's three-pillar machinery targets**, and MSS 2024 makes its open core *finite* (speeds ≤ `91^12`). The `p≡1 mod 4` cases are **easier** (Brouwer), a direction I initially had backwards.
- **Interpretive:** "tournament-space complexity by p tracks LRC difficulty by 2p" is a strong analogy grounded at p=3,7, not a theorem; the by-n difficulty here is about the *apex structure* (`p mod 4`, then pillars), not raw runner count.
- **Next:** extend Bridge 2 by computing the LRC champion Gram spectra at other `n=2p` (does the spectral gap track the pillar count?); test `n=46` (apex 23, tight-equiangular) for an equiangular-line certificate; verify the `p≡1 mod 4 ⇒ Brouwer/SOS` easy-regime conjecture (n=10) explicitly; and use `n=19` (the cleaner control) to check which parts of an LRC(14) argument are apex-7-specific.

— Related: `twentyeight-the-octonion-apex-and-the-three-pillars.md` (three pillars, {3,7} unique), `14-is-the-heptagon-dihedral-group-borsuk-ulam-not-brouwer-kps.md` (THE Brouwer-vs-Borsuk–Ulam / p-mod-4 certificate), HYP-3813 (mac-mini-S87, six merged-metagraph↔LRC bridges — the at-14 companion), HYP-3794 (MSS 2024 finiteness), `easier-frontier-targets-lrc-unit-distance-s678.md` (non-monotone, n=19 cleaner), `three-evens-two-poles-…` (equiangular/equioscillation, tight `{2,3,7,23}`), `even_degree_vs_even_parity_kps.py`/S11 (complement = odd anti-automorphism = the free ℤ₂), `the-SC-spine-is-the-half-tiling-…`, THM-580/581/582, the-determinant-lens (`d ⊥ H`). Script: `04-computation/lrc_by_n_apex_tournament_complexity_kps.py` (+ .out). Not a HYP reservation.
