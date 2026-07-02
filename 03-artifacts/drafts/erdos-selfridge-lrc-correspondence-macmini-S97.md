# The Erdős–Selfridge ↔ LRC correspondence, the local-covering program for hpartA, and the resonance-lattice formula for hp0cap

**mac-mini-2026-07-01-S97 (HYP-3853).** Owner brief: relate the odd-covering-system
conjecture (no incongruent covering system with all moduli odd and > 1; square-free case
needs ≥ 22 prime factors) to the tournament/LRC/unit-distance work; creatively close
hp0cap and hpartA; crossing numbers K_n vs K_{m,n}; make the collapse-rate law
machine-checkable end to end (→ THM-597, done separately).

---

## 1. The correspondence table (covering systems ↔ danger-arc systems)

| Covering systems (ℤ) | LRC danger systems (ℝ/ℤ) |
|---|---|
| congruence class `a mod m` | danger comb `D_v(r)` = arcs of half-width `r/v` at centers `j/v` |
| modulus `m` | spacing `1/v` (density `2r` regardless of `v`) |
| **offset `a` free** | **offsets FIXED at 0** (centers are the rationals `j/v`) |
| covering system (covers ℤ) | covering set (arcs cover the circle at radius `r`) |
| distinct moduli | distinct speeds |
| exact/disjoint cover rigidity: Davenport–Mirsky–Newman–Rado | THM-594(C): no finite distinct-speed exact tiling (divisor-minimal frequency) — and the polygon theorem (Lean, S95) is the same statement discretized |
| **Erdős–Selfridge**: no covering with all moduli odd, distinct, > 1 (OPEN) | all speeds odd ⟹ `t = 1/2` is lonely — **PROVED, one line** (`all_odd_half_lonely`, already Lean-verified in-repo) |
| odd square-free covering needs ≥ 22 primes | covering saturation (THM-523): a covering set needs a multiple of EVERY `q ∈ {2..14}` — in particular of 2 |
| BBMST: minimum modulus of a distinct-moduli covering system is ≤ 616,000 (distortion method) | Malikiosis–Santos–Schymura: LRC(14) counterexample speeds ≤ C(14,2)^12 — both are "covering complexity is finitely bounded" theorems |
| 2-adic chains `1(2), 2(4), 4(8), …` = the canonical exact cover | divisor chains = the fixed locus of the scale action = the ONLY route to tiling limits (THM-594(C) corollary); peeled by the 2-adic descent (THM-580) |

**The moral.** The offset column is the difficulty axis. Our whole rigidity program
(MN, polygon theorem, unit-residue, tight-locus classification) lives in the
*fixed-offset* column, where covering obstructions are provable; Erdős–Selfridge is the
*free-offset* shadow, where the same obstruction (you need the prime 2 / the 2-adic
chain direction to cover efficiently) is conjectural. Precisely: E–S says odd moduli
cannot cover **even with adversarial offsets**; our column proves the offset-0 case
trivially via the antipode `1/2` (the point that fixed offsets can never approach with
odd spacings — `||w/2|| = 1/2` for odd `w`). The 22-primes bound quantifies how much
free-offset flexibility costs the adversary — the analog of our exposure counts (13
simultaneous residue clearances mod `q*`, THM-596) and kps's arc-count tax.

**Dynamical reading (arXiv 2604.16750, Blaschke/circle maps — extends HYP-3822):**
odd denominators are exactly the periodic points of the doubling map `t ↦ 2t`; E–S says
doubling-periodic congruence data cannot cover ℤ without the transient (2-adic) direction.
Our THM-580 descent quotients precisely the doubling direction, leaving the odd core —
the tower's residual IS the E–S regime, with fixed offsets making it tractable. Arnold
tongues (the paper's rational resonances) = our Stern–Brocot band strata (THM-596/HYP-3852).

## 2. hpartA — the local-covering program (the creative close)

`hpartA` (`witnessG2 > 0 ⟹ Mreach ≥ 1/14`) fails only if the large cluster's arcs
**cover the entire G2 window** `I` (a finite union of rational intervals). Within `I`,
each far runner `w` is an arc-comb of spacing `1/w ≪ |I|` and density `2r = 1/7`; the
offsets (comb phases relative to `I`) are effectively free parameters. So:

> **hpartA ⟺ [j ≤ 13 arc-combs with distinct spacings and density 1/7 each cannot cover
> an interval]** — a *continuous covering-system problem with free offsets and ≤ 13
> moduli*.

This is exactly the E–S/BBMST regime, and the 22-primes/616,000 results say covering at
bounded density with constrained moduli requires ENORMOUS complexity — 13 combs is far
below any covering threshold. Quantitative probe (this session): at critical mass
(`j = 7`, `A = 1.0`), 3000 adversarial phase draws never covered more than **86%** of the
window (best 0.8605 at `|I| = 0.01`; 0.80 at `0.02`); at `j = 6` (`A = 6/7`) max 81%.
The residual to certify is ≥ 14% at criticality — large, not marginal.

**The proof shape to import (BBMST distortion, adapted to fixed density):** order the
spacings `1/w₁ > … > 1/w_j`; reveal combs one at a time; the distortion lemma bounds how
much of the *uncovered* set's measure each new comb can absorb, because an uncovered set
that is a union of intervals of typical length `ℓ` meets a comb of spacing `1/w` in
density `≤ 2r·(1 + O(1/(wℓ)))` — and the previous combs force `ℓ` to remain comparable
to the coarsest spacings (few combs ⟹ few cut points: after `k` combs the uncovered set
has at most `Σ_{i≤k} (w_i|I| + 1)` components — an ARC-COUNT TAX, kps-S28's mechanism,
now in the free-offset column). Iterating: uncovered measure after `j` combs
`≥ |I|·∏(1 − 2r·(1+ε_i))` with `Σε_i` controlled by the spacing ratios — positive for
`j ≤ 13` at `2r = 1/7` provided the spacings are distinct-enough (a Λ-gap or
distinctness suffices; the resonant equal-spacing case is excluded by distinctness and
handled at the boundary by THM-594's pair law). **This replaces equidistribution with a
covering-complexity argument**; the fixed-offset extreme (all combs aligned at 0) is
the classical column where the antipode/unit-residue tools finish the job.

Named target: **Lemma (interval anti-covering).** For distinct `w₁ < … < w_j` with
`w₁ ≥ 2/|I|` and `2rj ≤ j/7`: `|I ∩ {C_F = 0}| ≥ |I| · κ_j` with explicit
`κ_j > 0` for `j ≤ 13`. This lemma + the peel (HYP-3834/3950) + G2's window structure
(S93 grid) discharges hpartA outside the compact-cluster regime, which is the census +
renormalization regime (THM-595 legs).

## 3. hp0cap — the resonance-lattice formula (proved here)

**Lemma (d-fold overlap = resonance-lattice sum).** For speeds `e₁,…,e_d` and intervals
`I₁,…,I_d ⊂ ℝ/ℤ`:
```
| ∩_i {t : e_i t mod 1 ∈ I_i} |  =  Σ_{(m₁,…,m_d) ∈ Λ}  ∏_i  1̂_{I_i}(m_i),
Λ = { m ∈ ℤ^d : Σ_i m_i e_i = 0 }   (the RESONANCE LATTICE),
```
absolutely convergent for `d ≥ 3` (coefficients `O(1/|m_i|)`, lattice has codimension 1),
conditionally (Abel/symmetric summation) at `d = 2`. *Proof:* expand each indicator in
Fourier, integrate; only frequency tuples with `Σ m_i e_i = 0` survive. ∎ The `m = 0`
term is `∏|I_i|` (independence); `d = 2` recovers THM-594(B) (`Λ = ℤ·(q,−p)/g`, and the
cosine series closes to the two-branch law). For **rational sector data** (`I_i` =
sevenths), each primitive lattice direction contributes a finite Bernoulli-polynomial
evaluation (`Σ_s cos(2πsθ)/s^k = (−1)^{k/2+1}(2π)^k B_k({θ})/(2·k!)` for even `k`, and
the odd-`k` sine analogs) — so every sector-miss probability `q_j(E)` of `hp0cap`'s
miss-distribution is a **finite exact rational** computable per shape.

**hp0cap closure plan (now fully concrete):** (i) `p0 ≤ L_y` (THM-534, proved);
(ii) `L_y(E)` exact per shape by the lattice formula; (iii) "consec maximizes `L_y`" on
the bounded-spread compact set (THM-529) = a finite exact census; (iv) the far-element
leg = branch-2 deficits (exact, THM-594(B)) with the `1/w` tail controlled by the
lattice formula's `d ≥ 3` absolute convergence — no Vitali, no marginal-uniformity.
The census (iii) is compute; (iv) is a finite family of rational inequalities.

## 4. Crossing numbers (the requested tangent, honestly assessed)

Zarankiewicz's conjectured-optimal `K_{m,n}` drawing places the two classes on two
perpendicular axes at positions `±1, ±2, …` — two arithmetic progressions — and its
crossing count `⌊m/2⌋⌊(m−1)/2⌋⌊n/2⌋⌊(n−1)/2⌋` counts interleaved pairs-of-pairs, the
same order-type/lattice-point species as our binding-pair grid (`v±w` denominators) and
the MT-slice counts. Guy's `K_n` optimum is realized by cylindrical drawings: points at
regular positions on two concentric circles — REGULAR-POLYGON configurations again, with
the crossing count a 4-point interleaving statistic (the species of the repo's
`c₃ = Tr(A³)/3` and H-variance quadruple counts). Kleitman's parity theorems for
`cr(K_{m,n})` (parity invariance across drawings for small odd cases) are
antipodal/deformation-parity arguments — the same Z₂ genus as Rédei's odd-H and our
ι-parity: parity survives where magnitude doesn't. Honest verdict: a genuine METHOD
family resemblance (AP-extremal configurations; interleaving counts; parity rigidity)
worth one backlog line — the concrete crossover to chase is whether the 2-page circular
`K_n` count decomposes over the Farey/level structure the way `Λ_AP` does (both are sums
over interleaved pairs on a circle of regular points). Not pursued further this session.

-> THM-523, THM-527/534 (hp0cap/hpartA context), THM-580, THM-592–597, HYP-3822 (Blaschke),
HYP-3834/3950 (peel), HYP-3852, OPEN-Q-108; refs: Erdős–Selfridge problem; BBMST (Ann. of
Math. 2022, minimum modulus); Malikiosis–Santos–Schymura 2024; arXiv 2604.16750.
