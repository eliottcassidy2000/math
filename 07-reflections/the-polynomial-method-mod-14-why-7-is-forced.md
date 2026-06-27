# The polynomial method mod 14: why 7 is forced — units, null polynomials, and the c-lift threshold

*kind-pasteur-2026-06-27-S31ag, continuing the arXiv:2604.23906 synthesis. Having mapped LRC(14) onto
the paper's composite-`k+1` case, I dug into the exact algebraic mechanism by which the polynomial method
(Prop 4.1) fails at `k+1=14`, and found it gives a clean, quantitative reason for the apex prime 7 that
**matches both the project's "6 inner sectors" (φ(14)=6) AND mac-mini's gK8 "low-order moment-LP", from a
single source.** It also yields a one-parameter generalization of THM-573 (the c-lift sieve family) and a
sharp statement of why `c=7` is optimal.*

## 1. Three numbers, one wall

For the polynomial method to run over `ℤ_{k+1}` you need two things, and at `k+1=14` both fail by the
same prime 7:

- **All indices `1..k` invertible.** The construction `P(X)=∏_{i=1}^{k}(v_i i^{-1}+X)` needs `i^{-1}` mod
  `k+1`. Over `ℤ_14` only the **`φ(14)=6` units** `{1,3,5,9,11,13}` are invertible; the **7 non-units**
  `{2,4,6,7,8,10,12}` (six evens + the seven) are not. So the method can only form the indicator at
  `6` of the `13` positions. **Deficiency `13−6 = 7`.**
- **Interpolation: degree-`k` polys agreeing at `k+1` points are equal.** False over `ℤ_14` (not a
  field). The minimal monic **null polynomial** (vanishing as a function on all of `ℤ_14`) is
  `b_7(X)=∏_{j=0}^{6}(X−j)` — degree **7** — because `7! = 5040 = 14·360`, so a product of 7 consecutive
  integers is `≡0 mod 14`. Two degree-13 polynomials can agree as *functions* yet differ by `R(X)·b_7(X)`
  (`deg R ≤ 6`), so the **leading-coefficient relation `1 ≡ −|G|` is lost** above degree 7.

Both failures are the **same 7**: `(ℤ/14)^*` has order `6 = 7−1`, and the null polynomial has degree `7`.
This is the apex prime, seen at the level of the algebra the paper's proof runs on.

## 2. What survives: the low-order moments (= mac-mini's gK8)

The interpolation failure is *graded*. `P` and `Q` agree as functions ⟹ `P−Q ∈ (b_7, 14)`. Since
`deg b_7 = 7`, the **bottom `7` Mahler/finite-difference coordinates** `Δ^j(P−Q)(0)`, `j=0,…,6`, still
vanish: the method **does** pin down the elementary symmetric functions / moments of the data **up to
order 6**, and loses orders `7..13`. mac-mini's S60 finding — "the gK8 concentration is a **low-order
moment-LP led by pairwise `S2`**, on `S0..S4`" — is exactly this: the surviving low-order moments are the
ones the (broken) polynomial method still controls; the irreducible nut is the **top-order tail beyond
`b_7`**. The Clebsch cut-space `(ℤ/2)^4` and the 15 = C(6,2) pairs live on the 6-unit group; they *are*
the multiplicative structure that would have driven Prop 4.1 had 14 been prime. The two frames — algebra
(this note) and moment-LP (mac-mini) — are the **same decomposition**: `b_7` is the degree where control
is lost, order-2 (`S2`) is where it is strongest.

## 3. The sieve side: the c-lift threshold, and why `c=7` is uniquely optimal

The same 7 governs the *constructive* side (THM-573). Generalize the level-7 lift to any `c`:

> **THM-574 (candidate — the c-lift sieve family).** Fix a primitive 13-set `S`. Let `H_c = {s∈S : c|s}`.
> If `c ≤ 7` and `|H_c| ≥ 14 − c`, then `M(S) > 1/14`.
>
> *Proof sketch (same as THM-573).* `H_c = c·P`, `|P|=|H_c| ≤ 12` (primitive ⟹ not all divisible by
> `c`), so by LRC(≤13) there is a `P`-safe phase `v*` with `‖p v*‖ ≥ 1/13`. The `c` lifts
> `t_j=(v*+j)/c`, `j=0,…,c−1`, keep all of `H_c` safe (`‖c p·t_j‖=‖p v*‖`). A speed `w` coprime to `c`
> meets the `c` lift points at `c` equally-spaced phases (spacing `1/c`); the forbidden arc has length
> `2/14 = 1/7`, so it catches `⌈(1/7)/(1/c)⌉ = ⌈c/7⌉` of them. For `c ≤ 7`, `⌈c/7⌉ = 1`, so the
> `13−|H_c|` coprime speeds forbid `≤ 13−|H_c|` lifts; if `13−|H_c| < c`, a lift survives. ∎

The threshold is `|H_c| ≥ 14−c`. Among `c ≤ 7` it is **minimized at `c=7` (threshold 7)**. And `c=7` is
the **largest** `c` with `⌈c/7⌉ = 1`, i.e. the largest `c` whose lift-spacing `1/c` is still `≥` the
forbidden-arc width `1/7`. So:

> **7 is forced because `1/7` is simultaneously the forbidden-arc width *and* the finest lift-spacing that
> still guarantees one survivor per coprime speed.** `c=7` is the unique fixed point `lift-spacing =
> arc-width`. For `c>7` the spacing is finer than the arc, a coprime speed can kill `≥2` lifts, and the
> bound degrades; for `c<7` you need more multiples (`14−c > 7`).

This is the sieve-side twin of §1's algebra-side statement. Same prime, two mechanisms.

## 4. The residual and the CRT path

THM-573/574 (`c=7`) closes `|H_7| ≥ 7`. The residual is `≤6` multiples of 7. A `c=2` lift needs
`|H_2| ≥ 12` (threshold `14−2`) — much weaker, the "level-2 is useless alone" remark. So **no single
`c`-lift closes the residual**; this is structural, not a missed trick. The paper's composite-`k+1`
fallback is exactly to combine lifts at the **prime factors `2` and `7` via CRT** (their `c=2,3` for
`k=11`). For `k=13` the CRT object is the `2×7` grid; the `c=7` lift handles the mod-7 structure
(THM-573), and the mod-2 / 2-adic tower (HYP-2656/2661) must finish — but *jointly*, on the `ℤ/14`
grid, not as two independent single-`c` sieves. The analytic input (Node-3 / measure floor) supplies what
the CRT count cannot.

## 5. The measure floor is the right uniform object (S31ag computation)

`lrc14_largest_arc_witness_denom_kps.py` confirms the correct invariant. For the hardest unbounded-apex
family `{1..11,13,84m}` (`m=1..11`): non-tight, **minimal witness denominator grows `89→929`** (HYP-2866,
the `≈84m+5` ladder), **largest lonely arc → 0** (`Darc 654→1078`), **but lonely measure stays
`≥ 0.0054`** (plateau `≈0.0105 = (6/7)·meas(G_{1..11,13})`, the decoupling floor). For bounded cores
(`max ≤ 40`) the largest arc stays `≥ 0.0047` (`D ≤ 215`, matching the `V*` atlas ≈234). Conclusion:

- **Uniformly positive:** the lonely **measure** (`≥ m_P`). This is OPEN-Q-108 / the gK8 floor — the real
  crux, and it is uniform.
- **Not uniform:** largest arc, minimal witness denominator (both controlled only after peeling).

So **HYP-3088 is corrected:** the largest-arc/bounded-denominator-witness floor holds only for
bounded-core tuples; the genuinely uniform statement is the **measure floor**, and the unbounded apex is
handled by the bounded-speed reduction (Node-3 peeling), exactly as the paper assumes its speeds bounded
(Lemma 2.6, `∏u_i < B_k`). The literal "∀d≥D" Conjecture 7.1 is too strong for unbounded speeds
(aliasing `d|V` kills it); it is a *post-reduction* statement, and our witness route is its engine.

## 6. Net / next

- **New:** THM-574 (c-lift sieve family) with the sharp `c=7` optimality and the "`1/7` = arc width =
  finest one-survivor spacing" reason; the algebraic wall (6 units + degree-7 null poly `b_7`) unifying
  φ(14)=6, the apex 7, and mac-mini's low-order moment-LP; the corrected uniform object (measure floor).
- **Open (unchanged crux):** the measure floor `meas(lonely) ≥ c > 0` over bounded-core covering tuples
  (OPEN-Q-108 / gK8). The polynomial method tells us the *low-order* moments are controlled; the floor
  needs the **top-order tail** (beyond `b_7`) — the same nut as the reflection-Perron / signed
  Erdős–Turán residue tail.

→ THM-574 (new), THM-573, HYP-3087/3088 (3088 corrected), HYP-2656/2661 (2-adic), mac-mini-S60 (gK8),
OPEN-Q-108, arXiv:2604.23906 (Prop 4.1, Conj 7.1), [[lrc14-thread]].
