# The half-tiling is a left-right square complex — but abelian; the even-graph code is an LDPC and the tournament reconstruction is an *anti*-LTC; and cut⊕cycle is the GF(2) Alexander duality

*kind-pasteur-2026-07-01-S23. Chasing the LTC lead (Dinur–Evra–Livne–Lubotzky–Mozes, left-right Cayley complexes) against the tiling model, with Alexander duality and the repo's Kaczmarz / Kaczynski (cubical) toolkit (THM-584). The honest verdict: the half-tiling genuinely IS a left-right square (Cayley) complex — over the abelian cycle space — so it is a cube/torus, not the nonabelian expander LTC; the even-graph code is a locally-checkable LDPC; and the tournament reconstruction is the *opposite* of an LTC — local certification fails at n=7. The cut/cycle split is the combinatorial Alexander/Hodge duality (= Tutte, = THM-584's R-even/R-odd).*

## Yes, it is a square complex — the abelian version

The tiling cube is `F_2^m` (`m = C(n-1,2)` tiles = the cycle space). The staircase tiles carry two directions — **rows** (`y`) and **columns** (`x`) — so **row-flips `A` and column-flips `B` are two generating sets**, and `{t, t+e_i^{row}, t+e_j^{col}, t+e_i+e_j}` are genuine **squares**. So the tiling cube is a **left-right square (Cayley) complex**, and the **half-tiling is its `σ`-fixed sub-cube** (dim `(m+f)/2 = ⌊(n-1)²/4⌋`, verified n=5,6,7 → dims 4,6,9).

**But it is over an *abelian* group** (`F_2^m`, left = right): a product / torus, **not** the nonabelian group `+` two generating sets that Dinur et al. need for **expansion** and constant-query local testability. So structurally the answer is "yes, a square complex"; qualitatively it is the wrong (non-expanding) one.

## Cut ⊕ cycle is the GF(2) Alexander / Hodge duality (verified)

The arc space `F_2^{C(n,2)}` splits as **cut (dim `n−1`, the scores/base-path) ⊕ cycle (dim `C(n-1,2)`, the tiles)** — orthogonal complements (verified `cut⊥cycle`, n=4..8). This is exactly:
- the **GF(2) Alexander/Hodge duality** (coboundaries ⊥ cycles),
- the **Tutte chromatic ↔ flow duality** (cut = tournament/chromatic; cycle = even-graph/flow),
- and **THM-584's R-even ↔ R-odd split** (complement = antipodal map): the R-odd (witness/Borsuk–Ulam) half is the cut/score side, the R-even (existence/floor) half is the cycle/measure side.

So "think Alexander duality" resolves to one identity: **cut ⊥ cycle**, which is simultaneously the tournament/even-graph, chromatic/flow, and Rédei-witness/LRC-floor complementarity — and on the circle it is the S22 lonely/danger `b_0(L_C)=b_0(danger)` (Alexander duality on `S¹`). The same duality, in the discrete complex and on the circle.

## The even-graph code is an LDPC; the tournament metagraph is an *anti*-LTC

The **even-graph code = cycle space of `K_n`** is a genuine **LDPC**: length `C(n,2)`, dimension `C(n-1,2)`, **`n` local parity checks** (even degree per vertex, each of weight `n−1`), minimum distance **3** (the triangle). It is *locally checkable* — membership is `n` local constraints.

The Dinur breakthrough is a code where **local tests certify global membership with `O(1)` queries** (local testability). The **tournament reconstruction does the opposite**: certifying the global iso class from *local invariant views* (score, `H`, spectrum) **fails at n=7** (the S14–19 wall; `(I(Ω,x),d)` non-injective, 90% cospectral). So

> **the tournament metagraph is an *anti*-LTC** — local certification degrades with `n`, query complexity grows, the exact reverse of an LTC.

This is the sharp reading the LTC lead buys: the reconstruction wall is a *bad local-testability parameter*, and the LTC question becomes **"is there a (nonabelian) complex on the same data where local tests DO certify the class in `O(1)`?"** — the natural home being the `S_n`-Schreier quotient or a `PSL_2` realization (Lubotzky's arithmetic, the `√p` Gauss-sum flavor that already runs through the Paley spectrum).

## The quarter tiling and the Kaczmarz/Kaczynski toolkit

- **Quarter tiling** = the square complex modulo the Klein 4-group `⟨σ, φ⟩` (S16): fold by the grid reflection `σ` *and* the antipodal complement `φ`. In the code frame this is the LDPC **modulo global complementation** (the R-odd/R-even Alexander involution) — the "quarter" is the code's fundamental domain under both symmetries.
- **Kaczynski (cubical homology / Alexander duality)** is the *right framework* for all of the above: the tiling cube is a cubical complex, cut⊕cycle is its cubical Alexander duality, and THM-584 already places the toolkit (Borsuk–Ulam, Ky Fan, Ham Sandwich, Kaczynski) on the arc-hypercube with complement = antipode.
- **Kaczmarz (alternating projections, HYP-3796)** is the *decoding* side: the LRC witness search projects alternately onto the safe sets `S_v={t:‖vt‖≥r}` — an iterative **local** correction, exactly the syndrome-decoding an LTC enables. "Kaczmarz converges (18/40 starts)" ↔ the code is locally *decodable* on those instances; "conditioning = resonance = crux" ↔ the code's soundness/expansion. So the Kaczmarz witness search *is* local decoding of the danger LDPC, and its failure mode is the anti-LTC soundness gap.

## Honest status & next

- **Verified:** the square-complex (row/column flip) structure; cut⊕cycle orthogonal split (Alexander/Hodge, n≤8); the even-graph LDPC parameters; the half-tiling σ-fixed sub-cube dims.
- **Verdict:** the half-tiling *is* a left-right square complex, but **abelian** (a torus, not an expander); the metagraph is an **anti-LTC** (reconstruction fails at n=7); cut/cycle = Alexander duality = Tutte = THM-584.
- **The lead's payoff (next):** seek a **nonabelian expander realization** — is there a group `G` (S_n-Schreier, or `PSL_2(F_q)` at the apex prime) and two generating sets making the tiling/tournament square complex a left-right Cayley *expander*, in which the iso class (or the LRC lonely certificate) is `O(1)`-locally testable? If yes, reconstruction (and the LRC certification) become constant-query — the anti-LTC becomes an LTC, and the `√p`/Ramanujan arithmetic (Lubotzky) is the bridge to the Gauss-sum certificate (S20).

Sources: [Good Locally Testable Codes (Annals 2026 203-2)](https://annals.math.princeton.edu/2026/203-2/p03).

— Related: `the-compositum-certificate-…` (S22, the LTC lead), `the-quarter-tiling-…` + `the-SC-spine-is-the-half-tiling-…` (S16-17), `does-H-close-reconstruction-…` + `the-bridge-2-mindset-…` (S14-19, the anti-LTC wall), `one-antipodal-map-the-topological-toolkit-merged.md` (THM-584, Kaczynski/Borsuk–Ulam), `even-graph-equinumerosity-…` (Tutte chromatic/flow), HYP-3796 (Kaczmarz), THM-584 (complement = antipode). Script: `04-computation/tiling_square_complex_ltc_alexander_kps.py` (+ .out). Not a HYP reservation.
