# The master cycle-packing polynomial: the spectrum and the OCF are two specializations of ONE object

*monad-explorer-2026-06-15-S5. Builds directly on THM-505 and the reflection
`the-zeta-function-and-the-ocf-read-complementary-halves`. That reflection said the
spectrum and `H` "read complementary halves of the same closed-orbit data" and named the
shared object the **dynamical zeta function** (the trace / power-sum side). This reflection
sharpens "complementary halves" into "**two evaluations of one polynomial**," and shows the
structurally correct shared object is not the zeta function but a **disjoint-cycle packing
polynomial** — because `H` is itself a packing count, and the matching spectral reading is
the elementary-symmetric / Sachs Coefficient Theorem, not the power-sum / trace one. This
answers handoff (3) of the parent reflection: "read `H` as a special value of a tournament
zeta/`L`-type function in the fugacity variable."*

## The object

For a tournament `T` on `n` vertices with adjacency `A`, a **linear subdigraph** `L` is a
set of vertex-disjoint directed cycles (the empty set allowed). Directed cycles in a
tournament have length `≥ 3`. Define the **master cycle-packing polynomial**, one variable
`y_k` per cycle length `k`:

> **Φ(T; y₃, y₄, y₅, …) = Σ_{L linear subdigraph} ∏_{C ∈ L} y_{|C|}.**

`Φ` is the generating function of all vertex-disjoint cycle packings, graded by the
length-multiset of the packing. Grouping by length-multiset `λ`,
`Φ = Σ_λ N_λ ∏_i y_{λ_i}`, where `N_λ` is the number of packings with cycle-length
multiset `λ` (`N_∅ = 1`).

Two very different invariants of `T` are evaluations of this one polynomial.

### Specialization 1 — the spectrum (Sachs Coefficient Theorem)

Set `y_k = −x^{−k}` for **every** length `k`. Then

> **xⁿ · Φ(T; −x^{−k}) = det(xI − A).**

Reading off the coefficient of `x^{n−m}`:

> **[coeff of x^{n−m} in det(xI − A)] = e_m^{signed} := Σ_{L : |V(L)| = m} (−1)^{#cycles(L)}**,

the **signed, all-length** packing count of collections covering exactly `m` vertices.
These coefficients are `(±)` the elementary symmetric functions of the eigenvalues — so
`e_m^{signed}` is **purely spectral**. (This is the digraph Coefficient/Sachs Theorem;
verified here `200/200` at n=3..6 and `120/120` at n=7 by independent packing enumeration
vs. exact Faddeev–LeVerrier char-poly, `master_cycle_packing_monad.py`.)

### Specialization 2 — the OCF (Rédei `H`)

Set `y_k = 2` if `k` is odd, `0` if `k` is even. Then only odd-cycle packings survive and
each contributes `2^{#cycles}`:

> **Φ(T; 2·[k odd]) = Σ_{odd packings} 2^{#cycles} = I(Ω, 2) = H(T)**,

the number of directed Hamiltonian paths (OCF / THM-499). (Verified `H = I(Ω,2)` against
direct Hamiltonian-path enumeration, all samples, n=3..7.)

## The reframe: one polynomial, two loci

So the spectrum and `H` are not "two functions that happen to share closed-orbit data."
They are **the same polynomial `Φ` evaluated at two places**:

| invariant | locus in `Φ` | lengths | sign | weight |
|---|---|---|---|---|
| char. polynomial / spectrum | the **curve** `y_k = −x^{−k}` (all `x`) | all | signed `(−1)^{#cyc}` | length-graded `x^{−k}` |
| OCF `H` | the **point** `y_k = 2·[k odd]` | odd only | unsigned | fugacity `2` |

The "non-spectral defect" of THM-505 is exactly the **gap between this curve and this
point**. The spectrum probes `Φ` only along a one-parameter, all-length, *signed* curve; `H`
probes a single all-*unsigned*, *odd-only* point. Three differences separate them — (a) the
spectrum weights even and odd cycles identically (`y_k = −x^{−k}` for all `k`), (b) it carries
signs `(−1)^{#cyc}`, (c) it grades by total vertices `x^{−Σλ}` rather than by cycle count.
Each difference is a channel through which information visible to `Φ` becomes invisible to the
spectrum.

This is a cleaner statement than "complementary halves." The two are not halves of a
partition; they are two **specializations of a common parent**, and the parent is a
*packing* polynomial — the right category for `H`, which is a packing count, not an orbit
count. The zeta-function picture used the **power-sum** reading of the spectrum (traces
`tr A^k`, primitive-orbit counts `W_k`); the master-polynomial picture uses the
**elementary-symmetric** reading (char-poly coefficients `e_m = signed packing counts`). The
two readings are Newton-dual, but only the elementary-symmetric one puts the spectrum and
`H` in the **same units** (both are disjoint-cycle packing counts), which is why the
comparison becomes transparent.

## The spectral skeleton is Z-linear in the characteristic-polynomial coefficients

THM-505 wrote the spectral skeleton of `H` in the power-sum basis, e.g. at n=7
`S = 1 + 2c₃ + 2c₅ + 4·C(c₃,2) − 4W₆` with `W₆ = (tr A⁶ − tr A³)/6`. In the natural
(elementary-symmetric / Sachs) basis the skeleton collapses to a transparent **integer
combination of char-poly coefficients** `e_m`:

> **n=7:  H = (1 − 2e₃ − 2e₅ + 4e₆) + 4c₆ + 2c₇.**
> **n=8:  H = (1 − 2e₃ − 2e₅ + 4e₆ + 4e₈) + 4c₆ + 2c₇ + 4c₈ − 4·D₄₄.**

(`e_m` = coeff of `x^{n−m}` in `det(xI−A)`; `c_k` = number of `k`-cycles; `D₄₄` = number of
vertex-disjoint pairs of 4-cycles. Verified `2000/2000` at n=7, `1500/1500` at n=8;
the bracket is constant on every one of the 167 sampled cospectral classes at n=7.) The
bracket is *manifestly spectral* — it is a `Z`-combination of characteristic-polynomial
coefficients, no traces or rational `W_k` needed. "Spectral skeleton" now means literally
"`Z`-linear in the char poly."

And the **defect lock** of THM-505 is one Sachs coefficient:

> **e₆^{signed} = D₃₃ − c₆**  (verified n=6,7,8).

This says the signed packing count on 6 vertices fixes the *combination* `D₃₃ − c₆`
(disjoint triangle pairs minus 6-cycles), while `H` needs `D₃₃` and `c₆` *separately*
(`D₃₃` at independent-set level 2, weight `2² = 4`; `c₆` entering via the same pair-defect).
The carrier `c₆` is non-spectral precisely because it is **locked to `D₃₃` by a single
char-poly coefficient** — the carrier is the *unresolved even/odd split inside the Sachs
coefficient*. This is the whole THM-505 mechanism in one line: the spectrum fixes the
signed all-length Sachs coefficients `e_m`; the carriers of `H` are the resolutions of each
`e_m` into its even-vs-odd, signed-vs-unsigned constituents.

## Every face of Φ turns non-spectral at the same wall — n=6

`Φ` has more faces than just `H`. Three natural fugacity-2 packing counts and their signed
versions, probed against cospectral classes (does the invariant ever differ between two
cospectral tournaments?):

| face of `Φ` | first non-spectral at | split classes at n=7 |
|---|---|---|
| `H = I(Ω_odd, 2)` (odd packings, unsigned) | **n = 6** | 47 / 168 |
| `I(Ω_even, 2)` (even packings, unsigned) | **n = 6** | 46 / 168 |
| `I(Ω_all, 2)` (all packings, unsigned) | **n = 6** | 47 / 168 |
| signed odd packing `Σ(−1)^{#cyc}` | **n = 6** | **16** / 168 |
| signed even packing | **n = 6** | 46 / 168 |

(`master_cycle_packing_monad.py`, 6000–9000 samples per `n`; all faces spectral at n ≤ 5.)
Two readings:

1. **The wall is `n = 6` for everything**, because `n = 6 = 3 + 3` is the first size where
   two cycles can be vertex-disjoint at all — the disjoint-pair level `α₂` switches on, and
   `α₂` is the first packing datum the signed-all-length `e_m` cannot disentangle. The
   non-spectrality of `H` (THM-499), of `c₆` (THM-502), and of the even and all faces are
   the **same event**: the first disjoint pair.
2. **Signs partially restore spectrality.** The signed odd face splits only `16` classes vs
   `47` for unsigned `H`. Sign cancellation `(−1)^{#cyc}` moves the odd packing count *closer*
   to the spectral curve (which is fully signed), but does not reach it — the length-grading
   and the even-cycle channel still separate them. The unsigned fugacity-2 evaluation is the
   *most* non-spectral way to read the odd packings; that is exactly why `H` (Rédei's
   unsigned count) is the project's richest non-spectral invariant. **The choice of fugacity
   and sign is the dial that controls how much of `Φ` the spectrum can see.**

## The fugacity is a dial — and `x = −1` is the Euler characteristic of the packing complex

The cleanest evidence that fugacity selects *which* non-spectral functional you read comes
from comparing `x = 2` (the OCF) with `x = −1`. For any `x`, `I(Ω, x) = 1 + Σ_j x^j α_j`,
and the non-spectral content sits in the level-sums `α_j`. At `x = −1`,

> **I(Ω, −1) = Σ_j (−1)^j α_j = −χ̃(Ind Ω)**,

the reduced **Euler characteristic of the odd-cycle packing complex** `Ind(Ω)` (the
simplicial complex whose faces are the vertex-disjoint odd-cycle collections). It is the
`sgn_odd` face of `Φ`. At n=7 (two non-spectral levels, `α₁ ⊃ c₇`, `α₂ ≅ c₆`), the two
fugacities read two different directions of the same 2-D non-spectral content `(c₆, c₇)`:

> `H = I(Ω,2) = [1+e₃+e₅+e₆]·(spectral) + (4c₆ + 2c₇)`  — direction `(4,2)`,
> `χ-count = I(Ω,−1) = [1+e₃+e₅+e₆] + (c₆ − c₇)`  — direction `(1,−1)`.

(`I(Ω,−1) = (1+e₃+e₅+e₆) + (c₆−c₇)` verified `3000/3000`.) Within cospectral classes the
behaviour is fully explained by these directions: `H` splits a class **iff** `Δ(2c₆+c₇) ≠ 0`,
the Euler-characteristic count splits **iff** `Δ(c₆−c₇) ≠ 0` (both consistent with
`0` mismatches over 167 classes). The count splits 16, `H` splits 47, and **47 = 16 + 31**:
the 31 extra classes are exactly those where `(c₆,c₇)` varies but `c₆−c₇` does **not** —
i.e. where `c₆` and `c₇` **covary** (`Δc₆ ≈ Δc₇`). This is the project's own observed fact
that `c₆` and `c₇` co-vary within cospectral classes (46/47, parent reflection). The
alternating sign `x = −1` is *near-orthogonal to the covariation direction* `(1,1)` —
`(1,−1)·(1,1) = 0` — so the Euler characteristic cancels exactly the part of `H`'s
non-spectrality that lives along the dominant covariation, while `x = 2`'s direction
`(4,2)·(1,1) = 6 ≠ 0` amplifies it. **The fugacity is a dial choosing a linear functional of
the non-spectral level-sums; `x = 2` points into the covariation, `x = −1` points across it.**
The topology of the packing complex (its Euler characteristic) is "more spectral" than its
fugacity-2 size for a concrete linear-algebra reason, not a coincidence.

**Sign-damping holds but fades with `n`** (steering from the project owner, POKE-COORDINATION
2026-06-15). The signed `x = −1` count splits fewer cospectral classes than the unsigned `H`
at every `n` tested — n=7: **16 vs 47** (66% fewer), n=8: **656 vs 717** (8.5% fewer) — so
the effect is real but is a *small-`n`* phenomenon, decaying as the non-spectral content
gains dimension and the single alternating direction `((−1)^j)_j` can cancel only a shrinking
fraction of it. Crucially the *unsigned* counts at `x = 1` (total #packings: 47, 713) and
`x = 2` (`H`: 47, 717) split almost identically — **it is the alternating sign, not the
fugacity magnitude, that damps.** Sign-damping is the `n=7` shadow of the general fact that
`I(Ω, x)` is a *one-direction* probe of the `⌊n/3⌋`-D non-spectral level-sum space; only when
that space is 2-D (n=7) can a single sign-choice be nearly orthogonal to all of the
within-class variation.

## Why this matters / the transcendent line

The project's recurring slogan is "the spectrum is mean-field, the OCF is correlation."
The master polynomial says where both live: a tournament is a **gas of disjoint-cycle
packings**, and `Φ` is its full grand partition function with a separate fugacity per cycle
length. The eigenvalues are not a different object — they are `Φ` read along the single
all-length signed curve `y_k = −x^{−k}`, the locus where every length is tuned to one knob
and signs are forced. `H` is `Φ` read at the odd-only unsigned fugacity-2 point. Cospectral
tournaments are those that agree on the *entire spectral curve* yet can disagree at the OCF
point — and they first can at `n = 6`, the moment a second cycle finds room to stand apart
from the first. The "two dimensions" of THM-505 (the linear `⌊n/3⌋` that `H` reads vs. the
`A000009(n) − 3` of the full packing vector `N_λ`) are now two ways of *coordinatizing the
same `Φ`*: `H` collapses `Φ` onto its level-sum marginals (set all `y_k = 2`), while the
multivariate `Φ` itself — keeping the `y_k` independent — is the object whose coordinate
count is `A000009 − 3`. The spectrum, `H`, the even/all faces, the signed variants, and the
full fine packing vector are a single lattice of specializations of one polynomial; the
project has been studying its shadows one fugacity at a time.

## The general Sachs-basis skeleton (derived, verified n≤9)

The `w_m` pattern resolves cleanly. The canonical (simple-cycle-carrier) form is

> **H = S + Σ (carrier)·(±2^{level})**,  with  **S = 1 − 2e₃ − 2e₅ + 4·Σ_{m even, m≥6} e_m**,

and carriers = the simple cycle counts `c_k` (`k≥6`; odd `k` weight `2`, even `k` weight `4`)
plus the overlap/multi-cycle configs (`D₄₄` weight `−4`, `T₃₃₃` weight `+8`, …). Verified:

> n=7: `H = (1−2e₃−2e₅+4e₆) + 4c₆ + 2c₇`  (`2000/2000`)
> n=8: `H = (1−2e₃−2e₅+4e₆+4e₈) + 4c₆+2c₇+4c₈ − 4D₄₄`  (`1500/1500`)
> n=9: `H = (1−2e₃−2e₅+4e₆+4e₈) + 4c₆+2c₇+4c₈+2c₉ − 4D₄₄ + 8T₃₃₃`  (`800/800`)

The **conceptual reason** the skeleton picks up exactly `e₃, e₅` and the *even-m* coefficients:
`e₃ = −c₃`, `e₅ = −c₅` are the spectral level-1 odd cycles (THM-118), entering `H` at level 1
(weight `2`), so `2c₃ = −2e₃`, `2c₅ = −2e₅`. The **even-m** Sachs coefficient `e_{2t}` is the
spectral shadow of the disjoint-**odd-pair** level `α₂` (a disjoint odd pair has `ℓ=2`, sign
`+1`, and even vertex-count `2t`: `(3,3)→6, (3,5)→8, (5,5)/(3,7)→10`), so `α₂` enters `H` at
level 2 (weight `4 = 2²`), giving `+4 e_{2t}` in the skeleton plus the even-cycle carriers.
**Odd-m** coefficients `e_{2t+1}` (`m≥7`) contain a single odd cycle `c_{2t+1}` — itself a
non-spectral carrier — so the whole odd `e_m` goes into carriers, never the skeleton.
*Threshold caveat:* the bare `4·Σ_{even m}e_m` skeleton is exact only while `α₄ = 0`
(`n ≤ 11`); at `n = 12` the quadruple level `(3,3,3,3)→12` switches on and the even-m
coefficient `e₁₂` then mixes the `α₂` (weight 4) and `α₄` (weight 16) shadows, refining the
coefficient — the level structure of THM-505 reappears in the skeleton itself, one rung per
new packing level.

## Handoffs
2. **The even / all faces as new invariants.** `I(Ω_even, 2)` and `I(Ω_all, 2)` are
   non-spectral on the same schedule as `H` but have not been studied. Do they have a
   combinatorial meaning (analogue of Rédei's Hamiltonian-path count)? Is `I(Ω_all, 2)` a
   known graph polynomial of `T`? Their dimension growth law (the `H` analogue of
   `A000009 − 3`) is open.
3. **The signed packing count `Σ(−1)^{#cyc}` over odd packings.** It is "more spectral" than
   `H` (16 vs 47 split classes). Is the *signed* odd-packing count expressible from the
   spectrum plus a strictly smaller carrier set? If so it isolates exactly the part of `H`'s
   non-spectrality that survives sign cancellation — a finer wall than `n = 6`.
4. **`Φ` as a tournament invariant in its own right.** The multivariate `Φ(T; y₃, …, y_n)` is
   a complete-ish fingerprint (it sees all `N_λ`). How does it sit relative to the iso-class
   graph / the `H`-spectrum-as-fingerprint program? Two non-isomorphic tournaments with the
   same `Φ`?
