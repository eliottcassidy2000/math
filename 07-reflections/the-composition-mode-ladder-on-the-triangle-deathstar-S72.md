# The composition-mode ladder on the triangle: additive figurate, multiplicative Ham, factorial arborescence — and why {7,21} can only live in the middle

**death-star-2026-07-20-S72** (HYP-8555). Owner: relate the recent Ham-path/arborescence spectrum work
(S70–S71) to the extension frames the project has viewed triangular numbers through (polygonal,
figurate, Faulhaber, …); invent more; and *procedurally* analyze each along many axes. The organizing
find: **every extension frame is a composition operation applied to the same triangular substrate, and
the four modes `+, ×, !, ^` each carry their own representability/spectrum theorem.** The project
mapped the additive quadrant thoroughly; S70–S71 open the multiplicative and factorial ones — and a
density argument shows the `{7,21}` phenomenon can *only* occur in the multiplicative frame.

## The substrate and the four modes

The staircase `δ_{n−2}` has `m = C(n−1,2) = T_{n−2}` tiles — a **triangular number**. Every quantity in
the project is a composition operation applied to this substrate:

| mode | how it builds | representative frame | representability / spectrum theorem |
|---|---|---|---|
| **`+` additive** | sum of `k` terms | `k`-gonal `P(s,n)`; Faulhaber `Σk^p`; simplicial `C(n+d−1,d)` | **Fermat/Cauchy**: every `n` = sum of `≤s` `s`-gonal; **Pollock**: `≤5` tetrahedral |
| **`×` multiplicative** | product of terms | **Ham-path count** `H(T⊕S)=H(T)H(S)` | **spectrum = odds `\ {7,21}`** (S70) — a multiplicative monoid |
| **`!` factorial** | branching (parent-choices) | **arborescences** `Σa_r`, min `(n−1)!` (S71) | `(n−1)!`-bands; forbidden `{4,5,8,11–23,…}` |
| **`^` exponential** | `2^substrate` | tilings `2^{T_{n−2}}` | the hypercube `Q_{T_{n−2}}` on which the `×` and `!` invariants live |

## Procedurally generated grid (uniform axes)

A single analyzer (`04-computation/triangle_frames_procedural_deathstar_S72.py`) takes any frame's
sequence and reports: first terms, growth (finite-difference degree → polynomial deg, else exp/factorial),
and parity pattern. Output:

| frame | first terms | growth | parity |
|---|---|---|---|
| `T_n` triangular | 1,3,6,10,15,21 | poly deg 2 | period-4 `oo·ee` |
| squares `s=4` | 1,4,9,16,25 | poly deg 2 | alternating |
| pentagonal `s=5` | 1,5,12,22,35 | poly deg 2 | period-4 |
| `Σk²` (Faulhaber p=2) | 1,5,14,30,55 | poly deg **3** | period-4 |
| `Σk³=T_n²` (p=3) | 1,9,36,100,225 | poly deg **4** | period-4 |
| tetrahedral `d=3` | 1,4,10,20,35 | poly deg **3** | — |
| simplicial `d=4` | 1,5,15,35,70 | poly deg **4** | — |
| **Ham monoid** | 1,3,5,9,11,13 | (monoid) | **all odd** |
| **arborescence min `(n−1)!`** | 1,1,2,6,24 | factorial | — |
| tilings `2^{T}` | 1,1,2,8,64 | super-exp | — |
| Catalan (triangulations) | 1,1,2,5,14 | exp `~4ⁿ` | — |

The **additive frames are all polynomial** — Faulhaber `p` gives degree `p+1`, simplicial `d` gives
degree `d`: this *is* the project's "figurate degree ladder." The **`×/!/^` frames super-grow**. Only
the Ham frame is **all-odd** (Rédei's parity, the linchpin of `{7,21}`).

## The three sharpest connections

**1. `{7,21}` is the multiplicative analog of Fermat's polygonal representability.** Fermat/Cauchy asks
which `n` are *sums* of few polygonal numbers (additive basis); the Ham spectrum asks which `n` are
*products*/values of Ham counts (multiplicative monoid). Same "which integers does this frame realize"
question, transposed from `+` to `×`. The project had the additive side (Gauss `k=3`, Lagrange `k=4`,
Cauchy general); S70 supplies the multiplicative side.

**2. A density argument forces `{7,21}` into the multiplicative frame — it could not live anywhere else.**
Image densities (verified): triangular numbers `≤X` number `~√X` (density → 0); the Ham spectrum `≤X`
is `X/2 − 2` (density → ½); arborescence values are `~log X` bands (density → 0). A **finite, structural
forbidden set** (all-but-finitely-many realized) can only exist where the frame is **co-finite within
its allowed class**. The Ham frame is co-finite within the odds (`{7,21}` the only holes). The additive
figurate frames are polynomially *sparse* — their "forbidden set" is co-sparse (almost everything),
never finite; the factorial frame is band-sparse likewise. **So `{7,21}` is not an accident of tournaments
— it is the unique signature of the one composition mode dense enough to have a finite exception set.**

**3. The modes are nested, not parallel.** The additive substrate `T_{n−2}` is exponentiated to the
tiling hypercube `2^{T_{n−2}}` (`^`), and the multiplicative `H` and factorial `Σa_r` are *functions on
that hypercube*. So `+ → ^ → {×, !}`: the triangle is summed into a tile count, lifted to a hypercube,
and read by multiplicative/factorial invariants. The recent work sits two floors above the figurate one.

## More frames, procedurally generated (each with mode + representability)

`04-computation/triangle_frames_density_deathstar_S72.py`: centered triangular `1+3T_{n−1}` (ADD, deg 2);
**q-triangular** `[n+1,2]_q` at `q=2` = 1,**7**,35,155 (q-ADD, exp); **doubly-triangular** `T_{T_n}` =
1,6,**21**,55 (ADD², deg 4); **Gauss** `n=Σ₃ triangular` (verified `n≤40`, ADD-basis order 3); **Pollock**
`n=Σ_{≤5} tetrahedral` (DIM-basis); **tournament scores** sum to `T_{n−1}` (Landau's inequalities =
constrained partitions of a triangular number — a *partition* mode I add to the ladder).

*A cross-frame curiosity, flagged as likely coincidence (MISTAKE-197 discipline):* the two Ham-forbidden
values are figurate elsewhere — `21 = T_6` (triangular) and `7` = 2nd centered-hexagonal = q-triangular
`[3,2]_2`. Both are the "second nontrivial term" of another frame. I record it as a coincidence of
position, not a mechanism, unless a map is exhibited.

## Honest status
A synthesis + a reusable procedural analyzer + two verified quantitative axes (growth-degree ladder;
image-density ladder). The composition-mode organizing principle and the density-necessity of `{7,21}`
are the new content; the individual frames and their representability theorems are classical or already
in canon. Actionable spin-offs: is the **tournament-score frame** (partitions of `T_{n−1}`) a fifth
mode with its own forbidden phenomenon? and is there a **multiplicative Fermat** — "every odd `∉{7,21}`
is a product of Ham counts of *strong* tournaments" (the S70 monoid-generator question, restated)?

## Credit
death-star-S70/S71 (the Ham/arborescence spectra this generalizes); the project's polygonal/figurate/
Faulhaber/Zeckendorf/Goldbach additive-basis thread (kps `bsd-hodge-the-polygonal-ladder`, `goldbach-
polygonal-zeckendorf-*`, the shear/Proth catalog T1554), the "everything is the triangle" geometry;
classical: Gauss/Lagrange/Cauchy (polygonal), Pollock, Landau (scores), Rédei.

## Cross-links
S70/S71, `bsd-hodge-the-polygonal-ladder-and-what-merges-kps7`, `goldbach-polygonal-zeckendorf-s501`,
`everything-is-the-triangle`, `04-computation/triangle_frames_{procedural,density}_deathstar_S72.py`, HYP-8555.
