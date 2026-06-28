# The Joukowski bridge: tournament TRRT and the LRC cover bound are ONE Hermite–Biehler theorem

*kind-pasteur-2026-06-27-S31am. The owner asked to merge Joukowski, Hermite–Biehler, and Perron–Frobenius
into the LRC(14) effort. These were converging across the team in real time: oracle's TRRT strategy
(tournament real-rootedness via `I=A+xB` + Hermite–Biehler interlacing, Lemmas A/B verified to n=9),
mac-mini-S73c's port (Perron for the even half + Joukowski→Hermite–Biehler for the odd Worpitzky/Eulerian
leg), and my two-maps (HYP-3099: tournament = real axis, LRC = circle) + cyclotomic Joukowski (HYP-3162).
The merge is clean: **the LRC cover bound is the Joukowski image of the tournament TRRT — the same
Hermite–Biehler interlacing theorem conformally transported from the real axis to the circle — and
Perron–Frobenius is the leading mode on both.** (Note: HYP-3201 collided — mac-mini's Perron/HB and my
Carathéodory–Toeplitz both took it; this synthesis is HYP-3210, and the two angle-sets should be merged
under one id by whoever consolidates.)*

## The three tools, one picture
- **Hermite–Biehler** (1877): `P = A + xB` is Hurwitz/real-rooted ⟺ `A,B` real-rooted **and `B` interlaces
  `A`**. The interlacing IS the root-location certificate.
- **Joukowski** `w = z + 1/z`: the conformal 2-to-1 map sending the unit circle to `[−2,2]` (and `|z|=R`
  to `[−2,2]` after R-normalization). It turns "zeros on a circle" into "zeros on a real segment."
- **Perron–Frobenius**: a non-negative matrix has a simple positive dominant eigenvalue with a positive
  eigenvector — the leading collective mode.

## The two maps are two sides of Joukowski (the merge)
| | object | Hermite–Biehler lives on | leading mode |
|---|---|---|---|
| **Tournament (MAP 1)** | `I(Ω,x)=A(x)+xB(x)` | **real axis** (negative reals) — TRRT, oracle Lemmas A/B | spectral radius = deepest root |
| **LRC (MAP 2)** | coverage / the `L_y` dual | **circle** — via Joukowski → real `[−2,2]` | de Moivre ground angle `−1.8019` |

**Joukowski `w=z+1/z` is the bridge:** it carries the LRC circle-Hermite–Biehler onto the tournament-style
real-axis-Hermite–Biehler. The **de Moivre angles `{2cos(2πj/7)}`** are the image of the 7-fold ideal
(HYP-3162) — the fixed interlacing skeleton. So the LRC cover bound and the tournament TRRT are **the same
interlacing theorem**, read on the two sides of the Joukowski map. mac-mini-S73c's "Joukowski→Hermite–
Biehler port of the odd Worpitzky leg" is exactly this bridge applied to one leg; the picture says it is
not a port of a tool but a **conformal identity of the two problems**.

## What ports across (the proof value)
oracle's TRRT machinery is *verified to n=9* and reduced to two combinatorial lemmas. Under Joukowski it
becomes the LRC cover-bound machinery:
- **`I = A + xB` (deletion–contraction)** ↔ the LRC dual's **even ⊕ x·odd** split (mac-mini: biquadratic
  `A` ⊕ Eulerian/Worpitzky `xB`). Both legs are **exactly real-rooted** (the biquadratic by ℚ-collapse
  HYP-3132; the Eulerian by the classical Frobenius/Eulerian real-rootedness — this is why mac-mini's odd
  leg is rigorously real-rooted, sidestepping the coverage PGF's off-circle defect).
- **Lemma B (interlacing)** ↔ the LRC's **even/odd interlacing** = the cover bound. This is the shared
  crux: if the biquadratic and Eulerian legs interlace, the dual has the right root structure and consec
  is extremal.
- **Lemma A (existence of a good `C*`)** ↔ the existence of the reflection-`s↦6−s` decomposition that
  *produces* the even/odd legs (always available — the apex involution).

> **The single inherited target:** prove the **even (biquadratic) / odd (Eulerian–Worpitzky) interlacing**
> for the LRC k=8 dual — the Joukowski image of TRRT Lemma B. TRRT-B is verified to n=9 on the tournament
> side; its Joukowski'd LRC version is mac-mini's open interlacing. **One lemma now covers both problems.**

## Perron–Frobenius is the leading mode (and its honest limits)
- **Covariance side (mac-mini):** for consec the empty-sector covariance `C` is **entrywise non-negative**
  (ferromagnetic, all 15 `Cov>0`, HYP-3161) ⟹ **Perron–Frobenius applies**: a positive dominant
  eigenvector. VERIFIED: it is `≈99.7%` aligned with the **uniform vector `1`** (the ferromagnetic
  collective mode), `λ_max≈0.437`, `1ᵀC1 ≈ 6λ_max` (2.612 vs 2.623). HONEST: not *exactly* uniform — the
  **anchor** (the stationary runner, sector 0) breaks the 7-fold symmetry, so `1` is the Perron eigenvector
  only to `O(anchor)`. For antiferromagnetic configs (`k≤5`, some `Cov<0`) `C` is not entrywise
  non-negative and Perron–Frobenius does **not** apply — the spectral reason the ferromagnetic transition
  (HYP-3161) is the boundary of the Perron regime.
- **Root side:** the **spectral radius** (largest |zero|) of the tournament `I(Ω,x)` (the H-maximizer's
  zero hugging 0, HYP-3099) and the LRC coverage's dominant Joukowski mode (locked to the de Moivre ground
  angle `−1.8019`, HYP-3162) are the same "leading mode" on the two sides. **Perron–Frobenius eigenvalue =
  the deepest interlacing root = the dominant collective mode**, on both maps.

## The merged spectral toolbox for the node (consolidating the angles)
The k=8 node now has a coherent spectral attack, all on consec:
1. **Carathéodory–Toeplitz** (mine): consec maximizes `λ_min(T)` of the coverage-moment Toeplitz
   (EXHAUSTIVE over 3432 bounded k=8) — the moment-cone interior extremal; Szegő tools.
2. **Perron–Frobenius** (mac-mini): consec's covariance is ferromagnetic ⟹ uniform Perron mode; the
   transition is the Perron boundary.
3. **Hermite–Biehler / TRRT** (oracle + mac-mini, bridged here): even/odd interlacing = the cover bound =
   Joukowski'd TRRT-B.
4. **Joukowski / de Moivre** (mine): the bridge + the cyclotomic ideal skeleton.
These are four faces of one statement: **consec is the spectrally-extremal (Perron-dominant,
Toeplitz-interior, Hermite–Biehler-interlacing) configuration**, and the cover bound is the
real-axis↔circle (TRRT↔LRC) Joukowski identity.

## Net / honest status
- **MERGE (new):** Joukowski identifies the LRC cover bound with the tournament TRRT as one Hermite–Biehler
  interlacing theorem; the verified TRRT Lemmas A/B port to the LRC even/odd (biquadratic/Eulerian) legs;
  Perron–Frobenius is the shared leading mode (ferromagnetic uniform on `C`, de Moivre ground angle on the
  zeros).
- **CAVEATS:** the coverage PGF itself is only near-circle (HB is rigorous on the *dual* legs, which are
  exactly real-rooted); the Perron eigenvector is uniform only up to the anchor; HYP-3201 collided
  (consolidate the Perron/HB and Carathéodory/Toeplitz angle-sets).
- **TARGET:** the even/odd interlacing of the k=8 dual = Joukowski'd TRRT-B = the one lemma that closes
  both the tournament real-rootedness and the LRC cover bound.

→ HYP-3210 (this), HYP-3201 (collided: mac-mini Perron/HB + kps Carathéodory/Toeplitz), HYP-3099 (two
maps), HYP-3162 (Joukowski/de Moivre), HYP-3160/3161 (parity/ferromagnetic), TRRT strategy (oracle
hermite-biehler-trrt-strategy), mac-mini-S73c, Eulerian/Frobenius real-rootedness, Szegő/Carathéodory,
OPEN-Q-108.
