# The moment-problem merge: Carathéodory-Toeplitz rigidity, Griffiths/GKS, and Verblunsky all read the LRC ferromagnetic phase

*mac-mini-2026-06-27-S73d. The owner asked to merge in Carathéodory-Toeplitz moment rigidity and
antiferromagnetic Griffiths/GKS, plus a few self-generated concepts. They land on one object: the FM/AFM
split of the LRC coverage that the Lee-Yang circle (S72/S73) and the empty-sector covariance (S73c) already
found. Four classical lenses now read the SAME phase boundary. Builds on
[[two-creative-angles-for-the-k8-node-perron-spectral-and-joukowski-hermite-biehler]],
[[the-joukowski-de-moivre-bridge-lrc-circle-is-the-tournament-real-rooted-class]],
[[the-antiferromagnetic-tournament]], kps S31ai/HYP-3160.*

## The one object: the ferromagnetic/antiferromagnetic phase boundary
`X_j = 1[inner sector j empty]`, `q_t = P(N=t)` the miss-count law. consec (the AP) is the coherent/FM
config; dissociated sets are the frustrated/AFM ones. **Four classical lenses give the SAME split**
(`lrc_caratheodory_toeplitz_griffiths_macmini_S73.py`):

| lens | FM = consec (k≥8) | AFM = dissociated | what it is |
|---|---|---|---|
| **Lee-Yang circle** (S72) | zeros ON the circle | zeros OFF | Ising partition-function zeros |
| **Griffiths/GKS-I** | `Cov(X_i,X_j) ≥ 0` (✓ +0.010,+0.022) | some `Cov < 0` | ferromagnetic correlation ≥ 0 |
| **Carathéodory-Toeplitz** | `[q_{|i-j|}]` PSD (valid measure) | NOT PSD (−0.58) | trigonometric moment realizability |
| **Verblunsky/OPUC** | `\|α_k\| < 1` (interior) | `\|α_k\| > 1` (invalid) | Schur/Szegő reflection coeffs |

So the LRC covering-bound extremality is a **phase transition**, and the FM phase is exactly: positive
correlations (GKS) ⟺ a valid circle-moment sequence (Carathéodory-Toeplitz, PSD Toeplitz) ⟺ interior OPUC
(`|α_k|<1`) ⟺ Lee-Yang-on-circle. The AFM phase breaks all four at once.

## Carathéodory-Toeplitz moment rigidity (the owner's first cue)
The theorem: a sequence `(c_0,…,c_n)` is the moment sequence of a positive measure on the circle iff the
Toeplitz form `[c_{j−k}]` is PSD; and the **boundary** (rank-deficient Toeplitz) ⟺ a **unique, finitely-
supported (rigid)** measure. This is the moment skeleton of the Lee-Yang circle:
- consec's miss-law `q` gives a **PSD Toeplitz** (a valid circle measure) — the ordered/FM phase; dissociated
  `q` is **not** a valid moment sequence (Toeplitz indefinite). Verified.
- The **rigid boundary** = the AP's coverage condensed onto the **de Moivre / 7th-cyclotomic points** (the
  ideal of [[the-joukowski-de-moivre-bridge-lrc-circle-is-the-tournament-real-rooted-class]]). consec
  *approaches* the boundary as `k→8` (commensurate filling: `maxRefl 0.75` at k=8 vs `0.15` at k=13) but is
  **interior, not at, the boundary** — the rigidity is a limit, not attained. HONEST.

## Antiferromagnetic Griffiths/GKS (the owner's second cue)
GKS-I (ferromagnetic Ising): all correlations `⟨σ_A σ_B⟩ ≥ ⟨σ_A⟩⟨σ_B⟩`, i.e. `Cov ≥ 0`. Verified: consec
(k≥8) has **entrywise non-negative** empty-sector covariance — the FM regime — while AFM/dissociated configs
have negative entries. GKS-II (correlations **monotone increasing in the couplings**) is the route to the even
half: the AP is the **maximally-coupled** FM config (all runners in lockstep), so by GKS-II it should
**maximize** `Σ Cov(X_i,X_j)` — exactly kps's HYP-3160 even-half target and codex's HYP-3200 (consec max
`Σκ₂`, 0/3431). **So the even half = a GKS monotonicity statement.**

HONEST: the empty-sector events are not a *literal* Ising/FKG system, so GKS is here a **structural analogy /
target**, not a turnkey theorem — the covariance positivity is *verified*, and proving it *via* GKS/FKG needs
the ferromagnetic (lattice-supermodular) structure established. That is the remaining content.

## A few more concepts (self-generated)
- **Verblunsky / OPUC (Szegő recursion):** the reflection coefficients `α_k` parametrize the measure; `|α_k|≤1`
  with the boundary `|α_k|=1` the rigid (finite-support) case. Tested: FM ⟹ `|α_k|<1`, AFM ⟹ `|α_k|>1`
  (invalid). This is the *engine* under Carathéodory-Toeplitz, and it ties to the Lee-Yang radius (the α's
  encode the zero locations).
- **Pólya positive-definiteness:** a convex, decreasing, non-negative sequence is positive-definite. consec's
  `q` is decreasing (the ordered/FM shape); the AFM `q` is not monotone — Pólya is the elementary sufficient
  condition for the Toeplitz PSD-ness, the FM order parameter in moment language.
- **Fejér-Riesz factorization:** the non-negative coverage trig-polynomial `= |Schur poly|²`; the Lee-Yang
  zeros are the factor's roots — the half-degree (degree-3 resolvent / Joukowski) lives here.
- **Christoffel-Darboux kernel / Christoffel function:** the reproducing kernel of the moment problem; its
  equilibrium measure on the de Moivre points is the rigid AP measure; the Christoffel function detects the
  support — a quantitative "distance to rigidity."
- **FKG / Harris:** the lattice positive-association that would *prove* GKS-I for the empty-sector events if
  the occupation measure is shown log-supermodular — the cleanest provable route to the covariance positivity.

## Convergence with kps S31al (independent, same two cues — integrated)
The owner gave kps the same prompt; kps's S31al landed the same two merges in parallel, with two results that
SHARPEN this one:
- **Carathéodory-Toeplitz, sharper:** kps proved (exhaustively, 0/3432 bounded k=8 clusters) that **consec
  MAXIMIZES `lambda_min(T)`** — the Carathéodory PSD margin of the Toeplitz `T[j,k]=q_{|j−k|}`. That is far
  stronger than my PSD-yes/no: consec is the **most-interior** moment config (largest distance to the moment-
  cone boundary), and `q_0 = tr(T)/7` makes the cover bound a **spectral inequality on `T`** — Szegő/Schur-
  Cohn machinery, the proof route of choice. (Note: "most-interior `lambda_min`" and my Verblunsky "`|α_k|`
  smallest at large k" are the **same interiority**, two coordinates on it.)
- **Griffiths, honest failure mode:** kps found the naive greedy speed-path to consec is **non-monotone** in
  `Σκ₂` (dips through the antiferro phase, 3/60) — so **plain Griffiths fails on the speed-lattice** because
  the couplings are not free parameters. The right partial order is the coupling/coherence manifold via the
  **random-current representation (Aizenman)**. This is the precise upgrade of my "GKS is an analogy" caveat.
So the merged picture: my four-lens unification (Lee-Yang = GKS = Carathéodory = Verblunsky) + kps's two sharp
results (`lambda_min` extremality; random-current as the correct Griffiths order). My added value over S31al =
the **unification** (one phase, four lenses) + the **Verblunsky test** + Fejér-Riesz / Christoffel-Darboux / FKG.

## Honest status
- **VERIFIED:** the FM/AFM split is read identically by Lee-Yang (on/off circle), GKS-I (Cov sign),
  Carathéodory-Toeplitz (Toeplitz PSD), and Verblunsky (`|α_k|≶1`); consec (k≥8) is FM on all four, dissociated
  is AFM on all four; the FM regime begins at commensurate filling k=8.
- **MERGE/STRUCTURAL:** the even-half target (consec maximizes `Σ Cov`) = GKS-II monotonicity; the Lee-Yang
  circle = the Carathéodory-Toeplitz moment-cone boundary; the rigidity is the de Moivre/cyclotomic ideal.
- **NOT a proof, honest gaps:** GKS/FKG is an analogy until the ferromagnetic lattice structure is established;
  the rigidity is approached (interior `|α_k|<1`), not attained; the k=7 undercrowded edge fails all four. LRC(14)
  remains open. Net value: the even half now sits in two classical, heavily-developed theories (the
  trigonometric moment problem; ferromagnetic correlation inequalities), each with standard tools to try.

Related: HYP-3202 (this), kps-S31al (independent same-merge: λ_min extremality + random-current), HYP-3203 (my Perron/HB angles), HYP-3160 (kps even/odd), HYP-3200 (codex `Σκ₂` census),
HYP-3154 (Joukowski bridge), THM-577 (cap), [[the-antiferromagnetic-tournament]], OPEN-Q-108.
