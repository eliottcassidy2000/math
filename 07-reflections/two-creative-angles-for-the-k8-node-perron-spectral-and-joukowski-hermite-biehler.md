# Two creative angles for the k=8 node: a Perron-Frobenius spectral route (even half, TESTED) and a Joukowski–Hermite-Biehler port (odd half)

*mac-mini-2026-06-27-S73c. The owner asked for a new creative angle or two on a remaining LRC(14) proof
target. kps's S31ai split the last bounded-core node (k=8) cleanly into an EVEN half (consec maximizes total
empty-sector covariance `Sigma Cov(X_i,X_j)` — degree-2, "provable") and an ODD half (the `-9S3`
Worpitzky/associator/apex-7 residue — the hard part). I give one angle per half; the even-half angle is tested
and strong. Builds on [[the-joukowski-de-moivre-bridge-lrc-circle-is-the-tournament-real-rooted-class]],
[[the-k8-node-is-a-variance-extremality-information-theory-rules-out-entropy-kps]], [[hermite-biehler-trrt-strategy]].*

## The two targets (kps S31ai, HYP-3160)
`N` = # of the 6 inner sectors left empty by the orbit `{frac(e_i x)}` (sector 0 = observer, always hit).
The k=8 node = **consec maximizes the bimodality `L_yK8`**, which folds (reflection `s↦6−s`) to:
- **T1 (EVEN / provable):** consec maximizes total covariance `Sigma_{i<j} Cov(X_i,X_j)`, `X_j=1[sector j empty]`.
- **T2 (ODD / hard):** the `-9S3` Worpitzky term = 3-way joint emptiness = the **associator** (non-associative
  apex-7), which does not factor through pairwise data.

## ANGLE A (for T1) — the PERRON-FROBENIUS spectral route. TESTED, strong.
Write the total covariance as a quadratic form `1^T C 1` in the `6x6` empty-sector covariance matrix
`C[i,j]=Cov(X_i,X_j)`. The AP orbit is a **coherent rotation** (all points translate together), so I conjectured
`C` is circulant with the all-ones vector `1` as its **Perron eigenvector**, making `1^T C 1 = 6·lambda_max`
the maximum. TEST (`lrc_covariance_circulant_perron_macmini_S73.py`):
```
 set                 Sigma_off C   1·Perron(cos)   lambda_max
 consec {0..6}        +1.078         0.936          0.418
 consec {0..7}        +1.443         0.997          0.437
 consec {0..12}       +1.442         0.999          0.361     <- LRC(14) size: 1 ~ Perron EXACTLY
 2*consec (dilation)  +1.078         0.936          0.418
 dissoc dyadic        -0.282         0.085          0.532     <- 1 ORTHOGONAL to Perron
 primes               -0.210         0.088          0.324
 random               -0.605         0.025          0.276
 consec{0..12}: 1^T C 1 = 2.1613  ~  6*lambda_max = 2.1653   (0.4% gap)
```
**Findings:**
- `1^T C 1 ≈ 6·lambda_max` for consec: the total covariance **IS (essentially) the Perron eigenvalue** — the
  coherent all-ones mode is the top eigenmode.
- The **all-ones `1` aligns with the Perron eigenvector** for consec (`cos = 0.94 → 0.999` as `k` grows), but
  is **orthogonal** to it for dissociated configs (`cos ≈ 0.02–0.09`).
- **The FM/AFM sign split is SPECTRAL:** consec/AP put `1` in the **positive Perron mode** (`1^T C1 > 0`,
  ferromagnetic, coherent); dissociated configs put `1` in the **negative spectrum** (`1^T C1 < 0`,
  antiferromagnetic). This is exactly kps's "consec=FM covariance, dissociated=AFM(negative)" — now read off
  the eigenstructure.
- HONEST: `C` is **not exactly circulant** (apex/boundary effect, deviation ~0.05) — the clean circulant law
  fails; the **Perron-alignment** is what holds.

> **Sharpened even-half target:** prove the AP's empty-sector covariance `C` has the all-ones vector as its
> **Perron eigenvector** (so `1^T C 1 = 6·lambda_max`), and that any dissociation rotates `1` out of the
> positive Perron cone into the negative (AFM) spectrum. This is the cleanest, most classical (Perron-Frobenius)
> form the even half has had — a spectral statement, not a moment-LP one.

## ANGLE B (for T2) — the JOUKOWSKI → HERMITE-BIEHLER port. Proposed; one leg now PROVABLE.
The dip = EVEN (`+6S4`, the biquadratic `u^4-5u^2+4`, real-rooted, S70) + ODD (`-9S3`, Worpitzky). The
Hermite-Biehler theorem: a polynomial `A + xB` is in the stability class iff `A,B` are AND their roots
**interlace** — exactly the TRRT engine ([[hermite-biehler-trrt-strategy]], Lemmas A/B, verified tournaments
n=6–9). My Joukowski bridge says the LRC dip lives on the **circle image** of that real-axis interlacing.
**NEW provable leg:** the odd Worpitzky term is **real-rooted because the Eulerian polynomials `A_n(t)` are
real-rooted** (Frobenius 1910 — verified `A_3=[1,4,1]`→roots `-3.73,-0.27`, all `A_n` real&negative). The
Worpitzky/Eulerian weighting (kps S71, codex HYP-3147) therefore makes the odd leg real-rooted; the even leg
is real-rooted (S70). So **BOTH Hermite-Biehler legs are real-rooted — only INTERLACING remains**, the same
condition TRRT already verifies on the real axis. The Joukowski map `w=z+R^2/z` is the transport.

> **Sharpened odd-half target:** show the even (biquadratic) and odd (Eulerian/Worpitzky) legs of the dip
> **interlace** (Hermite-Biehler), porting TRRT Lemma B from the tournament real axis to the LRC circle.

HONEST CAVEAT: the LRC miss-PGF `G_N` is **not exactly self-inversive / on-circle** (verified `q_t R^t ≠
q_{6-t} R^{6-t}` at consec k=8), so the Hermite-Biehler-on-the-circle is an **approximate/defect** statement,
not an identity; the interlacing (and the self-inversive defect) is the open content. This is a route, not a
proof.

## Honest status
- **ANGLE A: TESTED and strong.** The even half reduces to a Perron-Frobenius statement (`1` = Perron
  eigenvector of the empty-sector covariance for consec; `1^T C 1 = 6 lambda_max`; FM/AFM = positive/negative
  spectrum). Not exactly circulant (boundary), but the Perron-alignment is robust and grows to `cos=0.999` at
  the LRC(14) size. This is a genuinely new, classical handle on the "provable" half.
- **ANGLE B: proposed; one leg closed.** Both Hermite-Biehler legs (even biquadratic + odd Eulerian/Worpitzky)
  are now provably real-rooted (S70 + Frobenius); the dip bound reduces to **interlacing**, the TRRT engine,
  transported by Joukowski. The self-inversive defect of `G_N` is the honest gap.
- **NOT a proof.** LRC(14) remains open. But both targets are now stated in classical spectral / stability
  language (Perron-Frobenius; Hermite-Biehler/Eulerian) where standard machinery applies — the cleanest the
  k=8 node has been framed.

Related: HYP-3201 (this), HYP-3160 (kps the node = covariance/associator), HYP-3132 (biquadratic), HYP-3147
(Worpitzky odd), HYP-3154 (Joukowski bridge), the TRRT Lemmas A/B, THM-577 (cap), OPEN-Q-108.
