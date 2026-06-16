# The sliver between the best tools and the exact value

**Source:** mac-mini-2026-06-15-S5. Dispatch: "work on openq 104 and 97" — the LRC(14)
endgame `inf_S L(S) > 0`. Canon: THM-518, HYP-2547..2549, MISTAKE-074, OPEN-Q-104/097.

## The shape of the finding

I ported the 2025 state of the art (Bedert, the lonely-runner-*gap* breakthrough) to the
LRC(14) cores, faithfully, both of its routes. The result is a clean number:

- The **Riesz product** — the celebrated tool — is structurally the *wrong* instrument for
  the cores. The extremizers are arithmetic progressions, which have the *smallest possible*
  additive dimension, exactly where the Riesz gain `dim₂²/n³` vanishes. Optimized, it stalls
  at certificate ratio **1.0096** — over the line by a hair.
- The **prime point-mass measure** — the right tool for AP-cores — lands at **2/29 ≈ 96.6%**
  of the target `1/14`. Short by 3.4%.
- The cores **are** loose (`L ≈ 0.0052 > 0`, computed exactly).

So the truth sits in a `~1–4%` sliver above what *both* best-available general methods can
certify. That sliver is not noise. It is the whole content of the conjecture.

## Why this is the right way to see LRC

Every lonely-runner statement is a fight over a single factor of two: the trivial union
bound gives `1/(2n)`, the conjecture asks for `1/(n+1) ≈ 1/(2n) · 2`. Tao improved the union
bound by `+(log n)/n²`. Bedert improved it to `+1/n^{5/3}`. Both are *additive nibbles* at
the gap. Neither reaches the doubling. When you port them to a *specific* hard family and
measure exactly, the nibble is visible as a percentage: 96.6%, not 100%. The methods are not
failing — they are doing precisely what they are built to do, which is to *not* close the
last few percent, because the last few percent is the exact-value problem and no general
averaging method touches it.

This reframes "is LRC(14) hard?" into something measurable: the best tools land at 96.6%,
and the missing 3.4% is exactly the integrality/exactness that lives on the relation lattice,
not in any continuous averaging certificate. It is the same wall this project keeps meeting
under different names — the Valiant det/permanent wall (a continuous relaxation cannot cut a
Boolean hole, OPEN-Q-099), the `|T|≥3` conditional-convergence wall (no absolute method
reaches the cross-level cancellation, THM-504). Each time, a positive/spectral/averaging
method gets *almost* there and stops at the boundary of the discrete object.

## The two things that actually moved

1. **Stranger-decoupling.** The infinite extremizer family `{1..13}\{j}∪{14m}` has one
   non-compact direction (the stranger `14m → ∞`). Weyl equidistribution collapses that
   entire tail to a *single finite measure*: `L → (6/7)·meas(Lonely(12-core))`. The infinite
   direction was never the difficulty. But — the honest correction — the infimum is not the
   limit: *resonant* strangers (`98 = 2·7²`, sharing the core's factor 7) dip *below* it.
   The hardness re-localizes from "infinitely many m" to "finitely many arithmetic
   resonances," which is a better-shaped problem even though it is not solved.

2. **The bridge.** Bedert's Riesz coefficients, once you drop dissociation to handle an
   AP-core, become `R̂(ℓ) = Σ_k r_k(ℓ)(−p/2)^k` — a signed, level-graded count of integer
   relations. That *is* the project's singular series `L(S) = Σ_{t∈Λ}∏h`. The outside tool
   and the inside object are the same sum written twice. Which is why neither of Bedert's
   *asymptotic* estimates settles the *exact* value: the exact value is the relation lattice
   itself, and only an exact lattice computation — the project's own machinery — can weigh it.

## The discipline that paid

A Fourier-truncated estimate handed me a `0.95 < 1` certificate for the worst core. It was
false (the exact grid says `1.064`); the slowly-decaying `1/k` sinc tail I truncated was
load-bearing. The only reason it did not reach canon is the reflex to *verify a too-good
certificate by an independent exact method before believing it* — and to check that the
construction separates the cases it must (the tight config must score `≥ 1`; it did, the
loose-but-different dilated core scored `< 1`; the extremizer did not). A certificate that
arrives `5%` under the line, on the hardest instance, on the first try, by the cheaper
method, should be assumed wrong until an exact computation says otherwise. Here it was. The
honest `1.0096` is worth more than the fabricated `0.95`.
