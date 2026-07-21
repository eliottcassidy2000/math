# NC2's resonance band is the Sheffer no-common-zero (S62/S64) — and cancellation requires mixed signs (a positivity crux)

**death-star-2026-07-21-S87** (HYP-8769). Owner: work to complete NC2, relate it to past work, follow the
insightful threads. NC2's frontier (after codex THM-2014/2017) is the **finite resonance band** (HYP-8766) of the
three-weight $P=Z^p a(s)+b(s)+\bar Z^q c(s)$: the $2r+1$ offsets $-r\le\deg h-r\deg b\le r$. This session connects
that frontier to three of my past threads and adds a positivity observation; **NC2 remains open** at the band.

## The resonance band IS my S62/S64 Sheffer no-common-zero (the "relates to past work")
For the single-straddle $P=Z\,A(s)+B(s)+\bar Z\,c_0$ ($p=q=1,r=2$), charge-0 forces $\#Z=\#\bar Z=i$, and
$$E[P^m]=\sum_{i=0}^{\lfloor m/2\rfloor}\binom{m}{i,\,i,\,m-2i}\,L\!\big(s^{\,i}A(s)^i c_0^{\,i}\,B(s)^{m-2i}\big),\qquad L(s^k)=k!,$$
a sum over **primitive-return channels** $i$ (codex's channels). This is *verbatim* the object of my S64: each
channel a factorial-weighted Laplace moment, the whole a Sheffer/Legendre-recursive sequence in $m$. codex's
named next proof — a **"differential/transseries tower with no common nonzero zero"** (the hyper-Bessel channels
share no common root) — is exactly my **S62 Hermite-no-common-root** (Lean-verified for the base $(1,1)$ case,
`GMC2HermiteNoCommonRoot.lean`) **generalized to the Sheffer hierarchy**, which my **S64** found does *not* close
by the top-term/domination argument once $\deg b\ge2$ — **precisely codex's resonance band** ($\deg h\approx r\deg b$,
the balanced-degree region). So the algebraic (Sheffer, mine) and analytic (hyper-Bessel, codex) attacks are the
*same* crux: **no common nonzero zero of the channel polynomials**. This unifies two independently-developed lines
onto one statement — the honest "relates to past work" the owner asked for.

## Cancellation requires mixed signs — a positivity crux (the S67 thread)
New observation (verified, 1280-example central-offset scan, **zero** two-sided nullcone elements): for a
**positive-coefficient** $P$ every channel $\binom{m}{i,i,m-2i}L(s^iA^ic_0^iB^{m-2i})$ is **positive**, so
$E[P^m]>0$ — it fires at once, NC2 holds *trivially by positivity*. Cancellation (a nullcone element) can only
happen with **mixed-sign coefficients**, where channels of opposite sign balance. That is the resonance band:
codex's "channels preserve factorial slope but destroy phase" is exactly "same magnitude, opposite sign."
**So the resonance band is a positivity-past-the-cancellation-wall problem** (my S67 / klein-S363): the obstruction
is signed cancellation among channels of comparable size, and the natural weapon is a **sum-of-squares /
Hankel-PD reformulation** of the channel sum, not another top-term estimate (consistent with codex ruling top-term
out). Concretely: is there a change of variable in which $E[P^m]$, or its generating function, is a manifestly
non-negative (Bargmann-type, my S84) form on the resonance band? If so, cancellation is impossible and NC2 closes.
This is the thread I'd pull next — it routes codex's analytic band through the repo's positivity toolkit.

## Threads followed, and the map
- **S86 pushforward/null-quadrature:** the whole band lives inside "$P_*(\text{Gaussian})$ has vanishing analytic
  moments"; the channel sum is the analytic-moment expansion.
- **S80 even/odd, toral/radial:** the channel index $i$ is the radial-height grading; the band is the *radial*
  (open, EMP) layer, the toral part being DvdK-closed.
- **THM-438 Wigner/Catalan:** the channels' factorial weights are the free-cumulant/moment weights; the
  "proportional-channel entropy saddle" (codex's central offset) smells like a free-probability rate function —
  a possible closed form for the central offset.
- **klein THM-1790 EMP floor / THM-1810 bosonic:** the band is the bosonic (permanent, no-cancellation) residue;
  the positivity crux above is the precise "why the bosonic side resists."

## Honest scope + computation
NC2 remains OPEN at the resonance band. Contribution: (1) identifying codex's resonance band with my S62/S64
Sheffer no-common-zero (two attacks, one crux); (2) the positivity observation — cancellation needs mixed signs,
so the band is an SOS/positivity problem (S67/S84), a fresh weapon; (3) computational confirmation (1280
central-offset examples, zero nullcone, depths 1–3). No new theorem. Cross-links: codex THM-2014/2017/HYP-8766,
my S62 (Hermite no-common-root), S64 (Sheffer non-closure), S67 (positivity past the wall), S84 (Bargmann),
S86 (null-quadrature), THM-438 (Wigner), THM-1790/1810 (EMP/bosonic), THM-1540/1775. Script
`nc2_resonance_band_sheffer_deathstar_S87.py` (+out). HYP-8769.
