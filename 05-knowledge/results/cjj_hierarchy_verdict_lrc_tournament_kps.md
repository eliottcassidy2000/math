# CJJ complete-LP-hierarchy verdict for LRC + tournament extremality: the subspace-vs-coset dichotomy

**kind-pasteur-2026-06-21** (the user's CJJ lead, arXiv:2211.01248/2112.09221). A consolidated,
honest verdict from this session's three-thread workflow + direct checks. Converges with mac-mini
(HYP-2750/2753) and codex (HYP-2751 signed-Tanner).

## The one structural fact

**CJJ completeness/integrality requires the optimizer to be an `F_q`-SUBSPACE** -- equivalently
(across all four CJJ views): closed under linear combination (view c) <=> higher Krawtchouk/
interaction moments are DETERMINED by lower ones (views a,d: the SoS/Mobius lift self-tightens
around the subspace, vacuous past a finite level). This is the property that makes the hierarchy
collapse to a closed level and the LP polytope integral.

## The two extremizers, split by it

| | LRC: consec = AP | Tournament: H-max = Paley |
|---|---|---|
| Object | additive **coset**, Freiman dim 1 | QR cyclic **code** = ideal in `F_p[x]/(x^p-1)` |
| Subspace? | **NO** (translate, not a subgroup) | **YES** (genuine `F_p`-subspace) |
| CJJ linear lift | **collapses** (HYP-2744, S16) | certificate-existence **closes** (MacWilliams `[7,4]->[7,3]` simplex) |
| The catch | the collapse is partly the **dilation symmetry** (HYP-2754); a per-`k` signed cut PROVES consec-max at `k=8` exactly (HYP-2752) but is **NON-UNIFORM** in `k` | the certifiable functional is **wrong**: `theta'`/Delsarte bounds code DISTANCE, but `H=I(Omega,2)` is **nonlinear** (deg-`m` elem-symm, THM-134); `theta_LP(Omega)` **anti-tracks** `H` (Paley 40 > non-Paley 29.5 yet Paley has higher `H`); and Paley is H-max only `p<=11`, `p=3 mod 4` |

## What this session established (verified)

- **HYP-2754 (k=8):** the moment-LP non-pinning of consec is EXACTLY the dilation tie `consec` vs
  `2*consec` (THM-531); consec's **dilation-orbit is the unique measS7-maximizer**. => the right
  hierarchy factors the **AFFINE group (translation+dilation)**, not just translation+linear-comb;
  after quotienting, AP (one affine generator) is the natural extremizer.
- **HYP-2755 (verified, corrected a buggy thread):** Paley=QR is the circulant H-maximizer
  `p=7 (189), p=11 (95095)`; the linear-code structure gives certificate-EXISTENCE -- but NOT a
  proof of H-extremality (the `theta'` functional is wrong for the nonlinear `H`).
- **HYP-2752 (thread B):** a low-level (`R<=3`) signed Boolean/type cut certifies consec-max at
  `k=8` (exact rational, 0 violations / 319 stratum + 3112 off-stratum) and `k=9`, but is
  non-uniform in `k` and its VALIDITY is non-structural (as hard as consec-max itself).

## The two walls (honest)

1. **LRC:** the consec-max signed cut EXISTS at every fixed `k` (low level) but the **uniform-in-`k`
   validity** is the wall -- the per-`k` LP-tuned cut's validity is itself the consec-max difficulty.
   The affine-hierarchy (HYP-2754) is the right frame; the residual is uniform AP-extremality on the
   full-residue stratum.
2. **Tournament:** `H` is a **nonlinear** spectral functional, not a code-distance enumerator, so
   CJJ/Delsarte (which certify distance) cannot bound it; `theta'` anti-tracks `H`. The only candidate
   is the fugacity-2 **Lasserre hardcore moment hierarchy** on `Omega(T)` -- and it can only work
   `p<=11` (where Paley is H-max).

## Net

The CJJ hierarchy is the correct LENS: it explains the subspace-vs-coset split, certificate-existence,
the affine reframing, and proves consec-max at fixed small `k`. It does NOT give a uniform proof of
either extremality -- the LRC by non-uniform validity (AP=coset), the tournament by the nonlinear
`H`-functional (Paley=subspace but `H` is not its distance). The two open extremalities, identified as
"the same `theta'` problem" (SESSION-LOG L9992), are split by linearity: Paley clears the
certificate-existence bar that AP fails, but neither clears the certificate-VALIDITY bar that the
nonlinearity (AP's affine-coset, `H`'s elementary-symmetric) imposes. -> HYP-2744/2749/2752/2754/2755,
HYP-2726/2738/2602, THM-126/134/531, OPEN-Q-108, the user's CJJ papers.
