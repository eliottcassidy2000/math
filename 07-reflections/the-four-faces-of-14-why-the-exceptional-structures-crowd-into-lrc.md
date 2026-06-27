# The Four Faces of 14: why the exceptional structures (Farey, two Fano planes, Clebsch) crowd into LRC(14)

*mac-mini-2026-06-27-S60. The owner asked: integrate incoming + past work, understand what REMAINS
in LRC(14), and "look to understand the deep underlying structure which causes these abstract features
to surface" — naming the Clebsch graph and the four Farey variations (the hyperoperation hierarchy on
(numerator, denominator)). This welds four prior threads into one answer: [[the-cut-side-is-classical-clebsch-and-the-permutohedron]]
(S40, Clebsch = cut-space of K5), [[multiplication-is-repeated-addition-the-lrc-hyperoperation-tower-s548]]
(the tower), [[the-proof-lives-at-the-exponential-variation-covering-bound-not-census]] (S59, the
redirect), and HYP-2829/HYP-3085 (the gK8 covering-moment bound). The answer is sharp: **the exceptional
structures are the four hyperoperation faces of the apex denominator 14 = 2·7, and LRC(14) is hard
because it is the first case forced to confront all four at once.***

## The thesis in one line
The torus `R/Z` carries a canonical tower of finite quotients `R/Z ⊇ (1/N)Z/Z ≅ Z/N`. For the apex
denominator `N = 14`, the **four levels of arithmetic on `Z/14`** are exactly the **four hyperoperation
variations on the pair `(numerator, denominator) = (1,14)`**, and each level has a known exceptional
combinatorial avatar. The proof must act on all four; the open part lives at the exponential.

## The four faces (the table that organizes everything)

| Hyperop on (a,b)=(1,14) | value | arithmetic of Z/14 | exceptional avatar | LRC role | status |
|---|---|---|---|---|---|
| **a+b** (H1, additive) | 1+14 = **15** | the GROUP — gaps, three-gap/Steinhaus, Farey/Stern-Brocot mediant | K5(5,10), K3,3(6,9) [a+b=15]; permutohedron walk; AP tight locus | the **census** (tight extremals) | OFF critical path (S59) |
| **a·b** (H2, multiplicative) | 1·14 = **14 = 2·7** | the RING — CRT `Z/14 = Z/2 × Z/7`; units/zero-divisors; multiples | **TWO Fano planes** (the 14 cyclic triples of Paley T7 = 2·PG(2,2), VERIFIED); QR{1,2,4} | the **apex/covering condition** (a runner ≡0 mod 14) | reduction proved (trivial dir.) |
| **b^a / periodicity** (H3, exponential) | `14`-periodicity | the DYNAMICS — a multiple of 14 is `1/14`-periodic; equidistribution | the `1/14`-lift; gamma-trick; Weyl/Erdős–Turán peel | the **covering bound = THE PROOF** | OPEN (CRUX 1 + Node-3) |
| **a^b / iterated** (H4, 2-adic tower) | `2^4 = 16` | the 2-PART — cut-space, Burnside `2^orbits`, Cayley–Dickson | **Clebsch graph** (= cut-space of K5 = folded 5-cube, 16=2^4); the support-six biplane | the descent `14→7→2`; the **S2 pairwise carrier** | the machinery |

The owner's *first* message (K5, K3,3, and `1/14` all have `a+b = 15`) was the **H1 face**. "Think
Clebsch and the four Farey variations" is the instruction to **climb the tower** — and the climb lands
the proof at H3 (exponential), with H4 (Clebsch) as the carrier of its binding term.

## Why these and not others: 7 is the first exceptional prime, 14 stacks the 2-adic tower on it
- **7 is the smallest prime with a non-degenerate finite projective geometry:** `PG(2,2)` = the Fano
  plane (7 points, 7 lines, `|PSL(2,7)| = 168`), which is *also* the octonion multiplication table and
  the Hamming `[7,4]` code. Below 7 the geometry is trivial; **at 7 the Fano/octonion/QR structure
  switches on.** This is the cycle-side, apex-7 obstruction (`H = I(Ω,2)`), and it is why the LRC
  reformulation has **7 sectors**. The "2 copies" in `T7`'s `14 = 2·7` cyclic triples is the CRT 2-part
  already visible inside the prime-7 object.
- **14 = 2·7 stacks the 2-adic tower on the Fano base.** The factor 2 brings the cut-space / hypercube
  / Cayley–Dickson / Clebsch structure (H4). The project's own `E(K_n) = Cut(n−1) ⊕ Cycle(C(n−1,2))`
  decomposition **IS this 2-vs-7 split refracted through the tournament**: the *cut* side (score,
  permutohedron → Clebsch folded cube) is the 2-adic / `a^b` face; the *cycle* side (odd cycles, H,
  apex-7, Fano) is the 7-multiplicative face.
- **LRC(14) is the FIRST OPEN case** (13 speeds). It is the smallest instance whose apex denominator is
  a 2× an exceptional prime — so it is the first forced to confront the 7-Fano cycle structure *and*
  its 2-adic doubling *simultaneously*. Earlier cases (≤7 runners, apex ≤8) never reach the 7-Fano
  threshold-with-a-double. **That confluence is why every small-number coincidence — Farey, the two
  Fano planes, K5/K3,3, Clebsch — crowds into exactly this case.** The exceptional graphs are not
  decoration; they are the combinatorial avatars of the two primes whose product is the apex.

## The weld (where the synthesis pays off concretely)
The open node is **CRUX 1 = the covering-moment bound** (`max_E L_y(E) ≤ scale·cap`, the bounded-core
positivity / OPEN-Q-108): the BOUND the covering proof needs — *not* the census classification (S59).
HYP-3085 (this session) localizes it: the Delsarte duals are **low-order moment functionals**
(`L_yK8 = 10S0 − 10S1 + 10S2 − 9S3 + 6S4`, supported on `S0..S4`), and the binding direction is the
**`+S2` term = the pairwise sector co-emptiness** (verified: S2 leads the consec−wide gap at k=9,10,11;
consec maximizes S2). And `S2` ranges over the `15 = C(6,2) = 2^4 − 1` inner-sector pairs — **the nonzero
vectors of the Clebsch cut-space `(Z/2)^4`** (S40). So:

> **The binding term of the covering-moment bound is an H4 (Clebsch / cut-space) object — the pairwise
> covariance — extremized by the H1 (consec / three-gap / additive) configuration.** CRUX 1 is literally
> an H1-extremal sitting on an H4-carrier: the cross-level weld S548 predicts LRC is made of.

This is the improved argument: prove CRUX 1 as a **low-order (≤4) moment-LP whose pairwise block is
certified by the Clebsch / `2-(16,6,2)` biplane design-Hodge decomposition** (the `4I + 2J` SRG
eigenstructure), with the bounded finite check + single-far periodicity (THM-563) + r≥2-lower (HYP-2829)
supplying the tail, and a finite `−9S3 + 6S4` correction at the tight k=8 row. Comfortable margins
(0.90–1.44) replace the razor-thin p0 dichotomy (0.13).

## What remains, mapped onto the four faces (the honest scope)
- **H1 (census):** the *classification* of tight configs is not logically required (S59) — but the tight
  configs are NOT off the critical path. **kps-S31af correction (HYP-3084):** the covering bound is *tight*
  (`2·{1..13} = {2,…,26}` is covering and `M = 1/14`), so there is **no free margin** — my earlier
  "strictly weaker" framing was wrong. The fix is illuminating: the **dilation `×2` that carries the AP
  into the covering case IS the H2 multiplicative face** (`14 = 2·7`), so the census tight locus reappears,
  dilated, as the equality locus of the *sharp* covering bound. The proof needs that sharp inequality
  (`p0 ≤ cap_k`, with `cap_k` the exact extremal), not the AP/GW dichotomy.
- **H2 (apex/CRT):** the trivial reduction is proved (no multiple of 14 ⟹ τ=1/14 works); the apex-MAJORITY
  branch is **proved** and now **sharpened by kps THM-573 to ≥7 multiples of 7** (the *level-7 lift sieve*,
  a single argument subsuming THM-570/571's ≥7-multiples-of-14 — the 14→7 two-Fano descent, H2→H3 in
  action), and the one-large-speed peel is **proved** (HYP-2906). The dilation `×2` (above) is the same
  H2 operation seen from the equality side.
- **H3 (exponential):** the live core. (i) **CRUX 1** — the bounded covering-core positivity (= HYP-3085's
  low-order/S2 bound); (ii) **Node-3** effective equidistribution (Erdős–Turán) for the unbounded peel
  (HYP-2900); (iii) the induction base (LRC ≤ 7, known).
- **H4 (2-adic):** not a separate obligation but the *carrier/method* — the Clebsch pairwise design-Hodge
  certificate that should close CRUX 1, and the prime-tower descent that drives THM-571.

## Rigorous vs suggestive (discipline)
- RIGOROUS: the covering reduction (H2); the two-Fano fact for T7 (14 cyclic triples, verified); the
  apex-majority/one-large proofs; Clebsch = cut-space of K5 (S40); the moment form and S2-driver
  localization (HYP-3085); `15 = 2^4 − 1` = Clebsch cut-space dimension count.
- SUGGESTIVE (a frame, not yet a theorem): that the *exact* design-Hodge eigenstructure of the 6-sector
  pairwise covariance is the biplane `4I+2J` and certifies the consec extremum (HYP-3085's "to verify");
  the `(numerator,denominator)`→hyperoperation→exceptional-graph dictionary as a *generative* principle
  rather than an a-posteriori organizing table. These are leads with teeth, flagged as such.

## The one-sentence answer to the owner's question
The abstract features surface because **LRC(14) is governed by the apex `14 = 2·7`, and the four
hyperoperation faces of `(1,14)` are the four arithmetic structures on `Z/14` — additive (Farey/three-gap,
the census), multiplicative (CRT into two Fano planes, the apex condition), exponential (periodicity, the
covering proof), and 2-adic (cut-space/Clebsch, the pairwise carrier)** — so the Farey mediants, the two
Fano planes, and the Clebsch graph are not analogies imported into the problem but the **native
combinatorial shadows of the prime 7 and its doubling**, and the proof is the act of moving the binding
obstruction up from the additive face (where it is tight/censused) to the exponential face (where a
multiple of 14 decouples and a safe sector survives).
