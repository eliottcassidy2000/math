---
id: HYP-2899
title: LRC(14) is TWO structures -- three-gap rigidity (FINITE, Node 2) + torus equidistribution (ANALYTIC, Node 3); no finite certificate (lcm family); the radical/prime handle; the Mobius/Legendre/Eisenstein modular frame (Eisenstein=E_2 where zeta(-1)=-1/12 lives)
status: SYNTHESIS + verified pieces (three-gap rigidity, lcm refutation, prime-witness); proof ROADMAP (owner's two-node split)
source: mac-mini-2026-06-22-S45 (owner: 3 parity-stratified recursion modes + the two-structure split + the lcm refutation)
related:
  - HYP-+2876  # REFUTED this session (finite certificate; lcm family)
  - THM-523    # q-witness (resonance-killing); kps S31p/r/s
  - the-three-modes-are-parity-stratified-and-lrc14-is-eisenstein-compose-legendre  # kps S31r recursion
  - HYP-2885   # additive-energy extremality (= Node 2)
---

# HYP-2899: LRC(14) is two structures, and no finite certificate exists

## 0. The recursion (building on kps S31r/s)
All three modes live on the half-tiling `h(n) = floor((n-1)^2/4)`. The **Legendre** recursion is the
3-set Venn (sizes A=B=h(n-1), C=D=h(n-2), E=F=h(n-3), G=h(n-4)): corners A,D,B; edges A+B-C=|A∪B|,
A+D-E=|A∪D|, B+D-F=|B∪D|; center A+B+D-C-E-F+G = |A∪B∪D| = h(n) [ODD]. **Eisenstein** (EVEN) is the A&B
edge |A∪B|=A+B-C alone (no D corner, no triple). **Mobius** (ALL) is the union IE skeleton. kps: LRC(14)
= **Eisenstein(even) ∘ Legendre(odd)** composed on the half-tiling.

## 1. The modular frame (interpretive, but the names are apt)
- **Mobius** mode (all sizes) <-> `zeta(s)^{-1} = Σ μ(n)/n^s` -- the multiplicative/IE skeleton (φ=μ*id, S44).
- **Legendre** mode (odd sizes) <-> the quadratic character / Legendre symbol `(·/p)` -> the QR structure,
  the Dirichlet L-function `L(s,χ)`, and the **apex-7** (QR(7)={1,2,4}=Fano).
- **Eisenstein** mode (even sizes) <-> the Eisenstein series `E_2(τ)=1-24 Σ σ_1(n) q^n`, whose
  normalization is `4/B_2=24` and `ζ(-1)=-B_2/2=-1/12`. **The owner's "1+2+3+…=ζ(-1)=-1/12" IS the
  Eisenstein constant.** So composing Eisenstein∘Legendre couples the -1/12 (even) with the QR/apex-7 (odd),
  on the Mobius skeleton -- the modular/L-function structure the LRC lives in.

## 2. The two structures (the owner's split = the two open nodes)
**Structure A -- three-gap (Steinhaus) rigidity = FINITE (Node 2).** Only the AP/consec has <=3 distinct
gap-lengths for ALL x (VERIFIED: consec=3, perturb=4, spread=13). Three-gap rigidity forces all-or-nothing
sector coverage => the AP's empty-sector count is most bimodal => the AP maximizes p0/L_y. So "consec
maximizes" (Node 2 / HYP-2885) has a finite/algebraic handle: **only APs are three-gap-rigid**, plausibly
closable by an AP-hull convex-order majorization.

**Structure B -- torus equidistribution = ANALYTIC (Node 3).** For a committed large speed, loneliness is
realized DYNAMICALLY: the orbit (v_1 t,…,v_k t) is a closed geodesic on T^k, and the committed comb removes
~1/7 by Weyl. Not a single tournament -- the equidistribution of the orbit. Effective Weyl / Erdős–Turán.

## 3. No finite certificate (HYP-2876 REFUTED this session)
S_X = {1,…,11,13, lcm(2..X)} is a primitive covering 13-set; its minimal witness denominator GROWS without
bound (12,12,14,14,14,41,41,… then exceeds 41). The committed speed is ≡0 mod every D<=X, so no small D
certifies it. Hence **a purely finite/combinatorial proof of LRC(14) is impossible** -- an analytic
equidistribution input (Node 3) is irreducibly required. HYP-2876 (mine) + HYP-2864 (codex) capture only
the BOUNDED covering sets.

## 4. The radical/prime handle (the effective hint for Node 3)
Refine the resonance-killing to PRIMES: at t=a/p (p prime), all runners are safe iff p divides NO runner
(nonzero residues mod p are >= 1/p >= 1/13 > 1/14). So **M(S) >= 1/(smallest surviving prime)**, and if any
prime p<=13 survives, M >= 1/13 > 1/14. A counterexample must therefore kill ALL primes {2,3,5,7,11,13}
(a runner divisible by each) AND b=14. The committed speed's **radical** (which primes divide it), NOT its
size, controls the witness: the witness denominator is the smallest prime not dividing the committed speed.
Only the lcm-type family (radical ⊇ {2,…,13}) escapes a small-prime witness -- and equidistribution then
shows the prime's NEIGHBORHOOD survives (the thin comb of the huge speed misses it). nextprime(X)-type
growth = the radical handle.

## 5. The proof roadmap (honest)
LRC(14) splits cleanly by bounded/unbounded:
- **Node 2 (bounded):** kill all primes<=13 + b=14 with bounded speeds => a bounded covering set =>
  dominated by the AP (three-gap rigidity, AP-hull majorization) => M >= 1/14.
- **Node 3 (unbounded):** a committed large speed => equidistribution (effective Weyl, radical handle) =>
  the surviving prime/neighborhood gives M >= 1/14.
Both are required (the lcm family forecloses Node-2-only). The modular frame (Mobius/Legendre/Eisenstein)
is the multiplicative spine; ζ(2) (Farey floor, coprime density 1/ζ(2)) is the Node-2 density; ζ(-1)=-1/12
(Eisenstein) is the even-mode constant. -> kps S31r (recursion), HYP-2885 (Node 2), THM-566 (Node 3).
