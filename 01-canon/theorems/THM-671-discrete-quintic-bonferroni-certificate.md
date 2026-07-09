---
id: THM-671
title: The discrete quintic Bonferroni certificate (B5-cert) — THM-604's depth-5 truncation ported to Z_q kill sets gives LM(q) ≥ B5(S,q) = Σ_{p≠0} f₅(C(p)), an EXACT-integer one-pass computation whose positivity is a decidable loneliness certificate (t = p/q via the pair-sum dispatch); the iid floor B₅(13) = 2052/7⁵ ≈ 0.1221 > 0 explains why the depth-1/2 union ledgers (C1, C4-Hunter) adversarially broke while the truth stayed fat (the cubic already fails at 13 constraints); VERIFIED: ≥1 certified modulus on 100% of tested covering instances at every scale V = 30..280, including the adversarial C1∪C4-killers, the 7-structured @91 cluster, and the covering equality extremal 2·{1..13}
status: PROVED (parts 1–3: the truncation identity is pointwise combinatorics; the dispatch is kps-S114's sorry-free Lean; the moment form is exact bookkeeping) + PROVED-CONDITIONAL (part 4: the iid floor under discrete depth-5 resolution, THM-604's hypothesis pattern) + VERIFIED (part 5: the census). OPEN (part 6): the a-priori supply — every covering 13-set with Vmax > V₀ has a resolved modulus in (Vmax, 2Vmax]; this is HYP-5732's aggregated resonance-dispersal + the E3-rigidity dichotomy, the ONE remaining item of the aggregated modular route.
source: klein-2026-07-09-S210 (HYP-5758; the owner-directed discrete moment-LP port)
depends_on:
  - THM-604   # the depth-5 truncation identity + iid values B_D(j) (the continuous template)
  - THM-668   # pair-sum/residue-band witness structure (mac-mini)
  - kps-S114 LRCPairSumDispatch  # mreach_ge_of_pairsum_band, sorry-free (the Lean consumer)
  - kps-S115 LRCCovering966      # the [1,18] native_decide base case
related:
  - HYP-5732  # the aggregated modular supply route (this is its crux instrument)
  - HYP-5731  # S209: supply abundant; C1/C4 adversarially breakable — the gap this closes
  - HYP-5730 (mac-mini)  # C1/C2/C3 live-ruler certificates (depth-1 ledger family)
  - THM-660/661          # the continuous moment-LP floors (the ported machinery)
  - LEM-015 + kps-S114 LRCSchurRigidity  # E3 ≤ C(k,2), equality iff dilated interval — the dichotomy dial for part 6
---

# THM-671 — the discrete quintic Bonferroni certificate

**Setting.** S a 13-speed set, q > Vmax a modulus, c = ⌈q/14⌉. Residue r is SAFE iff
14r ≥ q and 14(q−r) ≥ q (the closed middle band [c, q−c]; boundary = the M = 1/14
equality witnesses). Kill sets B_l = {p ∈ Z_q : v_l·p mod q unsafe}; merge l ~ k iff
v_l ≡ ±v_k (mod q) (identical kill sets — the THM-668 conjugate-pair mechanism);
C(p) = # killing classes at p. A p ≠ 0 with C(p) = 0 is a LIVE multiplier:
**t = p/q satisfies ‖v_l t‖ ≥ 1/14 for all l, i.e. M(S) ≥ 1/14 directly**
(`LRCPairSumDispatch.mreach_ge_of_pairsum_band`, sorry-free). LM(q) = #{p ≠ 0 : C(p) = 0}.

## 1. The certificate (proved; pointwise combinatorics)

For integers x ≥ 0 and odd D: 1_{x=0} ≥ f_D(x) := Σ_{d≤D} (−1)^d C(x,d) (THM-604's
truncation; the remainder is an alternating tail of fixed sign). Summing over p ≠ 0:

> **LM(q) ≥ B_D(S,q) := Σ_{p≠0} f_D(C(p))** — an exact integer, computable in ONE
> O(k·q) pass (mark each class's kill fibers, histogram the coverage, evaluate f_D).

Hence **B₅(S,q) > 0 is a decidable loneliness certificate**: it exhibits LM(q) > 0,
i.e. a live p, i.e. M(S) ≥ 1/14. `native_decide`-shaped: integers only.

## 2. The moment form (the LP connection)

B_D(S,q) = Σ_{d≤D} (−1)^d S_d(q), S_d(q) = Σ_{|T|=d} |∩_{l∈T} B_l \ {0}| — the exact
discrete joint-kill moments. The Bonferroni coefficients are the corner of the
moment-LP polytope (THM-660/661's degree-D minorant with c_d = (−1)^d binomials);
running the full LP over the exact S_d can only improve the floor. This is the
discrete port of the density-floor moment machinery: same minorant logic, exact
integer moments, no Farey integration needed.

## 3. The depth ladder explains the certificate history

iid values at 13 classes, kill density 1/7 (THM-604's table): B₁ = −0.857 (union
bound: hopeless — mac-mini's C1 fires only on gcd-structured moduli), B₃ = −0.099
(cubic FAILS at 13 — Hunter/C4-type depth-2 corrections cannot cross zero
generically; this is precisely why klein-S209's adversary broke C1∪C4 to zero while
exact LM stayed ≈ 0.135·q), **B₅ = 2052/7⁵ ≈ +0.1221** (positive, within 9% of the
truth (6/7)¹³ ≈ 0.1348), B₇ ≈ 0.1346. Depth 5 is the first truncation that clears at
13 constraints — the same reason THM-604 needed the quintic in the continuous frame.

## 4. The iid floor (proved conditional on discrete resolution)

If no d ≤ 5 subset of classes has a small resonance mod q (discrete depth-5
resolution: no n ∈ Z^13, 0 < ‖n‖₀ ≤ 5, small coefficients, with Σ n_l v_l ≡ 0 mod q
at heights that distort the d-fold fiber counts), then each S_d(q) =
C(ncl,d)·κ_q^d·(q−1)·(1+ε_d) with κ_q = (2c−1)/q ≤ 1/7 + 2/q and explicit ε_d,
giving B₅(S,q)/q ≥ 0.1221 − ε. Same hypothesis pattern as THM-604; on Z_q the
d-fold counts are exact lattice/CRT counts (no torus-band volumes needed).

## 5. Verification (exact; the S209 battlefield and beyond)

- **The adversarial C1∪C4-killers** (klein-S209, 0 union-certs): B5 certifies
  62–100% of ALL moduli q ∈ (V, 2V]; at the best modulus B5 = LM EXACTLY (e.g.
  V=120 worst: q=231, B5 = 34 = LM; V=280: q=541, B5 = 82 vs LM = 84).
- **7-structured @91**: 90/91 moduli certified; best q=117: B5 = 38 = LM.
- **Scale sweep V = 30..260** (random mid-band covering instances): EVERY instance
  at EVERY scale has ≥ 1 certified modulus; min best-B5/q = 0.113. The empirical
  V₀ of the B5 route is < 30; kps-S115's `LRCCovering966` native_decide covers
  Vmax ≤ 18. Remaining unchecked band: Vmax ∈ (18, 30) (finite; enumeration is the
  cost, the per-set check is trivial — flagged follow-up).
- **Structure dial**: E3 = 7 keeps 56% of moduli certified; the covering equality
  extremal 2·{1..13} (E3 = 42, M = 1/14 exactly) keeps 3.8% — certified at its
  boundary witnesses. Degradation is graceful and E3-monotone, feeding the part-6
  dichotomy.
- Forensics: rare B5-failures (LM > 0, B5 ≤ 0) show binomial-vs-actual coverage
  histograms with visible high-coverage resonance inflation — single moduli,
  irrelevant to certification (one good q suffices).

## 6. What remains (the ONE named item — part of HYP-5732)

**The a-priori supply**: every covering 13-set with Vmax > V₀ has some q ∈
(Vmax, 2Vmax] that is depth-5 resolved (or directly B5-positive). Route:
(i) ruler-specific relations Σ n_l v_l = mq (m ≠ 0, ‖n‖₀ ≤ 5, small height) pin
≤ ‖n‖₁ moduli each — a divisor-counting O(V) budget against ~V available moduli;
(ii) EXACT relations (m = 0) hit every modulus — their count is the E2/E3 additive
structure of S, bounded away from the AP extremum for covering sets
(LEM-015 + `LRCSchurRigidity`: E3 = C(13,2) iff dilated interval; near-max E3 ⟹
near-dilated-interval ⟹ the dilation-invariant boundary family, handled exactly);
(iii) the small-V floor: native_decide (kps-S115 pattern) up to V₀.
No analysis remains in (i)/(iii); (ii) is the E3-budget-vs-B5-penalty bookkeeping.

## Files

`04-computation/lrc14_discrete_quintic_bonferroni_klein_S210.py` (+ `.out`) — the
one-pass histogram implementation, battlefield census, scale sweep, structure dial,
forensics. Lean next step: part 1+2 as `LRCDiscreteBonferroni.lean` (f₅ over a
coverage histogram, integers; consumes into `mreach_ge_of_pairsum_band`) — decide-
shaped, no analysis.
