---
id: HYP-9331
title: "The Lucas-block detector finds rank-two attractor points. In a blind, pre-registered (B1') scan over 0 < |m| <= 10^4 of four Calabi-Yau threefold families, the survivors are exactly the two known rational rank-2 attractor points (Hulek-Verrill, m = -7; 4F3(1/3,1/3,2/3,2/3), m = -5832) together with the rational singular points. At ordinary primes the block defect has valuation 2 resp. 3, the moduli of Beukers' conjectured supercongruences there. Conjecture: for Lucas CY3 families, a nonsingular rational point survives (B1') iff it is a rank-2 attractor point."
status: >
  OPEN. FINITE-EXACT evidence: 4 families, 0 < |m| <= 10^4, filter primes 7..43, confirm primes 47..200.
  - Prediction P9 (pre-registered and pushed in 1846da553b before the run; outcome in ceb77fcd33) PASSED.
  - The attractor points are CITED:
    - Hulek-Verrill t_0 = -1/7 (AESZ 34): Candelas, de la Ossa, Elmi, van Straten, JHEP 10 (2020) 202,
      arXiv:1912.06146;
    - 4F3(1/3,1/3,2/3,2/3), t_0 = -1/8 in t = 729x: Boenisch, Elmi, Kashani-Poor, Klemm, "Time reversal and CP
      invariance in Calabi-Yau compactifications", JHEP 09 (2022) 019, arXiv:2204.06506, Table 1 (Beukers' [4]).
    - Both appear as Conjectures 1.27 (mod p^2) and 1.28 (mod p^3) in Beukers arXiv:2403.03301.
  - The conjectural "iff" goes beyond the data: only two attractor points and four families were tested.
  - The physics reading (BPS attractors, flux vacua) is an ANALOGY, as in S43's bridge note section 4.4; no mechanism
    is claimed.
session: kind-pasteur-2026-10-18-S2
source: 05-knowledge/results/lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md
scripts:
  - 04-computation/cm_detector_sporadic/weight3_attractor_scan.py (+ .out, .json)
  - 04-computation/cm_detector_sporadic/PREREGISTRATION.md (P9)
related: [HYP-9330, THM-4708, HYP-9310, HYP-9311, Beukers arXiv:2403.03301 Conj 1.27-1.28, Candelas-de la Ossa-Elmi-van Straten arXiv:1912.06146, Boenisch-Elmi-Kashani-Poor-Klemm arXiv:2204.06506, Moore hep-th/9807087, AESZ database arXiv:math/0507430]
---

# HYP-9331 — The Lucas-block detector finds rank-two attractor points

## Data (P9, pre-registered)

The test (B1') is that of HYP-9330: `sum_{k=p}^{2p-1} u_k m^{-k} == u_1 m^{-1} sum_{k<p} u_k m^{-k} (mod p^2)` at
every prime `7 <= p <= 200` not dividing `m`.

| family (coefficients `u_k`) | (B1') survivors, `0 < \|m\| <= 10^4` | attractor | singular points | `v_p` at ordinary primes |
|---|---|---|---|---|
| Hulek-Verrill `sum_{m_1+...+m_5 = k} (k!/prod m_i!)^2` (AESZ 34) | -7, 1, 9, 25 | -7 | 1, 9, 25 | 2 at -7 (37 of 41 primes; 3 at the other 4; one non-ordinary prime, also 2); 2 at 1, 9, 25 |
| `((3k)!/k!^3)^2` = 4F3(1/3,1/3,2/3,2/3), `t = 729x` | -5832, 729 | -5832 = -18^3 | 729 | 3 at -5832 (42 of 43 primes; 4 at one); 3 at 729 |
| `C(2k,k)^4` (control) | 256 | none known | 256 | 3 |
| `(5k)!/k!^5` (quintic, control) | 3125 | none known | 3125 | 3 |

- Generic `m = 5`: `v_p = 1` at ordinary primes, in every family.
- The modulus of Beukers' conjectured supercongruence is `p^2` at the Hulek-Verrill point (Conj. 1.27) and `p^3` at
  `t_0 = -1/8` (Conj. 1.28). These match the observed `v_p`.
- **Post hoc, not predicted.** The ordinary-restricted test (B1'_ord) of HYP-9330 also passes `m = -256` for
  `C(2k,k)^4`, i.e. `z = 256x = -1`, the second fixed point of the symmetry `z <-> 1/z` of the `(1/2)^4` operator.
  Its non-ordinary primes (`p <= 199`) are exactly `p == 3 (mod 4)`; there `v_p = 1`, and elsewhere `v_p = 2`. It
  behaves like a point with CM by `Q(i)`. (B1') rejects it, as it rejects weight-1 CM points (HYP-9330(B)).

## Conjecture

Let `L` be a Calabi-Yau threefold operator whose holomorphic period coefficients `u_k` have the Lucas property. Let
`x_0 = 1/m` be a rational nonsingular point.
1. `x_0` survives (B1') at all large `p` iff `x_0` is a rank-2 attractor point. Here "attractor" means that the
   local zeta factor of the fibre splits off a weight-2 piece over `Q`.
2. At such a point, `v_p(Delta_p)` at ordinary primes equals the depth of the supercongruence
   `F_p(x_0) == a_p(f_4) (mod p^r)`.

## Mechanism (heuristic, by analogy with HYP-9330(D))

- In weights 1 and 2, the split-prime half of (B1') is the truncated form of `F(t) == F_p(t) F(t^sigma)`, evaluated
  at a fixed point of the excellent Frobenius lift. The fixed points are the Serre-Tate canonical lifts, i.e. the CM
  points.
- For CY3 families, Beukers-Vlasenko construct excellent lifts (Dwork crystals III). The conjecture says their
  rational fixed points away from the singular locus are the rank-2 attractors. This would be a p-adic counterpart of
  the attractor mechanism (Moore 1998): attractor points are where `H^3` splits compatibly with the Hodge filtration.

## Falsifiers and next tests

- **More attractor families.** Scan further AESZ operators with known rational attractor points: Table 1 of
  Boenisch-Elmi-Kashani-Poor-Klemm, the Candelas-de la Ossa-Kuusela-McGovern Hulek-Verrill family, Elmi's thesis.
  Every known rational attractor should survive.
- **Families without attractors.** Scan families with no known attractor at moderate height. A nonsingular survivor
  there is either a new attractor candidate (to be checked by counting points mod `p`) or a counterexample.
- **Mixed points.** Points like `z = -1` for `(1/2)^4` test whether (B1'_ord) detects a wider class, namely
  local-zeta factorization with a CM factor.
