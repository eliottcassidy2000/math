---
id: THM-4916
title: "Local factors of the two-root-change multiplicity: (main, new) in the letter-count limit the class distribution of a word's value mod an odd prime p has L_p(0) = 1 + (3p-2)(p-1)/(2p^3) and L_p(r) = 1 - (3p-2)/(2p^3) for r != 0 (two root changes), and L_p(0) = 1 + (3p-4)/(4p^2), L_p(-1) = 1 - 1/p^2 + 1/(4p), other classes 1 - 1/p^2 (one root change), by a confinement cascade at the root changes; with the supporting identities (1a) seed confinement mod gcd(R,d) = gcd(D,U), (1b) root-change tails (m - Q(m))D + (m+1)U for odd m with gcd(Q(m)+1, m+1) | N, (1c) value(p J Z O^(j-1) Z) = (2^(j-1)+1)(D+2U) for j >= 2, (1d) the K-form expansion with a = D-E-U+d, (1e) the J-budget product family V V in F_<=2"
status: >
  PROVED. The main closed form is new for general p, in the letter-count model (all prefixes with n letters O/Z and
  at most two J, weighted as in M_2), as n -> infinity. The proof is the three-block confinement cascade (PROOF
  section 2.1) together with lemma (M) and the reduction to Pi J Pi J Pi delta_t. Lemma (M) and the reduction were
  proved by the independent audit (audit/AUDIT.md section 5); the reduction is Session223 section 5's argument.
  Checks:
  - exact rational equality at p = 3, 5, 7;
  - agreement within 1.2e-12 for every prime p <= 31, by two independent projector methods.
  The method and the case p = 3 (34/27, 47/54; one J 41/36, 8/9, 35/36) are Session223's (threeadic/PROOF.md,
  commit 1c524e45a9, 13:31:03), which also proves the 3-adic law at every precision.
  Priority of the supporting identities:
  - (1a) is KNOWN: Session206 weighted-parameter-capacity and exact-modular-image, THM-4888 (C);
  - (1b) is concurrent with Session223 cutoff/REPORT.md eq. (2), pushed 32 s before this packet's checkpoint. New
    here: the m, Q(m) form and the tail-confinement corollary;
  - (1c) is a one-line corollary of s.JZ = s.OJ (Session205, THM-4888 (D)1, Session199 hidden-structure Theorem 2)
    and of F((h,e,h,e).O^k) = (2^k+1)h, the mechanism of Session223's rectangle;
  - (1d) is Session199's: a = -f;
  - (1e) is Session199 physical-adjoint section 4 (K = hh' on balanced pairs). New: the J-budget form V V in F_<=2.
  NOT CLAIMED:
  - the limit of the census (height-weighted) factors, which is OPEN (pre-registered R2 FAILED);
  - any cutoff or coverage theorem. Chan-Pak OPEN.
  AUDIT: independent, PASS WITH CORRECTIONS (3 majors on priority, the cutoff paragraph and the Fermat-boost mechanism;
  14 minors; all applied; MISTAKE-734).
session: kind-pasteur-2026-11-07-S15
date: 2026-11-07
source: 05-knowledge/results/root_change_local_factors_20261107/PROOF.md
related:
  - HYP-9506
  - THM-4904
  - HYP-9494
  - HYP-9480
  - THM-4888
  - session199/physical-adjoint/PROOF.md
  - session200/carmichael-fibers/PROOF.md
  - 05-knowledge/results/session223_root_changes_and_local_factors_20261107/README.md
  - 04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py
script: 05-knowledge/results/root_change_local_factors_20261107/exact_limit.py; 05-knowledge/results/root_change_local_factors_20261107/resdp.c; 05-knowledge/results/root_change_local_factors_20261107/audit/cesaro_exact.py; 05-knowledge/results/root_change_local_factors_20261107/audit/cesaro_lazy.c
---

# THM-4916 - local factors of the two-root-change multiplicity

Proofs: [PROOF.md](../../05-knowledge/results/root_change_local_factors_20261107/PROOF.md), sections 1 and 2.1.
Audit: [AUDIT.md](../../05-knowledge/results/root_change_local_factors_20261107/audit/AUDIT.md).

## Main statement (new for general p)

Count every word over {O,Z,J} from t = (1,0,2,1) with n letters O/Z and at most K letters J. Weight it 2 if it has
fewer than K letters J and 1 otherwise, which is the weighting of M_K. Record the class mod p of its value, and let
L_n(r) be p times the weighted share of class r. As n -> infinity, for every odd prime p:

- **K = 2:** L(0) = 1 + (3p-2)(p-1)/(2p^3), and L(r) = 1 - (3p-2)/(2p^3) for r != 0. For example 34/27 and 47/54 at
  p = 3, and 151/125 and 237/250 at p = 5. (L(0) - 1)p tends to 3/2.
- **K = 1:** L(0) = 1 + (3p-4)/(4p^2), L(-1) = 1 - 1/p^2 + 1/(4p), and every other class has 1 - 1/p^2.

**Mechanism: a confinement cascade.** A long O/Z block with d != 0 mixes (D,U) to uniform. A root change whose input
has d = U = 0 (mod p) freezes the value class of the next block at 0 or at E. The confinement probability is 1/p^2 at
the first root change and (7p-4)/(4p^3) at the second.

## Supporting identities (credited)

1. **(1a)** H_s(m) = h - e(m mod 2) mod gcd(R,d); for a root-change seed, gcd(R,d) = gcd(D,U).
2. **(1b)** value(p J tail m) = (m - Q(m))D + (m+1)U for odd m, and (m+1-Q(m))D + (m+1)U - E for even m. So
   gcd(Q(m)+1, m+1) divides every odd-tail value.
3. **(1c)** value(p J Z O^(j-1) Z) = (2^(j-1)+1)(D+2U) for j >= 2.
4. **(1d)** K(v,w) = U_vU_w + E_vE_w - d_vd_w + a_vd_w + a_wd_v, with a = D-E-U+d (Session199's -f), and
   F_<=2 = {1} ∪ {K(t.alpha, t.beta) : #J(alpha) + #J(beta) <= 2}.
5. **(1e)** (x + 2^p y)(x' + 2^q y') is in F_<=2 for J-free states (x,.,y,1), (x',.,y',1) and p, q >= 1.
