---
id: THM-4916
title: "Exact identities behind the local factors of the two-root-change multiplicity: (1a) a seed's tail values are h - e (m mod 2) mod gcd(R,d), and gcd(R,d) = gcd(D,U) for a root-change seed; (1b) the value of p J (tail m) is (m - Q(m))D + (m+1)U for odd m and (m+1-Q(m))D + (m+1)U - E for even m, so gcd(Q(m)+1, m+1) divides every odd-tail value; (1c) value(p J Z O^(j-1) Z) = (2^(j-1)+1)(D + 2U) - Fermat-type multipliers 3, 5, 9, 17, 33, ...; (1d) the defect a = D - E - U + d is Z-invariant, J-reset and O-shifted by E - U, and K(v,w) = U_v U_w + E_v E_w - d_v d_w + a_v d_w + a_w d_v; (1e) F_<=2 contains the product set V V, V = {D + 2^p U : J-free states, p >= 1}"
status: >
  PROVED (new as stated): (1a)-(1e), by direct substitution into the physical moves and Session199's K-form. (1c) is
  checked on 889 words and (1e) on 3,000 random pairs, both with exact arithmetic. KNOWN and credited:
  - Session199's K-form and its source-word realisation;
  - Session200 Theorem 2, the family D + 2^(k+1) U;
  - seed confinement mod 3 as a mechanism for a mod-3 local factor: opus-S82, local_factor_mod3_prediction.py, from
    Session193 modular_sync and Session206's weighted capacities.
  Session223 (same owner prompt, concurrent) has the two-run product rectangle; (1e) is a different, two-state product
  family. The local-factor consequences are HYP-9506, with pre-registered tests R2 FAILED, R3-R5 CONFIRMED. NOT
  CLAIMED: any coverage theorem or cutoff; Chan-Pak OPEN. AUDIT: pending.
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
  - 04-computation/experiments/chanpak_mass_exploratory_20261106/local_factor_mod3_prediction.py
  - 05-knowledge/results/session223_root_changes_and_local_factors_20261107/README.md
script: 05-knowledge/results/root_change_local_factors_20261107/resdp.c; 05-knowledge/results/root_change_local_factors_20261107/decomp.c; 05-knowledge/results/root_change_local_factors_20261107/rc_lf.c
---

# THM-4916 - exact identities behind the local factors of the two-root-change multiplicity

Proofs: [PROOF.md](../../05-knowledge/results/root_change_local_factors_20261107/PROOF.md), section 1.

Words over {O,Z,J} act on t = (D,E,U,d) = (1,0,2,1) by O = (D,0,2U,d), Z = (D-E+U,d,2U-d,d) and
J = (D-E+U,D-E,D+U,D). A seed s = (h,R,e,d) has tail values H_s(m) = h + Rm - dQ(m) - e(m mod 2).

1. **(1a) Seed confinement.** H_s(m) = h - e(m mod 2) (mod gcd(R,d)). For s = J(p), gcd(R,d) = gcd(D_p, U_p).
2. **(1b) Root-change tail in pre-change coordinates.** For a state p = (D,E,U,d):
   - odd m: value(p J tail m) = (m - Q(m)) D + (m + 1) U;
   - even m: value(p J tail m) = (m + 1 - Q(m)) D + (m + 1) U - E.

   So gcd(Q(m) + 1, m + 1) divides every odd-tail value (tail confinement).
3. **(1c) Fermat-type multipliers.** For j >= 2, value(p J Z O^(j-1) Z) = (2^(j-1) + 1)(D + 2U). The multipliers are
   2^i + 1, and the Fermat primes 3, 5, 17, 257 are among them.
4. **(1d) Defect and K-form.**
   - a = D - E - U + d is fixed by Z, reset to 0 by J, and sent to a + E - U by O;
   - Session199's form satisfies K(v,w) = U_v U_w + E_v E_w - d_v d_w + a_v d_w + a_w d_v;
   - F_<=2 = {K(v,w) : #J(v) + #J(w) <= 2}.
5. **(1e) Product family.** For J-free states (x,eps,y,1) and (x',eps',y',1) and p, q >= 1,

       (x + 2^p y)(x' + 2^q y') = K((x,eps,y,1) J Z^p, (x',eps',y',1) J Z^q)   is in F_<=2.

Consequences, stated as conjectures in HYP-9506:
- (1a)-(1c) make N mod q a finite-automaton image, which gives divisor-type local factors heaviest at the
  Fermat-type primes;
- (1e) shows F_<=2 contains all products of two family-A values;
- every product or multiplier family misses all primes, which is the barrier for a cutoff proof.
