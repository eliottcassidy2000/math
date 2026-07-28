# M23 as a Galois group: exact status + in-repo Frobenius census (klein-S691, 2026-07-28)

**Prompt (owner puzzle bundle):** "Find a polynomial whose Galois group is
the Mathieu group M23." This question appears in Epoch AI's FrontierMath
open-problems set (https://epoch.ai/frontiermath/open-problems/inverse-galois).
The honest answer depends entirely on the base field, and the repo must not
blur that line.

## Status ledger (checked 2026-07-28)

- **Over Q: OPEN.** M23 is the unique sporadic simple group not known to be
  a Galois group over Q (all other 25 are settled; Malle–Matzat framework).
  The rigidity method fails so far: no known rationally rigid class
  vectors; the recent braid-orbit study (Barth–König-circle, arXiv:2202.08222)
  computed new braid invariants for M23 and still reports no realization.
  **Guardrail: never cite an M23/Q polynomial as existing; never claim the
  problem impossible either.**
- **Over number fields: PARTIAL.** M23 is realized over every number field
  K in which −1 is a sum of two squares (e.g. Q(i)); so the obstruction is
  genuinely arithmetic, not group-theoretic.
- **Over F₂(t): PROVED (Abhyankar, BAMS 27 (1992) / Israel J. Math 88
  (1994)).** The trinomial **f_t(x) = x²³ + t·x³ + 1** has
  Gal(f_t / F₂(t)) = M23; writing t = (x²³+1)/x³ gives dt/dx = x⁻⁴ ≠ 0, so
  this is his celebrated UNRAMIFIED M23-cover of the affine line in
  characteristic 2 — "nice equations for nice groups", tied to the binary
  Golay code [23,12,7] whose automorphism group is M23.
- **CORRECTION (census-caught, this session):** the session's first recall
  named x²³ + x + t as the M23 trinomial. The census REFUTED that: its
  factorization types over F_{2^k} are a clean A23 (218 alien types at
  89.4% of samples; exact witness: over GF(2), x²³+x+1 = (deg 2)(deg 8)
  (deg 13), and 13 is not an M23 cycle length). The wrong-recall family is
  kept in the script as the hostile control it turned out to be.
- **Char-2 ⇒ char-0 transfer is exactly what's missing:** wild ramification
  at infinity has no char-0 shadow; this is the same "the wall is the
  arithmetic of the base, not the group" shape as several repo walls.

## In-repo verification (FINITE-EXACT, not a new theorem)

`04-computation/m23_census.py` (pure Python, no deps): GF(2^k) arithmetic,
distinct-degree factorization over k = 5..13 (exhaustive t through 2¹²,
4000 random t at k=13; ≈24,300 factorizations, ~61 s); factorization type
of squarefree f_t = cycle type of Frobenius on the 23 roots; function-field
Chebotarev forces convergence to the class-fraction table. M23 (order
10 200 960) has exactly 12 cycle types on 23 points: 1²³, 1⁷2⁸, 1⁵3⁶,
1³2²4⁴, 1³5⁴, 1·2²3²6², 1²7³, 1·2·4·8², 1·11², 2·7·14, 3·5·15, 23.
Decisive negative control: A23-only types (1²⁰3, 1¹⁹2², 2·8·13, …) must
never occur. Note the 23-cycle frequency 2/23 does NOT separate M23 from
A23; separation is carried by the absent types.

**Results (frozen in the script's context block + session scratch):**
- x²³ + t·x³ + 1: **exactly the 12 M23 types, zero aliens, all 12 observed**
  (rules out proper transitive subgroups such as 23:11); pooled N = 12,160;
  best-field fit k=12 (exhaustive 4096): max deviation 0.0042, χ² ≈ 4.7 on
  ~11 dof. Small-field bias is genuine O(q^{-1/2}) arithmetic — sharpest at
  k=11 where 23 | 2¹¹−1 (that field supplied the lone totally-split 1²³
  sample). Non-squarefree t: none (unramified over A¹, as the theorem
  predicts).
- x²³ + x + t: A23 census (224 types, frequencies fitting A23's 2/z_λ with
  max |emp−pred| = 0.0129 over all 1255 partitions of 23); non-squarefree
  only at t = 0, exactly as forced by the derivative identity.
- Validation stack all green: Rabin irreducibility asserts for the 16
  hardcoded moduli, exp/log table checks, table-vs-bitpoly multiplication
  agreement, degree-1 counts vs direct root censuses (k ≤ 8), independent
  GF(2) DDF transport at t = 1.

## Why this stays in the repo

1. The Golay/M23 objects sit one street over from our tournament/code work
   (deletion-code/cogirth wheel THM-2069; Hamming-support THM-2082): a
   perfect code whose symmetry is the obstruction-free part, with the
   arithmetic base carrying the open difficulty — a useful mirror for
   LRC(14)'s "finite box ≠ discharge" guardrails.
2. If any fleet session is tempted by "M23 via specialization": the correct
   target is a REGULAR realization over Q(t); Hilbert irreducibility then
   gives infinitely many specializations. The known failure mode is that
   all candidate class vectors have nontrivial field of definition. Do not
   spend cycles on naive specialization hunts.

## Census results (klein-S691 run)

(appended after the census agent's run completes; script + raw output are
the source of truth)
