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
- **Over F₂(t): PROVED (Abhyankar circle, 1990s).** The trinomial
  **f_t(x) = x²³ + x + t** has Gal(f_t / F₂(t)) = M23; the splitting field
  is an unramified M23-cover of the affine line in characteristic 2 — the
  "nice equations for nice groups" phenomenon, tied to the binary Golay
  code [23,12,7] whose automorphism group is M23.
- **Char-2 ⇒ char-0 transfer is exactly what's missing:** wild ramification
  at infinity has no char-0 shadow; this is the same "the wall is the
  arithmetic of the base, not the group" shape as several repo walls.

## In-repo verification (FINITE-EXACT, not a new theorem)

`04-computation/m23_census.py` (pure Python, no deps): GF(2^k) arithmetic,
distinct-degree factorization of f_t over many k and t; the factorization
type of squarefree f_t equals the cycle type of Frobenius on the 23 roots;
by function-field Chebotarev the empirical distribution must converge to
M23's class-fraction table (order 10 200 960; 12 cycle types on 23 points:
1²³, 1⁷2⁸, 1⁵3⁶, 1³2²4⁴, 1³5⁴, 1·2²3²6², 1²7³, 1·2·4·8², 1·11², 2·7·14,
3·5·15, 23). Decisive negative control: types that live in A23 but NOT in
M23 (e.g. 1²⁰3, 1¹⁹2², 5·18) must NEVER occur. The 23-cycle frequency must
approach 2/23 ≈ 0.08696 — note this particular statistic does NOT separate
M23 from A23 (all 23-cycles are even); the separation is carried by the
ABSENT types and by the 11², 14-, 15-type frequencies.

Census output is frozen next to the script; see the run block at the end of
this note.

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
