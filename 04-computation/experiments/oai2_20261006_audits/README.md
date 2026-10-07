# Independent audits of the openai/math second reading (mac-mini-2026-10-06-oaimath2)

Two adversarial auditors checked the session's checkpoint `9258bf9d52` read-only, each with its own code. Their reports are `audit_A.md` and `audit_B.md`; every correction was applied before the final commit (MISTAKE-577 sharpened, MISTAKE-579).

* **Audit A** (`auditA_*`): algebra and number theory.
  * THM-4558: κ(q) by SAT, Hoffman bounds, decomposition/inertia groups of all primes ≤ 400, generalized spindles.
  * THM-4559: the idempotent certificate on Pálvölgyi's heptagon with `r = 3`.
  * THM-4561: the tower identities and bounds.
  * THM-4562: fibre counts over `GF(q)`, including `q = 25, 49, 121, 125`.
* **Audit B** (`auditB_*`): combinatorics.
  * THM-4560 / MISTAKE-577: the `n = 2` games with two solvers, plus exhaustive enumeration of all `2^24` gap sets (`auditB_gap_bruteforce`, `auditB_F2_bruteforce`), and the `n = 3` witnesses with all `1.7·10^8` subgrids.
  * Crossing counts and DDS offsets.
  * `u(21) = 57` and the 49-point set.
  * The Davenport law and 3-smooth thresholds.
  * The 14-vertex Barnette graph.
  * HYP-9215's `n = 6` one-move case, by a CEGAR hitting-set proof.

Third-party papers the auditors read (PDFs, texts) are not committed.
